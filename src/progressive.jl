"""Binary guide tree; leaves index the input sequences. Heights are UPGMA distances/2."""
struct GuideTree
    left::Union{Nothing,GuideTree}
    right::Union{Nothing,GuideTree}
    leaves::Vector{Int}
    height::Float64
end
GuideTree(i::Integer)=GuideTree(nothing,nothing,[Int(i)],0.0)
isleaf(t::GuideTree)=isnothing(t.left)

const DAYHOFF_CODES=let
    codes=zeros(UInt8,128)
    for (group,letters) in enumerate(("AGPST","DENQ","HKR","ILMV","FWY","C")),c in letters
        codes[Int(c)+1]=group
    end
    codes
end
const DNA_CODES=let
    codes=zeros(UInt8,128)
    for (group,c) in enumerate("ACGT")
        codes[Int(c)+1]=group
    end
    codes
end

@inline kmer_code(c::UInt8,alphabet)=alphabet==:protein ? DAYHOFF_CODES[Int(c)+1] : DNA_CODES[Int(c)+1]

function encoded_kmer_counts(s,k,alphabet)
    codes=codeunits(s)
    base=UInt64(alphabet==:protein ? 7 : 5)
    high=base^(k-1)
    counts=Dict{UInt64,Int}()
    key=UInt64(0); run=0
    @inbounds for i in eachindex(codes)
        value=UInt64(kmer_code(codes[i],alphabet))
        if iszero(value)
            key=0; run=0
        elseif run<k
            key=key*base+value
            run+=1
            run==k && (counts[key]=get(counts,key,0)+1)
        else
            old=UInt64(kmer_code(codes[i-k],alphabet))
            key=(key-old*high)*base+value
            counts[key]=get(counts,key,0)+1
        end
    end
    counts
end

function tuple_kmer_counts(s,k,alphabet)
    values=alphabet==:protein ? [Int(kmer_code(c,alphabet)) for c in codeunits(s)] :
        [Int(kmer_code(c,alphabet)) for c in codeunits(s)]
    counts=Dict{Tuple,Int}()
    for i in 1:length(values)-k+1
        word=Tuple(values[i:i+k-1])
        any(iszero,word) && continue
        counts[word]=get(counts,word,0)+1
    end
    counts
end

function distances_from_kmers(raw,counts,length_correction,workers)
    totals=[sum(values(c);init=0) for c in counts]
    n=length(raw); distances=zeros(n,n)
    parallel_for(max(n-1,0),workers) do i
        for j in i+1:n
            shared=sum(min(value,get(counts[j],word,0)) for (word,value) in counts[i];init=0)
            denom=min(totals[i],totals[j])
            d=raw[i]==raw[j] ? 0.0 : denom==0 ? 1.0 : 1-shared/denom
            if length_correction && denom>0
                x,y=max(length(raw[i]),length(raw[j])),min(length(raw[i]),length(raw[j]))
                d/=0.1*y/x+10000/(x+10000)+0.01
            end
            distances[i,j]=distances[j,i]=d
        end
    end
    distances
end

"""Initial distances from shared k-mer counts (Dayhoff six-group protein alphabet).
Use MAFFT's documented length correction by default; this is an approximate
similarity distance, not an evolutionary substitution model.
"""
function kmer_distances(records::AbstractVector{Sequence};alphabet=:protein,k=6,
                        length_correction=true,threads=1)
    alphabet in (:protein,:dna) || throw(ArgumentError("Invalid alphabet"))
    k>0 || throw(ArgumentError("k must be positive"))
    workers=resolve_threads(threads)
    raw=[replace(r.sequence,"-"=>"") for r in records]
    effective=isempty(raw) ? k : min(k,max(1,minimum(length.(raw))))
    # Base-7 protein words fit in UInt64 through k=22; DNA base-5 through k=27.
    encoded_limit=alphabet==:protein ? 22 : 27
    if effective<=encoded_limit
        counts=Vector{Dict{UInt64,Int}}(undef,length(raw))
        parallel_for(length(raw),workers) do i
            counts[i]=encoded_kmer_counts(raw[i],effective,alphabet)
        end
        distances_from_kmers(raw,counts,length_correction,workers)
    else
        counts=Vector{Dict{Tuple,Int}}(undef,length(raw))
        parallel_for(length(raw),workers) do i
            counts[i]=tuple_kmer_counts(raw[i],effective,alphabet)
        end
        distances_from_kmers(raw,counts,length_correction,workers)
    end
end

@inline function heap_less(a,b)
    a[1]<b[1] || (a[1]==b[1] && (a[2]<b[2] || (a[2]==b[2] && a[3]<b[3])))
end

function heap_push!(heap,item)
    push!(heap,item); child=length(heap)
    while child>1
        parent=child>>>1
        heap_less(heap[child],heap[parent]) || break
        heap[child],heap[parent]=heap[parent],heap[child]
        child=parent
    end
    heap
end

function heap_pop!(heap)
    root=first(heap); tail=pop!(heap)
    isempty(heap) && return root
    heap[1]=tail; parent=1
    while true
        left=parent<<1
        left>length(heap) && break
        right=left+1
        child=right<=length(heap) && heap_less(heap[right],heap[left]) ? right : left
        heap_less(heap[child],heap[parent]) || break
        heap[child],heap[parent]=heap[parent],heap[child]
        parent=child
    end
    root
end

"""Construct a deterministic UPGMA guide tree from a symmetric distance matrix.
`linkage=:mixed` blends minimum and size-weighted mean distances with
`mix` mean weight (default .1); `:upgma` uses the weighted mean alone.
"""
function guide_tree(distances::AbstractMatrix;linkage=:upgma,mix=0.1)
    n=size(distances,1)
    n>0 && size(distances,2)==n || throw(ArgumentError("Distance matrix must be nonempty and square"))
    all(x->isfinite(x) && x>=0,distances) && isapprox(distances,transpose(distances)) &&
        all(i->abs(distances[i,i])<1e-12,1:n) || throw(ArgumentError("Distances must be symmetric, finite, nonnegative with zero diagonal"))
    linkage in (:upgma,:mixed) && 0<=mix<=1 || throw(ArgumentError("Invalid linkage"))
    work=PackedSymmetric{Float64}(2*n-1)
    for a in 1:n-1,b in a+1:n
        work[a,b]=distances[a,b]
    end
    nodes=GuideTree[GuideTree(i) for i in 1:n]
    active=collect(1:n); sizes=ones(Int,2*n-1)
    alive=falses(2*n-1); alive[1:n].=true
    heap=Tuple{Float64,Int,Int}[]
    sizehint!(heap,n*(n-1)÷2)
    for a in 1:n-1,b in a+1:n
        heap_push!(heap,(work[a,b],a,b))
    end
    while length(active)>1
        best=heap_pop!(heap)
        while !alive[best[2]] || !alive[best[3]]
            best=heap_pop!(heap)
        end
        distance,a,b=best
        height=max(distance/2,nodes[a].height,nodes[b].height)
        push!(nodes,GuideTree(nodes[a],nodes[b],vcat(nodes[a].leaves,nodes[b].leaves),height))
        c=length(nodes); sizes[c]=sizes[a]+sizes[b]; alive[c]=true
        for d in active
            d in (a,b) && continue
            average=(sizes[a]/sizes[c])*work[a,d]+(sizes[b]/sizes[c])*work[b,d]
            value=linkage==:upgma ? average : mix*average+(1-mix)*min(work[a,d],work[b,d])
            work[c,d]=work[d,c]=value
            heap_push!(heap,(value,min(c,d),max(c,d)))
        end
        alive[a]=alive[b]=false
        filter!(d->d!=a && d!=b,active); push!(active,c)
    end
    last(nodes)
end

function tree_weights(tree,n)
    weights=zeros(n)
    function visit(node,parentheight)
        branch=max(parentheight-node.height,1e-8)/length(node.leaves)
        weights[node.leaves] .+= branch
        if !isleaf(node)
            visit(node.left,node.height); visit(node.right,node.height)
        end
    end
    visit(tree,tree.height)
    weights ./ sum(weights)
end

function reorder(a,ids,weights=nothing)
    lookup=Dict(r.id=>i for (i,r) in enumerate(a.records))
    indices=[lookup[id] for id in ids]
    Alignment(a.records[indices];weights=isnothing(weights) ? a.weights[indices] : weights)
end

function progressive(records,tree,weights;thread_count=1,kwargs...)
    if isleaf(tree)
        return Alignment([records[only(tree.leaves)]])
    end
    if thread_count>1
        if isleaf(tree.left)
            a=progressive(records,tree.left,weights;thread_count=1,kwargs...)
            b=progressive(records,tree.right,weights;thread_count=thread_count,kwargs...)
        elseif isleaf(tree.right)
            a=progressive(records,tree.left,weights;thread_count=thread_count,kwargs...)
            b=progressive(records,tree.right,weights;thread_count=1,kwargs...)
        else
            # Internal-node counts approximate the number of pending profile
            # merges. This keeps threads on the larger side of skewed trees.
            left_work=length(tree.left.leaves)-1
            right_work=length(tree.right.leaves)-1
            left_threads=clamp(round(Int,thread_count*left_work/(left_work+right_work)),
                1,thread_count-1)
            right_threads=thread_count-left_threads
            left_task=Threads.@spawn progressive(records,tree.left,weights;
                thread_count=left_threads,kwargs...)
            b=progressive(records,tree.right,weights;thread_count=right_threads,kwargs...)
            a=fetch(left_task)
        end
    else
        a=progressive(records,tree.left,weights;thread_count=1,kwargs...)
        b=progressive(records,tree.right,weights;thread_count=1,kwargs...)
    end
    a=Alignment(a.records;weights=weights[tree.left.leaves])
    b=Alignment(b.records;weights=weights[tree.right.leaves])
    result=merge_profiles(a,b;threads=thread_count,kwargs...)
    # Child profile records follow tree leaf order, so absolute tree weights
    # can be restored at every merge without changing their relative values.
    Alignment(result.records;weights=weights[tree.leaves])
end

function aligned_distances(a;threads=1)
    workers=resolve_threads(threads)
    n=length(a); matrix=zeros(n,n)
    rows=codeunits.(getfield.(a.records,:sequence))
    gap=UInt8('-')
    spans=[(something(findfirst(!=(gap),row),1),something(findlast(!=(gap),row),0)) for row in rows]
    parallel_for(max(n-1,0),workers) do i
        for j in i+1:n
            ar,br=spans[i],spans[j]
            distance=1.0
            if first(ar)<=last(ar) && first(br)<=last(br)
                lo,hi=max(first(ar),first(br)),min(last(ar),last(br))
                compared=different=0
                @inbounds for k in lo:hi
                    rows[i][k]==gap && rows[j][k]==gap && continue
                    compared+=1
                    different+=rows[i][k]!=rows[j][k]
                end
                compared>0 && (distance=different/compared)
            end
            matrix[i,j]=matrix[j,i]=distance
        end
    end
    matrix
end

function tree_partitions(tree,n)
    partitions=Vector{Int}[]
    function visit(node)
        length(node.leaves)<n && push!(partitions,sort(node.leaves))
        if !isleaf(node)
            visit(node.left); visit(node.right)
        end
    end
    visit(tree)
    # A cut and its complement describe the same bipartition.
    seen=Set{Tuple}(); out=Vector{Int}[]
    for p in partitions
        canonical=1 in p ? p : setdiff(1:n,p)
        key=Tuple(canonical)
        key in seen && continue
        push!(seen,key); push!(out,p)
    end
    out
end

"""Refine an MSA along guide-tree bipartitions, removing group-only gap columns.
Accept realignments only when the exact affine weighted sum-of-pairs score
increases. Return an Alignment, or diagnostics with `return_trace=true`.
The same fixed weights are used for every candidate and score comparison.
"""
function refine_alignment(a::Alignment;model=ScoringModel(),tree=nothing,max_iterations=2,
                          anchors=:none,anchor_window=12,max_cells=100_000_000,
                          threads=:auto,return_trace=false)
    max_iterations>=0 || throw(ArgumentError("max_iterations must be nonnegative"))
    model.mode==:query && throw(ArgumentError(":query is for pairwise/profile insertion, not MSA refinement"))
    workers=resolve_threads(threads)
    current=trim_columns(a)
    ids=getfield.(a.records,:id); weights=copy(a.weights)
    pair_cache,score=pair_scores(current,model;threads=workers)
    scores=[score]; accepted=0; cycles=0
    if length(a)>1 && max_iterations>0
        isnothing(tree) && (tree=guide_tree(aligned_distances(current;threads=workers)))
        sort(tree.leaves)==collect(1:length(a)) || throw(ArgumentError("Guide-tree leaves must index this alignment"))
        partitions=[(left,setdiff(1:length(a),left)) for left in tree_partitions(tree,length(a))]
        for cycle in 1:max_iterations
            cycles=cycle; changed=false
            for (left,right) in partitions
                candidate=merge_profiles(subset(current,left),subset(current,right);
                    model,anchors,anchor_window,max_cells,threads=workers)
                candidate=reorder(candidate,ids,weights)
                candidate_score=score-cross_score(pair_cache,weights,left,right)+
                    cross_score(candidate,left,right,model;threads=workers)
                if candidate_score>score+1e-9
                    current=candidate
                    # Rebuild only after an accepted move. This also records an
                    # independently summed exact score instead of an accumulated delta.
                    pair_cache,score=pair_scores(current,model;threads=workers)
                    changed=true; accepted+=1
                    push!(scores,score)
                end
            end
            changed || break
        end
    end
    return_trace ? (alignment=current,scores=scores,accepted=accepted,cycles=cycles) : current
end

"""Multiple sequence alignment entirely in Julia.
Pipeline: duplicate collapse → reduced-alphabet k-mer distances → UPGMA tree
→ progressive affine profile merges → optional guide-tree rebuilding →
tree-partition iterative refinement. Original IDs/order/duplicates are restored.
`strategy=:accurate` defaults to two tree passes, two refinement cycles and full DP.
`:balanced` uses one tree pass, one refinement cycle and FFT anchors; `:large`
uses one tree pass, no refinement and FFT anchors. Explicit `retree`,
`max_iterations` or `anchors` override the selected preset. `:auto` selects
these presets from the number of distinct sequences and `strategy_thresholds`.
`anchors=:fft` enables heuristic FFT anchor acceleration. `return_trace=true`
returns the alignment plus guide tree, scores and refinement diagnostics.
"""
function align_sequences(records::AbstractVector{Sequence};model=ScoringModel(),kmer=6,
                          strategy=:accurate,retree=nothing,max_iterations=nothing,
                          linkage=:upgma,mix=0.1,sequence_weighting=:tree,anchors=nothing,anchor_window=12,
                          strategy_thresholds=(100,500),max_cells=100_000_000,
                          threads=:auto,return_trace=false)
    requested_strategy=alignment_strategy(strategy)
    sequence_weighting in (:tree,:uniform) || throw(ArgumentError("Invalid sequence_weighting"))
    model.mode==:query && throw(ArgumentError(":query is for pairwise/profile insertion"))
    workers=resolve_threads(threads)
    validate_ids(records)
    raw=[Sequence(r.id,replace(r.sequence,"-"=>"")) for r in records]
    for r in raw,c in r.sequence
        residue_index(model,c)
    end
    reps=unique(r->r.sequence,raw)
    selected_strategy=select_strategy(requested_strategy,length(reps),strategy_thresholds)
    defaults=strategy_defaults(selected_strategy)
    retree=isnothing(retree) ? defaults.retree : retree
    max_iterations=isnothing(max_iterations) ? defaults.max_iterations : max_iterations
    anchors=isnothing(anchors) ? defaults.anchors : anchors
    retree>=1 && max_iterations>=0 || throw(ArgumentError("Invalid iteration counts"))
    if isempty(raw)
        empty=Alignment(Sequence[])
        return return_trace ? (alignment=empty,tree=nothing,representatives=String[],scores=[0.0],accepted=0,cycles=0,strategy=strategy_name(selected_strategy)) : empty
    end
    tree=guide_tree(kmer_distances(reps;alphabet=model.alphabet,k=kmer,threads=workers);linkage,mix)
    weights=sequence_weighting==:tree ? tree_weights(tree,length(reps)) : fill(1/length(reps),length(reps))
    ids=getfield.(reps,:id)
    current=reorder(progressive(reps,tree,weights;model,anchors,anchor_window,max_cells,
        thread_count=workers),ids,weights)
    scores=[alignment_score(current;model,threads=workers)]
    for pass in 2:retree
        tree=guide_tree(aligned_distances(current;threads=workers);linkage,mix)
        candidate=reorder(progressive(reps,tree,weights;model,anchors,anchor_window,max_cells,
            thread_count=workers),ids,weights)
        score=alignment_score(candidate;model,threads=workers)
        if score>last(scores)+1e-9
            current=candidate; push!(scores,score)
        end
    end
    refined=refine_alignment(current;model,tree,max_iterations,anchors,anchor_window,max_cells,
        threads=workers,return_trace=true)
    append!(scores,refined.scores[2:end])
    lookup=Dict(reps[i].sequence=>(refined.alignment[i].sequence,weights[i]) for i in eachindex(reps))
    multiplicity=Dict{String,Int}()
    for r in raw
        multiplicity[r.sequence]=get(multiplicity,r.sequence,0)+1
    end
    result=Alignment([Sequence(r.id,lookup[r.sequence][1]) for r in raw];
        weights=[lookup[r.sequence][2]/multiplicity[r.sequence] for r in raw])
    # Restoring duplicates adds only within-representative pairs. Compute that
    # constant from multiplicities instead of rescoring every expanded pair.
    duplicate_constant=0.0
    for (i,r) in enumerate(refined.alignment.records)
        count=multiplicity[replace(r.sequence,"-"=>"")]
        if count>1
            duplicate_constant+=refined.alignment.weights[i]^2*(count-1)/(2count)*
                pair_score(r.sequence,r.sequence,model)
        end
    end
    scores .+= duplicate_constant
    return_trace ? (alignment=result,tree=tree,representatives=ids,scores=scores,
        accepted=refined.accepted,cycles=refined.cycles,strategy=strategy_name(selected_strategy)) : result
end
align_sequences(strings::AbstractVector{<:AbstractString};kwargs...)=
    align_sequences([Sequence("seq_$i",s) for (i,s) in enumerate(strings)];kwargs...)
align_sequences(sequences::AbstractVector{<:BioSequences.BioSequence};kwargs...)=
    align_sequences([Sequence("seq_$i",s) for (i,s) in enumerate(sequences)];kwargs...)

"""Align independent sequence sets with a bounded Julia task pool.
Multiple sets use outer parallelism and one thread inside each alignment,
avoiding nested oversubscription. A single set receives the full thread budget.
Result order matches input order.
"""
function align_sequence_sets(groups::AbstractVector;threads=:auto,kwargs...)
    workers=resolve_threads(threads)
    results=Vector{Alignment}(undef,length(groups))
    if length(groups)==1
        results[1]=align_sequences(groups[1];threads=workers,kwargs...)
    else
        parallel_for(length(groups),workers) do i
            results[i]=align_sequences(groups[i];threads=1,kwargs...)
        end
    end
    results
end
