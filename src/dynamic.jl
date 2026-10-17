struct PairwiseResult
    alignment::Alignment
    score::Float64
end

@inline gap_cost(model::ScoringModel{GlobalMode},state,i,j,n,m)=
    (model.gap_open,model.gap_extend)
@inline function gap_cost(model::ScoringModel{OverlapMode},state,i,j,n,m)
    terminal=state==2 ? (j==0 || j==m) : (i==0 || i==n)
    terminal ? (model.terminal_gap,model.terminal_gap) : (model.gap_open,model.gap_extend)
end
@inline function gap_cost(model::ScoringModel{QueryMode},state,i,j,n,m)
    terminal=state!=2 && (i==0 || i==n)
    terminal ? (model.terminal_gap,model.terminal_gap) : (model.gap_open,model.gap_extend)
end

function column_dot(a,transformed,i,j)
    total=0.0
    @inbounds for k in axes(a,1)
        total+=a[k,i]*transformed[k,j]
    end
    total
end

# Exact affine DP for single sequences. For profiles, residue substitution
# scores are weighted averages; gap costs are scaled by column occupancy.
# This is an approximation to the exact affine weighted sum-of-pairs problem.
function affine_path_profiles(ap,bp,model;free_left=true,free_right=true)
    n,m=size(ap,2),size(bp,2)
    transformed=model.matrix*bp
    occupancy_a=vec(sum(ap;dims=1)); occupancy_b=vec(sum(bp;dims=1))
    # Six score vectors instead of three full floating-point matrices.
    previous=fill(-Inf,3,n+1); current=similar(previous)
    traceback=zeros(UInt8,3,n+1,m+1)
    previous[1,1]=0.0
    for i in 1:n
        op,ex=gap_cost(model,2,i,0,n,m)
        !free_left && ((op,ex)=(model.gap_open,model.gap_extend))
        choices=(previous[1,i]+op*occupancy_a[i],previous[2,i]+ex*occupancy_a[i],previous[3,i]+op*occupancy_a[i])
        k=argmax(choices); previous[2,i+1]=choices[k]; traceback[2,i+1,1]=k
    end
    for j in 1:m
        fill!(current,-Inf)
        op,ex=gap_cost(model,3,0,j,n,m)
        !free_left && ((op,ex)=(model.gap_open,model.gap_extend))
        choices=(previous[1,1]+op*occupancy_b[j],previous[2,1]+op*occupancy_b[j],previous[3,1]+ex*occupancy_b[j])
        k=argmax(choices); current[3,1]=choices[k]; traceback[3,1,j+1]=k
        @inbounds for i in 1:n
            diagonal=(previous[1,i],previous[2,i],previous[3,i])
            k=argmax(diagonal)
            current[1,i+1]=diagonal[k]+column_dot(ap,transformed,i,j)
            traceback[1,i+1,j+1]=k
            op,ex=gap_cost(model,2,i,j,n,m)
            j==m && !free_right && ((op,ex)=(model.gap_open,model.gap_extend))
            choices=(current[1,i]+op*occupancy_a[i],current[2,i]+ex*occupancy_a[i],current[3,i]+op*occupancy_a[i])
            k=argmax(choices); current[2,i+1]=choices[k]; traceback[2,i+1,j+1]=k
            op,ex=gap_cost(model,3,i,j,n,m)
            i==n && !free_right && ((op,ex)=(model.gap_open,model.gap_extend))
            choices=(previous[1,i+1]+op*occupancy_b[j],previous[2,i+1]+op*occupancy_b[j],previous[3,i+1]+ex*occupancy_b[j])
            k=argmax(choices); current[3,i+1]=choices[k]; traceback[3,i+1,j+1]=k
        end
        previous,current=current,previous
    end
    last=(previous[1,n+1],previous[2,n+1],previous[3,n+1])
    state=argmax(last); score=last[state]
    path=UInt8[]; i,j=n,m
    while i>0 || j>0
        push!(path,state)
        next=traceback[state,i+1,j+1]
        state in (1,2) && (i-=1)
        state in (1,3) && (j-=1)
        state=next
    end
    reverse!(path),score
end

function affine_path(a::Alignment,b::Alignment,model;free_left=true,free_right=true)
    affine_path_profiles(profiles(a,model),profiles(b,model),model;free_left,free_right)
end

function apply_path(a,b,path;threads=1)
    input=vcat(a.records,b.records)
    records=Vector{Sequence}(undef,length(input))
    workers=length(input)>=8 ? resolve_threads(threads) : 1
    parallel_for(length(input),workers) do k
        output=Vector{UInt8}(undef,length(path))
        position=0
        from_a=k<=length(a)
        record=input[k]
        source=codeunits(record.sequence)
        @inbounds for (column,state) in enumerate(path)
            consumes=from_a ? state in (1,2) : state in (1,3)
            if consumes
                position+=1
                output[column]=source[position]
            else
                output[column]=UInt8('-')
            end
        end
        records[k]=Sequence(record.id,String(output))
    end
    # Preserve relative within-profile weights while assigning each group its
    # number-of-sequences share in the combined profile.
    weights=vcat(a.weights .* length(a),b.weights .* length(b))
    Alignment(records;weights)
end

function anchored_path(a,b,model,anchors;threads=1)
    isempty(anchors) && return first(affine_path(a,b,model))
    path=UInt8[]; ia,ib=0,0
    n,m=alignment_width(a),alignment_width(b)
    ap,bp=profiles(a,model),profiles(b,model)
    for (index,(i,j)) in enumerate(anchors)
        piece,_=affine_path_profiles(view(ap,:,ia+1:i-1),view(bp,:,ib+1:j-1),model;
            free_left=index==1,free_right=false)
        append!(path,piece); push!(path,1)
        ia,ib=i,j
    end
    piece,_=affine_path_profiles(view(ap,:,ia+1:n),view(bp,:,ib+1:m),model;
        free_left=false,free_right=true)
    append!(path,piece)
end

"""Align two existing profiles, preserving their columns and relative sequence weights.
`anchors=:none` uses full affine DP. `anchors=:fft` forces high-identity FFT
anchor matches and aligns intervening blocks; this optional acceleration is heuristic.
Gap costs for profiles use residue occupancy and approximate the affine WSP objective.
"""
function merge_profiles(a::Alignment,b::Alignment;model=ScoringModel(),anchors=:none,
                         anchor_window=12,max_cells=100_000_000,threads=1)
    anchors in (:none,:fft) || throw(ArgumentError("anchors must be :none or :fft"))
    max_cells>0 || throw(ArgumentError("max_cells must be positive"))
    workers=resolve_threads(threads)
    validate_ids(vcat(a.records,b.records))
    validate_residues(a,model); validate_residues(b,model)
    isempty(a) && return b
    isempty(b) && return a
    cells=(Int128(alignment_width(a))+1)*(Int128(alignment_width(b))+1)
    if anchors==:fft
        selected=fft_anchors(a,b;model,window=anchor_window)
        # Guard the largest block, not the unsegmented full matrix.
        borders=vcat([(0,0)],selected,[(alignment_width(a)+1,alignment_width(b)+1)])
        largest=maximum(Int128(borders[i+1][1]-borders[i][1])*Int128(borders[i+1][2]-borders[i][2]) for i in 1:length(borders)-1)
        largest<=max_cells || throw(ArgumentError("DP block exceeds max_cells; increase limit or use smaller profiles"))
        return apply_path(a,b,anchored_path(a,b,model,selected;threads=workers);threads=workers)
    end
    cells<=max_cells || throw(ArgumentError("Alignment exceeds max_cells=$max_cells"))
    path,_=affine_path(a,b,model)
    apply_path(a,b,path;threads=workers)
end
merge_profiles(a::AbstractVector{Sequence},b::AbstractVector{Sequence};kwargs...)=merge_profiles(Alignment(a),Alignment(b);kwargs...)

"""Optimal two-sequence affine alignment. Return `.alignment` and unweighted `.score`.
Existing gaps are removed. `model.mode=:query` fits the first sequence into the second.
"""
function pairalign(a::Sequence,b::Sequence;model=ScoringModel(),max_cells=100_000_000,
                   threads=:auto)
    aa=Alignment([Sequence(a.id,replace(a.sequence,"-"=>""))])
    bb=Alignment([Sequence(b.id,replace(b.sequence,"-"=>""))])
    result=merge_profiles(aa,bb;model,max_cells,threads)
    score=pair_score(result[1].sequence,result[2].sequence,model)
    PairwiseResult(result,score)
end
const RawSequence=Union{AbstractString,BioSequences.BioSequence}
pairalign(a::RawSequence,b::RawSequence;kwargs...)=
    pairalign(Sequence("query",a),Sequence("reference",b);kwargs...)
