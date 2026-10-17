const PROTEIN_ALPHABET=collect("ARNDCQEGHILKMFPSTWYVBZX*")
const BLOSUM62=let
    lines=readlines(joinpath(@__DIR__,"..","data","BLOSUM62.tsv"))
    reduce(vcat,[permutedims(parse.(Float64,split(line)[2:end])) for line in lines[2:end]])
end

abstract type TerminalMode end
struct GlobalMode <: TerminalMode end
struct OverlapMode <: TerminalMode end
struct QueryMode <: TerminalMode end
terminal_mode(mode)=mode==:global ? GlobalMode() : mode==:overlap ? OverlapMode() : QueryMode()

"""Substitution scores and affine gap penalties.
Protein defaults: BLOSUM62, gap_open=-10, gap_extend=-1. DNA uses match=2,
mismatch=-3. A gap of length L costs gap_open+(L-1)*gap_extend.
Modes: :global (all ends penalized), :overlap (both sequences' ends free),
:query (first argument must fit inside second). terminal_gap defaults to 0.
Supply a finite square `matrix` and unique ASCII `symbols` for custom scores.
"""
struct ScoringModel{M<:TerminalMode}
    alphabet::Symbol
    symbols::Vector{Char}
    matrix::Matrix{Float64}
    indices::Vector{Int}
    gap_open::Float64
    gap_extend::Float64
    terminal_gap::Float64
    mode::Symbol
    function ScoringModel(;alphabet=:protein,matrix=nothing,symbols=nothing,
                           match=nothing,mismatch=nothing,gap_open=-10.0,gap_extend=-1.0,
                           terminal_gap=0.0,mode=:overlap)
        alphabet in (:protein,:dna) || throw(ArgumentError("alphabet must be :protein or :dna"))
        mode in (:global,:overlap,:query) || throw(ArgumentError("Invalid alignment mode"))
        all(x->isfinite(x) && x<=0,(gap_open,gap_extend,terminal_gap)) ||
            throw(ArgumentError("Gap penalties must be finite and nonpositive"))
        chars=isnothing(symbols) ? (alphabet==:protein ? copy(PROTEIN_ALPHABET) : collect("ACGTRYSWKMBDHVN")) : collect(symbols)
        length(unique(chars))==length(chars) && !isempty(chars) &&
            all(c->isascii(c) && c!='-',chars) || throw(ArgumentError("Invalid scoring alphabet"))
        if !isnothing(matrix)
            (!isnothing(match) || !isnothing(mismatch)) && throw(ArgumentError("Choose matrix or match/mismatch"))
            scores=Matrix{Float64}(matrix)
        elseif !isnothing(match) || !isnothing(mismatch) || alphabet==:dna
            ma=isnothing(match) ? 2.0 : Float64(match)
            mi=isnothing(mismatch) ? -3.0 : Float64(mismatch)
            scores=[alphabet==:dna && (a=='N' || b=='N') ? 0.0 : a==b ? ma : mi for a in chars,b in chars]
        else
            chars==PROTEIN_ALPHABET || throw(ArgumentError("Custom symbols require a matrix"))
            scores=copy(BLOSUM62)
        end
        size(scores)==(length(chars),length(chars)) && all(isfinite,scores) || throw(ArgumentError("Invalid substitution matrix"))
        isapprox(scores,transpose(scores)) || throw(ArgumentError("Substitution matrix must be symmetric"))
        indices=zeros(Int,128)
        for (i,c) in enumerate(chars)
            indices[Int(c)+1]=i
        end
        selected=terminal_mode(mode)
        new{typeof(selected)}(alphabet,chars,scores,indices,gap_open,gap_extend,terminal_gap,mode)
    end
end

function residue_index(model,c)
    i=isascii(c) ? model.indices[Int(c)+1] : 0
    i>0 || throw(ArgumentError("Unsupported $(model.alphabet) residue: $c"))
    i
end

function validate_residues(a::Alignment,model)
    for r in a.records,c in r.sequence
        c=='-' || residue_index(model,c)
    end
    nothing
end

function profiles(a::Alignment,model)
    validate_residues(a,model)
    counts=zeros(length(model.symbols),alignment_width(a))
    gap=UInt8('-')
    for (i,r) in enumerate(a.records),(j,c) in enumerate(codeunits(r.sequence))
        c==gap || (counts[byte_residue_index(model,c),j]+=a.weights[i])
    end
    counts
end

@inline function byte_residue_index(model,c::UInt8)
    i=model.indices[Int(c)+1]
    i>0 || throw(ArgumentError("Unsupported $(model.alphabet) residue: $(Char(c))"))
    i
end

function pair_score(a::AbstractString,b::AbstractString,model)
    # All Sequence values are validated ASCII, so code units are residues. Work
    # directly on them: double-gap columns are skipped without allocating a
    # compressed copy, while their relative position does not affect end gaps.
    aa,bb=codeunits(a),codeunits(b)
    length(aa)==length(bb) || throw(ArgumentError("Aligned sequences must have equal lengths"))
    gap=UInt8('-'); n=length(aa)
    alo=nextind(aa,n); ahi=0; blo=nextind(bb,n); bhi=0
    @inbounds for j in 1:n
        if aa[j]!=gap
            alo=min(alo,j); ahi=j
        end
        if bb[j]!=gap
            blo=min(blo,j); bhi=j
        end
    end
    total=0.0; previous=0
    @inbounds for j in 1:n
        ca,cb=aa[j],bb[j]
        ca==gap && cb==gap && continue
        if ca!=gap && cb!=gap
            total+=model.matrix[byte_residue_index(model,ca),byte_residue_index(model,cb)]
            previous=0
        else
            state=ca==gap ? 2 : 1
            terminal=pair_terminal(model,state,j,alo,ahi,blo,bhi)
            total+=terminal ? model.terminal_gap : state==previous ? model.gap_extend : model.gap_open
            previous=state
        end
    end
    total
end

@inline pair_terminal(::ScoringModel{GlobalMode},state,j,alo,ahi,blo,bhi)=false
@inline pair_terminal(::ScoringModel{OverlapMode},state,j,alo,ahi,blo,bhi)=
    j<alo || j>ahi || j<blo || j>bhi
@inline pair_terminal(::ScoringModel{QueryMode},state,j,alo,ahi,blo,bhi)=
    state==2 && (j<alo || j>ahi)

# Internal compatibility for callers that already hold residue vectors.
pair_score(a,b,model)=pair_score(String(a),String(b),model)

function pair_scores(a::Alignment,model;threads=1)
    validate_residues(a,model)
    workers=resolve_threads(threads)
    n=length(a); scores=PackedSymmetric{Float64}(n)
    parallel_for(max(n-1,0),workers) do i
        @inbounds for j in i+1:n
            score=pair_score(a[i].sequence,a[j].sequence,model)
            scores[i,j]=scores[j,i]=score
        end
    end
    total=0.0
    @inbounds for i in 1:n-1,j in i+1:n
        total+=a.weights[i]*a.weights[j]*scores[i,j]
    end
    scores,total
end

function cross_score(a::Alignment,left,right,model;threads=1)
    workers=resolve_threads(threads)
    partial=zeros(length(left))
    parallel_for(length(left),workers) do k
        i=left[k]; subtotal=0.0
        @inbounds for j in right
            subtotal+=a.weights[i]*a.weights[j]*pair_score(a[i].sequence,a[j].sequence,model)
        end
        partial[k]=subtotal
    end
    sum(partial)
end

function cross_score(scores::AbstractMatrix,weights,left,right)
    total=0.0
    @inbounds for i in left,j in right
        total+=weights[i]*weights[j]*scores[i,j]
    end
    total
end

"""Weighted sum-of-pairs score with affine gaps; double-gap columns are ignored.
In :query mode, each pair uses its first sequence as the query.
"""
function alignment_score(a::Alignment;model=ScoringModel(),threads=1)
    last(pair_scores(a,model;threads))
end
