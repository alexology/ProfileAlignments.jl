"""Named sequence, normalized to uppercase ASCII. Alignment gaps use `-`."""
struct Sequence
    id::String
    sequence::String
    function Sequence(id::AbstractString, sequence::AbstractString)
        isempty(id) && throw(ArgumentError("Sequence ID cannot be empty"))
        occursin(r"[\r\n]",id) && throw(ArgumentError("Sequence ID cannot contain newlines"))
        s=uppercase(replace(sequence,r"\s"=>""))
        isascii(s) || throw(ArgumentError("Sequence must be ASCII"))
        new(String(id),s)
    end
end
Sequence(id::AbstractString, sequence::BioSequences.BioSequence)=Sequence(id,String(sequence))

"""Alignment effort preset. Presets only select public algorithm parameters."""
abstract type AlignmentStrategy end
struct AccurateStrategy <: AlignmentStrategy end
struct BalancedStrategy <: AlignmentStrategy end
struct LargeStrategy <: AlignmentStrategy end
struct AutoStrategy <: AlignmentStrategy end

alignment_strategy(strategy::AlignmentStrategy)=strategy
alignment_strategy(strategy::Symbol)=strategy==:accurate ? AccurateStrategy() :
    strategy==:balanced ? BalancedStrategy() : strategy==:large ? LargeStrategy() :
    strategy==:auto ? AutoStrategy() :
    throw(ArgumentError("strategy must be :accurate, :balanced, :large or :auto"))
strategy_name(::AccurateStrategy)=:accurate
strategy_name(::BalancedStrategy)=:balanced
strategy_name(::LargeStrategy)=:large
strategy_name(::AutoStrategy)=:auto
strategy_defaults(::AccurateStrategy)=(retree=2,max_iterations=2,anchors=:none)
strategy_defaults(::BalancedStrategy)=(retree=1,max_iterations=1,anchors=:fft)
strategy_defaults(::LargeStrategy)=(retree=1,max_iterations=0,anchors=:fft)
function select_strategy(strategy::AlignmentStrategy,unique_sequences,thresholds)
    length(thresholds)==2 || throw(ArgumentError("strategy_thresholds must contain two integers"))
    accurate_max,balanced_max=thresholds
    accurate_max isa Integer && balanced_max isa Integer && 0<=accurate_max<=balanced_max ||
        throw(ArgumentError("strategy_thresholds must be ordered nonnegative integers"))
    strategy isa AutoStrategy || return strategy
    unique_sequences<=accurate_max ? AccurateStrategy() :
        unique_sequences<=balanced_max ? BalancedStrategy() : LargeStrategy()
end

"""Packed symmetric matrix with an implicit zero diagonal."""
struct PackedSymmetric{T} <: AbstractMatrix{T}
    dimension::Int
    values::Vector{T}
end
PackedSymmetric{T}(dimension::Integer) where T=
    PackedSymmetric{T}(dimension,zeros(T,dimension*(dimension-1)÷2))
Base.size(a::PackedSymmetric)=(a.dimension,a.dimension)
Base.IndexStyle(::Type{<:PackedSymmetric})=IndexCartesian()
@inline function packed_index(i,j)
    i==j && return 0
    i>j && ((i,j)=(j,i))
    (j-1)*(j-2)÷2+i
end
@inline Base.getindex(a::PackedSymmetric{T},i::Int,j::Int) where T=
    i==j ? zero(T) : @inbounds a.values[packed_index(i,j)]
@inline function Base.setindex!(a::PackedSymmetric,value,i::Int,j::Int)
    i==j || (@inbounds a.values[packed_index(i,j)]=value)
    value
end

function validate_ids(records)
    length(Set(r.id for r in records))==length(records) || throw(ArgumentError("Sequence IDs must be unique"))
    nothing
end

"""Rectangular alignment with names, order and normalized positive sequence weights.
`records` are Sequence objects. Empty alignments are supported.
"""
struct Alignment
    records::Vector{Sequence}
    weights::Vector{Float64}
    function Alignment(records::AbstractVector{Sequence};weights=ones(length(records)))
        validate_ids(records)
        rows=collect(records)
        isempty(rows) || all(r->length(r.sequence)==length(first(rows).sequence),rows) ||
            throw(ArgumentError("Aligned sequences must have equal lengths"))
        w=Float64.(weights)
        length(w)==length(rows) || throw(ArgumentError("Weight count must equal sequence count"))
        all(x->isfinite(x) && x>0,w) || throw(ArgumentError("Weights must be finite and positive"))
        if !isempty(w)
            w ./= maximum(w)
            w ./= sum(w)
            all(>(0),w) || throw(ArgumentError("Weight ratios exceed floating-point range"))
        end
        new(rows,w)
    end
end
Base.length(a::Alignment)=length(a.records)
Base.isempty(a::Alignment)=isempty(a.records)
Base.getindex(a::Alignment,i::Integer)=a.records[i]
Base.iterate(a::Alignment,args...)=iterate(a.records,args...)
alignment_width(a::Alignment)=isempty(a) ? 0 : length(first(a.records).sequence)

function resolve_threads(threads)
    threads === :auto && return Threads.nthreads()
    threads isa Integer || throw(ArgumentError("threads must be :auto or a positive integer"))
    threads > 0 || throw(ArgumentError("threads must be positive"))
    min(Int(threads),Threads.nthreads())
end

# Run independent indexed jobs with bounded task creation. The caller owns the
# output slots, so results stay deterministic even when jobs finish out of order.
function parallel_for(f,n::Integer,workers::Integer)
    n <= 0 && return nothing
    workers=min(n,workers)
    if workers <= 1
        for i in 1:n
            f(i)
        end
        return nothing
    end
    next=Threads.Atomic{Int}(1)
    @sync for _ in 1:workers
        Threads.@spawn begin
            while true
                i=Threads.atomic_add!(next,1)
                i>n && break
                f(i)
            end
        end
    end
    nothing
end

function alignment_matrix(a::Alignment)
    mat=Matrix{UInt8}(undef,length(a),alignment_width(a))
    for (i,r) in enumerate(a.records)
        copyto!(view(mat,i,:),codeunits(r.sequence))
    end
    mat
end

function trim_columns(a::Alignment)
    isempty(a) && return a
    mat=alignment_matrix(a)
    gap=UInt8('-')
    keep=[any(!=(gap),col) for col in eachcol(mat)]
    Alignment([Sequence(r.id,String(mat[i,keep])) for (i,r) in enumerate(a.records)];weights=a.weights)
end

function subset(a::Alignment,indices)
    trim_columns(Alignment(a.records[indices];weights=a.weights[indices]))
end

"""Read FASTA and preserve full headers, record order and gap columns."""
function read_fasta(path::AbstractString)
    records=Sequence[]
    id=nothing
    buffer=IOBuffer()
    open(path) do io
        for line in eachline(io)
            line=strip(line)
            isempty(line) && continue
            if startswith(line,">")
                isnothing(id) || push!(records,Sequence(id,String(take!(buffer))))
                id=strip(line[2:end])
            else
                isnothing(id) && throw(ArgumentError("FASTA sequence precedes header in $path"))
                print(buffer,line)
            end
        end
    end
    isnothing(id) || push!(records,Sequence(id,String(take!(buffer))))
    validate_ids(records)
    records
end

"""Write Sequence records or an Alignment to FASTA; preserve full identifiers."""
function write_fasta(records::AbstractVector{Sequence},path::AbstractString;width=80)
    width>0 || throw(ArgumentError("width must be positive"))
    validate_ids(records)
    mkpath(dirname(abspath(path)))
    open(path,"w") do io
        for r in records
            println(io,">",r.id)
            for i in 1:width:length(r.sequence)
                println(io,r.sequence[i:min(i+width-1,end)])
            end
        end
    end
    path
end
write_fasta(a::Alignment,path::AbstractString;kwargs...)=write_fasta(a.records,path;kwargs...)
