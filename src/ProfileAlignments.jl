module ProfileAlignments

import BioSequences

export Sequence, Alignment, ScoringModel, GuideTree, AlignmentStrategy,
    AccurateStrategy, BalancedStrategy, LargeStrategy, AutoStrategy, read_fasta, write_fasta,
    pairalign, merge_profiles, align_sequences, alignment_score, guide_tree,
    align_sequence_sets, kmer_distances, refine_alignment, fft_lags, fft_anchors, align_fasta

include("types.jl")
include("scoring.jl")
include("fft.jl")
include("dynamic.jl")
include("progressive.jl")

"""Read input FASTA, align in Julia, and write output FASTA. Return the Alignment."""
function align_fasta(input::AbstractString, output::AbstractString; kwargs...)
    result=align_sequences(read_fasta(input);kwargs...)
    write_fasta(result,output)
    result
end

end
