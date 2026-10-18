# ProfileAlignments.jl

Standalone package for sequence and profile alignment, written in Julia.
Julia ≥1.10; it uses BioSequences.jl for interoperability with biological
types, with no external FFT library, MAFFT executable, R or remote service.
AASV uses it as its alignment engine.

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/ProfileAlignments.jl")
using ProfileAlignments, BioSequences

records = [Sequence("reference", aa"ACDEFGHIKLMNPQ"),
           Sequence("fragment", aa"DEFGHIK"),
           Sequence("insertion", aa"ACDEFGWHIKLMNPQ")]
result = align_sequences(records)
write_fasta(result, "aligned.fasta")
```

From this directory, the package can be used and tested on its own:

```powershell
julia --project=. examples/protein_alignment.jl
julia --project=. -e 'using Pkg; Pkg.test()'
julia -t auto --project=. benchmark/threading.jl
```

## Algorithm

The design follows the family of progressive and iterative strategies
described in the [MAFFT algorithms](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html).
It is a new implementation, with its own choices and parameters.

1. Collapses identical sequences, keeping the original names.
2. Estimates distances from k-mer counts in the protein alphabet reduced to six
   Dayhoff groups, using the length correction documented by MAFFT.
3. Builds a deterministic UPGMA guide tree. `linkage=:mixed` is also
   available, mixing minimum and average distance (`mix=0.1`).
4. Aligns the profiles progressively along the tree. Scores are weighted
   BLOSUM62 averages and gap costs are affine.
5. Optionally rebuilds the tree from the distances in the MSA and repeats the
   alignment, keeping the new MSA if it improves the score.
6. Refines the result by splitting the groups at the tree branches and
   realigning them. Only improvements of the exact weighted sum-of-pairs score
   are accepted.
7. Restores duplicates, identifiers and the original order. Duplicates share
   the weight of their representative, so identical copies are not overweighted.

`sequence_weighting=:tree` uses weights derived from branch lengths divided
among the descendant leaves; `:uniform` gives every representative the same
weight. Weights stay fixed while candidates are compared and during refinement.

## API

| Function / type | Use |
|---|---|
| `Sequence(id, sequence)` | sequence with an identifier |
| `Alignment(records; weights=...)` | rectangular profile with normalised positive weights |
| `ScoringModel(...)` | matrix, penalties and end handling |
| `pairalign(a, b; model=...)` | affine alignment of two sequences, result `.alignment` and `.score` |
| `merge_profiles(a, b; model=...)` | merge of two profiles preserving their columns |
| `align_sequences(records; ...)` | progressive and refined MSA |
| `align_sequence_sets(groups; ...)` | parallel alignment of independent sets |
| `refine_alignment(alignment; ...)` | refinement of an existing MSA |
| `alignment_score(alignment; model=...)` | weighted affine sum-of-pairs score |
| `kmer_distances(records; ...)` | approximate initial distances |
| `guide_tree(distances; ...)` | guide tree |
| `fft_lags(a, b; ...)`, `fft_anchors(a, b; ...)` | candidate diagonals and anchor chain |
| `read_fasta`, `write_fasta`, `align_fasta` | FASTA input/output and complete FASTA workflow |

`Sequence`, `pairalign` and `align_sequences` accept both strings and
`BioSequence` values, including `LongAA`, `LongDNA{2}` and `LongDNA{4}`. At the
engine boundary, symbols are converted once into the ASCII representation used
by the optimised DP loops; results therefore keep the historical string API.
The BioSequences alphabet also recognises `O`, `U` and `J`, but the default
BLOSUM62 model still rejects them because it has no defined scores for them.

```julia
# Short sequence inserted into a reference: only the query's terminal gaps
# are free. The query sequence is always the first argument.
model = ScoringModel(mode=:query)
pair = pairalign("DEFGHIK", "ACDEFGHIKLMNPQ"; model)

# DNA with scalar match/mismatch, no automatic reverse complement.
dna = align_sequences(["ACGTACGT", "CGTAC", "ACGTTACGT"];
                      model=ScoringModel(alphabet=:dna))

# Trace of the scores of accepted candidates only.
trace = align_sequences(records; retree=2, max_iterations=5,
                        return_trace=true)
trace.scores
trace.alignment
```

### Computation strategies

`align_sequences` accepts explicit presets. The default remains `:accurate`,
so upgrading does not change the algorithm automatically:

```julia
accurate = align_sequences(records; strategy=:accurate)
balanced = align_sequences(records; strategy=:balanced)
large    = align_sequences(records; strategy=:large)
automatic = align_sequences(records; strategy=:auto)
```

- `:accurate`: two tree passes, two refinement cycles, full DP;
- `:balanced`: one tree, one refinement cycle, FFT anchors;
- `:large`: one tree, no global refinement, FFT anchors;
- `:auto`: picks among the three presets at 100 and 500 distinct sequences; the
  thresholds are changed with `strategy_thresholds=(accurate_max, balanced_max)`.

`retree`, `max_iterations` and `anchors`, when given, take precedence over the
preset. The FFT strategies are heuristic and may produce a different alignment.
The types `AccurateStrategy`, `BalancedStrategy` and `LargeStrategy` select the
same behaviour through dispatch.

For many independent sets, `align_sequence_sets` distributes the sets across
Julia tasks while preserving the order of the results. With several sets it
assigns one thread to each alignment, avoiding nested parallelism; with a
single set it uses the whole thread budget.

By default the computation uses the available Julia threads
(`threads=:auto`). Start Julia with `julia -t auto` or set `JULIA_NUM_THREADS`
to enable more than one. `threads=1` forces serial execution; a positive
integer limits the number of threads used by the package. K-mer distances,
distances derived from the MSA and pairwise scores are distributed
dynamically. The two independent branches of the guide tree are aligned in
parallel with a bounded task budget. The order of reductions and refinement
decisions stays fixed, so the same input produces the same alignment and the
same scores with `threads=1` and `threads=:auto`.

The gain depends on the shape of the tree and on the length and number of the
sequences. The last merges of the progressive stage have inherently serial
dependencies, so do not expect a speed-up proportional to the number of cores.
`benchmark/threading.jl [number_of_sequences] [length]` compares the two modes
on the current machine without additional benchmark dependencies.
`benchmark/scaling.jl [sizes...]` compares the three presets and reports time,
allocations and final width.

## Representation and memory

The public API keeps `Sequence` and FASTA as ASCII strings. The internal
kernels read the `UInt8` `codeunits` directly, without copying the residues
into `Char` matrices. Outputs and filtered columns are built in byte buffers of
exact size.

Common k-mers are encoded without collisions in `UInt64` (base 7 for the six
Dayhoff groups, base 5 for DNA); values of `k` that are too large use the
generic tuple path. UPGMA construction and the score cache use packed
symmetric matrices with an implicit diagonal. The DP end modes are type
parameters of `ScoringModel`, so Julia compiles separate kernels for global,
overlap and query alignment. The radix-2 FFT remains pure Julia and reuses
buffers across channels.

`trace.tree` indexes the unique representatives listed in
`trace.representatives`; the result keeps all the original sequences.
`trace.scores` reports the scores of the complete result, including the
constant term from the duplicates.

For `ScoringModel`, the modes are `:global` (penalises all ends), `:overlap`
(free ends on both sides, the default) and `:query`. The latter applies to
`pairalign` and `merge_profiles`, not to MSA/refinement.
A gap of length L costs **gap_open + (L−1) × gap_extend**. Protein defaults:
BLOSUM62, `gap_open=-10`, `gap_extend=-1`, `terminal_gap=0`. DNA:
`match=2`, `mismatch=-3`; `N` scores zero. Other DNA ambiguity codes follow
the scalar scores; they are not expanded into IUPAC probabilities.
Custom `matrix` and `symbols` can be passed; the matrix must be finite and
symmetric. Proteins: the 20 standard residues, `B`, `Z`, `X` and `*`.

## Optional FFT

```julia
fast = align_sequences(records; anchors=:fft, anchor_window=12)
```

A radix-2 Julia FFT computes the correlations of the normalised
volume/polarity properties for proteins, following the principle of the
[2002 MAFFT paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC135756/).
For DNA, four indicator channels are used. High-identity windows are searched
between the correlation peaks; a monotone chain of their centres splits the
matrix into smaller blocks. When there are no anchors, full DP is used.

This acceleration **forces the anchor matches** and may miss better
alignments. It is disabled by default. It does not reproduce MAFFT's segment
selection method or its exact scoring system. The parameters are not
interchangeable with `--op`, `--ep` and `--maxiterate`.

## Correctness and limitations

- Pairwise affine alignment is checked against an independent enumeration of
  all alignments of small inputs, in all three end modes. The suite also checks
  the FFT against direct transforms and correlations, profile columns,
  duplicates, determinism and monotonicity of refinement.
- Profile-to-profile DP uses the product of weighted frequencies for
  substitutions and column occupancy for gaps. This approximates the affine
  sum-of-pairs problem: it does not guarantee the best global profile
  alignment. Refinement instead uses the actual affine pair score, ignoring
  columns in which both sequences have a gap.
- With full DP, time is O(L₁×L₂) and traceback memory is O(L₁×L₂);
  intermediate scores use O(L₁) vectors, in addition to the profiles.
  `max_cells=100_000_000` caps the cells of the comparison, or of the largest
  block when anchors are used. Tree and distances require O(N²) memory;
  the distance queue makes UPGMA construction O(N² log N). The cost of
  refinement depends on the bipartitions and on the pairs crossing them;
  `max_iterations=0` disables it.
- It does not implement PartTree, L/G/E-INS-i consistency scores, SIMD or
  structural RNA alignment. It is not a complete clone of MAFFT, and neither
  biological nor performance equivalence with MAFFT has been demonstrated.
- The current checks are algorithmic and use synthetic data. Before final
  scientific use, benchmarks on real alignments and an assessment of how
  alignment choices affect AASV are needed.

Licensed under GPL ≥2, text in `COPYING`. Independent implementation based on
published descriptions; no MAFFT source code was copied. The BLOSUM62 matrix is
the amino acid substitution data from the Biostrings version used in the R
comparison, in the order ARNDCQEGHILKMFPSTWYVBZX*. `B`, `Z` and `X` therefore
follow the Biostrings convention; `X` against any residue scores −1. See
`THIRD_PARTY_NOTICES.md` for provenance and references.
