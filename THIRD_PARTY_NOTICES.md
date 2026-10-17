# Statistical and biochemical data

`data/BLOSUM62.tsv` contains the BLOSUM62 statistical substitution matrix
exported from the Biostrings data set during development. It includes
Biostrings' ambiguity-symbol conventions for B, Z and X. No Biostrings
implementation code is included or loaded at runtime.

Biostrings source and documentation:
https://bioconductor.org/packages/Biostrings/

Matrix reference: Henikoff S, Henikoff JG (1992), Amino acid substitution
matrices from protein blocks. PNAS 89:10915–10919.
https://doi.org/10.1073/pnas.89.22.10915

Volume/polarity constants in `src/fft.jl` are biochemical data from Grantham
(1974), Amino acid difference formula to help explain protein evolution.
https://doi.org/10.1126/science.185.4154.862

# Algorithm references

The implementation is original Julia code inspired by published algorithms;
no MAFFT code has been copied. Architecture references:

- Katoh et al. (2002), https://doi.org/10.1093/nar/gkf436
- Katoh et al. (2005), https://doi.org/10.1093/nar/gki198
- MAFFT algorithm descriptions:
  https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html

Default gap penalties and profile occupancy approximations are this package's
own choices and do not reproduce MAFFT's normalized substitution scores,
gap treatment or FFT segment selection.
