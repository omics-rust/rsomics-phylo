# rsomics-phylo

`rsomics-phylo` is one product for alignment trimming, sequence distances,
distance-based tree inference, tree comparison, tree measures, and the Hommola
cospeciation test. It replaces operation-sized crates with one installable
binary and one library built around shared alignment, matrix, and tree models.

The 0.2 source adds WPGMA. The registry still provides 0.1.0; that immutable
archive does not contain WPGMA. Native distribution is prepared separately.

```console
cargo install rsomics-phylo
rsomics-phylo trim -i aligned.fa -o trimmed.fa --maximum-gap-fraction 0.2
rsomics-phylo distance -i aligned.fa -o distances.tsv --metric k2p --threads 8
rsomics-phylo infer -i distances.tsv -o tree.nwk --method nj
rsomics-phylo compare \
  --reference expected.nwk \
  --comparison observed.nwk \
  --metric wrf \
  --output comparison.tsv
```

The command surface is:

- `trim`: remove columns whose gap fraction exceeds a threshold;
- `distance`: literal Hamming, JC69, or K2P sequence distances;
- `infer`: neighbor joining, UPGMA, or WPGMA from a labeled distance matrix;
- `compare`: RF, weighted RF, Kuhner-Felsenstein branch score, directional
  false-positive/false-negative counts, or cophenetic distance;
- `tip-distance`: all named-tip path distances by branch length or edge count;
- `summarize`: selected balance and timing measures, node ages, and
  coalescence ages;
- `cospeciation`: observed Hommola correlation and an optional permutation
  test with an optional null table.

FASTA records retain their complete headers, while the first whitespace token
is the unique sequence identity used in matrices. Alignments must be nonempty,
rectangular, and uniquely identified. `trim` counts both `-` and `.` as gaps,
never retains an all-gap column, and fails if no column survives. A column map
uses zero-based retained and original coordinates.

Distance matrices are labeled square TSV tables. The first header field is
`id`; row identities must exactly match the remaining header identities in the
same order. Values must be symmetric, diagonal zero, and finite and
nonnegative where required by the operation.

Hamming compares every byte literally. JC69 and K2P use pairwise deletion and
accept uppercase `A`, `C`, `G`, and `T` or `U`; every other byte is missing.
Saturated corrected distances are emitted as `NaN`. Neighbor joining emits an
unrooted three-way final join and clamps negative inferred branches by default;
`--negative-branches allow` retains them. UPGMA and WPGMA emit rooted
ultrametric trees. UPGMA weights merged-cluster distances by their number of
taxa; WPGMA weights the two clusters equally. WPGMA requires at least three
taxa and rejects the NJ-only negative-branch override.

Tree comparisons require the exact same named taxa and report any mismatch.
Newick underscores remain underscores unless `--convert-underscores` is set.
The root's optional length is not an edge and is excluded from branch metrics.
Weighted comparisons and branch-length paths fail when a required length is
missing or negative. Rooted split semantics are opt-in with `--rooted`.

The Hommola statistic follows scikit-bio's observed correlation. Permutations
use a deterministic SplitMix64 stream identified in the result; the stream is
not claimed to reproduce NumPy's random sequence. Version 0.1 runs this test
on one thread. `--threads` is accepted only by `distance`, so unsupported
parallelism cannot be requested silently.

Named outputs are written transactionally. Paired primary and secondary
outputs are committed together. Standard input is selected with `-`; only one
input per command may use it. Results go to standard output unless a path is
named. `--json` reserves standard output for the shared rsomics result or error
envelope and therefore requires a named data output.

Compatibility is checked against trimAl 1.5.1, scikit-bio 0.7.3, and DendroPy
5.0.8. The tracked oracle includes randomized distance and inference cases and
real CLI checks. Representative performance, fixture checksums, and raw timing
series are in `benchmarks/`.
