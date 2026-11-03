# WPGMA resource evidence

Measured November3,2026 on Apple M2,8GiB, macOS27.2 arm64, Rust1.91.0.
Source implementation784cef7cd42640880a1d24f35954edea54008fdb used a local
release binary SHA256
4a53cffab9503280fbcab521137a8088ed3c849050b43e97f9cf1b79f971f6e3.
This is not a measurement of every subsequently packaged native binary.
Earlier operation measurements remain separately identified in `benchmarks/`.

Both implementations read the same119918537B labeled2500-taxon distance
matrix, generated from eight-dimensional LCG(seed77) points. SHA256
2bf32f3e62c6e44216e5aad79202d1aa12dce8e3533a36a2a31383c9537beb8b.
The reference uses scikit-bio0.7.3 `upgma(weighted=True)`, SciPy1.18.0,
NumPy2.5.1 and Python3.12.11. BLAS/OpenMP controls are1. There is no native
parallel-inference claim and the shared host is not exclusive.

One retained warmup per implementation preceded five alternating measured
pairs. All twelve complete trees agree in rooted topology, tip identities
and branch lengths; an independent clade-map replay found maximum absolute
branch difference0.0. No sampled-output checksum substitutes for this check.

| Measurement | Native | Reference |
|---|---:|---:|
| Median maximum-process RSS |224673792B|275070976B|
| Median complete-process elapsed |0.290964458s|2.372631208s|

The accepted resource advantage is50397184B/18.32%lower median RSS. RSS comes
from macOS `/usr/bin/time -l`, not an aggregate process-tree measurement.
Complete-process elapsed includes startup, input, inference and output;
Python import is included. Its8.15x ratio is not an algorithm-throughput claim.
Reference inference-only times have a different scope from native CLI time.
No cold-cache, affinity, cross-machine or universal input-size claim is made.

Raw measured vectors in execution order:

```text
native wall: 0.308947833 0.292775250 0.287726709 0.289310833 0.290964458
source wall: 2.376230208 2.372631208 2.366805250 2.401765500 2.369910375
native RSS: 224673792 224673792 224706560 224657408 224641024
source RSS: 275070976 271941632 280182784 279576576 269762560
source core: 0.070483291 0.067601750 0.073021250 0.075880750 0.071606250
warmup native: 0.339880583s 224657408B
warmup source: 2.578822500s 271302656B
```

Original commands, outputs, timing records, provenance and input remain in
the external evidence packet
`phylo-wpgma-performance-20261103.sEuLzo9O`. Canonical control-plane receipt:
`rsomics-world/.autopilot/oracles/phylo-wpgma-performance-2026-11-03.json`.
