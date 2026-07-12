# Benchmark methodology

## Question

Can barcode-only marginalization plus exact radius-one hash lookup reduce the
MAPseq collapse step by at least one order of magnitude while preserving the R
mutscan barcode representative map?

## Compared workflows

1. R `mutscan::groupSimilarSequences()` on complete barcode/UMI/index feature
   keys from a fixed one-million-pair read slice.
2. R `groupSimilarSequences()` on barcode-only marginal scores.
3. Python barcode-only grouping on identical barcodes and scores.
4. Python grouping plus SQLite mapping persistence and re-aggregation.

Parameters were Hamming radius 1, minimum representative score 2, minimum score
ratio 0, and no singleton-combination filter.

## Correctness checks

- Compare sequence-to-representative maps, not only representative counts.
- Reconcile retained-read totals across raw, barcode-collapsed,
  UMI-collapsed, and molecule-summary tables.
- Compare the full retained total with the existing RDS.
- Validate SQLite integrity and mapping coverage.

## Results

| Measurement | Result |
|---|---:|
| Input read pairs in timing slice | 1,000,000 |
| Retained read pairs | 943,816 |
| Composite keys | 302,674 |
| Unique barcodes | 43,915 |
| Composite R grouping | 442.7 s |
| Barcode-only R grouping | 3.856 s |
| Python grouping | 0.717 s |
| Python full collapse/re-aggregation | 1.725 s |
| Python/R barcode-map differences | 0 |

## Caveat

The early part of the composite R timing experienced CPU contention from two
unrelated RDS inspection processes. Treat 256.6× as an observed workflow timing,
not a hardware-pure language microbenchmark. The cardinality reduction,
uncontended barcode-only timings, and exact map equality independently support
the conclusion that the improvement safely exceeds 10×.

