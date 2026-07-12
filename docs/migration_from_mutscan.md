# Migrating a MAPseq workflow from mutscan

## Conceptual mapping

| R mutscan concept | pymutscan equivalent |
|---|---|
| `digestFastqs()` for contiguous MAPseq fields | `pymutscan digest` / `digest_fastqs()` |
| `summaryTable` | SQLite `raw_counts` |
| `groupSimilarSequences()` on barcodes | `barcode_mapping` via `pymutscan collapse` |
| sample index parsed from a composite name | explicit `sample_index` plus `sample_index_mapping` |
| unique/collapsed UMI count | `umi_mapping`, `umi_collapsed_counts`, `molecule_counts` |
| `SummarizedExperiment` RDS | SQLite/TSV primary output plus optional R exporter |

## Recommended migration

1. Record the exact R1/R2 layout and the RT-index whitelist.
2. Run pymutscan on a bounded read subset.
3. Reconcile retained-read totals with `digestFastqs()`.
4. Aggregate scores by barcode and compare the complete Python representative
   map with `mutscan::groupSimilarSequences()` under identical thresholds.
5. Inspect exact, corrected, ambiguous, and unassigned RT-index read totals.
6. Compare read-count and molecule-count matrices by sample.
7. Only then scale to the full dataset.

## Identifier change

Legacy notebooks often parse a row name such as:

```text
BARCODE_UMIINDEX
```

The compatibility RDS exporter instead provides explicit row metadata:

```text
barcode | sampleIndex | umi | sequence
```

Prefer those fields over parsing row names. A composite `sequence` value remains
available for software that requires a unique row identifier.

## Threshold equivalence

For exact representative-map comparison, keep these aligned:

- barcode score definition (normally total reads across samples/UMIs);
- `collapseMaxDist` / `collapse_max_dist`;
- `collapseMinScore` / `collapse_min_score`;
- `collapseMinRatio` / `collapse_min_ratio`;
- singleton-combination filtering before barcode marginalization.

Changing the clustering grain intentionally changes the biological correction:
barcode-only grouping is not expected to reproduce a composite
barcode/UMI/index map.

