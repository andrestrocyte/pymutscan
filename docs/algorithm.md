# MAPseq data model and optimized collapse

## Experimental fields

The FASTQ header contains the Illumina UDI used to identify the sequenced
library. It is distinct from the reverse-transcription sample index embedded in
R2, which identifies one of the samples pooled before amplification.

Observed layouts in this repository are:

| Experiment | R1 | R2 |
| --- | --- | --- |
| 3994 | barcode (30) + virus constant (7) + skipped tail | UMI (16) + RT sample index (6) + skipped tail |
| 2502 | barcode (30) + virus constant (7) + skipped tail | UMI (18) + RT sample index (14) + skipped tail |

The four accepted 7-base virus constants are `CCGTACT`, `CTGTACT`, `TCGTACT`,
and `TTGTACT`. Spike-ins use a separate 24-base barcode and 13-base constant
layout and must be processed with a separate configuration.

## Why the former representation is expensive

With `elementsReverse = "VVS"`, mutscan concatenates both reverse variable
segments and then combines them with the forward variable sequence. For 3994,
the clustering key therefore represents 30 barcode bases, 16 UMI bases, and 6
sample-index bases. A million read pairs in the measured column-2 slice yielded
302,674 distinct composite keys but only 43,915 distinct barcodes.

This representation also changes the correction rule: two one-error barcodes
cannot merge unless their UMI and index strings are identical or sufficiently
close at the same time. That is not barcode-only error correction.

## Exact radius-one algorithm

For fixed length `L` over DNA, every string at Hamming distance at most one is
the sequence itself or one of `3L` single-base substitutions. Store all observed
barcodes in a hash table, order them by decreasing marginal read score with
lexicographic tie-breaking, and for every still-unassigned representative probe
those `1 + 3L` keys directly.

For the 30-base MAPseq barcode, each representative performs at most 91 hash
lookups. This gives expected `O(nL)` time and `O(n)` storage at radius one. It
exactly implements mutscan's greedy representative semantics, including
`collapseMinScore` and `collapseMinRatio`. The Python and R maps were identical
for every one of 43,915 barcodes in the real-data validation slice.

## Information-preserving workflow

1. Stream reads to `raw_counts(barcode, sample_index, umi, read_count)`.
2. Optionally remove low-support combinations before clustering.
3. Compute `barcode_scores` by summing across sample index and UMI.
4. Derive `barcode_mapping(barcode, representative)` from barcode strings only.
5. Join the mapping to `raw_counts`.
6. Aggregate by representative, sample index, and UMI into `collapsed_counts`.

The raw count table and the mapping are retained. Filtering is therefore an
explicit, reversible parameter rather than an irreversible loss of provenance.

## UMI treatment

UMIs are retained verbatim in the initial optimized workflow. Barcode collapse
and UMI deduplication answer different questions and should not share a distance
calculation. If UMI error correction is enabled later, it should run within each
`(representative barcode, sample index)` stratum and preserve both raw UMI and
UMI representative columns.

