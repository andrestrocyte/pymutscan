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

## Exact radius-one and radius-two algorithms

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

At radius two, the complete Hamming ball has size:

$$
1 + 3L + 9\binom{L}{2}.
$$

For $L=30$, exactly 4,006 keys are possible. pymutscan encodes DNA with two
bits per base and caches the corresponding XOR masks once per sequence length
and radius. Candidate lookup therefore allocates no neighbor strings and has
expected $O(nL^2)$ time with $O(n+L^2)$ storage. Exhaustive-reference equality
is asserted in the public benchmark before performance is reported.

For fewer than 384 unique sequences, an adaptive exact exhaustive loop avoids
the fixed 4,006-probe cost. The threshold changes only the search strategy, not
ordering, score rules, or representative assignments.

## Information-preserving workflow

1. Stream reads to `raw_counts(barcode, sample_index, umi, read_count)` and,
   when a library identity is supplied, `library_counts`.
2. Optionally remove low-support combinations before clustering.
3. Compute `barcode_scores` by summing across sample index and UMI.
4. Derive `barcode_mapping(barcode, representative)` from barcode strings only.
5. Join the mapping to `raw_counts`.
6. Aggregate by representative, sample index, and UMI into `collapsed_counts`.

The raw count table and the mapping are retained. Filtering is therefore an
explicit, reversible parameter rather than an irreversible loss of provenance.

## UMI treatment

Barcode collapse and UMI deduplication answer different questions and never
share a distance calculation. UMI correction runs within each `(representative
barcode, sample index)` stratum and preserves both raw UMI and representative.
The compatibility method gives every UMI equal score. The optional directional
method uses the edge rule $n_p \ge 2n_c-1$ and traverses directed paths, so a
high-support root can absorb a low-support descendant through an intermediate
UMI only when every edge is abundance-consistent.

## Composition and distance boundaries

Read elements `S/U/C/V/P` follow mutscan terminology; `I` denotes the embedded
MAPseq sample index. Elements can occur on either read and a final `-1` consumes
the remainder. Hamming radii one and two are the high-throughput substitution
model. Levenshtein distance is available for indels, but uses an exhaustive
fallback because an indel-aware candidate index would change the performance
and memory guarantees.
