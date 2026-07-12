#!/usr/bin/env python3
"""Generate deterministic, public-safe MAPseq example FASTQs and truth tables."""

from __future__ import annotations

import csv
import gzip
import json
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "synthetic"
CONSTANT = "CCGTACT"
INDEX_A = "CGTGAT"
INDEX_B = "ACATCG"
INDEX_A_ERROR = "CGTGAC"

BARCODE_A = "ACGTACGTACGTACGTACGTACGTACGTAC"
BARCODE_A_ERROR = "TCGTACGTACGTACGTACGTACGTACGTAC"
BARCODE_B = "TGCATGCATGCATGCATGCATGCATGCATG"
BARCODE_C = "GATTACAGATTACAGATTACAGATTACAGA"

UMI_A = "AAAACCCCGGGGTTTT"
UMI_A_ERROR = "CAAACCCCGGGGTTTT"
UMI_B = "TTTTGGGGCCCCAAAA"
UMI_SWITCH = "ACACACACGTGTGTGT"


def fastq_record(name: str, sequence: str, quality: str) -> str:
    return f"@{name}\n{sequence}\n+\n{quality}\n"


def add_reads(
    records: list[tuple[str, str, str, str, str]],
    label: str,
    barcode: str,
    umi: str,
    index: str,
    count: int,
    *,
    constant: str = CONSTANT,
    quality_char: str = "I",
) -> None:
    for replicate in range(count):
        name = f"synthetic:{label}:{replicate:04d} 1:N:0:AAAAAA+CCCCCC"
        r1 = barcode + constant + "A"
        r2 = umi + index + "T"
        records.append((name, r1, quality_char * len(r1), r2, quality_char * len(r2)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str, str, str, str]] = []

    # Dominant biological molecules.
    add_reads(records, "a-main", BARCODE_A, UMI_A, INDEX_A, 40)
    add_reads(records, "a-second-umi", BARCODE_A, UMI_B, INDEX_A, 15)
    add_reads(records, "b-main", BARCODE_B, UMI_B, INDEX_B, 32)
    add_reads(records, "c-main", BARCODE_C, UMI_A, INDEX_B, 12)

    # One-base barcode, UMI, and sample-index errors.
    add_reads(records, "a-barcode-error", BARCODE_A_ERROR, UMI_A, INDEX_A, 4)
    add_reads(records, "a-umi-error", BARCODE_A, UMI_A_ERROR, INDEX_A, 3)
    add_reads(records, "a-index-error", BARCODE_A, UMI_A, INDEX_A_ERROR, 5)

    # A shared UMI across two distant barcodes mimics a template-switch candidate.
    add_reads(records, "switch-a", BARCODE_A, UMI_SWITCH, INDEX_A, 6)
    add_reads(records, "switch-b", BARCODE_B, UMI_SWITCH, INDEX_A, 5)

    # Reads expected to fail distinct QC filters.
    add_reads(records, "bad-constant", BARCODE_A, UMI_A, INDEX_A, 3, constant="AAAAAAA")
    add_reads(records, "low-quality", BARCODE_B, UMI_A, INDEX_B, 2, quality_char="!")
    add_reads(records, "ambiguous-barcode", "N" + BARCODE_A[1:], UMI_A, INDEX_A, 2)
    add_reads(records, "ambiguous-umi", BARCODE_C, "N" + UMI_B[1:], INDEX_B, 2)

    random.Random(20260712).shuffle(records)
    r1_path = OUT / "example_R1.fastq.gz"
    r2_path = OUT / "example_R2.fastq.gz"
    with gzip.open(r1_path, "wt", encoding="ascii") as r1_handle, gzip.open(
        r2_path, "wt", encoding="ascii"
    ) as r2_handle:
        for number, (_, r1, q1, r2, q2) in enumerate(records):
            pair_name = f"synthetic:{number:06d} 1:N:0:AAAAAA+CCCCCC"
            r1_handle.write(fastq_record(pair_name, r1, q1))
            r2_handle.write(fastq_record(pair_name, r2, q2))

    with (OUT / "sample_indices.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_name", "sample_index"])
        writer.writerow(["sample_A", INDEX_A])
        writer.writerow(["sample_B", INDEX_B])

    truth = {
        "description": "Synthetic MAPseq-like paired reads; no biological data.",
        "seed": 20260712,
        "total_read_pairs": len(records),
        "layout": {
            "R1": "barcode[30] + constant[7] + skipped_tail[1]",
            "R2": "UMI[16] + sample_index[6] + skipped_tail[1]",
        },
        "expected_filters": {
            "constant_mismatch": 3,
            "low_variable_quality": 2,
            "ambiguous_barcode": 2,
            "ambiguous_umi": 2,
            "retained_reads": len(records) - 9,
        },
        "canonical_barcodes": [BARCODE_A, BARCODE_B, BARCODE_C],
        "barcode_error": {BARCODE_A_ERROR: BARCODE_A},
        "sample_index_error": {INDEX_A_ERROR: INDEX_A},
        "umi_error": {UMI_A_ERROR: UMI_A},
        "template_switch_candidate": {
            "umi": UMI_SWITCH,
            "barcodes": [BARCODE_A, BARCODE_B],
            "sample_index": INDEX_A,
        },
    }
    (OUT / "truth.json").write_text(json.dumps(truth, indent=2) + "\n", encoding="utf-8")

    composition = Counter(item[0].split()[0].split(":")[1] for item in records)
    print(f"Wrote {len(records)} paired reads to {OUT}")
    print("Synthetic scenario labels:", dict(sorted(composition.items())))


if __name__ == "__main__":
    main()

