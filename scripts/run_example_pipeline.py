#!/usr/bin/env python3
"""Run the complete pymutscan workflow on the bundled synthetic dataset."""

from __future__ import annotations

import json
from pathlib import Path

from pymutscan import collapse_database, collapse_umis, digest_fastqs, map_sample_indices
from pymutscan.pipeline import export_table

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "synthetic"
OUTPUT = ROOT / "examples" / "output"


def main() -> None:
    if not (DATA / "example_R1.fastq.gz").exists():
        raise SystemExit("Run `python scripts/generate_synthetic_mapseq.py` first.")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    database = OUTPUT / "example.sqlite"
    if database.exists():
        database.unlink()
    for suffix in ("-wal", "-shm"):
        sidecar = Path(str(database) + suffix)
        if sidecar.exists():
            sidecar.unlink()

    results = {
        "digest": digest_fastqs(
            DATA / "example_R1.fastq.gz",
            DATA / "example_R2.fastq.gz",
            database,
        ),
        "sample_indices": map_sample_indices(database, ["CGTGAT", "ACATCG"]),
        "barcodes": collapse_database(database, collapse_min_score=2),
        "umis": collapse_umis(database),
    }

    for table in (
        "raw_counts",
        "barcode_mapping",
        "sample_index_mapping",
        "collapsed_counts",
        "umi_mapping",
        "umi_collapsed_counts",
        "molecule_counts",
        "qc",
    ):
        export_table(database, table, OUTPUT / f"{table}.tsv.gz")

    (OUTPUT / "run_summary.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(results, indent=2, sort_keys=True))
    print(f"Outputs: {OUTPUT}")


if __name__ == "__main__":
    main()
