import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from pymutscan import collapse_database, collapse_umis, digest_fastqs, map_sample_indices

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "synthetic"


class SyntheticExampleTests(unittest.TestCase):
    def test_example_matches_documented_truth(self):
        truth = json.loads((DATA / "truth.json").read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as tmp:
            database = Path(tmp) / "example.sqlite"
            qc = digest_fastqs(
                DATA / "example_R1.fastq.gz",
                DATA / "example_R2.fastq.gz",
                database,
            )
            for key, expected in truth["expected_filters"].items():
                self.assertEqual(qc[key], expected)

            index_result = map_sample_indices(database, ["CGTGAT", "ACATCG"])
            self.assertEqual(index_result["reads_corrected"], 5)
            barcode_result = collapse_database(database, collapse_min_score=2)
            self.assertEqual(barcode_result["input_barcodes"], 4)
            self.assertEqual(barcode_result["representatives"], 3)
            umi_result = collapse_umis(database)
            self.assertLess(umi_result["representative_umis"], umi_result["input_umis"])

            con = sqlite3.connect(database)
            raw_total = con.execute("SELECT sum(read_count) FROM raw_counts").fetchone()[0]
            collapsed_total = con.execute(
                "SELECT sum(read_count) FROM collapsed_counts"
            ).fetchone()[0]
            umi_total = con.execute(
                "SELECT sum(read_count) FROM umi_collapsed_counts"
            ).fetchone()[0]
            barcode_error = next(iter(truth["barcode_error"]))
            expected_rep = truth["barcode_error"][barcode_error]
            actual_rep = con.execute(
                "SELECT representative FROM barcode_mapping WHERE barcode = ?",
                (barcode_error,),
            ).fetchone()[0]
            switch_umi = truth["template_switch_candidate"]["umi"]
            switch_barcodes = con.execute(
                """
                SELECT count(DISTINCT barcode)
                FROM collapsed_counts
                WHERE umi = ? AND sample_index = ?
                """,
                (switch_umi, truth["template_switch_candidate"]["sample_index"]),
            ).fetchone()[0]
            con.close()

            self.assertEqual(raw_total, truth["expected_filters"]["retained_reads"])
            self.assertEqual(raw_total, collapsed_total)
            self.assertEqual(raw_total, umi_total)
            self.assertEqual(actual_rep, expected_rep)
            self.assertEqual(switch_barcodes, 2)


if __name__ == "__main__":
    unittest.main()
