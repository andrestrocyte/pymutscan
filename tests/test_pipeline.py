import gzip
import sqlite3
import tempfile
import unittest
from pathlib import Path

from pymutscan.pipeline import (
    MapSeqConfig,
    collapse_database,
    collapse_umis,
    digest_fastqs,
    map_sample_indices,
)


def write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'F' * len(sequence)}\n")


class PipelineTests(unittest.TestCase):
    def test_digest_and_barcode_only_collapse(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1, r2, db = tmp / "r1.fq.gz", tmp / "r2.fq.gz", tmp / "counts.sqlite"
            barcode = "A" * 30
            neighbor = "C" + "A" * 29
            constant = "CCGTACT"
            reads1 = [("r1 1:N:0:X", barcode + constant), ("r2 1:N:0:X", neighbor + constant)]
            reads2 = [("r1 2:N:0:X", "G" * 16 + "AACCGG"), ("r2 2:N:0:X", "T" * 16 + "AACCGG")]
            write_fastq(r1, reads1)
            write_fastq(r2, reads2)
            qc = digest_fastqs(r1, r2, db, config=MapSeqConfig())
            self.assertEqual(qc["retained_reads"], 2)
            result = collapse_database(db, collapse_min_score=0)
            self.assertEqual(result["input_barcodes"], 2)
            self.assertEqual(result["representatives"], 1)
            con = sqlite3.connect(db)
            rows = con.execute(
                "SELECT barcode, sample_index, umi, read_count FROM collapsed_counts ORDER BY umi"
            ).fetchall()
            con.close()
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(row[0] == barcode for row in rows))
            self.assertEqual({row[1] for row in rows}, {"AACCGG"})
            self.assertEqual({row[2] for row in rows}, {"G" * 16, "T" * 16})

    def test_sample_index_mapping_is_independent_and_auditable(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "counts.sqlite"
            con = sqlite3.connect(db)
            con.execute("CREATE TABLE raw_counts(barcode, sample_index, umi, read_count)")
            con.executemany(
                "INSERT INTO raw_counts VALUES (?, ?, ?, ?)",
                [("A" * 30, "AACCGG", "G" * 16, 4), ("A" * 30, "AACCGT", "T" * 16, 2)],
            )
            con.execute("CREATE TABLE run_metadata(key TEXT PRIMARY KEY, value TEXT)")
            con.commit()
            con.close()
            result = map_sample_indices(db, ["AACCGG", "TTTTTT"])
            self.assertEqual(result["reads_exact"], 4)
            self.assertEqual(result["reads_corrected"], 2)
            con = sqlite3.connect(db)
            row = con.execute(
                "SELECT sample_index, distance, status FROM sample_index_mapping WHERE observed_index='AACCGT'"
            ).fetchone()
            con.close()
            self.assertEqual(row, ("AACCGG", 1, "corrected"))

    def test_umi_collapse_stays_within_barcode_and_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "counts.sqlite"
            con = sqlite3.connect(db)
            con.execute("CREATE TABLE collapsed_counts(barcode, sample_index, umi, read_count)")
            con.executemany(
                "INSERT INTO collapsed_counts VALUES (?, ?, ?, ?)",
                [
                    ("A" * 30, "INDEX1", "AAAAAAAAAAAAAAAA", 3),
                    ("A" * 30, "INDEX1", "CAAAAAAAAAAAAAAA", 2),
                    ("A" * 30, "INDEX2", "CAAAAAAAAAAAAAAA", 5),
                ],
            )
            con.execute("CREATE TABLE run_metadata(key TEXT PRIMARY KEY, value TEXT)")
            con.commit()
            con.close()
            result = collapse_umis(db)
            self.assertEqual(result["input_umis"], 3)
            self.assertEqual(result["representative_umis"], 2)
            con = sqlite3.connect(db)
            rows = con.execute(
                "SELECT sample_index, umi, read_count FROM umi_collapsed_counts ORDER BY sample_index"
            ).fetchall()
            molecules = con.execute(
                "SELECT sample_index, molecule_count, read_count FROM molecule_counts ORDER BY sample_index"
            ).fetchall()
            con.close()
            self.assertEqual(rows, [("INDEX1", "AAAAAAAAAAAAAAAA", 5), ("INDEX2", "CAAAAAAAAAAAAAAA", 5)])
            self.assertEqual(molecules, [("INDEX1", 1, 5), ("INDEX2", 1, 5)])


if __name__ == "__main__":
    unittest.main()
