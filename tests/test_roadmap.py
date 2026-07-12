import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from pymutscan import (
    MapSeqConfig,
    call_template_switch_evidence,
    collapse_database,
    collapse_umis,
    digest_fastqs,
    digest_libraries,
    edit_distance,
    export_sparse_matrix,
    group_directional_sequences,
    group_similar_sequences,
    import_sample_metadata,
    load_config,
    merge_databases,
)

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / "examples/synthetic/example_R1.fastq.gz"
R2 = ROOT / "examples/synthetic/example_R2.fastq.gz"


class RoadmapTests(unittest.TestCase):
    def test_exact_radius_two_and_levenshtein(self):
        sequences = ["AAAA", "AACC", "TTTT", "AAA"]
        mapping = group_similar_sequences(
            sequences[:3], [100, 2, 10], collapse_max_dist=2, collapse_min_score=0
        )
        self.assertEqual(mapping["AACC"], "AAAA")
        self.assertEqual(mapping["TTTT"], "TTTT")
        indel_mapping = group_similar_sequences(
            ["AAAA", "AAA"],
            [10, 1],
            collapse_max_dist=1,
            distance_metric="levenshtein",
        )
        self.assertEqual(indel_mapping["AAA"], "AAAA")
        self.assertEqual(edit_distance("ACGT", "AGT"), 1)

    def test_directional_umi_rule(self):
        directional = group_directional_sequences(
            ["AAAA", "AAAT", "AATT"],
            [10, 6, 2],
            collapse_max_dist=1,
        )
        self.assertEqual(directional["AAAT"], "AAAT")
        self.assertEqual(directional["AATT"], "AAAT")
        transitive = group_directional_sequences(
            ["AAAA", "AAAT", "AATT"], [10, 5, 3], collapse_max_dist=1
        )
        self.assertEqual(transitive["AATT"], "AAAA")
        radius_two = group_directional_sequences(
            ["AAAA", "AACC", "TTTT"], [10, 2, 5], collapse_max_dist=2
        )
        self.assertEqual(radius_two["AACC"], "AAAA")

    def test_native_json_and_toml_configuration(self):
        legacy_positional = MapSeqConfig(30, 16, 6, ("CCGTACT",), 25.0, 1, 2, 3)
        self.assertEqual(legacy_positional.min_average_phred, 25.0)
        self.assertEqual(legacy_positional.max_ambiguous_sample_index, 3)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "config.json").write_text(
                json.dumps({"config": {"barcode_length": 24, "umi_length": 12}}),
                encoding="utf-8",
            )
            (tmp / "config.toml").write_text(
                "[presets.custom]\nbarcode_length = 27\numi_length = 10\n",
                encoding="utf-8",
            )
            self.assertEqual(load_config(tmp / "config.json").barcode_length, 24)
            self.assertEqual(load_config(tmp / "config.toml", preset="custom").barcode_length, 27)

    def test_single_end_composition(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fastq = tmp / "single.fastq"
            sequence = "G" + "TTTT" + "AAA" + "AACCGG" + "CC" + "ACGT"
            fastq.write_text(f"@read\n{sequence}\n+\n{'I' * len(sequence)}\n", encoding="ascii")
            config = MapSeqConfig(
                constant_forward=("AAA",),
                elements_forward="SUCVPI",
                element_lengths_forward=(1, 4, 3, 6, 2, 4),
                primer_forward=("CC",),
            )
            database = tmp / "single.sqlite"
            result = digest_fastqs(fastq, None, database, config=config)
            self.assertEqual(result["retained_reads"], 1)
            row = sqlite3.connect(database).execute("SELECT * FROM raw_counts").fetchone()
            self.assertEqual(row, ("AACCGG", "ACGT", "TTTT", 1))

            desired = "AACCGG" + "TTTT" + "ACGT"
            reverse_complement = desired.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            reverse_fastq = tmp / "reverse.fastq"
            reverse_fastq.write_text(
                f"@read\n{reverse_complement}\n+\n{'I' * len(reverse_complement)}\n",
                encoding="ascii",
            )
            reverse_database = tmp / "reverse.sqlite"
            digest_fastqs(
                reverse_fastq,
                None,
                reverse_database,
                config=MapSeqConfig(
                    constant_forward=(),
                    elements_forward="VUI",
                    element_lengths_forward=(6, 4, 4),
                    reverse_complement_forward=True,
                ),
            )
            reverse_row = sqlite3.connect(reverse_database).execute(
                "SELECT * FROM raw_counts"
            ).fetchone()
            self.assertEqual(reverse_row, ("AACCGG", "ACGT", "TTTT", 1))

    def test_multilibrary_merge_metadata_sparse_and_switch_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            combined = tmp / "combined.sqlite"
            result = digest_libraries(
                [
                    {"library_id": "lane-1", "lane": "L001", "r1": str(R1), "r2": str(R2)},
                    {"library_id": "lane-2", "lane": "L002", "r1": str(R1), "r2": str(R2)},
                ],
                combined,
            )
            self.assertEqual(result["libraries"], 2)
            con = sqlite3.connect(combined)
            self.assertEqual(con.execute("SELECT count(*) FROM libraries").fetchone()[0], 2)
            self.assertEqual(
                con.execute("SELECT sum(read_count) FROM library_counts").fetchone()[0], 244
            )
            con.close()

            first = tmp / "first.sqlite"
            second = tmp / "second.sqlite"
            digest_fastqs(R1, R2, first, library_id="first")
            digest_fastqs(R1, R2, second, library_id="second")
            merged = tmp / "merged.sqlite"
            self.assertEqual(merge_databases([first, second], merged)["libraries"], 2)
            self.assertEqual(
                sqlite3.connect(merged).execute("SELECT sum(read_count) FROM raw_counts").fetchone()[0],
                244,
            )
            legacy = tmp / "legacy.sqlite"
            legacy_con = sqlite3.connect(legacy)
            legacy_con.execute(
                "CREATE TABLE raw_counts(barcode TEXT, sample_index TEXT, umi TEXT, read_count INTEGER, PRIMARY KEY(barcode, sample_index, umi)) WITHOUT ROWID"
            )
            legacy_con.execute("INSERT INTO raw_counts VALUES ('AAAA', 'CCCC', 'GGGG', 3)")
            legacy_con.commit()
            legacy_con.close()
            legacy_merged = tmp / "legacy-merged.sqlite"
            self.assertEqual(merge_databases([legacy], legacy_merged)["libraries"], 1)

            collapse_database(combined)
            collapse_umis(combined, method="directional")
            con = sqlite3.connect(combined)
            self.assertEqual(
                con.execute("SELECT sum(read_count) FROM library_umi_collapsed_counts").fetchone()[0],
                244,
            )
            self.assertEqual(
                con.execute("SELECT count(DISTINCT library_id) FROM library_molecule_counts").fetchone()[0],
                2,
            )
            con.close()
            metadata = tmp / "metadata.tsv"
            metadata.write_text(
                "sample_index\tregion\tanimal\nCGTGAT\tbrain\tA1\nACATCG\tbrain\tA1\n",
                encoding="utf-8",
            )
            self.assertEqual(import_sample_metadata(combined, metadata), 2)
            sparse = export_sparse_matrix(combined, tmp / "matrix")
            self.assertGreater(sparse["nonzero"], 0)
            self.assertTrue((tmp / "matrix/matrix.mtx").read_text().startswith("%%MatrixMarket"))
            evidence = call_template_switch_evidence(combined, min_reads_per_barcode=2)
            self.assertGreaterEqual(evidence["candidates"], 1)


if __name__ == "__main__":
    unittest.main()
