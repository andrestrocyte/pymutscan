"""Streaming MAPseq digestion, persistence, collapsing, and export."""

from __future__ import annotations

import csv
import gzip
import json
import sqlite3
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

from .collapse import group_similar_sequences, hamming_distance
from .fastq import paired_records


@dataclass(frozen=True)
class MapSeqConfig:
    barcode_length: int = 30
    umi_length: int = 16
    sample_index_length: int = 6
    constant_forward: tuple[str, ...] = ("CCGTACT", "CTGTACT", "TCGTACT", "TTGTACT")
    min_average_phred: float = 20.0
    max_ambiguous_barcode: int = 0
    max_ambiguous_umi: int = 0
    max_ambiguous_sample_index: int = 0

    @property
    def constant_length(self) -> int:
        lengths = {len(x) for x in self.constant_forward}
        if len(lengths) != 1:
            raise ValueError("all constant sequences must have equal length")
        return next(iter(lengths))


def _average_phred(quality: str) -> float:
    return sum(ord(char) - 33 for char in quality) / len(quality) if quality else 0.0


def _connect(path: str | Path) -> sqlite3.Connection:
    con = sqlite3.connect(path)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=MEMORY")
    return con


def _initialize(con: sqlite3.Connection) -> None:
    con.executescript(
        """
        CREATE TABLE IF NOT EXISTS raw_counts (
            barcode TEXT NOT NULL,
            sample_index TEXT NOT NULL,
            umi TEXT NOT NULL,
            read_count INTEGER NOT NULL,
            PRIMARY KEY (barcode, sample_index, umi)
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS qc (key TEXT PRIMARY KEY, value INTEGER NOT NULL);
        CREATE TABLE IF NOT EXISTS run_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        """
    )


def _flush_counts(con: sqlite3.Connection, counts: Counter[tuple[str, str, str]]) -> None:
    con.executemany(
        """
        INSERT INTO raw_counts(barcode, sample_index, umi, read_count)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(barcode, sample_index, umi)
        DO UPDATE SET read_count = read_count + excluded.read_count
        """,
        ((b, i, u, n) for (b, i, u), n in counts.items()),
    )
    con.commit()
    counts.clear()


def digest_fastqs(
    r1_path: str | Path,
    r2_path: str | Path,
    database: str | Path,
    *,
    config: MapSeqConfig = MapSeqConfig(),
    max_reads: int | None = None,
    flush_unique: int = 500_000,
) -> dict[str, int | float]:
    """Stream paired MAPseq reads into a normalized SQLite count table."""
    database = Path(database)
    database.parent.mkdir(parents=True, exist_ok=True)
    con = _connect(database)
    _initialize(con)
    counts: Counter[tuple[str, str, str]] = Counter()
    qc: Counter[str] = Counter()
    started = perf_counter()
    c0 = config.barcode_length
    c1 = c0 + config.constant_length
    u0 = config.umi_length
    u1 = u0 + config.sample_index_length
    constants = set(x.upper() for x in config.constant_forward)

    for _, seq1, qual1, seq2, qual2 in paired_records(r1_path, r2_path):
        if max_reads is not None and qc["total_reads"] >= max_reads:
            break
        qc["total_reads"] += 1
        if len(seq1) < c1 or len(seq2) < u1:
            qc["wrong_length"] += 1
            continue
        barcode = seq1[:c0]
        constant = seq1[c0:c1]
        umi = seq2[:u0]
        sample_index = seq2[u0:u1]
        if constant not in constants:
            qc["constant_mismatch"] += 1
            continue
        if barcode.count("N") > config.max_ambiguous_barcode:
            qc["ambiguous_barcode"] += 1
            continue
        if umi.count("N") > config.max_ambiguous_umi:
            qc["ambiguous_umi"] += 1
            continue
        if sample_index.count("N") > config.max_ambiguous_sample_index:
            qc["ambiguous_sample_index"] += 1
            continue
        if (
            _average_phred(qual1[:c0]) < config.min_average_phred
            or _average_phred(qual2[:u1]) < config.min_average_phred
        ):
            qc["low_variable_quality"] += 1
            continue
        counts[(barcode, sample_index, umi)] += 1
        qc["retained_reads"] += 1
        if len(counts) >= flush_unique:
            _flush_counts(con, counts)
    if counts:
        _flush_counts(con, counts)

    con.executemany(
        "INSERT OR REPLACE INTO qc(key, value) VALUES (?, ?)", sorted(qc.items())
    )
    metadata = {
        "r1_path": str(Path(r1_path).resolve()),
        "r2_path": str(Path(r2_path).resolve()),
        "config": json.dumps(asdict(config), sort_keys=True),
        "format_version": "1",
    }
    con.executemany(
        "INSERT OR REPLACE INTO run_metadata(key, value) VALUES (?, ?)", metadata.items()
    )
    con.commit()
    unique_combinations = con.execute("SELECT count(*) FROM raw_counts").fetchone()[0]
    con.close()
    return {
        **qc,
        "unique_combinations": unique_combinations,
        "elapsed_seconds": perf_counter() - started,
    }


def collapse_database(
    database: str | Path,
    *,
    collapse_max_dist: float = 1,
    collapse_min_score: float = 2,
    collapse_min_ratio: float = 0,
    min_combo_reads: int = 1,
) -> dict[str, int | float]:
    """Collapse barcodes only and retain sample-index/UMI stratification."""
    con = _connect(database)
    started = perf_counter()
    has_index_mapping = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sample_index_mapping'"
    ).fetchone() is not None
    con.executescript(
        """
        DROP TABLE IF EXISTS barcode_scores;
        DROP TABLE IF EXISTS barcode_mapping;
        DROP TABLE IF EXISTS collapsed_counts;
        CREATE TABLE barcode_scores (
            barcode TEXT PRIMARY KEY,
            score INTEGER NOT NULL
        ) WITHOUT ROWID;
        """
    )
    con.execute(
        """
        INSERT INTO barcode_scores(barcode, score)
        SELECT barcode, sum(read_count)
        FROM raw_counts
        WHERE read_count >= ?
        GROUP BY barcode
        """,
        (min_combo_reads,),
    )
    rows = con.execute("SELECT barcode, score FROM barcode_scores").fetchall()
    sequences = [row[0] for row in rows]
    scores = {row[0]: row[1] for row in rows}
    mapping = group_similar_sequences(
        sequences,
        scores,
        collapse_max_dist=collapse_max_dist,
        collapse_min_score=collapse_min_score,
        collapse_min_ratio=collapse_min_ratio,
    )
    con.execute(
        """
        CREATE TABLE barcode_mapping (
            barcode TEXT PRIMARY KEY,
            representative TEXT NOT NULL
        ) WITHOUT ROWID
        """
    )
    con.executemany(
        "INSERT INTO barcode_mapping(barcode, representative) VALUES (?, ?)",
        mapping.items(),
    )
    index_join = "LEFT JOIN sample_index_mapping i ON i.observed_index = r.sample_index" if has_index_mapping else ""
    index_value = "COALESCE(i.sample_index, r.sample_index)" if has_index_mapping else "r.sample_index"
    con.executescript(
        f"""
        CREATE TABLE collapsed_counts (
            barcode TEXT NOT NULL,
            sample_index TEXT NOT NULL,
            umi TEXT NOT NULL,
            read_count INTEGER NOT NULL,
            PRIMARY KEY (barcode, sample_index, umi)
        ) WITHOUT ROWID;
        INSERT INTO collapsed_counts(barcode, sample_index, umi, read_count)
        SELECT m.representative, {index_value}, r.umi, sum(r.read_count)
        FROM raw_counts r
        JOIN barcode_mapping m USING (barcode)
        {index_join}
        WHERE r.read_count >= {int(min_combo_reads)}
        GROUP BY m.representative, {index_value}, r.umi;
        """
    )
    parameters = {
        "collapse_max_dist": collapse_max_dist,
        "collapse_min_score": collapse_min_score,
        "collapse_min_ratio": collapse_min_ratio,
        "min_combo_reads": min_combo_reads,
    }
    con.execute(
        "INSERT OR REPLACE INTO run_metadata(key, value) VALUES ('collapse_parameters', ?)",
        (json.dumps(parameters, sort_keys=True),),
    )
    con.commit()
    n_representatives = con.execute(
        "SELECT count(DISTINCT representative) FROM barcode_mapping"
    ).fetchone()[0]
    n_collapsed_combinations = con.execute(
        "SELECT count(*) FROM collapsed_counts"
    ).fetchone()[0]
    con.close()
    return {
        "input_barcodes": len(sequences),
        "representatives": n_representatives,
        "collapsed_combinations": n_collapsed_combinations,
        "elapsed_seconds": perf_counter() - started,
    }


def map_sample_indices(
    database: str | Path,
    sample_indices: list[str] | tuple[str, ...],
    *,
    max_distance: int = 1,
) -> dict[str, int]:
    """Map observed RT sample indices to a known whitelist.

    A non-exact string is corrected only when it has one unique closest
    whitelist entry within ``max_distance``. Raw strings remain in
    ``raw_counts`` and every decision is stored in ``sample_index_mapping``.
    """
    whitelist = tuple(dict.fromkeys(x.upper() for x in sample_indices))
    if not whitelist:
        raise ValueError("sample_indices must not be empty")
    lengths = {len(x) for x in whitelist}
    if len(lengths) != 1:
        raise ValueError("all sample indices must have equal length")
    con = _connect(database)
    observed = [row[0] for row in con.execute("SELECT DISTINCT sample_index FROM raw_counts")]
    rows: list[tuple[str, str | None, int | None, str]] = []
    status_counts: Counter[str] = Counter()
    for index in observed:
        if len(index) != next(iter(lengths)):
            row = (index, None, None, "unassigned")
        else:
            distances = [(hamming_distance(index, expected), expected) for expected in whitelist]
            best_distance = min(distance for distance, _ in distances)
            best = [expected for distance, expected in distances if distance == best_distance]
            if best_distance > max_distance:
                row = (index, None, best_distance, "unassigned")
            elif len(best) > 1:
                row = (index, None, best_distance, "ambiguous")
            elif best_distance == 0:
                row = (index, best[0], 0, "exact")
            else:
                row = (index, best[0], best_distance, "corrected")
        rows.append(row)
        status_counts[row[3]] += 1
    con.executescript(
        """
        DROP TABLE IF EXISTS sample_index_mapping;
        CREATE TABLE sample_index_mapping (
            observed_index TEXT PRIMARY KEY,
            sample_index TEXT,
            distance INTEGER,
            status TEXT NOT NULL
        ) WITHOUT ROWID;
        """
    )
    con.executemany(
        "INSERT INTO sample_index_mapping(observed_index, sample_index, distance, status) VALUES (?, ?, ?, ?)",
        rows,
    )
    con.execute(
        "INSERT OR REPLACE INTO run_metadata(key, value) VALUES ('sample_index_parameters', ?)",
        (json.dumps({"whitelist": whitelist, "max_distance": max_distance}),),
    )
    con.commit()
    read_counts = dict(
        con.execute(
            """
            SELECT i.status, sum(r.read_count)
            FROM raw_counts r JOIN sample_index_mapping i
            ON i.observed_index = r.sample_index
            GROUP BY i.status
            """
        ).fetchall()
    )
    con.close()
    return {
        **{f"indices_{key}": value for key, value in status_counts.items()},
        **{f"reads_{key}": value for key, value in read_counts.items()},
    }


def collapse_umis(
    database: str | Path,
    *,
    collapse_max_dist: float = 1,
) -> dict[str, int | float]:
    """Collapse UMIs within each barcode/sample-index stratum.

    All UMI scores are equal, reproducing mutscan's lexicographic greedy UMI
    representative choice. The input ``collapsed_counts`` and a complete
    ``umi_mapping`` remain available for audit.
    """
    con = _connect(database)
    if con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='collapsed_counts'"
    ).fetchone() is None:
        con.close()
        raise ValueError("collapsed_counts is required; run barcode collapse first")
    started = perf_counter()
    con.executescript(
        """
        DROP TABLE IF EXISTS umi_mapping;
        DROP TABLE IF EXISTS umi_collapsed_counts;
        DROP TABLE IF EXISTS molecule_counts;
        CREATE TABLE umi_mapping (
            barcode TEXT NOT NULL,
            sample_index TEXT NOT NULL,
            umi TEXT NOT NULL,
            representative TEXT NOT NULL,
            PRIMARY KEY (barcode, sample_index, umi)
        ) WITHOUT ROWID;
        CREATE TABLE umi_collapsed_counts (
            barcode TEXT NOT NULL,
            sample_index TEXT NOT NULL,
            umi TEXT NOT NULL,
            read_count INTEGER NOT NULL,
            PRIMARY KEY (barcode, sample_index, umi)
        ) WITHOUT ROWID;
        """
    )
    cursor = con.execute(
        "SELECT barcode, sample_index, umi, read_count FROM collapsed_counts ORDER BY barcode, sample_index"
    )
    current: tuple[str, str] | None = None
    group: list[tuple[str, int]] = []
    mapping_rows: list[tuple[str, str, str, str]] = []
    collapsed_rows: list[tuple[str, str, str, int]] = []
    n_strata = n_input = n_representatives = 0

    def flush_group(key: tuple[str, str] | None, values: list[tuple[str, int]]) -> None:
        nonlocal n_strata, n_input, n_representatives
        if key is None or not values:
            return
        barcode, sample_index = key
        umis = [item[0] for item in values]
        mapping = group_similar_sequences(
            umis, [1.0] * len(umis), collapse_max_dist=collapse_max_dist
        )
        sums: Counter[str] = Counter()
        for umi, reads in values:
            representative = mapping[umi]
            mapping_rows.append((barcode, sample_index, umi, representative))
            sums[representative] += reads
        collapsed_rows.extend(
            (barcode, sample_index, representative, reads)
            for representative, reads in sums.items()
        )
        n_strata += 1
        n_input += len(values)
        n_representatives += len(sums)

    for barcode, sample_index, umi, read_count in cursor:
        key = (barcode, sample_index)
        if current is not None and key != current:
            flush_group(current, group)
            group = []
            if len(mapping_rows) >= 100_000:
                con.executemany("INSERT INTO umi_mapping VALUES (?, ?, ?, ?)", mapping_rows)
                con.executemany("INSERT INTO umi_collapsed_counts VALUES (?, ?, ?, ?)", collapsed_rows)
                con.commit()
                mapping_rows.clear()
                collapsed_rows.clear()
        current = key
        group.append((umi, int(read_count)))
    flush_group(current, group)
    if mapping_rows:
        con.executemany("INSERT INTO umi_mapping VALUES (?, ?, ?, ?)", mapping_rows)
        con.executemany("INSERT INTO umi_collapsed_counts VALUES (?, ?, ?, ?)", collapsed_rows)
    con.executescript(
        """
        CREATE TABLE molecule_counts (
            barcode TEXT NOT NULL,
            sample_index TEXT NOT NULL,
            molecule_count INTEGER NOT NULL,
            read_count INTEGER NOT NULL,
            PRIMARY KEY (barcode, sample_index)
        ) WITHOUT ROWID;
        INSERT INTO molecule_counts(barcode, sample_index, molecule_count, read_count)
        SELECT barcode, sample_index, count(*), sum(read_count)
        FROM umi_collapsed_counts
        GROUP BY barcode, sample_index;
        """
    )
    con.execute(
        "INSERT OR REPLACE INTO run_metadata(key, value) VALUES ('umi_collapse_parameters', ?)",
        (json.dumps({"collapse_max_dist": collapse_max_dist, "score_method": "equal_lexicographic"}),),
    )
    con.commit()
    con.close()
    return {
        "strata": n_strata,
        "input_umis": n_input,
        "representative_umis": n_representatives,
        "elapsed_seconds": perf_counter() - started,
    }


def export_table(database: str | Path, table: str, output: str | Path) -> int:
    """Export an allowed database table as TSV or TSV.GZ."""
    allowed = {"raw_counts", "barcode_scores", "barcode_mapping", "sample_index_mapping", "collapsed_counts", "umi_mapping", "umi_collapsed_counts", "molecule_counts", "qc", "run_metadata"}
    if table not in allowed:
        raise ValueError(f"table must be one of {sorted(allowed)}")
    con = sqlite3.connect(database)
    cursor = con.execute(f"SELECT * FROM {table}")
    columns = [item[0] for item in cursor.description]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if output.suffix == ".gz" else open
    count = 0
    with opener(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in cursor:
            writer.writerow(row)
            count += 1
    con.close()
    return count
