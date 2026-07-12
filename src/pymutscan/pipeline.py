"""Streaming MAPseq digestion, persistence, collapsing, and export."""

from __future__ import annotations

import csv
import gzip
import json
import sqlite3
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from time import perf_counter

from .collapse import group_directional_sequences, group_similar_sequences, hamming_distance
from .config import MapSeqConfig
from .fastq import read_sets


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
        CREATE TABLE IF NOT EXISTS libraries (
            library_id TEXT PRIMARY KEY,
            r1_path TEXT NOT NULL,
            r2_path TEXT,
            lane TEXT,
            metadata_json TEXT NOT NULL
        ) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS library_counts (
            library_id TEXT NOT NULL,
            barcode TEXT NOT NULL,
            sample_index TEXT NOT NULL,
            umi TEXT NOT NULL,
            read_count INTEGER NOT NULL,
            PRIMARY KEY (library_id, barcode, sample_index, umi)
        ) WITHOUT ROWID;
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


def _flush_library_counts(
    con: sqlite3.Connection,
    library_id: str,
    counts: Counter[tuple[str, str, str]],
) -> None:
    con.executemany(
        """
        INSERT INTO library_counts(library_id, barcode, sample_index, umi, read_count)
        VALUES (?, ?, ?, ?, ?)
        ON CONFLICT(library_id, barcode, sample_index, umi)
        DO UPDATE SET read_count = read_count + excluded.read_count
        """,
        ((library_id, barcode, index, umi, count) for (barcode, index, umi), count in counts.items()),
    )


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _extract_composition(
    sequence: str,
    quality: str,
    elements: str,
    lengths: tuple[int, ...],
) -> tuple[dict[str, str], dict[str, str]] | None:
    offset = 0
    sequences: dict[str, list[str]] = {key: [] for key in "VUCPIS"}
    qualities: dict[str, list[str]] = {key: [] for key in "VUCPIS"}
    for element, length in zip(elements.upper(), lengths):
        end = len(sequence) if length == -1 else offset + length
        if end > len(sequence):
            return None
        sequences[element].append(sequence[offset:end])
        qualities[element].append(quality[offset:end])
        offset = end
    return (
        {key: "".join(parts) for key, parts in sequences.items()},
        {key: "".join(parts) for key, parts in qualities.items()},
    )


def _extract_fields(
    seq1: str,
    qual1: str,
    seq2: str | None,
    qual2: str | None,
    config: MapSeqConfig,
) -> tuple[str, str, str, str, str, str, str, str] | None:
    if config.reverse_complement_forward:
        seq1, qual1 = _reverse_complement(seq1), qual1[::-1]
    if seq2 is not None and qual2 is not None and config.reverse_complement_reverse:
        seq2, qual2 = _reverse_complement(seq2), qual2[::-1]
    if config.elements_forward is None and config.elements_reverse is None:
        if seq2 is None or qual2 is None:
            return None
        c0 = config.barcode_length
        c1 = c0 + config.constant_length
        u0 = config.umi_length
        u1 = u0 + config.sample_index_length
        if len(seq1) < c1 or len(seq2) < u1:
            return None
        return (
            seq1[:c0],
            seq2[:u0],
            seq2[u0:u1],
            seq1[c0:c1],
            "",
            qual1[:c0],
            qual2[:u0],
            qual2[u0:u1],
        )
    forward = _extract_composition(
        seq1,
        qual1,
        config.elements_forward or "S",
        config.element_lengths_forward or (-1,),
    )
    reverse = (
        _extract_composition(
            seq2,
            qual2 or "",
            config.elements_reverse or "S",
            config.element_lengths_reverse or (-1,),
        )
        if seq2 is not None
        else ({key: "" for key in "VUCPIS"}, {key: "" for key in "VUCPIS"})
    )
    if forward is None or reverse is None:
        return None
    fseq, fqual = forward
    rseq, rqual = reverse
    return (
        fseq["V"] + rseq["V"],
        fseq["U"] + rseq["U"],
        fseq["I"] + rseq["I"],
        fseq["C"],
        rseq["C"],
        fqual["V"] + rqual["V"],
        fqual["U"] + rqual["U"],
        fqual["I"] + rqual["I"],
    )


def digest_fastqs(
    r1_path: str | Path,
    r2_path: str | Path | None,
    database: str | Path,
    *,
    config: MapSeqConfig = MapSeqConfig(),
    max_reads: int | None = None,
    flush_unique: int = 500_000,
    library_id: str = "library-1",
    lane: str | None = None,
    library_metadata: dict[str, object] | None = None,
) -> dict[str, int | float]:
    """Stream single- or paired-end reads into aggregate and library tables."""
    database = Path(database)
    database.parent.mkdir(parents=True, exist_ok=True)
    con = _connect(database)
    _initialize(con)
    counts: Counter[tuple[str, str, str]] = Counter()
    qc: Counter[str] = Counter()
    started = perf_counter()
    constants_forward = set(x.upper() for x in config.constant_forward)
    constants_reverse = set(x.upper() for x in config.constant_reverse)
    primers_forward = set(x.upper() for x in config.primer_forward)
    primers_reverse = set(x.upper() for x in config.primer_reverse)
    con.execute(
        """
        INSERT OR REPLACE INTO libraries(library_id, r1_path, r2_path, lane, metadata_json)
        VALUES (?, ?, ?, ?, ?)
        """,
        (
            library_id,
            str(Path(r1_path).resolve()),
            str(Path(r2_path).resolve()) if r2_path is not None else None,
            lane,
            json.dumps(library_metadata or {}, sort_keys=True),
        ),
    )

    for _, seq1, qual1, seq2, qual2 in read_sets(r1_path, r2_path):
        if max_reads is not None and qc["total_reads"] >= max_reads:
            break
        qc["total_reads"] += 1
        fields = _extract_fields(seq1, qual1, seq2, qual2, config)
        if fields is None:
            qc["wrong_length"] += 1
            continue
        barcode, umi, sample_index, constant_fwd, constant_rev, barcode_q, umi_q, index_q = fields
        if (
            constants_forward
            and (config.elements_forward is None or "C" in config.elements_forward)
            and constant_fwd not in constants_forward
        ):
            qc["constant_mismatch"] += 1
            continue
        if (
            constants_reverse
            and config.elements_reverse is not None
            and "C" in config.elements_reverse
            and constant_rev not in constants_reverse
        ):
            qc["constant_mismatch_reverse"] += 1
            continue
        if config.elements_forward and "P" in config.elements_forward:
            extracted = _extract_composition(
                _reverse_complement(seq1) if config.reverse_complement_forward else seq1,
                qual1[::-1] if config.reverse_complement_forward else qual1,
                config.elements_forward,
                config.element_lengths_forward,
            )
            if extracted is None or (primers_forward and extracted[0]["P"] not in primers_forward):
                qc["primer_mismatch"] += 1
                continue
        if seq2 is not None and config.elements_reverse and "P" in config.elements_reverse:
            extracted = _extract_composition(
                _reverse_complement(seq2) if config.reverse_complement_reverse else seq2,
                (qual2 or "")[::-1] if config.reverse_complement_reverse else (qual2 or ""),
                config.elements_reverse,
                config.element_lengths_reverse,
            )
            if extracted is None or (primers_reverse and extracted[0]["P"] not in primers_reverse):
                qc["primer_mismatch_reverse"] += 1
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
            (barcode_q and _average_phred(barcode_q) < config.min_average_phred)
            or (
                (umi_q or index_q)
                and _average_phred(umi_q + index_q) < config.min_average_phred
            )
        ):
            qc["low_variable_quality"] += 1
            continue
        counts[(barcode, sample_index, umi)] += 1
        qc["retained_reads"] += 1
        if len(counts) >= flush_unique:
            _flush_library_counts(con, library_id, counts)
            _flush_counts(con, counts)
    if counts:
        _flush_library_counts(con, library_id, counts)
        _flush_counts(con, counts)

    con.executemany(
        "INSERT OR REPLACE INTO qc(key, value) VALUES (?, ?)", sorted(qc.items())
    )
    metadata = {
        "r1_path": str(Path(r1_path).resolve()),
        "r2_path": str(Path(r2_path).resolve()) if r2_path is not None else "",
        "config": json.dumps(asdict(config), sort_keys=True),
        "format_version": "2",
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


def digest_libraries(
    libraries: list[dict[str, object]],
    database: str | Path,
    *,
    config: MapSeqConfig = MapSeqConfig(),
    max_reads_per_library: int | None = None,
) -> dict[str, int | float]:
    """Ingest multiple lanes/libraries without concatenating FASTQ files."""
    if not libraries:
        raise ValueError("libraries must not be empty")
    identifiers = [str(item.get("library_id", "")) for item in libraries]
    if any(not identifier for identifier in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError("every library requires a unique non-empty library_id")
    totals: Counter[str] = Counter()
    elapsed = 0.0
    for library in libraries:
        result = digest_fastqs(
            str(library["r1"]),
            str(library["r2"]) if library.get("r2") else None,
            database,
            config=config,
            max_reads=max_reads_per_library,
            library_id=str(library["library_id"]),
            lane=str(library["lane"]) if library.get("lane") is not None else None,
            library_metadata=dict(library.get("metadata", {})),
        )
        elapsed += float(result["elapsed_seconds"])
        totals.update(
            {key: int(value) for key, value in result.items() if key != "elapsed_seconds"}
        )
    con = _connect(database)
    con.execute("DELETE FROM qc")
    con.executemany(
        "INSERT INTO qc(key, value) VALUES (?, ?)",
        sorted((key, value) for key, value in totals.items() if key != "unique_combinations"),
    )
    unique = con.execute("SELECT count(*) FROM raw_counts").fetchone()[0]
    con.execute(
        "INSERT OR REPLACE INTO run_metadata(key, value) VALUES ('library_manifest', ?)",
        (json.dumps(libraries, sort_keys=True),),
    )
    con.commit()
    con.close()
    return {**totals, "unique_combinations": unique, "libraries": len(libraries), "elapsed_seconds": elapsed}


def merge_databases(
    sources: list[str | Path],
    output: str | Path,
) -> dict[str, int]:
    """Merge independently digested databases with collision-safe library IDs."""
    if not sources:
        raise ValueError("sources must not be empty")
    output_path = Path(output).resolve()
    if any(Path(source).resolve() == output_path for source in sources):
        raise ValueError("output database must differ from every source")
    con = _connect(output_path)
    _initialize(con)
    seen_libraries = {row[0] for row in con.execute("SELECT library_id FROM libraries")}
    merged_qc: Counter[str] = Counter()
    for source in sources:
        source_con = sqlite3.connect(source)
        legacy_id = f"legacy-{Path(source).stem}"
        has_libraries = source_con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='libraries'"
        ).fetchone() is not None
        if has_libraries:
            library_rows = source_con.execute(
                "SELECT library_id, r1_path, r2_path, lane, metadata_json FROM libraries"
            ).fetchall()
        else:
            library_rows = [
                (
                    legacy_id,
                    str(Path(source).resolve()),
                    None,
                    None,
                    json.dumps({"source_format": "legacy-aggregate"}),
                )
            ]
        collisions = seen_libraries.intersection(row[0] for row in library_rows)
        if collisions:
            source_con.close()
            con.close()
            raise ValueError(f"duplicate library identities: {sorted(collisions)}")
        con.executemany("INSERT INTO libraries VALUES (?, ?, ?, ?, ?)", library_rows)
        seen_libraries.update(row[0] for row in library_rows)
        has_library_counts = source_con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='library_counts'"
        ).fetchone() is not None
        if not has_library_counts and len(library_rows) != 1:
            source_con.close()
            con.close()
            raise ValueError("database has multiple libraries but no library_counts table")
        fallback_library_id = library_rows[0][0]
        library_counts = (
            source_con.execute(
                "SELECT library_id, barcode, sample_index, umi, read_count FROM library_counts"
            ).fetchall()
            if has_library_counts
            else [
                (fallback_library_id, barcode, index, umi, count)
                for barcode, index, umi, count in source_con.execute(
                    "SELECT barcode, sample_index, umi, read_count FROM raw_counts"
                )
            ]
        )
        con.executemany("INSERT INTO library_counts VALUES (?, ?, ?, ?, ?)", library_counts)
        con.executemany(
            """
            INSERT INTO raw_counts VALUES (?, ?, ?, ?)
            ON CONFLICT(barcode, sample_index, umi)
            DO UPDATE SET read_count = read_count + excluded.read_count
            """,
            source_con.execute("SELECT barcode, sample_index, umi, read_count FROM raw_counts"),
        )
        if source_con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='qc'"
        ).fetchone():
            merged_qc.update(dict(source_con.execute("SELECT key, value FROM qc")))
        source_con.close()
    con.execute("DELETE FROM qc")
    merged_qc["retained_reads"] = con.execute(
        "SELECT COALESCE(sum(read_count), 0) FROM raw_counts"
    ).fetchone()[0]
    con.executemany(
        "INSERT INTO qc(key, value) VALUES (?, ?)", sorted(merged_qc.items())
    )
    con.execute(
        "INSERT OR REPLACE INTO run_metadata VALUES ('merged_sources', ?)",
        (json.dumps([str(Path(source).resolve()) for source in sources]),),
    )
    con.commit()
    combinations = con.execute("SELECT count(*) FROM raw_counts").fetchone()[0]
    libraries = con.execute("SELECT count(*) FROM libraries").fetchone()[0]
    con.close()
    return {"sources": len(sources), "libraries": libraries, "unique_combinations": combinations}


def collapse_database(
    database: str | Path,
    *,
    collapse_max_dist: float = 1,
    collapse_min_score: float = 2,
    collapse_min_ratio: float = 0,
    min_combo_reads: int = 1,
    distance_metric: str = "hamming",
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
        distance_metric=distance_metric,
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
    has_library_counts = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='library_counts'"
    ).fetchone() is not None
    if has_library_counts:
        library_index_join = (
            "LEFT JOIN sample_index_mapping i ON i.observed_index = r.sample_index"
            if has_index_mapping
            else ""
        )
        library_index_value = (
            "COALESCE(i.sample_index, r.sample_index)" if has_index_mapping else "r.sample_index"
        )
        con.executescript(
            f"""
            DROP TABLE IF EXISTS library_collapsed_counts;
            CREATE TABLE library_collapsed_counts (
                library_id TEXT NOT NULL,
                barcode TEXT NOT NULL,
                sample_index TEXT NOT NULL,
                umi TEXT NOT NULL,
                read_count INTEGER NOT NULL,
                PRIMARY KEY (library_id, barcode, sample_index, umi)
            ) WITHOUT ROWID;
            INSERT INTO library_collapsed_counts
            SELECT r.library_id, m.representative, {library_index_value}, r.umi, sum(r.read_count)
            FROM library_counts r
            JOIN barcode_mapping m USING (barcode)
            {library_index_join}
            WHERE r.read_count >= {int(min_combo_reads)}
            GROUP BY r.library_id, m.representative, {library_index_value}, r.umi;
            """
        )
    parameters = {
        "collapse_max_dist": collapse_max_dist,
        "collapse_min_score": collapse_min_score,
        "collapse_min_ratio": collapse_min_ratio,
        "min_combo_reads": min_combo_reads,
        "distance_metric": distance_metric,
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
    method: str = "equal",
) -> dict[str, int | float]:
    """Collapse UMIs within each barcode/sample-index stratum.

    All UMI scores are equal, reproducing mutscan's lexicographic greedy UMI
    representative choice. The input ``collapsed_counts`` and a complete
    ``umi_mapping`` remain available for audit.
    """
    if method not in {"equal", "directional"}:
        raise ValueError("method must be 'equal' or 'directional'")
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
        if method == "equal":
            mapping = group_similar_sequences(
                umis, [1.0] * len(umis), collapse_max_dist=collapse_max_dist
            )
        else:
            mapping = group_directional_sequences(
                umis,
                {umi: reads for umi, reads in values},
                collapse_max_dist=collapse_max_dist,
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
    has_library_collapsed = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='library_collapsed_counts'"
    ).fetchone() is not None
    if has_library_collapsed:
        con.executescript(
            """
            DROP TABLE IF EXISTS library_umi_collapsed_counts;
            DROP TABLE IF EXISTS library_molecule_counts;
            CREATE TABLE library_umi_collapsed_counts (
                library_id TEXT NOT NULL,
                barcode TEXT NOT NULL,
                sample_index TEXT NOT NULL,
                umi TEXT NOT NULL,
                read_count INTEGER NOT NULL,
                PRIMARY KEY (library_id, barcode, sample_index, umi)
            ) WITHOUT ROWID;
            INSERT INTO library_umi_collapsed_counts
            SELECT l.library_id, l.barcode, l.sample_index, m.representative, sum(l.read_count)
            FROM library_collapsed_counts l
            JOIN umi_mapping m USING (barcode, sample_index, umi)
            GROUP BY l.library_id, l.barcode, l.sample_index, m.representative;
            CREATE TABLE library_molecule_counts (
                library_id TEXT NOT NULL,
                barcode TEXT NOT NULL,
                sample_index TEXT NOT NULL,
                molecule_count INTEGER NOT NULL,
                read_count INTEGER NOT NULL,
                PRIMARY KEY (library_id, barcode, sample_index)
            ) WITHOUT ROWID;
            INSERT INTO library_molecule_counts
            SELECT library_id, barcode, sample_index, count(*), sum(read_count)
            FROM library_umi_collapsed_counts
            GROUP BY library_id, barcode, sample_index;
            """
        )
    con.execute(
        "INSERT OR REPLACE INTO run_metadata(key, value) VALUES ('umi_collapse_parameters', ?)",
        (
            json.dumps(
                {
                    "collapse_max_dist": collapse_max_dist,
                    "score_method": (
                        "equal_lexicographic" if method == "equal" else "directional_2n_minus_1"
                    ),
                }
            ),
        ),
    )
    con.commit()
    con.close()
    return {
        "strata": n_strata,
        "input_umis": n_input,
        "representative_umis": n_representatives,
        "elapsed_seconds": perf_counter() - started,
    }


def import_sample_metadata(
    database: str | Path,
    metadata_path: str | Path,
    *,
    key_column: str = "sample_index",
) -> int:
    """Import arbitrary sample annotations while retaining the original fields."""
    with Path(metadata_path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or key_column not in reader.fieldnames:
            raise ValueError(f"metadata must contain a {key_column!r} column")
        rows = []
        for row in reader:
            key = (row.get(key_column) or "").strip()
            if not key:
                raise ValueError("sample metadata keys must be non-empty")
            rows.append((key, json.dumps(row, sort_keys=True)))
    if len({row[0] for row in rows}) != len(rows):
        raise ValueError("sample metadata keys must be unique")
    con = _connect(database)
    con.executescript(
        """
        CREATE TABLE IF NOT EXISTS sample_metadata (
            sample_index TEXT PRIMARY KEY,
            metadata_json TEXT NOT NULL
        ) WITHOUT ROWID;
        """
    )
    con.executemany(
        "INSERT OR REPLACE INTO sample_metadata(sample_index, metadata_json) VALUES (?, ?)", rows
    )
    con.commit()
    con.close()
    return len(rows)


def export_sparse_matrix(
    database: str | Path,
    output_directory: str | Path,
    *,
    value: str = "molecule_count",
) -> dict[str, int]:
    """Export a barcode-by-sample matrix in Matrix Market coordinate format."""
    if value not in {"molecule_count", "read_count"}:
        raise ValueError("value must be 'molecule_count' or 'read_count'")
    con = sqlite3.connect(database)
    if con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='molecule_counts'"
    ).fetchone() is None:
        con.close()
        raise ValueError("molecule_counts is required; run UMI collapse first")
    barcodes = [row[0] for row in con.execute("SELECT DISTINCT barcode FROM molecule_counts ORDER BY barcode")]
    samples = [row[0] for row in con.execute("SELECT DISTINCT sample_index FROM molecule_counts ORDER BY sample_index")]
    barcode_index = {barcode: index + 1 for index, barcode in enumerate(barcodes)}
    sample_index = {sample: index + 1 for index, sample in enumerate(samples)}
    entries = [
        (barcode_index[barcode], sample_index[sample], int(count))
        for barcode, sample, count in con.execute(
            f"SELECT barcode, sample_index, {value} FROM molecule_counts WHERE {value} != 0"
        )
    ]
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "matrix.mtx").open("w", encoding="ascii", newline="") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% rows=barcodes columns=samples generated by pymutscan\n")
        handle.write(f"{len(barcodes)} {len(samples)} {len(entries)}\n")
        for row, column, count in entries:
            handle.write(f"{row} {column} {count}\n")
    with (output / "barcodes.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("barcode\n" + "\n".join(barcodes) + "\n")
    metadata_exists = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sample_metadata'"
    ).fetchone() is not None
    metadata = (
        dict(con.execute("SELECT sample_index, metadata_json FROM sample_metadata"))
        if metadata_exists
        else {}
    )
    with (output / "samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_index", "metadata_json"])
        for sample in samples:
            writer.writerow([sample, metadata.get(sample, "{}")])
    con.close()
    return {"barcodes": len(barcodes), "samples": len(samples), "nonzero": len(entries)}


def call_template_switch_evidence(
    database: str | Path,
    *,
    min_reads_per_barcode: int = 2,
    min_barcodes: int = 2,
) -> dict[str, int]:
    """Materialize shared-UMI evidence without assigning a causal mechanism."""
    if min_reads_per_barcode < 1 or min_barcodes < 2:
        raise ValueError("support thresholds must be positive and min_barcodes at least two")
    con = _connect(database)
    if con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='collapsed_counts'"
    ).fetchone() is None:
        con.close()
        raise ValueError("collapsed_counts is required")
    groups: dict[tuple[str, str], list[tuple[str, int]]] = {}
    for barcode, sample, umi, reads in con.execute(
        "SELECT barcode, sample_index, umi, read_count FROM collapsed_counts"
    ):
        if reads >= min_reads_per_barcode:
            groups.setdefault((sample, umi), []).append((barcode, int(reads)))
    rows = []
    for (sample, umi), supports in groups.items():
        if len(supports) < min_barcodes:
            continue
        counts = [count for _, count in supports]
        rows.append(
            (
                sample,
                umi,
                len(supports),
                sum(counts),
                min(counts),
                max(counts),
                "shared_umi_high_support" if min(counts) >= 10 else "shared_umi_evidence",
                json.dumps(sorted(supports), separators=(",", ":")),
            )
        )
    con.executescript(
        """
        DROP TABLE IF EXISTS template_switch_evidence;
        CREATE TABLE template_switch_evidence (
            sample_index TEXT NOT NULL,
            umi TEXT NOT NULL,
            n_barcodes INTEGER NOT NULL,
            total_reads INTEGER NOT NULL,
            min_barcode_reads INTEGER NOT NULL,
            max_barcode_reads INTEGER NOT NULL,
            evidence_level TEXT NOT NULL,
            barcode_support_json TEXT NOT NULL,
            PRIMARY KEY (sample_index, umi)
        ) WITHOUT ROWID;
        """
    )
    con.executemany("INSERT INTO template_switch_evidence VALUES (?, ?, ?, ?, ?, ?, ?, ?)", rows)
    con.execute(
        "INSERT OR REPLACE INTO run_metadata VALUES ('template_switch_parameters', ?)",
        (json.dumps({"min_reads_per_barcode": min_reads_per_barcode, "min_barcodes": min_barcodes}),),
    )
    con.commit()
    con.close()
    return {"candidates": len(rows)}


def export_table(database: str | Path, table: str, output: str | Path) -> int:
    """Export an allowed database table as TSV or TSV.GZ."""
    allowed = {"raw_counts", "library_counts", "libraries", "barcode_scores", "barcode_mapping", "sample_index_mapping", "collapsed_counts", "library_collapsed_counts", "umi_mapping", "umi_collapsed_counts", "library_umi_collapsed_counts", "molecule_counts", "library_molecule_counts", "sample_metadata", "template_switch_evidence", "qc", "run_metadata"}
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
