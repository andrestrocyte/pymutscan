#!/usr/bin/env python3
"""Reproducible search benchmarks with correctness and peak-memory checks."""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import random
import statistics
import tracemalloc
from pathlib import Path
from time import perf_counter

from pymutscan import group_similar_sequences
from pymutscan.collapse import hamming_distance

DNA = "ACGT"


def dataset(size: int, length: int = 30, seed: int = 20260712) -> tuple[list[str], list[int]]:
    rng = random.Random(seed + size)
    sequences: list[str] = []
    seen: set[str] = set()
    while len(sequences) < size:
        sequence = "".join(rng.choice(DNA) for _ in range(length))
        if sequence not in seen:
            seen.add(sequence)
            sequences.append(sequence)
    # Add exact known radius-two relationships without changing cardinality.
    for index in range(1, min(size, 100), 10):
        parent = sequences[index - 1]
        child = list(parent)
        child[0] = DNA[(DNA.index(child[0]) + 1) % 4]
        child[7] = DNA[(DNA.index(child[7]) + 1) % 4]
        candidate = "".join(child)
        if candidate not in seen:
            seen.remove(sequences[index])
            sequences[index] = candidate
            seen.add(candidate)
    scores = list(range(size, 0, -1))
    return sequences, scores


def exhaustive(sequences: list[str], scores: list[int], radius: int) -> dict[str, str]:
    score_map = dict(zip(sequences, scores))
    active = set(sequences)
    mapping: dict[str, str] = {}
    for query in sorted(active, key=lambda sequence: (-score_map[sequence], sequence)):
        if query not in active:
            continue
        claimed = [
            candidate
            for candidate in active
            if hamming_distance(query, candidate) <= radius
        ]
        for candidate in claimed:
            mapping[candidate] = query
            active.remove(candidate)
    return mapping


def timed(callable_, repeats: int) -> tuple[float, object]:
    values = []
    result = None
    for _ in range(repeats):
        started = perf_counter()
        result = callable_()
        values.append(perf_counter() - started)
    return statistics.median(values), result


def peak_memory(callable_) -> int:
    tracemalloc.start()
    callable_()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", type=int, nargs="+", default=[128, 256, 384, 500, 2000, 5000])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--profile-size", type=int, default=20000)
    parser.add_argument("--output", type=Path, default=Path("benchmarks_public/search_benchmark.json"))
    parser.add_argument(
        "--dataset-output",
        type=Path,
        default=Path("benchmarks_public/radius2_validation_500.tsv.gz"),
    )
    args = parser.parse_args()
    records = []
    for size in args.sizes:
        sequences, scores = dataset(size)
        optimized_seconds, optimized_mapping = timed(
            lambda: group_similar_sequences(
                sequences, scores, collapse_max_dist=2, collapse_min_score=0
            ),
            args.repeats,
        )
        exhaustive_seconds, exhaustive_mapping = timed(
            lambda: exhaustive(sequences, scores, 2), args.repeats
        )
        if optimized_mapping != exhaustive_mapping:
            raise RuntimeError("optimized and exhaustive mappings differ")
        optimized_peak = peak_memory(
            lambda: group_similar_sequences(
                sequences, scores, collapse_max_dist=2, collapse_min_score=0
            )
        )
        records.append(
            {
                "unique_sequences": size,
                "length": len(sequences[0]),
                "radius": 2,
                "search_strategy": "adaptive_exhaustive" if size < 384 else "packed_xor",
                "optimized_seconds_median": optimized_seconds,
                "exhaustive_seconds_median": exhaustive_seconds,
                "speedup": exhaustive_seconds / optimized_seconds,
                "optimized_peak_bytes": optimized_peak,
                "mappings_identical": True,
            }
        )
        print(json.dumps(records[-1], sort_keys=True))
    profile_sequences, profile_scores = dataset(args.profile_size)
    profile_seconds, _ = timed(
        lambda: group_similar_sequences(
            profile_sequences, profile_scores, collapse_max_dist=2, collapse_min_score=0
        ),
        1,
    )
    profile_peak = peak_memory(
        lambda: group_similar_sequences(
            profile_sequences, profile_scores, collapse_max_dist=2, collapse_min_score=0
        )
    )
    large_scale_profile = {
        "unique_sequences": args.profile_size,
        "length": len(profile_sequences[0]),
        "radius": 2,
        "optimized_seconds": profile_seconds,
        "optimized_peak_bytes": profile_peak,
        "exhaustive_skipped": True,
        "reason": "quadratic reference is retained at validation scales only",
    }
    print(json.dumps(large_scale_profile, sort_keys=True))
    payload = {
        "schema_version": 1,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "seed": 20260712,
        "repeats": args.repeats,
        "results": records,
        "large_scale_profile": large_scale_profile,
        "notes": (
            "Synthetic fixed-length DNA; exact mapping equality is asserted for every validation "
            "scale. tracemalloc reports incremental grouping allocations after inputs and the "
            "immutable mask cache are constructed."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    validation_sequences, validation_scores = dataset(500)
    args.dataset_output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.dataset_output, "wt", encoding="ascii", newline="") as handle:
        handle.write("sequence\tscore\n")
        for sequence, score in zip(validation_sequences, validation_scores):
            handle.write(f"{sequence}\t{score}\n")


if __name__ == "__main__":
    main()
