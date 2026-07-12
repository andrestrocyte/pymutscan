"""Sequence grouping algorithms compatible with mutscan's greedy semantics."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from math import floor

DNA = "ACGT"


def hamming_distance(left: str, right: str) -> int:
    """Return Hamming distance for two equal-length strings."""
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(a != b for a, b in zip(left, right))


def _integer_tolerance(max_distance: float, length: int) -> int:
    # The C++ implementation passes a double tolerance into a function taking
    # int, so fractional values are truncated rather than rounded.
    return int(max_distance if max_distance >= 1 else floor(max_distance * length))


def _radius_one_candidates(sequence: str) -> Iterable[str]:
    yield sequence
    chars = list(sequence)
    for pos, original in enumerate(chars):
        for base in DNA:
            if base != original:
                chars[pos] = base
                yield "".join(chars)
        chars[pos] = original


def group_similar_sequences(
    sequences: Sequence[str],
    scores: Sequence[float] | Mapping[str, float],
    *,
    collapse_max_dist: float = 0.0,
    collapse_min_score: float = 0.0,
    collapse_min_ratio: float = 0.0,
) -> dict[str, str]:
    """Map sequences to greedy abundance-ordered representatives.

    Radius zero and one use exact hash lookup. Larger radii intentionally use a
    correctness-first fallback and are not the optimized MAPseq path.
    """
    if not sequences:
        return {}
    normalized = [str(s).upper() for s in sequences]
    length = len(normalized[0])
    if any(len(s) != length for s in normalized):
        raise ValueError("all sequences must have the same length")
    if any(set(s) - set(DNA) for s in normalized):
        raise ValueError("sequences must contain only A, C, G, and T")

    if isinstance(scores, Mapping):
        score_by_sequence = {s.upper(): float(v) for s, v in scores.items()}
    else:
        if len(scores) != len(normalized):
            raise ValueError("sequences and scores must have equal length")
        # mutscan's std::map::insert keeps the first score for duplicate keys.
        score_by_sequence: dict[str, float] = {}
        for sequence, score in zip(normalized, scores):
            score_by_sequence.setdefault(sequence, float(score))

    unique = set(normalized)
    missing = unique.difference(score_by_sequence)
    if missing:
        raise ValueError(f"missing scores for {len(missing)} sequences")
    ordered = sorted(unique, key=lambda s: (-score_by_sequence[s], s))
    tolerance = _integer_tolerance(collapse_max_dist, length)
    active = set(unique)
    representative: dict[str, str] = {}

    for query in ordered:
        if query not in active:
            continue
        if collapse_min_score > 0 and score_by_sequence[query] < collapse_min_score:
            for sequence in active:
                representative[sequence] = sequence
            break

        if tolerance == 0:
            candidates = (query,)
        elif tolerance == 1:
            candidates = _radius_one_candidates(query)
        else:
            candidates = (
                sequence for sequence in active
                if hamming_distance(query, sequence) <= tolerance
            )

        claimed: list[str] = []
        for candidate in candidates:
            if candidate not in active:
                continue
            if (
                candidate == query
                or score_by_sequence[query]
                >= collapse_min_ratio * score_by_sequence[candidate]
            ):
                representative[candidate] = query
                claimed.append(candidate)
        for candidate in claimed:
            active.remove(candidate)

    return {sequence: representative[sequence] for sequence in normalized}

