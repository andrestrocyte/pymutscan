"""Sequence grouping algorithms compatible with mutscan's greedy semantics."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from functools import lru_cache
from math import floor

DNA = "ACGT"
_DNA_CODE = {base: code for code, base in enumerate(DNA)}
RADIUS_TWO_PACKED_MIN_SEQUENCES = 384


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


def _radius_two_candidates(sequence: str) -> Iterable[str]:
    """Yield every DNA string at Hamming distance at most two exactly once."""
    yield from _radius_one_candidates(sequence)
    chars = list(sequence)
    for left in range(len(chars) - 1):
        original_left = chars[left]
        for right in range(left + 1, len(chars)):
            original_right = chars[right]
            for left_base in DNA:
                if left_base == original_left:
                    continue
                chars[left] = left_base
                for right_base in DNA:
                    if right_base == original_right:
                        continue
                    chars[right] = right_base
                    yield "".join(chars)
                chars[right] = original_right
            chars[left] = original_left


def _encode_dna(sequence: str) -> int:
    value = 0
    for base in sequence:
        value = (value << 2) | _DNA_CODE[base]
    return value


@lru_cache(maxsize=None)
def _hamming_masks(length: int, radius: int) -> tuple[int, ...]:
    singles = tuple(
        delta << (2 * (length - position - 1))
        for position in range(length)
        for delta in (1, 2, 3)
    )
    if radius == 1:
        return (0, *singles)
    doubles = tuple(
        (left_delta << (2 * (length - left - 1)))
        ^ (right_delta << (2 * (length - right - 1)))
        for left in range(length - 1)
        for right in range(left + 1, length)
        for left_delta in (1, 2, 3)
        for right_delta in (1, 2, 3)
    )
    return (0, *singles, *doubles)


def _observed_neighbors(
    sequence: str,
    radius: int,
    observed_by_code: Mapping[int, str],
) -> Iterable[str]:
    """Probe radius-one/two neighbors as packed integers without string allocation."""
    encoded = _encode_dna(sequence)
    for mask in _hamming_masks(len(sequence), radius):
        candidate = observed_by_code.get(encoded ^ mask)
        if candidate is not None:
            yield candidate


def edit_distance(left: str, right: str, max_distance: int | None = None) -> int:
    """Return Levenshtein distance, optionally stopping beyond a threshold."""
    if max_distance is not None and abs(len(left) - len(right)) > max_distance:
        return max_distance + 1
    if len(left) > len(right):
        left, right = right, left
    previous = list(range(len(left) + 1))
    for row, right_base in enumerate(right, 1):
        current = [row]
        row_minimum = row
        for column, left_base in enumerate(left, 1):
            value = min(
                current[-1] + 1,
                previous[column] + 1,
                previous[column - 1] + (left_base != right_base),
            )
            current.append(value)
            row_minimum = min(row_minimum, value)
        if max_distance is not None and row_minimum > max_distance:
            return max_distance + 1
        previous = current
    return previous[-1]


def group_similar_sequences(
    sequences: Sequence[str],
    scores: Sequence[float] | Mapping[str, float],
    *,
    collapse_max_dist: float = 0.0,
    collapse_min_score: float = 0.0,
    collapse_min_ratio: float = 0.0,
    distance_metric: str = "hamming",
) -> dict[str, str]:
    """Map sequences to greedy abundance-ordered representatives.

    Hamming radii zero, one, and two use exact packed-key lookup. Larger radii
    and Levenshtein distance use a correctness-first exhaustive fallback.
    """
    if not sequences:
        return {}
    if distance_metric not in {"hamming", "levenshtein"}:
        raise ValueError("distance_metric must be 'hamming' or 'levenshtein'")
    normalized = [str(s).upper() for s in sequences]
    length = len(normalized[0])
    if distance_metric == "hamming" and any(len(s) != length for s in normalized):
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
    packed_lookup = distance_metric == "hamming" and (
        tolerance == 1
        or (tolerance == 2 and len(unique) >= RADIUS_TWO_PACKED_MIN_SEQUENCES)
    )
    observed_by_code = (
        {_encode_dna(sequence): sequence for sequence in unique}
        if packed_lookup
        else {}
    )
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
        elif packed_lookup:
            candidates = _observed_neighbors(query, tolerance, observed_by_code)
        else:
            candidates = (
                sequence
                for sequence in active
                if (
                    hamming_distance(query, sequence)
                    if distance_metric == "hamming"
                    else edit_distance(query, sequence, tolerance)
                )
                <= tolerance
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


def group_directional_sequences(
    sequences: Sequence[str],
    scores: Sequence[float] | Mapping[str, float],
    *,
    collapse_max_dist: float = 1,
) -> dict[str, str]:
    """Collapse lower-count sequences using the directional UMI criterion.

    A representative may claim a neighbor only when ``parent >= 2*child - 1``,
    the directional adjacency rule popularized for UMI error correction.
    """
    if isinstance(scores, Mapping):
        score_map = {str(key).upper(): float(value) for key, value in scores.items()}
    else:
        if len(sequences) != len(scores):
            raise ValueError("sequences and scores must have equal length")
        score_map = {}
        for sequence, score in zip(sequences, scores):
            score_map.setdefault(str(sequence).upper(), float(score))
    normalized = [str(sequence).upper() for sequence in sequences]
    if not normalized:
        return {}
    length = len(normalized[0])
    if any(len(sequence) != length for sequence in normalized):
        raise ValueError("all sequences must have the same length")
    if any(set(sequence) - set(DNA) for sequence in normalized):
        raise ValueError("sequences must contain only A, C, G, and T")
    missing = set(normalized).difference(score_map)
    if missing:
        raise ValueError(f"missing scores for {len(missing)} sequences")
    tolerance = _integer_tolerance(collapse_max_dist, length)
    if tolerance > 2:
        raise ValueError("directional correction supports Hamming radius zero, one, or two")
    active = set(normalized)
    packed_lookup = tolerance == 1 or (
        tolerance == 2 and len(active) >= RADIUS_TWO_PACKED_MIN_SEQUENCES
    )
    observed_by_code = (
        {_encode_dna(sequence): sequence for sequence in active} if packed_lookup else {}
    )
    result: dict[str, str] = {}
    for query in sorted(active, key=lambda sequence: (-score_map[sequence], sequence)):
        if query not in active:
            continue
        result[query] = query
        active.remove(query)
        frontier = [query]
        while frontier:
            parent = frontier.pop()
            candidates = (
                ()
                if tolerance == 0
                else _observed_neighbors(parent, tolerance, observed_by_code)
                if packed_lookup
                else (
                    candidate
                    for candidate in tuple(active)
                    if hamming_distance(parent, candidate) <= tolerance
                )
            )
            for candidate in candidates:
                if candidate not in active:
                    continue
                if score_map[parent] >= 2 * score_map[candidate] - 1:
                    result[candidate] = query
                    active.remove(candidate)
                    frontier.append(candidate)
    return {sequence: result[sequence] for sequence in normalized}
