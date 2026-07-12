"""Small, allocation-conscious paired FASTQ reader."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO


@contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            yield handle
    else:
        with path.open("rt", encoding="ascii", newline="") as handle:
            yield handle


def records(path: str | Path) -> Iterator[tuple[str, str, str]]:
    with open_text(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ record in {path}")
            if not name.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record in {path}")
            if len(sequence) != len(quality):
                raise ValueError(f"sequence/quality length mismatch in {path}")
            yield name.rstrip("\r\n"), sequence.upper(), quality


def paired_records(
    r1_path: str | Path, r2_path: str | Path
) -> Iterator[tuple[str, str, str, str, str]]:
    r1_iter = records(r1_path)
    r2_iter = records(r2_path)
    sentinel = object()
    while True:
        r1 = next(r1_iter, sentinel)
        r2 = next(r2_iter, sentinel)
        if r1 is sentinel and r2 is sentinel:
            return
        if r1 is sentinel or r2 is sentinel:
            raise ValueError("paired FASTQ files contain different record counts")
        name1, seq1, qual1 = r1
        name2, seq2, qual2 = r2
        if name1.split()[0] != name2.split()[0]:
            raise ValueError(f"paired read names differ: {name1} / {name2}")
        yield name1, seq1, qual1, seq2, qual2


def read_sets(
    r1_path: str | Path,
    r2_path: str | Path | None = None,
) -> Iterator[tuple[str, str, str, str | None, str | None]]:
    """Yield single-end records or synchronized paired-end records."""
    if r2_path is None:
        for name, sequence, quality in records(r1_path):
            yield name, sequence, quality, None, None
        return
    yield from paired_records(r1_path, r2_path)
