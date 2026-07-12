"""Validated experiment configuration and named MAPseq presets."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10 CI
    import tomli as tomllib


@dataclass(frozen=True)
class MapSeqConfig:
    """Read layout and filtering configuration.

    The legacy length fields remain authoritative unless an explicit element
    composition is supplied. Composition symbols follow mutscan: ``V``
    variable/barcode, ``U`` UMI, ``C`` constant, ``P`` primer, and ``S`` skip.
    pymutscan adds ``I`` for a sample-index segment.
    """

    barcode_length: int = 30
    umi_length: int = 16
    sample_index_length: int = 6
    constant_forward: tuple[str, ...] = ("CCGTACT", "CTGTACT", "TCGTACT", "TTGTACT")
    min_average_phred: float = 20.0
    max_ambiguous_barcode: int = 0
    max_ambiguous_umi: int = 0
    max_ambiguous_sample_index: int = 0
    constant_reverse: tuple[str, ...] = ()
    primer_forward: tuple[str, ...] = ()
    primer_reverse: tuple[str, ...] = ()
    elements_forward: str | None = None
    element_lengths_forward: tuple[int, ...] = ()
    elements_reverse: str | None = None
    element_lengths_reverse: tuple[int, ...] = ()
    reverse_complement_forward: bool = False
    reverse_complement_reverse: bool = False

    def __post_init__(self) -> None:
        for name in ("barcode_length", "umi_length", "sample_index_length"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative")
        self._validate_composition(self.elements_forward, self.element_lengths_forward, "forward")
        self._validate_composition(self.elements_reverse, self.element_lengths_reverse, "reverse")
        for values in (
            self.constant_forward,
            self.constant_reverse,
            self.primer_forward,
            self.primer_reverse,
        ):
            if any(set(value.upper()) - set("ACGT") for value in values):
                raise ValueError("constant and primer sequences must contain only A, C, G, and T")

    @staticmethod
    def _validate_composition(elements: str | None, lengths: tuple[int, ...], read: str) -> None:
        if elements is None:
            if lengths:
                raise ValueError(f"{read} element lengths require an element string")
            return
        if set(elements.upper()) - set("VUCPSI"):
            raise ValueError("composition elements must use V, U, C, P, S, or I")
        if len(elements) != len(lengths):
            raise ValueError(f"{read} elements and lengths must have equal length")
        if any(length < -1 for length in lengths) or lengths.count(-1) > 1:
            raise ValueError("element lengths must be non-negative, with at most one trailing -1")
        if -1 in lengths and lengths[-1] != -1:
            raise ValueError("a -1 remainder element must be last")

    @property
    def constant_length(self) -> int:
        lengths = {len(value) for value in self.constant_forward}
        if len(lengths) != 1:
            raise ValueError("all forward constant sequences must have equal length")
        return next(iter(lengths))

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, values: dict[str, Any]) -> MapSeqConfig:
        allowed = {field.name for field in fields(cls)}
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"unknown configuration keys: {sorted(unknown)}")
        tuple_fields = {
            "constant_forward",
            "constant_reverse",
            "primer_forward",
            "primer_reverse",
            "element_lengths_forward",
            "element_lengths_reverse",
        }
        normalized = {
            key: tuple(value) if key in tuple_fields and value is not None else value
            for key, value in values.items()
        }
        return cls(**normalized)

    @classmethod
    def from_file(cls, path: str | Path, *, preset: str | None = None) -> MapSeqConfig:
        path = Path(path)
        with path.open("rb") as handle:
            data = json.load(handle) if path.suffix.lower() == ".json" else tomllib.load(handle)
        if preset is not None:
            try:
                data = data["presets"][preset]
            except KeyError as error:
                raise ValueError(f"preset {preset!r} not found in {path}") from error
        elif "config" in data:
            data = data["config"]
        return cls.from_dict(data)


PRESETS: dict[str, MapSeqConfig] = {
    "mapseq-default": MapSeqConfig(),
    "mapseq-30-18-14": MapSeqConfig(umi_length=18, sample_index_length=14),
    "mapseq-template-switch": MapSeqConfig(
        elements_forward="VCS",
        element_lengths_forward=(30, 7, -1),
        elements_reverse="UIS",
        element_lengths_reverse=(18, 14, -1),
    ),
}


def get_preset(name: str) -> MapSeqConfig:
    try:
        return PRESETS[name]
    except KeyError as error:
        raise ValueError(f"unknown preset {name!r}; choose from {sorted(PRESETS)}") from error


def load_config(path: str | Path | None = None, *, preset: str | None = None) -> MapSeqConfig:
    if path is not None:
        return MapSeqConfig.from_file(path, preset=preset)
    return get_preset(preset) if preset else MapSeqConfig()


def load_library_manifest(path: str | Path) -> list[dict[str, object]]:
    """Load ``libraries`` from a JSON or TOML ingestion manifest."""
    path = Path(path)
    with path.open("rb") as handle:
        data = json.load(handle) if path.suffix.lower() == ".json" else tomllib.load(handle)
    libraries = data.get("libraries")
    if not isinstance(libraries, list) or not all(isinstance(item, dict) for item in libraries):
        raise ValueError("manifest must contain a list named 'libraries'")
    base = path.parent
    normalized = []
    for item in libraries:
        library = dict(item)
        if "r1" not in library:
            raise ValueError("each manifest library requires r1")
        for key in ("r1", "r2"):
            if library.get(key) and not Path(str(library[key])).is_absolute():
                library[key] = str((base / str(library[key])).resolve())
        normalized.append(library)
    return normalized
