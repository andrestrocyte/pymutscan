"""Fast MAPseq preprocessing and exact radius-1 barcode collapsing."""

from .collapse import group_similar_sequences
from .pipeline import (
    MapSeqConfig,
    collapse_database,
    collapse_umis,
    digest_fastqs,
    export_table,
    map_sample_indices,
)

__all__ = [
    "MapSeqConfig",
    "collapse_database",
    "collapse_umis",
    "digest_fastqs",
    "export_table",
    "group_similar_sequences",
    "map_sample_indices",
]

__version__ = "0.2.0"
