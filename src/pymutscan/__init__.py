"""Fast MAPseq preprocessing and exact radius-one/two barcode collapsing."""

from .collapse import edit_distance, group_directional_sequences, group_similar_sequences
from .config import MapSeqConfig, get_preset, load_config, load_library_manifest
from .pipeline import (
    call_template_switch_evidence,
    collapse_database,
    collapse_umis,
    digest_fastqs,
    digest_libraries,
    export_sparse_matrix,
    export_table,
    import_sample_metadata,
    map_sample_indices,
    merge_databases,
)

__all__ = [
    "MapSeqConfig",
    "call_template_switch_evidence",
    "collapse_database",
    "collapse_umis",
    "digest_fastqs",
    "digest_libraries",
    "edit_distance",
    "export_sparse_matrix",
    "export_table",
    "get_preset",
    "group_directional_sequences",
    "group_similar_sequences",
    "import_sample_metadata",
    "load_config",
    "load_library_manifest",
    "map_sample_indices",
    "merge_databases",
]

__version__ = "0.3.0"
