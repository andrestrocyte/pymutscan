"""Command-line interface for pymutscan."""

from __future__ import annotations

import argparse
import json

from .config import PRESETS, load_config, load_library_manifest
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


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pymutscan")
    commands = parser.add_subparsers(dest="command", required=True)
    digest = commands.add_parser("digest", help="extract barcode/index/UMI counts")
    digest.add_argument("--r1", required=True)
    digest.add_argument("--r2")
    digest.add_argument("--database", required=True)
    digest.add_argument("--config")
    digest.add_argument("--preset", choices=sorted(PRESETS))
    digest.add_argument("--library-id", default="library-1")
    digest.add_argument("--lane")
    digest.add_argument("--barcode-length", type=int, default=30)
    digest.add_argument("--umi-length", type=int, default=16)
    digest.add_argument("--sample-index-length", type=int, default=6)
    digest.add_argument("--constant-forward", action="append")
    digest.add_argument("--min-average-phred", type=float, default=20)
    digest.add_argument("--max-reads", type=int)

    ingest = commands.add_parser("ingest", help="ingest a JSON/TOML multi-library manifest")
    ingest.add_argument("--manifest", required=True)
    ingest.add_argument("--database", required=True)
    ingest.add_argument("--config")
    ingest.add_argument("--preset", choices=sorted(PRESETS))
    ingest.add_argument("--max-reads-per-library", type=int)

    presets = commands.add_parser("presets", help="list named experiment presets")
    presets.add_argument("--json", action="store_true")

    collapse = commands.add_parser("collapse", help="collapse barcodes and retain index/UMI")
    collapse.add_argument("--database", required=True)
    collapse.add_argument("--max-distance", type=float, default=1)
    collapse.add_argument("--min-score", type=float, default=2)
    collapse.add_argument("--min-ratio", type=float, default=0)
    collapse.add_argument("--min-combo-reads", type=int, default=1)
    collapse.add_argument("--distance-metric", choices=["hamming", "levenshtein"], default="hamming")

    indices = commands.add_parser("map-indices", help="correct RT sample indices against a whitelist")
    indices.add_argument("--database", required=True)
    indices.add_argument("--sample-index", action="append", required=True)
    indices.add_argument("--max-distance", type=int, default=1)

    umis = commands.add_parser("collapse-umis", help="collapse UMIs within barcode/index strata")
    umis.add_argument("--database", required=True)
    umis.add_argument("--max-distance", type=float, default=1)
    umis.add_argument("--method", choices=["equal", "directional"], default="equal")

    export = commands.add_parser("export", help="export a database table as TSV")
    export.add_argument("--database", required=True)
    export.add_argument("--table", required=True)
    export.add_argument("--output", required=True)

    merge = commands.add_parser("merge", help="merge independently digested SQLite databases")
    merge.add_argument("--source", action="append", required=True)
    merge.add_argument("--output", required=True)

    metadata = commands.add_parser("import-metadata", help="import sample metadata TSV")
    metadata.add_argument("--database", required=True)
    metadata.add_argument("--input", required=True)
    metadata.add_argument("--key-column", default="sample_index")

    sparse = commands.add_parser("export-sparse", help="export Matrix Market barcode-by-sample data")
    sparse.add_argument("--database", required=True)
    sparse.add_argument("--output-directory", required=True)
    sparse.add_argument("--value", choices=["molecule_count", "read_count"], default="molecule_count")

    switches = commands.add_parser("template-switch", help="materialize non-causal shared-UMI evidence")
    switches.add_argument("--database", required=True)
    switches.add_argument("--min-reads-per-barcode", type=int, default=2)
    switches.add_argument("--min-barcodes", type=int, default=2)
    return parser


def main() -> None:
    args = _parser().parse_args()
    if args.command == "digest":
        config = load_config(args.config, preset=args.preset)
        if not args.config and not args.preset:
            values = config.to_dict()
            values.update(
                barcode_length=args.barcode_length,
                umi_length=args.umi_length,
                sample_index_length=args.sample_index_length,
                min_average_phred=args.min_average_phred,
            )
            if args.constant_forward:
                values["constant_forward"] = tuple(args.constant_forward)
            config = type(config).from_dict(values)
        result = digest_fastqs(
            args.r1,
            args.r2,
            args.database,
            config=config,
            max_reads=args.max_reads,
            library_id=args.library_id,
            lane=args.lane,
        )
    elif args.command == "ingest":
        result = digest_libraries(
            load_library_manifest(args.manifest),
            args.database,
            config=load_config(args.config, preset=args.preset),
            max_reads_per_library=args.max_reads_per_library,
        )
    elif args.command == "presets":
        payload = {name: config.to_dict() for name, config in PRESETS.items()}
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            print("\n".join(sorted(payload)))
        return
    elif args.command == "map-indices":
        result = map_sample_indices(
            args.database, args.sample_index, max_distance=args.max_distance
        )
    elif args.command == "collapse":
        result = collapse_database(
            args.database,
            collapse_max_dist=args.max_distance,
            collapse_min_score=args.min_score,
            collapse_min_ratio=args.min_ratio,
            min_combo_reads=args.min_combo_reads,
            distance_metric=args.distance_metric,
        )
    elif args.command == "collapse-umis":
        result = collapse_umis(
            args.database, collapse_max_dist=args.max_distance, method=args.method
        )
    elif args.command == "export":
        result = {"rows_exported": export_table(args.database, args.table, args.output)}
    elif args.command == "merge":
        result = merge_databases(args.source, args.output)
    elif args.command == "import-metadata":
        result = {
            "samples_imported": import_sample_metadata(
                args.database, args.input, key_column=args.key_column
            )
        }
    elif args.command == "export-sparse":
        result = export_sparse_matrix(
            args.database, args.output_directory, value=args.value
        )
    else:
        result = call_template_switch_evidence(
            args.database,
            min_reads_per_barcode=args.min_reads_per_barcode,
            min_barcodes=args.min_barcodes,
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
