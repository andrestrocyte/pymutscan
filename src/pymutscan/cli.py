"""Command-line interface for pymutscan."""

from __future__ import annotations

import argparse
import json

from .pipeline import (
    MapSeqConfig,
    collapse_database,
    collapse_umis,
    digest_fastqs,
    export_table,
    map_sample_indices,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pymutscan")
    commands = parser.add_subparsers(dest="command", required=True)
    digest = commands.add_parser("digest", help="extract barcode/index/UMI counts")
    digest.add_argument("--r1", required=True)
    digest.add_argument("--r2", required=True)
    digest.add_argument("--database", required=True)
    digest.add_argument("--barcode-length", type=int, default=30)
    digest.add_argument("--umi-length", type=int, default=16)
    digest.add_argument("--sample-index-length", type=int, default=6)
    digest.add_argument("--constant-forward", action="append")
    digest.add_argument("--min-average-phred", type=float, default=20)
    digest.add_argument("--max-reads", type=int)

    collapse = commands.add_parser("collapse", help="collapse barcodes and retain index/UMI")
    collapse.add_argument("--database", required=True)
    collapse.add_argument("--max-distance", type=float, default=1)
    collapse.add_argument("--min-score", type=float, default=2)
    collapse.add_argument("--min-ratio", type=float, default=0)
    collapse.add_argument("--min-combo-reads", type=int, default=1)

    indices = commands.add_parser("map-indices", help="correct RT sample indices against a whitelist")
    indices.add_argument("--database", required=True)
    indices.add_argument("--sample-index", action="append", required=True)
    indices.add_argument("--max-distance", type=int, default=1)

    umis = commands.add_parser("collapse-umis", help="collapse UMIs within barcode/index strata")
    umis.add_argument("--database", required=True)
    umis.add_argument("--max-distance", type=float, default=1)

    export = commands.add_parser("export", help="export a database table as TSV")
    export.add_argument("--database", required=True)
    export.add_argument("--table", required=True)
    export.add_argument("--output", required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    if args.command == "digest":
        constants = tuple(args.constant_forward) if args.constant_forward else MapSeqConfig().constant_forward
        config = MapSeqConfig(
            barcode_length=args.barcode_length,
            umi_length=args.umi_length,
            sample_index_length=args.sample_index_length,
            constant_forward=constants,
            min_average_phred=args.min_average_phred,
        )
        result = digest_fastqs(
            args.r1, args.r2, args.database, config=config, max_reads=args.max_reads
        )
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
        )
    elif args.command == "collapse-umis":
        result = collapse_umis(args.database, collapse_max_dist=args.max_distance)
    else:
        result = {"rows_exported": export_table(args.database, args.table, args.output)}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
