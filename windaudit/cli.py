"""Command line interface: ``windaudit run --inputs DIR --out DIR``."""

from __future__ import annotations

import argparse
import json
import os
import sys

from .report import run_audit

DEFAULT_NAMES = {
    "relative": "relative_windings.json",
    "same": "same_windings.json",
    "absolute": "abs_winding.json",
    "umbilicus": "umbilicus.json",
}


def _resolve_inputs(args) -> dict:
    paths = {}
    for key, default_name in DEFAULT_NAMES.items():
        explicit = getattr(args, key)
        if explicit:
            path = explicit
        elif args.inputs:
            path = os.path.join(args.inputs, default_name)
        else:
            raise SystemExit(f"windaudit: missing --{key} (or --inputs directory)")
        if not os.path.isfile(path):
            raise SystemExit(f"windaudit: input not found: {path}")
        paths[key] = path
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="windaudit",
        description="Pre-fit consistency auditing of winding annotations.",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run_p = sub.add_parser("run", help="run the full audit and write reports")
    run_p.add_argument("--inputs", help="directory with relative_windings.json, "
                       "same_windings.json, abs_winding.json and umbilicus.json")
    run_p.add_argument("--relative", help="relative-winding point collections")
    run_p.add_argument("--same", help="same-winding point collections")
    run_p.add_argument("--absolute", help="absolute-winding point collections")
    run_p.add_argument("--umbilicus", help="umbilicus control points")
    run_p.add_argument("--out", required=True, help="output directory")
    run_p.add_argument("--spiral-sense", default="CW", choices=["CW", "ACW"],
                       help="upstream spiral_outward_sense of the scan (default CW)")
    run_p.add_argument("--seed", type=int, default=20260917)
    run_p.add_argument("--shuffle-rounds", type=int, default=20)
    run_p.add_argument("--null-rounds", type=int, default=10)
    run_p.add_argument("--skip-planted", action="store_true",
                       help="skip planted-defect validation")
    run_p.add_argument("--skip-sweep", action="store_true",
                       help="skip the sensitivity sweep")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    paths = _resolve_inputs(args)
    report = run_audit(
        relative_path=paths["relative"],
        same_path=paths["same"],
        absolute_path=paths["absolute"],
        umbilicus_path=paths["umbilicus"],
        out_dir=args.out,
        spiral_sense=args.spiral_sense,
        seed=args.seed,
        n_shuffle_rounds=args.shuffle_rounds,
        n_null_rounds=args.null_rounds,
        run_planted=not args.skip_planted,
        run_sweep=not args.skip_sweep,
    )
    corpus = report["corpus"]
    summary = {
        "points_audited": sum(corpus[k]["n_points"] for k in ("relative", "same", "absolute")),
        "order_pairs_tested": report["tier1"]["order_audit"]["n_pairs_tested"],
        "order_violations": report["tier1"]["order_audit"]["n_violations"],
        "same_sheet_links": report["matching"]["n_links"],
        "conflicting_pairs": report["matching"]["n_conflict_pairs"],
        "effective_cycles": report["graph"]["n_effective_cycles"],
        "failing_cycles": report["graph"]["n_effective_nonzero_holonomy"],
        "min_edge_cut": report["repair"]["min_edge_cut"]["cut_edges"],
        "certified_links": report["certified_links"]["n_links"],
        "runtime_seconds": report["runtime_seconds"],
    }
    json.dump(summary, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
