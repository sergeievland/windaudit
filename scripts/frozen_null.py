"""Measure the false-alarm rate of the frozen configuration with a large null.

The published audit runs ten noise rounds as part of every report, which bounds
the round-level false-alarm probability only loosely. This script runs the same
null — every relative and same-winding point jittered by an independent
1-voxel Gaussian, the audit rerun, and any newly implicated collection counted
as an alarm — for many more rounds, with the configuration exactly as frozen
and a seed that played no part in selecting it.

    python scripts/frozen_null.py data/paris4 results/frozen_null.json --rounds 200
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from windaudit.calibrate import calibrate_from_absolute_anchors  # noqa: E402
from windaudit.controls import JITTER_SIGMA, _jitter, audit_once, implicated_nodes  # noqa: E402
from windaudit.frame import Umbilicus, build_nodes, node_tag  # noqa: E402
from windaudit.pcl import load_point_collections  # noqa: E402

SEED = 20260921  # fixed after the configuration was frozen; not used in selection


def clopper_pearson(k, n, level=0.95):
    lo = 0.0 if k == 0 else float(beta.ppf((1 - level) / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - (1 - level) / 2, k + 1, n - k))
    return lo, hi


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs")
    parser.add_argument("out")
    parser.add_argument("--rounds", type=int, default=200)
    parser.add_argument("--spiral-sense", default="CW")
    args = parser.parse_args(argv)

    d = Path(args.inputs)
    relative = load_point_collections(str(d / "relative_windings.json"))
    same = load_point_collections(str(d / "same_windings.json"))
    absolute = load_point_collections(str(d / "abs_winding.json"))
    umbilicus = Umbilicus.load(str(d / "umbilicus.json"))
    nodes = build_nodes(relative, same, absolute, umbilicus, spiral_sense=args.spiral_sense)
    cfg = calibrate_from_absolute_anchors(nodes, spiral_sense=args.spiral_sense)

    _, base_match, base_repair = audit_once(relative, same, absolute, umbilicus, cfg)
    base = implicated_nodes(base_match, base_repair)

    rng = np.random.default_rng(SEED)
    rounds, started = [], time.perf_counter()
    for i in range(args.rounds):
        rel_j = _jitter(relative, rng, JITTER_SIGMA)
        same_j = _jitter(same, rng, JITTER_SIGMA)
        _, match_j, repair_j = audit_once(rel_j, same_j, absolute, umbilicus, cfg)
        new = sorted(node_tag(k) for k in implicated_nodes(match_j, repair_j) - base)
        rounds.append(new)
        if (i + 1) % 25 == 0:
            print(f"{i + 1}/{args.rounds} rounds, {sum(bool(r) for r in rounds)} with an alarm", flush=True)

    k, n = sum(bool(r) for r in rounds), len(rounds)
    lo, hi = clopper_pearson(k, n)
    counts = {}
    for r in rounds:
        for tag in r:
            counts[tag] = counts.get(tag, 0) + 1
    report = {
        "seed": SEED,
        "sigma_voxels": JITTER_SIGMA,
        "config_sha256": cfg.config_hash,
        "rounds": n,
        "rounds_with_alarm": k,
        "alarm_rate": k / n,
        "clopper_pearson_95": [lo, hi],
        "mean_new_implicated_per_round": float(np.mean([len(r) for r in rounds])),
        "collections_ever_implicated": dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))),
        "per_round_new_implicated": rounds,
        "runtime_seconds": round(time.perf_counter() - started, 1),
    }
    Path(args.out).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"{k}/{n} rounds with an alarm; 95% interval {lo:.3f}-{hi:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
