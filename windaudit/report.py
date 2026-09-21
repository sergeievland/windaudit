"""Full audit pipeline: load, calibrate, audit, repair, localize, export."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import time
from dataclasses import replace
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy

from . import __version__
from .audit import order_audit, same_collection_step_audit, shuffled_order_control
from .calibrate import FrozenConfig, calibrate_from_absolute_anchors
from .controls import implicated_nodes, planted_defect_validation
from .frame import (
    KIND_ABSOLUTE,
    KIND_RELATIVE,
    KIND_SAME,
    NodeKey,
    Umbilicus,
    build_nodes,
    node_tag,
)
from .graph import WindingGraph, constant_label_nodes
from .localize import build_review_queue
from .matching import MatchResult, build_links
from .pcl import Collection, Point, load_point_collections, save_point_collections
from .repair import (RepairResult, edges_consistent, quarantine_membership,
                     solve_min_edge_cut, solve_min_quarantine)

SWEEP_R_SAME_FACTORS = (0.25, 0.35, 0.45)
SWEEP_DTHETA_DEG = (3.0, 4.0, 6.0)
SWEEP_DZ = (30.0, 40.0, 60.0)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sensitivity_sweep(nodes, cfg: FrozenConfig) -> Tuple[Dict[NodeKey, float], List[dict]]:
    """Rerun matching and repair over the 27-point protocol grid; per
    collection, the fraction of settings that implicate it."""
    counts: Dict[NodeKey, int] = {}
    rows: List[dict] = []
    constant = constant_label_nodes(nodes)
    for f in SWEEP_R_SAME_FACTORS:
        for dth in SWEEP_DTHETA_DEG:
            for dz in SWEEP_DZ:
                variant = replace(
                    cfg,
                    r_same=f * cfg.wrap_spacing_median,
                    sector_dtheta_rad=math.radians(dth),
                    sector_dz=dz,
                )
                match = build_links(nodes, variant)
                repair = solve_min_quarantine(match.edges)
                stats = WindingGraph(match, constant).stats()
                for key in implicated_nodes(match, repair):
                    counts[key] = counts.get(key, 0) + 1
                rows.append({
                    "r_same_factor": f,
                    "sector_dtheta_deg": dth,
                    "sector_dz": dz,
                    "n_links": len(match.links),
                    "n_edges": len(match.edges),
                    "n_conflict_pairs": len(match.conflicts),
                    "n_effective_cycles": stats["n_effective_cycles"],
                    "n_effective_nonzero_holonomy": stats["n_effective_nonzero_holonomy"],
                    "quarantine": sorted(node_tag(k) for k in repair.removed),
                    "quarantine_optimal": repair.optimal,
                })
    n = len(rows)
    return {k: v / n for k, v in counts.items()}, rows


def certified_links(match: MatchResult, repair: RepairResult,
                    out_path: Optional[str]) -> dict:
    """Links in the consistent remainder, exported as two-point same-winding
    collections in the upstream point-collection format."""
    if not repair.optimal:
        raise RuntimeError("refusing link export: repair was not proven optimal")
    removed = set(repair.removed)
    cut = set(repair.cut_edges)
    remainder = {pair: value for pair, value in match.edges.items()
                 if not (set(pair) & removed) and pair not in cut}
    if not edges_consistent(remainder):
        raise RuntimeError("refusing link export: retained graph is inconsistent")
    kept = []
    for pair in sorted(match.edges):
        if pair[0] in removed or pair[1] in removed or pair in cut:
            continue
        kept.extend(match.edge_links[pair])
    if out_path is not None:
        cols: Dict[int, Collection] = {}
        for i, link in enumerate(kept, start=1):
            cols[i] = Collection(
                id=i,
                name=f"windaudit_{node_tag(link.a)}#{link.pid_a}_{node_tag(link.b)}#{link.pid_b}",
                points={
                    1: Point(id=1, xyz=list(link.xyz_a), wind_a=None),
                    2: Point(id=2, xyz=list(link.xyz_b), wind_a=None),
                },
                metadata={"windaudit_version": __version__,
                          "evidence_type": "geometry_inferred_graph_consistent",
                          "ct_verified": False},
            )
        save_point_collections(out_path, cols)
    linked = {k for pair in match.edges for k in pair
              if pair[0] not in removed and pair[1] not in removed and pair not in cut}
    return {"n_links": len(kept), "n_collections": len(linked),
            "path": os.path.basename(out_path) if out_path else None}


def run_audit(
    relative_path: str,
    same_path: str,
    absolute_path: str,
    umbilicus_path: str,
    out_dir: Optional[str] = None,
    spiral_sense: str = "CW",
    seed: int = 20260917,
    n_shuffle_rounds: int = 20,
    n_null_rounds: int = 10,
    run_planted: bool = True,
    run_sweep: bool = True,
) -> dict:
    t0 = time.time()
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
    inputs = {"relative": relative_path, "same": same_path,
              "absolute": absolute_path, "umbilicus": umbilicus_path}

    relative = load_point_collections(relative_path)
    same = load_point_collections(same_path)
    absolute = load_point_collections(absolute_path)
    umbilicus = Umbilicus.load(umbilicus_path)

    nodes = build_nodes(relative, same, absolute, umbilicus, spiral_sense=spiral_sense)
    cfg = calibrate_from_absolute_anchors(nodes, spiral_sense=spiral_sense)
    constant = constant_label_nodes(nodes)

    def corpus_stats(cols: Dict[int, Collection]) -> dict:
        return {
            "n_collections": len(cols),
            "n_nonempty_collections": sum(1 for c in cols.values() if c.points),
            "n_points": sum(len(c.points) for c in cols.values()),
        }

    corpus = {
        "relative": corpus_stats(relative),
        "same": corpus_stats(same),
        "absolute": corpus_stats(absolute),
        "n_nodes": len(nodes),
        "n_seam_unsafe_nodes": sum(1 for n in nodes.values() if not n.seam_safe),
        "n_nodes_crossing_cut": sum(
            1 for n in nodes.values() if any(q.crossings for q in n.points)
        ),
    }

    order = order_audit(nodes, cfg)
    steps = same_collection_step_audit(nodes, cfg)
    shuffled = shuffled_order_control(nodes, cfg, n_shuffle_rounds, seed)

    match = build_links(nodes, cfg)
    graph = WindingGraph(match, constant)
    stats = graph.stats()
    cycles = graph.cycles()
    failing = [c for c in cycles if c.holonomy != 0]
    quarantine = solve_min_quarantine(match.edges)
    edge_cut = solve_min_edge_cut(match.edges)

    robustness: Dict[NodeKey, float] = {}
    sweep_rows: List[dict] = []
    if run_sweep:
        robustness, sweep_rows = sensitivity_sweep(nodes, cfg)

    queue = build_review_queue(nodes, match, quarantine, robustness, failing,
                               steps.violations)

    linked_labeled = {
        k for e in match.edges for k in e if k[0] in (KIND_RELATIVE, KIND_ABSOLUTE)
    }
    on_effective_cycle = set()
    for c in cycles:
        if c.effective:
            on_effective_cycle.update(c.tree_path)
    coverage = {
        "relative_chains": sum(1 for k in nodes if k[0] == KIND_RELATIVE),
        "relative_chains_linked": sum(1 for k in linked_labeled if k[0] == KIND_RELATIVE),
        "same_traces": sum(1 for k in nodes if k[0] == KIND_SAME),
        "same_traces_linked": len({k for e in match.edges for k in e if k[0] == KIND_SAME}),
        "collections_on_effective_cycles": len(on_effective_cycle),
    }

    cert = certified_links(
        match, edge_cut,
        os.path.join(out_dir, "certified_links.json") if out_dir else None,
    )

    planted = None
    if run_planted:
        planted = planted_defect_validation(
            relative, same, absolute, umbilicus, cfg, seed=seed,
            n_null_rounds=n_null_rounds,
        )

    report = {
        "tool": {"name": "windaudit", "version": __version__},
        "environment": {"python": platform.python_version(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
        "inputs": {k: os.path.basename(v) for k, v in inputs.items()},
        "input_sha256": {k: sha256_file(v) for k, v in inputs.items()},
        "corpus": corpus,
        "calibration": cfg.as_dict(),
        "config_sha256": cfg.config_hash,
        "tier1": {
            "order_audit": {
                "n_pairs_tested": order.n_pairs_tested,
                "n_violations": order.n_violations,
                "flagged_collections": sorted(node_tag(k) for k in order.flagged_nodes),
            },
            "shuffled_control": shuffled,
            "same_winding_step_audit": {
                "n_collections_audited": steps.n_collections_audited,
                "n_steps_tested": steps.n_steps_tested,
                "n_violations": steps.n_violations,
                "flagged_collections": sorted(node_tag(k) for k in steps.flagged_nodes),
                "violations": [
                    {"collection": node_tag(v.node), "pid_a": v.pid_a, "pid_b": v.pid_b,
                     "xyz_a": v.xyz_a, "xyz_b": v.xyz_b, "residual": round(v.residual, 2)}
                    for v in steps.violations
                ],
            },
        },
        "matching": {
            "n_links": len(match.links),
            "n_edges": len(match.edges),
            "n_conflict_pairs": len(match.conflicts),
            "n_links_across_cut": match.n_cut_links,
            "n_unsupported_links": match.n_unsupported_links,
            "n_cells": match.n_cells,
            "n_vetoed_cells": match.n_vetoed_cells,
            "n_points_considered": match.n_points_considered,
            "n_points_with_slope": match.n_points_with_slope,
            "conflicts": [
                {"a": node_tag(c.a), "b": node_tag(c.b), "offsets": c.offsets,
                 "n_links": len(c.links)}
                for c in match.conflicts
            ],
        },
        "graph": stats,
        "coverage": coverage,
        "failing_cycles": [
            {"closing_edge": [node_tag(c.u), node_tag(c.v)], "holonomy": c.holonomy,
             "closing_support": c.support, "cycle": [node_tag(k) for k in c.tree_path]}
            for c in failing
        ],
        "repair": {
            "min_quarantine": {
                "removed": sorted(node_tag(k) for k in quarantine.removed),
                "optimal": quarantine.optimal,
                "status": quarantine.status,
                "links_removed": quarantine.support_removed,
                "membership": {
                    node_tag(k): v for k, v in sorted(
                        quarantine_membership(
                            match.edges,
                            sorted(set(quarantine.removed)
                                   | {c.a for c in match.conflicts}
                                   | {c.b for c in match.conflicts}),
                        ).items()
                    )
                },
            },
            "min_edge_cut": {
                "cut_edges": [[node_tag(a), node_tag(b)] for a, b in edge_cut.cut_edges],
                "optimal": edge_cut.optimal,
                "status": edge_cut.status,
                "links_removed": edge_cut.support_removed,
            },
            "links_retained_per_link_removed": (
                (len(match.links) - edge_cut.support_removed) / edge_cut.support_removed
                if edge_cut.support_removed else None
            ),
        },
        "certified_links": cert,
        "sensitivity_sweep": sweep_rows,
        "robustness": {node_tag(k): round(v, 4) for k, v in sorted(robustness.items())},
        "review_queue": queue,
        "planted_defect_validation": planted,
        "seed": seed,
        "runtime_seconds": round(time.time() - t0, 1),
    }

    if out_dir is not None:
        with open(os.path.join(out_dir, "audit_report.json"), "w", encoding="utf-8") as f:
            json.dump(report, f, indent=1)
        with open(os.path.join(out_dir, "review_queue.json"), "w", encoding="utf-8") as f:
            json.dump(queue, f, indent=1)
    return report
