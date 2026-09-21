"""Planted-defect validation and a measurement-noise null control.

Two defect models mirror the two real failure modes of winding annotation:

* skipped wrap - a relative ladder misses (or double-counts) one wrap: every
  ``wind_a`` after an interior boundary shifts by one;
* sheet switch - a same-winding trace drifts onto the adjacent wrap: every
  point after an interior boundary moves one wrap spacing radially.

A defect is only visible to any annotation-level audit if independent
evidence spans the planted boundary. The validation therefore enumerates, on
the unmodified corpus, every boundary where it does ("sites"), plants both
defect directions at every site, reruns the whole pipeline, and reports:

* coverage - how many sites the corpus offers, out of all boundaries;
* alarm - the audit implicates a collection it did not implicate before;
* localized - the planted collection itself is newly implicated;
* shortlisted - the planted collection is newly on a failing cycle.

The null control jitters every point by isotropic measurement noise and
reports how often the audit raises an alarm when nothing is wrong.
"""

from __future__ import annotations

import copy
import math
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .calibrate import FrozenConfig
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
from .matching import MatchResult, build_links
from .pcl import Collection
from .repair import quarantine_membership, solve_min_quarantine

Site = Tuple[int, int, List[int]]
JITTER_SIGMA = 1.0


def implicated_nodes(match: MatchResult, repair) -> Set[NodeKey]:
    """Collections the audit names: conflicting pairs plus the quarantine."""
    out: Set[NodeKey] = set()
    for c in match.conflicts:
        out.add(c.a)
        out.add(c.b)
    out.update(repair.removed)
    return out


def failing_cycle_nodes(match: MatchResult, nodes) -> Set[NodeKey]:
    """Collections on a nonzero-holonomy cycle or in a conflicting pair."""
    out: Set[NodeKey] = set()
    for c in match.conflicts:
        out.add(c.a)
        out.add(c.b)
    graph = WindingGraph(match, constant_label_nodes(nodes))
    for cyc in graph.cycles():
        if cyc.holonomy != 0:
            out.update(cyc.tree_path)
    return out


def audit_once(relative, same, absolute, umbilicus, cfg: FrozenConfig):
    nodes = build_nodes(relative, same, absolute, umbilicus,
                        spiral_sense=cfg.spiral_sense,
                        max_step_rad=cfg.max_chain_step_rad)
    match = build_links(nodes, cfg)
    repair = solve_min_quarantine(match.edges)
    return nodes, match, repair


def _attachments(match: MatchResult, key: NodeKey) -> List[Tuple[int, NodeKey]]:
    out = []
    for link in match.links:
        if link.a == key:
            out.append((link.pid_a, link.b))
        elif link.b == key:
            out.append((link.pid_b, link.a))
    return out


def is_reachable(match: MatchResult, key: NodeKey, pids: List[int], boundary: int,
                 min_support: int) -> bool:
    """True when evidence independent of the planted collection connects the
    two sides of the boundary, each side carrying at least ``min_support``
    links to its partner (the matcher never accepts less). For a same-winding
    collection the connecting evidence must include a label-bearing
    collection: a path made only of same-winding traces asserts nothing but
    zero offsets."""
    left, right = set(pids[:boundary]), set(pids[boundary:])
    left_count: Dict[NodeKey, int] = {}
    right_count: Dict[NodeKey, int] = {}
    for pid, partner in _attachments(match, key):
        if pid in left:
            left_count[partner] = left_count.get(partner, 0) + 1
        elif pid in right:
            right_count[partner] = right_count.get(partner, 0) + 1
    left_p = {p for p, c in left_count.items() if c >= min_support}
    right_p = {p for p, c in right_count.items() if c >= min_support}
    if not left_p or not right_p:
        return False
    labeled = key[0] in (KIND_RELATIVE, KIND_ABSOLUTE)
    if labeled and left_p & right_p:
        return True

    parent: Dict[NodeKey, NodeKey] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in match.edges:
        if key not in (a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    shared = {find(p) for p in left_p} & {find(p) for p in right_p}
    if not shared:
        return False
    if labeled:
        return True
    label_roots = {
        find(k) for e in match.edges for k in e
        if k != key and k[0] in (KIND_RELATIVE, KIND_ABSOLUTE)
    }
    return bool(shared & label_roots)


def _ordered_pids(col: Collection, kind: str) -> List[int]:
    if kind == KIND_RELATIVE:
        return [p.id for p in col.ordered_points() if p.has_winding]
    return [p.id for p in col.ordered_points()]


def enumerate_sites(match: MatchResult, collections: Dict[int, Collection],
                    kind: str, min_support: int,
                    min_points: int = 4) -> Tuple[List[Site], int]:
    """(reachable sites, total interior boundaries) for one defect family."""
    sites: List[Site] = []
    total = 0
    for cid in sorted(collections):
        pids = _ordered_pids(collections[cid], kind)
        if len(pids) < min_points:
            continue
        for boundary in range(1, len(pids)):
            total += 1
            if is_reachable(match, (kind, cid), pids, boundary, min_support):
                sites.append((cid, boundary, pids))
    return sites, total


def plant_step_defect(relative: Dict[int, Collection], site: Site, shift: int):
    cid, boundary, pids = site
    mutated = dict(relative)
    col = copy.deepcopy(relative[cid])
    for pid in pids[boundary:]:
        col.points[pid].wind_a += shift
    mutated[cid] = col
    return mutated


def plant_sheet_switch(same: Dict[int, Collection], umbilicus: Umbilicus,
                       spacing: float, site: Site, direction: float):
    cid, boundary, pids = site
    mutated = dict(same)
    col = copy.deepcopy(same[cid])
    for pid in pids[boundary:]:
        p = col.points[pid]
        cx, cy = umbilicus.axis_xy(p.xyz[2])
        dx, dy = p.xyz[0] - cx, p.xyz[1] - cy
        r = math.hypot(dx, dy)
        if r < 1e-9:
            continue
        scale = (r + direction * spacing) / r
        p.xyz[0] = cx + dx * scale
        p.xyz[1] = cy + dy * scale
    mutated[cid] = col
    return mutated


@dataclass
class PlantedTrial:
    defect: str
    collection: str
    boundary: int
    direction: int
    reachable_after: bool
    alarm: bool
    localized: bool
    shortlisted: bool
    n_new_implicated: int
    n_new_shortlisted: int
    membership: str
    certainly_named: bool


def _jitter(cols: Dict[int, Collection], rng: np.random.Generator, sigma: float):
    out = {}
    for cid in sorted(cols):
        col = copy.deepcopy(cols[cid])
        for pid in sorted(col.points):
            p = col.points[pid]
            p.xyz[0] += float(rng.normal(0.0, sigma))
            p.xyz[1] += float(rng.normal(0.0, sigma))
            p.xyz[2] += float(rng.normal(0.0, sigma))
        out[cid] = col
    return out


def planted_defect_validation(
    relative: Dict[int, Collection],
    same: Dict[int, Collection],
    absolute: Dict[int, Collection],
    umbilicus: Umbilicus,
    cfg: FrozenConfig,
    seed: int,
    n_null_rounds: int = 10,
    max_sites_per_type: Optional[int] = None,
) -> dict:
    rng = np.random.default_rng(seed)
    base_nodes, base_match, base_repair = audit_once(relative, same, absolute, umbilicus, cfg)
    base_implicated = implicated_nodes(base_match, base_repair)
    base_failing = failing_cycle_nodes(base_match, base_nodes)

    min_support = cfg.min_edge_support
    step_sites, step_total = enumerate_sites(base_match, relative, KIND_RELATIVE, min_support)
    switch_sites, switch_total = enumerate_sites(base_match, same, KIND_SAME, min_support)
    if max_sites_per_type is not None:
        for sites in (step_sites, switch_sites):
            if len(sites) > max_sites_per_type:
                keep = sorted(rng.choice(len(sites), max_sites_per_type, replace=False))
                sites[:] = [sites[i] for i in keep]

    trials: List[PlantedTrial] = []
    for defect, sites in (("skipped_wrap", step_sites), ("sheet_switch", switch_sites)):
        for site in sites:
            cid, boundary, pids = site
            for direction in (-1, 1):
                if defect == "skipped_wrap":
                    key: NodeKey = (KIND_RELATIVE, cid)
                    rel_m = plant_step_defect(relative, site, direction)
                    nodes, match, repair = audit_once(rel_m, same, absolute, umbilicus, cfg)
                else:
                    key = (KIND_SAME, cid)
                    same_m = plant_sheet_switch(same, umbilicus, cfg.wrap_spacing_median,
                                                site, float(direction))
                    nodes, match, repair = audit_once(relative, same_m, absolute, umbilicus, cfg)
                new_imp = implicated_nodes(match, repair) - base_implicated
                new_fail = failing_cycle_nodes(match, nodes) - base_failing
                status = quarantine_membership(match.edges, [key])[key]
                trials.append(PlantedTrial(
                    defect=defect,
                    collection=node_tag(key),
                    boundary=boundary,
                    direction=direction,
                    reachable_after=is_reachable(match, key, pids, boundary, min_support),
                    alarm=bool(new_imp),
                    localized=key in new_imp,
                    shortlisted=key in new_fail,
                    n_new_implicated=len(new_imp),
                    n_new_shortlisted=len(new_fail),
                    membership=status,
                    certainly_named=(status == "certain" and key not in base_implicated),
                ))

    def summary(name: str, n_sites: int, n_total: int, chains_total: int) -> dict:
        sub = [t for t in trials if t.defect == name]
        n = len(sub)

        def rate(flag):
            return sum(1 for t in sub if getattr(t, flag)) / n if n else None

        listed = [t for t in sub if t.shortlisted]
        return {
            "boundaries_total": n_total,
            "reachable_sites": n_sites,
            "site_coverage": n_sites / n_total if n_total else None,
            "collections_with_sites": len({t.collection for t in sub}),
            "collections_total": chains_total,
            "n_trials": n,
            "alarm_rate": rate("alarm"),
            "localization_rate": rate("localized"),
            "certain_localization_rate": rate("certainly_named"),
            "membership_counts": {
                value: sum(1 for t in sub if t.membership == value)
                for value in ("certain", "possible", "excluded", "consistent", "unknown")
            },
            "shortlist_rate": rate("shortlisted"),
            "mean_shortlist_size": (
                float(np.mean([t.n_new_shortlisted for t in listed])) if listed else None
            ),
            "still_reachable_after_plant": rate("reachable_after"),
        }

    null_alarms = []
    for _ in range(n_null_rounds):
        rel_j = _jitter(relative, rng, JITTER_SIGMA)
        same_j = _jitter(same, rng, JITTER_SIGMA)
        _, match_j, repair_j = audit_once(rel_j, same_j, absolute, umbilicus, cfg)
        null_alarms.append(len(implicated_nodes(match_j, repair_j) - base_implicated))

    def n_chains(cols, kind):
        return sum(1 for c in cols.values() if len(_ordered_pids(c, kind)) >= 4)

    return {
        "seed": seed,
        "baseline_implicated": sorted(node_tag(n) for n in base_implicated),
        "skipped_wrap": summary("skipped_wrap", len(step_sites), step_total,
                                n_chains(relative, KIND_RELATIVE)),
        "sheet_switch": summary("sheet_switch", len(switch_sites), switch_total,
                                n_chains(same, KIND_SAME)),
        "null_jitter": {
            "sigma_voxels": JITTER_SIGMA,
            "n_rounds": n_null_rounds,
            "rounds_with_alarm": sum(1 for a in null_alarms if a),
            "mean_new_implicated": float(np.mean(null_alarms)) if null_alarms else None,
        },
        "trials": [asdict(t) for t in trials],
    }
