"""Turn audit findings into a coordinate-addressed review queue."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional

from .frame import Node, NodeKey, node_tag
from .graph import CycleFinding
from .matching import EdgeConflict, MatchResult
from .repair import RepairResult


def _conflict_bracket(node: Node, conflict: EdgeConflict, key: NodeKey) -> Optional[dict]:
    """Group this collection's matched points by implied offset; two adjacent
    offsets bracket the stretch where the disagreement arises."""
    by_pid = {q.pid: q for q in node.points}
    per_offset: Dict[int, List[tuple]] = defaultdict(list)
    for link in conflict.links:
        if link.a == key:
            pid, offset = link.pid_a, link.offset
        elif link.b == key:
            pid, offset = link.pid_b, -link.offset
        else:
            continue
        q = by_pid.get(pid)
        if q is not None:
            per_offset[offset].append((q.r, list(q.xyz)))
    offsets = sorted(per_offset)
    for oa, ob in zip(offsets, offsets[1:]):
        lo = max(per_offset[oa], key=lambda t: t[0])
        hi = min(per_offset[ob], key=lambda t: t[0])
        return {
            "partner": node_tag(conflict.b if conflict.a == key else conflict.a),
            "offsets": [oa, ob],
            "xyz_at_first_offset": lo[1],
            "xyz_at_second_offset": hi[1],
        }
    return None


def _cycle_witness(match: MatchResult, finding: CycleFinding) -> List[dict]:
    """The links that realise the closing edge of a failing cycle."""
    key = (finding.u, finding.v)
    return [
        {"xyz_a": l.xyz_a, "xyz_b": l.xyz_b, "radial_gap": round(l.gap, 2)}
        for l in match.edge_links.get(key, [])
    ]


def build_review_queue(
    nodes: Dict[NodeKey, Node],
    match: MatchResult,
    repair: RepairResult,
    robustness: Dict[NodeKey, float],
    failing_cycles: List[CycleFinding],
    step_violations=(),
) -> List[dict]:
    evidence: Dict[NodeKey, List[str]] = defaultdict(list)
    locations: Dict[NodeKey, List[dict]] = defaultdict(list)
    for v in step_violations:
        evidence[v.node].append("wrap-scale radial step against the trace's local trend")
        locations[v.node].append({
            "step_from_xyz": v.xyz_a,
            "step_to_xyz": v.xyz_b,
            "residual_voxels": round(v.residual, 2),
        })
    for c in match.conflicts:
        for key, other in ((c.a, c.b), (c.b, c.a)):
            evidence[key].append(f"links to {node_tag(other)} imply offsets {c.offsets}")
            if key in nodes:
                bracket = _conflict_bracket(nodes[key], c, key)
                if bracket is not None:
                    locations[key].append(bracket)
    for cyc in failing_cycles:
        members = ", ".join(node_tag(k) for k in cyc.tree_path)
        for key in cyc.tree_path:
            evidence[key].append(f"on a cycle with holonomy {cyc.holonomy:+d} ({members})")
        witness = _cycle_witness(match, cyc)
        for key in (cyc.u, cyc.v):
            locations[key].append({"closing_edge_links": witness})
    removed = set(repair.removed)
    for key in removed:
        evidence[key].append("member of the exact minimum quarantine")

    queue = []
    for key in sorted(evidence):
        node = nodes.get(key)
        if node is None:
            continue
        queue.append({
            "collection": node_tag(key),
            "name": node.name,
            "n_points": node.n,
            "evidence": sorted(set(evidence[key])),
            "robustness": round(robustness.get(key, 0.0), 4),
            "in_quarantine": key in removed,
            "locations": locations.get(key, []),
            "example_xyz": list(node.points[0].xyz),
        })
    queue.sort(key=lambda e: (-e["robustness"], not e["in_quarantine"], e["collection"]))
    return queue
