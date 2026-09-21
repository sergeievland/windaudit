"""Tier-1 audits: consistency inside single collections.

These tests rest on one physical fact - papyrus wraps do not cross - so at a
common z inside a small sector, radial order must agree with winding order.
Points are projected to a common z along their own wrap's radial slope (the
same projection the matcher uses), and label differences are taken in one
branch frame. No cross-collection matching and no model are involved.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from .calibrate import FrozenConfig
from .frame import (
    KIND_ABSOLUTE,
    KIND_RELATIVE,
    KIND_SAME,
    Node,
    NodeKey,
    label_diff,
    wrap_angle,
)
from .matching import point_slopes, projected_gap


@dataclass
class OrderViolation:
    node: NodeKey
    pid_a: int
    pid_b: int
    xyz_a: List[float]
    xyz_b: List[float]
    d_label: int
    gap: float
    kind: str


@dataclass
class OrderAuditResult:
    n_pairs_tested: int = 0
    violations: List[OrderViolation] = field(default_factory=list)
    flagged_nodes: Dict[NodeKey, int] = field(default_factory=dict)

    @property
    def n_violations(self) -> int:
        return len(self.violations)

    @property
    def violation_rate(self) -> float:
        return self.n_violations / self.n_pairs_tested if self.n_pairs_tested else 0.0


def order_audit(
    nodes: Dict[NodeKey, Node],
    cfg: FrozenConfig,
    kinds: Tuple[str, ...] = (KIND_RELATIVE, KIND_ABSOLUTE),
) -> OrderAuditResult:
    """Within one collection, same-sector pairs must order the same way by
    projected radius and by winding label ("sign" violations); a same-label
    pair further apart than the crowding floor is a "split"."""
    sense = cfg.sense
    result = OrderAuditResult()
    for key in sorted(nodes):
        node = nodes[key]
        if node.kind not in kinds or not node.seam_safe:
            continue
        slopes = point_slopes(node, sense)
        pts = node.points
        count = 0
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                a, b = pts[i], pts[j]
                if abs(a.z - b.z) > cfg.sector_dz:
                    continue
                if abs(wrap_angle(a.theta - b.theta)) > cfg.sector_dtheta_rad:
                    continue
                gap = projected_gap(a, slopes[i], b, slopes[j])
                if gap is None:
                    continue
                result.n_pairs_tested += 1
                d_label = label_diff(a, b, sense)
                kind = None
                if d_label == 0:
                    if abs(gap) > cfg.crowd_min_spacing:
                        kind = "split"
                elif d_label * gap < 0 and abs(gap) > cfg.r_same:
                    kind = "sign"
                if kind is not None:
                    result.violations.append(
                        OrderViolation(key, a.pid, b.pid, a.xyz, b.xyz, d_label, gap, kind)
                    )
                    count += 1
        if count:
            result.flagged_nodes[key] = count
    return result


@dataclass
class StepViolation:
    node: NodeKey
    pid_a: int
    pid_b: int
    xyz_a: List[float]
    xyz_b: List[float]
    residual: float


@dataclass
class StepAuditResult:
    n_steps_tested: int = 0
    n_collections_audited: int = 0
    violations: List[StepViolation] = field(default_factory=list)
    flagged_nodes: Dict[NodeKey, int] = field(default_factory=dict)

    @property
    def n_violations(self) -> int:
        return len(self.violations)


STEP_TREND_HALF_WINDOW = 3


def same_collection_step_audit(
    nodes: Dict[NodeKey, Node], cfg: FrozenConfig
) -> StepAuditResult:
    """Risk figure for the same-winding hypothesis. A same-winding trace that
    takes a wrap-scale radial step against its own local trend asserts one
    wrap while its geometry jumps to another. Each consecutive step's radial
    move is compared with the median slope of up to three neighbouring steps
    on either side; a residual above the crowding floor is flagged."""
    result = StepAuditResult()
    for key in sorted(nodes):
        node = nodes[key]
        if node.kind != KIND_SAME:
            continue
        result.n_collections_audited += 1
        pts = node.points
        if len(pts) < 4:
            continue
        steps = [(b.r - a.r, b.z - a.z) for a, b in zip(pts, pts[1:])]
        count = 0
        for i, (a, b) in enumerate(zip(pts, pts[1:])):
            dr, dz = steps[i]
            if abs(dz) < 1.0:
                continue
            lo, hi = max(0, i - STEP_TREND_HALF_WINDOW), i + STEP_TREND_HALF_WINDOW
            neighbours = [
                s_dr / s_dz
                for j, (s_dr, s_dz) in enumerate(steps)
                if j != i and lo <= j <= hi and abs(s_dz) >= 1.0
            ]
            if not neighbours:
                continue
            residual = dr - float(np.median(neighbours)) * dz
            result.n_steps_tested += 1
            if abs(residual) > cfg.crowd_min_spacing:
                result.violations.append(
                    StepViolation(key, a.pid, b.pid, a.xyz, b.xyz, residual)
                )
                count += 1
        if count:
            result.flagged_nodes[key] = count
    return result


def shuffled_order_control(
    nodes: Dict[NodeKey, Node],
    cfg: FrozenConfig,
    n_rounds: int,
    seed: int,
) -> dict:
    """Permute winding labels within each relative chain and repeat the order
    audit on the same pairs. The gap between the real and shuffled violation
    rates is the audit's power on this corpus."""
    rng = np.random.default_rng(seed)
    rel_keys = [k for k in sorted(nodes) if nodes[k].kind == KIND_RELATIVE]
    rates = []
    for _ in range(n_rounds):
        shuffled: Dict[NodeKey, Node] = {}
        for key in rel_keys:
            node = nodes[key]
            labels = np.array([q.label for q in node.points])
            rng.shuffle(labels)
            clone = Node(kind=node.kind, cid=node.cid, name=node.name,
                         source_cid=node.source_cid, points=[],
                         seam_safe=node.seam_safe)
            for q, lab in zip(node.points, labels):
                q2 = type(q)(**q.__dict__)
                q2.label = int(lab)
                clone.points.append(q2)
            shuffled[key] = clone
        rates.append(order_audit(shuffled, cfg, kinds=(KIND_RELATIVE,)).violation_rate)
    real = order_audit({k: nodes[k] for k in rel_keys}, cfg, kinds=(KIND_RELATIVE,))
    return {
        "n_pairs": real.n_pairs_tested,
        "real_violation_rate": real.violation_rate,
        "shuffled_mean_rate": float(np.mean(rates)) if rates else None,
        "shuffled_min_rate": float(np.min(rates)) if rates else None,
        "shuffled_max_rate": float(np.max(rates)) if rates else None,
        "n_rounds": n_rounds,
        "seed": seed,
    }
