"""Same-sheet link extraction between collections, tilt- and cut-aware.

Two points from different collections are linked as same-wrap evidence when,
inside one sector, their radii agree within ``r_same`` after projecting each
point along its own wrap's local radial slope to a common reference z. Wraps
tilt, and on real scrolls the radial drift across a sector's height is
comparable to the wrap spacing, so a tilt-blind comparison manufactures links
to the neighbouring wrap. Points whose wrap provides no local slope
(single-slice ladders) are never projected; the partner is projected to their
z instead, and slope-less pairs are compared only within ``STRICT_DZ``.

Sectors where any single collection places two radially adjacent windings
closer than ``crowd_min_spacing`` are vetoed: where wraps crowd below the
matcher's resolution, no links are emitted.

Each link asserts one integer equation between node gauges,
``gauge_b - gauge_a = label_a - label_b + sense * cut_crossing(a, b)``.
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .calibrate import FrozenConfig
from .frame import Node, NodeKey, NodePoint, TWO_PI, label_diff, wrap_angle

SLOPE_WINDOW_DZ = 60.0
SLOPE_WINDOW_DTHETA_RAD = math.radians(6.0)
SLOPE_MIN_ZSPAN = 10.0
STRICT_DZ = 10.0

EdgeKey = Tuple[NodeKey, NodeKey]


@dataclass
class Link:
    a: NodeKey
    b: NodeKey
    pid_a: int
    pid_b: int
    xyz_a: List[float]
    xyz_b: List[float]
    offset: int
    gap: float


@dataclass
class EdgeConflict:
    a: NodeKey
    b: NodeKey
    offsets: List[int]
    links: List[Link]


@dataclass
class MatchResult:
    links: List[Link]
    edges: Dict[EdgeKey, Tuple[int, int]]
    edge_links: Dict[EdgeKey, List[Link]]
    conflicts: List[EdgeConflict]
    n_cells: int
    n_vetoed_cells: int
    n_points_considered: int
    n_points_with_slope: int
    n_cut_links: int
    n_unsupported_links: int


def point_slopes(node: Node, sense: int) -> List[Optional[float]]:
    """Per-point local radial slope dr/dz, fitted on the point's own wrap
    (same label after cut correction) inside a small (z, theta) window; None
    when that wrap is too flat in z there to constrain a slope."""
    pts = node.points
    slopes: List[Optional[float]] = []
    for q in pts:
        zs, rs = [], []
        for p in pts:
            if abs(p.z - q.z) > SLOPE_WINDOW_DZ:
                continue
            if abs(wrap_angle(p.theta - q.theta)) > SLOPE_WINDOW_DTHETA_RAD:
                continue
            if label_diff(q, p, sense) != 0:
                continue
            zs.append(p.z)
            rs.append(p.r)
        if len(zs) >= 2 and (max(zs) - min(zs)) >= SLOPE_MIN_ZSPAN:
            z_arr = np.asarray(zs)
            r_arr = np.asarray(rs)
            zc = z_arr - z_arr.mean()
            slopes.append(float(np.dot(zc, r_arr - r_arr.mean()) / np.dot(zc, zc)))
        else:
            slopes.append(None)
    return slopes


def projected_gap(
    qa: NodePoint, sa: Optional[float], qb: NodePoint, sb: Optional[float]
) -> Optional[float]:
    """Radial gap ``r_a - r_b`` at a common z, or None when not comparable."""
    if sa is None and sb is None:
        if abs(qb.z - qa.z) > STRICT_DZ:
            return None
        return qa.r - qb.r
    if sa is None:
        return qa.r - (qb.r + sb * (qa.z - qb.z))
    if sb is None:
        return (qa.r + sa * (qb.z - qa.z)) - qb.r
    z_ref = 0.5 * (qa.z + qb.z)
    return (qa.r + sa * (z_ref - qa.z)) - (qb.r + sb * (z_ref - qb.z))


def _cell_of(point: NodePoint, cfg: FrozenConfig, ntheta: int) -> Tuple[int, int]:
    zi = int(math.floor(point.z / cfg.sector_dz))
    ti = int(point.theta / TWO_PI * ntheta) % ntheta
    return zi, ti


def build_links(nodes: Dict[NodeKey, Node], cfg: FrozenConfig) -> MatchResult:
    """Extract same-sheet links under the frozen configuration. Seam-unsafe
    nodes are excluded."""
    sense = cfg.sense
    ntheta = max(3, int(round(TWO_PI / cfg.sector_dtheta_rad)))
    slope_of: Dict[Tuple[NodeKey, int], Optional[float]] = {}
    cells: Dict[Tuple[int, int], List[Tuple[NodeKey, NodePoint]]] = defaultdict(list)
    n_points = 0
    for key in sorted(nodes):
        node = nodes[key]
        if not node.seam_safe:
            continue
        for q, s in zip(node.points, point_slopes(node, sense)):
            slope_of[(key, q.pid)] = s
            cells[_cell_of(q, cfg, ntheta)].append((key, q))
            n_points += 1

    pair_links: Dict[EdgeKey, List[Link]] = defaultdict(list)
    seen = set()
    n_vetoed = 0
    n_cut_links = 0

    for cell in sorted(cells):
        zi, ti = cell
        group = list(cells[cell])
        for dz, dt in ((0, 1), (1, -1), (1, 0), (1, 1)):
            neighbour = (zi + dz, (ti + dt) % ntheta)
            if neighbour != cell:
                group += cells.get(neighbour, [])

        # Crowding veto, evaluated in the branch of the cell's first point and
        # at the neighbourhood's median z.
        ref = group[0][1]
        z_ref = float(np.median([q.z for _, q in group]))
        by_node: Dict[NodeKey, Dict[int, List[float]]] = defaultdict(lambda: defaultdict(list))
        for key, q in group:
            s = slope_of[(key, q.pid)]
            r_proj = q.r if s is None else q.r + s * (z_ref - q.z)
            by_node[key][label_diff(q, ref, sense) + ref.label].append(r_proj)
        crowded = False
        for by_label in by_node.values():
            labels = sorted(by_label)
            for la, lb in zip(labels, labels[1:]):
                if lb - la == 1:
                    gap = min(by_label[lb]) - max(by_label[la])
                    if abs(gap) < cfg.crowd_min_spacing:
                        crowded = True
                        break
            if crowded:
                break
        if crowded:
            n_vetoed += 1
            continue

        for (ka, qa), (kb, qb) in itertools.combinations(group, 2):
            if ka == kb:
                continue
            if abs(qa.z - qb.z) > cfg.sector_dz:
                continue
            if abs(wrap_angle(qa.theta - qb.theta)) > cfg.sector_dtheta_rad:
                continue
            if ka > kb:
                ka, qa, kb, qb = kb, qb, ka, qa
            pair_id = (ka, qa.pid, kb, qb.pid)
            if pair_id in seen:
                continue
            seen.add(pair_id)
            gap = projected_gap(qa, slope_of[(ka, qa.pid)], qb, slope_of[(kb, qb.pid)])
            if gap is None or abs(gap) > cfg.r_same:
                continue
            link = Link(
                a=ka, b=kb, pid_a=qa.pid, pid_b=qb.pid,
                xyz_a=list(qa.xyz), xyz_b=list(qb.xyz),
                offset=label_diff(qa, qb, sense), gap=gap,
            )
            if abs(qa.theta - qb.theta) > math.pi:
                n_cut_links += 1
            pair_links[(ka, kb)].append(link)

    edges: Dict[EdgeKey, Tuple[int, int]] = {}
    edge_links: Dict[EdgeKey, List[Link]] = {}
    conflicts: List[EdgeConflict] = []
    n_unsupported = 0
    for pair in sorted(pair_links):
        by_offset: Dict[int, List[Link]] = defaultdict(list)
        for link in pair_links[pair]:
            by_offset[link.offset].append(link)
        supported = {o: ls for o, ls in by_offset.items()
                     if len(ls) >= cfg.min_edge_support}
        n_unsupported += sum(len(ls) for o, ls in by_offset.items() if o not in supported)
        if len(supported) == 1:
            (offset, links), = supported.items()
            edges[pair] = (offset, len(links))
            edge_links[pair] = links
        elif len(supported) > 1:
            links = [l for o in sorted(supported) for l in supported[o]]
            conflicts.append(EdgeConflict(pair[0], pair[1], sorted(supported), links))

    kept = [l for pair in sorted(edge_links) for l in edge_links[pair]]
    kept += [l for c in conflicts for l in c.links]
    return MatchResult(
        links=kept,
        edges=edges,
        edge_links=edge_links,
        conflicts=conflicts,
        n_cells=len(cells),
        n_vetoed_cells=n_vetoed,
        n_points_considered=n_points,
        n_points_with_slope=sum(1 for s in slope_of.values() if s is not None),
        n_cut_links=n_cut_links,
        n_unsupported_links=n_unsupported,
    )
