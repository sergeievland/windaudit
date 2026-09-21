"""Cylindrical scroll frame from the umbilicus, and audit nodes.

Angles follow the upstream spiral-fitting convention
(``sample_spiral.get_theta``): ``theta = atan2(y, x)`` around the umbilicus,
taken in ``[0, 2*pi)``, with the integer winding branch cut at ``theta = 0``.
For a clockwise-outward scroll (upstream ``spiral_outward_sense = "CW"``)
radius grows with theta, so a sheet followed forward across the cut gains one
winding; for ``"ACW"`` it loses one.

Labels are integer windings expressed in that branch frame:

* relative-winding chains: ``label = wind_a + sense * crossings`` where
  ``crossings`` is the signed number of cut crossings walked along the chain
  (the annotation-level counterpart of the theta=0 branch transport in
  upstream ``find_inconsistent_windings.py``);
* same-winding collections: ``label = sense * crossings``;
* absolute-winding collections carry branch-frame windings already
  (``label = wind_a``). Because the scan-space cut only approximates the
  fitted spiral's cut, absolute collections are split into angular clusters
  and clusters close to the cut are marked seam-unsafe.

Every node carries one unknown integer gauge; only label differences inside a
node, corrected for the cut, carry information.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .pcl import Collection

TWO_PI = 2.0 * math.pi

KIND_RELATIVE = "relative"
KIND_SAME = "same"
KIND_ABSOLUTE = "absolute"

SPIRAL_SENSES = {"CW": 1, "ACW": -1}

ABSOLUTE_CLUSTER_GAP_RAD = math.radians(8.0)
ABSOLUTE_CUT_GUARD_RAD = math.radians(10.0)

NodeKey = Tuple[str, Union[int, str]]


def sense_sign(spiral_sense: str) -> int:
    try:
        return SPIRAL_SENSES[spiral_sense.upper()]
    except KeyError:
        raise ValueError("spiral sense must be 'CW' or 'ACW'") from None


class Umbilicus:
    """Piecewise-linear umbilicus: z -> (x, y) axis position."""

    def __init__(self, control_points_xyz: np.ndarray):
        cp = np.asarray(control_points_xyz, dtype=np.float64)
        if cp.ndim != 2 or cp.shape[1] != 3 or len(cp) < 2:
            raise ValueError("umbilicus needs an (N, 3) array of [x, y, z], N >= 2")
        order = np.argsort(cp[:, 2], kind="stable")
        self._z = cp[order, 2]
        self._x = cp[order, 0]
        self._y = cp[order, 1]

    @classmethod
    def load(cls, path: str) -> "Umbilicus":
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        cp = np.array(
            [[c["x"], c["y"], c["z"]] for c in data["control_points"]],
            dtype=np.float64,
        )
        return cls(cp)

    @property
    def z_range(self) -> Tuple[float, float]:
        return float(self._z[0]), float(self._z[-1])

    def axis_xy(self, z: float) -> Tuple[float, float]:
        return (
            float(np.interp(z, self._z, self._x)),
            float(np.interp(z, self._z, self._y)),
        )

    def cylindrical(self, xyz) -> Tuple[float, float]:
        """(r, theta) of a scan-voxel point, theta in [0, 2*pi)."""
        x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
        cx, cy = self.axis_xy(z)
        theta = math.atan2(y - cy, x - cx) % TWO_PI
        if theta >= TWO_PI:
            theta = 0.0
        return math.hypot(x - cx, y - cy), theta


def wrap_angle(d: float) -> float:
    """Wrap an angle difference into (-pi, pi]."""
    d = math.fmod(d, TWO_PI)
    if d > math.pi:
        d -= TWO_PI
    elif d <= -math.pi:
        d += TWO_PI
    return d


def cut_crossing(theta_from: float, theta_to: float) -> int:
    """Signed crossings of the theta = 0 cut on the short arc from
    ``theta_from`` to ``theta_to`` (both in [0, 2*pi)): +1 forward, -1
    backward, 0 when the arc stays inside the branch."""
    diff = theta_to - theta_from
    if diff < -math.pi:
        return 1
    if diff > math.pi:
        return -1
    return 0


@dataclass
class NodePoint:
    pid: int
    xyz: List[float]
    r: float
    theta: float
    z: float
    wind_a: float
    crossings: int
    label: int


def label_diff(qa: NodePoint, qb: NodePoint, sense: int) -> int:
    """Winding of ``qa`` minus winding of ``qb``, both expressed in the branch
    of ``qa``, for two points of one node (or, with gauges added, of two
    nodes). Zero means the points sit on the same sheet."""
    return qa.label - qb.label + sense * cut_crossing(qa.theta, qb.theta)


@dataclass
class Node:
    """One collection (or angular cluster) with a single unknown gauge."""

    kind: str
    cid: Union[int, str]
    name: str
    source_cid: int
    points: List[NodePoint] = field(default_factory=list)
    seam_safe: bool = True

    @property
    def key(self) -> NodeKey:
        return (self.kind, self.cid)

    @property
    def n(self) -> int:
        return len(self.points)


def _chain_points(
    collection: Collection,
    umbilicus: Umbilicus,
    use_wind: bool,
    sense: int,
    count_crossings: bool,
    max_step_rad: float,
) -> Tuple[List[NodePoint], bool]:
    pts: List[NodePoint] = []
    crossings = 0
    prev_theta: Optional[float] = None
    seam_safe = True
    for p in collection.ordered_points():
        if use_wind and not p.has_winding:
            continue
        r, theta = umbilicus.cylindrical(p.xyz)
        if prev_theta is not None:
            if abs(wrap_angle(theta - prev_theta)) > max_step_rad:
                seam_safe = False
            crossings += cut_crossing(prev_theta, theta)
        prev_theta = theta
        wind = float(p.wind_a) if use_wind else 0.0
        if use_wind and abs(wind - round(wind)) > 1e-6:
            raise ValueError(
                f"collection {collection.id} point {p.id}: winding annotation "
                f"{wind} is not an integer"
            )
        c = crossings if count_crossings else 0
        pts.append(
            NodePoint(
                pid=p.id,
                xyz=list(p.xyz),
                r=r,
                theta=theta,
                z=float(p.xyz[2]),
                wind_a=wind,
                crossings=c,
                label=int(round(wind)) + sense * c,
            )
        )
    return pts, seam_safe


def build_chain_node(
    collection: Collection,
    umbilicus: Umbilicus,
    kind: str,
    sense: int,
    max_step_rad: float = math.radians(120.0),
) -> Optional[Node]:
    """Node for a relative-winding or same-winding collection."""
    if kind not in (KIND_RELATIVE, KIND_SAME):
        raise ValueError(f"chain node kind must be relative or same, got {kind}")
    pts, safe = _chain_points(
        collection, umbilicus, use_wind=(kind == KIND_RELATIVE), sense=sense,
        count_crossings=True, max_step_rad=max_step_rad,
    )
    if not pts:
        return None
    return Node(kind=kind, cid=collection.id, name=collection.name,
                source_cid=collection.id, points=pts, seam_safe=safe)


def build_absolute_nodes(
    collection: Collection, umbilicus: Umbilicus
) -> List[Node]:
    """Angular clusters of one absolute-winding collection."""
    pts, _ = _chain_points(
        collection, umbilicus, use_wind=True, sense=1,
        count_crossings=False, max_step_rad=math.pi,
    )
    if not pts:
        return []
    order = sorted(range(len(pts)), key=lambda i: pts[i].theta)
    thetas = [pts[i].theta for i in order]
    # Circular clustering: open the circle at its widest gap.
    gaps = [
        (thetas[(k + 1) % len(thetas)] - thetas[k]) % TWO_PI
        for k in range(len(thetas))
    ]
    start = (int(np.argmax(gaps)) + 1) % len(thetas) if len(thetas) > 1 else 0
    rotated = order[start:] + order[:start]
    clusters: List[List[int]] = [[rotated[0]]]
    for prev, cur in zip(rotated, rotated[1:]):
        if (pts[cur].theta - pts[prev].theta) % TWO_PI > ABSOLUTE_CLUSTER_GAP_RAD:
            clusters.append([cur])
        else:
            clusters[-1].append(cur)
    nodes = []
    for k, members in enumerate(clusters):
        members = sorted(members, key=lambda i: pts[i].pid)
        cpts = [pts[i] for i in members]
        near_cut = any(
            min(q.theta, TWO_PI - q.theta) < ABSOLUTE_CUT_GUARD_RAD for q in cpts
        )
        nodes.append(
            Node(
                kind=KIND_ABSOLUTE,
                cid=f"{collection.id}.{k}",
                name=f"{collection.name}#{k}",
                source_cid=collection.id,
                points=cpts,
                seam_safe=not near_cut,
            )
        )
    return nodes


def classify(collection: Collection) -> str:
    """Upstream classification (``classify_pcl``): absolute by metadata,
    otherwise relative when any point carries ``wind_a``, else same."""
    if collection.is_absolute:
        return KIND_ABSOLUTE
    if collection.has_winding_annotations:
        return KIND_RELATIVE
    return KIND_SAME


def build_nodes(
    relative: Dict[int, Collection],
    same: Dict[int, Collection],
    absolute: Dict[int, Collection],
    umbilicus: Umbilicus,
    spiral_sense: str = "CW",
    max_step_rad: float = math.radians(120.0),
) -> Dict[NodeKey, Node]:
    """Build all audit nodes. Each input file is checked against the upstream
    content classification; a collection whose content contradicts the file
    it came from is an input error. Empty collections are dropped."""
    sense = sense_sign(spiral_sense)
    nodes: Dict[NodeKey, Node] = {}
    for expected, cols in (
        (KIND_RELATIVE, relative),
        (KIND_SAME, same),
        (KIND_ABSOLUTE, absolute),
    ):
        for cid in sorted(cols):
            col = cols[cid]
            if not col.points:
                continue
            kind = classify(col)
            if kind != expected:
                raise ValueError(
                    f"collection {cid} ({col.name!r}) in the {expected} input "
                    f"classifies as {kind} by content"
                )
            if kind == KIND_ABSOLUTE:
                for node in build_absolute_nodes(col, umbilicus):
                    nodes[node.key] = node
            else:
                node = build_chain_node(col, umbilicus, kind, sense, max_step_rad)
                if node is not None:
                    nodes[node.key] = node
    return nodes


def node_tag(key: NodeKey) -> str:
    return f"{key[0]}:{key[1]}"
