"""Frozen audit configuration, calibrated from absolute anchors only.

Every radial threshold is a fixed multiple of the wrap spacing implied by the
absolute-winding anchors. Absolute anchors never enter the cross-collection
evidence that the thresholds gate on this corpus, so calibration and findings
draw on separate annotation families. Protocol constants are module-level and
the resulting configuration is hashed into every report.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np

from .frame import KIND_ABSOLUTE, Node, NodeKey, NodePoint, sense_sign, wrap_angle

SECTOR_DZ = 40.0
SECTOR_DTHETA_RAD = math.radians(4.0)
CALIB_DZ = 200.0
CALIB_DTHETA_RAD = math.radians(8.0)
MAX_CHAIN_STEP_RAD = math.radians(120.0)
R_SAME_FACTOR = 0.35
CROWD_FACTOR = 0.70
MIN_EDGE_SUPPORT = 2
MIN_CALIBRATION_PAIRS = 10


@dataclass(frozen=True)
class FrozenConfig:
    wrap_spacing_median: float
    wrap_spacing_p10: float
    wrap_spacing_p90: float
    n_calibration_pairs: int
    anchor_order_agreement: float
    r_same: float
    crowd_min_spacing: float
    sector_dz: float
    sector_dtheta_rad: float
    max_chain_step_rad: float
    min_edge_support: int
    spiral_sense: str

    @property
    def sense(self) -> int:
        return sense_sign(self.spiral_sense)

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def config_hash(self) -> str:
        payload = json.dumps(self.as_dict(), sort_keys=True).encode()
        return hashlib.sha256(payload).hexdigest()


def _anchor_pairs(nodes: Dict[NodeKey, Node]):
    anchors: List[NodePoint] = []
    for key in sorted(nodes, key=str):
        if nodes[key].kind == KIND_ABSOLUTE:
            anchors.extend(nodes[key].points)
    for i in range(len(anchors)):
        for j in range(i + 1, len(anchors)):
            a, b = anchors[i], anchors[j]
            if abs(a.z - b.z) > CALIB_DZ:
                continue
            if abs(wrap_angle(a.theta - b.theta)) > CALIB_DTHETA_RAD:
                continue
            dw = a.wind_a - b.wind_a
            if dw == 0:
                continue
            yield a, b, dw


def calibrate_from_absolute_anchors(
    nodes: Dict[NodeKey, Node], spiral_sense: str = "CW"
) -> FrozenConfig:
    samples = []
    agree = 0
    for a, b, dw in _anchor_pairs(nodes):
        samples.append(abs(a.r - b.r) / abs(dw))
        if (a.r - b.r) * dw > 0:
            agree += 1
    if len(samples) < MIN_CALIBRATION_PAIRS:
        raise ValueError(
            f"calibration needs >= {MIN_CALIBRATION_PAIRS} same-sector "
            f"absolute-anchor pairs, got {len(samples)}"
        )
    arr = np.asarray(samples, dtype=np.float64)
    med = float(np.median(arr))
    sense_sign(spiral_sense)
    return FrozenConfig(
        wrap_spacing_median=med,
        wrap_spacing_p10=float(np.percentile(arr, 10)),
        wrap_spacing_p90=float(np.percentile(arr, 90)),
        n_calibration_pairs=int(len(arr)),
        anchor_order_agreement=agree / len(arr),
        r_same=R_SAME_FACTOR * med,
        crowd_min_spacing=CROWD_FACTOR * med,
        sector_dz=SECTOR_DZ,
        sector_dtheta_rad=SECTOR_DTHETA_RAD,
        max_chain_step_rad=MAX_CHAIN_STEP_RAD,
        min_edge_support=MIN_EDGE_SUPPORT,
        spiral_sense=spiral_sense.upper(),
    )
