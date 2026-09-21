import math

import numpy as np
import pytest

from windaudit.calibrate import FrozenConfig
from windaudit.frame import Umbilicus
from windaudit.pcl import Collection, Point

SPACING = 20.0


@pytest.fixture
def umbilicus():
    return Umbilicus(np.array([[0.0, 0.0, -1000.0], [0.0, 0.0, 20000.0]]))


@pytest.fixture
def config():
    return FrozenConfig(
        wrap_spacing_median=SPACING,
        wrap_spacing_p10=0.9 * SPACING,
        wrap_spacing_p90=1.1 * SPACING,
        n_calibration_pairs=100,
        anchor_order_agreement=1.0,
        r_same=0.35 * SPACING,
        crowd_min_spacing=0.70 * SPACING,
        sector_dz=40.0,
        sector_dtheta_rad=math.radians(4.0),
        max_chain_step_rad=math.radians(120.0),
        min_edge_support=2,
        spiral_sense="CW",
    )


def wrap_xyz(wrap: float, theta: float, z: float, tilt: float = 0.0,
             spacing: float = SPACING, r0: float = 100.0):
    """A point on an ideal clockwise-outward spiral (upstream "CW"): with
    theta in [0, 2*pi), radius grows by one spacing per integer winding and
    continuously with theta, plus an optional radial drift with z."""
    theta = theta % (2.0 * math.pi)
    r = r0 + spacing * (wrap + theta / (2.0 * math.pi)) + tilt * z
    return [r * math.cos(theta), r * math.sin(theta), z]


def make_collection(cid, rows, absolute=False, name=None):
    """rows: list of (x, y, z, wind_a-or-None)."""
    points = {}
    for i, (x, y, z, w) in enumerate(rows, start=1):
        points[i] = Point(id=i, xyz=[float(x), float(y), float(z)],
                          wind_a=None if w is None else float(w))
    metadata = {"winding_is_absolute": True} if absolute else {}
    return Collection(id=cid, name=name or f"col{cid}", points=points,
                      metadata=metadata)


def ladder(cid, theta, z, wraps, tilt=0.0):
    """A radial relative-winding ladder: one point per wrap at fixed theta, z."""
    rows = [tuple(wrap_xyz(w, theta, z, tilt)) + (w,) for w in wraps]
    return make_collection(cid, rows)


def wrap_trace(cid, wrap, theta, zs, tilt=0.0):
    """A same-winding trace running up one wrap at fixed theta."""
    rows = [tuple(wrap_xyz(wrap, theta, z, tilt)) + (None,) for z in zs]
    return make_collection(cid, rows)
