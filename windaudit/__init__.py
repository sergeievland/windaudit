"""windaudit: pre-fit certification of winding annotations for scroll unwrapping.

Audits VC3D point-collection winding annotations (absolute, relative and
same-winding) against the papyrus non-crossing invariant and exact integer
cycle consistency, using only the annotations and the umbilicus: no trained
model, no fitted checkpoint, no GPU.
"""

__version__ = "0.2.1"

from .calibrate import FrozenConfig, calibrate_from_absolute_anchors
from .frame import Umbilicus, build_nodes
from .graph import WindingGraph
from .matching import build_links
from .pcl import load_point_collections, save_point_collections
from .repair import solve_min_edge_cut, solve_min_quarantine
from .report import run_audit

__all__ = [
    "__version__",
    "FrozenConfig",
    "Umbilicus",
    "WindingGraph",
    "build_links",
    "build_nodes",
    "calibrate_from_absolute_anchors",
    "load_point_collections",
    "run_audit",
    "save_point_collections",
    "solve_min_edge_cut",
    "solve_min_quarantine",
]
