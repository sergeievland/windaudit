"""Reader/writer for the VC3D point-collection JSON format (version "1").

Field semantics mirror the upstream loader in
``ScrollPrize/villa/spiral-fitting/point_collection.py``:

* collections and points are keyed by stringified integer ids;
* each point carries ``p = [x, y, z]`` in scan voxels and an optional
  integer-valued ``wind_a`` winding annotation;
* chain order is point-id order;
* ``metadata.winding_is_absolute`` marks absolute-winding collections;
* a collection whose points carry no ``wind_a`` at all is a same-winding
  collection (all points assert one shared winding).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

SUPPORTED_VERSION = "1"


@dataclass
class Point:
    id: int
    xyz: List[float]  # [x, y, z] scan voxels
    wind_a: Optional[float]  # None when the source carried no annotation

    @property
    def has_winding(self) -> bool:
        return self.wind_a is not None and math.isfinite(self.wind_a)


@dataclass
class Collection:
    id: int
    name: str
    points: Dict[int, Point]
    metadata: dict = field(default_factory=dict)
    color: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    source_file: str = ""

    @property
    def is_absolute(self) -> bool:
        return bool(self.metadata.get("winding_is_absolute", False))

    @property
    def has_winding_annotations(self) -> bool:
        """True when at least one point carried ``wind_a`` in the source json.

        Matches upstream ``normalise_pcl_winding_annotations``: this flag is the
        only signal separating a relative-winding collection from a same-winding
        one once missing annotations are zero-filled downstream.
        """
        return any(p.has_winding for p in self.points.values())

    def ordered_points(self) -> List[Point]:
        """Points in chain order (point-id order, matching upstream)."""
        return [self.points[pid] for pid in sorted(self.points)]


def load_point_collections(path: str) -> Dict[int, Collection]:
    """Load a point-collection file; raises ValueError on version mismatch."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    version = data.get("vc_pointcollections_json_version")
    if version != SUPPORTED_VERSION:
        raise ValueError(
            f"unsupported vc_pointcollections_json_version {version!r} in {path}"
        )
    collections: Dict[int, Collection] = {}
    for cid_str, col in data.get("collections", {}).items():
        cid = int(cid_str)
        if cid in collections:
            raise ValueError(f"duplicate normalised collection id {cid} in {path}")
        points: Dict[int, Point] = {}
        for pid_str, p in col.get("points", {}).items():
            pid = int(pid_str)
            if pid in points:
                raise ValueError(f"duplicate normalised point id {pid} in collection {cid}")
            xyz = [float(v) for v in p["p"]]
            if len(xyz) != 3 or not all(math.isfinite(v) for v in xyz):
                raise ValueError(f"point {cid}:{pid} must have three finite coordinates")
            wind = p.get("wind_a")
            if wind is not None and not math.isfinite(float(wind)):
                raise ValueError(f"point {cid}:{pid} has non-finite wind_a")
            points[pid] = Point(
                id=pid,
                xyz=xyz,
                wind_a=float(wind) if wind is not None else None,
            )
        collections[cid] = Collection(
            id=cid,
            name=col.get("name", str(cid)),
            points=points,
            metadata=col.get("metadata", {}) or {},
            color=col.get("color", [0.0, 0.0, 0.0]),
            source_file=path,
        )
    return collections


def save_point_collections(path: str, collections: Dict[int, Collection]) -> None:
    """Write collections back in the upstream JSON layout."""
    out = {"vc_pointcollections_json_version": SUPPORTED_VERSION, "collections": {}}
    for cid in sorted(collections):
        col = collections[cid]
        points = {}
        for pid in sorted(col.points):
            p = col.points[pid]
            entry: dict = {"p": list(p.xyz)}
            if p.has_winding:
                entry["wind_a"] = float(p.wind_a)
            points[str(pid)] = entry
        out["collections"][str(cid)] = {
            "name": col.name,
            "metadata": col.metadata,
            "color": col.color,
            "points": points,
        }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
