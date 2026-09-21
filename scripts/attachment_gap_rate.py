"""How often would a branch ray split an annotation point from its attachment?

For every endpoint of every equation in a patch graph, the attachment gap is
the segment between the annotation point and the sample where its strip
actually ends on the surface. If the branch ray passes through that segment,
the equation is off by one winding unless the gap is transported. For a ray
placed uniformly at random around the axis, the expected number of straddled
endpoints is the sum of the gaps' angular widths divided by 2*pi. This gives a
frame-free estimate of how exposed a graph is to the defect, for upstream's
fitted frame as much as for the scan frame.

    python scripts/attachment_gap_rate.py results/wide_corrected
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_erosion

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from windaudit.scanspace import PatchGraph, ScanPatch, load_umbilicus, polyline_ijs, theta_at  # noqa: E402


def load_patch(pid):
    p = ScanPatch.load(ROOT / "data/scanspace_wide" / pid)
    valid = np.any(p.xyz != -1, axis=-1)
    cells = int(p.metadata.get("spiral_patch_erode_cells", 1))
    if cells > 0:
        p.xyz[~binary_erosion(valid, iterations=cells, border_value=0)] = -1
        p.__post_init__()
    return p


def main(run):
    run = Path(run)
    umb = load_umbilicus(ROOT / "data/paris4/umbilicus.json")
    reached = json.loads((run / "reached.json").read_text())
    edges = json.loads((run / "measured_edges.json").read_text())
    patches, graphs, widths, seen = {}, {}, [], set()
    for e in edges:
        for pid, ij, zyx, pt in ((e["P"], e["from_ij"], e["from_zyx"], e["from_point_id"]),
                                 (e["R"], e["to_ij"], e["to_zyx"], e["to_point_id"])):
            key = (e["pcl_id"], pt, pid)
            if key in seen:
                continue
            seen.add(key)
            if pid not in patches:
                patches[pid] = load_patch(pid)
                graphs[pid] = PatchGraph(patches[pid])
            entry = np.asarray(reached[pid]["entry_ij"], np.float32)
            centres = graphs[pid].route(entry, np.asarray(ij, np.float32))
            xyz, valid = patches[pid].lift(polyline_ijs(np.vstack([entry, centres, ij]), 1.0))
            end = xyz[valid][-1]
            point = np.asarray(zyx, np.float64)[::-1]
            a, b = float(theta_at(end, umb)), float(theta_at(point, umb))
            widths.append(abs((a - b + np.pi) % (2 * np.pi) - np.pi))
    w = np.asarray(widths)
    out = {"endpoints": int(len(w)), "mean_gap_radians": float(w.mean()), "max_gap_radians": float(w.max()),
           "expected_straddled_endpoints_uniform_ray": float(w.sum() / (2 * np.pi)),
           "probability_at_least_one_poisson": float(1 - np.exp(-w.sum() / (2 * np.pi)))}
    (run / "attachment_gap_rate.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
