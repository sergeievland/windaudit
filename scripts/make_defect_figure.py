"""Render the figure for the upstream defects shown on real data.

    python scripts/make_defect_figure.py results/wide_corrected/upstream_defects.png

Needs the fetched patch surfaces (the left panel draws the real strip samples).
"""

import io
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.ndimage import binary_erosion  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from windaudit.scanspace import PatchGraph, ScanPatch, load_umbilicus, polyline_ijs, theta_at  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e6e5e0"
NEUTRAL = "#9a9892"
BLUE = "#2a78d6"
ORANGE = "#eb6834"


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=0)


def gap_geometry():
    umb = load_umbilicus(ROOT / "data/paris4/umbilicus.json")
    reached = json.loads((ROOT / "results/wide/reached.json").read_text())
    edge = next(e for e in json.loads((ROOT / "results/wide/measured_edges.json").read_text())
                if e["pcl_id"] == 82 and e["from_point_id"] == 1561 and e["to_point_id"] == 1560)
    pid = edge["R"]
    patch = ScanPatch.load(ROOT / "data/scanspace_wide" / pid)
    valid = np.any(patch.xyz != -1, axis=-1)
    cells = int(patch.metadata.get("spiral_patch_erode_cells", 1))
    if cells > 0:
        patch.xyz[~binary_erosion(valid, iterations=cells, border_value=0)] = -1
        patch.__post_init__()
    entry = np.asarray(reached[pid]["entry_ij"], np.float32)
    to = np.asarray(edge["to_ij"], np.float32)
    centres = PatchGraph(patch).route(entry, to)
    xyz, ok = patch.lift(polyline_ijs(np.vstack([entry, centres, to]), 1.0))
    xyz = xyz[ok][-4:]
    point = np.asarray(edge["to_zyx"], float)[::-1]

    def rel(q):
        q = np.asarray(q, float)
        u = umb(q[..., 2])
        return np.stack([q[..., 0] - u[..., 1], q[..., 1] - u[..., 0]], -1)

    return rel(xyz), rel(point), float(theta_at(point, umb)), float(theta_at(xyz[-1], umb))


def main(out_path):
    path, point, th_point, th_end = gap_geometry()
    fig = plt.figure(figsize=(12.5, 4.6), facecolor=SURFACE)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28)

    ax = fig.add_subplot(grid[0])
    style(ax)
    x0 = point[0]
    ax.axhline(0, color=NEUTRAL, lw=1.4, zorder=1)
    ax.text(-3.3, 0.12, "branch ray  θ = 0", color=INK_2, fontsize=9, ha="left", va="bottom")
    ax.plot(path[:, 0] - x0, path[:, 1], color=BLUE, lw=2, zorder=2)
    ax.scatter(path[:-1, 0] - x0, path[:-1, 1], s=14, color=BLUE, edgecolors=SURFACE, linewidths=1.2, zorder=3)
    end = path[-1]
    ax.scatter([end[0] - x0], [end[1]], s=90, color=BLUE, edgecolors=SURFACE, linewidths=2, zorder=4)
    ax.scatter([0], [point[1]], s=90, color=ORANGE, edgecolors=SURFACE, linewidths=2, zorder=4)
    ax.plot([end[0] - x0, 0], [end[1], point[1]], color=ORANGE, lw=2, ls=(0, (3, 2)), zorder=2)
    ax.annotate(f"strip on the surface ends here\nθ = {th_end:.5f}", xy=(end[0] - x0, end[1]),
                xytext=(1.9, 1.25), color=INK, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8))
    ax.annotate(f"annotation point col109 #1560\nθ = {th_point:.5f}", xy=(0, point[1]),
                xytext=(-3.2, -2.6), color=INK, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8))
    ax.text(-3.3, -1.25, "1.37 voxels across the ray,\nnever transported", color=ORANGE,
            fontsize=9, ha="left", va="center", fontweight="bold")
    ax.text(3.3, -2.1, "strip samples", color=BLUE, fontsize=9, ha="left", va="center")
    ax.set_xlabel("radial offset from the annotation point (voxels)", color=INK_2, fontsize=9)
    ax.set_ylabel("distance from the branch ray (voxels)", color=INK_2, fontsize=9)
    ax.set_title("The attachment gap: one step no count includes", loc="left", color=INK,
                 fontsize=11, fontweight="bold", pad=10)
    ax.set_xlim(-3.5, 6.5)
    ax.set_ylim(-3.2, 2.1)
    ax.set_aspect("equal", adjustable="datalim")

    ax = fig.add_subplot(grid[1])
    style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    rows = [("As upstream builds it", 7, 645), ("With the gap transported", 2, 645)]
    ys = np.arange(len(rows))[::-1] * 1.0
    h = 0.32
    for y, (label, up, fixed) in zip(ys, rows):
        ax.barh(y + h / 2 + 0.02, fixed, height=h, color=BLUE, edgecolor=SURFACE, lw=2)
        ax.barh(y - h / 2 - 0.02, max(up, 4), height=h, color=ORANGE, edgecolor=SURFACE, lw=2)
        ax.text(fixed + 10, y + h / 2 + 0.02, f"{fixed} of 645", va="center", color=INK, fontsize=10, fontweight="bold")
        ax.text(max(up, 4) + 10, y - h / 2 - 0.02, f"{up} of 645", va="center", color=INK, fontsize=10, fontweight="bold")
        ax.text(4, y + h + 0.1, label, va="bottom", ha="left", color=INK_2, fontsize=9)
    ax.set_yticks([])
    ax.set_xlim(0, 760)
    ax.set_ylim(-0.55, 1.62)
    ax.set_xlabel("real equations kept in the repair model", color=INK_2, fontsize=9)
    ax.set_title("What the repair solver keeps", loc="left", color=INK, fontsize=11, fontweight="bold", pad=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    ax.legend(handles, ["patched solver", "upstream solver"], loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=2, frameon=False, fontsize=9, labelcolor=INK)

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=160, bbox_inches="tight", facecolor=SURFACE)
    buffer.seek(0)
    image = Image.open(buffer).convert("RGB")
    image.quantize(colors=192, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(out_path, optimize=True)


if __name__ == "__main__":
    main(sys.argv[1])
