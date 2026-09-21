"""Render the summary figure from a committed audit report.

usage: python scripts/make_figure.py results/paris4/audit_report.json data/paris4 results/paris4/summary.png

The noise bar reads results/frozen_null.json (500 frozen rounds), next to the report.
"""

import io
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from windaudit.frame import Umbilicus, build_nodes  # noqa: E402
from windaudit.pcl import load_point_collections  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e6e5e0"
NEUTRAL = "#9a9892"
SERIES = {"relative": "#2a78d6", "same": "#eb6834", "absolute": "#1baf7a"}


def style_axis(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=0)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def bars(ax, labels, values, colors, notes, title):
    style_axis(ax)
    xs = range(len(labels))
    ax.bar(xs, values, width=0.56, color=colors, edgecolor=SURFACE, linewidth=2)
    for x, v, note in zip(xs, values, notes):
        ax.text(x, v + 2.5, f"{v:.0f}%", ha="center", va="bottom", color=INK,
                fontsize=11, fontweight="bold")
        ax.text(x, -9, note, ha="center", va="top", color=INK_2, fontsize=8)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(labels, color=INK, fontsize=9)
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100%"])
    ax.set_title(title, loc="left", color=INK, fontsize=11, fontweight="bold", pad=10)


def main(report_path, data_dir, out_path):
    with open(report_path, encoding="utf-8") as f:
        report = json.load(f)
    shuffled = report["tier1"]["shuffled_control"]
    planted = report["planted_defect_validation"]

    rel = load_point_collections(os.path.join(data_dir, "relative_windings.json"))
    same = load_point_collections(os.path.join(data_dir, "same_windings.json"))
    ab = load_point_collections(os.path.join(data_dir, "abs_winding.json"))
    umb = Umbilicus.load(os.path.join(data_dir, "umbilicus.json"))
    nodes = build_nodes(rel, same, ab, umb, spiral_sense=report["calibration"]["spiral_sense"])

    fig = plt.figure(figsize=(13, 4.4), facecolor=SURFACE)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.6], wspace=0.32)

    ax = fig.add_subplot(grid[0])
    bars(
        ax,
        ["annotations", "labels shuffled"],
        [100 * shuffled["real_violation_rate"], 100 * shuffled["shuffled_mean_rate"]],
        [SERIES["relative"], NEUTRAL],
        [f"{report['tier1']['order_audit']['n_pairs_tested']:,} pairs",
         f"{shuffled['n_rounds']} rounds"],
        "Order violations",
    )

    ax = fig.add_subplot(grid[1])
    sw, ss = planted["skipped_wrap"], planted["sheet_switch"]
    null_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(report_path))), "frozen_null.json")
    with open(null_path, encoding="utf-8") as f:
        null = json.load(f)
    bars(
        ax,
        ["skipped wrap", "sheet switch", "noise only"],
        [100 * sw["alarm_rate"], 100 * ss["alarm_rate"],
         100 * null["rounds_with_alarm"] / null["rounds"]],
        [SERIES["relative"], SERIES["relative"], NEUTRAL],
        [f"{sw['n_trials']} plants", f"{ss['n_trials']} plants",
         f"{null['rounds']} rounds"],
        "Alarms: planted defects vs noise",
    )

    ax = fig.add_subplot(grid[2])
    style_axis(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    for kind, label in (("same", "same-winding"), ("relative", "relative-winding"),
                        ("absolute", "absolute-winding")):
        xs, zs = [], []
        for node in nodes.values():
            if node.kind != kind:
                continue
            for q in node.points:
                xs.append(q.theta * 180.0 / 3.141592653589793)
                zs.append(q.z)
        ax.scatter(xs, zs, s=9 if kind != "absolute" else 26, color=SERIES[kind],
                   edgecolors=SURFACE, linewidths=0.4, label=f"{label} ({len(xs):,})",
                   zorder=3 if kind == "absolute" else 2)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel("theta around the umbilicus (degrees)", color=INK_2, fontsize=9)
    ax.set_ylabel("z (scan voxels)", color=INK_2, fontsize=9)
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8,
              frameon=False, labelcolor=INK, handletextpad=0.3, columnspacing=1.2)
    n_points = sum(report['corpus'][k]['n_points'] for k in ('relative', 'same', 'absolute'))
    ax.set_title(f"PHercParis4 annotations audited: {n_points:,} points",
                 loc="left", color=INK, fontsize=11, fontweight="bold", pad=10)

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=160, bbox_inches="tight", facecolor=SURFACE)
    buffer.seek(0)
    image = Image.open(buffer).convert("RGB")
    palette = image.quantize(colors=192, method=Image.Quantize.MEDIANCUT,
                             dither=Image.Dither.NONE)
    palette.save(out_path, optimize=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
