"""Check that the patch-graph results do not depend on the angular frame.

The scan-space transport measures θ around the umbilicus; upstream measures it
after a fitted spiral transform. SCANSPACE_TRANSPORT.md proves that any
admissible transform changes every graph equation only by a difference of
per-patch potentials, so cycle sums, the inconsistent cycles and the solver's
optimal edit sets are frame-independent. This script tests that claim on the
real patch graph: it re-runs the whole wide audit under several synthetic
admissible transforms — rotations about the umbilicus by a smooth,
position-dependent angle that moves many points across the branch ray — and
compares the outputs with the untransformed run.

    python scripts/frame_invariance.py results/wide out/frame_invariance

Needs the public patch surfaces (scripts/fetch_data.py). Surface projection is
frame-independent and is reused from the audit cache; every angle-dependent
quantity is recomputed under each transform.
"""

from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import windaudit.scanspace as scanspace  # noqa: E402
import wide_patch_audit  # noqa: E402

UNTWISTED = scanspace.theta_at

# Rotation about the umbilicus by alpha(p). For fixed radius and z the map
# theta -> theta + alpha is a bijection of the circle because |d alpha / d theta| < 1,
# so each transform is a homeomorphism that fixes the umbilicus and carries a
# small loop around it to a loop winding once: admissible in the sense of the
# proposition.
TRANSFORMS = {
    "twist_z": dict(a=1.1, b=0.0, c=0.0),
    "twist_r": dict(a=0.0, b=1.3, c=0.0),
    "twist_mixed": dict(a=0.9, b=0.7, c=0.45),
}


def make_theta(a, b, c):
    def theta(xyz, umbilicus):
        base = UNTWISTED(xyz, umbilicus)
        xyz = np.asarray(xyz, dtype=np.float64)
        yx = xyz[..., 1::-1] - umbilicus(xyz[..., 2])
        r = np.linalg.norm(yx, axis=-1)
        z = xyz[..., 2]
        alpha = a * np.sin(z / 700.0 + 0.3) + b * np.cos(r / 150.0) + c * np.sin(base)
        return np.mod(base + alpha, 2 * np.pi)
    return theta


def use(theta):
    scanspace.theta_at = theta
    wide_patch_audit.theta_at = theta


def load(outdir):
    out = Path(outdir)
    return {
        "edges": json.loads((out / "measured_edges.json").read_text()),
        "initial": json.loads((out / "initial_graph.json").read_text()),
        "before": json.loads((out / "before_full.json").read_text()),
        "after": json.loads((out / "after_full.json").read_text()),
    }


def edge_id(e):
    return (e["pcl_id"], e["from_point_id"], e["to_point_id"])


def potential_residuals(pairs):
    """Edges (P, R, value); return how many violate value == phi(R) - phi(P)."""
    adj = collections.defaultdict(list)
    for p, r, v in pairs:
        adj[p].append((r, v))
        adj[r].append((p, -v))
    phi, bad = {}, 0
    for start in adj:
        if start in phi:
            continue
        phi[start] = 0
        stack = [start]
        while stack:
            u = stack.pop()
            for w, v in adj[u]:
                if w in phi:
                    bad += phi[w] != phi[u] + v
                else:
                    phi[w] = phi[u] + v
                    stack.append(w)
    return bad // 2


def holonomy_signature(edges):
    """Multiset of nonzero fundamental-cycle sums over a fixed spanning forest."""
    order = sorted(edges, key=edge_id)
    adj = collections.defaultdict(list)
    for e in order:
        adj[e["P"]].append((e["R"], e["D"], edge_id(e)))
        adj[e["R"]].append((e["P"], -e["D"], edge_id(e)))
    phi, tree = {}, set()
    for start in sorted(adj):
        if start in phi:
            continue
        phi[start] = 0
        queue = collections.deque([start])
        while queue:
            u = queue.popleft()
            for w, v, k in sorted(adj[u], key=lambda t: (str(t[0]), t[2])):
                if w not in phi:
                    phi[w] = phi[u] + v
                    tree.add(k)
                    queue.append(w)
    sums = {}
    for e in order:
        if edge_id(e) not in tree:
            s = phi[e["P"]] + e["D"] - phi[e["R"]]
            if s:
                sums[edge_id(e)] = s
    return sums


def edits(result):
    return sorted((e["rel_pcl_id"], e["from_point_id"], e["to_point_id"]) for e in result["edges"])


def main(reference_dir, outroot):
    reference = load(reference_dir)
    ref_edges = {edge_id(e): e for e in reference["edges"]}
    ref_holonomy = holonomy_signature(reference["edges"])
    before_py = ROOT / "upstream/villa/find_inconsistent_windings.py"
    after_py = ROOT / "out/wide_solver/spiral-fitting/find_inconsistent_windings.py"
    report = {"reference": str(reference_dir), "transforms": {}}
    for name, params in TRANSFORMS.items():
        use(make_theta(**params))
        outdir = Path(outroot) / name
        wide_patch_audit.main(ROOT / "data/scanspace_wide", before_py, after_py, outdir)
        use(UNTWISTED)
        got = load(outdir)
        got_edges = {edge_id(e): e for e in got["edges"]}
        same_edges = set(got_edges) == set(ref_edges)
        changed = sum(got_edges[k]["D"] != ref_edges[k]["D"] for k in ref_edges if k in got_edges)
        gauge_violations = potential_residuals(
            [(e["P"], e["R"], got_edges[k]["D"] - e["D"]) for k, e in ref_edges.items() if k in got_edges])
        holonomy = holonomy_signature(got["edges"])
        row = {
            "parameters": params,
            "same_edge_set": same_edges,
            "edges": len(got_edges),
            "edges_whose_value_changed": int(changed),
            "value_changes_not_explained_by_patch_potentials": int(gauge_violations),
            "inconsistent_fundamental_cycles": got["initial"]["inconsistent_non_tree_edges"],
            "fundamental_cycle_sums_identical": holonomy == ref_holonomy,
        }
        for side in ("before", "after"):
            row[side] = {
                "equations_in_model": got[side]["equations_in_model"],
                "num_edges_changed": got[side]["num_edges_changed"],
                "same_edges_edited_as_untransformed": edits(got[side]) == edits(reference[side]),
                "retained_graph_consistent": got[side]["independent_retained_graph_check"]["consistent"],
            }
        report["transforms"][name] = row
        print(name, json.dumps(row), flush=True)
    report["invariant"] = all(
        r["same_edge_set"] and r["value_changes_not_explained_by_patch_potentials"] == 0
        and r["fundamental_cycle_sums_identical"]
        and all(r[s]["equations_in_model"] == reference[s]["equations_in_model"]
                and r[s]["num_edges_changed"] == reference[s]["num_edges_changed"] for s in ("before", "after"))
        for r in report["transforms"].values())
    Path(outroot).mkdir(parents=True, exist_ok=True)
    (Path(outroot) / "frame_invariance.json").write_text(json.dumps(report, indent=2) + "\n")
    print("INVARIANT" if report["invariant"] else "NOT INVARIANT")
    return 0 if report["invariant"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
