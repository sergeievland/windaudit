"""Write every inconsistent cycle of a patch graph as a certificate a reader can add up.

Each graph equation says that the winding gauge of patch R minus that of patch
P equals an integer D. A certificate is a closed walk through such equations.
Walking an equation along its stored direction contributes +D; walking it
against its stored direction contributes -D. The certificate lists that signed
contribution for every step, so checking it needs nothing but addition: if
the contributions do not sum to zero, the equations on the walk cannot all
hold at once.

    python scripts/make_certificates.py results/wide_corrected
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict, deque
from pathlib import Path

RULE = ("Each step walks one equation 'gauge(to_patch) - gauge(from_patch) = contribution'. "
        "Add the contributions. A nonzero total proves that the equations on this walk "
        "cannot all hold simultaneously.")


def attachment_distances(directory):
    hits_path = directory / "point_surface_hits.json"
    if not hits_path.exists():
        return {}
    out = {}
    for collection in json.loads(hits_path.read_text()):
        for point in collection["points"]:
            if point["hits"]:
                best = sorted(point["hits"], key=lambda h: (-h["area"], h["distance"], h["id"]))[0]
                out[(collection["source_file"], int(point["id"]))] = (best["id"], round(float(best["distance"]), 3))
    return out


def main(directory):
    directory = Path(directory)
    edges = json.loads((directory / "measured_edges.json").read_text())
    distances = attachment_distances(directory)
    adjacency = defaultdict(list)
    for i, e in enumerate(edges):
        adjacency[e["P"]].append((e["R"], e["D"], i, "stored"))
        adjacency[e["R"]].append((e["P"], -e["D"], i, "reversed"))
    potential, tree = {}, defaultdict(list)
    for start in sorted(adjacency):
        if start in potential:
            continue
        potential[start] = 0
        queue = deque([start])
        while queue:
            a = queue.popleft()
            for b, d, i, way in adjacency[a]:
                if b not in potential:
                    potential[b] = potential[a] + d
                    queue.append(b)
                    tree[a].append((b, d, i, way))
                    tree[b].append((a, -d, i, "reversed" if way == "stored" else "stored"))

    def step(index, way, frm, to):
        e = edges[index]
        pts = [(e["source_file"], int(e["from_point_id"])), (e["source_file"], int(e["to_point_id"]))]
        return {
            "annotation": e.get("pcl_name"),
            "source_file": e["source_file"],
            "point_ids": [int(e["from_point_id"]), int(e["to_point_id"])],
            "from_patch": frm,
            "to_patch": to,
            "walked": way,
            "contribution": int(e["D"] if way == "stored" else -e["D"]),
            "annotation_winding_difference": int(e["raw_winding_delta"]) * (1 if way == "stored" else -1),
            "point_xyz": [list(reversed(e["from_zyx"])), list(reversed(e["to_zyx"]))],
            "attachment": [
                {"patch": distances[p][0], "distance_voxels": distances[p][1]} if p in distances else None
                for p in pts],
        }

    certificates = []
    for i, e in enumerate(edges):
        if potential[e["R"]] - potential[e["P"]] == e["D"]:
            continue
        parent = {e["R"]: None}
        queue = deque([e["R"]])
        while e["P"] not in parent:
            a = queue.popleft()
            for b, d, j, way in tree[a]:
                if b not in parent:
                    parent[b] = (a, j, way)
                    queue.append(b)
        back, node = [], e["P"]
        while parent[node] is not None:
            a, j, way = parent[node]
            back.append(step(j, way, a, node))
            node = a
        walk = [step(i, "stored", e["P"], e["R"])] + list(reversed(back))
        assert all(x["to_patch"] == y["from_patch"] for x, y in zip(walk, walk[1:] + walk[:1]))
        total = sum(s["contribution"] for s in walk)
        assert total != 0
        certificates.append({"sum_of_contributions": total, "steps": walk})
    out = {"rule": RULE, "graph_equations": len(edges), "certificates": certificates}
    (directory / "certificates.json").write_text(json.dumps(out, indent=2) + "\n")
    for n, c in enumerate(certificates, 1):
        print(f"certificate {n}: " + " + ".join(f"({s['contribution']})" for s in c["steps"])
              + f" = {c['sum_of_contributions']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
