"""Winding gauge graph: components, potentials, cycle holonomy.

Nodes are collections (one unknown integer gauge each); edges are internally
consistent same-sheet offsets. Around every independent cycle the signed
offsets must sum to zero; a nonzero sum (holonomy) is an integer proof that
the annotations along that cycle cannot all be right.

A cycle is tautological - it cannot fail - when every node on it carries a
constant label and every offset on it is zero (for example, a loop of
same-winding traces that never crosses the branch cut). Tautological cycles
are counted separately and excluded from the effective-cycle statistics.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .frame import KIND_SAME, Node, NodeKey
from .matching import EdgeKey, MatchResult


@dataclass
class CycleFinding:
    u: NodeKey
    v: NodeKey
    holonomy: int
    support: int
    tree_path: List[NodeKey]
    effective: bool


def constant_label_nodes(nodes: Dict[NodeKey, Node]) -> Set[NodeKey]:
    return {
        key for key, node in nodes.items()
        if len({q.label for q in node.points}) <= 1
    }


class WindingGraph:
    def __init__(self, match: MatchResult, constant_nodes: Optional[Iterable[NodeKey]] = None):
        self.edges: Dict[EdgeKey, Tuple[int, int]] = dict(match.edges)
        if constant_nodes is None:
            constant = {k for e in self.edges for k in e if k[0] == KIND_SAME}
        else:
            constant = set(constant_nodes)
        self._constant = constant
        self.adj: Dict[NodeKey, List[Tuple[NodeKey, int, int]]] = defaultdict(list)
        for (a, b), (off, support) in self.edges.items():
            self.adj[a].append((b, off, support))
            self.adj[b].append((a, -off, support))
        for k in self.adj:
            self.adj[k].sort()
        self.potential: Dict[NodeKey, int] = {}
        self.components: List[List[NodeKey]] = []
        self.tree_parent: Dict[NodeKey, NodeKey] = {}
        self._tree_edges: Set[frozenset] = set()
        self._build()

    def _build(self) -> None:
        for start in sorted(self.adj):
            if start in self.potential:
                continue
            comp = [start]
            self.potential[start] = 0
            stack = [start]
            while stack:
                u = stack.pop()
                for v, off, _ in self.adj[u]:
                    if v in self.potential:
                        continue
                    self.potential[v] = self.potential[u] + off
                    self.tree_parent[v] = u
                    self._tree_edges.add(frozenset((u, v)))
                    comp.append(v)
                    stack.append(v)
            self.components.append(sorted(comp))

    def _path_to_root(self, k: NodeKey) -> List[NodeKey]:
        path = [k]
        while path[-1] in self.tree_parent:
            path.append(self.tree_parent[path[-1]])
        return path

    def tree_path(self, u: NodeKey, v: NodeKey) -> List[NodeKey]:
        pu, pv = self._path_to_root(u), self._path_to_root(v)
        index_v = {k: i for i, k in enumerate(pv)}
        for i, k in enumerate(pu):
            if k in index_v:
                return pu[: i + 1] + pv[: index_v[k]][::-1]
        raise ValueError("nodes are not in one component")

    def _edge_offset(self, a: NodeKey, b: NodeKey) -> int:
        if (a, b) in self.edges:
            return self.edges[(a, b)][0]
        return -self.edges[(b, a)][0]

    def cycles(self) -> List[CycleFinding]:
        """All fundamental cycles (one per non-tree edge), with holonomy."""
        out: List[CycleFinding] = []
        for (a, b), (off, support) in sorted(self.edges.items()):
            if frozenset((a, b)) in self._tree_edges:
                continue
            holonomy = self.potential[b] - self.potential[a] - off
            path = self.tree_path(a, b)
            offsets = [off] + [self._edge_offset(x, y) for x, y in zip(path, path[1:])]
            effective = (
                any(k not in self._constant for k in path)
                or any(o != 0 for o in offsets)
            )
            out.append(CycleFinding(a, b, holonomy, support, path, effective))
        return out

    def stats(self) -> dict:
        cycles = self.cycles()
        effective = [c for c in cycles if c.effective]
        return {
            "n_nodes": len(self.adj),
            "n_edges": len(self.edges),
            "n_components": len(self.components),
            "largest_component": max((len(c) for c in self.components), default=0),
            "n_cycles": len(cycles),
            "n_effective_cycles": len(effective),
            "n_nonzero_holonomy": sum(1 for c in cycles if c.holonomy != 0),
            "n_effective_nonzero_holonomy": sum(1 for c in effective if c.holonomy != 0),
        }
