"""Exact, deterministic repair of an inconsistent winding graph.

Both repairs act on the same integer-potential model: every node carries an
integer gauge, and a kept edge ``(a, b, offset)`` requires
``gauge_b - gauge_a == offset``.

* Minimum quarantine removes collections. Its objective is lexicographic:
  fewest collections, then least link support lost, then a canonical choice
  by node order.
* Minimum edge cut removes edges. Its objective is lexicographic: least link
  support lost, then a canonical choice by edge order.

Each stage is an integer program solved with HiGHS (``scipy.optimize.milp``)
with a zero optimality gap, and every stage's optimum is fixed as a
constraint before the next. The canonical stage walks the candidates in a
fixed order and keeps a candidate out of the removal whenever an optimal
repair without it exists. The result is therefore unique for a given graph,
independent of solver version or tie-breaking inside the solver. A repair is
reported as optimal only if every stage succeeds and an independent traversal
confirms that the remainder is consistent.
"""

from __future__ import annotations

import io
import itertools
import os
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from numbers import Integral
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix, vstack

from .frame import NodeKey

Edges = Dict[Tuple[NodeKey, NodeKey], Tuple[int, int]]

SOLVER_OPTIONS = {"mip_rel_gap": 0.0}


@dataclass
class RepairResult:
    removed: List[NodeKey]
    cut_edges: List[Tuple[NodeKey, NodeKey]]
    objective: Optional[float]
    optimal: bool
    status: str
    support_removed: int


def _validate_edges(edges: Edges) -> None:
    for (a, b), (offset, support) in edges.items():
        if not isinstance(offset, Integral) or isinstance(offset, bool):
            raise ValueError("edge offsets must be integers")
        if not isinstance(support, Integral) or isinstance(support, bool) or support < 0:
            raise ValueError("edge support must be a non-negative integer")
    if sum(abs(off) + sup for off, sup in edges.values()) > 2**40:
        raise ValueError("graph exceeds supported numerical range")


def edges_consistent(edges: Edges) -> bool:
    """True when integer gauges exist that satisfy every edge equation."""
    adj: Dict[NodeKey, List[Tuple[NodeKey, int]]] = {}
    for (a, b), (off, _) in edges.items():
        adj.setdefault(a, []).append((b, off))
        adj.setdefault(b, []).append((a, -off))
    pot: Dict[NodeKey, int] = {}
    for start in adj:
        if start in pot:
            continue
        pot[start] = 0
        stack = [start]
        while stack:
            u = stack.pop()
            for v, off in adj[u]:
                if v in pot:
                    if pot[v] != pot[u] + off:
                        return False
                else:
                    pot[v] = pot[u] + off
                    stack.append(v)
    return True


def _potential_bound(edges: Edges) -> float:
    """Any consistent sub-graph admits integer potentials within +/- this
    bound (one gauge per component fixed at zero): every potential is a sum of
    offsets along a simple path."""
    n_nodes = len({k for e in edges for k in e})
    max_off = max((abs(off) for off, _ in edges.values()), default=0)
    total = sum(abs(off) for off, _ in edges.values())
    return float(min(total, n_nodes * max_off) + 1)


@contextmanager
def _quiet_native_output():
    """Silence diagnostics that the native solver writes straight to the
    process's stdout/stderr file descriptors."""
    try:
        fds = [sys.stdout.fileno(), sys.stderr.fileno()]
    except (AttributeError, ValueError, io.UnsupportedOperation):
        yield
        return
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(fd) for fd in fds]
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        for fd in fds:
            os.dup2(devnull, fd)
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        for fd, keep in zip(fds, saved):
            os.dup2(keep, fd)
            os.close(keep)
        os.close(devnull)


class _LexProgram:
    """Integer program with a block of binary selection variables, solved
    lexicographically and canonicalised."""

    def __init__(self, A: csr_matrix, ub: np.ndarray, lb_vars: np.ndarray,
                 ub_vars: np.ndarray, select: Sequence[int]):
        self.n_vars = A.shape[1]
        self.rows = [A]
        self.row_lo = [np.full(A.shape[0], -np.inf)]
        self.row_hi = [ub]
        self.lb = lb_vars.astype(float).copy()
        self.hb = ub_vars.astype(float).copy()
        self.select = list(select)
        self.ok = True
        self.messages: List[str] = []

    def _solve(self, cost: np.ndarray):
        A = vstack(self.rows).tocsr()
        with _quiet_native_output():
            res = milp(
                cost,
                constraints=LinearConstraint(A, np.concatenate(self.row_lo),
                                             np.concatenate(self.row_hi)),
                integrality=np.ones(self.n_vars),
                bounds=Bounds(self.lb, self.hb),
                options=SOLVER_OPTIONS,
            )
        return res

    def minimise_and_fix(self, weights: np.ndarray):
        """Minimise ``weights . select`` and pin the optimum as an equality."""
        cost = np.zeros(self.n_vars)
        cost[self.select] = weights
        res = self._solve(cost)
        self.messages.append(res.message)
        if res.status != 0 or res.x is None:
            self.ok = False
            return None
        value = float(np.round(np.dot(weights, np.round(res.x[self.select]))))
        row = np.zeros((1, self.n_vars))
        row[0, self.select] = weights
        self.rows.append(csr_matrix(row))
        self.row_lo.append(np.array([value]))
        self.row_hi.append(np.array([value]))
        return res.x, value

    def value_range(self, var: int) -> Optional[Tuple[int, int]]:
        """Smallest and largest value ``var`` takes over the feasible set.

        Both directions are solved as optimisations of a feasible program
        rather than as feasibility tests of a pinned one: the program already
        contains a known solution, so an infeasible sub-problem — which some
        solver builds report as a generic error rather than as infeasibility —
        never has to be distinguished from a failure. None if a solve did not
        reach optimality.
        """
        out: List[int] = []
        for sign in (1.0, -1.0):
            cost = np.zeros(self.n_vars)
            cost[var] = sign
            res = self._solve(cost)
            if res.status != 0 or res.x is None:
                self.messages.append("membership test failed: " + str(res.message))
                return None
            out.append(int(round(float(res.x[var]))))
        return out[0], out[1]

    def canonicalise(self, x: np.ndarray) -> np.ndarray:
        """Walk candidates in order; keep each out of the selection whenever
        an optimal solution without it exists."""
        current = np.round(x).astype(int)
        zero_cost = np.zeros(self.n_vars)
        for var in self.select:
            if current[var] == 0:
                self.hb[var] = 0.0
                continue
            saved = self.hb[var]
            self.hb[var] = 0.0
            res = self._solve(zero_cost)
            if res.status == 0 and res.x is not None:
                current = np.round(res.x).astype(int)
            elif res.status == 2:  # Proven infeasible, not a timeout or solver error.
                self.hb[var] = saved
                self.lb[var] = 1.0
            else:
                self.hb[var] = saved
                self.ok = False
                self.messages.append("canonicalisation failed: " + str(res.message))
                break
        return current


def _edge_rows(edge_list, nidx, n, extra_cols):
    """Big-M rows: |u_b - u_a - offset| <= M * (sum of the edge's release vars)."""
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    ub: List[float] = []
    bound = _potential_bound(dict(edge_list))
    big_m = 2.0 * bound + max(abs(off) for _, (off, _) in edge_list) + 1.0
    for i, ((a, b), (off, _)) in enumerate(edge_list):
        ia, ib = nidx[a], nidx[b]
        release = extra_cols(i, a, b)
        for sign, rhs in ((1.0, float(off)), (-1.0, float(-off))):
            r = len(ub)
            rows += [r, r] + [r] * len(release)
            cols += [ib, ia] + release
            vals += [sign, -sign] + [-big_m] * len(release)
            ub.append(rhs)
    return rows, cols, vals, np.array(ub), bound


def _quarantine_program(edges: Edges):
    """Build the quarantine program with both objective optima pinned.

    Returns ``(prog, node_list, n, stage)``; ``stage`` is ``None`` when a stage
    did not solve. Every optimal quarantine of ``edges`` — and no other node
    set — satisfies the returned program, so membership questions about the
    set of optima can be answered by feasibility alone.
    """
    node_list = sorted({k for e in edges for k in e})
    nidx = {k: i for i, k in enumerate(node_list)}
    edge_list = sorted(edges.items())
    n, m = len(node_list), len(edge_list)
    # y_e = OR(x_a, x_b): every lost edge is charged exactly once.
    rows, cols, vals, ub, bound = _edge_rows(
        edge_list, nidx, n, lambda i, a, b: [2 * n + i]
    )
    upper = list(ub)
    for i, ((a, b), _) in enumerate(edge_list):
        xa, xb, y = n + nidx[a], n + nidx[b], 2 * n + i
        for cc, vv in (([xa, y], [1., -1.]),
                       ([xb, y], [1., -1.]),
                       ([y, xa, xb], [1., -1., -1.])):
            rows.extend([len(upper)] * len(cc))
            cols.extend(cc)
            vals.extend(vv)
            upper.append(0.)
    A = coo_matrix((vals, (rows, cols)), shape=(len(upper), 2 * n + m)).tocsr()
    lb = np.concatenate([np.full(n, -bound), np.zeros(n + m)])
    hb = np.concatenate([np.full(n, bound), np.ones(n + m)])
    prog = _LexProgram(A, np.array(upper), lb, hb, list(range(n, 2 * n + m)))
    stage = prog.minimise_and_fix(np.concatenate([np.ones(n), np.zeros(m)]))
    count = stage[1] if stage else None
    if stage:
        support = np.array([sup for _, (_, sup) in edge_list], dtype=float)
        stage = prog.minimise_and_fix(np.concatenate([np.zeros(n), support]))
    # Once both objectives are fixed, node order defines the canonical answer.
    prog.select = list(range(n, 2 * n))
    return prog, node_list, n, count, stage


def solve_min_quarantine(edges: Edges) -> RepairResult:
    """Fewest collections to remove so the remaining edges are consistent;
    ties broken by least link support lost, then by node order."""
    _validate_edges(edges)
    if edges_consistent(edges):
        return RepairResult([], [], 0.0, True, "already consistent", 0)
    edge_list = sorted(edges.items())
    prog, node_list, n, count, stage = _quarantine_program(edges)
    removed: List[NodeKey] = []
    if stage:
        final = prog.canonicalise(stage[0])
        removed = [node_list[i] for i in range(n) if final[n + i] == 1]
    hidden = set(removed)
    verified = bool(removed) and edges_consistent(
        {e: v for e, v in edges.items() if e[0] not in hidden and e[1] not in hidden}
    )
    return RepairResult(
        removed=removed,
        cut_edges=[],
        objective=count,
        optimal=prog.ok and verified,
        status="optimal" if prog.ok and verified else "; ".join(prog.messages),
        support_removed=int(sum(sup for (a, b), (_, sup) in edge_list
                                if a in hidden or b in hidden)),
    )


def quarantine_membership(edges: Edges,
                          targets: Optional[Sequence[NodeKey]] = None
                          ) -> Dict[NodeKey, str]:
    """Classify collections by their membership in the set of optimal repairs.

    For each target collection the answer is one of

    * ``"certain"`` — it is removed by *every* optimal quarantine, so the
      evidence names it without relying on any tie-break convention;
    * ``"possible"`` — some optimal quarantine removes it and some does not;
    * ``"excluded"`` — no optimal quarantine removes it;
    * ``"consistent"`` — the graph has no contradiction to repair;
    * ``"unknown"`` — the solver did not settle the question.

    Each answer costs two solves on the program whose feasible set is exactly
    the set of optimal quarantines, so it is independent of which optimum
    ``solve_min_quarantine`` happens to return, and of any tie-break.
    """
    _validate_edges(edges)
    node_list_all = sorted({k for e in edges for k in e})
    wanted = list(node_list_all if targets is None else targets)
    if edges_consistent(edges):
        return {k: "consistent" for k in wanted}
    prog, node_list, n, _count, stage = _quarantine_program(edges)
    nidx = {k: i for i, k in enumerate(node_list)}
    out: Dict[NodeKey, str] = {}
    for key in wanted:
        if stage is None or key not in nidx:
            out[key] = "excluded" if (stage is not None) else "unknown"
            continue
        span = prog.value_range(n + nidx[key])
        if span is None:
            out[key] = "unknown"
        elif span[0] == 1:
            out[key] = "certain"
        elif span[1] == 0:
            out[key] = "excluded"
        else:
            out[key] = "possible"
    return out


def solve_min_edge_cut(edges: Edges) -> RepairResult:
    """Least total link support to cut so the remaining edges are consistent;
    ties broken by edge order."""
    _validate_edges(edges)
    if edges_consistent(edges):
        return RepairResult([], [], 0.0, True, "already consistent", 0)
    node_list = sorted({k for e in edges for k in e})
    nidx = {k: i for i, k in enumerate(node_list)}
    edge_list = sorted(edges.items())
    n, m = len(node_list), len(edge_list)
    rows, cols, vals, ub, bound = _edge_rows(
        edge_list, nidx, n, lambda i, a, b: [n + i]
    )
    A = coo_matrix((vals, (rows, cols)), shape=(len(ub), n + m)).tocsr()
    lb = np.concatenate([np.full(n, -bound), np.zeros(m)])
    hb = np.concatenate([np.full(n, bound), np.ones(m)])
    prog = _LexProgram(A, ub, lb, hb, list(range(n, n + m)))

    support = np.array([float(sup) for _, (_, sup) in edge_list])
    stage = prog.minimise_and_fix(support)
    cut: List[Tuple[NodeKey, NodeKey]] = []
    value = None
    if stage:
        value = stage[1]
        final = prog.canonicalise(stage[0])
        cut = [edge_list[i][0] for i in range(m) if final[n + i] == 1]
    cut_set = set(cut)
    verified = bool(cut) and edges_consistent(
        {e: v for e, v in edges.items() if e not in cut_set}
    )
    return RepairResult(
        removed=[],
        cut_edges=cut,
        objective=value,
        optimal=prog.ok and verified,
        status="optimal" if prog.ok and verified else "; ".join(prog.messages),
        support_removed=int(sum(edges[e][1] for e in cut)),
    )


def brute_force_min_quarantine(edges: Edges, max_size: int = 4) -> Optional[List[NodeKey]]:
    """Exhaustive reference for tests: smallest node subset (up to max_size)
    whose removal leaves a consistent edge set; None if none exists."""
    if edges_consistent(edges):
        return []
    nodes = sorted({k for e in edges for k in e})
    for size in range(1, max_size + 1):
        for subset in itertools.combinations(nodes, size):
            hidden = set(subset)
            rest = {
                e: v for e, v in edges.items() if e[0] not in hidden and e[1] not in hidden
            }
            if edges_consistent(rest):
                return list(subset)
    return None


# Retained for callers that used the pre-0.2.1 private name.
_edges_consistent = edges_consistent
