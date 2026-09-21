import numpy as np
import pytest

from windaudit.repair import (
    brute_force_min_quarantine,
    solve_min_edge_cut,
    solve_min_quarantine,
    _edges_consistent,
)


def _key(i):
    return ("n", i)


def test_consistent_graph_needs_nothing():
    edges = {(_key(1), _key(2)): (1, 3), (_key(2), _key(3)): (2, 3),
             (_key(1), _key(3)): (3, 3)}
    assert _edges_consistent(edges)
    assert solve_min_quarantine(edges).removed == []
    assert solve_min_edge_cut(edges).cut_edges == []


def test_single_bad_edge_is_cut():
    edges = {(_key(1), _key(2)): (1, 5), (_key(2), _key(3)): (2, 5),
             (_key(1), _key(3)): (4, 1)}  # should be 3
    cut = solve_min_edge_cut(edges)
    assert cut.optimal
    assert cut.cut_edges == [(_key(1), _key(3))]
    q = solve_min_quarantine(edges)
    assert q.optimal
    assert len(q.removed) == 1


def test_edge_cut_prefers_low_support():
    edges = {(_key(1), _key(2)): (1, 9), (_key(2), _key(3)): (2, 9),
             (_key(1), _key(3)): (4, 1)}
    cut = solve_min_edge_cut(edges)
    assert cut.cut_edges == [(_key(1), _key(3))]
    assert cut.support_removed == 1


def test_quarantine_hub_explains_many_cycles():
    """One node whose every edge is shifted by +1: quarantining it is cheaper
    than cutting all its edges."""
    edges = {}
    for i in range(1, 5):
        for j in range(i + 1, 5):
            edges[(_key(i), _key(j))] = (0, 2)
    hub = _key(9)
    for i in range(1, 5):
        edges[(_key(i), hub)] = (i, 2)  # inconsistent by design
    q = solve_min_quarantine(edges)
    assert q.optimal
    assert q.removed == [hub]


@pytest.mark.parametrize("seed", range(8))
def test_quarantine_matches_brute_force(seed):
    rng = np.random.default_rng(seed)
    n = 7
    truth = {i: int(rng.integers(-5, 6)) for i in range(n)}
    edges = {}
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.55:
                edges[(_key(i), _key(j))] = (truth[j] - truth[i],
                                             int(rng.integers(1, 4)))
    # corrupt up to two edges
    keys = sorted(edges)
    for k in rng.choice(len(keys), size=min(2, len(keys)), replace=False):
        off, sup = edges[keys[k]]
        edges[keys[k]] = (off + int(rng.integers(1, 3)), sup)
    got = solve_min_quarantine(edges)
    ref = brute_force_min_quarantine(edges, max_size=4)
    assert got.optimal
    assert ref is not None
    assert len(got.removed) == len(ref)
    # verify the solver's answer really restores consistency
    hidden = set(got.removed)
    rest = {e: v for e, v in edges.items()
            if e[0] not in hidden and e[1] not in hidden}
    assert _edges_consistent(rest)


def test_large_offsets_do_not_break_bounds():
    edges = {(_key(1), _key(2)): (500, 1), (_key(2), _key(3)): (-499, 1),
             (_key(1), _key(3)): (0, 1)}  # inconsistent by 1
    q = solve_min_quarantine(edges)
    assert q.optimal
    assert len(q.removed) == 1
    cut = solve_min_edge_cut(edges)
    assert cut.optimal
    assert len(cut.cut_edges) == 1


def test_quarantine_ties_resolve_canonically():
    """Every node of an inconsistent triangle with equal support is an
    optimal quarantine; the canonical choice is the last node in order."""
    edges = {(_key(1), _key(2)): (1, 2), (_key(2), _key(3)): (1, 2),
             (_key(1), _key(3)): (5, 2)}
    first = solve_min_quarantine(edges)
    assert first.optimal
    assert first.removed == [_key(3)]
    reordered = dict(reversed(list(edges.items())))
    assert solve_min_quarantine(reordered).removed == first.removed


def test_quarantine_prefers_least_support_before_order():
    edges = {(_key(1), _key(2)): (1, 1), (_key(2), _key(3)): (1, 9),
             (_key(1), _key(3)): (5, 2)}
    q = solve_min_quarantine(edges)
    assert q.optimal
    assert q.removed == [_key(1)]
    assert q.support_removed == 3


def test_edge_cut_ties_resolve_canonically():
    edges = {(_key(1), _key(2)): (1, 3), (_key(2), _key(3)): (1, 3),
             (_key(1), _key(3)): (5, 3)}
    cut = solve_min_edge_cut(edges)
    assert cut.optimal
    assert cut.cut_edges == [(_key(2), _key(3))]
