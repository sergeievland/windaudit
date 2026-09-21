"""Independent exhaustive objectives and failure-path regression tests."""
import itertools
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from windaudit.repair import (_LexProgram, _edges_consistent, quarantine_membership,
                              solve_min_quarantine, solve_min_edge_cut, RepairResult)
from windaudit.report import certified_links


def oracle(edges, quarantine):
    candidates = sorted({v for e in edges for v in e}) if quarantine else sorted(edges)
    feasible = []
    for bits in itertools.product((0, 1), repeat=len(candidates)):
        chosen = {k for k, bit in zip(candidates, bits) if bit}
        lost = {e for e in edges if (bool(set(e) & chosen) if quarantine else e in chosen)}
        if _edges_consistent({e: v for e, v in edges.items() if e not in lost}):
            support = sum(edges[e][1] for e in lost)
            score = (len(chosen), support, bits) if quarantine else (support, bits)
            feasible.append((score, chosen))
    return min(feasible, key=lambda x: x[0])[1]


@pytest.mark.parametrize('seed', range(24))
def test_both_full_lexicographic_objectives_against_exhaustive(seed):
    rng = np.random.default_rng(seed)
    nodes = [('relative', i) for i in range(5)]
    edges = {(a, b): (int(rng.integers(-2, 3)), int(rng.integers(0, 9)))
             for a, b in itertools.combinations(nodes, 2) if rng.random() < .8}
    for quarantine, solve in [(True, solve_min_quarantine), (False, solve_min_edge_cut)]:
        result = solve(edges)
        assert result.optimal
        assert set(result.removed if quarantine else result.cut_edges) == oracle(edges, quarantine)


@pytest.mark.parametrize('status', [1, 3, 4])
def test_canonical_solver_failure_cannot_claim_optimal(monkeypatch, status):
    prog = _LexProgram(csr_matrix([[1.]]), np.array([1.]), np.zeros(1), np.ones(1), [0])
    monkeypatch.setattr(prog, '_solve', lambda _: SimpleNamespace(status=status, x=None, message='injected failure'))
    prog.canonicalise(np.ones(1))
    assert not prog.ok


def test_export_refuses_unproven_repair(tmp_path):
    path = tmp_path/'links.json'
    with pytest.raises(RuntimeError, match='not proven optimal'):
        certified_links(None, RepairResult([], [], None, False, 'timeout', 0), str(path))
    assert not path.exists()


def test_export_independently_checks_remainder(tmp_path):
    a, b, c = [('relative', i) for i in range(3)]
    match = SimpleNamespace(edges={(a,b):(0,1), (b,c):(0,1), (a,c):(1,1)})
    with pytest.raises(RuntimeError, match='inconsistent'):
        certified_links(match, RepairResult([], [], 0, True, 'optimal', 0), str(tmp_path/'links.json'))


@pytest.mark.parametrize('value', [(0,-1), (0,1.5), (float('nan'),1), (1.5,1)])
def test_bad_graph_values_rejected(value):
    for solve in [solve_min_quarantine, solve_min_edge_cut]:
        with pytest.raises(ValueError):
            solve({(('r',1),('r',2)):value})


def membership_oracle(edges):
    """Every optimal quarantine, by exhaustive enumeration."""
    nodes = sorted({v for e in edges for v in e})
    best, optima = None, []
    for bits in itertools.product((0, 1), repeat=len(nodes)):
        chosen = {k for k, bit in zip(nodes, bits) if bit}
        lost = {e for e in edges if set(e) & chosen}
        if _edges_consistent({e: v for e, v in edges.items() if e not in lost}):
            score = (len(chosen), sum(edges[e][1] for e in lost))
            if best is None or score < best:
                best, optima = score, [chosen]
            elif score == best:
                optima.append(chosen)
    return optima


@pytest.mark.parametrize('seed', range(24))
def test_membership_matches_exhaustive_set_of_optima(seed):
    rng = np.random.default_rng(seed)
    nodes = [('relative', i) for i in range(5)]
    edges = {(a, b): (int(rng.integers(-2, 3)), int(rng.integers(0, 9)))
             for a, b in itertools.combinations(nodes, 2) if rng.random() < .8}
    if not edges or _edges_consistent(edges):
        pytest.skip('graph has nothing to repair')
    optima = membership_oracle(edges)
    got = quarantine_membership(edges)
    for key in sorted({v for e in edges for v in e}):
        expected = ('certain' if all(key in o for o in optima)
                    else 'possible' if any(key in o for o in optima) else 'excluded')
        assert got[key] == expected, (key, got[key], expected)
    # The canonical repair must itself be one of the optima.
    assert set(solve_min_quarantine(edges).removed) in [set(o) for o in optima]


def test_membership_on_consistent_graph_reports_nothing_to_repair():
    a, b, c = [('relative', i) for i in range(3)]
    edges = {(a, b): (1, 2), (b, c): (1, 2), (a, c): (2, 2)}
    assert set(quarantine_membership(edges).values()) == {'consistent'}
