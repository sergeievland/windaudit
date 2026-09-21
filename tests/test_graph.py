import math

from windaudit.frame import build_nodes
from windaudit.graph import WindingGraph, constant_label_nodes
from windaudit.matching import build_links
from tests.conftest import ladder, make_collection, wrap_trace


def _nodes(umbilicus, relative=None, same=None, absolute=None):
    return build_nodes(relative or {}, same or {}, absolute or {}, umbilicus)


def _consistent_scene(umbilicus):
    """Two ladders 300 z apart, bridged by three wrap traces: real cycles,
    all consistent."""
    rel = {
        1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 8)),
        2: ladder(2, theta=1.0, z=800.0, wraps=range(1, 8)),
    }
    same = {
        w: wrap_trace(w, wrap=w, theta=1.0, zs=range(470, 840, 20))
        for w in (3, 4, 5)
    }
    return rel, same


def test_consistent_scene_has_zero_holonomy(umbilicus, config):
    rel, same = _consistent_scene(umbilicus)
    match = build_links(_nodes(umbilicus, rel, same), config)
    graph = WindingGraph(match)
    stats = graph.stats()
    assert stats["n_effective_cycles"] >= 2
    assert stats["n_nonzero_holonomy"] == 0


def test_planted_skip_creates_holonomy(umbilicus, config):
    rel, same = _consistent_scene(umbilicus)
    # Ladder 2 skips a wrap: annotations 5.. shift by +1.
    bad = ladder(2, theta=1.0, z=800.0, wraps=range(1, 8))
    for pid, p in bad.points.items():
        if p.wind_a >= 5:
            p.wind_a += 1
    rel[2] = bad
    match = build_links(_nodes(umbilicus, rel, same), config)
    graph = WindingGraph(match)
    cycles = graph.cycles()
    assert any(c.holonomy != 0 for c in cycles)


def test_tautological_cycles_are_classified(umbilicus, config):
    """Two wrap traces on one wrap linked in two places close a cycle whose
    edges are all zero offsets between unlabeled nodes: it cannot fail and
    must not count as effective."""
    same = {
        1: wrap_trace(1, wrap=3, theta=1.0, zs=range(400, 700, 20)),
        2: wrap_trace(2, wrap=3, theta=1.0, zs=range(400, 700, 20)),
        3: wrap_trace(3, wrap=3, theta=1.0, zs=range(400, 700, 20)),
    }
    nodes = _nodes(umbilicus, same=same)
    match = build_links(nodes, config)
    graph = WindingGraph(match, constant_label_nodes(nodes))
    stats = graph.stats()
    assert stats["n_cycles"] >= 1
    assert stats["n_effective_cycles"] == 0
    assert stats["n_nonzero_holonomy"] == 0


def test_potentials_reproduce_offsets(umbilicus, config):
    rel, same = _consistent_scene(umbilicus)
    match = build_links(_nodes(umbilicus, rel, same), config)
    graph = WindingGraph(match)
    for (a, b), (off, _) in match.edges.items():
        assert graph.potential[b] - graph.potential[a] == off


def test_ring_claimed_as_one_sheet_is_caught(umbilicus, config):
    """Overlapping same-winding traces that close a ring at constant radius
    claim a sheet meets itself one full turn later: the branch-corrected
    offsets around the ring sum to one winding."""
    same = {}
    for k in range(12):
        centre = k * 30 + 15
        thetas = [math.radians(centre + d) for d in range(-20, 21, 2)]
        rows = [(300 * math.cos(t), 300 * math.sin(t), 500.0, None) for t in thetas]
        same[k + 1] = make_collection(k + 1, rows)
    nodes = _nodes(umbilicus, same=same)
    match = build_links(nodes, config)
    graph = WindingGraph(match, constant_label_nodes(nodes))
    failing = [c for c in graph.cycles() if c.holonomy != 0]
    assert len(failing) == 1
    assert abs(failing[0].holonomy) == 1
    assert failing[0].effective
