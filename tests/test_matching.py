import math

from windaudit.frame import build_nodes
from windaudit.matching import build_links, point_slopes
from windaudit.frame import build_chain_node, KIND_RELATIVE, KIND_SAME
from tests.conftest import ladder, make_collection, wrap_trace, wrap_xyz


def _nodes(umbilicus, relative=None, same=None, absolute=None):
    return build_nodes(relative or {}, same or {}, absolute or {}, umbilicus)


def test_ladder_and_trace_link_with_correct_offset(umbilicus, config):
    rel = {1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 8))}
    same = {1: wrap_trace(1, wrap=4, theta=1.0, zs=[470, 490, 510, 530])}
    match = build_links(_nodes(umbilicus, rel, same), config)
    edge = match.edges[(("relative", 1), ("same", 1))]
    # offset = gauge_same - gauge_rel = label_rel - label_same = 4 - 0
    assert edge[0] == 4
    assert edge[1] >= 1
    assert not match.conflicts


def test_no_link_across_wraps(umbilicus, config):
    rel = {1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 8))}
    same = {1: wrap_trace(1, wrap=4, theta=1.0, zs=[470, 490, 510, 530])}
    match = build_links(_nodes(umbilicus, rel, same), config)
    offsets = {match.edges[e][0] for e in match.edges}
    assert offsets == {4}


def test_tilt_does_not_create_cross_wrap_links(umbilicus, config):
    """With wraps drifting 0.4 voxel of radius per voxel of z, a tilt-blind
    matcher would pair a vertical trace with the neighbouring wrap of a
    single-slice ladder; the projection must keep the offset unique."""
    tilt = 0.4
    rel = {1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 10), tilt=tilt)}
    same = {
        1: wrap_trace(1, wrap=4, theta=1.0,
                      zs=[440, 460, 480, 500, 520, 540, 560], tilt=tilt)
    }
    match = build_links(_nodes(umbilicus, rel, same), config)
    pair = (("relative", 1), ("same", 1))
    assert pair in match.edges
    assert match.edges[pair][0] == 4
    assert not match.conflicts


def test_crowding_veto_blocks_compressed_regions(umbilicus, config):
    """When wrap spacing collapses below the crowding floor, no links are
    emitted there."""
    spacing = 0.5 * config.crowd_min_spacing
    rel_rows = [tuple(wrap_xyz(w, 1.0, 500.0, spacing=spacing)) + (w,)
                for w in range(1, 8)]
    rel = {1: make_collection(1, rel_rows)}
    same_rows = [tuple(wrap_xyz(4, 1.0, z, spacing=spacing)) + (None,)
                 for z in [480, 500, 520]]
    same = {1: make_collection(1, same_rows)}
    match = build_links(_nodes(umbilicus, rel, same), config)
    assert not match.edges
    assert match.n_vetoed_cells > 0


def test_conflicting_offsets_reported_not_merged(umbilicus, config):
    """A trace that sits on wrap 4 where it crosses one slice of a two-slice
    ladder and on wrap 5 where it crosses the other must surface as an edge
    conflict, never be merged into one offset."""
    rows = [tuple(wrap_xyz(w, 1.0, z)) + (w,) for z in (400.0, 600.0) for w in range(1, 10)]
    rel = {1: make_collection(1, rows)}
    trace = ([tuple(wrap_xyz(4, 1.0, z)) + (None,) for z in (380, 400, 420)]
             + [tuple(wrap_xyz(5, 1.0, z)) + (None,) for z in (580, 600, 620)])
    same = {1: make_collection(1, trace)}
    match = build_links(_nodes(umbilicus, rel, same), config)
    assert len(match.conflicts) == 1
    assert match.conflicts[0].offsets == [4, 5]
    assert match.edges == {}


def test_slopes_are_per_label(umbilicus):
    """A diagonal chain stepping one wrap per z step must not mix wraps into
    one slope estimate: single-point labels give no slope."""
    rows = [tuple(wrap_xyz(w, 1.0, 500.0 + 20 * w)) + (w,) for w in range(1, 6)]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_RELATIVE, sense=1)
    assert point_slopes(node, 1) == [None] * 5


def test_slopes_follow_the_wrap(umbilicus):
    tilt = 0.3
    node = build_chain_node(wrap_trace(1, wrap=2, theta=1.0,
                                       zs=[400, 420, 440, 460], tilt=tilt),
                            umbilicus, KIND_SAME, sense=1)
    slopes = point_slopes(node, 1)
    assert all(s is not None for s in slopes)
    for s in slopes:
        assert abs(s - tilt) < 0.02


def test_seam_unsafe_nodes_are_excluded(umbilicus, config):
    rel = {1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 6))}
    rows = [(100, 0, 500, None), (-100, 10, 500, None), (100, 5, 500, None)]
    same = {1: make_collection(1, rows)}
    nodes = _nodes(umbilicus, rel, same)
    assert not nodes[("same", 1)].seam_safe
    match = build_links(nodes, config)
    assert all(("same", 1) not in pair for pair in match.edges)


def test_link_across_the_cut_has_branch_corrected_offset(umbilicus, config):
    """A ladder just before theta = 0 and a trace just after it on the same
    wrap: the trace's integer winding is one higher, and the offset must
    absorb exactly that."""
    eps = 0.01
    rel = {1: ladder(1, theta=2 * math.pi - eps, z=500.0, wraps=range(1, 8))}
    rows = [tuple(wrap_xyz(5, eps, z)) + (None,) for z in (480, 500, 520)]
    same = {1: make_collection(1, rows)}
    match = build_links(_nodes(umbilicus, rel, same), config)
    pair = (("relative", 1), ("same", 1))
    assert pair in match.edges
    # gauge_same = W(trace) - 0 = 5; gauge_rel = W - label = 0; offset = 5.
    assert match.edges[pair][0] == 5
    assert match.n_cut_links >= 1


def test_single_link_is_not_an_edge(umbilicus, config):
    rel = {1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 8))}
    rows = [tuple(wrap_xyz(4, 1.0, 500.0)) + (None,),
            tuple(wrap_xyz(4, 1.0, 700.0)) + (None,),
            tuple(wrap_xyz(4, 1.0, 900.0)) + (None,),
            tuple(wrap_xyz(4, 1.0, 1100.0)) + (None,)]
    same = {1: make_collection(1, rows)}
    match = build_links(_nodes(umbilicus, rel, same), config)
    assert match.edges == {}
    assert match.n_unsupported_links == 1
