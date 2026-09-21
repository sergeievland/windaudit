import math

import numpy as np
import pytest

from windaudit.frame import (
    KIND_ABSOLUTE,
    KIND_RELATIVE,
    KIND_SAME,
    Umbilicus,
    build_absolute_nodes,
    build_chain_node,
    build_nodes,
    cut_crossing,
    label_diff,
    wrap_angle,
)
from tests.conftest import make_collection, wrap_xyz


def test_cylindrical_frame_matches_upstream_convention(umbilicus):
    r, theta = umbilicus.cylindrical([10.0, 0.0, 5.0])
    assert r == pytest.approx(10.0)
    assert theta == pytest.approx(0.0)
    r, theta = umbilicus.cylindrical([0.0, -7.0, 5.0])
    assert r == pytest.approx(7.0)
    assert theta == pytest.approx(1.5 * math.pi)


def test_umbilicus_interpolates_axis():
    umb = Umbilicus(np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 100.0]]))
    r, _ = umb.cylindrical([60.0, 0.0, 50.0])
    assert r == pytest.approx(10.0)


def test_wrap_angle():
    assert wrap_angle(3 * math.pi) == pytest.approx(math.pi)
    assert wrap_angle(-3 * math.pi) == pytest.approx(math.pi)
    assert wrap_angle(0.3) == pytest.approx(0.3)


def test_cut_crossing_direction():
    assert cut_crossing(2 * math.pi - 0.1, 0.1) == 1
    assert cut_crossing(0.1, 2 * math.pi - 0.1) == -1
    assert cut_crossing(1.0, 1.2) == 0
    assert cut_crossing(math.pi - 0.1, math.pi + 0.1) == 0


def test_trace_across_cut_changes_branch_label(umbilicus):
    """A same-winding trace walking forward through theta = 0 gains one
    branch winding (CW) and keeps label differences equal to zero."""
    thetas = [2 * math.pi - 0.2, 2 * math.pi - 0.05, 0.05, 0.2]
    rows = [tuple(wrap_xyz(3, t, 50.0)) + (None,) for t in thetas]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_SAME, sense=1)
    assert node.seam_safe
    assert [q.label for q in node.points] == [0, 0, 1, 1]
    for a in node.points:
        for b in node.points:
            assert label_diff(a, b, 1) == 0


def test_acw_sense_flips_crossing(umbilicus):
    thetas = [2 * math.pi - 0.1, 0.1]
    rows = [(100 * math.cos(t), 100 * math.sin(t), 50.0, None) for t in thetas]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_SAME, sense=-1)
    assert [q.label for q in node.points] == [0, -1]


def test_half_turn_without_cut_keeps_label(umbilicus):
    """Accumulated angle alone never changes the integer winding: only the
    branch cut does."""
    thetas = [0.2 + k * 0.5 for k in range(7)]  # 0.2 .. 3.2 rad, no cut
    rows = [tuple(wrap_xyz(2, t, 50.0)) + (None,) for t in thetas]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_SAME, sense=1)
    assert {q.label for q in node.points} == {0}


def test_relative_ladder_labels(umbilicus):
    rows = [tuple(wrap_xyz(w, 1.0, 50.0)) + (w,) for w in (4, 5, 6)]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_RELATIVE, sense=1)
    assert [q.label for q in node.points] == [4, 5, 6]


def test_seam_unsafe_flag(umbilicus):
    rows = [(100, 0, 0, 1.0), (-100, 10, 0, 2.0)]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_RELATIVE, sense=1)
    assert not node.seam_safe


def test_relative_node_skips_unannotated_points(umbilicus):
    rows = [(100, 0, 0, 1.0), (120, 0, 0, None), (140, 0, 0, 2.0)]
    node = build_chain_node(make_collection(1, rows), umbilicus, KIND_RELATIVE, sense=1)
    assert [q.pid for q in node.points] == [1, 3]


def test_non_integer_winding_is_rejected(umbilicus):
    rows = [(100, 0, 0, 1.5)]
    with pytest.raises(ValueError):
        build_chain_node(make_collection(1, rows), umbilicus, KIND_RELATIVE, sense=1)


def test_absolute_clusters_and_cut_guard(umbilicus):
    rows = [tuple(wrap_xyz(10, t, 50.0)) + (10,) for t in (1.0, 1.02, 1.04)]
    rows += [tuple(wrap_xyz(11, t, 50.0)) + (11,) for t in (3.0, 3.02)]
    rows += [tuple(wrap_xyz(12, 0.05, 50.0)) + (12,)]
    nodes = build_absolute_nodes(make_collection(1, rows, absolute=True), umbilicus)
    sizes = sorted(n.n for n in nodes)
    assert sizes == [1, 2, 3]
    near_cut = [n for n in nodes if not n.seam_safe]
    assert len(near_cut) == 1 and near_cut[0].n == 1
    assert all(n.kind == KIND_ABSOLUTE for n in nodes)


def test_content_classification_is_enforced(umbilicus):
    annotated = make_collection(1, [(100, 0, 0, 1.0)])
    with pytest.raises(ValueError):
        build_nodes({}, {1: annotated}, {}, umbilicus)


def test_empty_collections_are_dropped(umbilicus):
    nodes = build_nodes({}, {1: make_collection(1, [])}, {}, umbilicus)
    assert nodes == {}
