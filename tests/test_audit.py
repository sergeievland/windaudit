from windaudit.audit import (
    order_audit,
    same_collection_step_audit,
    shuffled_order_control,
)
from windaudit.frame import build_nodes
from tests.conftest import ladder, make_collection, wrap_trace, wrap_xyz


def _nodes(umbilicus, relative=None, same=None, absolute=None):
    return build_nodes(relative or {}, same or {}, absolute or {}, umbilicus)


def test_clean_ladder_passes(umbilicus, config):
    nodes = _nodes(umbilicus, relative={1: ladder(1, 1.0, 500.0, range(1, 12))})
    result = order_audit(nodes, config)
    assert result.n_pairs_tested > 0
    assert result.n_violations == 0


def test_swapped_labels_are_caught(umbilicus, config):
    bad = ladder(1, 1.0, 500.0, range(1, 12))
    bad.points[3].wind_a, bad.points[8].wind_a = (
        bad.points[8].wind_a,
        bad.points[3].wind_a,
    )
    nodes = _nodes(umbilicus, relative={1: bad})
    result = order_audit(nodes, config)
    assert result.n_violations > 0
    assert ("relative", 1) in result.flagged_nodes


def test_tilted_trace_does_not_false_flag(umbilicus, config):
    nodes = _nodes(
        umbilicus,
        same={1: wrap_trace(1, wrap=3, theta=1.0,
                            zs=range(400, 900, 20), tilt=0.4)},
    )
    step = same_collection_step_audit(nodes, config)
    assert step.n_steps_tested > 0
    assert step.n_violations == 0


def test_wrap_jump_in_trace_is_caught(umbilicus, config):
    rows = ([tuple(wrap_xyz(3, 1.0, z)) + (None,) for z in range(400, 500, 20)]
            + [tuple(wrap_xyz(4, 1.0, z)) + (None,) for z in range(500, 600, 20)])
    nodes = _nodes(umbilicus, same={1: make_collection(1, rows)})
    step = same_collection_step_audit(nodes, config)
    assert step.n_violations >= 1
    assert ("same", 1) in step.flagged_nodes


def test_shuffled_control_has_power(umbilicus, config):
    nodes = _nodes(
        umbilicus,
        relative={i: ladder(i, 1.0, 400.0 + 60 * i, range(1, 14))
                  for i in range(1, 7)},
    )
    control = shuffled_order_control(nodes, config, n_rounds=5, seed=1)
    assert control["real_violation_rate"] == 0.0
    assert control["shuffled_mean_rate"] > 0.2
