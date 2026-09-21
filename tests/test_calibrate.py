import pytest

from windaudit.calibrate import calibrate_from_absolute_anchors
from windaudit.frame import build_nodes
from tests.conftest import SPACING, make_collection, wrap_xyz


def _anchor_collection(cid, wraps, theta=1.0, z=500.0):
    rows = [tuple(wrap_xyz(w, theta, z)) + (w,) for w in wraps]
    return make_collection(cid, rows, absolute=True)


def test_calibration_recovers_spacing(umbilicus):
    absolute = {
        1: _anchor_collection(1, range(1, 12)),
        2: _anchor_collection(2, range(1, 12), z=600.0),
    }
    nodes = build_nodes({}, {}, absolute, umbilicus)
    cfg = calibrate_from_absolute_anchors(nodes)
    assert cfg.wrap_spacing_median == pytest.approx(SPACING, rel=0.05)
    assert cfg.r_same == pytest.approx(0.35 * cfg.wrap_spacing_median)
    assert cfg.crowd_min_spacing == pytest.approx(0.70 * cfg.wrap_spacing_median)
    assert cfg.anchor_order_agreement == 1.0
    assert cfg.min_edge_support == 2
    assert cfg.sense == 1


def test_acw_config(umbilicus):
    nodes = build_nodes({}, {}, {1: _anchor_collection(1, range(1, 12))}, umbilicus,
                        spiral_sense="ACW")
    cfg = calibrate_from_absolute_anchors(nodes, spiral_sense="acw")
    assert cfg.spiral_sense == "ACW"
    assert cfg.sense == -1


def test_invalid_sense_is_rejected(umbilicus):
    nodes = build_nodes({}, {}, {1: _anchor_collection(1, range(1, 12))}, umbilicus)
    with pytest.raises(ValueError):
        calibrate_from_absolute_anchors(nodes, spiral_sense="left")


def test_config_hash_is_stable(umbilicus):
    nodes = build_nodes({}, {}, {1: _anchor_collection(1, range(1, 12))}, umbilicus)
    assert (calibrate_from_absolute_anchors(nodes).config_hash
            == calibrate_from_absolute_anchors(nodes).config_hash)


def test_too_few_anchors_is_an_error(umbilicus):
    nodes = build_nodes({}, {}, {1: _anchor_collection(1, [1, 2])}, umbilicus)
    with pytest.raises(ValueError):
        calibrate_from_absolute_anchors(nodes)
