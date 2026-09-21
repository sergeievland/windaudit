import math

from windaudit.controls import (
    enumerate_sites,
    plant_sheet_switch,
    plant_step_defect,
    planted_defect_validation,
)
from windaudit.frame import KIND_RELATIVE, KIND_SAME, build_nodes
from windaudit.matching import build_links
from tests.conftest import ladder, wrap_trace


def _scene():
    rel = {
        1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 8)),
        2: ladder(2, theta=1.0, z=800.0, wraps=range(1, 8)),
    }
    same = {
        w: wrap_trace(w, wrap=w, theta=1.0, zs=range(470, 840, 20))
        for w in (3, 4, 5)
    }
    return rel, same


def test_sites_exist_in_bridged_scene(umbilicus, config):
    rel, same = _scene()
    match = build_links(build_nodes(rel, same, {}, umbilicus), config)
    step_sites, step_total = enumerate_sites(match, rel, KIND_RELATIVE, 2)
    switch_sites, switch_total = enumerate_sites(match, same, KIND_SAME, 2)
    assert step_sites and switch_sites
    assert step_total == 2 * 6
    assert switch_total == 3 * 18
    assert len(step_sites) <= step_total


def test_plant_step_defect_shifts_tail_only(umbilicus):
    rel, _ = _scene()
    pids = [p.id for p in rel[1].ordered_points()]
    mutated = plant_step_defect(rel, (1, 3, pids), shift=1)
    for pid in pids[:3]:
        assert mutated[1].points[pid].wind_a == rel[1].points[pid].wind_a
    for pid in pids[3:]:
        assert mutated[1].points[pid].wind_a == rel[1].points[pid].wind_a + 1
    assert mutated[2] is rel[2]


def test_plant_sheet_switch_moves_tail_radially(umbilicus, config):
    _, same = _scene()
    pids = [p.id for p in same[3].ordered_points()]
    mutated = plant_sheet_switch(same, umbilicus, config.wrap_spacing_median,
                                 (3, 5, pids), direction=1.0)
    for pid in pids[:5]:
        assert mutated[3].points[pid].xyz == same[3].points[pid].xyz
    r_new, t_new = umbilicus.cylindrical(mutated[3].points[pids[5]].xyz)
    r_old, t_old = umbilicus.cylindrical(same[3].points[pids[5]].xyz)
    assert math.isclose(r_new - r_old, config.wrap_spacing_median, abs_tol=1e-6)
    assert math.isclose(t_new, t_old, abs_tol=1e-9)


def test_planted_validation_on_bridged_scene(umbilicus, config):
    rel, same = _scene()
    result = planted_defect_validation(rel, same, {}, umbilicus, config,
                                       seed=3, n_null_rounds=3)
    assert result["baseline_implicated"] == []
    assert result["skipped_wrap"]["n_trials"] > 0
    assert result["skipped_wrap"]["alarm_rate"] == 1.0
    assert result["sheet_switch"]["alarm_rate"] > 0.5
    assert result["null_jitter"]["rounds_with_alarm"] == 0
