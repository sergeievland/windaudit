"""Golden regression on the bundled PHercParis4 corpus: a fresh audit must
reproduce the committed report exactly (fast sections; the sensitivity sweep
and planted-defect validation are reproduced by scripts/reproduce.sh)."""

import json
import os

import pytest

from windaudit.report import run_audit

ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA = os.path.join(ROOT, "data", "paris4")
GOLDEN = os.path.join(ROOT, "results", "paris4", "audit_report.json")


@pytest.fixture(scope="module")
def report():
    return run_audit(
        relative_path=os.path.join(DATA, "relative_windings.json"),
        same_path=os.path.join(DATA, "same_windings.json"),
        absolute_path=os.path.join(DATA, "abs_winding.json"),
        umbilicus_path=os.path.join(DATA, "umbilicus.json"),
        out_dir=None,
        run_planted=False,
        run_sweep=False,
    )


@pytest.fixture(scope="module")
def golden():
    with open(GOLDEN, encoding="utf-8") as f:
        return json.load(f)


def _plain(obj):
    return json.loads(json.dumps(obj))


def test_inputs_are_the_pinned_files(report, golden):
    assert report["input_sha256"] == golden["input_sha256"]


@pytest.mark.parametrize("section", [
    "corpus", "calibration", "config_sha256", "matching", "graph", "coverage",
    "failing_cycles",
])
def test_section_reproduces(report, golden, section):
    assert _plain(report[section]) == golden[section]


def test_tier1_reproduces(report, golden):
    live = _plain(report["tier1"])
    assert live["order_audit"] == golden["tier1"]["order_audit"]
    assert live["same_winding_step_audit"] == golden["tier1"]["same_winding_step_audit"]
    assert live["shuffled_control"] == golden["tier1"]["shuffled_control"]


def test_repair_and_certified_links_reproduce(report, golden):
    live = _plain(report["repair"])
    for key in ("min_quarantine", "min_edge_cut"):
        assert live[key]["removed" if key == "min_quarantine" else "cut_edges"] == \
            golden["repair"][key]["removed" if key == "min_quarantine" else "cut_edges"]
        assert live[key]["optimal"] == golden["repair"][key]["optimal"]
    assert report["certified_links"]["n_links"] == golden["certified_links"]["n_links"]


def test_headline_numbers(golden):
    assert golden["tier1"]["order_audit"]["n_pairs_tested"] == 11447
    assert golden["tier1"]["order_audit"]["n_violations"] == 0
    assert golden["matching"]["n_links"] == 2371
    assert golden["matching"]["n_conflict_pairs"] == 0
    assert golden["graph"]["n_effective_nonzero_holonomy"] == 0
    planted = golden["planted_defect_validation"]
    assert planted["skipped_wrap"]["alarm_rate"] == 1.0
    assert planted["null_jitter"]["rounds_with_alarm"] == 0
