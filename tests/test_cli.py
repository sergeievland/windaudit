import json

import pytest

from windaudit.cli import main
from windaudit.pcl import load_point_collections, save_point_collections
from tests.conftest import ladder, make_collection, wrap_trace, wrap_xyz


def _write_scene(tmp_path):
    rel = {
        1: ladder(1, theta=1.0, z=500.0, wraps=range(1, 8)),
        2: ladder(2, theta=1.0, z=800.0, wraps=range(1, 8)),
    }
    same = {w: wrap_trace(w, wrap=w, theta=1.0, zs=range(470, 840, 20))
            for w in (3, 4, 5)}
    rows = [tuple(wrap_xyz(w, 1.0 + 0.001 * w, z)) + (w,)
            for z in (500.0, 560.0) for w in range(1, 10)]
    absolute = {1: make_collection(1, rows, absolute=True)}
    save_point_collections(str(tmp_path / "relative_windings.json"), rel)
    save_point_collections(str(tmp_path / "same_windings.json"), same)
    save_point_collections(str(tmp_path / "abs_winding.json"), absolute)
    umb = {"control_points": [{"x": 0, "y": 0, "z": -1000}, {"x": 0, "y": 0, "z": 20000}]}
    (tmp_path / "umbilicus.json").write_text(json.dumps(umb))


def test_cli_end_to_end(tmp_path, capsys):
    _write_scene(tmp_path)
    out = tmp_path / "out"
    assert main(["run", "--inputs", str(tmp_path), "--out", str(out),
                 "--null-rounds", "2", "--shuffle-rounds", "2"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["failing_cycles"] == 0
    assert summary["certified_links"] > 0
    report = json.loads((out / "audit_report.json").read_text())
    assert report["planted_defect_validation"]["skipped_wrap"]["alarm_rate"] == 1.0
    certified = load_point_collections(str(out / "certified_links.json"))
    assert len(certified) == summary["certified_links"]
    assert all(len(c.points) == 2 for c in certified.values())


def test_cli_reports_missing_input(tmp_path):
    with pytest.raises(SystemExit):
        main(["run", "--inputs", str(tmp_path), "--out", str(tmp_path / "o")])
