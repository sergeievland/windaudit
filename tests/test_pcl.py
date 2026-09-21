import json

import pytest

from windaudit.pcl import load_point_collections, save_point_collections
from tests.conftest import make_collection


def test_round_trip(tmp_path):
    cols = {
        1: make_collection(1, [(1, 2, 3, 5.0), (4, 5, 6, 6.0)]),
        2: make_collection(2, [(7, 8, 9, None)], absolute=False),
        3: make_collection(3, [(1, 1, 1, 12.0)], absolute=True),
    }
    path = str(tmp_path / "pcl.json")
    save_point_collections(path, cols)
    loaded = load_point_collections(path)
    assert sorted(loaded) == [1, 2, 3]
    assert loaded[1].points[1].wind_a == 5.0
    assert loaded[1].points[2].xyz == [4.0, 5.0, 6.0]
    assert loaded[2].points[1].wind_a is None
    assert not loaded[2].has_winding_annotations
    assert loaded[3].is_absolute
    assert loaded[3].has_winding_annotations


def test_version_check(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"vc_pointcollections_json_version": "2",
                                "collections": {}}))
    with pytest.raises(ValueError):
        load_point_collections(str(path))


def test_chain_order_is_point_id_order(tmp_path):
    data = {
        "vc_pointcollections_json_version": "1",
        "collections": {
            "7": {
                "name": "c",
                "points": {
                    "10": {"p": [0, 0, 0], "wind_a": 1.0},
                    "2": {"p": [1, 0, 0], "wind_a": 2.0},
                    "30": {"p": [2, 0, 0], "wind_a": 3.0},
                },
            }
        },
    }
    path = tmp_path / "pcl.json"
    path.write_text(json.dumps(data))
    loaded = load_point_collections(str(path))
    assert [p.id for p in loaded[7].ordered_points()] == [2, 10, 30]

@pytest.mark.parametrize('point', [
    {'p':[1,2]}, {'p':[1,2,3,4]}, {'p':[1,float('nan'),3]},
    {'p':[1,2,3],'wind_a':float('inf')},
])
def test_reject_invalid_numeric_inputs(tmp_path, point):
    path=tmp_path/'bad.json'
    path.write_text(json.dumps({'vc_pointcollections_json_version':'1',
                               'collections':{'1':{'points':{'1':point}}}}))
    with pytest.raises(ValueError):
        load_point_collections(str(path))


def test_reject_normalised_id_collision(tmp_path):
    path=tmp_path/'bad.json'
    path.write_text(json.dumps({'vc_pointcollections_json_version':'1',
                               'collections':{'1':{'points':{}},'01':{'points':{}}}}))
    with pytest.raises(ValueError,match='duplicate'):
        load_point_collections(str(path))
