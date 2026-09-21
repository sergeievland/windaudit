"""Execute the genuine pinned upstream loader, isolated from GPU imports.

No reimplemented parser and no mocked parser: compile the original function
AST. This tests file-format compatibility, not upstream surface linking.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Optional, Dict, Any
ROOT=Path(__file__).resolve().parents[1]

def test_export_with_pinned_upstream_loader(tmp_path):
    path=ROOT/'upstream/villa/point_collection.py'
    provenance=json.loads((path.parent/'provenance.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==provenance['sha256']
    module=ast.parse(path.read_text())
    f=next(n for n in module.body if isinstance(n,ast.FunctionDef) and n.name=='load_point_collection')
    ns=dict(json=json,Optional=Optional,Dict=Dict,Any=Any)
    exec(compile(ast.Module(body=[f],type_ignores=[]),str(path),'exec'),ns)
    source=ROOT/'results/paris4/certified_links.json'
    loaded=ns['load_point_collection'](str(source))
    raw=json.loads(source.read_text())['collections']
    assert len(loaded)==2371
    for cid,col in loaded.items():
        assert len(col['points'])==2
        for pid,p in col['points'].items():
            assert p['p']==raw[str(cid)]['points'][str(pid)]['p']
            assert math.isnan(p['winding_annotation'])
    bad=tmp_path/'bad.json';bad.write_text('{"vc_pointcollections_json_version":"99"}')
    assert ns['load_point_collection'](str(bad)) is None
