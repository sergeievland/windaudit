"""CPU regression harness: extracts the actual function without GPU imports.
Usage: python check_solver.py /path/to/find_inconsistent_windings.py
Two synthetic regressions, one per solver defect; exit 1 before the patch, 0 after.
"""
import ast
import json
import sys
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


def load(path):
    tree = ast.parse(open(path, encoding='utf-8').read())
    func = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'solve_min_edge_fix')
    ns = dict(np=np, Bounds=Bounds, LinearConstraint=LinearConstraint, milp=milp,
              coo_matrix=coo_matrix, _to_py=lambda x: np.asarray(x).tolist())
    exec(compile(ast.Module(body=[func], type_ignores=[]), path, 'exec'), ns)
    return ns['solve_min_edge_fix']


def fixture(rows):
    reached = {p: {'entry_ij': [0,0], 'acc': 0} for a,b,d in rows for p in [a,b]}
    adj = {p: [] for p in reached}
    for i,(a,b,d) in enumerate(rows):
        adj[a].append(dict(neighbor=b, pcl_id=i, pcl_name=str(i), source_file='test.json',
                           from_point_id=1, to_point_id=2, from_ij=[0,0],to_ij=[0,0],
                           from_zyx=[0,0,0],to_zyx=[0,0,0],winding_delta=d,kind='relative'))
    return reached, adj


def check(path):
    solve = load(path)
    r,a = fixture([('A','B',1000),('B','C',1000),('A','C',2000)])
    result = solve(r,a,lambda *args: 0)
    bounds_ok = result['objective'] == 0 and result['num_edges_changed'] == 0
    # Uneditable A-B=0 and B-C=0 must force the sole editable A-C=1 to change.
    r,a = fixture([('A','B',0),('B','C',0),('A','C',1)])
    result = solve(r,a,lambda *args: 0,allowed_edge_keys={(2,frozenset((1,2)))})
    retained_ok = result['num_edges_considered'] == 3 and result['num_edges_changed'] == 1
    print(json.dumps({'path':path,'large_offset_consistent_graph':bounds_ok,
                      'retain_uneditable_equations':retained_ok},indent=2))
    return bounds_ok and retained_ok

if __name__ == '__main__':
    raise SystemExit(0 if check(sys.argv[1]) else 1)
