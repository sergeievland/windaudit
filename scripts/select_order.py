"""Recorded strict-null selection. No rejected margin is silently frozen."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
import numpy as np
from windaudit.pcl import load_point_collections
from windaudit.frame import Umbilicus,build_nodes,KIND_RELATIVE,KIND_SAME
from windaudit.calibrate import calibrate_from_absolute_anchors
from windaudit.controls import _jitter
from windaudit.matching import build_links
from windaudit.order_constraints import *
from windaudit.order_controls import sites,planted,conflict_nodes

GRID=[0.5,0.75,1.0,1.25,1.5,2.0,3.0]

def main(out):
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    report=dict(grid=GRID,seed=777,n_null_rounds=20,sigma_voxels=1.0,
                selection_rule='Reject any margin with a new canonical-quarantine alarm in the 20 one-voxel Gaussian null rounds; rank eligible margins by the sum of the two exhaustive planted alarm rates, tie by smaller margin.',
                protocol_difference='Legacy selection tolerated 2/20 null alarms; this extension requires zero. Null draws and mutation families otherwise follow the legacy protocol.',
                baseline_alarm_definition='Canonical quarantine members plus both endpoints of conflicting equality pairs, matching the legacy alarm convention; null alarm means a newly implicated member. Existing inconsistencies are recorded separately, never labelled confirmed annotation errors.',
                selected=None,rows=[])
    def save():out.write_text(json.dumps(report,indent=1)+'\n')
    save() # record grid and rule before evaluating it
    p=Path('data/paris4');rel,same,absolute=[load_point_collections(str(p/f)) for f in ['relative_windings.json','same_windings.json','abs_winding.json']];u=Umbilicus.load(str(p/'umbilicus.json'))
    nodes=build_nodes(rel,same,absolute,u);cfg=calibrate_from_absolute_anchors(nodes)
    match=build_links(nodes,cfg);pairs=candidates(nodes,cfg)
    baselines=[]
    for m in GRID:
        arcs,ev=build_constraints(nodes,cfg,m,match,pairs)
        rep=repair(arcs,targets=[],method='cycles')
        if not rep['optimal']:raise RuntimeError('unverified baseline repair')
        row=dict(margin=m,configuration=configuration(cfg,m),n_collections=len({k for a in arcs for k in (a.u,a.v)}),n_arcs=len(arcs),n_relation_groups=len({a.group for a in arcs}),n_point_pairs=len({(a,pa,b,pb) for a,pa,b,pb,g in ev}),baseline_repair=rep,certificate=cycle_certificate(arcs),null_rounds=[],eligible=None)
        row['coverage']={}
        for name,cols,kind in [('skipped_wrap',rel,KIND_RELATIVE),('sheet_switch',same,KIND_SAME)]:
            ss,total=sites(arcs,ev,cols,kind,cfg.min_edge_support)
            row['coverage'][name]=dict(reachable_sites=len(ss),boundaries_total=total,site_coverage=len(ss)/total if total else None)
        report['rows'].append(row);baselines.append((nodes,arcs,ev,rep));save();print('baseline',m,row['n_collections'],row['coverage'],flush=True)
    rng=np.random.default_rng(777)
    for i in range(20):
        r=_jitter(rel,rng,1.);s=_jitter(same,rng,1.)
        ns=build_nodes(r,s,absolute,u);mt=build_links(ns,cfg);ps=candidates(ns,cfg)
        for row in report['rows']:
            aa,_=build_constraints(ns,cfg,row['margin'],mt,ps)
            rr=repair(aa,targets=[],method='cycles')
            if not rr['optimal']:raise RuntimeError('unverified null repair')
            new=sorted((set(rr['removed'])|conflict_nodes(aa))-set(row['baseline_repair']['removed']))
            row['null_rounds'].append(dict(round=i,new_implicated=new,alarm=bool(new)))
            save()
        print('null',i,[sum(t['alarm'] for t in r['null_rounds']) for r in report['rows']],flush=True)
    for row,base in zip(report['rows'],baselines):
        row['null_rounds_with_alarm']=sum(t['alarm'] for t in row['null_rounds'])
        row['eligible']=row['null_rounds_with_alarm']==0
        if row['eligible']:
            row['planted_defect_validation']=planted(rel,same,absolute,u,cfg,row['margin'],base)
        else:
            row['planted_defect_validation']=None
            row['rejection_reason']='Non-zero strict null false-alarm count; not eligible for planted-based selection.'
        save()
    eligible=[r for r in report['rows'] if r['eligible']]
    if eligible:
        best=max(eligible,key=lambda r:(sum(r['planted_defect_validation'][k]['alarm_rate'] or 0 for k in ['skipped_wrap','sheet_switch']),-r['margin']))
        report['selected']=best['margin'];report['frozen_configuration']=best['configuration']
    else:
        report['status']='no_admissible_margin'
        report['frozen_configuration']=None
    save()
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'results/selection_grid_order.json')
