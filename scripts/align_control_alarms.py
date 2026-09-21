"""Include equality-conflict endpoints in recorded incremental alarm counts.

The first exploratory runs stored only canonical quarantines. This transparent
replay uses the identical mutations, adds the legacy conflict-endpoint rule,
and never changes a solver result or an all-optima membership.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import concurrent.futures
import json
import numpy as np
from windaudit.pcl import load_point_collections
from windaudit.frame import Umbilicus,build_nodes,node_tag,KIND_RELATIVE,KIND_SAME
from windaudit.calibrate import calibrate_from_absolute_anchors
from windaudit.matching import build_links
from windaudit.controls import _jitter,_ordered_pids,plant_step_defect,plant_sheet_switch

def initialise():
    global R,S,A,U,CFG
    p=Path('data/paris4');R,S,A=[load_point_collections(str(p/f)) for f in ['relative_windings.json','same_windings.json','abs_winding.json']];U=Umbilicus.load(str(p/'umbilicus.json'));CFG=calibrate_from_absolute_anchors(build_nodes(R,S,A,U))
def conflicts(r,s):
    match=build_links(build_nodes(r,s,A,U),CFG)
    return sorted({node_tag(k) for c in match.conflicts for k in [c.a,c.b]})
def trial(t):
    kind,cid=t['collection'].split(':');cid=int(cid)
    col=(R if kind==KIND_RELATIVE else S)[cid]
    site=(cid,t['boundary'],_ordered_pids(col,kind))
    r=plant_step_defect(R,site,t['direction']) if kind==KIND_RELATIVE else R
    s=plant_sheet_switch(S,U,CFG.wrap_spacing_median,site,t['direction']) if kind==KIND_SAME else S
    return conflicts(r,s)
def summaries(result):
    result['summary']={}
    for cohort in ['all_legacy_sites','sample_new_only_sites']:
        result['summary'][cohort]={}
        for name in ['skipped_wrap','sheet_switch']:
            rr=[r for r in result['trials'] if r['cohort']==cohort and r['defect']==name]
            resolved=[r for r in rr if r['alarm'] is not None]
            result['summary'][cohort][name]=dict(n_trials=len(rr),n_resolved=len(resolved),n_unknown=len(rr)-len(resolved),alarm_rate_resolved=sum(r['alarm'] for r in resolved)/len(resolved) if resolved else None,certain_localization_rate_resolved=sum(r['certainly_named'] for r in resolved)/len(resolved) if resolved else None,alarm_rate_all_trials_bounds=[sum(r['alarm'] for r in resolved)/len(rr),(sum(r['alarm'] for r in resolved)+len(rr)-len(resolved))/len(rr)] if rr else None,certain_rate_all_trials_bounds=[sum(r['certainly_named'] for r in resolved)/len(rr),(sum(r['certainly_named'] for r in resolved)+len(rr)-len(resolved))/len(rr)] if rr else None,membership_counts={v:sum(r['membership']==v for r in rr) for v in ['certain','possible','excluded','consistent','unknown']})
def main():
    initialise();p=Path('results/selection_grid_order.json');grid=json.loads(p.read_text())
    assert all(len(row['null_rounds'])==20 for row in grid['rows'])
    rng=np.random.default_rng(777)
    for i in range(20):
        extra=conflicts(_jitter(R,rng,1.),_jitter(S,rng,1.))
        for row in grid['rows']:
            t=row['null_rounds'][i];new=set(t['new_implicated'])|(set(extra)-set(row['baseline_repair']['removed']))
            t['new_implicated']=sorted(new);t['alarm']=bool(new)
    for row in grid['rows']:row['null_rounds_with_alarm']=sum(t['alarm'] for t in row['null_rounds'])
    grid['baseline_alarm_definition']='Canonical quarantine members plus both endpoints of conflicting equality pairs, matching the legacy alarm convention; null alarm means a newly implicated member. Existing inconsistencies are recorded separately, never labelled confirmed annotation errors.'
    p.write_text(json.dumps(grid,indent=1)+'\n')
    p=Path('results/paris4/order_diagnostic.json');result=json.loads(p.read_text());assert len(result['trials'])==238
    with concurrent.futures.ProcessPoolExecutor(max_workers=4,initializer=initialise) as pool:
        for t,extra in zip(result['trials'],pool.map(trial,result['trials'])):
            extra=set(extra)-set(result['baseline_removed'])
            if t['alarm'] is not None:
                t['alarm']=bool(t['alarm'] or extra)
                t['localized']=bool(t['localized'] or t['collection'] in extra)
    summaries(result);p.write_text(json.dumps(result,indent=1)+'\n')
    print('Replayed 20 identical null mutations and 238 identical planted mutations for the legacy conflict-endpoint alarm convention.')
if __name__=='__main__':main()
