"""Diagnostic evaluation of a preregistered reference margin, not selection.

Evaluate every legacy reachable boundary plus a fixed random sample of newly
covered boundaries. Full new coverage is enumerated before sampling. A rejected
margin never becomes deployable merely because its sample detects defects.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import concurrent.futures
import multiprocessing
import json
import numpy as np
from windaudit.pcl import load_point_collections
from windaudit.frame import Umbilicus,build_nodes,KIND_RELATIVE,KIND_SAME,node_tag
from windaudit.calibrate import calibrate_from_absolute_anchors
from windaudit.matching import build_links
from windaudit.controls import enumerate_sites,plant_step_defect,plant_sheet_switch
from windaudit.order_constraints import build_constraints,repair,cycle_certificate,configuration
from windaudit.order_controls import sites,audit,conflict_nodes
MARGIN=1.0

def initialise():
    global REL,SAME,ABS,U,CFG,BASE
    p=Path('data/paris4');REL,SAME,ABS=[load_point_collections(str(p/f)) for f in ['relative_windings.json','same_windings.json','abs_winding.json']];U=Umbilicus.load(str(p/'umbilicus.json'))
    ns=build_nodes(REL,SAME,ABS,U);CFG=calibrate_from_absolute_anchors(ns)
    aa,_=build_constraints(ns,CFG,MARGIN)
    BASE=set(repair(aa,targets=[],method='cycles')['removed'])

def trial(job):
    cohort,name,site,direction=job
    kind=KIND_RELATIVE if name=='skipped_wrap' else KIND_SAME
    key=(kind,site[0]);tag=node_tag(key)
    r=plant_step_defect(REL,site,direction) if name=='skipped_wrap' else REL
    s=plant_sheet_switch(SAME,U,CFG.wrap_spacing_median,site,direction) if name=='sheet_switch' else SAME
    _,aa,ev,rep=audit(r,s,ABS,U,CFG,MARGIN,[key],time_budget=60.)
    if not rep['optimal']:
        return dict(cohort=cohort,defect=name,collection=tag,boundary=site[1],direction=direction,alarm=None,localized=None,membership='unknown',certainly_named=None,baseline_already_quarantined=tag in BASE,solver_status=rep.get('status','unresolved'))
    new=(set(rep['removed'])|conflict_nodes(aa))-BASE;membership=rep['membership'][tag]
    return dict(cohort=cohort,defect=name,collection=tag,boundary=site[1],direction=direction,
                alarm=bool(new),localized=tag in new,membership=membership,certainly_named=membership=='certain' and tag not in BASE,
                baseline_already_quarantined=tag in BASE)

def main(out="results/paris4/order_diagnostic.json",resume=True):
    initialise();ns=build_nodes(REL,SAME,ABS,U);match=build_links(ns,CFG);arcs,ev=build_constraints(ns,CFG,MARGIN,match)
    rng=np.random.default_rng(20260918);jobs=[];counts={}
    for name,cols,kind in [('skipped_wrap',REL,KIND_RELATIVE),('sheet_switch',SAME,KIND_SAME)]:
        old,total=enumerate_sites(match,cols,kind,CFG.min_edge_support)
        new,_=sites(arcs,ev,cols,kind,CFG.min_edge_support)
        oldids={(s[0],s[1]) for s in old};fresh=[s for s in new if (s[0],s[1]) not in oldids]
        ids=sorted(rng.choice(len(fresh),min(32,len(fresh)),replace=False))
        sample=[fresh[i] for i in ids]
        counts[name]=dict(boundaries_total=total,legacy_sites=len(old),extended_structural_sites=len(new),new_only_sites=len(fresh),sampled_new_only_sites=len(sample))
        jobs.extend((cohort,name,site,d) for cohort,ss in [('all_legacy_sites',old),('sample_new_only_sites',sample)] for site in ss for d in [-1,1])
    path=Path(out);path.parent.mkdir(parents=True,exist_ok=True)
    result=dict(margin=MARGIN,margin_status='diagnostic_only_not_selected',sample_seed=20260918,configuration=configuration(CFG,MARGIN),coverage=counts,baseline_removed=sorted(BASE),trials=[],summary={})
    if resume and path.exists():
        previous=json.loads(path.read_text())
        if previous['configuration']==result['configuration']:
            result['trials']=previous['trials']
    result['repair_time_budget_seconds']=60
    def key(row):return (row['cohort'],row['defect'],row['collection'],row['boundary'],row['direction'])
    done={key(row) for row in result['trials']}
    pending=[j for j in jobs if (j[0],j[1],node_tag((KIND_RELATIVE if j[1]=='skipped_wrap' else KIND_SAME,j[2][0])),j[2][1],j[3]) not in done]
    path.write_text(json.dumps(result,indent=1)+'\n')
    with concurrent.futures.ProcessPoolExecutor(max_workers=4,initializer=initialise,mp_context=multiprocessing.get_context("spawn")) as pool:
        futures=[pool.submit(trial,j) for j in pending]
        for i,future in enumerate(concurrent.futures.as_completed(futures)):
            row=future.result()
            result['trials'].append(row)
            result['trials'].sort(key=key)
            path.write_text(json.dumps(result,indent=1)+'\n')
            print(i+1,len(jobs),row['cohort'],row['defect'],flush=True)
    for cohort in ['all_legacy_sites','sample_new_only_sites']:
        result['summary'][cohort]={}
        for name in counts:
            rr=[r for r in result['trials'] if r['cohort']==cohort and r['defect']==name]
            resolved=[r for r in rr if r['alarm'] is not None]
            result['summary'][cohort][name]=dict(n_trials=len(rr),n_resolved=len(resolved),n_unknown=len(rr)-len(resolved),alarm_rate_resolved=sum(r['alarm'] for r in resolved)/len(resolved) if resolved else None,certain_localization_rate_resolved=sum(r['certainly_named'] for r in resolved)/len(resolved) if resolved else None,alarm_rate_all_trials_bounds=[sum(r['alarm'] for r in resolved)/len(rr),(sum(r['alarm'] for r in resolved)+len(rr)-len(resolved))/len(rr)] if rr else None,certain_rate_all_trials_bounds=[sum(r['certainly_named'] for r in resolved)/len(rr),(sum(r['certainly_named'] for r in resolved)+len(rr)-len(resolved))/len(rr)] if rr else None,membership_counts={v:sum(r['membership']==v for r in rr) for v in ['certain','possible','excluded','consistent','unknown']})
    path.write_text(json.dumps(result,indent=1)+'\n')
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',default='results/paris4/order_diagnostic.json');ap.add_argument('--fresh',action='store_true');args=ap.parse_args()
    main(args.out,not args.fresh)
