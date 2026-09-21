"""Order-only additions to the frozen legacy measurement pipeline."""
from collections import defaultdict
from pathlib import Path
import json
import numpy as np
from .controls import _ordered_pids, _jitter, plant_step_defect, plant_sheet_switch
from .frame import build_nodes, KIND_RELATIVE, KIND_SAME, node_tag
from .matching import build_links
from .order_constraints import build_constraints, candidates, repair, cycle_certificate


def conflict_nodes(arcs):
    """Legacy alarm convention: both endpoints of conflicting equalities."""
    relations=defaultdict(set)
    for a in arcs:
        if a.group.startswith('equality:'):
            relations[tuple(sorted((a.u,a.v)))].add(a.group)
    return {node_tag(k) for pair,groups in relations.items() if len(groups)>1 for k in pair}


def sites(arcs, evidence, collections, kind, min_support):
    """Directed structural coverage, not a claim of sensitivity to a 1-wrap step.

    Two sides must close a directed cycle through evidence independent of the
    target. Each attachment direction needs the frozen minimum support.
    Inequality slack can still hide the planted defect.
    """
    by_group=defaultdict(list)
    for a in arcs: by_group[a.group].append(a)
    attachments=defaultdict(list)
    for a,pa,b,pb,g in evidence:
        for arc in by_group[g]:
            attachments[a].append((pa,b,'out' if arc.u==a else 'in'))
            attachments[b].append((pb,a,'out' if arc.u==b else 'in'))
    result=[];total=0
    for cid,col in sorted(collections.items()):
        pids=_ordered_pids(col,kind)
        if len(pids)<4:continue
        total+=len(pids)-1
        k=(kind,cid); ats=attachments[k]
        if not ats:continue
        adj=defaultdict(set)
        for a in arcs:
            if k not in (a.u,a.v):adj[a.u].add(a.v)
        reach={}
        for partner in {p for _,p,_ in ats}:
            visited={partner};stack=[partner]
            while stack:
                v=stack.pop()
                for w in adj[v]-visited:visited.add(w);stack.append(w)
            reach[partner]=visited
        idx={pid:i for i,pid in enumerate(pids)}
        for boundary in range(1,len(pids)):
            counts=defaultdict(int)
            for pid,p,direction in ats:
                if pid in idx:counts[(idx[pid]<boundary,direction,p)]+=1
            supported={side:{d:{p for (s,di,p),v in counts.items() if s==side and di==d and v>=min_support} for d in ['in','out']} for side in [True,False]}
            ok=any(reach[a]&supported[not side]['in'] for side in [True,False] for a in supported[side]['out'])
            if ok:result.append((cid,boundary,pids))
    return result,total


def audit(relative,same,absolute,umbilicus,cfg,margin,targets=(),time_budget=None):
    nodes=build_nodes(relative,same,absolute,umbilicus,spiral_sense=cfg.spiral_sense,max_step_rad=cfg.max_chain_step_rad)
    arcs,ev=build_constraints(nodes,cfg,margin)
    result=repair(arcs,targets=list(targets),method='cycles',time_budget=time_budget)
    if not result['optimal'] and time_budget is None:raise RuntimeError('order repair did not reach a verified optimum')
    return nodes,arcs,ev,result


def planted(relative,same,absolute,umbilicus,cfg,margin,baseline, site_lists=None,progress=None):
    nodes,arcs,ev,base=baseline
    base_removed=set(base['removed'])|conflict_nodes(arcs)
    if site_lists is None:
        site_lists={name:sites(arcs,ev,cols,kind,cfg.min_edge_support) for name,cols,kind in [('skipped_wrap',relative,KIND_RELATIVE),('sheet_switch',same,KIND_SAME)]}
    rows=[];summary={}
    for name,cols,kind in [('skipped_wrap',relative,KIND_RELATIVE),('sheet_switch',same,KIND_SAME)]:
        ss,total=site_lists[name]
        for cid,boundary,pids in ss:
            key=(kind,cid);tag=node_tag(key)
            for direction in [-1,1]:
                r=plant_step_defect(relative,(cid,boundary,pids),direction) if name=='skipped_wrap' else relative
                s=plant_sheet_switch(same,umbilicus,cfg.wrap_spacing_median,(cid,boundary,pids),direction) if name=='sheet_switch' else same
                _,aa,ee,result=audit(r,s,absolute,umbilicus,cfg,margin,[key])
                new=(set(result['removed'])|conflict_nodes(aa))-base_removed
                status=result['membership'][tag]
                row=dict(defect=name,collection=tag,boundary=boundary,direction=direction,
                         alarm=bool(new),localized=tag in new,membership=status,
                         certainly_named=status=='certain' and tag not in base_removed)
                rows.append(row)
                if progress:progress(row)
        sub=[r for r in rows if r['defect']==name]
        summary[name]=dict(boundaries_total=total,reachable_sites=len(ss),site_coverage=len(ss)/total if total else None,n_trials=len(sub),
                           alarm_rate=sum(r['alarm'] for r in sub)/len(sub) if sub else None,
                           certain_localization_rate=sum(r['certainly_named'] for r in sub)/len(sub) if sub else None,
                           localization_rate=sum(r['localized'] for r in sub)/len(sub) if sub else None,
                           membership_counts={v:sum(r['membership']==v for r in sub) for v in ['certain','possible','excluded','consistent','unknown']})
    return dict(**summary,trials=rows)
