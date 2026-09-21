"""Explicit opt-in entrypoint for experimental order constraints."""
import argparse
import json
import math
from dataclasses import replace
from pathlib import Path
from .pcl import load_point_collections
from .frame import build_nodes,Umbilicus,KIND_RELATIVE,KIND_SAME,node_tag
from .calibrate import calibrate_from_absolute_anchors
from .matching import build_links
from .order_constraints import candidates,build_constraints,configuration,repair,cycle_certificate,minimum_edge_cut
from .order_controls import sites
from .report import sha256_file,SWEEP_R_SAME_FACTORS,SWEEP_DTHETA_DEG,SWEEP_DZ

def run(inputs,margin,diagnostic=False,sweep=False):
    if not diagnostic:
        raise ValueError('No validated default margin is distributed. Use --diagnostic for explicitly unvalidated order evidence.')
    p=Path(inputs)
    paths=[p/f for f in ['relative_windings.json','same_windings.json','abs_winding.json','umbilicus.json']]
    r,s,a=[load_point_collections(str(f)) for f in paths[:3]];u=Umbilicus.load(str(paths[3]))
    n=build_nodes(r,s,a,u);cfg=calibrate_from_absolute_anchors(n)
    match=build_links(n,cfg);pairs=candidates(n,cfg);arcs,ev=build_constraints(n,cfg,margin,match,pairs)
    cert=cycle_certificate(arcs)
    rr=repair(arcs,targets=[],method='cycles')
    wanted=set(rr['removed'])|set(cert['collections'])
    targets=[k for k in n if node_tag(k) in wanted]
    rr=repair(arcs,targets=targets,method='cycles',time_budget=120.)
    cut=minimum_edge_cut(arcs,time_budget=120.)
    coverage={}
    for name,cols,kind in [('skipped_wrap',r,KIND_RELATIVE),('sheet_switch',s,KIND_SAME)]:
        ss,total=sites(arcs,ev,cols,kind,cfg.min_edge_support)
        coverage[name]=dict(reachable_sites=len(ss),boundaries_total=total,site_coverage=len(ss)/total if total else None)
    sweeps=[];groups={x['group'] for x in cert['arcs']}
    if sweep:
        for factor in SWEEP_R_SAME_FACTORS:
            for theta in SWEEP_DTHETA_DEG:
                for dz in SWEEP_DZ:
                    c=replace(cfg,r_same=factor*cfg.wrap_spacing_median,sector_dtheta_rad=math.radians(theta),sector_dz=dz)
                    aa,_=build_constraints(n,c,margin)
                    g={x.group for x in aa}
                    sweeps.append(dict(r_same_factor=factor,sector_dtheta_deg=theta,sector_dz=dz,
                                       consistent=cycle_certificate(aa)['consistent'],reference_cycle_relations_present=bool(groups and groups<=g)))
    return dict(extension_version='0.3.0-experimental',status='diagnostic_unvalidated_margin',margin=margin,
                input_sha256={f.name:sha256_file(str(f)) for f in paths},configuration=configuration(cfg,margin),
                assumptions=['Radially ordered local sheets about the chosen umbilicus','Radial projection is accurate enough in z AND angular variation is bounded; non-intersection alone does not establish this'],
                counts=dict(legacy_graph_collections=len({k for e in match.edges for k in e}),extended_graph_collections=len({k for x in arcs for k in (x.u,x.v)}),
                            candidate_comparable_pairs=len(pairs),order_and_equality_point_pairs=len({(a,pa,b,pb) for a,pa,b,pb,g in ev}),arcs=len(arcs),relation_groups=len({x.group for x in arcs})),
                coverage=coverage,certificate=cert,quarantine=rr,edge_cut=cut,
                certificate_membership={k:rr['membership'][k] for k in cert['collections']},
                sensitivity_sweep=sweeps,reference_cycle_robustness=sum(x['reference_cycle_relations_present'] for x in sweeps)/len(sweeps) if sweeps else None,
                ct_verified=False,confirmed_annotation_error=False)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--inputs',default='data/paris4');ap.add_argument('--out',required=True);ap.add_argument('--margin',type=float,required=True);ap.add_argument('--diagnostic',action='store_true');ap.add_argument('--sweep',action='store_true');args=ap.parse_args()
    report=run(args.inputs,args.margin,args.diagnostic,args.sweep);p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(report,indent=1)+'\n')
    print(json.dumps(dict(status=report['status'],counts=report['counts'],certificate=report['certificate']['collections'])))
if __name__=='__main__':main()
