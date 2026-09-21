"""Rebuild deterministic order graph and independently check saved evidence.

This validates saved measurements structurally; expensive null and planted
remeasurement is explicitly separate in reproduce_order.sh.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
from windaudit.pcl import load_point_collections
from windaudit.frame import build_nodes,Umbilicus,KIND_RELATIVE,KIND_SAME
from windaudit.calibrate import calibrate_from_absolute_anchors
from windaudit.matching import build_links
from windaudit.order_constraints import candidates,build_constraints,configuration,cycle_certificate
from windaudit.order_controls import sites

def main():
    p=Path('data/paris4');r,s,a=[load_point_collections(str(p/f)) for f in ['relative_windings.json','same_windings.json','abs_winding.json']];u=Umbilicus.load(str(p/'umbilicus.json'))
    n=build_nodes(r,s,a,u);cfg=calibrate_from_absolute_anchors(n);match=build_links(n,cfg);pairs=candidates(n,cfg)
    grid=json.loads(Path('results/selection_grid_order.json').read_text())
    assert grid['selected'] is None and grid['status']=='no_admissible_margin'
    for row in grid['rows']:
        arcs,ev=build_constraints(n,cfg,row['margin'],match,pairs)
        assert row['configuration']==configuration(cfg,row['margin'])
        assert row['certificate']==cycle_certificate(arcs)
        assert row['n_arcs']==len(arcs)
        assert row['n_collections']==len({k for x in arcs for k in (x.u,x.v)})
        assert len(row['null_rounds'])==20
        assert row['null_rounds_with_alarm']==sum(t['alarm'] for t in row['null_rounds'])>0
        for trial in row['null_rounds']: assert trial['alarm']==bool(trial['new_implicated'])
        for name,cols,kind in [('skipped_wrap',r,KIND_RELATIVE),('sheet_switch',s,KIND_SAME)]:
            ss,total=sites(arcs,ev,cols,kind,cfg.min_edge_support)
            assert row['coverage'][name]==dict(reachable_sites=len(ss),boundaries_total=total,site_coverage=len(ss)/total)
    report=json.loads(Path('results/paris4/order_report.json').read_text())
    arcs,ev=build_constraints(n,cfg,report['margin'],match,pairs)
    assert report['configuration']==configuration(cfg,report['margin'])
    assert report['certificate']==cycle_certificate(arcs)
    cert=report['certificate'];assert cert['total_bound']<0
    assert sum(x['bound'] for x in cert['arcs'])==cert['total_bound']
    for i,arc in enumerate(cert['arcs']):assert arc['v']==cert['arcs'][(i+1)%len(cert['arcs'])]['u']
    diag=json.loads(Path('results/paris4/order_diagnostic.json').read_text())
    identities=[(t['cohort'],t['defect'],t['collection'],t['boundary'],t['direction']) for t in diag['trials']]
    assert len(identities)==len(set(identities))==238
    for cohort,families in diag['summary'].items():
        for name,summary in families.items():
            rr=[t for t in diag['trials'] if t['cohort']==cohort and t['defect']==name]
            assert len(rr)==summary['n_trials']
            assert summary['n_unknown']==sum(t['membership']=='unknown' for t in rr)
            assert summary['n_unknown']+summary['n_resolved']==len(rr)
    print('Order graph, coverage, hashes, negative-cycle witness, 140 null records and 238 diagnostic identities verified.')
if __name__=='__main__':main()
