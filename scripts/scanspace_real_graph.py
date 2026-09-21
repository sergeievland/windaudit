"""Build a CPU real-patch graph and compare the pinned upstream solver.

Usage: python scripts/scanspace_real_graph.py PATCHES BEFORE.py AFTER.py OUTDIR
Only explicitly present patch folders are read. No network access.
"""
import ast
from collections import deque
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.ndimage import binary_erosion
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from windaudit.scanspace import ScanPatch, PatchGraph, load_umbilicus, theta_at, strip_winding_delta, audit_edge_integrability
from windaudit.scanspace_surface import Surface

ROOT=Path(__file__).resolve().parents[1]
def serial(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,set):return sorted(x)
    raise TypeError(type(x).__name__)

def extract(path,names,extra=None):
    tree=ast.parse(Path(path).read_text())
    body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert {n.name for n in body}==set(names)
    ns=dict(np=np,Bounds=Bounds,LinearConstraint=LinearConstraint,milp=milp,
            coo_matrix=coo_matrix,_to_py=lambda x:np.asarray(x).tolist())
    ns.update(extra or {})
    exec(compile(ast.Module(body=body,type_ignores=[]),str(path),'exec'),ns)
    return ns

class Chain:
    def __init__(self,p):self.p=p
    def iter_chain(self):return iter(sorted(self.p['points'].values(),key=lambda p:p['id']))

def edge_key(e):return (e['pcl_id'],frozenset((e['from_point_id'],e['to_point_id'])))
def check(edges,removed=()):
    adj={}
    for k,e in enumerate(edges):
        if k in removed:continue
        a,b,d=e['P'],e['R'],e['D'];adj.setdefault(a,[]).append((b,d,k));adj.setdefault(b,[]).append((a,-d,k))
    pot={};parent={};bad=[];allowed=set()
    for start in sorted(adj):
        if start in pot:continue
        pot[start]=0;parent[start]=None;q=deque([start])
        while q:
            a=q.popleft()
            for b,d,k in adj[a]:
                if b not in pot:pot[b]=pot[a]+d;parent[b]=(a,k);q.append(b)
                elif pot[b]!=pot[a]+d:
                    if k in bad:continue
                    bad.append(k);allowed.add(k)
                    ancestors={a:[]};n=a;path=[]
                    while parent[n] is not None:
                        n,eid=parent[n];path=path+[eid];ancestors[n]=path
                    n=b;path=[]
                    while n not in ancestors:
                        n,eid=parent[n];path.append(eid)
                    allowed.update(path);allowed.update(ancestors[n])
    return {'consistent':not bad,'inconsistent_non_tree_edges':len(bad),'potentials':pot},allowed

def main(directory,before,after,outdir,inputs=None):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    inputs=Path(inputs) if inputs is not None else ROOT/'data/paris4'
    umbilicus=load_umbilicus(inputs/'umbilicus.json')
    patches={};surfaces={};areas={};patch_report=[];rejected=[]
    for folder in sorted(Path(directory).iterdir()):
        if not folder.is_dir():continue
        patch=ScanPatch.load(folder);valid=np.any(patch.xyz!=-1,axis=-1)
        cells=int(patch.metadata.get('spiral_patch_erode_cells',1))
        if cells>0:patch.xyz[~binary_erosion(valid,iterations=cells,border_value=0)]=-1;patch.__post_init__()
        z=patch.xyz[...,2];in_band=(z>=12417)&(z<12419)&np.any(patch.xyz!=-1,axis=-1)
        # Select by actual continuous quad extent as well as vertex hits: a
        # two-voxel window can fall between grid rows. Record this explicitly.
        v=patch.valid_quad
        corners=np.stack([z[:-1,:-1],z[1:,:-1],z[:-1,1:],z[1:,1:]])
        crosses=v&(corners.min(axis=0)<12419)&(corners.max(axis=0)>=12417)
        if not crosses.any():
            rejected.append({'id':folder.name,'reason':'No valid quad intersects the band after mask handling and erosion.'})
            continue
        patches[folder.name]=patch;surfaces[folder.name]=Surface(patch)
        areas[folder.name]=float(v.sum()/np.prod(patch.metadata['scale']))
        patch_report.append({'id':folder.name,'valid_quads':int(v.sum()),'band_crossing_quads':int(crosses.sum()),'has_vertex_in_band':bool(in_band.any()),'erosion_cells':cells})
    print('loaded surfaces',len(patches),flush=True)
    gate=[]
    for name,patch in patches.items():
        r={'patch':name,**audit_edge_integrability(patch,umbilicus)};gate.append(r)
        print('edge gate',name,r['inconsistent_edges'],flush=True)
    (outdir/'band_transport_gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    if any(r['inconsistent_edges'] for r in gate):
        raise SystemExit('STOP: scan-space path independence failed on selected real patches')
    collections=[];points=[]
    for fn in ['abs_winding.json','relative_windings.json','same_windings.json']:
        raw=json.loads((inputs/fn).read_text())['collections']
        for cid,pcl in sorted(raw.items(),key=lambda item:int(item[0])):
            pp={}
            for pid,p in sorted(pcl['points'].items(),key=lambda item:int(item[0])):
                pp[int(pid)]={**p,'id':int(pid),'zyx':np.asarray(p['p'],np.float32)[::-1],
                              'winding_annotation':float(p['wind_a']) if p.get('wind_a') is not None else float('nan')}
            c={**pcl,'points':pp,'source_file':fn,'source_collection_id':cid}
            c['metadata']=dict(c.get('metadata',{}))
            if fn=='abs_winding.json':c['metadata']['winding_is_absolute']=True
            collections.append(c);points.extend(pp.values())
    all_hits={};start=time.perf_counter()
    for n,p in enumerate(points):
        hits=[]
        for name,surface in surfaces.items():
            hit=surface.project(p['p'],2.5)
            if hit is not None:hits.append({**hit,'id':name,'area':areas[name]})
        if hits:
            hits.sort(key=lambda h:(-h['area'],h['distance'],h['id']))
            p['on_patch']=hits[0];all_hits[id(p)]=hits
        if (n+1)%1000==0:print('projected',n+1,flush=True)
    linking_seconds=time.perf_counter()-start
    cross=[]
    for c in collections:
        attached=sum('on_patch' in p for p in c['points'].values())
        if attached<2 and not c['metadata'].get('winding_is_absolute'):continue
        annotated=[p for p in c['points'].values() if np.isfinite(p['winding_annotation'])]
        c['has_winding_annotations']=bool(annotated)
        if annotated:c['points']={k:p for k,p in c['points'].items() if np.isfinite(p['winding_annotation'])}
        else:
            for p in c['points'].values():p['winding_annotation']=0.0
        c['chain']=Chain(c);cross.append(c)
    def tour_adjustments(transform,dr,tour):
        th=theta_at(np.asarray([p['p'] for p in tour]),umbilicus);d=np.diff(th)
        return np.r_[0,np.cumsum((d>np.pi).astype(int)-(d < -np.pi).astype(int))]
    ns=extract(before,['classify_pcl','build_rel_adjacency'],{'_tour_unwrap_adjustments':tour_adjustments})
    # Upstream already uses list-position identities; there is no collision fix.
    adjacency=ns['build_rel_adjacency'](dict(enumerate(cross)),patches,None,None)
    graphs={};memo={}
    def strip(pid,a,b):
        key=(pid,tuple(a),tuple(b))
        if key not in memo:
            if pid not in graphs:graphs[pid]=PatchGraph(patches[pid])
            r=strip_winding_delta(graphs[pid],a,b,umbilicus)
            memo[key]=None if r is None else r['delta_windings']
        return memo[key]
    reached={};component_roots=[]
    for seed in sorted(adjacency):
        if seed in reached:continue
        component_roots.append(seed);entry=adjacency[seed][0]['from_ij']
        reached[seed]={'acc':0,'entry_ij':entry,'hops':0};q=deque([seed])
        while q:
            a=q.popleft();state=reached[a]
            for e in adjacency[a]:
                b=e['neighbor']
                if b in reached:continue
                d=strip(a,state['entry_ij'],e['from_ij'])
                if d is None:continue
                reached[b]={'acc':state['acc']-d-e['winding_delta'],'entry_ij':e['to_ij'],'hops':state['hops']+1};q.append(b)
    edges=[];seen=set();unmeasurable=0
    for a,links in adjacency.items():
        for e in links:
            b=e['neighbor'];key=edge_key(e)
            if key in seen or a not in reached or b not in reached:continue
            seen.add(key);sa=strip(a,reached[a]['entry_ij'],e['from_ij']);sb=strip(b,reached[b]['entry_ij'],e['to_ij'])
            if sa is None or sb is None:unmeasurable+=1;continue
            edges.append({**e,'P':a,'R':b,'D':int(e['winding_delta']+sa-sb)})
    initial,allowed_indices=check(edges);allowed={edge_key(edges[k]) for k in allowed_indices}
    results={}
    unrestricted={}
    for label,path in [('before',before),('after',after)]:
        solve=extract(path,['solve_min_edge_fix'])['solve_min_edge_fix']
        started=time.perf_counter();result=solve(reached,adjacency,strip,allowed_edge_keys=allowed,time_limit=120.)
        result['wall_seconds']=time.perf_counter()-started
        removed_keys={(e['rel_pcl_id'],frozenset((e['from_point_id'],e['to_point_id']))) for e in result['edges']}
        removed={k for k,e in enumerate(edges) if edge_key(e) in removed_keys}
        result['independent_retained_graph_check']=check(edges,removed)[0]
        result['retained_edges']=len(edges)-len(removed)
        results[label]=result
        print(label,result['solver_status'],result['num_edges_changed'],result['independent_retained_graph_check']['consistent'],flush=True)
        started=time.perf_counter()
        full=solve(reached,adjacency,strip,allowed_edge_keys=None,time_limit=120.)
        full['wall_seconds']=time.perf_counter()-started
        changed={(e['rel_pcl_id'],frozenset((e['from_point_id'],e['to_point_id']))) for e in full['edges']}
        removed={k for k,e in enumerate(edges) if edge_key(e) in changed}
        full['independent_retained_graph_check']=check(edges,removed)[0]
        unrestricted[label]=full
    target=[]
    for c in collections:
        if c['source_file']=='same_windings.json' and c['source_collection_id'] in ['149','150','151']:
            target.append({'collection':'same:'+c['source_collection_id'],'points':[
                {'id':p['id'],'xyz':p['p'],'hits':all_hits.get(id(p),[])} for p in c['points'].values()]})
    report={'schema_version':1,'band':[12417,12419],'patches':patch_report,'rejected_patches':rejected,'projection':'CPU villa Patch.project triangle fallback; largest area then nearest, tolerance 2.5, erosion 1 or metadata override; not native surface-index bit parity.',
        'band_filter':'Continuous quad intersection, not vertex-only ROI prefilter.',
        'linking_seconds':linking_seconds,'total_source_points':len(points),'attached_points':sum('on_patch' in p for p in points),
        'cross_collections':len(cross),'reached_patches':len(reached),'component_roots':component_roots,'measurable_edges':len(edges),'unmeasurable_edges':unmeasurable,'initial_graph':initial,
        'editable_edges':len(allowed),'solver_comparison':results,'unrestricted_solver_comparison':unrestricted,'target_traces':target,
        'real_data_improvement_demonstrated':results['after']['independent_retained_graph_check']['consistent'] and not results['before']['independent_retained_graph_check']['consistent']}
    for name,data in [('real_graph.json',report),('reached.json',reached),('rel_adjacency.json',adjacency),('measured_edges.json',edges)]:
        (outdir/name).write_text(json.dumps(data,indent=2,default=serial)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
