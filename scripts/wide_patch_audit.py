"""Build a CPU real-patch graph and compare the pinned upstream solver.

Usage: python scripts/wide_patch_audit.py PATCHES BEFORE.py AFTER.py OUTDIR
Only explicitly present patch folders are read. No network access.
"""
import ast
import contextlib
import hashlib
import io
import os
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
from windaudit.wide_transport_gate import audit as fast_gate

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

def main(directory,before,after,outdir,inputs=None,zlo=6000,zhi=18000):
    zlo,zhi=int(zlo),int(zhi)
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    inputs=Path(inputs) if inputs is not None else ROOT/'data/paris4'
    umbilicus=load_umbilicus(inputs/'umbilicus.json')
    cache=ROOT/'out/wide_cache';cache.mkdir(parents=True,exist_ok=True)
    code_files=[Path(__file__),ROOT/'windaudit/scanspace.py',ROOT/'windaudit/scanspace_surface.py',ROOT/'windaudit/wide_transport_gate.py']
    prefix=hashlib.sha256(b''.join(p.read_bytes() for p in code_files)+b''.join((inputs/n).read_bytes() for n in ['umbilicus.json','abs_winding.json','relative_windings.json','same_windings.json'])).hexdigest()
    fingerprints={}
    patches={};surfaces={};areas={};patch_report=[];rejected=[]
    for folder in sorted(Path(directory).iterdir()):
        if not folder.is_dir():continue
        patch=ScanPatch.load(folder);valid=np.any(patch.xyz!=-1,axis=-1)
        cells=int(patch.metadata.get('spiral_patch_erode_cells',1))
        if cells>0:patch.xyz[~binary_erosion(valid,iterations=cells,border_value=0)]=-1;patch.__post_init__()
        z=patch.xyz[...,2];in_band=(z>=zlo)&(z<zhi)&np.any(patch.xyz!=-1,axis=-1)
        # Select by actual continuous quad extent as well as vertex hits: a
        # two-voxel window can fall between grid rows. Record this explicitly.
        v=patch.valid_quad
        corners=np.stack([z[:-1,:-1],z[1:,:-1],z[:-1,1:],z[1:,1:]])
        crosses=v&(corners.min(axis=0)<zhi)&(corners.max(axis=0)>=zlo)
        if not crosses.any():
            rejected.append({'id':folder.name,'reason':'No valid quad intersects the band after mask handling and erosion.'})
            continue
        patches[folder.name]=patch
        fingerprints[folder.name]=hashlib.sha256(prefix.encode()+patch.xyz.tobytes()).hexdigest()
        areas[folder.name]=float(v.sum()/np.prod(patch.metadata['scale']))
        patch_report.append({'id':folder.name,'valid_quads':int(v.sum()),'band_crossing_quads':int(crosses.sum()),'has_vertex_in_band':bool(in_band.any()),'erosion_cells':cells})
    print('loaded surfaces',len(patches),flush=True)
    gate=[]
    for name,patch in patches.items():
        cached=cache/(fingerprints[name]+'.gate.json')
        if cached.exists() and not os.environ.get('WIDE_FORCE_RECOMPUTE'):r=json.loads(cached.read_text())
        else:
            r={'patch':name,**fast_gate(patch,umbilicus)}
            cached.write_text(json.dumps(r))
        gate.append(r)
        print('edge gate',name,r['inconsistent_edges'],flush=True)
        (outdir/'band_transport_gate.json').write_text(json.dumps(gate,indent=2)+'\n')
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
    point_xyz=np.asarray([p['p'] for p in points],np.float32)
    for n,(name,patch) in enumerate(patches.items()):
        cached=cache/(fingerprints[name]+'.hits.json')
        if cached.exists() and not os.environ.get('WIDE_FORCE_RECOMPUTE'):patch_hits=json.loads(cached.read_text())
        else:
            patch_hits=[];surface=Surface(patch)
            lower=surface.triangles.min(axis=(0,1));upper=surface.triangles.max(axis=(0,1))
            candidates=np.flatnonzero(np.all((point_xyz>=lower-2.501)&(point_xyz<=upper+2.501),axis=1))
            for pi in candidates:
                hit=surface.project(points[pi]['p'],2.5)
                if hit is not None:patch_hits.append([int(pi),{**hit,'id':name,'area':areas[name]}])
            del surface
            cached.write_text(json.dumps(patch_hits,default=serial))
        for pi,hit in patch_hits:all_hits.setdefault(id(points[pi]),[]).append(hit)
        if (n+1)%25==0:print('projected surfaces',n+1,'attached',len(all_hits),flush=True)
    for p in points:
        hits=all_hits.get(id(p),[])
        if hits:
            hits.sort(key=lambda h:(-h['area'],h['distance'],h['id']))
            p['on_patch']=hits[0];all_hits[id(p)]=hits
    linking_seconds=time.perf_counter()-start
    point_inventory=[{'source_file':c['source_file'],'collection_id':c['source_collection_id'],
        'points':[{'id':p['id'],'xyz':p['p'],'hits':all_hits.get(id(p),[])} for p in c['points'].values()]} for c in collections]
    (outdir/'point_surface_hits.json').write_text(json.dumps(point_inventory,indent=2,default=serial)+'\n')
    # Triangle indexes are no longer needed for graph construction.
    surfaces.clear()
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
            if pid not in graphs:
                if len(graphs)>=8:del graphs[next(iter(graphs))]
                graphs[pid]=PatchGraph(patches[pid])
            r=strip_winding_delta(graphs[pid],a,b,umbilicus)
            memo[key]=None if r is None else r['delta_windings']
            if len(graphs[pid].cache)>8:
                del graphs[pid].cache[next(iter(graphs[pid].cache))]
        return memo[key]
    def bfs():
        reached={};roots=[]
        for seed in sorted(adjacency):
            if seed in reached:continue
            roots.append(seed);entry=adjacency[seed][0]['from_ij']
            reached[seed]={'acc':0,'entry_ij':entry,'hops':0};q=deque([seed])
            while q:
                a=q.popleft();state=reached[a]
                for e in adjacency[a]:
                    b=e['neighbor']
                    if b in reached:continue
                    d=strip(a,state['entry_ij'],e['from_ij'])
                    if d is None:continue
                    reached[b]={'acc':state['acc']-d-e['winding_delta'],'entry_ij':e['to_ij'],'hops':state['hops']+1};q.append(b)
        return reached,roots
    reached,component_roots=bfs()
    # Optional attachment-gap transport (WIDE_ATTACHMENT_GAP=1). Upstream counts
    # branch crossings along the annotation tour at the annotation points' own
    # coordinates, but each within-patch strip ends at the point's attachment on
    # the surface -- precisely, at the strip's last valid sample. The segment
    # between that sample and the annotation point is transported by neither,
    # so a point lying across the branch ray from its attachment shifts the edge
    # by one winding. Here that segment is transported and added to the edge's
    # winding delta, in both directions. BFS entries do not depend on winding
    # values, so the tree is unchanged; accumulated windings are recomputed.
    gap_records=[]
    gap_enabled=os.environ.get('WIDE_ATTACHMENT_GAP')=='1'
    if gap_enabled:
        from windaudit.scanspace import crossing_delta, polyline_ijs
        ends={}
        def strip_end(pid,a,b):
            key=(pid,tuple(a),tuple(b))
            if key not in ends:
                if pid not in graphs:
                    if len(graphs)>=8:del graphs[next(iter(graphs))]
                    graphs[pid]=PatchGraph(patches[pid])
                centres=graphs[pid].route(a,b)
                if centres is None:ends[key]=None
                else:
                    xyz,valid=patches[pid].lift(polyline_ijs(np.vstack([a,centres,b]),1.0))
                    ends[key]=xyz[valid][-1] if valid.sum()>=2 else None
            return ends[key]
        for a,links in adjacency.items():
            if a not in reached:continue
            for e in links:
                b=e['neighbor']
                if b not in reached:continue
                fa=strip_end(a,reached[a]['entry_ij'],e['from_ij']);fb=strip_end(b,reached[b]['entry_ij'],e['to_ij'])
                if fa is None or fb is None:continue
                pa=np.asarray(e['from_zyx'],np.float64)[::-1];pb=np.asarray(e['to_zyx'],np.float64)[::-1]
                ga=crossing_delta(theta_at(np.stack([fa,pa]),umbilicus))
                gb=crossing_delta(theta_at(np.stack([pb,fb]),umbilicus))
                e['attachment_gap_delta']=int(ga+gb)
                if ga or gb:
                    gap_records.append({'pcl_id':e['pcl_id'],'pcl_name':e.get('pcl_name'),'from_patch':a,'to_patch':b,
                        'from_point_id':e['from_point_id'],'to_point_id':e['to_point_id'],'gap_from':int(ga),'gap_to':int(gb),
                        'from_point_xyz':pa.tolist(),'from_strip_end_xyz':np.asarray(fa,np.float64).tolist(),
                        'to_point_xyz':pb.tolist(),'to_strip_end_xyz':np.asarray(fb,np.float64).tolist(),
                        'winding_delta_uncorrected':int(e['winding_delta'])})
                e['winding_delta']=int(e['winding_delta']+ga+gb)
        reached,component_roots=bfs()
    (outdir/'attachment_gaps.json').write_text(json.dumps({'enabled':gap_enabled,
        'directed_edges_corrected':len(gap_records),'records':gap_records},indent=2,default=serial)+'\n')
    edges=[];seen=set();unmeasurable=0
    for a,links in adjacency.items():
        for e in links:
            b=e['neighbor'];key=edge_key(e)
            if key in seen or a not in reached or b not in reached:continue
            seen.add(key);sa=strip(a,reached[a]['entry_ij'],e['from_ij']);sb=strip(b,reached[b]['entry_ij'],e['to_ij'])
            if sa is None or sb is None:unmeasurable+=1;continue
            edges.append({**e,'P':a,'R':b,'D':int(e['winding_delta']+sa-sb)})
    initial,allowed_indices=check(edges);allowed={edge_key(edges[k]) for k in allowed_indices}
    for name,data in [('reached.json',reached),('rel_adjacency.json',adjacency),('measured_edges.json',edges),('initial_graph.json',initial)]:
        (outdir/name).write_text(json.dumps(data,indent=2,default=serial)+'\n')
    print('GRAPH',len(reached),len(edges),initial['inconsistent_non_tree_edges'],flush=True)
    results={}
    unrestricted={}
    for label,path in [('before',before),('after',after)]:
        solve=extract(path,['solve_min_edge_fix'])['solve_min_edge_fix']
        capture=io.StringIO();started=time.perf_counter()
        with contextlib.redirect_stdout(capture),contextlib.redirect_stderr(capture):
            result=solve(reached,adjacency,strip,allowed_edge_keys=allowed,time_limit=120.)
        result['wall_seconds']=time.perf_counter()-started
        removed_keys={(e['rel_pcl_id'],frozenset((e['from_point_id'],e['to_point_id']))) for e in result['edges']}
        removed={k for k,e in enumerate(edges) if edge_key(e) in removed_keys}
        result['independent_retained_graph_check']=check(edges,removed)[0]
        result['retained_edges']=len(edges)-len(removed)
        result['equations_in_model']=result['num_edges_considered']
        result['milp_inequality_rows']=2*result['num_edges_considered']
        (outdir/(label+'_full.json')).write_text(json.dumps(result,indent=2,default=serial)+'\n')
        (outdir/(label+'_stdout.txt')).write_text(capture.getvalue()+json.dumps(result,indent=2,default=serial)+'\n')
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
    geometry=[]
    target_by_id={t['collection']:t for t in target}
    for left,right in [('same:149','same:150'),('same:150','same:151'),('same:149','same:151')]:
        witnesses=[]
        for a in target_by_id[left]['points']:
            for b in target_by_id[right]['points']:
                for ha in a['hits']:
                    for hb in b['hits']:
                        if ha['id']!=hb['id']:continue
                        d=strip(ha['id'],ha['ij'],hb['ij'])
                        witnesses.append({'patch':ha['id'],'left_point':a['id'],'right_point':b['id'],
                            'left_ij':ha['ij'],'right_ij':hb['ij'],'left_distance':ha['distance'],'right_distance':hb['distance'],'transport':d})
        connected=[w for w in witnesses if w['transport'] is not None]
        values=sorted({w['transport'] for w in connected})
        verdict='unresolved'
        if values==[0]:verdict='same_local_winding_supported_on_shared_connected_surfaces'
        elif values and 0 not in values:verdict='different_local_windings_supported_on_shared_connected_surfaces'
        elif values:verdict='ambiguous_surface_attachments'
        geometry.append({'left':left,'right':right,'verdict':verdict,'transport_values':values,'witnesses':witnesses,
                         'connected_left_points':len({w['left_point'] for w in connected}),
                         'connected_right_points':len({w['right_point'] for w in connected}),
                         'left_total_points':len(target_by_id[left]['points']),
                         'right_total_points':len(target_by_id[right]['points']),
                         'scope':'Reconstructed surface evidence at fixed 2.5-voxel point attachment tolerance; no annotation equality edges used for these paths.'})
    (outdir/'target_geometry.json').write_text(json.dumps(geometry,indent=2,default=serial)+'\n')
    before_keys={(e['rel_pcl_id'],e['from_point_id'],e['to_point_id'],e['residual']) for e in results['before']['edges']}
    after_keys={(e['rel_pcl_id'],e['from_point_id'],e['to_point_id'],e['residual']) for e in results['after']['edges']}
    report={'schema_version':1,'band':[zlo,zhi],'patches':patch_report,'rejected_patches':rejected,'projection':'CPU villa Patch.project triangle fallback; largest area then nearest, tolerance 2.5, erosion 1 or metadata override; not native surface-index bit parity.',
        'band_filter':'Continuous quad intersection, not vertex-only ROI prefilter.',
        'linking_seconds':linking_seconds,'total_source_points':len(points),'attached_points':sum('on_patch' in p for p in points),
        'cross_collections':len(cross),'reached_patches':len(reached),'component_roots':component_roots,'measurable_edges':len(edges),'unmeasurable_edges':unmeasurable,'initial_graph':initial,
        'editable_edges':len(allowed),'solver_comparison':results,'unrestricted_solver_comparison':unrestricted,'target_traces':target,
        'repair_difference':{'only_before':sorted(before_keys-after_keys),'only_after':sorted(after_keys-before_keys)},
        'behaviour':{'real_equations':len(edges),'before_equations':results['before']['num_edges_considered'],'after_equations':results['after']['num_edges_considered']},
        'real_data_improvement_demonstrated':results['after']['independent_retained_graph_check']['consistent'] and not results['before']['independent_retained_graph_check']['consistent']}
    for name,data in [('real_graph.json',report),('reached.json',reached),('rel_adjacency.json',adjacency),('measured_edges.json',edges)]:
        (outdir/name).write_text(json.dumps(data,indent=2,default=serial)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
