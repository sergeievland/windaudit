"""Check saved wide inputs, retained graph consistency and repeat equality."""
from collections import defaultdict, deque
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]

# Which of several equally optimal edits a MILP returns depends on the solver
# build (HiGHS tie-breaking differs between SciPy 1.17.0 and 1.17.1 on this
# graph). Repeat comparisons therefore exclude the chosen edits and the
# potentials that follow from them, and check optimality and consistency of
# each run's own choice instead.
TIE_DEPENDENT={'linking_seconds','wall_seconds','edges','potentials','repair_difference'}

def canonical(value):
    if isinstance(value,dict):return {k:canonical(v) for k,v in value.items() if k not in TIE_DEPENDENT}
    if isinstance(value,list):return [canonical(v) for v in value]
    return value

def consistent(edges,changes):
    removed={(e['rel_pcl_id'],frozenset([e['from_point_id'],e['to_point_id']])) for e in changes}
    graph=defaultdict(list)
    for e in edges:
        if (e['pcl_id'],frozenset([e['from_point_id'],e['to_point_id']])) in removed:continue
        graph[e['P']].append((e['R'],e['D']));graph[e['R']].append((e['P'],-e['D']))
    seen={}
    for start in graph:
        if start in seen:continue
        seen[start]=0;queue=deque([start])
        while queue:
            a=queue.popleft()
            for b,delta in graph[a]:
                expected=seen[a]+delta
                if b in seen:
                    if seen[b]!=expected:return False
                else:seen[b]=expected;queue.append(b)
    return True

def main(saved,repeated=None):
    saved=Path(saved)
    manifest_path=saved/'dataset_manifest.json'
    if not manifest_path.exists():manifest_path=ROOT/'results/wide/dataset_manifest.json'  # both constructions use one dataset
    manifest=json.loads(manifest_path.read_text())
    total=0
    for patch in manifest['patches']:
        for filename,record in patch['files'].items():
            data=(ROOT/'data/scanspace_wide'/patch['id']/filename).read_bytes()
            assert len(data)==record['bytes']
            assert hashlib.sha256(data).hexdigest()==record['sha256']
            total+=len(data)
    assert total==manifest['total_bytes']<=manifest['announced_upper_bound']<=1_000_000_000
    edges=json.loads((saved/'measured_edges.json').read_text())
    report=json.loads((saved/'real_graph.json').read_text())
    assert consistent(edges,[])==report['initial_graph']['consistent']
    for label in ['before','after']:
        result=json.loads((saved/(label+'_full.json')).read_text())
        assert consistent(edges,result['edges'])==result['independent_retained_graph_check']['consistent']
        assert result['num_edges_changed']==len(result['edges'])
    assert report['solver_comparison']['after']['equations_in_model']==len(edges)
    print('Dataset hashes and budget verified; both retained graphs checked independently.')
    if repeated:
        names=['real_graph.json','reached.json','rel_adjacency.json','measured_edges.json','initial_graph.json','band_transport_gate.json','point_surface_hits.json','target_geometry.json','before_full.json','after_full.json']
        if (saved/'attachment_gaps.json').exists():names.append('attachment_gaps.json')
        for name in names:
            a=json.loads((saved/name).read_text());b=json.loads((Path(repeated)/name).read_text())
            assert canonical(a)==canonical(b),name
        repeat_edges=json.loads((Path(repeated)/'measured_edges.json').read_text())
        for label in ['before','after']:
            a=json.loads((saved/(label+'_full.json')).read_text());b=json.loads((Path(repeated)/(label+'_full.json')).read_text())
            assert a['num_edges_changed']==b['num_edges_changed'] and a['objective']==b['objective'],label
            assert consistent(repeat_edges,b['edges']),label
        print(f'Repeat: all {len(names)} compared JSON outputs identical after excluding runtimes and the choice among equally optimal edits;')
        print('each run\'s own edits are optimal and leave the graph consistent.')

if __name__=='__main__':main(*sys.argv[1:])
