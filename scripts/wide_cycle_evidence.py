"""Extract real contradiction certificates and check proposed delta values."""
from collections import defaultdict,deque
import json
from pathlib import Path
import sys
from verify_wide import consistent

def main(directory):
    directory=Path(directory)
    edges=json.loads((directory/'measured_edges.json').read_text())
    adjacency=defaultdict(list)
    for i,e in enumerate(edges):
        adjacency[e['P']].append((e['R'],e['D'],i,1))
        adjacency[e['R']].append((e['P'],-e['D'],i,-1))
    potentials={};tree=defaultdict(list);tree_ids=set();components=0
    for start in sorted(adjacency):
        if start in potentials:continue
        components+=1;potentials[start]=0;queue=deque([start])
        while queue:
            a=queue.popleft()
            for b,d,i,sign in adjacency[a]:
                if b in potentials:continue
                potentials[b]=potentials[a]+d;queue.append(b);tree_ids.add(i)
                tree[a].append((b,d,i,sign));tree[b].append((a,-d,i,-sign))
    witnesses=[]
    for i,e in enumerate(edges):
        if potentials[e['R']]-potentials[e['P']]==e['D']:continue
        # Close the non-tree edge P -> R by the unique tree route R -> P.
        queue=deque([e['R']]);parent={e['R']:None}
        while e['P'] not in parent:
            a=queue.popleft()
            for b,d,j,sign in tree[a]:
                if b not in parent:parent[b]=(a,d,j,sign);queue.append(b)
        backwards=[];node=e['P']
        while parent[node] is not None:
            a,d,j,sign=parent[node];backwards.append({'edge_index':j,'direction':sign,'delta':d,'from':a,'to':node});node=a
        cycle=[{'edge_index':i,'direction':1,'delta':e['D'],'from':e['P'],'to':e['R']}]+list(reversed(backwards))
        assert all(a['to']==b['from'] for a,b in zip(cycle,cycle[1:]+cycle[:1]))
        residual=sum(c['delta'] for c in cycle);assert residual!=0
        xyz=[list(reversed(p)) for c in cycle for p in [edges[c['edge_index']]['from_zyx'],edges[c['edge_index']]['to_zyx']]]
        witnesses.append({'cycle':cycle,'integer_sum':residual,
            'endpoint_bbox_xyz':[[min(p[k] for p in xyz) for k in range(3)],[max(p[k] for p in xyz) for k in range(3)]],
            'source_edges':[{**edges[c['edge_index']],'edge_index':c['edge_index']} for c in cycle]})
    comparisons={}
    for label in ['before','after']:
        result=json.loads((directory/(label+'_full.json')).read_text())
        proposed={(x['rel_pcl_id'],frozenset([x['from_point_id'],x['to_point_id']])):x for x in result['edges']}
        corrected=[];mechanical=True
        for e in edges:
            x=dict(e);change=proposed.get((e['pcl_id'],frozenset([e['from_point_id'],e['to_point_id']])))
            if change:
                if change['action']!='change_winding_delta':mechanical=False
                else:
                    assert change['suggested_winding_delta']-change['current_winding_delta']==change['residual']
                    direct=(e['P'],e['R'])==(change['from_patch'],change['to_patch'])
                    assert direct or (e['R'],e['P'])==(change['from_patch'],change['to_patch'])
                    x['D']+=(1 if direct else -1)*change['residual']
            corrected.append(x)
        comparisons[label]={'retained_graph_consistent':consistent(edges,result['edges']),
                            'all_proposals_are_numeric_delta_edits':mechanical,
                            'full_graph_after_proposed_delta_edits_consistent':consistent(corrected,[]) if mechanical else None}
    output={'vertices_in_measured_graph':len(adjacency),'edges':len(edges),'components':components,
            'fundamental_cycles':len(edges)-len(tree_ids),'inconsistent_fundamental_cycles':len(witnesses),
            'contradiction_witnesses':witnesses,'independent_checks':comparisons,
            'observed_quality_improvement':comparisons['after']['retained_graph_consistent'] and not comparisons['before']['retained_graph_consistent'],
            'scope':'Real sampled surfaces and original annotations; no planted defects. Both solver proposals are checked against every measured edge.'}
    (directory/'cycle_evidence.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='contradiction_witnesses'},indent=2))

if __name__=='__main__':main(sys.argv[1])
