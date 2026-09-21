"""Independent Floyd-Warshall and exhaustive repair references."""
import itertools
import random
import numpy as np
import pytest
from windaudit.order_constraints import Arc, negative_cycle, repair, minimum_edge_cut

def feasible(arcs):
    nodes=sorted({x for a in arcs for x in (a.u,a.v)})
    d={(u,v):0 if u==v else float('inf') for u in nodes for v in nodes}
    for a in arcs:d[a.u,a.v]=min(d[a.u,a.v],a.c)
    for k in nodes:
        for i in nodes:
            for j in nodes:d[i,j]=min(d[i,j],d[i,k]+d[k,j])
    return all(d[k,k]>=0 for k in nodes)

def reference(arcs):
    nodes=sorted({x for a in arcs for x in (a.u,a.v)})
    best=None; opts=[]
    for bits in itertools.product([0,1],repeat=len(nodes)):
        removed={k for k,b in zip(nodes,bits) if b}
        if not feasible([a for a in arcs if a.u not in removed and a.v not in removed]):continue
        lost={a.group for a in arcs if a.u in removed or a.v in removed}
        obj=(len(removed),len(lost))
        if best is None or obj<best:best=obj;opts=[removed]
        elif obj==best:opts.append(removed)
    return best,opts

def graphs():
    rng=random.Random(2819)
    for _ in range(35):
        nodes=[('same',i) for i in range(rng.randint(2,5))]
        arcs=[]
        for i in range(rng.randint(2,9)):
            u,v=rng.sample(nodes,2);c=rng.randint(-3,3)
            arcs.append(Arc(u,v,c,str(i)))
            if rng.random()<.3:arcs.append(Arc(v,u,-c,str(i)))
        yield arcs

@pytest.mark.parametrize('method',['big_m','cycles'])
@pytest.mark.parametrize('arcs',list(graphs()))
def test_bellman_ford_and_repairs_against_independent_references(arcs,method):
    assert bool(negative_cycle(arcs)) == (not feasible(arcs))
    result=repair(arcs,method=method)
    assert result['optimal']
    if feasible(arcs):return
    best,opts=reference(arcs)
    assert (result['objective'],result['groups_lost'])==best
    for k in {x for a in arcs for x in (a.u,a.v)}:
        counts=sum(k in s for s in opts)
        expected='certain' if counts==len(opts) else 'possible' if counts else 'excluded'
        assert result['membership'][f'{k[0]}:{k[1]}']==expected
    cut=repair(arcs,edge_cut=True,method=method)
    gs=sorted({a.group for a in arcs})
    optimum=min(sum(bs) for bs in itertools.product([0,1],repeat=len(gs))
                if feasible([a for a in arcs if not dict(zip(gs,bs))[a.group]]))
    assert cut['optimal'] and cut['objective']==optimum
    decomposed=minimum_edge_cut(arcs)
    assert decomposed['optimal'] and decomposed['objective']==optimum
    assert decomposed['cut_groups']==cut['cut_groups']

def test_self_loop_and_empty():
    a=Arc(('same',1),('same',1),-1,'x')
    assert negative_cycle([a])==[0]
    assert repair([a])['membership']=={'same:1':'certain'}
    assert not negative_cycle([])

def test_integer_bounds_required():
    with pytest.raises(ValueError):negative_cycle([Arc(('s',1),('s',2),.5,'x')])

def test_expired_budget_never_claims_a_repair():
    a=Arc(('same',1),('same',1),-1,'x')
    result=repair([a],method='cycles',time_budget=0.)
    assert not result['optimal']
    assert result['membership']=={'same:1':'unknown'}
