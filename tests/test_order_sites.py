from windaudit.order_constraints import Arc,repair
from windaudit.order_controls import sites
from windaudit.pcl import Collection,Point

def test_directed_paths_not_undirected_connectivity():
    a,b,c=('relative',1),('relative',2),('relative',3)
    col=Collection(id=1,name='target',points={i:Point(id=i,xyz=[0,0,0],wind_a=i) for i in range(1,5)})
    ev=[(a,i,b,i,'ab') for i in [1,2]]+[(a,i,c,i,'ac') for i in [3,4]]
    arcs=[Arc(a,b,0,'ab'),Arc(a,c,0,'ac'),Arc(b,c,0,'bc')]
    assert sites(arcs,ev,{1:col},'relative',2)[0]==[]
    arcs=[Arc(a,b,0,'ab'),Arc(c,a,0,'ac'),Arc(b,c,0,'bc')]
    ss,total=sites(arcs,ev,{1:col},'relative',2)
    assert total==3 and [(x[0],x[1]) for x in ss]==[(1,2)]

def test_duplicate_support_cannot_buy_repair_votes():
    a,b,c=('s',1),('s',2),('s',3)
    arcs=[Arc(a,b,-1,'ab'),Arc(b,c,-1,'bc'),Arc(c,a,-1,'ca')]
    original=repair(arcs,method='cycles')
    duplicated=repair(arcs+[arcs[0]]*40,method='cycles')
    assert original==duplicated
