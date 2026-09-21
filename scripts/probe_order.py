import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from windaudit.pcl import load_point_collections
from windaudit.frame import build_nodes,Umbilicus
from windaudit.calibrate import calibrate_from_absolute_anchors
from windaudit.order_constraints import *
p=Path('data/paris4'); r,s,a=[load_point_collections(str(p/f)) for f in ['relative_windings.json','same_windings.json','abs_winding.json']]; u=Umbilicus.load(str(p/'umbilicus.json'))
n=build_nodes(r,s,a,u); c=calibrate_from_absolute_anchors(n); pairs=candidates(n,c); match=build_links(n,c)
print('pairs',len(pairs),flush=True)
for m in [0.5,0.75,1.,1.25,1.5,2.,3.]:
 arcs,ev=build_constraints(n,c,m,match,pairs)
 cert=cycle_certificate(arcs)
 print(m,len(arcs),len({k for x in arcs for k in (x.u,x.v)}),cert,flush=True)
