"""Compare real graph evidence; measured wall-clock durations are volatile."""
import json
from pathlib import Path
import sys
def stable(x):
    if isinstance(x,dict):return {k:stable(v) for k,v in x.items() if k not in ['wall_seconds','linking_seconds']}
    if isinstance(x,list):return [stable(v) for v in x]
    return x
def main(a,b):
    for name in ['real_graph.json','band_transport_gate.json','reached.json','rel_adjacency.json','measured_edges.json']:
        x=(Path(a)/name).read_bytes();y=(Path(b)/name).read_bytes()
        if name=='real_graph.json':assert stable(json.loads(x))==stable(json.loads(y)),name
        else:assert x==y,name
        print(name+': reproduced')
if __name__=='__main__':main(*sys.argv[1:])
