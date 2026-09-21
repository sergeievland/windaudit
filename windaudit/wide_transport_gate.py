"""Equivalent edge-integrability audit with bulk sparse parent-edge lookup."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order
from .scanspace import PatchGraph, theta_at

def audit(patch, umbilicus):
    graph=PatchGraph(patch);centres=graph.centres
    xyz,valid=patch.lift(centres);assert valid.all()
    theta=theta_at(xyz,umbilicus)
    row,col=graph.csr.nonzero();keep=row<col;row=row[keep];col=col[keep]
    d=theta[col]-theta[row]
    delta=(d < -np.pi).astype(int)-(d > np.pi).astype(int)
    indices=np.flatnonzero(np.all(centres[row]!=centres[col],axis=1))
    xyz,valid=patch.lift((centres[row[indices]]+centres[col[indices]])/2)
    k=indices[valid];tm=theta_at(xyz[valid],umbilicus)
    d1=tm-theta[row[k]];d2=theta[col[k]]-tm
    delta[k]=(d1 < -np.pi).astype(int)-(d1 > np.pi).astype(int)+(d2 < -np.pi).astype(int)-(d2 > np.pi).astype(int)
    dg=coo_matrix((np.r_[delta,-delta],(np.r_[row,col],np.r_[col,row])),shape=(graph.n,graph.n)).tocsr()
    sentinel=np.iinfo(np.int64).max;potential=np.full(graph.n,sentinel,dtype=np.int64)
    components=0
    for start in range(graph.n):
        if potential[start]!=sentinel:continue
        components+=1;order,parent=breadth_first_order(graph.csr,start,directed=False)
        potential[start]=0;nodes=order[1:]
        incoming=np.asarray(dg[parent[nodes],nodes]).ravel()
        for node,change in zip(nodes,incoming):potential[node]=potential[parent[node]]+int(change)
    residual=potential[col]-potential[row]-delta
    return {'valid_quads':graph.n,'components':components,'edges':len(row),
            'inconsistent_edges':int(np.count_nonzero(residual)),
            'max_abs_residual':int(np.max(np.abs(residual))) if len(residual) else 0,
            'invalid_diagonal_midpoints':int((~valid).sum())}
