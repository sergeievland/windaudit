import gate
import numpy as np,json
tifffile=gate.tifffile
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order
out=[]
for folder in sorted(gate.PATCHES.iterdir()):
 a=np.stack([tifffile.imread(folder/(c+'.tif')) for c in 'xyz'],axis=-1).astype(np.float32);valid=np.any(a!=-1,axis=-1);v=valid[:-1,:-1]&valid[1:,:-1]&valid[:-1,1:]&valid[1:,1:];ii,jj=np.nonzero(v);n=len(ii);ids=np.full(v.shape,-1);ids[v]=np.arange(n);centres=np.stack([ii,jj],axis=1).astype(np.float32)+.5
 def theta(ij):
  xyz,bad=gate.lift(a,v,ij);assert bad==0;yx=xyz[:,1::-1]-gate.u(xyz[:,2]);return np.arctan2(yx[:,0],yx[:,1])%(2*np.pi)
 th=theta(centres);rows=[];cols=[];ds=[];invalid=0
 for di,dj in [(0,1),(1,0),(1,1),(1,-1)]:
  ni=ii+di;nj=jj+dj;ok=(ni>=0)&(nj>=0)&(ni<v.shape[0])&(nj<v.shape[1]);src=np.where(ok)[0];dst=ids[ni[ok],nj[ok]];ok=dst>=0;src=src[ok];dst=dst[ok]
  if di and dj:
   mid=(centres[src]+centres[dst])/2;mi,mj=np.floor(mid).astype(int).T;mv=v[mi,mj];invalid+=int((~mv).sum());d=th[dst]-th[src];delta=(d < -np.pi).astype(int)-(d>np.pi).astype(int);tm=theta(mid[mv]);d1=tm-th[src[mv]];d2=th[dst[mv]]-tm;delta[mv]=(d1 < -np.pi).astype(int)-(d1>np.pi).astype(int)+(d2 < -np.pi).astype(int)-(d2>np.pi).astype(int)
  else:
   d=th[dst]-th[src];delta=(d < -np.pi).astype(int)-(d>np.pi).astype(int)
  rows.extend(src);cols.extend(dst);ds.extend(delta)
 rows=np.array(rows);cols=np.array(cols);ds=np.array(ds);g=csr_matrix((np.ones(2*len(rows)),(np.r_[rows,cols],np.r_[cols,rows])),shape=(n,n));dg=csr_matrix((np.r_[ds,-ds],(np.r_[rows,cols],np.r_[cols,rows])),shape=(n,n));pot=np.full(n,999999,dtype=int)
 for s in range(n):
  if pot[s]!=999999:continue
  order,pred=breadth_first_order(g,s,directed=False);pot[s]=0
  for t in order[1:]:pot[t]=pot[pred[t]]+int(dg[pred[t],t])
 residual=pot[cols]-pot[rows]-ds;record={'patch':folder.name,'edges':len(rows),'inconsistent_edges':int(np.count_nonzero(residual)),'max_abs_residual':int(np.max(np.abs(residual))),'invalid_diagonal_midpoints':invalid};out.append(record);print(record,flush=True)
(gate.ROOT/'cycle_probe.json').write_text(json.dumps(out,indent=2))
