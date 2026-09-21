"""CPU scan-space path-independence gate, fixed upstream sampling step 1."""
import json, pathlib, hashlib, time, sys
import numpy as np
from scipy import sparse
from scipy.ndimage import distance_transform_edt
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.interpolate import interp1d
from PIL import Image
class tifffile:
 @staticmethod
 def imread(path): return np.asarray(Image.open(path))
REPO=pathlib.Path(__file__).resolve().parents[2]
ROOT=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else REPO/'out/scanspace_gate'
ROOT.mkdir(parents=True,exist_ok=True)
PATCHES=REPO/'data/scanspace_probe'
p=json.loads((REPO/'data/paris4/umbilicus.json').read_text())['control_points'];p=sorted(p,key=lambda p:p['z'])
u=interp1d(np.array([q['z'] for q in p],np.float32),np.array([[q['y'],q['x']] for q in p],np.float32),axis=0,fill_value='extrapolate')
def sample(points):
 segs=[]
 for a,b in zip(points,points[1:]):
  n=max(1,int(np.ceil(float(np.linalg.norm(b-a)))))+1;t=np.linspace(0,1,n,dtype=np.float32)[:,None];s=a[None]*(1-t)+b[None]*t
  segs.append(s[1:] if segs else s)
 return np.concatenate(segs).astype(np.float32)
def lift(a,v,ij):
 i,j=np.floor(ij).astype(int).T;valid=(i>=0)&(j>=0)&(i<v.shape[0])&(j<v.shape[1]);valid[valid]&=v[i[valid],j[valid]];ij=ij[valid];i,j=np.floor(ij).astype(int).T;di,dj=(ij-np.floor(ij)).T;tl=a[i,j];tr=a[i,j+1];bl=a[i+1,j];br=a[i+1,j+1];top=tl+(tr-tl)*dj[:,None];bottom=bl+(br-bl)*dj[:,None];return top+(bottom-top)*di[:,None],int((~valid).sum())
def count(a,v,ij):
 xyz,bad=lift(a,v,sample(ij));yx=xyz[:,1::-1]-u(xyz[:,2]);theta=np.arctan2(yx[:,0],yx[:,1])%(2*np.pi);d=np.diff(theta);return int(np.sum(d < -np.pi)-np.sum(d>np.pi)),bad

def run(folder):
 a=np.stack([tifffile.imread(folder/(c+'.tif')) for c in 'xyz'],axis=-1).astype(np.float32);valid=np.any(a!=-1,axis=-1);v=valid[:-1,:-1]&valid[1:,:-1]&valid[:-1,1:]&valid[1:,1:];ii,jj=np.nonzero(v);n=len(ii);ids=np.full(v.shape,-1);ids[v]=np.arange(n);centres=np.stack([ii,jj],axis=1).astype(np.float32)+.5;dt=distance_transform_edt(np.pad(v,1))[1:-1,1:-1];rows=[];cols=[];base=[]
 for di in [-1,0,1]:
  for dj in [-1,0,1]:
   if not (di or dj):continue
   ni=ii+di;nj=jj+dj;ok=(ni>=0)&(nj>=0)&(ni<v.shape[0])&(nj<v.shape[1]);src=np.where(ok)[0];dst=ids[ni[ok],nj[ok]];ok=dst>=0;rows.extend(src[ok]);cols.extend(dst[ok]);base.extend([np.hypot(di,dj)]*int(ok.sum()))
 rows=np.array(rows);cols=np.array(cols);base=np.array(base);g=sparse.csr_matrix((base,(rows,cols)),shape=(n,n));nc,labels=connected_components(g);largest=np.where(labels==np.argmax(np.bincount(labels)))[0];rng=np.random.default_rng(20260920);starts=rng.choice(largest,min(20,len(largest)),replace=False);ends=rng.choice(largest,min(40,len(largest)),replace=False);paths={};records=[];bad_total=0
 for mode in range(5):
  weights=base*(1+[4,0,20,4,4][mode]/np.maximum(dt[ii[cols],jj[cols]],1))
  if mode>=3:weights*=np.exp(rng.uniform(-2,2,len(weights)))
  gg=sparse.csr_matrix((weights,(rows,cols)),shape=(n,n))
  for si,s in enumerate(starts):
   _,pred=dijkstra(gg,directed=True,indices=int(s),return_predecessors=True)
   for ti,t in enumerate(ends):
    if s==t:continue
    path=[int(t)]
    while path[-1]!=s:
     k=int(pred[path[-1]])
     if k<0:raise RuntimeError('unreachable')
     path.append(k)
    path=path[::-1];key=(si,ti);ij=centres[path];delta,bad=count(a,v,np.concatenate([ij[:1],ij,ij[-1:]]));bad_total+=bad
    if mode==0:paths[key]=path;records.append({'s':int(s),'t':int(t),'deltas':[delta],'distinct_routes':1,'loop_deltas':[]})
    else:
     rec=next(r for r in records if r['s']==s and r['t']==t);rec['deltas'].append(delta);rec['distinct_routes']+=int(path!=paths[key]);loop=path+paths[key][-2::-1];ld,lb=count(a,v,centres[loop]);rec['loop_deltas'].append(ld);bad_total+=lb
  print(folder.name,'mode',mode,'done',flush=True)
 mismatch=[r for r in records if len(set(r['deltas']))>1];loops=[x for r in records for x in r['loop_deltas']];result={'patch':folder.name,'shape':list(a.shape),'valid_quads':n,'components':int(nc),'pairs':len(records),'routes_per_pair':5,'distinct_alternative_routes':sum(r['distinct_routes']-1 for r in records),'disagreeing_pairs':len(mismatch),'pair_disagreement_fraction':len(mismatch)/len(records),'closed_contours':len(loops),'nonzero_closed_contours':sum(x!=0 for x in loops),'invalid_samples_dropped':bad_total,'examples':mismatch[:10], 'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir())}}
 (ROOT/(folder.name+'_gate.json')).write_text(json.dumps(result,indent=2));return result
if __name__=='__main__':
 results=[]
 for folder in sorted(PATCHES.iterdir()):
  results.append(run(folder));(ROOT/'gate_partial.json').write_text(json.dumps(results,indent=2))
