"""CPU projection onto actual patch faces, with conservative spatial pruning.

Uses villa Patch.project's two triangles per valid quad and its barycentric
region tests. This is the documented brute-force fallback geometry, not a
claim of bit equivalence to the optional vc_spiral native surface index.
"""
import numpy as np
from scipy.spatial import cKDTree


class Surface:
    def __init__(self, patch):
        i,j=np.nonzero(patch.valid_quad)
        xyz=patch.xyz
        self.triangles=np.concatenate([
            np.stack([xyz[i+1,j],xyz[i,j],xyz[i,j+1]],axis=1),
            np.stack([xyz[i+1,j],xyz[i,j+1],xyz[i+1,j+1]],axis=1)])
        self.bases=np.tile(np.column_stack([i,j]),(2,1))
        self.types=np.repeat([0,1],len(i))
        self.centres=self.triangles.mean(axis=1)
        self.radii=np.linalg.norm(self.triangles-self.centres[:,None],axis=2).max(axis=1)
        self.tree=cKDTree(self.centres)
        self.max_radius=float(self.radii.max())

    def project(self, point, tolerance=2.5):
        point=np.asarray(point,dtype=np.float32)
        ids=np.asarray(sorted(self.tree.query_ball_point(point,self.max_radius+tolerance+1e-3)),dtype=int)
        if not len(ids):return None
        ok=np.linalg.norm(self.centres[ids]-point,axis=1)<=self.radii[ids]+tolerance+1e-3
        ids=ids[ok]
        if not len(ids):return None
        a,b,c=self.triangles[ids].transpose(1,0,2);ab=b-a;ac=c-a
        ap=point-a;bp=point-b;cp=point-c
        dot=lambda x,y:np.sum(x*y,axis=-1)
        d1=dot(ap,ab);d2=dot(ap,ac);d3=dot(bp,ab);d4=dot(bp,ac);d5=dot(cp,ab);d6=dot(cp,ac)
        vc=d1*d4-d3*d2;vb=d5*d2-d1*d6;va=d3*d6-d5*d4
        ma=(d1<=0)&(d2<=0);mb=(d3>=0)&(d4<=d3);mab=(vc<=0)&(d1>=0)&(d3<=0)
        mc=(d6>=0)&(d5<=d6);mac=(vb<=0)&(d2>=0)&(d6<=0);mbc=(va<=0)&((d4-d3)>=0)&((d5-d6)>=0)
        mf=~(ma|mb|mab|mc|mac|mbc);bary=np.zeros((len(ids),3),np.float32);eps=1e-8
        bary[ma,0]=1;bary[mb,1]=1
        v=d1/(d1-d3+eps);bary[mab,0]=1-v[mab];bary[mab,1]=v[mab]
        bary[mc,2]=1
        w=d2/(d2-d6+eps);bary[mac,0]=1-w[mac];bary[mac,2]=w[mac]
        w=(d4-d3)/((d4-d3)+(d5-d6)+eps);bary[mbc,1]=1-w[mbc];bary[mbc,2]=w[mbc]
        denom=va+vb+vc+eps;v=vb/denom;w=vc/denom
        bary[mf,1]=v[mf];bary[mf,2]=w[mf];bary[mf,0]=1-v[mf]-w[mf]
        nearest=bary[:,0,None]*a+bary[:,1,None]*b+bary[:,2,None]*c
        dist=np.linalg.norm(nearest-point,axis=1);k=int(np.argmin(dist))
        if dist[k]>tolerance:return None
        face=ids[k];u,v,w=bary[k]
        offset=[u,w] if self.types[face]==0 else [u+w,v+w]
        ij=self.bases[face]+offset
        return {'ij':np.asarray(ij,np.float32).tolist(),'distance':float(dist[k]),
                'face_id':int(face)}
