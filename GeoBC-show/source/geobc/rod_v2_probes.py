"""Exact TET10 field evaluation via physical inverse mapping, no nearest-node U."""
import numpy as np
from scipy.spatial import cKDTree

EDGES=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)]
DL=np.array([[-1.,-1.,-1.],[1,0,0],[0,1,0],[0,0,1]])


def shape(r):
    l=np.column_stack([1-r.sum(axis=1),r])
    N=np.empty((len(r),10));D=np.empty((len(r),10,3))
    N[:,:4]=l*(2*l-1);D[:,:4]=(4*l-1)[:,:,None]*DL[None,:,:]
    for k,(a,b) in enumerate(EDGES,4):
        N[:,k]=4*l[:,a]*l[:,b]
        D[:,k]=4*(l[:,a,None]*DL[b]+l[:,b,None]*DL[a])
    return N,D,l


def locate(xyz,cells,points):
    centers=xyz[cells[:,:4]].mean(axis=1);tree=cKDTree(centers)
    chosen=np.full(len(points),-1,dtype=int);natural=np.zeros((len(points),3));err=np.full(len(points),np.inf)
    for width in [16,64,256]:
        missing=np.flatnonzero(chosen<0)
        if not len(missing):break
        _,near=tree.query(points[missing],k=min(width,len(cells)))
        near=np.asarray(near).reshape(len(missing),-1)
        for j in range(near.shape[1]):
            todo=np.flatnonzero(chosen[missing]<0)
            if not len(todo):break
            source=missing[todo];element=near[todo,j];v=xyz[cells[element]];target=points[source]
            A=(v[:,1:4]-v[:,0,None,:]).transpose(0,2,1)
            r=np.linalg.solve(A,(target-v[:,0])[...,None])[...,0]
            # Clipped initial guesses avoid divergence on distant candidates only.
            r=np.clip(r,-.3,1.3)
            for _ in range(10):
                N,D,l=shape(r);value=np.einsum('ni,nij->nj',N,v);delta=value-target
                J=np.einsum('nik,nij->njk',D,v)
                try:step=np.linalg.solve(J,delta[...,None])[...,0]
                except np.linalg.LinAlgError:break
                r-=np.clip(step,-1.,1.)
                if np.max(np.abs(step))<1e-10:break
            N,D,l=shape(r);error=np.linalg.norm(np.einsum('ni,nij->nj',N,v)-target,axis=1)
            good=(l.min(axis=1)>=-1e-7)&(l.max(axis=1)<=1+1e-7)&(error<1e-7)
            chosen[source[good]]=element[good];natural[source[good]]=r[good];err[source[good]]=error[good]
    if np.any(chosen<0):
        raise ValueError(f'{(chosen<0).sum()} physical probes not found; no silent exclusion')
    N=shape(natural)[0]
    return chosen,N,float(err.max())


def evaluate(cells,U,element,N):
    return np.einsum('ni,nij->nj',N,U[cells[element]])


def selfcheck():
    # Curved manufactured isoparametric element; physical linear U is exact.
    ref=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1], [.5,0,0],[.5,.5,0],[0,.5,0],[0,0,.5],[0,.5,.5],[.5,0,.5]])
    xyz=ref.copy();xyz[4:,2]+=.05*ref[4:,0]
    cells=np.arange(10).reshape(1,10);r=np.array([[.2,.1,.3],[.1,.3,.2]])
    N=shape(r)[0];points=N@xyz;el,weights,error=locate(xyz,cells,points)
    U=np.column_stack([xyz[:,0]+2*xyz[:,1],xyz[:,2]*3,xyz[:,0]-xyz[:,2]])
    expected=np.column_stack([points[:,0]+2*points[:,1],3*points[:,2],points[:,0]-points[:,2]])
    result=evaluate(cells,U,el,weights)
    assert np.max(np.abs(result-expected))<1e-10
    return dict(passed=True,max_error=float(np.max(np.abs(result-expected))),mapping_residual_mm=error)
