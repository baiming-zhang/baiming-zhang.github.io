"""Signed inverse condition number and exact coordinate derivatives for TET10."""
import numpy as np
IDEAL=np.array([[1.,.5,.5],[0.,np.sqrt(3)/2,np.sqrt(3)/6],[0.,0.,np.sqrt(2/3)]])
IDEAL_INV=np.linalg.inv(IDEAL)
EDGES=[(0,1),(1,2),(0,2),(0,3),(2,3),(1,3)]
def sicn_and_gradient(a):
    cofactor=np.stack([np.cross(a[...,1],a[...,2],axis=-1),np.cross(a[...,2],a[...,0],axis=-1),np.cross(a[...,0],a[...,1],axis=-1)],axis=-1) if False else np.stack([np.cross(a[..., :,1],a[..., :,2]),np.cross(a[..., :,2],a[..., :,0]),np.cross(a[..., :,0],a[..., :,1])],axis=-1)
    det=np.sum(a[..., :,0]*cofactor[..., :,0],axis=-1);s=np.maximum(np.sum(a*a,axis=(-2,-1)),1e-30);t=np.maximum(np.sum(cofactor*cofactor,axis=(-2,-1)),1e-30);den=np.sqrt(s*t);q=3*det/den;b=np.swapaxes(a,-1,-2)@a;g=3*cofactor/den[...,None,None]-q[...,None,None]*a/s[...,None,None]-q[...,None,None]*(s[...,None,None]*a-a@b)/t[...,None,None];return q,g,det,cofactor
def reference_derivatives(points):
    points=np.asarray(points);lam=np.c_[1-points.sum(1),points];grad=np.array([[-1.,-1.,-1.],[1,0,0],[0,1,0],[0,0,1]]);d=np.zeros((len(points),10,3));d[:,:4]=(4*lam[:,:,None]-1)*grad[None]
    for k,(a,b) in enumerate(EDGES,4):d[:,k]=4*(lam[:,a,None]*grad[b]+lam[:,b,None]*grad[a])
    return d
def objective_gradient(x,base,cells,free,dN,target,reference_det,bounds_scale,movement_weight=1e-5):
    coords=base.copy();delta=np.asarray(x).reshape(-1,3);coords[free]+=delta;gradient=np.zeros_like(coords);loss=0.;normalizer=len(cells)*len(dN)
    for first in range(0,len(cells),2000):
        c=cells[first:first+2000];jac=np.einsum('eni,qnj->eqij',coords[c],dN);a=jac@IDEAL_INV;q,dq,det,cof=sicn_and_gradient(a);bad=np.minimum(q-target,0.);loss+=float(np.sum(bad*bad))/normalizer;ga=(2*bad)[...,None,None]*dq/normalizer;scale=reference_det[first:first+len(c),None];violation=np.minimum(det/scale-.1,0.);loss+=10*float(np.sum(violation*violation))/normalizer;ga+=20*violation[...,None,None]*cof/scale[...,None,None]/normalizer;gj=ga@IDEAL_INV.T;gx=np.einsum('eqij,qnj->eni',gj,dN);np.add.at(gradient,c.ravel(),gx.reshape(-1,3))
    loss+=movement_weight*float(np.mean((delta/bounds_scale[:,None])**2));g=gradient[free]+2*movement_weight*delta/bounds_scale[:,None]**2/delta.size
    return loss,g.ravel()
