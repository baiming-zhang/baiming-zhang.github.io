"""Solenoidal Legendre Galerkin reference for symmetric linear free-film modes."""
import numpy as np
from numpy.polynomial import legendre as lg
from scipy.linalg import eigh,expm

def build(N,h,rho,mu,gamma,E,k):
    x,weights=lg.leggauss(4*N+8);H=h/2;weights=weights*H
    phi=[];psi=[];dp=[];top=[];dtop=[];ddtop=[]
    for j in range(N):
        poly=lg.Legendre.basis(2*j);integ=poly.integ()
        phi.append(poly(x));psi.append(H*(integ(x)-integ(0)))
        dp.append(poly.deriv()(x)/H);top.append(poly(1))
        dtop.append(poly.deriv()(1)/H);ddtop.append(poly.deriv(2)(1)/H**2)
    phi,psi,dp=np.array(phi),np.array(psi),np.array(dp)
    M=rho*((phi*weights)@phi.T+k*k*(psi*weights)@psi.T)
    shear=dp+k*k*psi;D=mu*(4*k*k*(phi*weights)@phi.T+(shear*weights)@shear.T)
    lam,v=eigh(M);assert min(lam)>0;inv=(v/np.sqrt(lam))@v.T
    q=np.zeros(N);q[0]=h;top=np.array(top)
    B=-1j*k*np.array([q,top]) if E>0 else -1j*k*q[None,:]
    kval=np.array([gamma*k*k/2,2*E]) if E>0 else np.array([gamma*k*k/2])
    ks=np.sqrt(kval);C=ks[:,None]*B@inv;L=inv@D@inv
    A=np.block([[-L,-C.conj().T],[C,np.zeros((len(ks),len(ks)))]])
    return dict(A=A,L=L,M=M,D=D,inv=inv,ks=ks,B=B,q=q,top=top,
        dtop=np.array(dtop),ddtop=np.array(ddtop),N=N,h=h,rho=rho,mu=mu,gamma=gamma,E=E,k=k)

def observe(model,z):
    N=model['N'];a=model['inv']@z[:N];zdot=model['A']@z;adot=model['inv']@zdot[:N]
    y=z[N:]/model['ks'];u=model['top']@a;ut=model['top']@adot;k=model['k'];h=model['h'];mu=model['mu'];rho=model['rho'];eta=y[0];g=y[1]
    q=model['q']@a;w=-1j*k*q/2
    shear=mu*(model['dtop']@a+1j*k*w);expected=-1j*k*model['E']*g
    pressure=1j/k*(rho*ut-mu*(model['ddtop']@a-k*k*u));visc=2j*mu*k*u;cap=model['gamma']*k*k*eta/2
    return dict(eta=eta,g=g,u_surface=u,u_mean=q/h,flux=q,tangential_residual=shear-expected,tangential_scale=abs(shear)+abs(expected),normal_residual=pressure+visc-cap,normal_scale=abs(pressure)+abs(visc)+abs(cap),energy=.5*float(np.vdot(z,z).real),dissipation=float(np.vdot(z[:N],model['L']@z[:N]).real))

def trajectory(model,times,g0):
    z=np.zeros(model['N']+2,dtype=complex);z[-1]=model['ks'][-1]*g0;initial=z.copy();states=[];rows=[];before=0.
    for t in times:
        if t>before:z=expm(model['A']*(t-before))@z
        states.append(z.copy());rows.append(observe(model,z));before=t
    return initial,np.array(states),rows
