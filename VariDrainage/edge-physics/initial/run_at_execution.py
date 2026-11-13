"""Independent SI static-preparation feasibility; NOT native moving validation.

Minimize constant-activity capillary grand potential on the directional
dual-face area already used by C1. Independent measured/model inputs only.
Fixed h traces at x=240/320um are initialization constraints, released in any
future moving assay. No time evolution or empirical drainage law here.
"""
from pathlib import Path
import hashlib,json,subprocess,time
import numpy as np
import shapely
from shapely import MultiPoint,box,STRtree
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parents[3]
O=R/'results/initial_edge_20261113'
W,Y,CUT,H0,RADIUS=320e-6,160e-6,240e-6,6e-6,500e-6

def geometry(mode,dx):
    points=[]
    for lo,hi,d in [(0,CUT,dx*(2 if mode=='local' else 1)),(CUT,W,dx)]:
        nx,ny=round((hi-lo)/d),round(Y/d)
        points.extend(((lo+(i+.5)*(hi-lo)/nx,(j+.5)*Y/ny) for j in range(ny) for i in range(nx)))
    for x in [CUT,W]:points.extend((x,(j+.5)*dx) for j in range(round(Y/dx)))
    xy=np.array(points);domain=box(0,0,W,Y)
    poly=shapely.intersection(np.array(shapely.get_parts(shapely.voronoi_polygons(MultiPoint(xy),extend_to=domain,ordered=True))),domain)
    i,j=STRtree(poly).query(poly,predicate='intersects');mask=i<j;i,j=i[mask],j[mask]
    length=shapely.length(shapely.intersection(shapely.boundary(poly[i]),shapely.boundary(poly[j])));mask=length>1e-15
    return xy,shapely.area(poly),i[mask],j[mask],length[mask]

def solve(mode,dx,gamma,viscosity):
    xy,a,i,j,length=geometry(mode,dx);d=np.linalg.norm(xy[i]-xy[j],axis=1)
    free=(xy[:,0]>CUT+1e-12)&(xy[:,0]<W-1e-12);bath=xy[:,0]>W-1e-12
    trace=np.abs(xy[:,0]-CUT)<1e-12
    far=H0+2*(W-CUT)**2/(RADIUS+np.sqrt(RADIUS**2-(W-CUT)**2));mu=-gamma/RADIUS
    h=np.full(len(a),H0);h[bath]=far
    # A linear interpolation is only Newton's start; no interior profile is held.
    h[free]=H0+(far-H0)*(xy[free,0]-CUT)/(W-CUT)
    ids=np.flatnonzero(free);scale=gamma*H0**2
    def fields(v,second=False):
        z=(v[i]-v[j])/(2*d);f=np.hypot(1,z);extra=length*d*z*z/(f+1)
        e=2*gamma*extra.sum()-mu*np.dot(a[free],v[free]);g=np.zeros(len(a));force=gamma*length*z/f
        np.add.at(g,i,force);np.add.at(g,j,-force);g[free]-=mu*a[free]
        if not second:return e,g
        k=gamma*length/(2*d*f**3);row=np.r_[i,j,i,j];col=np.r_[i,j,j,i]
        H=coo_matrix((np.r_[k,k,-k,-k],(row,col)),shape=(len(a),len(a))).tocsr()
        return e,g,H
    history=[]
    for iteration in range(50):
        e,g,H=fields(h,True);res=float(np.max(abs(g[free]/a[free]))/abs(mu));history.append(res)
        if res<1e-9:break
        step=spsolve(H[ids][:,ids],-g[free]);alpha=1.
        for _ in range(40):
            trial=h.copy();trial[free]+=alpha*step
            if min(trial)>0 and fields(trial)[0]<=e+1e-4*alpha*np.dot(g[free],step):h=trial;break
            alpha*=.5
        else:raise RuntimeError('static line search failed')
    e,g=fields(h);pressure=g/a;pressure[free]+=mu
    # After release of the initial trace: native-like liquid effort, bath fixed.
    pressure[bath]=mu
    hf=.5*(h[i]+h[j]);hf=np.where(bath[i]&~bath[j],h[j],np.where(bath[j]&~bath[i],h[i],hf))
    conductance=hf**3/(12*viscosity)*length/d;conductance[bath[i]&bath[j]]=0
    q=conductance*(pressure[i]-pressure[j]);edge=xy[:,0]>CUT+1e-12
    inward=float(np.sum((edge[j].astype(float)-edge[i].astype(float))*q))
    bath_gain=float(np.sum((bath[j].astype(float)-bath[i].astype(float))*q))
    name=f'{mode}{round(dx*1e6)}'
    np.savetxt(O/(name+'.csv'),np.c_[xy,a,h,pressure,free,trace,bath],delimiter=',',header='x_m,y_m,A_m2,h_m,pressure_Pa,free_initial,trace_initial,bath',comments='')
    result=dict(name=name,N=len(a),free_nodes=len(ids),iterations=iteration,pressure_residual=history[-1],history=history,grand_potential_J=e,initial_patch_V_m3=float(np.dot(a[edge&~bath],h[edge&~bath])),interface_into_edge_m3_s=inward,bath_gain_m3_s=bath_gain,minimum_h_m=float(min(h)),far_height_m=far,bath_pressure_Pa=mu,initial_trace_reaction_max_Pa=float(max(abs(pressure[trace]))),pressure_gate=history[-1]<1e-8)
    return result,xy,h,pressure,free

def main():
    O.mkdir(exist_ok=True);assert not (O/'protocol.json').exists();start=time.monotonic()
    old=json.loads((R/'results/edge_preparation_20261113/counter/status.json').read_text())
    m=next(c['summary'] for c in old['cases'] if c['name']=='c6_balanced_local')
    gamma=m['clean_tension_N_m']+m['RT_J_mol']*m['Gamma_inf_mol_m2']*np.log1p(-m['theta0'])
    protocol=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),gamma_N_m=float(gamma),viscosity_Pa_s=m['viscosity_Pa_s'],material_source='counter/status.json c6_balanced_local; unchanged SDS equilibrium activity',W_m=W,Y_m=Y,cut_m=CUT,h0_m=H0,radius_m=RADIUS,cases=['local20','uniform20','local10','uniform10'],pressure_gate=1e-8,scope='Independent static constant-activity preparation; regular sites and exact physical trace columns; no native dynamics, conserved-S preparation or convergence claim',inventory='New preparation has S=theta0*Gamma_inf*curved_area; not matched total S to flat-jump case',no_rate_or_exponent_fit=True)
    (O/'protocol.json').write_text(json.dumps(protocol,indent=2))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True);results=[]
    for mode,dx,color,ls in [('local',20e-6,'#3c5488','-'),('uniform',20e-6,'#e64b35','--'),('local',10e-6,'#3c5488',':'),('uniform',10e-6,'#e64b35','-.')]:
        result,xy,h,p,free=solve(mode,dx,gamma,m['viscosity_Pa_s']);results.append(result)
        # Central row is selected before observing solution, for profile display.
        band=np.abs(xy[:,1]-Y/2)<=dx/2+1e-12;xx=xy[band,0];order=np.argsort(xx)
        axs[0].plot(xx[order]*1e6,h[band][order]*1e6,ls,color=color,label=result['name'])
        axs[1].plot(np.arange(len(result['history'])),result['history'],ls,color=color,label=result['name'])
    axs[0].axvline(CUT*1e6,color='.6',lw=.8);axs[0].set(xlabel='x (µm)',ylabel='Initial thickness (µm)',title='Variational static edge, fixed physical traces')
    axs[0].legend(fontsize=8);axs[1].set(yscale='log',xlabel='Newton iteration',ylabel='Relative pressure residual',title='Static feasibility only; no moving validation');axs[1].axhline(1e-8,color='.6',lw=.8)
    for ext in ['png','pdf']:fig.savefig(O/('initial_edge.'+ext),dpi=200)
    plt.close(fig)
    summary=dict(cases=results,wall_s=time.monotonic()-start,all_pressure_gates=all(c['pressure_gate'] for c in results),warning='No time integration; pressure at released trace can launch a transient; bath gain near zero is equilibrium, not proven sustained drainage. Native variable-S mechanics, normal flux and regional power remain untested.')
    (O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
