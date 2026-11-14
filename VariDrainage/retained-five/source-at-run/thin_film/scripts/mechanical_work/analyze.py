"""Independent SI inertia, viscosity, clipped regional action and cut traction."""
from pathlib import Path
import importlib.util,json,numpy as np,shapely
from shapely import Polygon,LineString,box
from shapely.geometry.polygon import orient
import run as cfg
R,O,RAW,save,sha=cfg.R,cfg.O,cfg.RAW,cfg.save,cfg.sha
spec=importlib.util.spec_from_file_location('curved_saved_audit',Path(__file__).resolve().parents[1]/'curved_edge/analyze.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
CUT=240e-6;DOMAIN=box(0,0,320e-6,160e-6);RIGHT=box(CUT,0,320e-6,160e-6)

def viscosity(a,mid,faces,nodes,mu,dt,pyramid_observer=None):
    xy=np.c_[a['x_m'],a['y_m']];u=np.c_[mid['ux_m_s'],mid['uy_m_s']].ravel();n=len(xy);D=2*n;free=nodes['fixed']==0;bath=a['bath']>0;poly=old.polys(a);adj=[[] for _ in xy];eps=1e-10*np.hypot(320e-6,160e-6)
    for i,j in zip(faces['i'].astype(int),faces['j'].astype(int)):adj[i].append(j);adj[j].append(i)
    edges=[[] for _ in xy]
    for i,p in enumerate(poly):
        vertices=np.array(orient(p,sign=1).exterior.coords)[:-1]
        for A,B in zip(vertices,np.roll(vertices,-1,axis=0)):
            delta=B-A;length=np.linalg.norm(delta)
            if length<=eps:continue
            midpoint=.5*(A+B);normal=np.array([delta[1],-delta[0]])/length
            wall=any(abs(A[k]-c)<=eps and abs(B[k]-c)<=eps for k,c in [(0,0),(0,320e-6),(1,0),(1,160e-6)])
            if wall:j=-1
            else:
                js=np.array(adj[i]);d=xy[js]-xy[i];err=np.abs(np.sum((midpoint-.5*(xy[js]+xy[i]))*d,axis=1))/np.linalg.norm(d,axis=1);j=int(js[np.argmin(err)]);assert min(err)<=eps
            edges[i].append((j,midpoint,normal,length,A,B))
    def row(i,axis):
        b=np.zeros(D)
        if i>=0 and free[i]:b[2*i+axis]=1
        return b
    G=np.zeros((n,2,2,D))
    for i in range(n):
        if bath[i]:continue
        gram=np.zeros((2,2));samples=[]
        for j,midpoint,normal,length,A,B in edges[i]:
            fixed=j<0 or not free[j] or bath[j];d=midpoint-xy[i] if fixed else xy[j]-xy[i];w=length/np.linalg.norm(d);gram+=w*np.outer(d,d);samples.append((-1 if fixed else j,d,w))
        inverse=np.linalg.inv(gram)
        for j,d,w in samples:
            coefficient=inverse@(w*d)
            for axis in range(2):G[i,axis]+=np.outer(coefficient,row(j,axis)-row(i,axis))
    force=np.zeros(D);globalD=0.;rightD=0.;traction=0.;cut_length=0.;line=LineString([(CUT,0),(CUT,160e-6)])
    for i in range(n):
        if bath[i]:continue
        h=a['V_m3'][i]/a['A_m2'][i];density=2*mu*h
        for j,midpoint,normal,length,A,B in edges[i]:
            dr=midpoint-xy[i];distance=dr@normal;assert distance>0;H=G[i].copy()
            for axis in range(2):
                if j>=0 and free[j] and not bath[j]:
                    skew=midpoint-.5*(xy[i]+xy[j]);res=.5*(row(j,axis)-row(i,axis))+.5*(skew[0]*(G[i,axis,0]+G[j,axis,0])+skew[1]*(G[i,axis,1]+G[j,axis,1]))-dr[0]*G[i,axis,0]-dr[1]*G[i,axis,1]
                else:res=-row(i,axis)-dr[0]*G[i,axis,0]-dr[1]*G[i,axis,1]
                H[axis]+=np.outer(normal/distance,res)
            basis=[H[0,0],H[1,1],H[0,1]+H[1,0],H[0,0]+H[1,1]];factors=[1,1,.5,1];values=np.array([b@u for b in basis]);area=.5*length*distance
            for b,value,factor in zip(basis,values,factors):force+=density*area*factor*b*value
            strain=float(np.dot(factors,values**2));tri=Polygon([xy[i],A,B]);globalD+=dt*density*area*strain;rightD+=dt*density*tri.intersection(RIGHT).area*strain
            if pyramid_observer is not None:pyramid_observer(i,xy[i],u[2*i:2*i+2],H,u,tri,density,area,dt)
            xs=np.array([xy[i,0],A[0],B[0]])
            if np.any(xs<CUT) and np.any(xs>=CUT):
                segment=tri.intersection(line);length_cut=segment.length
                if length_cut>0:
                    yy=.5*(segment.bounds[1]+segment.bounds[3]);grad=np.einsum('abd,d->ab',H,u);velocity=u[2*i:2*i+2]+grad@(np.array([CUT,yy])-xy[i]);div=np.trace(grad);tauXX=density*(grad[0,0]+div);tauYX=.5*density*(grad[0,1]+grad[1,0]);traction-=dt*length_cut*(velocity[0]*tauXX+velocity[1]*tauYX);cut_length+=length_cut
    return force.reshape(n,2),globalD,rightD,traction,cut_length,poly

def relative(a,b):return float(np.max(np.abs(a-b))/max(1e-40,float(np.max(np.abs(b)))))
def main():
    status=json.loads((O/'status.json').read_text());assert status['state']=='completed';status.update(sampled_steps=True,native_checks=24);save(O/'status.json',status);old.O=O;old.RAW=RAW;old.main();save(O/'saved_formula_audit.json',json.loads((O/'summary.json').read_text()))
    audits={};diagnostics={};regional_rows={}
    for rec in status['cases']:
        if not rec['parameters']['observer'] or not rec['completed']:continue
        root=RAW/rec['name'];m=rec['summary'];dt=m['dt_s'];stats=old.read(root/'mechanics.csv');nodes_all=old.read(root/'mechanical_nodes.csv');region_all=old.read(root/'regional_mechanics.csv');states=old.read(root/'states.csv');faces_all=old.read(root/'faces.csv');ratios=stats['g1_N']/stats['target_N']
        diagnostics[rec['name']]=dict(accepted_steps=len(stats),above_native_force_target=int(sum(ratios>1.000001)),maximum_force_target_ratio=float(max(ratios)),maximum_Newton_iterations=int(max(stats['iters'])),maximum_CG_iterations=int(max(stats['cg_iters'])),negative_curvature_hits=int(sum(stats['negative_curvature'])),native_viscous_action_J=float(sum(stats['visc_J'])),native_inertial_BE_loss_J=float(sum(stats['inertia_loss_J'])))
        results=[];rows=[]
        for step in sorted(set(map(int,nodes_all['step']))):
            a=states[(states['step']==step)&(states['stage']==0)];mid=states[(states['step']==step)&(states['stage']==1)];nodes=nodes_all[nodes_all['step']==step];regions=region_all[region_all['step']==step];faces=faces_all[(faces_all['step']==step)&(faces_all['stage']==0)];native=stats[step-1];mask=nodes['fixed']==0;v0=np.c_[nodes['u0x_m_s'],nodes['u0y_m_s']];v1=np.c_[nodes['u1x_m_s'],nodes['u1y_m_s']];dx=dt*v1
            I=m['rho_kg_m3']*nodes['V_m3'][:,None]*(v1-v0)/dt;I[~mask]=0;V,D,DR,traction,cut_length,polys=viscosity(a,mid,faces,nodes,m['viscosity_Pa_s'],dt)
            Ic=np.c_[nodes['gix_N'],nodes['giy_N']];Vc=np.c_[nodes['gvx_N'],nodes['gvy_N']];C=np.c_[nodes['gcx_N'],nodes['gcy_N']];Cr=np.c_[nodes['grx_N'],nodes['gry_N']]
            w0=shapely.area(shapely.intersection(polys,RIGHT))/shapely.area(polys);p1=old.polys(mid);w1=shapely.area(shapely.intersection(p1,RIGHT))/shapely.area(p1)
            K0=.5*m['rho_kg_m3']*nodes['V_m3']*np.sum(v0*v0,axis=1);K1=.5*m['rho_kg_m3']*nodes['V_m3']*np.sum(v1*v1,axis=1);BE=.5*m['rho_kg_m3']*nodes['V_m3']*np.sum((v1-v0)**2,axis=1);residual=(I+V+C)*dx
            errors=dict(inertial_force=relative(I,Ic),viscous_force=relative(V,Vc),viscous_action=abs(D-native['visc_J'])/max(1e-40,abs(native['visc_J'])),cut_length=abs(cut_length/160e-6-1),weights=float(max(max(abs(w0-nodes['w0'])),max(abs(w1-nodes['w1'])))),regional=0.)
            for region in [0,1]:
                r=regions[region];wa=w0 if region else 1-w0;wb=w1 if region else 1-w1;clocal=Cr if region else C-Cr;inertial_work=float(np.sum(wa[:,None]*I*dx));viscous_work=float(np.sum(wa[:,None]*V*dx));capillary_work=float(np.sum(wa[:,None]*C*dx));caplocal=float(np.sum(clocal*dx));dE=float(np.sum(wb*nodes['e1_J']-wa*nodes['e0_J']));dK=float(np.sum(wb*K1-wa*K0));ownE=float(np.sum((wb-wa)*nodes['e1_J']));ownK=float(np.sum((wb-wa)*K1));bepart=float(np.sum(wa*BE));diss=DR if region else D-DR;caprem=caplocal-dE+ownE;capport=caplocal-capillary_work;viscport=diss-viscous_work;balance=dE+dK+bepart+diss+caprem-ownE-ownK-capport-viscport;force_work=float(np.sum(wa[:,None]*residual));scale=sum(abs(x) for x in [dE,dK,bepart,diss,caprem,ownE,ownK,capport,viscport])+1e-40
                predicted=[dE,dK,ownE,ownK,bepart,diss,inertial_work,viscous_work,capillary_work,caplocal,caprem,capport,viscport,balance];keys=['dE_J','dK_J','ownershipE_J','ownershipK_J','inertia_loss_J','viscous_D_J','inertial_work_J','viscous_node_J','capillary_node_J','capillary_local_J','capillary_remainder_J','capillary_port_J','viscous_port_J','residual_J'];errors['regional']=max(errors['regional'],max(abs(x-r[k]) for x,k in zip(predicted,keys))/scale,abs(balance-force_work)/scale)
                recovered=traction if region else -traction;errors['regional']=max(errors['regional'],abs(recovered-r['recovered_viscous_port_J'])/scale)
                rows.append(dict(step=step,t_s=step*dt,region=region,discrete_viscous_port_J=viscport,recovered_viscous_port_J=recovered,capillary_port_J=capport,viscous_dissipation_J=diss,endpoint_capillary_remainder_J=caprem,force_residual_work_J=force_work,balance_J=balance,scale_J=scale))
            errors['passed']=max(errors[k] for k in ['inertial_force','viscous_force','viscous_action','regional'])<=1e-7 and max(errors['cut_length'],errors['weights'])<=1e-8;results.append(dict(step=step,**errors))
        audits[rec['name']]=results;regional_rows[rec['name']]=rows
        print(rec['name'],'independent mechanical saved-step audits',sum(r['passed'] for r in results),'/',len(results),'force ratio max',max(ratios),flush=True)
    save(O/'mechanical_audit.json',audits);save(O/'mechanical_regions.json',regional_rows);save(O/'summary.json',dict(attempts=len(status['cases']),completed=sum(r['completed'] for r in status['cases']),native_checks=24,observer_diagnostics=diagnostics,independent_mechanical_audits=audits,full_C1_complete=False,scope='Sampled independent SI inertia/viscosity/clipped regional work and recovered traction; exact discrete identity is not continuum traction validation. Native dynamics unchanged.'))
if __name__=='__main__':main()
