"""Frozen-state cut-work diagnosis: bilinear commutator and broken Green identity.

No native simulation, material, trajectory or previous result file is changed.
"""
from pathlib import Path
import sys,json,hashlib,subprocess,time
import numpy as np
from shapely import box,Polygon
from shapely.geometry.polygon import orient
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'mechanical_work'))
import analyze as prior
R=prior.R;SOURCE=prior.O;RAW=prior.RAW;O=R/'results/cut_work_20261113'
CUT=prior.CUT;LEFT=box(0,0,CUT,160e-6)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,data):p.write_text(json.dumps(data,indent=2),encoding='utf-8')
def strain(H):return np.array([H[0,0],H[1,1],H[0,1]+H[1,0],np.trace(H)])
FACTOR=np.array([1.,1.,.5,1.])
def stress(H,density):return density*np.array([[H[0,0]+np.trace(H),.5*(H[0,1]+H[1,0])],[.5*(H[0,1]+H[1,0]),H[1,1]+np.trace(H)]])
def boundary_power(poly,xi,ui,H,density,dt):
    """Exact midpoint edge quadrature: affine velocity, constant stress."""
    totals=dict(cut=0.,wall=0.,internal=0.)
    if poly.is_empty or poly.area==0:return totals
    vertices=np.asarray(orient(poly,sign=1).exterior.coords)[:-1];sigma=stress(H,density);eps=2e-14
    for A,B in zip(vertices,np.roll(vertices,-1,axis=0)):
        edge=B-A;length=np.linalg.norm(edge)
        if length==0:continue
        normal=np.array([edge[1],-edge[0]])/length;mid=.5*(A+B);velocity=ui+H@(mid-xi);power=dt*length*velocity@sigma@normal
        if abs(A[0]-CUT)<eps and abs(B[0]-CUT)<eps:kind='cut'
        elif any(abs(A[k]-c)<eps and abs(B[k]-c)<eps for k,c in [(0,0),(0,320e-6),(1,0),(1,160e-6)]):kind='wall'
        else:kind='internal'
        totals[kind]+=power
    return totals
def affine_guards():
    rng=np.random.default_rng(11303);errors=[]
    for k in range(48):
        # Rotated random affine field; independent polygon integral vs volume.
        xy=np.array([[.2,.1],[1.2,.25],[.65,1.3]])*1e-4;poly=Polygon(xy).intersection(box(0,0,240e-6,160e-6));H=rng.normal(size=(2,2))*20;ui=rng.normal(size=2)*.001;density=.003*6e-6;dt=.5e-6
        b=boundary_power(poly,xy[0],ui,H,density,dt);expected=dt*density*poly.area*np.dot(FACTOR,strain(H)**2);actual=sum(b.values());errors.append(abs(actual-expected)/max(1e-40,abs(expected)))
    assert max(errors)<1e-10
    return dict(cases=len(errors),maximum_relative_error=max(errors),purpose='Independent divergence-theorem quadrature/constitutive-factor guard; not native affine reproduction with nonzero wall data.')
def main():
    start=time.monotonic();O.mkdir(exist_ok=True);assert not (O/'protocol.json').exists(),'Frozen campaign exists; do not overwrite.'
    status=json.loads((SOURCE/'status.json').read_text());assert status['state']=='completed'
    manifest=json.loads((SOURCE/'raw_manifest.json').read_text())
    paths=[HERE/'analyze.py',HERE.parent/'mechanical_work/analyze.py',R/'thin_film/include/vf/thinfilm/MovingVoronoiFilm.h']
    protocol=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),sources={str(p.relative_to(R)):sha(p) for p in paths},input_manifest_sha256=sha(SOURCE/'raw_manifest.json'),input_campaign='mechanical_work_fixed_20261113',scope='Read-only 64 saved mechanical steps. Native bilinear form, clipped face-pyramid dissipation, nodal weighted work, broken Green boundary integrals; no physics changes.',gates=dict(native_reconstruction=1e-7,commutator_identity=1e-7,broken_green=1e-8,cut_trace=1e-7),seed=11303,timeout_s=600)
    save(O/'protocol.json',protocol);(O/'analyze_at_execution.py').write_bytes((HERE/'analyze.py').read_bytes());save(O/'status.json',dict(state='running'))
    output={};guards=affine_guards()
    for rec in status['cases']:
        if not rec['parameters']['observer']:continue
        for f in rec['files']:assert sha(R/f['path'])==f['sha256'],f['path']
        root=RAW/rec['name'];m=rec['summary'];dt=m['dt_s'];nodes_all=prior.old.read(root/'mechanical_nodes.csv');states=prior.old.read(root/'states.csv');faces_all=prior.old.read(root/'faces.csv');region_all=prior.old.read(root/'regional_mechanics.csv');rows=[]
        for step in sorted(set(map(int,nodes_all['step']))):
            a=states[(states['step']==step)&(states['stage']==0)];mid=states[(states['step']==step)&(states['stage']==1)];nodes=nodes_all[nodes_all['step']==step];faces=faces_all[(faces_all['step']==step)&(faces_all['stage']==0)];native=region_all[(region_all['step']==step)&(region_all['region']==0)][0]
            weight=1-nodes['w0'];wu=(weight[:,None]*np.c_[mid['ux_m_s'],mid['uy_m_s']]).ravel();acc=dict(geometric=0.,commutator=0.,bilinear=0.,cut=0.,wall=0.,internal=0.,D=0.);max_green=0.
            def observe(i,xi,ui,H,u,tri,density,area,dt):
                nonlocal max_green
                Hu=np.einsum('abd,d->ab',H,u);Hw=np.einsum('abd,d->ab',H,wu);s=strain(Hu);sw=strain(Hw);sq=float(np.dot(FACTOR,s*s));bilinear=float(np.dot(FACTOR,s*sw));part=tri.intersection(LEFT);partial=part.area;D=dt*density*partial*sq
                acc['D']+=D;acc['geometric']+=dt*density*(partial-weight[i]*area)*sq;acc['commutator']+=dt*density*area*(weight[i]*sq-bilinear);acc['bilinear']+=dt*density*area*bilinear
                b=boundary_power(part,xi,ui,Hu,density,dt)
                for key in b:acc[key]+=b[key]
                max_green=max(max_green,abs(sum(b.values())-D))
            force,D,DR,traction,length,polys=prior.viscosity(a,mid,faces,nodes,m['viscosity_Pa_s'],dt,observe)
            nodework=dt*np.sum(weight[:,None]*force*np.c_[mid['ux_m_s'],mid['uy_m_s']]);port=acc['geometric']+acc['commutator'];cut=acc['cut'];gap=port-cut;bridge=acc['internal']+acc['wall']-acc['bilinear'];scale=sum(abs(x) for x in [acc['D'],acc['bilinear'],port,cut,acc['internal'],acc['wall']])+1e-40
            errors=dict(native_force=prior.relative(force,np.c_[nodes['gvx_N'],nodes['gvy_N']]),native_dissipation=abs(acc['D']-native['viscous_D_J'])/scale,bilinear_nodework=abs(acc['bilinear']-nodework)/scale,commutator_port=abs(port-native['viscous_port_J'])/scale,broken_green=abs(acc['D']-cut-acc['wall']-acc['internal'])/scale,per_pyramid_green=max_green/scale,cut_trace=abs(cut-native['recovered_viscous_port_J'])/scale,gap_bridge=abs(gap-bridge)/scale)
            passed=max(errors.values())<1e-7 and errors['broken_green']<1e-8
            rows.append(dict(step=step,t_s=step*dt,dt_s=dt,**{k+'_J':v for k,v in acc.items()},port_J=port,gap_J=gap,bridge_J=bridge,errors=errors,passed=passed))
        output[rec['name']]=rows;print(rec['name'],sum(r['passed'] for r in rows),'/',len(rows),flush=True)
    summary={}
    for name,rows in output.items():
        norm=lambda key:float(np.linalg.norm([r[key] for r in rows]));den=max(1e-40,norm('port_J'))
        summary[name]=dict(saved_steps=len(rows),passed=sum(r['passed'] for r in rows),relative_gap=norm('gap_J')/den,geometric_over_port=norm('geometric_J')/den,commutator_over_port=norm('commutator_J')/den,wall_over_port=norm('wall_J')/den,internal_minus_nodal_over_port=float(np.linalg.norm([r['internal_J']-r['bilinear_J'] for r in rows]))/den,maximum_error=max(max(r['errors'].values()) for r in rows))
    save(O/'rows.json',output);save(O/'summary.json',dict(cases=summary,affine_guards=guards,interpretation='Exact discrete port = clipped-vs-cell ownership term + weighted-velocity reconstruction commutator. Exact gap = internal broken-field boundary work + outer-wall work - nodal weighted work. This explains a definition/reconstruction mismatch, not physical convergence. Fields are piecewise affine and generally nonconforming.',full_C1_complete=False))
    passed=all(r['passed'] for rows in output.values() for r in rows);save(O/'status.json',dict(state='completed',passed=passed,saved_steps=sum(map(len,output.values())),wall_s=time.monotonic()-start));assert passed
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
