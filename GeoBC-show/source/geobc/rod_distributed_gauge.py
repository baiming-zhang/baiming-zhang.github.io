"""Area-integral tangential/axial bore restraints, a separately approved BC model."""
import numpy as np
from skfem import LinearForm,asm,ElementVector,ElementTetP2
from geobc.rod_fem import quadratic_mesh,boundary_facets
from geobc.tet_reference_faces import reference_facet_basis

def coefficients(m,bc,L):
    mesh,to_original,to_vertex=quadratic_mesh(m['xyz'],m['cells'])
    facets=boundary_facets(mesh,m['faces'],m['face_group_ids'],to_vertex,1)
    fb=reference_facet_basis(mesh,ElementVector(ElementTetP2()),facets,intorder=8)
    @LinearForm
    def tangent(v,w):
        x=w.x[0]-L;y=w.x[1];r=np.sqrt(x*x+y*y)
        return -y/r*v[0]+x/r*v[1]
    @LinearForm
    def axial(v,w):return v[2]
    area=float(fb.dx.sum());split=fb.split_indices()
    result=[]
    for form in [tangent,axial]:
        flat=asm(form,fb);ordered=np.column_stack([flat[d] for d in split]);native=np.empty_like(ordered);native[to_original]=ordered;result.append(native/area)
    ct,cz=result;normal=bc['radial_normal'].astype(float);members=np.linalg.norm(normal,axis=1)>.5
    assert np.max(np.abs(cz.sum(0)-[0,0,1]))<1e-10
    assert np.max(np.abs(ct[~members]))<1e-12 and np.max(np.abs(cz[~members]))<1e-12
    assert np.linalg.norm(ct.sum(0))<1e-4
    # Choose radial pivots of largest magnitude; integral equations use only free DOFs.
    radial=[];spc=[];reduced=ct.copy();dependent=set()
    for i in np.flatnonzero(members):
        n=normal[i,:2];p=int(np.argmax(np.abs(n)));q=1-p
        dependent.add((int(i),p));reduced[i,q]-=ct[i,p]*n[q]/n[p];reduced[i,p]=0
        if abs(n[q])<1e-9:spc.append((int(i),p))
        else:radial.append([(int(i),p,float(n[p])),(int(i),q,float(n[q]))])
    equations=[]
    for row in [reduced,cz]:
        terms=[(int(i),int(j),float(row[i,j])) for i,j in np.argwhere(np.abs(row)>1e-15)]
        pivot=max(range(len(terms)),key=lambda k:abs(terms[k][2]));terms.insert(0,terms.pop(pivot));assert terms[0][:2] not in dependent
        dependent.add(terms[0][:2]);assert len({(i,j) for i,j,c in terms})==len(terms)
        scale=abs(terms[0][2]);equations.append([(i,j,c/scale) for i,j,c in terms])
    # Test the algebraic radial substitution with a deterministic admissible field.
    rng=np.random.default_rng(2026);u=rng.normal(size=m['xyz'].shape)
    for i in np.flatnonzero(members):
        n=normal[i,:2];p=int(np.argmax(np.abs(n)));q=1-p;u[i,p]=-n[q]/n[p]*u[i,q]
    assert abs(np.sum(ct*u)-np.sum(reduced*u))<1e-12
    return dict(tangent=ct,axial=cz,radial=radial,spc=spc,integrals=equations,area_mm2=area,members=members)

def replace_boundary(deck,m,c):
    text=deck.read_text();a=text.index('*BOUNDARY\n');b=text.index('*STEP,',a);ids=m['node_ids']
    lines=['*BOUNDARY']+[f'{int(ids[i])},{j+1},{j+1},0' for i,j in c['spc']]
    for terms in c['radial']+c['integrals']:
        lines+=['*EQUATION',str(len(terms))]
        for start in range(0,len(terms),4):lines.append(','.join(f'{int(ids[i])},{j+1},{v:.14g}' for i,j,v in terms[start:start+4]))
    deck.write_text(text[:a]+'\n'.join(lines)+'\n'+text[b:],encoding='ascii')

def equilibrium(folder,m,bc,c,parameters):
    from tools.check_rod_v3_pilot import diagnostics
    from geobc.rod_v3_ccx import read_dat
    from geobc.rod_io import write_json
    report=diagnostics(folder,m,bc,parameters=parameters)
    f=np.load(folder/'fields.npz');force=f['force_peak_N'];blocks=read_dat(folder/'rod.dat')
    rf=next(x['rows'] for x in blocks if x['kind']=='RF' and abs(x['time']-1)<1e-7)
    mapping={int(n):i for i,n in enumerate(m['node_ids'])};r=np.zeros_like(force)
    for row in rf:r[mapping[int(row[0])]]=row[1:]
    r-=force;n=bc['radial_normal'].astype(float);n/=np.maximum(np.linalg.norm(n,axis=1,keepdims=True),1e-30)
    def remove_radial(v):return v-n*np.sum(n*v,axis=1)[:,None]
    residual=remove_radial(r);basis=np.column_stack([remove_radial(c['tangent']).ravel(),remove_radial(c['axial']).ravel()])
    multipliers=np.linalg.lstsq(basis,residual.ravel(),rcond=None)[0]
    free=float(np.linalg.norm(residual.ravel()-basis@multipliers)/max(np.linalg.norm(force),1))
    t=float(np.sum(c['tangent']*f['U_peak_mm']));z=float(np.sum(c['axial']*f['U_peak_mm']))
    report.update(original_point_gauge_free_residual_not_applicable=report['free_residual_relative'],free_residual_relative=free,mean_tangential_U_mm=t,mean_axial_U_mm=z,boundary_model='area_integral_reference_bore_tangent_and_axial',passed=report['force_balance_relative']<1e-5 and report['moment_balance_deformed_relative']<1e-5 and free<1e-4 and report['max_reference_radial_constraint_violation_mm']<1e-5 and max(abs(t),abs(z))<1e-5)
    write_json(folder/'equilibrium_review.json',report)
    return report
