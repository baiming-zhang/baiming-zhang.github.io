"""Reported-elasticity liner / pure-hoop orthotropic annulus verification.

The pure-hoop material-axis field rotates at every integration point. This
is a plane-strain code benchmark, not a manufactured layup or final COPV.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from skfem import Basis,InteriorFacetBasis,ElementVector,ElementTetP2
from geobc.rod_fem import quadratic_mesh,boundary_facets
from geobc.copv_material import engineering_stiffness,rotate_stiffness,rotate_stress,stress_to_tensor
from geobc.copv_fem import solve_orthotropic,assemble_stiffness,engineering_strain,apply_stiffness
from geobc.copv_layered_reference import PARAMETERS,generate_mesh


def material_input():
    repository=Path(__file__).resolve().parents[1]
    path=repository/'configs/copv_material_candidates.json'
    config=json.loads(path.read_text(encoding='utf-8'))
    candidate=config['elastic_control_candidate']
    if not candidate['usable_as_numerical_elastic_input']:
        raise ValueError('Material candidate not accepted for elastic benchmark')
    materials={1:candidate['liner']['elastic_3d_MPa'],2:candidate['overwrap']['elastic_3d_MPa']}
    sources={}
    for role in ('liner','overwrap'):
        src=config['sources'][candidate[role]['source']]
        pdf=(repository/src['local_file']).resolve()
        actual=hashlib.sha256(pdf.read_bytes()).hexdigest()
        if actual!=src['sha256']:raise ValueError(f'Source PDF checksum mismatch: {pdf}')
        sources[role]={**src,'actual_sha256':actual,'sha256_verified':True}
    provenance=dict(candidate_id=candidate['id'],material_json_path=str(path),material_json_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    units=config['units'],constants_MPa=materials,sources=sources,
                    status='Reported elastic benchmark inputs, not certified material lots; pure-hoop axis assumption is numerical verification only')
    return materials,provenance


def local_stiffness():
    materials,_=material_input()
    return np.array([engineering_stiffness(**materials[k]) for k in (1,2)])


def hoop_axes_from_theta(theta):
    theta=np.asarray(theta,float)
    Q=np.zeros(theta.shape+(3,3))
    Q[...,0,0]=-np.sin(theta);Q[...,1,0]=np.cos(theta)
    Q[...,2,1]=1.
    Q[...,0,2]=np.cos(theta);Q[...,1,2]=np.sin(theta)
    return Q # [e_theta,e_z,e_r], right handed


def hoop_axes(x):
    x=np.asarray(x,float)
    return hoop_axes_from_theta(np.arctan2(x[...,1],x[...,0]))


def stiffness_field(x,material_ids,matrices=None):
    x=np.asarray(x,float);ids=np.broadcast_to(material_ids,x.shape[:-1])
    matrices=local_stiffness() if matrices is None else matrices
    return rotate_stiffness(matrices[ids-1],hoop_axes(x))


def cylindrical_normals(matrices=None):
    """Normal C block in order [r,theta,z], independently permuted from local."""
    matrices=local_stiffness() if matrices is None else matrices
    return matrices[:,[2,0,1]][:,:,[2,0,1]]


def analytic_constants(matrices=None):
    # Normalize radial coordinate rho=r/c to avoid ill-scaled r^(+/-k).
    cn=cylindrical_normals(matrices);k=np.sqrt(cn[:,1,1]/cn[:,0,0])
    a,c,b=PARAMETERS['inner_radius_mm'],PARAMETERS['interface_radius_mm'],PARAMETERS['outer_radius_mm']
    def terms(r,j):
        powers=np.array([(r/c)**k[j],(r/c)**(-k[j])])
        return powers, np.array([cn[j,0,0]*k[j]+cn[j,0,1],-cn[j,0,0]*k[j]+cn[j,0,1]])*powers/r
    _,sa=terms(a,0);_,sb=terms(b,1)
    uc1,sc1=terms(c,0);uc2,sc2=terms(c,1)
    M=np.array([[*sa,0,0],[0,0,*sb],[*sc1,*(-sc2)],[*uc1,*(-uc2)]])
    rhs=np.array([-PARAMETERS['pressure_MPa'],0,0,0])
    scale=np.max(np.abs(M),axis=1)
    coefficients=np.linalg.solve(M/scale[:,None],rhs/scale).reshape(2,2)
    return coefficients,k,cn


def analytic_fields(x,material_ids=None,matrices=None):
    x=np.asarray(x,float);r=np.linalg.norm(x[...,:2],axis=-1)
    ids=np.where(r<=PARAMETERS['interface_radius_mm'],1,2) if material_ids is None else np.broadcast_to(material_ids,r.shape)
    ab,k,cn=analytic_constants(matrices)
    j=ids-1;rho=r/PARAMETERS['interface_radius_mm']
    positive,negative=ab[j,0]*rho**k[j],ab[j,1]*rho**(-k[j])
    ur=positive+negative;dudr=k[j]*(positive-negative)/r
    strain_r,strain_t=dudr,ur/r
    sr=cn[j,0,0]*strain_r+cn[j,0,1]*strain_t
    st=cn[j,1,0]*strain_r+cn[j,1,1]*strain_t
    sz=cn[j,2,0]*strain_r+cn[j,2,1]*strain_t
    er=np.zeros_like(x);er[...,:2]=x[...,:2]/r[...,None]
    et=np.stack([-er[...,1],er[...,0],np.zeros_like(r)],axis=-1)
    sigma=sr[...,None,None]*er[...,:,None]*er[...,None,:]+st[...,None,None]*et[...,:,None]*et[...,None,:]
    sigma[...,2,2]=sz
    return ur[...,None]*er,sigma


def analytic_end_axial_resultants(matrices=None):
    ab,k,cn=analytic_constants(matrices)
    radii=[PARAMETERS['inner_radius_mm'],PARAMETERS['interface_radius_mm'],PARAMETERS['outer_radius_mm']]
    c=PARAMETERS['interface_radius_mm'];result=[]
    for j in range(2):
        ri,ro=radii[j],radii[j+1]
        def radial_integral(power):
            if abs(power+1)<1e-10:return c*np.log(ro/ri)
            return c*((ro/c)**(power+1)-(ri/c)**(power+1))/(power+1)
        value=2*np.pi*((cn[j,2,0]*k[j]+cn[j,2,1])*ab[j,0]*radial_integral(k[j])
                       +(-cn[j,2,0]*k[j]+cn[j,2,1])*ab[j,1]*radial_integral(-k[j]))
        result.append(float(value))
    return result


def solve_reference(directory):
    out=Path(directory);started=time.perf_counter();p=PARAMETERS
    materials,provenance=material_input();matrices=local_stiffness()
    with np.load(out/'mesh.npz') as data:source={k:data[k] for k in data.files}
    mesh,to_original,to_vertex=quadratic_mesh(source['xyz'],source['cells'])
    if not np.array_equal(to_vertex[source['cells'][:,:4]].T,mesh.t):raise ValueError('Changed element order')
    mids=source['material_ids'];basis=Basis(mesh,ElementVector(ElementTetP2()),intorder=4)
    qxyz=np.moveaxis(np.asarray(basis.global_coordinates()),0,-1)
    axes=hoop_axes(qxyz)
    C=stiffness_field(qxyz,mids[:,None],matrices) # TRUE GP direction, not centroid direction
    fixed=boundary_facets(mesh,source['faces'],source['face_group_ids'],to_vertex,1)
    loaded=boundary_facets(mesh,source['faces'],source['face_group_ids'],to_vertex,2)
    interface=boundary_facets(mesh,source['faces'],source['face_group_ids'],to_vertex,4)
    if np.any(mesh.f2t[:,interface]<0) or np.any(mids[mesh.f2t[0,interface]]==mids[mesh.f2t[1,interface]]):raise ValueError('Invalid interface')
    def displacement(x):return analytic_fields(np.asarray(x).T,matrices=matrices)[0].T
    def pressure(x,n):return -p['pressure_MPa']*n
    arrays,diag=solve_orthotropic(mesh,C,fixed,loaded,displacement,pressure)
    exact_u,_=analytic_fields(arrays['xyz'],matrices=matrices)
    _,exact_sigma=analytic_fields(qxyz,mids[:,None],matrices)
    sigma=stress_to_tensor(arrays['stress']);w=np.asarray(arrays['quadrature_weight'])
    u_l2=float(np.linalg.norm(arrays['displacement']-exact_u)/np.linalg.norm(exact_u))
    stress_errors={str(j):float(np.sqrt(np.sum((sigma[mids==j]-exact_sigma[mids==j])**2*w[mids==j,...,None,None])/np.sum(exact_sigma[mids==j]**2*w[mids==j,...,None,None]))) for j in (1,2)}
    coef=np.zeros(basis.N);f=np.zeros(basis.N)
    for i,d in enumerate(basis.split_indices()):coef[d]=arrays['displacement'][:,i];f[d]=arrays['nodal_load'][:,i]
    end0=np.isclose(arrays['xyz'][:,2],0.,atol=1e-9);endL=np.isclose(arrays['xyz'][:,2],p['length_mm'],atol=1e-9)
    residuals={};energies={};reactions={};end_exact=analytic_end_axial_resultants(matrices)
    for layer in (1,2):
        selected=np.flatnonzero(mids==layer)
        Kl=assemble_stiffness(basis.with_elements(selected),C[selected])
        energies[str(layer)]=float(.5*coef@(Kl@coef))
        residual=Kl@coef-(f if layer==1 else 0)
        nodal=np.column_stack([residual[d] for d in basis.split_indices()]);residuals[layer]=nodal
        bottom,top=nodal[end0].sum(0),nodal[endL].sum(0)
        expected=end_exact[layer-1]
        error=float(max(abs(bottom[2]/(-expected)-1),abs(top[2]/expected-1)))
        reactions[str(layer)]=dict(z0_N=bottom.tolist(),zL_N=top.tolist(),analytic_zL_axial_N=expected,axial_relative_error=error)
    ifnodes=np.isclose(np.linalg.norm(arrays['xyz'][:,:2],axis=1),p['interface_radius_mm'],atol=1e-7)&~(end0|endL)
    weak=float(np.linalg.norm((residuals[1]+residuals[2])[ifnodes])/np.linalg.norm(residuals[1][ifnodes]))
    traces=[];values=[];trace_mids=[];coords=[]
    for side in (0,1):
        fb=InteriorFacetBasis(mesh,basis.elem,facets=interface,intorder=5,side=side)
        x=np.moveaxis(np.asarray(fb.global_coordinates()),0,-1)
        field=fb.interpolate(coef);cell_mids=mids[fb.tind]
        face_C=stiffness_field(x,cell_mids[:,None],matrices)
        traces.append(stress_to_tensor(np.moveaxis(apply_stiffness(face_C,engineering_strain(field)),0,-1)))
        values.append(np.moveaxis(np.asarray(field),0,-1));trace_mids.append(cell_mids);coords.append(x)
    if not np.allclose(coords[0],coords[1],atol=1e-10):raise ValueError('Interface point mismatch')
    x=coords[0];er=x.copy();er[...,2]=0;er/=np.linalg.norm(er,axis=-1)[...,None]
    et=np.stack([-er[...,1],er[...,0],np.zeros(er.shape[:2])],axis=-1)
    s1=np.where((trace_mids[0]==1)[:,None,None,None],traces[0],traces[1])
    s2=np.where((trace_mids[0]==2)[:,None,None,None],traces[0],traces[1])
    t1=np.einsum('fqij,fqj->fqi',s1,er);t2=np.einsum('fqij,fqj->fqi',s2,er)
    qw=np.asarray(fb.dx)
    ifpoint=np.array([[p['interface_radius_mm'],0,0]])
    sif1=analytic_fields(ifpoint,np.array([1]),matrices)[1][0]
    sif2=analytic_fields(ifpoint,np.array([2]),matrices)[1][0]
    pressure_if=-float(sif1[0,0])
    def mean(value):return float(np.sum(value*qw)/qw.sum())
    rr1,rr2=np.sum(t1*er,axis=-1),np.sum(t2*er,axis=-1)
    hoop1=np.einsum('fqi,fqij,fqj->fq',et,s1,et);hoop2=np.einsum('fqi,fqij,fqj->fq',et,s2,et)
    interface_result=dict(displacement_jump_max_mm=float(np.linalg.norm(values[0]-values[1],axis=-1).max()),
        weak_interface_force_relative_residual=weak,
        traction_trace_jump_relative_surface_L2=float(np.sqrt(np.sum(np.sum((t1-t2)**2,axis=-1)*qw)/qw.sum())/pressure_if),
        radial_stress_layer1_mean_MPa=mean(rr1),radial_stress_layer2_mean_MPa=mean(rr2),analytic_radial_stress_MPa=-pressure_if,
        radial_trace_mean_max_relative_error=max(abs(mean(rr1)/(-pressure_if)-1),abs(mean(rr2)/(-pressure_if)-1)),
        hoop_stress_layer1_mean_MPa=mean(hoop1),hoop_stress_layer2_mean_MPa=mean(hoop2),
        analytic_hoop_layer1_MPa=float(sif1[1,1]),analytic_hoop_layer2_MPa=float(sif2[1,1]))
    local_sigma=rotate_stress(arrays['stress'],axes,to_global=False)
    # [fiber=theta, transverse=z, normal=r]; independent reference has no shear.
    exact_local=np.einsum('...ki,...kl,...lj->...ij',axes,exact_sigma,axes)
    actual_local=stress_to_tensor(local_sigma)
    local_error=float(np.sqrt(np.sum((actual_local-exact_local)**2*w[...,None,None])/np.sum(exact_local**2*w[...,None,None])))
    inner_ur=float(analytic_fields(np.array([[p['inner_radius_mm'],0,0]]),matrices=matrices)[0][0,0])
    exact_energy=float(.5*p['pressure_MPa']*2*np.pi*p['inner_radius_mm']*p['length_mm']*inner_ur)
    work=diag['traction_work_Nmm']+diag['prescribed_boundary_reaction_work_Nmm']
    row=dict(parameters=p,material_provenance=provenance,n_nodes=diag['n_nodes'],n_elements=diag['n_elements'],
        seconds=time.perf_counter()-started,displacement_relative_L2=u_l2,stress_tensor_relative_volume_L2_by_layer=stress_errors,
        rotate_back_local_stress_relative_volume_L2=local_error,interface=interface_result,
        end_reactions_by_layer=reactions,max_layer_axial_reaction_relative_error=max(r['axial_relative_error'] for r in reactions.values()),
        strain_energy_by_layer_Nmm=energies,strain_energy_analytic_Nmm=exact_energy,
        energy_analytic_relative_error=abs(diag['strain_energy_Nmm']/exact_energy-1),
        boundary_work_energy_relative_error=abs(work/(2*diag['strain_energy_Nmm'])-1),
        axes_max_within_element_component_range=float(np.ptp(axes[mids==2],axis=1).max()),
        stiffness_max_within_element_relative_range=float(np.max(np.ptp(C[mids==2],axis=1))/np.max(np.abs(C[mids==2]))),
        stiffness_association='actual C[E,Q,6,6] evaluated and rotated at every volume GP; face traces likewise rotate at every interface GP',
        diagnostics=diag)
    U=np.empty_like(arrays['displacement']);U[to_original]=arrays['displacement']
    np.savez_compressed(out/'solution.npz',**source,displacement=U,stress=arrays['stress'],stress_local=local_sigma,
                        material_axes=axes,quadrature_xyz=qxyz,quadrature_weight=w)
    (out/'solution.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
    return row
