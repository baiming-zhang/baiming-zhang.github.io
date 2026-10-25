"""Three-component axisymmetric vessel verification, not production labels.

Only geometry/protocol is read from the 3D source to construct the 2D solve.
Displacements/stresses from 3D are loaded by the comparison tool afterwards.
Quad9 physical geometry uses all normal-offset nodes and exact reference edges.
"""
from __future__ import annotations
from copy import copy
import json
from pathlib import Path
import time
import numpy as np
from skfem import Basis,FacetBasis,ElementLineP1,ElementQuad2,LinearForm,asm
from skfem.quadrature import get_quadrature
from geobc.copv_axisymmetric import quad9_meridian,vector_basis,material_axes_rtz,assemble_axisymmetric,solve_axial_mpc,strain_axisymmetric,apply_C
from geobc.copv_structured_mesh import analytic_segments
from geobc.copv_orthotropic_reference import local_stiffness,material_input
from geobc.copv_material import rotate_stiffness,rotate_stress,rotate_strain

ANGLES=np.array([55.,-55.,-55.,55.])


def meridian_mesh(parameters,counts,liner_elements,wrap_per_ply,source_mesh_path=None):
    ns=sum(counts);nw=4*wrap_per_ply;nt=liner_elements+nw
    a=quad9_meridian(nt,ns,1,2,1)
    rz,normal=analytic_segments(parameters,counts)
    offsets=np.r_[np.linspace(0,parameters['liner_thickness_mm'],2*liner_elements+1),
                   parameters['liner_thickness_mm']+np.linspace(0,parameters['overwrap_thickness_mm'],2*nw+1)[1:]]
    ij=a['node_parametric_indices'];k,s=ij.T
    xyz=rz[s]+normal[s]*offsets[k,None]
    if source_mesh_path is not None:
        # No solution/label arrays are accessed here.
        with np.load(source_mesh_path) as source:
            si=source['structured_node_ijk'];selected=np.flatnonzero(si[:,0]==0)
            lookup={(int(si[i,2]),int(si[i,1])):i for i in selected}
            rows=np.array([lookup[tuple(v)] for v in ij]);original=source['xyz'][rows]
            actual=np.c_[np.linalg.norm(original[:,:2],axis=1),original[:,2]]
            if np.max(np.abs(actual-xyz))>1e-10:raise ValueError('Actual phi=0 geometry differs from analytic parameter grid')
            xyz=actual
    a['xyz']=xyz;a['mesh'].doflocs=xyz[a['sk_to_original']].T.copy()
    if hasattr(a['mesh'],'_cached_mapping'):delattr(a['mesh'],'_cached_mapping')
    a.update(parameters_mm=parameters,meridian_segment_elements=list(counts),liner_elements=liner_elements,
             wrap_per_ply=wrap_per_ply,offsets_mm=offsets,meridian_inner_rz=rz,meridian_normal_rz=normal)
    return a


def exact_edge_basis(mesh,element,facets,intorder=5):
    find=np.asarray(facets,int);tind=mesh.f2t[0,find]
    if np.any(mesh.f2t[1,find]>=0):raise ValueError('Exterior meridian edges expected')
    X,W=get_quadrature(mesh.brefdom,intorder);ref=ElementQuad2().doflocs
    corner=mesh.facets[:,find];match=corner[:,None,:]==mesh.t[None,:,tind]
    if not np.all(match.sum(1)==1):raise ValueError('Ambiguous reference edge corners')
    local=match.argmax(1);rc=ref[local].transpose(2,1,0);line=ElementLineP1()
    phi=np.stack([line.lbasis(X,j)[0] for j in range(2)]);Y=np.einsum('afj,jq->afq',rc,phi)
    original=mesh.mapping();mapping=copy(original);physical=original.G(X,find=find)
    if np.max(abs(original.F(Y,tind=tind)-physical))>1e-10:raise ValueError('Reference edge geometry mismatch')
    if original.detDF(Y,tind=tind).min()<=0:raise ValueError('Nonpositive parent edge Jacobian')
    def inverse(x,tind=None,**kwargs):
        if not np.array_equal(tind,parents) or np.shape(x)!=physical.shape or np.max(abs(np.asarray(x)-physical))>1e-10:raise ValueError('Wrong exact-edge request')
        return Y
    parents=tind.copy();mapping.invF=inverse
    fb=FacetBasis(mesh,element,mapping=mapping,facets=find,quadrature=(X,W));fb.reference_volume_coordinates=Y
    return fb


def edge_load(basis,facets,traction,intorder=5):
    fb=exact_edge_basis(basis.mesh,basis.elem,facets,intorder)
    @LinearForm
    def loading(v,w):
        n=np.array([w.n[0],np.zeros_like(w.n[0]),w.n[1]])
        return 2*np.pi*w.x[0]*np.sum(traction(w.x,n)*v,axis=0)
    return asm(loading,fb),fb


def frames(basis,alpha,reference=None,elements=None):
    reference=basis.X if reference is None else reference
    DF=basis.mapping.DF(reference,tind=elements)
    ds=np.moveaxis(DF[:,1],0,-1)
    if np.any(ds[...,1]<=0):raise ValueError('Meridian is not increasing z')
    Q=material_axes_rtz(ds,alpha)
    t=ds/np.linalg.norm(ds,axis=-1)[...,None]
    n=np.stack([t[...,1],np.zeros(t.shape[:-1]),-t[...,0]],axis=-1)
    dt=np.moveaxis(DF[:,0],0,-1)
    cos=np.sum(dt*n[...,[0,2]],axis=-1)/np.linalg.norm(dt,axis=-1)
    if cos.min()<=0:raise ValueError('Normal-offset orientation mismatch')
    return Q,dict(min_offset_normal_cosine=float(cos.min()),min_material_frame_det=float(np.linalg.det(Q).min()))


def solve_vessel_2d(output,parameters,counts,liner_elements=2,wrap_per_ply=1,case='reported',pressure=1.,source_mesh_path=None):
    out=Path(output);out.mkdir(parents=True,exist_ok=True);started=time.perf_counter();cpu_started=time.process_time()
    a=meridian_mesh(parameters,counts,liner_elements,wrap_per_ply,source_mesh_path)
    mesh_done=time.perf_counter()
    b=vector_basis(a,intorder=5);m=a['mesh'];indices=a['cell_parametric_indices'];offset_cell=indices[:,0]
    min_jac=float(np.min(b.mapping.detDF(b.X)))
    if not np.isfinite(min_jac) or min_jac<=0:
        raise ValueError(f'Nonpositive vessel meridian volume-GP Jacobian: {min_jac}')
    ply=np.where(offset_cell<liner_elements,0,1+(offset_cell-liner_elements)//wrap_per_ply).astype(int)
    alpha=np.zeros(len(ply));alpha[ply>0]=ANGLES[ply[ply>0]-1]
    Q,frame_info=frames(b,np.deg2rad(alpha[:,None]));matrices=local_stiffness()
    if case=='isotropic':matrices[1]=matrices[0]
    elif case!='reported':raise ValueError('Unsupported material case')
    C=rotate_stiffness(matrices[(ply>0).astype(int)][:,None],Q);K=assemble_axisymmetric(b,C)
    Fpressure,pb=edge_load(b,a['facets']['inner'],lambda x,n:-pressure*n)
    cap=pressure*np.pi*parameters['boss_bore_radius_mm']**2;Fcap=np.zeros(b.N);capareas={}
    for name,sign in [('bottom',-1),('top',1)]:
        fb=exact_edge_basis(m,b.elem,a['facets'][name]);x=np.asarray(fb.global_coordinates())
        area=float(np.sum(2*np.pi*x[0]*fb.dx));capareas[name]=area
        def traction(x,n):
            v=np.zeros_like(n);v[2]=sign*cap/area;return v
        Fcap+=edge_load(b,a['facets'][name],traction)[0]
    loads_done=time.perf_counter()
    data,diag=solve_axial_mpc(a,b,K,Fpressure+Fcap,C)
    solve_done=time.perf_counter()
    scalar=data['displacement_rtz'][a['sk_to_original']];coef=np.zeros(b.N)
    for i,d in enumerate(b.split_indices()):coef[d]=scalar[:,i]
    pwork=float(coef@Fpressure);capwork=float(coef@Fcap)
    px=np.asarray(pb.global_coordinates());pn=np.asarray(pb.normals);pw=np.asarray(pb.dx)*2*np.pi*px[0]
    half=float(np.sum((-pressure*pn[1]*pw)[px[1]>0])+cap);expected_half=pressure*np.pi*parameters['inner_radius_mm']**2
    # Explicit local Quad9 reference coordinates, not Gmsh numbering assumptions.
    sk_top=np.flatnonzero(ElementQuad2().doflocs[:,1]==1)
    original_top=a['sk_to_original'][m.dofs.element_dofs[sk_top]].T
    cut_cells=np.flatnonzero(np.all(np.isclose(a['xyz'][original_top,1],0,atol=1e-8),axis=1))
    cut=None
    if len(cut_cells):
        g,w=np.polynomial.legendre.leggauss(3);g=(g+1)/2;w=w/2;ref=np.stack([g,np.ones(3)])
        cb=Basis(m,b.elem,quadrature=(ref,w),elements=cut_cells)
        qc=np.asarray(cb.global_coordinates());cq,_=frames(cb,np.deg2rad(alpha[cut_cells,None]),reference=ref,elements=cut_cells)
        cc=rotate_stiffness(matrices[(ply[cut_cells]>0).astype(int)][:,None],cq)
        sigma=apply_C(cc,strain_axisymmetric(cb.interpolate(coef),qc[0]))
        df=cb.mapping.DF(ref,tind=cut_cells);line_jac=np.linalg.norm(df[:,0],axis=0)
        areaweight=2*np.pi*qc[0]*line_jac*w[None,:]
        resultant=float(np.sum(sigma[2]*areaweight));cut=dict(axial_force_N=resultant,expected_N=expected_half,relative_error=abs(resultant/expected_half-1))
    _,provenance=material_input()
    qcount=data['stress_rtz_voigt'].shape[1]
    stats_data=dict(stress_local=rotate_stress(data['stress_rtz_voigt'],Q,to_global=False),
        strain_local=rotate_strain(data['strain_rtz_engineering'],Q,to_global=False),
        quadrature_weight=data['quadrature_volume_weight'],gp_ply_ids=np.broadcast_to(ply[:,None],(len(ply),qcount)))
    from tools.analyze_copv_convergence import layer_statistics
    stats=layer_statistics(stats_data)
    diag.update(pressure_work_Nmm=pwork,cap_work_Nmm=capwork,pressure_area_mm2=float(pw.sum()),
        cap_resultant_N=cap,cap_area_mm2=capareas,half_axial_load_N=half,expected_half_axial_load_N=expected_half,
        half_body_axial_load_relative_error=abs(half/expected_half-1),midsection_stress_integral=cut,
        boundary_work_energy_relative_error=abs((pwork+capwork)/(2*diag['strain_energy_Nmm'])-1),
        cap_work_vs_force_extension_relative_error=abs(capwork/(cap*diag['end_extension_mm'])-1))
    assert max(diag['reduced_relative_residual'],diag['boundary_work_energy_relative_error'],diag['cap_work_vs_force_extension_relative_error'])<1e-8
    assert diag['half_body_axial_load_relative_error']<.01
    manifest=dict(parameters_mm=parameters,meridian_segment_elements=list(counts),liner_elements=liner_elements,wrap_per_ply=wrap_per_ply)
    report=dict(case=case,status='independent_axisymmetric_vessel_validation_not_production',production_ready=False,
        manifest=manifest,n_nodes=len(a['xyz']),n_quad9=len(ply),seconds=time.perf_counter()-started,pressure_MPa=pressure,
        wrap_ply_angles_degrees=ANGLES.tolist(),material_provenance=provenance,frame_checks=frame_info,
        max_abs_utheta_mm=float(abs(data['displacement_rtz'][:,1]).max()),max_abs_ur_mm=float(abs(data['displacement_rtz'][:,0]).max()),
        min_meridian_quadrature_Jacobian=min_jac,diagnostics=diag,physical_ply_statistics=stats,
        limitations=['Exact revolution of P2 meridian differs slightly from 3D P2 circumference approximation',
                     'Axisymmetric three-component prototype; no nonaxisymmetric BC or defect representation',
                     'Reported numerical plies and sliding end resultants, not manufactured COPV certification'])
    metrics_done=time.perf_counter()
    np.savez_compressed(out/'solution.npz',**data,**stats_data,material_axes_rtz=Q,cell_ply_ids=ply,
        cell_parametric_indices=indices,node_parametric_indices=a['node_parametric_indices'],offsets_mm=a['offsets_mm'],
        original_quad9_cells=a['cells'])
    npz_done=time.perf_counter()
    (out/'solution.json').write_text(json.dumps(report,indent=2)+'\n')
    json_done=time.perf_counter()
    timing=dict(mesh_wall_seconds=mesh_done-started,assembly_and_load_wall_seconds=loads_done-mesh_done,
        solve_and_field_extraction_wall_seconds=solve_done-loads_done,metrics_wall_seconds=metrics_done-solve_done,
        npz_compression_and_write_wall_seconds=npz_done-metrics_done,json_export_wall_seconds=json_done-npz_done,
        pipeline_wall_seconds=json_done-started,process_CPU_seconds=time.process_time()-cpu_started,
        scope='Whole solve_vessel_2d call through NPZ+JSON export; tiny timing-sidecar write excluded. Does not include later 3D comparison/plotting.')
    (out/'timing.json').write_text(json.dumps(timing,indent=2)+'\n');report['timing']=timing
    return report
