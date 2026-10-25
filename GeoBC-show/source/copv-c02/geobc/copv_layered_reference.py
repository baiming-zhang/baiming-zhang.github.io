"""Bonded two-isotropic-layer 3-D annulus verification, synthetic materials.

Plane strain with analytic end displacements; no anisotropic overwrap,
closed-end COPV, production material, damage, or burst claim.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json
import time
import numpy as np
from skfem import Basis,FacetBasis,InteriorFacetBasis,ElementVector,ElementTetP2
from geobc.rod_fem import quadratic_mesh,boundary_facets
from geobc.copv_material import engineering_stiffness,stress_to_tensor
from geobc.copv_fem import solve_orthotropic,assemble_stiffness,engineering_strain,apply_stiffness

PARAMETERS=dict(inner_radius_mm=20.,interface_radius_mm=25.,outer_radius_mm=30.,length_mm=30.,pressure_MPa=1.)
MATERIALS={1:dict(E_MPa=70000.,nu=.33),2:dict(E_MPa=135000.,nu=.28)}


def analytic_constants(parameters=None,materials=None):
    p=dict(PARAMETERS);p.update(parameters or {})
    m=materials or MATERIALS
    a,c,b=p['inner_radius_mm'],p['interface_radius_mm'],p['outer_radius_mm']
    f1=(1+m[1]['nu'])/m[1]['E_MPa'];f2=(1+m[2]['nu'])/m[2]['E_MPa']
    matrix=np.array([[1,-1/a**2,0,0],[0,0,1,-1/b**2],
                     [1,-1/c**2,-1,1/c**2],
                     [f1*(1-2*m[1]['nu'])*c,f1/c,-f2*(1-2*m[2]['nu'])*c,-f2/c]])
    result=np.linalg.solve(matrix,[-p['pressure_MPa'],0,0,0])
    return result.reshape(2,2) # [[A1,B1],[A2,B2]]


def analytic_fields(x,material_ids=None,parameters=None,materials=None):
    p=dict(PARAMETERS);p.update(parameters or {})
    m=materials or MATERIALS
    constants=analytic_constants(p,m)
    x=np.asarray(x,float);r=np.linalg.norm(x[...,:2],axis=-1)
    ids=np.where(r<=p['interface_radius_mm'],1,2) if material_ids is None else np.broadcast_to(material_ids,r.shape)
    A,B=constants[ids-1,0],constants[ids-1,1]
    E=np.where(ids==1,m[1]['E_MPa'],m[2]['E_MPa'])
    nu=np.where(ids==1,m[1]['nu'],m[2]['nu'])
    er=np.zeros_like(x);er[...,:2]=x[...,:2]/r[...,None]
    et=np.stack([-er[...,1],er[...,0],np.zeros_like(r)],axis=-1)
    ur=(1+nu)/E*((1-2*nu)*A*r+B/r)
    sr,st=A-B/r**2,A+B/r**2
    sigma=sr[...,None,None]*er[...,:,None]*er[...,None,:]+st[...,None,None]*et[...,:,None]*et[...,None,:]
    sigma[...,2,2]=2*nu*A
    return ur[...,None]*er,sigma


def generate_mesh(output_dir,mesh_size_mm,parameters=None):
    import gmsh
    p=dict(PARAMETERS);p.update(parameters or {})
    a,c,b,L=(p[k] for k in ('inner_radius_mm','interface_radius_mm','outer_radius_mm','length_mm'))
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        gmsh.option.setNumber('Mesh.RandomSeed',2026)
        gmsh.model.add('bonded_two_layer_plane_strain_reference')
        occ=gmsh.model.occ
        def annulus(ri,ro):
            outer,inner=occ.addCircle(0,0,0,ro),occ.addCircle(0,0,0,ri)
            f=occ.addPlaneSurface([occ.addCurveLoop([outer]),occ.addCurveLoop([inner])])
            return [dt for dt in occ.extrude([(2,f)],0,0,L) if dt[0]==3]
        v1,v2=annulus(a,c),annulus(c,b)
        pieces,_=occ.fragment(v1,v2)
        occ.synchronize()
        volumes=[dt for dt in pieces if dt[0]==3]
        if len(volumes)!=2:
            raise ValueError(f'Expected two material volumes, got {volumes}')
        volume_layer={}
        for _,tag in volumes:
            if gmsh.model.isInside(3,tag,[(a+c)/2,0,L/2]):volume_layer[tag]=1
            elif gmsh.model.isInside(3,tag,[(b+c)/2,0,L/2]):volume_layer[tag]=2
            else:raise ValueError('Cannot assign material volume')
        if sorted(volume_layer.values())!=[1,2]:raise ValueError('Material partition incorrect')
        surface_counts=Counter(t for d,t in gmsh.model.getBoundary(volumes,combined=False,oriented=False) if d==2)
        groups={1:[],2:[],3:[],4:[]}
        for tag,count in surface_counts.items():
            kind=gmsh.model.getType(2,tag);area=occ.getMass(2,tag)
            if count==2:
                if kind!='Cylinder' or abs(area-2*np.pi*c*L)>1e-6:raise ValueError('Unexpected shared interface')
                groups[4].append(tag)
            elif kind=='Plane':groups[1].append(tag)
            elif abs(area-2*np.pi*a*L)<1e-6:groups[2].append(tag)
            elif abs(area-2*np.pi*b*L)<1e-6:groups[3].append(tag)
            else:raise ValueError(f'Unidentified surface {tag}')
        if [len(groups[k]) for k in (1,2,3,4)]!=[4,1,1,1]:raise ValueError('Shared material interface not conformal')
        for k,tags in groups.items():gmsh.model.addPhysicalGroup(2,tags,k)
        for tag,layer in volume_layer.items():gmsh.model.addPhysicalGroup(3,[tag],10+layer)
        gmsh.option.setNumber('Mesh.MeshSizeMin',mesh_size_mm)
        gmsh.option.setNumber('Mesh.MeshSizeMax',mesh_size_mm)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0)
        gmsh.option.setNumber('Mesh.Algorithm3D',1)
        gmsh.option.setNumber('Mesh.SecondOrderLinear',0)
        gmsh.model.mesh.generate(3);gmsh.model.mesh.optimize('Netgen')
        gmsh.model.mesh.setOrder(2);gmsh.model.mesh.optimize('HighOrder')
        node_ids,coords,_=gmsh.model.mesh.getNodes()
        node_ids,xyz=np.asarray(node_ids,np.int64),np.asarray(coords).reshape(-1,3)
        order=np.argsort(node_ids);node_ids,xyz=node_ids[order],xyz[order]
        all_cells,all_cell_ids,material_ids=[],[],[]
        for volume,layer in sorted(volume_layer.items(),key=lambda kv:kv[1]):
            ids,nodes=gmsh.model.mesh.getElementsByType(11,volume)
            cells=np.asarray(nodes,np.int64).reshape(-1,10)
            all_cells.append(cells);all_cell_ids.extend(ids);material_ids.extend([layer]*len(cells))
        conn=np.concatenate(all_cells);cell_ids=np.array(all_cell_ids,np.int64)
        active=np.isin(node_ids,np.unique(conn));node_ids,xyz=node_ids[active],xyz[active]
        cells=np.searchsorted(node_ids,conn)
        faces,face_groups=[],[]
        for k,tags in groups.items():
            for tag in tags:
                _,nodes=gmsh.model.mesh.getElementsByType(9,tag)
                nodes=np.asarray(nodes,np.int64).reshape(-1,6)
                assert np.isin(nodes,node_ids).all()
                faces.append(np.searchsorted(node_ids,nodes));face_groups.extend([k]*len(nodes))
        quality=np.asarray(gmsh.model.mesh.getElementQualities(cell_ids,'minSICN'))
        if quality.min()<=0:raise ValueError('Nonpositive quadratic element quality')
        faces,face_groups=np.concatenate(faces),np.array(face_groups,np.int32)
        used1=np.unique(cells[np.array(material_ids)==1]);used2=np.unique(cells[np.array(material_ids)==2])
        shared=np.intersect1d(used1,used2);interface=np.unique(faces[face_groups==4])
        if not np.array_equal(shared,interface):raise ValueError('Shared volume nodes do not match interface nodes')
        np.savez_compressed(out/'mesh.npz',xyz=xyz,cells=cells,cell_ids=cell_ids,material_ids=np.array(material_ids,np.int32),
                            node_ids=node_ids,faces=faces,face_group_ids=face_groups)
        gmsh.write(str(out/'geometry.step'))
        info=dict(parameters=p,mesh_size_mm=mesh_size_mm,n_nodes=len(xyz),n_elements=len(cells),
                  min_sicn=float(quality.min()),n_shared_quadratic_interface_nodes=len(shared),
                  conformal_interface=True,face_groups={1:'end planes',2:'inner pressure',3:'outer free',4:'bonded material interface'},
                  material_volumes=volume_layer,gmsh_version=gmsh.__version__)
        (out/'mesh_manifest.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
        return info
    finally:gmsh.finalize()


def material_stiffness():
    matrices=[]
    for k in (1,2):
        E,nu=MATERIALS[k]['E_MPa'],MATERIALS[k]['nu'];G=E/(2*(1+nu))
        matrices.append(engineering_stiffness(E1=E,E2=E,E3=E,nu12=nu,nu13=nu,nu23=nu,G12=G,G13=G,G23=G))
    return np.array(matrices)


def solve_reference(directory):
    out=Path(directory);started=time.perf_counter();p=PARAMETERS
    with np.load(out/'mesh.npz') as data:source={k:data[k] for k in data.files}
    mesh,to_original,to_vertex=quadratic_mesh(source['xyz'],source['cells'])
    if not np.array_equal(to_vertex[source['cells'][:,:4]].T,mesh.t):
        raise ValueError('Element ordering changed during mesh conversion')
    mids=source['material_ids'];C=material_stiffness()[mids-1]
    fixed=boundary_facets(mesh,source['faces'],source['face_group_ids'],to_vertex,1)
    loaded=boundary_facets(mesh,source['faces'],source['face_group_ids'],to_vertex,2)
    interface=boundary_facets(mesh,source['faces'],source['face_group_ids'],to_vertex,4)
    if np.any(mesh.f2t[:,interface]<0) or np.any(mids[mesh.f2t[0,interface]]==mids[mesh.f2t[1,interface]]):
        raise ValueError('Interface facets do not separate two materials')
    def displacement(x):return analytic_fields(np.asarray(x).T)[0].T
    def pressure(x,n):return -p['pressure_MPa']*n
    arrays,diag=solve_orthotropic(mesh,C,fixed,loaded,displacement,pressure)
    exact_u,_=analytic_fields(arrays['xyz'])
    _,exact_sigma=analytic_fields(arrays['quadrature_xyz'],mids[:,None])
    sigma=stress_to_tensor(arrays['stress']);w=arrays['quadrature_weight']
    u_l2=float(np.linalg.norm(arrays['displacement']-exact_u)/np.linalg.norm(exact_u))
    def stress_error(mask):
        return float(np.sqrt(np.sum((sigma[mask]-exact_sigma[mask])**2*w[mask,...,None,None])/np.sum(exact_sigma[mask]**2*w[mask,...,None,None])))
    stress_l2={str(k):stress_error(mids==k) for k in (1,2)}
    basis=Basis(mesh,ElementVector(ElementTetP2()),intorder=4)
    coef=np.zeros(basis.N)
    for k,dofs in enumerate(basis.split_indices()):coef[dofs]=arrays['displacement'][:,k]
    # Separate element assembly partitions reactions without averaging interface stresses.
    layer_residual={};layer_energy={};end_reactions={}
    end0=np.isclose(arrays['xyz'][:,2],0.,atol=1e-9);endL=np.isclose(arrays['xyz'][:,2],p['length_mm'],atol=1e-9)
    radii=[p['inner_radius_mm'],p['interface_radius_mm'],p['outer_radius_mm']]
    constants=analytic_constants();end_errors=[]
    nodal_load=np.zeros(basis.N)
    for k,dofs in enumerate(basis.split_indices()):nodal_load[dofs]=arrays['nodal_load'][:,k]
    for layer in (1,2):
        selected=np.flatnonzero(mids==layer)
        Kl=assemble_stiffness(basis.with_elements(selected),C[selected])
        layer_energy[str(layer)]=float(.5*coef@(Kl@coef))
        residual=Kl@coef-(nodal_load if layer==1 else 0)
        nodal=np.column_stack([residual[d] for d in basis.split_indices()])
        layer_residual[layer]=nodal
        bottom,top=nodal[end0].sum(0),nodal[endL].sum(0)
        expected=float(2*MATERIALS[layer]['nu']*constants[layer-1,0]*np.pi*(radii[layer]**2-radii[layer-1]**2))
        error=float(max(abs(bottom[2]/(-expected)-1),abs(top[2]/expected-1)))
        end_errors.append(error)
        end_reactions[str(layer)]=dict(z0_N=bottom.tolist(),zL_N=top.tolist(),analytic_zL_axial_N=expected,axial_relative_error=error)
    # Interface equilibrium is assessed away from prescribed end nodes.
    interface_node=np.isclose(np.linalg.norm(arrays['xyz'][:,:2],axis=1),p['interface_radius_mm'],atol=1e-7)&~(end0|endL)
    weak_jump=float(np.linalg.norm((layer_residual[1]+layer_residual[2])[interface_node])/np.linalg.norm(layer_residual[1][interface_node]))
    traces=[];Utraces=[];mattraces=[];coords=[]
    for side in (0,1):
        fb=InteriorFacetBasis(mesh,basis.elem,facets=interface,intorder=5,side=side)
        field=fb.interpolate(coef)
        traces.append(stress_to_tensor(np.moveaxis(apply_stiffness(C[fb.tind],engineering_strain(field)),0,-1)))
        Utraces.append(np.moveaxis(np.asarray(field),0,-1));mattraces.append(mids[fb.tind]);coords.append(np.moveaxis(np.asarray(fb.global_coordinates()),0,-1))
    if not np.allclose(coords[0],coords[1],atol=1e-10):raise ValueError('Interface quadrature point mismatch')
    x=coords[0];er=x.copy();er[...,2]=0;er/=np.linalg.norm(er,axis=-1)[...,None]
    et=np.stack([-er[...,1],er[...,0],np.zeros(er.shape[:2])],axis=-1)
    trace1=np.where((mattraces[0]==1)[:,None,None,None],traces[0],traces[1])
    trace2=np.where((mattraces[0]==2)[:,None,None,None],traces[0],traces[1])
    traction1=np.einsum('fqij,fqj->fqi',trace1,er);traction2=np.einsum('fqij,fqj->fqi',trace2,er)
    qweight=np.asarray(fb.dx)
    exact_if1=analytic_fields(np.array([[p['interface_radius_mm'],0,0]]),np.array([1]))[1][0]
    exact_if2=analytic_fields(np.array([[p['interface_radius_mm'],0,0]]),np.array([2]))[1][0]
    interface_pressure=-float(exact_if1[0,0])
    def surface_average(value):return float(np.sum(value*qweight)/qweight.sum())
    rr1=np.sum(traction1*er,axis=-1);rr2=np.sum(traction2*er,axis=-1)
    hoop1=np.einsum('fqi,fqij,fqj->fq',et,trace1,et);hoop2=np.einsum('fqi,fqij,fqj->fq',et,trace2,et)
    trace_jump=float(np.sqrt(np.sum(np.sum((traction1-traction2)**2,axis=-1)*qweight)/qweight.sum())/interface_pressure)
    trace_summary=dict(displacement_jump_max_mm=float(np.max(np.linalg.norm(Utraces[0]-Utraces[1],axis=-1))),
        weak_interface_force_relative_residual=weak_jump,
        traction_trace_jump_relative_surface_L2=trace_jump,
        radial_stress_layer1_mean_MPa=surface_average(rr1),radial_stress_layer2_mean_MPa=surface_average(rr2),
        analytic_radial_stress_MPa=-interface_pressure,
        radial_trace_mean_max_relative_error=max(abs(surface_average(rr1)/(-interface_pressure)-1),abs(surface_average(rr2)/(-interface_pressure)-1)),
        hoop_stress_layer1_mean_MPa=surface_average(hoop1),hoop_stress_layer2_mean_MPa=surface_average(hoop2),
        analytic_hoop_layer1_MPa=float(exact_if1[1,1]),analytic_hoop_layer2_MPa=float(exact_if2[1,1]),
        note='Separate one-sided FE traces; strong pointwise traction is not enforced by conforming displacement FEM; weak force equilibrium is checked independently')
    inner_ur=float(analytic_fields(np.array([[p['inner_radius_mm'],0,0]]))[0][0,0])
    exact_energy=float(.5*p['pressure_MPa']*2*np.pi*p['inner_radius_mm']*p['length_mm']*inner_ur)
    work=diag['traction_work_Nmm']+diag['prescribed_boundary_reaction_work_Nmm']
    row=dict(parameters=p,synthetic_materials=MATERIALS,analytic_A_B=constants.tolist(),
        n_nodes=diag['n_nodes'],n_elements=diag['n_elements'],seconds=time.perf_counter()-started,
        displacement_relative_L2=u_l2,stress_tensor_relative_volume_L2_by_layer=stress_l2,
        interface=trace_summary,end_reactions_by_layer=end_reactions,max_layer_axial_reaction_relative_error=max(end_errors),
        strain_energy_by_layer_Nmm=layer_energy,strain_energy_analytic_Nmm=exact_energy,
        energy_analytic_relative_error=abs(diag['strain_energy_Nmm']/exact_energy-1),
        boundary_work_energy_relative_error=abs(work/(2*diag['strain_energy_Nmm'])-1),diagnostics=diag)
    # All exported nodal fields retain original Gmsh order. Stress remains per GP.
    U=np.empty_like(arrays['displacement']);U[to_original]=arrays['displacement']
    np.savez_compressed(out/'solution.npz',**source,displacement=U,stress=arrays['stress'],quadrature_xyz=arrays['quadrature_xyz'],quadrature_weight=w)
    (out/'solution.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
    return row
