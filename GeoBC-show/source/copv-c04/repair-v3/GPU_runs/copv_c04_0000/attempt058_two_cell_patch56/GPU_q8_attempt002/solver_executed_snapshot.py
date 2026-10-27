"""C03 detailed structural mesh/order/closure study: 3-D orthotropic elasticity on PETSc CUDA.

The accepted configuration defines ideal bonds, excluded seals and equivalent
caps explicitly. Raw stress is stored by original element/GP, never smoothed
across materials. Linear load/unload states use exact field superposition. The default detailed
benchmark is P1 on the accepted straight geometry; original TET10 midside
output values are affine interpolation, not independent P2 unknowns.
"""
import argparse,gc,hashlib,json,os,sysconfig,time,subprocess
from pathlib import Path
os.environ.setdefault('PETSC_OPTIONS','-use_gpu_aware_mpi 0')
_cfg=[sysconfig.get_config_vars()]
try:
    from setuptools._distutils import sysconfig as _dc
    _cfg.append(_dc.get_config_vars())
except ImportError:pass
for cfg in _cfg:cfg.update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from mpi4py import MPI
from petsc4py import PETSc
import basix,basix.ufl,ufl,dolfinx
from dolfinx import mesh,fem
from dolfinx.fem import petsc as fp
from geobc.copv_material import engineering_stiffness
from geobc.copv_structured_mesh import analytic_segments
from geobc.copv_geometry import DEFAULT_PARAMETERS

def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');tmp.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stiffness(material):
    if 'E_MPa' in material:
        E=material['E_MPa'];nu=material['nu'];G=E/(2*(1+nu));return engineering_stiffness(E1=E,E2=E,E3=E,G12=G,G13=G,G23=G,nu12=nu,nu13=nu,nu23=nu)
    return engineering_stiffness(**{k:material[k] for k in ['E1','E2','E3','G12','G13','G23','nu12','nu13','nu23']})
def eng(e):return ufl.as_vector([e[0,0],e[1,1],e[2,2],2*e[1,2],2*e[0,2],2*e[0,1]])

def compact_material_frame(domain,tangent_field,frame_parameters):
    x=ufl.SpatialCoordinate(domain);r=ufl.sqrt(x[0]*x[0]+x[1]*x[1]);er=ufl.as_vector([x[0]/r,x[1]/r,0]);et=ufl.as_vector([-x[1]/r,x[0]/r,0]);m=-tangent_field[0]*er-ufl.as_vector([0,0,tangent_field[1]]);n=ufl.cross(m,et);mask,ca,sa=frame_parameters[0],frame_parameters[1],frame_parameters[2];a=ca*m+sa*et;b=-sa*m+ca*et;rotation=ufl.as_matrix([[a[i],b[i],n[i]] for i in range(3)]);return mask*rotation+(1-mask)*ufl.Identity(3)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--assembly-quadrature',type=int,choices=[4,6,8],default=6);ap.add_argument('--repair-preflight',type=Path,required=True);ap.add_argument('--mesh',type=Path,required=True);ap.add_argument('--surface',type=Path);ap.add_argument('--review',type=Path,required=True);ap.add_argument('--contract',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--degree',type=int,choices=[1,2],default=2);ap.add_argument('--preconditioner',choices=['full-gamg','p1-gamg'],default='full-gamg');ap.add_argument('--cap-convention',choices=['circular','discrete'],default='discrete');ap.add_argument('--save-raw-gp',action='store_true');ap.add_argument('--save-gp-local',action='store_true');ap.add_argument('--quarantined-research',action='store_true');ap.add_argument('--research-validation-only',action='store_true');ap.add_argument('--validate-native-schema',action='store_true');args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();phase={}
    def status(name,**extra):
        phase.update({'pid':os.getpid(),'phase':name,'seconds':time.perf_counter()-start,**extra});write(out/'progress.json',phase);print(json.dumps(phase),flush=True)
    if not PETSc.Sys.hasExternalPackage('cuda'):raise RuntimeError('CUDA required; CPU fallback refused')
    if not args.quarantined_research:raise RuntimeError('This research-only executable requires --quarantined-research')
    if args.degree!=2 or args.preconditioner!='p1-gamg' or not args.save_gp_local:raise RuntimeError('Frozen research pilot requires P2 / exact CUDA P1 coarse space and original local GP outputs')
    args.research_validation_only=True;args.validate_native_schema=True
    contract=json.loads(args.contract.read_text());review=json.loads(args.review.read_text());parts={p['part_id']:p for p in review['parts'] if p['structural']};source=np.load(args.mesh);xyz=source['xyz'];cells=source['cells'];pids=source['cell_part_ids'];nc=len(cells)
    params=contract['geometry_parameters_mm'];end=params['straight_length_mm']/2+params['dome_depth_mm']+params['neck_length_mm'];liner_end=end-3.;upper_z=end+22.;lower_z=-(end+15.);upper_bore=contract['load']['upper_bore_radius_mm'];lower_bore=contract['load']['lower_bore_radius_mm']
    if not args.research_validation_only:
        if not review.get('C04_production_gate',{}).get('passed',False):raise RuntimeError('C04 production mesh gate not passed; no FE label allowed')
        if not args.save_gp_local:raise RuntimeError('C04 accepted target schema requires original-cell local GP stress')
        if args.degree!=2 or args.preconditioner!='p1-gamg':raise RuntimeError('C04 accepted targets require the frozen P2 / exact embedded CUDA P1 coarse protocol')
        if review.get('parameters_mm')!=params:raise RuntimeError('C04 geometry/contract parameter identity mismatch')
    if contract.get('mesh_sha256') and sha(args.mesh)!=contract['mesh_sha256']:raise RuntimeError('C04 mesh/contract identity mismatch')
    repair_preflight=json.loads(args.repair_preflight.read_text())
    if not repair_preflight['candidate_FE_validation_allowed'] or repair_preflight['production_eligible']:raise RuntimeError('Repair candidate has no geometry/topology/load validation permission')
    if repair_preflight['mesh_sha256']!=sha(args.mesh):raise RuntimeError('Repair preflight mesh identity mismatch')
    candidate=args.repair_preflight.parent.parent
    geometry_path=candidate/'geometry_check/review.json'
    geometry=json.loads(geometry_path.read_text())
    if (candidate/'local_geometry_hold.json').exists() or not geometry['passed']:raise RuntimeError('Actual CAD geometry hold')
    if sha(geometry_path)!=repair_preflight['local_geometry_review_sha256'] or geometry['mesh_sha256']!=sha(args.mesh):raise RuntimeError('Geometry review identity mismatch')
    if review.get('parameters_mm')!=params:raise RuntimeError('Geometry/contract parameter identity mismatch')
    def resources(stage,host_required=0,gpu_required=24576):
        info={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemAvailable:','MemTotal:'))}
        data=subprocess.check_output(['nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip().splitlines()
        if len(data)!=1:raise RuntimeError('Expected the frozen single GPU')
        free,util=[int(v.strip()) for v in data[0].split(',')];record=dict(stage=stage,host_available_bytes=info['MemAvailable'],host_required_bytes=int(host_required),GPU_free_MiB=free,GPU_utilization_percent=util,GPU_required_MiB=int(gpu_required))
        write(out/('resources_'+stage+'.json'),record)
        if info['MemAvailable']<host_required or free<gpu_required or util>10:raise RuntimeError('Resources occupied or insufficient; preserve candidate and retry in a NEW run after resources become available')
    # CPU guard before mesh/form allocation, then actual sparse storage before CUDA.
    resources('before_domain',host_required=8*1024**3+nc*args.assembly_quadrature*1024)
    native_path=args.repair_preflight.parent/'native_boundary_vectors.npz'
    closure_path=args.repair_preflight.parent/'quadratic_load_closure.json'
    if sha(native_path)!=repair_preflight['native_vectors_sha256'] or sha(closure_path)!=repair_preflight['closure_sha256']:raise RuntimeError('Quadratic load preflight identity mismatch')
    independent_native=np.load(native_path);independent_boundary=json.loads(closure_path.read_text())
    if not independent_boundary['passed'] or not np.array_equal(independent_native['node_ids'],source['node_ids']):raise RuntimeError('Quadratic load preflight/node identity failed')
    if review['connected_components']!=1 or source['min_det_jac'].min()<=0:raise RuntimeError('Unaccepted mesh')
    # Native Gmsh to Basix TET10 reference mapping, as in the R02 solver.
    edges=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)];top=basix.topology(basix.CellType.tetrahedron);perm=[0,1,2,3]+[4+next(i for i,e in enumerate(edges) if set(e)==set(edge)) for edge in top[1]]
    element=basix.ufl.element('Lagrange','tetrahedron',2,shape=(3,),lagrange_variant=basix.LagrangeVariant.equispaced);domain=mesh.create_mesh(MPI.COMM_SELF,cells[:,perm].astype(np.int64),ufl.Mesh(element),xyz);V=fem.functionspace(domain,basix.ufl.element('Lagrange','tetrahedron',args.degree,shape=(3,),lagrange_variant=basix.LagrangeVariant.equispaced));coords=V.tabulate_dof_coordinates();dist,node_order=cKDTree(xyz).query(coords)
    if dist.max()>1e-8 or len(np.unique(node_order))!=len(coords):raise RuntimeError('Original node bijection failed')
    cell_order=domain.topology.original_cell_index;part_local=pids[cell_order];qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,args.assembly_quadrature);nq=len(qp);scalar=basix.create_element(basix.ElementFamily.P,basix.CellType.tetrahedron,2,basix.LagrangeVariant.equispaced);tab=scalar.tabulate(1,qp);geommap=domain.geometry.dofmap
    Tspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(2,),degree=args.assembly_quadrature,scheme='default'));T=fem.Function(Tspace);FrameSpace=fem.functionspace(domain,basix.ufl.element('DG','tetrahedron',0,shape=(3,)));FrameParameters=fem.Function(FrameSpace);Q=compact_material_frame(domain,T,FrameParameters);Cspace=fem.functionspace(domain,basix.ufl.element('DG','tetrahedron',0,shape=(6,6)));C=fem.Function(Cspace);matrices={name:stiffness(value) for name,value in contract['materials'].items()}
    rz,normal=analytic_segments(params,(4,2500,16,2500,4));tree=cKDTree(rz);tangent=np.c_[-normal[:,1],normal[:,0]];p_by_id={pid:matrices[p['material']] for pid,p in parts.items()};angle_by_id={pid:float(p.get('angle',0)) for pid,p in parts.items()};carbon=set(pid for pid,p in parts.items() if p['material']=='carbon');frame_error=0.;status('material_rotation',nodes=len(coords),geometry_nodes=len(xyz),cells=nc,dofs=3*len(coords),displacement_degree=args.degree,quadrature_points=nq)
    def chunk_geometry(ids):
        nodes=domain.geometry.x[geommap[ids]];position=np.einsum('qi,eij->eqj',tab[0,:,:,0],nodes);deriv=np.stack([np.einsum('qi,eij->eqj',tab[j,:,:,0],nodes) for j in [1,2,3]],axis=-1);return position,np.linalg.det(deriv)
    def material_axes(positions,ids):
        axes=np.broadcast_to(np.eye(3),(len(ids),positions.shape[1],3,3)).copy();sel=np.isin(part_local[ids],list(carbon))
        if sel.any():
            x=positions[sel];r=np.linalg.norm(x[...,:2],axis=-1);_,near=tree.query(np.stack([r,x[...,2]],axis=-1));t=tangent[near];er=np.concatenate([x[...,:2]/r[...,None],np.zeros(r.shape+(1,))],axis=-1);et=np.stack([-er[...,1],er[...,0],np.zeros_like(r)],axis=-1);m=-t[...,0,None]*er;m[...,2]-=t[...,1];n=np.cross(m,et);angle=np.deg2rad([angle_by_id[int(v)] for v in part_local[ids][sel]])[:,None,None];ca=np.cos(angle);sa=np.sin(angle);axes[sel]=np.stack([ca*m+sa*et,-sa*m+ca*et,n],axis=-1)
        return axes
    for first in range(0,nc,4000):
        ids=np.arange(first,min(first+4000,nc),dtype=np.int32);positions,_=chunk_geometry(ids);axes=material_axes(positions,ids)
        frame_error=max(frame_error,float(abs(np.linalg.det(axes)-1).max()))
        radii=np.linalg.norm(positions[...,:2],axis=-1);_,nearest=tree.query(np.stack([radii,positions[...,2]],axis=-1));tangent_values=tangent[nearest]
        for j,e in enumerate(ids):
            pid=int(part_local[e]);angle=np.deg2rad(angle_by_id[pid]);T.x.array.reshape(-1,2)[Tspace.dofmap.cell_dofs(e)]=tangent_values[j];FrameParameters.x.array.reshape(-1,3)[FrameSpace.dofmap.cell_dofs(e)]=np.array([[float(pid in carbon),np.cos(angle),np.sin(angle)]]);C.x.array.reshape(-1,36)[Cspace.dofmap.cell_dofs(e)]=p_by_id[pid].reshape(1,36)
    if frame_error>1e-10:raise RuntimeError('Invalid GP material frame')
    # Exterior facet tags use the documented pressure cavity and actual geometry.
    domain.topology.create_connectivity(2,3);boundary=mesh.exterior_facet_indices(domain.topology);fg=mesh.entities_to_geometry(domain,2,boundary,False);fx=domain.geometry.x[fg];centroid=fx.mean(axis=1)
    if fx.shape[1]==6:centroid=(-fx[:,:3].sum(1)+4*fx[:,3:].sum(1))/9
    r=np.linalg.norm(centroid[:,:2],axis=1);z=centroid[:,2];inner=np.interp(abs(z),rz[rz[:,1]>=0,1],rz[rz[:,1]>=0,0]);inner=np.where(z>liner_end,upper_bore,np.where(z<-liner_end,lower_bore,inner));values=np.full(len(boundary),4,np.int32);values[r<=inner+.30]=1;values[np.isclose(z,upper_z,atol=1e-5)]=2;values[np.isclose(z,lower_z,atol=1e-5)]=3
    for value in [1,2,3]:
        if not np.any(values==value):raise RuntimeError('Missing pressure/cap boundary '+str(value))
    order=np.argsort(boundary);facet_tags=mesh.meshtags(domain,2,boundary[order].astype(np.int32),values[order]);dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':args.assembly_quadrature});ds=ufl.Measure('ds',domain=domain,subdomain_data=facet_tags,metadata={'quadrature_degree':4});normal_ufl=ufl.FacetNormal(domain)
    cap_areas={name:float(fem.assemble_scalar(fem.form(1*ds(tag)))) for name,tag in [('top',2),('bottom',3)]};forces={name:float(np.pi*contract['load'][name+'_bore_radius_mm']**2) for name in ['upper','lower']}
    circular_forces=forces.copy();polygon_forces={}
    corners=np.unique(cells[:,:4])
    for name,zend in [('upper',upper_z),('lower',lower_z)]:
        radius=contract['load'][name+'_bore_radius_mm'];xx=xyz[corners];selected=xx[np.isclose(xx[:,2],zend,atol=1e-7)&np.isclose(np.linalg.norm(xx[:,:2],axis=1),radius,atol=1e-7),:2]
        if len(selected)<3:raise RuntimeError('Cannot identify the actual discrete bore polygon')
        selected=selected[np.argsort(np.arctan2(selected[:,1],selected[:,0]))];polygon_forces[name]=float(.5*abs(np.sum(selected[:,0]*np.roll(selected[:,1],-1)-selected[:,1]*np.roll(selected[:,0],-1))))
    if args.cap_convention=='discrete':forces=polygon_forces.copy()
    forces=independent_boundary['cap_forces_N'].copy()
    for cap in ['top','bottom']:
        if abs(cap_areas[cap]/independent_boundary['cap_surface_areas_mm2'][cap]-1)>1e-9:raise RuntimeError('Independent curved cap-area identity failed')
    write(out/'cap_definition.json',dict(selected='quadratic_opening_edge_integral',circular_peak_N=circular_forces,corner_polygon_peak_N=polygon_forces,quadratic_opening_peak_N=forces))
    u=ufl.TrialFunction(V);v=ufl.TestFunction(V)
    def local_strain(w):return eng(Q.T*ufl.sym(ufl.grad(w))*Q)
    bilinear=fem.form(ufl.dot(local_strain(v),C*local_strain(u))*dx);linear=fem.form(-ufl.dot(normal_ufl,v)*ds(1)+(forces['upper']/cap_areas['top'])*v[2]*ds(2)-(forces['lower']/cap_areas['bottom'])*v[2]*ds(3))
    status('assembling',material_frame_det_max_error=frame_error,cap_areas_mm2=cap_areas,cap_forces_N=forces)
    A=fp.create_matrix(bilinear);preallocated=A.getInfo();status('matrix_preallocated',matrix_info=preallocated)
    estimated_sparse_bytes=int(preallocated['nz_allocated'])*(8+PETSc.IntType().itemsize)
    # Sparse storage is already preallocated. Reserve the actual coefficient
    # pack plus assembly scratch, rather than another complete sparse matrix.
    estimated_coefficient_pack_bytes=nc*(39+2*nq)*8
    resources('before_assembly',host_required=estimated_coefficient_pack_bytes+2*1024**3)
    write(out/'material_frame_storage.json',dict(method='Two quadrature tangent components and three DG0 mask/cos/sin parameters; spatial radial basis reconstructed exactly in UFL',assembly_points=nq,full_frame_value_bytes=nc*nq*9*8,compact_frame_value_bytes=nc*(nq*2+3)*8,estimated_coefficient_pack_bytes=estimated_coefficient_pack_bytes,preallocated_sparse_bytes=estimated_sparse_bytes,scope='Same axisymmetric material-frame values at the same quadrature points, not an interpolated or changed material model.'))
    fp.assemble_matrix(A,bilinear);A.assemble();status('matrix_assembled',matrix_info=A.getInfo());bhost=fp.assemble_vector(linear);bhost.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,mode=PETSc.ScatterMode.REVERSE);F=bhost.array.copy();status('loads_assembled')
    csr_i,csr_j,csr_v=A.getValuesCSR();actual_sparse_bytes=csr_i.nbytes+csr_j.nbytes+csr_v.nbytes;del csr_i,csr_j,csr_v
    resources('before_CUDA',host_required=2*1024**3,gpu_required=max(24576,int((2.0*actual_sparse_bytes+3*1024**3)/1024**2)))
    write(out/'quadrature_protocol.json',dict(assembly_degree=args.assembly_quadrature,output_degree=4,output_points=14,load_degree=4,actual_CSR_bytes=actual_sparse_bytes,scope='Assembly sensitivity candidate; original GP/node IDs and 14-point field layout preserved.'))
    corner_ids=np.unique(cells[:,:4]);corner_map=np.full(len(xyz),-1,np.int64);corner_map[corner_ids]=np.arange(len(corner_ids));bottom=np.flatnonzero(np.isclose(coords[:,2],lower_z,atol=1e-7)&(corner_map[node_order]>=0));radius=np.linalg.norm(coords[bottom,:2],axis=1).max();chosen=[int(bottom[np.argmin(np.linalg.norm(coords[bottom]-point,axis=1))]) for point in [[radius,0,lower_z],[-radius,0,lower_z],[0,radius,lower_z]]];a,b,c=chosen
    if len(set(chosen))!=3:raise RuntimeError('Non-distinct gauge nodes')
    fixed=np.array([3*a,3*a+1,3*a+2,3*b+1,3*b+2,3*c+2],np.int32);gaugecoords=coords[chosen];centerpoint=coords.mean(0);raw=[]
    for axis in np.eye(3):raw.append(np.tile(axis,(len(coords),1)).ravel())
    for axis in np.eye(3):raw.append(np.cross(np.broadcast_to(axis,coords.shape),coords-centerpoint).ravel())
    if np.linalg.matrix_rank(np.column_stack([v0[fixed] for v0 in raw]))!=6:raise RuntimeError('3-2-1 gauge does not remove exactly six rigid motions')
    gauge_rows=[]
    for dof in fixed:
        ind,val=A.getRow(int(dof));gauge_rows.append((ind.copy(),val.copy()))
    status('applying_gauge',fixed_dofs=fixed.tolist());A.zeroRowsColumns(fixed,diag=1.);bhost.array[fixed]=0.;A.setBlockSize(3);near=[]
    for vec in raw:
        vec[fixed]=0
        for previous in near:vec-=np.dot(previous,vec)*previous
        vec/=np.linalg.norm(vec);near.append(vec)
    if args.degree==1 and args.preconditioner=='p1-gamg':raise RuntimeError('Embedded P1 coarse path is only for a P2 solve')
    if args.preconditioner == 'p1-gamg':
        # Explicit p-coarsening avoids building a dense six-mode aggregation
        # hierarchy directly on four million P2 displacement unknowns. The fine
        # operator is unchanged; the coarse space is the exact embedded P1 space.
        status('building_P1_preconditioner')
        endpoints=np.full((len(xyz),2),-1,np.int64)
        for j,(i0,i1) in enumerate(edges):endpoints[cells[:,4+j]]=cells[:,[i0,i1]]
        endpoints[corner_ids]=corner_ids[:,None];ep=corner_map[endpoints[node_order]]
        if np.any(ep<0):raise RuntimeError('Incomplete P2-to-P1 edge map')
        coarse_fixed=np.array([3*corner_map[node_order[int(d)//3]]+int(d)%3 for d in fixed],dtype=PETSc.IntType)
        rr=np.repeat(np.arange(3*len(xyz)),2);cc=(3*np.repeat(ep,3,axis=0)+np.tile(np.arange(3),len(xyz))[:,None]).ravel();vv=np.full(len(rr),.5);keep=~np.isin(cc,coarse_fixed);transfer=sp.coo_matrix((vv[keep],(rr[keep],cc[keep])),shape=(3*len(xyz),3*len(corner_ids))).tocsr();transfer.sum_duplicates();Pcpu=PETSc.Mat().createAIJ(size=transfer.shape,csr=(transfer.indptr.astype(PETSc.IntType),transfer.indices.astype(PETSc.IntType),transfer.data),comm=MPI.COMM_SELF)
        coarse=A.PtAP(Pcpu);coarse.zeroRowsColumns(coarse_fixed,diag=1.);coarse.setBlockSize(3);coarse.setVecType('cuda');coarse_info=coarse.getInfo();status('P1_preconditioner_projected',coarse_dofs=coarse.getSize()[0],coarse_nonzeros=coarse_info['nz_used'],coarse_matrix_type=coarse.getType(),coarse_backend='CUDA P1 GAMG and CUDA transfer vectors');Pcuda=Pcpu.convert('aijcusparse');Kcoarse=coarse.convert('aijcusparse')
        if Pcuda.handle!=Pcpu.handle:Pcpu.destroy()
        del transfer,rr,cc,vv,keep,endpoints,ep;gc.collect()
        coarse_modes=[];coarse_xyz=xyz[corner_ids];coarse_center=coarse_xyz.mean(0)
        for i in range(6):
            axis=np.eye(3)[i%3];values=np.broadcast_to(axis,coarse_xyz.shape).copy() if i<3 else np.cross(np.broadcast_to(axis,coarse_xyz.shape),coarse_xyz-coarse_center);values=values.ravel();values[coarse_fixed]=0
            for old in coarse_modes:values-=np.dot(old,values)*old
            values/=np.linalg.norm(values);coarse_modes.append(values)
        coarse_vecs=[]
        for values in coarse_modes:
            pv=PETSc.Vec().createSeq(len(values),comm=MPI.COMM_SELF);pv.setType('cuda');pv.array[:]=values;coarse_vecs.append(pv)
        coarse_ns=PETSc.NullSpace().create(vectors=coarse_vecs,comm=MPI.COMM_SELF);Kcoarse.setNearNullSpace(coarse_ns);del coarse_modes;gc.collect()
    status('converting_CUDA',source_matrix_type=A.getType());gpu=A.convert('aijcusparse');matrix_info=A.getInfo()
    # PETSc may return the same underlying Mat when the source is already the
    # requested CUDA type. Destroying that alias invalidates the solver matrix.
    if gpu.handle!=A.handle:A.destroy()
    rhs=PETSc.Vec().createSeq(len(F),comm=MPI.COMM_SELF);rhs.setType('cuda');rhs.array[:]=bhost.array;answer=rhs.duplicate();answer.set(0);bhost.destroy();nearvecs=[]
    for v0 in near:
        pv=rhs.duplicate();pv.array[:]=v0;nearvecs.append(pv)
    ns=PETSc.NullSpace().create(vectors=nearvecs,comm=MPI.COMM_SELF);gpu.setNearNullSpace(ns);del near,raw;gc.collect()
    ksp=PETSc.KSP().create(MPI.COMM_SELF);ksp.setOptionsPrefix('copv_');ksp.setOperators(gpu);ksp.setType('fgmres');ksp.setGMRESRestart(60);pc=ksp.getPC();opts=PETSc.Options()
    if args.preconditioner == 'p1-gamg':
        pc.setType('mg');pc.setMGLevels(2);pc.setMGInterpolation(1,Pcuda);pc.setMGRestriction(1,Pcuda)
        smooth=pc.getMGSmoother(1);smooth.setOperators(gpu);smooth.setType('chebyshev');smooth.getPC().setType('jacobi');smooth.setTolerances(max_it=5)
        cksp=pc.getMGCoarseSolve();cksp.setOptionsPrefix('copv_p1_');cksp.setOperators(Kcoarse);cksp.setType('preonly');cksp.getPC().setType('gamg')
        for key,value in {'pc_gamg_type':'agg','pc_gamg_agg_nsmooths':0,'pc_gamg_threshold':0.,'pc_gamg_aggressive_coarsening':1,'pc_gamg_aggressive_square_graph':False,'mg_levels_ksp_type':'chebyshev','mg_levels_pc_type':'jacobi','mg_levels_ksp_max_it':3,'mg_coarse_ksp_type':'preonly','mg_coarse_pc_type':'svd'}.items():opts['copv_p1_'+key]=value
        cksp.setFromOptions()
        backend_description='CPU FFCx assembly; CUDA P2 FGMRES/Jacobi-Chebyshev and exact P1 transfer; CUDA P1 GAMG coarse correction; tiny CPU SVD at final coarse level, explicit in KSP view'

    else:
        pc.setType('gamg')
        opts['copv_pc_gamg_type']='agg';opts['copv_pc_gamg_agg_nsmooths']=0;opts['copv_pc_gamg_threshold']=0.;opts['copv_pc_gamg_aggressive_coarsening']=1;opts['copv_pc_gamg_aggressive_square_graph']=False
        opts['copv_mg_levels_ksp_type']='chebyshev';opts['copv_mg_levels_pc_type']='jacobi';opts['copv_mg_levels_ksp_max_it']=3;opts['copv_mg_coarse_ksp_type']='preonly';opts['copv_mg_coarse_pc_type']='svd'
        backend_description='CPU FFCx assembly; PETSc CUDA FGMRES/GAMG on full displacement operator, unsmoothed aggregation and Jacobi/Chebyshev levels; small coarsest SVD on CPU, hierarchy recorded in KSP view'
    ksp.setTolerances(rtol=1e-10,max_it=20000);ksp.setFromOptions();iteration=[]
    def monitor(solver,it,norm):
        if it%100==0:
            rec=dict(iteration=int(it),reported_residual=float(norm),seconds=time.perf_counter()-start);iteration.append(rec);status('solving_peak',**rec)
        if it and it%100==0:
            checkpoint=out/'peak_iterate.npy';tmp=out/'peak_iterate.tmp.npy';np.save(tmp,solver.getSolution().array);tmp.replace(checkpoint)
            write(out/'checkpoint.json',dict(mesh_sha256=sha(args.mesh),contract_sha256=sha(args.contract),iteration=int(it),reported_residual=float(norm)))
    checkpoint=out/'peak_iterate.npy'
    if checkpoint.exists() and (out/'checkpoint.json').exists():
        record=json.loads((out/'checkpoint.json').read_text())
        if record['mesh_sha256']==sha(args.mesh) and record['contract_sha256']==sha(args.contract):answer.array[:]=np.load(checkpoint);ksp.setInitialGuessNonzero(True);status('restored_iterate',previous_iteration=record['iteration'])
    ksp.setMonitor(monitor);status('solving_peak',matrix_nonzeros=matrix_info['nz_used'],matrix_type=gpu.getType(),vector_type=rhs.getType());solve_start=time.perf_counter();ksp.solve(rhs,answer);solve_seconds=time.perf_counter()-solve_start
    if ksp.getConvergedReason()<=0:raise RuntimeError('CUDA peak solve failed '+str(ksp.getConvergedReason()))
    viewer=PETSc.Viewer().createASCII(str(out/'ksp_view.txt'),comm=MPI.COMM_SELF);ksp.view(viewer);viewer.destroy();write(out/'petsc_options.json',PETSc.Options().getAll())
    primary_iterations=int(ksp.getIterationNumber());refinement_iterations=0;refinement=[]
    # Recompute the actual algebraic residual. Restarted correction solves
    # address loss of recursive-CG residual accuracy on thin/tiny elements.
    for cycle in range(4):
        actual=rhs.duplicate();gpu.mult(answer,actual);actual.axpy(-1,rhs);norm=float(actual.norm()/rhs.norm());refinement.append(dict(cycle=cycle,true_relative_residual=norm));status('checking_true_residual',refinement_cycle=cycle,true_relative_residual=norm)
        if norm<=1e-7 or cycle==3:actual.destroy();break
        actual.scale(-1);correction=rhs.duplicate();correction.set(0);ksp.cancelMonitor();ksp.setInitialGuessNonzero(False);ksp.setTolerances(rtol=1e-9);ksp.solve(actual,correction)
        if ksp.getConvergedReason()<=0:raise RuntimeError('Residual-correction solve failed')
        refinement_iterations+=int(ksp.getIterationNumber());answer.axpy(1,correction);actual.destroy();correction.destroy()
    solve_seconds=time.perf_counter()-solve_start;write(out/'residual_refinement.json',refinement)
    peak_iterations=primary_iterations+refinement_iterations;peak_reason=int(ksp.getConvergedReason());peak=answer.array.copy();residual=rhs.duplicate();gpu.mult(answer,residual);energy=float(.5*answer.dot(residual));residual.axpy(-1,rhs);true_residual=float(residual.norm()/rhs.norm());work=float(answer.dot(rhs))
    if true_residual>1e-7:raise RuntimeError('True residual gate failed '+str(true_residual))
    reactions=np.array([np.dot(value,peak[index])-F[dof] for dof,(index,value) in zip(fixed,gauge_rows)]);totalforce=F.reshape(-1,3).sum(0);totalmoment=np.cross(coords,F.reshape(-1,3)).sum(0);reaction_vectors=np.zeros((len(coords),3));reaction_vectors.reshape(-1)[fixed]=reactions;force_balance=totalforce+reaction_vectors.sum(0);moment_balance=totalmoment+np.cross(coords,reaction_vectors).sum(0);scale=np.pi*params['inner_radius_mm']**2
    if np.linalg.norm(force_balance)/scale>1e-5 or np.linalg.norm(moment_balance)/(scale*(upper_z-lower_z))>1e-5:raise RuntimeError('Global force/moment gate failed')
    if abs(2*energy/work-1)>1e-7:raise RuntimeError('Energy/work gate failed')
    write(out/'peak_algebra.json',dict(true_relative_residual=true_residual,iterations=peak_iterations,reason=int(ksp.getConvergedReason()),strain_energy_Nmm=energy,external_work_Nmm=work,solve_seconds=solve_seconds,matrix_nonzeros=matrix_info['nz_used'],force_balance_N=force_balance.tolist(),moment_balance_Nmm=moment_balance.tolist()))
    original=np.empty((len(xyz),3));original[node_order]=peak.reshape(-1,3)
    if args.degree==1:
        for j,(i0,i1) in enumerate(edges):original[cells[:,4+j]]=.5*(original[cells[:,i0]]+original[cells[:,i1]])
    np.savez_compressed(out/'nodal_fields.npz',U_mm=original,node_ids=source['node_ids'],gauge_dofs_FEniCSx=fixed,gauge_original_nodes=node_order[chosen],gauge_coordinates_mm=gaugecoords,gauge_reactions_N=reactions)
    if not args.research_validation_only or args.validate_native_schema:
        if args.degree!=2:raise RuntimeError('Native label schema validation requires P2')
        # Preserve consistent peak loads in the original geometry node order.
        # Independent native vectors integrate the true quadratic face map.
        loads=[independent_native[key] for key in ['pressure_load_N','upper_cap_load_N','lower_cap_load_N']]
        actual_load=np.zeros_like(original);actual_load[node_order]=F.reshape(-1,3);actual_reaction=np.zeros_like(original);actual_reaction[node_order]=reaction_vectors;locked=np.zeros(original.shape,bool)
        for dof in fixed:locked[node_order[dof//3],dof%3]=True
        native_load_error=float(np.linalg.norm(sum(loads)-actual_load)/np.linalg.norm(actual_load))
        if native_load_error>1e-10:raise RuntimeError('Native consistent-load identity failed '+str(native_load_error))
        np.savez_compressed(out/'native_load_and_bc.npz',node_ids=source['node_ids'],consistent_load_N=actual_load,pressure_load_N=loads[0],upper_cap_load_N=loads[1],lower_cap_load_N=loads[2],reaction_N=actual_reaction,locked_components=locked)
        write(out/'native_load_identity.json',dict(passed=True,relative_difference=native_load_error,node_order='original Gmsh node_ids',pressure_MPa=1.,units='mm N MPa',basis='Independent native quadratic Tri6 coordinate/basis quadrature; actual quadratic bore-edge cap area; compared with true FEniCSx-assembled RHS'))
    # Independent half-load solve checks linearity; unloading states are exact
    # linear superposition, not claimed as nonlinear plastic path integrations.
    ksp.cancelMonitor();rhs.scale(.5);answer.set(0);ksp.setInitialGuessNonzero(False);ksp.solve(rhs,answer)
    if ksp.getConvergedReason()<=0:raise RuntimeError('Half-load solve failed')
    half_error=float(np.linalg.norm(answer.array-.5*peak)/np.linalg.norm(.5*peak));rhs.scale(0);answer.set(0);ksp.solve(rhs,answer);zero_norm=float(answer.norm())
    if half_error>1e-6 or zero_norm>1e-12:raise RuntimeError('Half/zero-load verification failed')
    write(out/'linearity.json',dict(half_load_U_relative_difference=half_error,zero_load_U_norm=zero_norm,state_generation='Exact linear superposition from peak, verified by independent half/zero RHS solves; no nonlinear/hysteretic unloading claim'))
    # Rebuild only the output material-frame quadrature at the frozen 14 points.
    # Assembly frames live at degree6/8; their quadrature coefficients cannot be
    # evaluated as if they were an interpolating field at a different point set.
    del bilinear,linear,Q,T,FrameParameters;gc.collect()
    qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,4);nq=len(qp);assert nq==14;tab=scalar.tabulate(1,qp)
    Qspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(3,3),degree=4,scheme='default'));Q=fem.Function(Qspace)
    field=fem.Function(V);field.x.array[:]=peak;local_expression=fem.Expression(C*local_strain(field),qp);status('extracting_raw_GP',peak_solve_seconds=solve_seconds,true_relative_residual=true_residual)
    if args.save_raw_gp:
        raw_local=np.lib.format.open_memmap(out/'GP_stress_local_MPa.npy',mode='w+',dtype='float64',shape=(nc,nq,6));raw_global=np.lib.format.open_memmap(out/'GP_stress_xyz_MPa.npy',mode='w+',dtype='float64',shape=(nc,nq,6));raw_xyz=np.lib.format.open_memmap(out/'GP_xyz_mm.npy',mode='w+',dtype='float64',shape=(nc,nq,3));raw_weights=np.lib.format.open_memmap(out/'GP_volume_weights_mm3.npy',mode='w+',dtype='float64',shape=(nc,nq))
    else:
        raw_local=np.lib.format.open_memmap(out/'GP_stress_local_MPa.npy',mode='w+',dtype='float64',shape=(nc,nq,6)) if args.save_gp_local else None
        raw_global=raw_xyz=raw_weights=None
    cell_mean=np.zeros((nc,6));cell_local_mean=np.zeros((nc,6));cell_vm_mean=np.zeros(nc);stats={pid:dict(material=parts[pid]['material'],name=parts[pid]['name'],max_GP_von_Mises_MPa=0.,max_abs_fiber_sigma11_MPa=0.,volume_mm3=0.) for pid in parts}
    pairs=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)];minJ=float('inf')
    for first in range(0,nc,4000):
        ids=np.arange(first,min(first+4000,nc),dtype=np.int32);xyzq,jac=chunk_geometry(ids);minJ=min(minJ,float(jac.min()));axes=material_axes(xyzq,ids)
        for local,e in enumerate(ids):Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(int(e))]=axes[local].reshape(-1,9)
        stress=local_expression.eval(domain,ids).reshape(len(ids),nq,6);tensor=np.zeros((len(ids),nq,3,3))
        for j,(a0,a1) in enumerate(pairs):tensor[...,a0,a1]=stress[...,j];tensor[...,a1,a0]=stress[...,j]
        transformed=axes@tensor@np.swapaxes(axes,-1,-2);globalstress=np.stack([transformed[...,i,j] for i,j in pairs],axis=-1);w=jac*qw[None,:];orig=cell_order[ids]
        if raw_local is not None:raw_local[orig]=stress
        if args.save_raw_gp:raw_global[orig]=globalstress;raw_xyz[orig]=xyzq;raw_weights[orig]=w
        cell_mean[orig]=np.sum(globalstress*w[...,None],axis=1)/w.sum(1)[:,None];cell_local_mean[orig]=np.sum(stress*w[...,None],axis=1)/w.sum(1)[:,None];s=globalstress;vm=np.sqrt(.5*((s[...,0]-s[...,1])**2+(s[...,1]-s[...,2])**2+(s[...,2]-s[...,0])**2)+3*np.sum(s[...,3:]**2,axis=-1));cell_vm_mean[orig]=np.sum(vm*w,axis=1)/w.sum(1)
        for pid in np.unique(part_local[ids]):
            mask=part_local[ids]==pid;st=stats[int(pid)];st['max_GP_von_Mises_MPa']=max(st['max_GP_von_Mises_MPa'],float(vm[mask].max()));st['max_abs_fiber_sigma11_MPa']=max(st['max_abs_fiber_sigma11_MPa'],float(abs(stress[mask,...,0]).max()));st['volume_mm3']+=float(w[mask].sum())
    for mm in [raw_local,raw_global,raw_xyz,raw_weights]:
        if mm is not None:mm.flush()
    if minJ<=0:raise RuntimeError('Nonpositive integration Jacobian')
    np.savez_compressed(out/'cell_display_fields.npz',mean_stress_xyz_MPa=cell_mean,mean_stress_local_MPa=cell_local_mean,mean_GP_von_Mises_MPa=cell_vm_mean,cell_part_ids=pids,cell_ids=source['cell_ids'])
    maxU=float(np.linalg.norm(original,axis=1).max());states=[dict(state=i+1,pressure_MPa=p,max_U_mm=p*maxU,strain_energy_Nmm=p*p*energy,derivation='exact linear superposition') for i,p in enumerate(contract['load']['path_pressure_MPa'])];write(out/'load_unload_states.json',states)
    report=dict(status='completed',case_id=contract['case_id'],backend=backend_description,preconditioner=args.preconditioner,dolfinx=dolfinx.__version__,petsc=list(PETSc.Sys.getVersion()),matrix_type=gpu.getType(),vector_type=answer.getType(),geometry_mapping=review['geometry_mapping'],nodes=len(coords),geometry_nodes=len(xyz),cells=nc,dofs=3*len(coords),displacement_degree=args.degree,cap_convention='quadratic-opening-edge',geometry_parameters_mm=params,research_pilot=True,dataset_release=contract.get('dataset_release','quarantined_not_for_production_training'),research_validation_only=args.research_validation_only,native_schema_saved=bool(not args.research_validation_only or args.validate_native_schema),production_eligible=not args.research_validation_only,raw_GP_saved=args.save_raw_gp,local_GP_saved=bool(args.save_raw_gp or args.save_gp_local),structural_parts=12,display_only_seals=6,pressure_MPa=1.,true_relative_residual=true_residual,peak_iterations=peak_iterations,peak_reason=peak_reason,peak_solve_seconds=solve_seconds,max_U_mm=maxU,strain_energy_Nmm=energy,external_work_Nmm=work,energy_work_relative_error=abs(2*energy/work-1),force_balance_N=force_balance.tolist(),moment_balance_Nmm=moment_balance.tolist(),load_resultant_N=totalforce.tolist(),gauge_reactions_N=reactions.tolist(),cap_forces_N=forces,cap_areas_mm2=cap_areas,min_GP_Jacobian_mm3=minJ,material_frame_det_max_error=frame_error,half_load_relative_error=half_error,zero_load_U_norm=zero_norm,part_statistics=stats,mesh_sha256=sha(args.mesh),contract_sha256=sha(args.contract),solver_sha256=sha(Path(__file__)),total_seconds=time.perf_counter()-start,limitations=['Linear elastic ideal bonds; no seal compression/contact/friction/damage/yield',review['geometry_mapping']+'; repair candidate only, new same-geometry sensitivity not yet certified','von Mises for composites is a diagnostic, not a failure criterion','20 displayed states use exact linear superposition, validated at half and zero load; not 20 nonlinear integrations'])
    report['assembly_frame_storage']='compact_exact_axisymmetric_frame';report['assembly_quadrature_degree']=args.assembly_quadrature;report['output_quadrature_degree']=4;report['repair_preflight']=repair_preflight;report['candidate_verification_only']=True;report['original_shape_target_unresolved']=not repair_preflight['original_volume_SICN_target_passed'];report['production_eligible']=False
    write(out/'summary.json',report);status('completed',total_seconds=report['total_seconds']);print(json.dumps(report),flush=True)

if __name__=='__main__':main()
