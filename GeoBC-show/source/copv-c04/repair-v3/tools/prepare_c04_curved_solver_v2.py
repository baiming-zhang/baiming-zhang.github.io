"""Add explicit assembly quadrature while preserving the original 14 output GPs."""
from pathlib import Path
import runpy
W=Path(__file__).resolve().parents[1];runpy.run_path(str(W/'notes/prepare_c04_curved_solver.py'));p=W/'publication/GeoBC-Neural-Operator/tools/copv_fenicsx_c04_repair.py';s=p.read_text(encoding='utf-8')
s=s.replace('import argparse,gc,hashlib,json,os,sysconfig,time','import argparse,gc,hashlib,json,os,sysconfig,time,subprocess')
s=s.replace("ap.add_argument('--repair-preflight'","ap.add_argument('--assembly-quadrature',type=int,choices=[4,6,8],default=6);ap.add_argument('--repair-preflight'",1)
anchor="    native_path=args.repair_preflight.parent/'native_boundary_vectors.npz'"
addition="""    candidate=args.repair_preflight.parent.parent
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
"""
assert anchor in s;s=s.replace(anchor,addition+anchor,1)
s=s.replace("qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,4)","qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,args.assembly_quadrature)",1)
s=s.replace("value_shape=(3,3),degree=4,scheme='default'","value_shape=(3,3),degree=args.assembly_quadrature,scheme='default'",1)
s=s.replace("metadata={'quadrature_degree':4});ds=", "metadata={'quadrature_degree':args.assembly_quadrature});ds=",1)
begin=s.index('    for first in range(0,nc,4000):');end=s.index('    if frame_error>1e-10:',begin)
old=s[begin:end];axes_begin=old.index('axes=np.broadcast_to');axes_end=old.index('        frame_error=max')
axes=old[axes_begin:axes_end]
function='    def material_axes(positions,ids):\n        '+axes+'        return axes\n'
function=function.replace('nq,3,3','positions.shape[1],3,3')
newloop="""    for first in range(0,nc,4000):
        ids=np.arange(first,min(first+4000,nc),dtype=np.int32);positions,_=chunk_geometry(ids);axes=material_axes(positions,ids)
        frame_error=max(frame_error,float(abs(np.linalg.det(axes)-1).max()))
        for j,e in enumerate(ids):Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(e)]=axes[j].reshape(-1,9);C.x.array.reshape(-1,36)[Cspace.dofmap.cell_dofs(e)]=p_by_id[int(part_local[e])].reshape(1,36)
"""
s=s[:begin]+function+newloop+s[end:]
anchor="    A=fp.create_matrix(bilinear);status('matrix_preallocated',matrix_info=A.getInfo());fp.assemble_matrix(A,bilinear);"
replace="""    A=fp.create_matrix(bilinear);preallocated=A.getInfo();status('matrix_preallocated',matrix_info=preallocated)
    estimated_sparse_bytes=int(preallocated['nz_allocated'])*(8+PETSc.IntType().itemsize)
    resources('before_assembly',host_required=estimated_sparse_bytes+3*1024**3)
    fp.assemble_matrix(A,bilinear);"""
assert anchor in s;s=s.replace(anchor,replace,1)
# Check actual CPU CSR storage before any CUDA matrices/vectors are created.
anchor="    corner_ids=np.unique(cells[:,:4]);"
replace="""    csr_i,csr_j,csr_v=A.getValuesCSR();actual_sparse_bytes=csr_i.nbytes+csr_j.nbytes+csr_v.nbytes;del csr_i,csr_j,csr_v
    resources('before_CUDA',host_required=2*1024**3,gpu_required=max(24576,int((2.0*actual_sparse_bytes+3*1024**3)/1024**2)))
    write(out/'quadrature_protocol.json',dict(assembly_degree=args.assembly_quadrature,output_degree=4,output_points=14,load_degree=4,actual_CSR_bytes=actual_sparse_bytes,scope='Assembly sensitivity candidate; original GP/node IDs and 14-point field layout preserved.'))
    corner_ids=np.unique(cells[:,:4]);"""
assert anchor in s;s=s.replace(anchor,replace,1)
anchor="    field=fem.Function(V);field.x.array[:]=peak;local_expression="
replace="""    # Rebuild only the output material-frame quadrature at the frozen 14 points.
    # Assembly frames live at degree6/8; their quadrature coefficients cannot be
    # evaluated as if they were an interpolating field at a different point set.
    del bilinear,linear,Q;gc.collect()
    qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,4);nq=len(qp);assert nq==14;tab=scalar.tabulate(1,qp)
    Qspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(3,3),degree=4,scheme='default'));Q=fem.Function(Qspace)
    field=fem.Function(V);field.x.array[:]=peak;local_expression="""
assert anchor in s;s=s.replace(anchor,replace,1)
anchor="stress=local_expression.eval(domain,ids).reshape(len(ids),nq,6);axes=np.stack([Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(int(e))].reshape(nq,3,3) for e in ids]);tensor="
replace="""axes=material_axes(xyzq,ids)
        for local,e in enumerate(ids):Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(int(e))]=axes[local].reshape(-1,9)
        stress=local_expression.eval(domain,ids).reshape(len(ids),nq,6);tensor="""
assert anchor in s;s=s.replace(anchor,replace,1)
s=s.replace("report['repair_preflight']=repair_preflight;", "report['assembly_quadrature_degree']=args.assembly_quadrature;report['output_quadrature_degree']=4;report['repair_preflight']=repair_preflight;",1)
compile(s,str(p),'exec');p.write_text(s,encoding='utf-8');print('Prepared guarded degree4/6/8 assembly with frozen 14-GP output; runtime validation pending')
