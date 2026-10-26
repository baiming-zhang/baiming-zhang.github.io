"""Store the same axisymmetric material frames in two quadrature components."""
from pathlib import Path
import runpy
W=Path(__file__).resolve().parents[1];runpy.run_path(str(W/'notes/prepare_c04_curved_solver_v2.py'));p=W/'publication/GeoBC-Neural-Operator/tools/copv_fenicsx_c04_repair.py';s=p.read_text(encoding='utf-8')
helper="""
def compact_material_frame(domain,tangent_field,frame_parameters):
    x=ufl.SpatialCoordinate(domain);r=ufl.sqrt(x[0]*x[0]+x[1]*x[1]);er=ufl.as_vector([x[0]/r,x[1]/r,0]);et=ufl.as_vector([-x[1]/r,x[0]/r,0]);m=-tangent_field[0]*er-ufl.as_vector([0,0,tangent_field[1]]);n=ufl.cross(m,et);mask,ca,sa=frame_parameters[0],frame_parameters[1],frame_parameters[2];a=ca*m+sa*et;b=-sa*m+ca*et;rotation=ufl.as_matrix([[a[i],b[i],n[i]] for i in range(3)]);return mask*rotation+(1-mask)*ufl.Identity(3)

"""
s=s.replace('\ndef main():',helper+'\ndef main():',1)
old="Qspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(3,3),degree=args.assembly_quadrature,scheme='default'));Q=fem.Function(Qspace);Cspace="
new="Tspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(2,),degree=args.assembly_quadrature,scheme='default'));T=fem.Function(Tspace);FrameSpace=fem.functionspace(domain,basix.ufl.element('DG','tetrahedron',0,shape=(3,)));FrameParameters=fem.Function(FrameSpace);Q=compact_material_frame(domain,T,FrameParameters);Cspace="
assert old in s;s=s.replace(old,new,1)
old="        for j,e in enumerate(ids):Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(e)]=axes[j].reshape(-1,9);C.x.array.reshape(-1,36)[Cspace.dofmap.cell_dofs(e)]=p_by_id[int(part_local[e])].reshape(1,36)"
new="""        radii=np.linalg.norm(positions[...,:2],axis=-1);_,nearest=tree.query(np.stack([radii,positions[...,2]],axis=-1));tangent_values=tangent[nearest]
        for j,e in enumerate(ids):
            pid=int(part_local[e]);angle=np.deg2rad(angle_by_id[pid]);T.x.array.reshape(-1,2)[Tspace.dofmap.cell_dofs(e)]=tangent_values[j];FrameParameters.x.array.reshape(-1,3)[FrameSpace.dofmap.cell_dofs(e)]=np.array([[float(pid in carbon),np.cos(angle),np.sin(angle)]]);C.x.array.reshape(-1,36)[Cspace.dofmap.cell_dofs(e)]=p_by_id[pid].reshape(1,36)"""
assert old in s;s=s.replace(old,new,1)
old="resources('before_assembly',host_required=estimated_sparse_bytes+3*1024**3)"
new="""# Sparse storage is already preallocated. Reserve the actual coefficient
    # pack plus assembly scratch, rather than another complete sparse matrix.
    estimated_coefficient_pack_bytes=nc*(39+2*nq)*8
    resources('before_assembly',host_required=estimated_coefficient_pack_bytes+2*1024**3)
    write(out/'material_frame_storage.json',dict(method='Two quadrature tangent components and three DG0 mask/cos/sin parameters; spatial radial basis reconstructed exactly in UFL',assembly_points=nq,full_frame_value_bytes=nc*nq*9*8,compact_frame_value_bytes=nc*(nq*2+3)*8,estimated_coefficient_pack_bytes=estimated_coefficient_pack_bytes,preallocated_sparse_bytes=estimated_sparse_bytes,scope='Same axisymmetric material-frame values at the same quadrature points, not an interpolated or changed material model.'))"""
assert old in s;s=s.replace(old,new,1)
s=s.replace('del bilinear,linear,Q;gc.collect()','del bilinear,linear,Q,T,FrameParameters;gc.collect()',1)
s=s.replace("report['assembly_quadrature_degree']=args.assembly_quadrature;", "report['assembly_frame_storage']='compact_exact_axisymmetric_frame';report['assembly_quadrature_degree']=args.assembly_quadrature;",1)
compile(s,str(p),'exec');p.write_text(s,encoding='utf-8');print('Prepared exact compact-frame candidate; CPU equivalence validation required before GPU')
