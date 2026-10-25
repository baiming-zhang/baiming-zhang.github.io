"""Geometry-only periodic structured Hex27 COPV pilot, not a simulation dataset.

Native Gmsh reference coordinates determine all element connectivity. Physical
nodes are evaluated on the analytic normal-offset surfaces (not linear midpoints).
Current tetra10-only rod solver cannot consume these cells without adaptation.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import time
from geobc.copv_geometry import DEFAULT_PARAMETERS

SURFACES={1:'wetted_inner',2:'external_outer',3:'bonded_interface',
          4:'left_neck_liner',5:'right_neck_liner',6:'left_neck_overwrap',7:'right_neck_overwrap'}


def analytic_segments(parameters, counts):
    import numpy as np
    R=parameters['inner_radius_mm']; L=parameters['straight_length_mm']; D=parameters['dome_depth_mm']
    b=parameters['boss_bore_radius_mm']; neck=parameters['neck_length_mm']
    ctrl=np.array([[R,L/2],[R,L/2+.4*D],[b,L/2+.6*D],[b,L/2+D]])
    def dome(n):
        t=np.linspace(0,1,2*n+1)[:,None]
        rz=(1-t)**3*ctrl[0]+3*(1-t)**2*t*ctrl[1]+3*(1-t)*t*t*ctrl[2]+t**3*ctrl[3]
        dr=3*((1-t)**2*(ctrl[1]-ctrl[0])+2*t*(1-t)*(ctrl[2]-ctrl[1])+t*t*(ctrl[3]-ctrl[2]))
        normal=np.stack([dr[:,1],-dr[:,0]],axis=1);normal/=np.linalg.norm(normal,axis=1)[:,None]
        return rz,normal
    dl,nl=dome(counts[1]); dr,nr=dome(counts[3])
    pieces=[np.column_stack([np.full(2*counts[0]+1,b),np.linspace(-L/2-D-neck,-L/2-D,2*counts[0]+1)]),
            dl[::-1]*[1,-1],np.column_stack([np.full(2*counts[2]+1,R),np.linspace(-L/2,L/2,2*counts[2]+1)]),
            dr,np.column_stack([np.full(2*counts[4]+1,b),np.linspace(L/2+D,L/2+D+neck,2*counts[4]+1)])]
    normals=[np.tile([1.,0.],(len(pieces[0]),1)),nl[::-1]*[1,-1],np.tile([1.,0.],(len(pieces[2]),1)),nr,np.tile([1.,0.],(len(pieces[4]),1))]
    join=lambda a:np.concatenate([v if i==0 else v[1:] for i,v in enumerate(a)])
    return join(pieces),join(normals)


def generate_copv_structured_mesh(output_dir, nphi=24, segment_counts=(2,8,20,8,2),
                                  parameters=None, liner_elements=1, overwrap_elements=2):
    import numpy as np
    import gmsh
    start=time.perf_counter();p=dict(DEFAULT_PARAMETERS)
    if parameters:
        if set(parameters)-set(p):raise ValueError('Unknown COPV parameter')
        p.update(parameters)
    if nphi<12 or len(segment_counts)!=5 or min(segment_counts)<1 or min(liner_elements,overwrap_elements)<1:
        raise ValueError('Invalid pilot mesh counts')
    if any(not np.isfinite(v) or v<=0 for v in p.values()):raise ValueError('Positive finite geometry required')
    if p['boss_bore_radius_mm']>=p['inner_radius_mm'] or p['liner_thickness_mm']+p['overwrap_thickness_mm']>=.1*p['inner_radius_mm']:
        raise ValueError('Invalid pilot geometry wall')
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    rz,normal=analytic_segments(p,segment_counts)
    offsets=np.concatenate([np.linspace(0,p['liner_thickness_mm'],2*liner_elements+1),
                            p['liner_thickness_mm']+np.linspace(0,p['overwrap_thickness_mm'],2*overwrap_elements+1)[1:]])
    nphi2=2*nphi; ns=sum(segment_counts); nt=liner_elements+overwrap_elements; kinterface=2*liner_elements
    shape=(nphi2,2*ns+1,2*nt+1)
    ijk=np.indices(shape).reshape(3,-1).T
    phi=ijk[:,0]*(2*np.pi/nphi2); s=ijk[:,1]; k=ijk[:,2]
    rr=rz[s,0]+offsets[k]*normal[s,0]; zz=rz[s,1]+offsets[k]*normal[s,1]
    xyz=np.column_stack([rr*np.cos(phi),rr*np.sin(phi),zz])
    node_ids=np.arange(1,len(xyz)+1,dtype=np.int64)
    grid=node_ids.reshape(shape)
    gmsh.initialize()
    try:
        gmsh.model.add('copv_periodic_hex27_geometry_pilot')
        for entity in SURFACES:gmsh.model.addDiscreteEntity(2,entity)
        gmsh.model.addDiscreteEntity(3,10,[1,3,4,5]);gmsh.model.addDiscreteEntity(3,11,[2,3,6,7])
        # Geometric classification of nodes; interface is one shared node set.
        classification=np.where(k<kinterface,10,11)
        classification=np.where(s==0,np.where(k<kinterface,4,6),classification)
        classification=np.where(s==2*ns,np.where(k<kinterface,5,7),classification)
        classification=np.where(k==0,1,classification);classification=np.where(k==2*nt,2,classification)
        classification=np.where(k==kinterface,3,classification)
        for entity in np.unique(classification):
            mask=classification==entity
            gmsh.model.mesh.addNodes(3 if entity>=10 else 2,int(entity),node_ids[mask],xyz[mask].ravel())
        hlocal=np.asarray(gmsh.model.mesh.getElementProperties(12)[4]).reshape(27,3)
        qlocal=np.asarray(gmsh.model.mesh.getElementProperties(10)[4]).reshape(9,2)
        hoff=np.rint(hlocal+1).astype(int);qoff=np.rint(qlocal+1).astype(int)
        elements=np.indices((nphi,ns,nt)).reshape(3,-1).T
        inds=2*elements[:,None,:]+hoff[None,:,:]
        cells=grid[inds[:,:,0]%nphi2,inds[:,:,1],inds[:,:,2]]-1
        cell_ids=np.arange(1,len(cells)+1,dtype=np.int64)
        material=np.where(elements[:,2]<liner_elements,10,11)
        for entity,name in [(10,'liner_volume'),(11,'overwrap_volume')]:
            mask=material==entity
            gmsh.model.mesh.addElementsByType(entity,12,cell_ids[mask],(cells[mask]+1).ravel())
            gmsh.model.addPhysicalGroup(3,[entity],entity);gmsh.model.setPhysicalName(3,entity,name)
        faces=[];face_group=[]
        def add_faces(group,axis,fixed,range1,range2,reverse=False):
            for u in range1:
                for v in range2:
                    refs=np.column_stack([np.full(9,fixed),2*u+qoff[:,0],2*v+qoff[:,1]])
                    if axis==2:refs=refs[:,[1,2,0]] # (phi,s) at fixed offset
                    elif axis==1:refs=refs[:,[1,0,2]] # (phi,offset) at fixed s
                    # Native quad orientation is documented in manifest; solver must use outward solid normals.
                    if reverse:refs=refs[[0,3,2,1,7,6,5,4,8]]
                    f=grid[refs[:,0]%nphi2,refs[:,1],refs[:,2]]-1
                    faces.append(f);face_group.append(group)
        add_faces(1,2,0,range(nphi),range(ns),True)
        add_faces(2,2,2*nt,range(nphi),range(ns))
        add_faces(3,2,kinterface,range(nphi),range(ns))
        for boundary in [0,2*ns]:
            add_faces(4 if boundary==0 else 5,1,boundary,range(nphi),range(liner_elements),boundary!=0)
            add_faces(6 if boundary==0 else 7,1,boundary,range(nphi),range(liner_elements,nt),boundary!=0)
        faces=np.array(faces,dtype=np.int64);face_group=np.array(face_group,dtype=np.int32)
        face_ids=np.arange(len(cells)+1,len(cells)+1+len(faces),dtype=np.int64)
        for group,name in SURFACES.items():
            mask=face_group==group;gmsh.model.mesh.addElementsByType(group,10,face_ids[mask],(faces[mask]+1).ravel())
            gmsh.model.addPhysicalGroup(2,[group],group);gmsh.model.setPhysicalName(2,group,name)
        local,weights=gmsh.model.mesh.getIntegrationPoints(12,'Gauss4')
        _,det,_=gmsh.model.mesh.getJacobians(12,local)
        det=np.asarray(det)
        quality=np.asarray(gmsh.model.mesh.getElementQualities(cell_ids,'minSICN'))
        min_det=np.asarray(gmsh.model.mesh.getElementQualities(cell_ids,'minDetJac'))
        if not np.isfinite(det).all() or min(det.min(),quality.min(),min_det.min())<=0:
            raise RuntimeError(f'Nonpositive Jacobian/quality: {det.min()}, {quality.min()}, {min_det.min()}')
        interface_nodes=np.unique(faces[face_group==3])
        liner_nodes=np.unique(cells[material==10]);wrap_nodes=np.unique(cells[material==11])
        if not np.array_equal(np.intersect1d(liner_nodes,wrap_nodes),interface_nodes):
            raise RuntimeError('Interface is not exactly shared by both materials')
        if len(np.unique(cells))!=len(xyz):raise RuntimeError('Unused nodes')
        gmsh.write(str(out/'mesh.msh'))
        np.savez_compressed(out/'mesh.npz',xyz=xyz,cells=cells,faces=faces,node_ids=node_ids,cell_ids=cell_ids,
                            cell_material_ids=material,face_ids=face_ids,face_group_ids=face_group,
                            gmsh_hex27_reference_nodes=hlocal,gmsh_quad9_reference_nodes=qlocal,
                            structured_node_ijk=ijk,meridian_inner_rz=rz,meridian_normal_rz=normal,
                            thickness_offsets_mm=offsets,min_sicn=quality,min_det_jac=min_det)
        report=dict(schema_version=1,status='geometry_mesh_pilot_only',simulation_labels_generated=False,
                    production_ready=False,parameters_mm=p,nphi=nphi,meridian_segment_elements=list(segment_counts),
                    liner_elements=liner_elements,overwrap_elements=overwrap_elements,
                    node_count=len(xyz),cell_count=len(cells),quad9_count=len(faces),gmsh_version=gmsh.__version__,
                    gmsh_cell_type=12,gmsh_face_type=10,reference_axes=['phi','increasing-z meridian','outward normal offset'],
                    min_sampled_jacobian_mm3=float(det.min()),min_det_jac_mm3=float(min_det.min()),min_sicn=float(quality.min()),
                    shared_interface_node_count=len(interface_nodes),periodic_phi_seam='modulo indexing; no duplicated end nodes',
                    all_nodes_used=True,interface_exact_shared_nodes=True,boundary_groups=SURFACES,
                    layer_normal_offsets_mm=offsets.tolist(),generation_seconds=time.perf_counter()-start,
                    limitations=['Not compatible with current tetra10 rod FE solver',
                                 'Future MeshHex2 adapter must map reference-node coordinates, not assume identical ordering',
                                 'Curved physical nodes lie on analytic normal offsets; cell geometry is triquadratic interpolation',
                                 'No mesh convergence or FE solve yet; double-open neck requires reviewed pressure end resultants',
                                 'No layup/material/pressure freeze; no winding manufacturability claim'])
        report['artifacts']={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['mesh.msh','mesh.npz']}
        (out/'mesh_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
        return report
    finally:gmsh.finalize()


def preview_structured_mesh(output_dir):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    out=Path(output_dir);a=np.load(out/'mesh.npz');x=a['xyz'];f=a['faces'];g=a['face_group_ids']
    fig=plt.figure(figsize=(10,6));ax=fig.add_subplot(111,projection='3d',proj_type='ortho')
    mask=g!=3
    color=np.array(['#76a398' if v==2 else '#dba85b' if v in [4,5] else '#9cb4c6' for v in g])[mask]
    ax.add_collection3d(Poly3DCollection(x[f[mask,:4]],facecolors=color,edgecolors='#36564d',linewidths=.22,alpha=.92))
    lo=x.min(0);hi=x.max(0);ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_zlim(lo[2],hi[2]);ax.set_box_aspect(hi-lo)
    ax.view_init(20,-65);ax.set_axis_off();fig.suptitle(f'COPV periodic Hex27 mesh | {len(x):,} nodes | geometry only')
    fig.text(.5,.05,'Curved midside/face/body nodes retained in mesh; preview draws corner-edge outlines.\nShared liner / overwrap interface; open neck ends; no simulation labels.',ha='center',fontsize=9)
    fig.savefig(out/'mesh_preview.png',dpi=180,bbox_inches='tight');plt.close(fig)
