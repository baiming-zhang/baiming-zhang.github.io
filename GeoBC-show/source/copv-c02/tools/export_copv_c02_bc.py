"""Export the executed C02 boundary/load locations without changing the solve.

Facet classification reproduces the executed solver on the accepted affine
geometry. Masks and reactions come from saved actual constraint DOFs.
"""
from pathlib import Path
import argparse, base64, gzip, hashlib, json
import numpy as np
from geobc.copv_geometry import DEFAULT_PARAMETERS
from geobc.copv_structured_mesh import analytic_segments

def encode(a, dtype):
    return base64.b64encode(np.asarray(a, dtype=dtype).tobytes()).decode()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    root=args.root;out=args.output;out.mkdir(parents=True,exist_ok=True)
    report=json.loads((root/'fenicsx/summary.json').read_text(encoding='utf-8'))
    meshfile=root/'mesh_h6_netgen/mesh.npz'
    assert hashlib.sha256(meshfile.read_bytes()).hexdigest()==report['mesh_sha256']
    m=np.load(meshfile);s=np.load(root/'surface/surface.npz');a=np.load(root/'fenicsx/nodal_fields.npz')
    x=m['xyz'];u=a['U_mm'];cells=m['cells'];faces=s['faces'][s['is_exterior']];owners=s['owner_cells'][s['is_exterior']]
    center=x[faces[:,:3]].mean(1);area=.5*np.cross(x[faces[:,1]]-x[faces[:,0]],x[faces[:,2]]-x[faces[:,0]])
    opposite=np.array([next(int(v) for v in cells[e,:4] if v not in set(f[:3])) for e,f in zip(owners,faces)])
    flip=np.einsum('ij,ij->i',area,x[opposite]-center)>0;area[flip]*=-1
    magnitude=np.linalg.norm(area,axis=1);normals=area/magnitude[:,None]
    rz,_=analytic_segments(DEFAULT_PARAMETERS,(4,2500,16,2500,4));z=center[:,2];r=np.linalg.norm(center[:,:2],axis=1)
    inner=np.interp(abs(z),rz[rz[:,1]>=0,1],rz[rz[:,1]>=0,0]);inner=np.where(z>277,10.2,np.where(z<-277,8.25,inner))
    tags=np.full(len(faces),4);tags[r<=inner+.30]=1;tags[np.isclose(z,302,atol=1e-5)]=2;tags[np.isclose(z,-295,atol=1e-5)]=3
    old=json.loads((root/'surface/load_closure_review.json').read_text(encoding='utf-8'))
    counts={str(i):int((tags==i).sum()) for i in [1,2,3]}
    assert [counts[str(i)] for i in [1,2,3]]==[old['pressure_boundary_facets'],old['top_cap_facets'],old['bottom_cap_facets']]
    constraints=[];fixed=a['gauge_dofs_FEniCSx'];reaction=a['gauge_reactions_N'];starts=[(0,3),(3,5),(5,6)]
    for k,(start,stop) in enumerate(starts):
        node=int(a['gauge_original_nodes'][k]);assert np.allclose(x[node],a['gauge_coordinates_mm'][k],atol=1e-10)
        axes=(fixed[start:stop]%3).astype(int).tolist();R=np.zeros(3)
        for axis,value in zip(axes,reaction[start:stop]):R[axis]=value
        constraints.append(dict(id='ABC'[k],type='nodal_translation_gauge',original_node_index=node,gmsh_node_id=int(m['node_ids'][node]),position_mm=x[node].tolist(),U_peak_mm=u[node].tolist(),locked_axes=axes,reaction_peak_N=R.tolist(),values_mm=[0.]*len(axes)))
    assert [v['locked_axes'] for v in constraints]==[[0,1,2],[1,2],[2]]
    groups=[];caps=[]
    for tag,name,color in [(1,'pressure',[.96,.34,.08,.19]),(2,'upper_cap',[.85,.12,.32,.65]),(3,'lower_cap',[.85,.12,.32,.65])]:
        select=np.flatnonzero(tags==tag);node_ids,inverse=np.unique(faces[select,:3],return_inverse=True)
        groups.append(dict(name=name,kind=5,xyz=encode(x[node_ids],'<f4'),u=encode(u[node_ids],'<f4'),tri=encode(inverse.reshape(-1,3),'<u4'),color=color))
        if tag>1:
            weights=magnitude[select];pos=np.average(center[select],axis=0,weights=weights);disp=np.average(u[faces[select,:3]].mean(1),axis=0,weights=weights)
            key='upper' if tag==2 else 'lower';direction=[0.,0.,1. if tag==2 else -1.]
            caps.append(dict(id=key,position_mm=pos.tolist(),U_peak_mm=disp.tolist(),direction=direction,force_peak_N=report['cap_forces_N'][key],annular_area_mm2=float(weights.sum()),facet_count=len(select)))
    # Direction samples are actual pressure facets, not an invented wall/cap.
    candidates=np.flatnonzero(tags==1);chosen=set()
    for targetz in [-240.,-160.,-80.,0.,80.,160.,240.]:
        near=candidates[abs(z[candidates]-targetz)<12]
        if not len(near):near=candidates
        for phi in [0.,np.pi/2,np.pi,3*np.pi/2]:
            radius=float(np.interp(abs(targetz),rz[rz[:,1]>=0,1],rz[rz[:,1]>=0,0]));target=np.array([radius*np.cos(phi),radius*np.sin(phi),targetz]);chosen.add(int(near[np.argmin(np.linalg.norm(center[near]-target,axis=1))]))
    arrows=[]
    for f in sorted(chosen):arrows.append(dict(position_mm=center[f].tolist(),U_peak_mm=u[faces[f,:3]].mean(0).tolist(),direction=(-normals[f]).tolist(),length_mm=min(24.,float(r[f]*.65))))
    resultant=-area[tags==1].sum(0)+np.array([0.,0.,report['cap_forces_N']['upper']-report['cap_forces_N']['lower']])
    assert np.linalg.norm(resultant-np.array(old['pressure_plus_caps_force_resultant_N']))<1e-8
    metadata=dict(schema='copv_bc_load_visualization_v1',mesh_sha256=report['mesh_sha256'],contract_sha256=report['contract_sha256'],peak_pressure_MPa=report['pressure_MPa'],constraints=constraints,pressure=dict(type='surface_traction',formula='t=-p*n_solid',facet_count=counts['1'],arrows=arrows,external_pressure_MPa=0.),caps=caps,wall_supports=[],hinges=[],connections='Perfectly bonded shared structural interfaces',pressure_plus_caps_resultant_peak_N=resultant.tolist(),gauge_reaction_sum_peak_N=np.sum([v['reaction_peak_N'] for v in constraints],axis=0).tolist(),notes=['Constraint masks are the actual saved six translational DOFs; no whole end face is clamped.','Pressure highlight contains every selected wetted facet, including port interiors/shoulders; arrows are direction samples.','Caps are equivalent resultants distributed on the actual annular ends; no cap solid or hinge is added.','Symbols use reference-load directions and follow saved displacement for display; arrow lengths are schematic, not force magnitudes.','The small axial resultant is the recorded circular-cap/polygon-bore discretization remainder, carried by the numerical gauge.'])
    (out/'bc-loads.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    payload=json.dumps(dict(parts=groups),separators=(',',':')).encode();(out/'bc-load-surfaces.json.gz').write_bytes(gzip.compress(payload,mtime=0))
    print(json.dumps(dict(facet_counts=counts,locked_components=6,pressure_arrows=len(arrows),resultant_N=resultant.tolist(),surface_bytes=(out/'bc-load-surfaces.json.gz').stat().st_size)))

if __name__=='__main__':main()
