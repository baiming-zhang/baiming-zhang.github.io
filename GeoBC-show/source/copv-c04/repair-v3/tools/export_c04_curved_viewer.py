"""Export completed curved native fields and actual Tri6 BC/loads to PRIVATE staging."""
from pathlib import Path
import argparse,base64,gzip,hashlib,json,os,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np
from geobc.copv_quadratic_boundary import native_boundary_loads,triangle_rule
from geobc.copv_structured_mesh import analytic_segments
W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def enc(a,dtype):return base64.b64encode(np.asarray(a,dtype=dtype).tobytes()).decode()
def save(p,v):p.write_bytes(gzip.compress(json.dumps(v,ensure_ascii=False,separators=(',',':')).encode(),compresslevel=3,mtime=0))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def quant(a):
 scale=max(float(np.max(abs(a)))/32760,1e-30);q=np.rint(a/scale).astype('<i2');return enc(q,'<i2'),scale,float(np.max(abs(q.astype(float)*scale-a)))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();c=a.candidate.resolve();run=a.run.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);r=read(run/'summary.json');assert r['status']=='completed' and r['candidate_verification_only'] and not r['production_eligible'];assert sha(c/'mesh/mesh.npz')==r['mesh_sha256'];m=np.load(c/'mesh/mesh.npz');n=np.load(run/'nodal_fields.npz');v=np.load(run/'native_load_and_bc.npz');s=np.load(run/'cell_display_fields.npz');sf=np.load(c/'preflight/surface.npz');review=read(c/'mesh/mesh_review_curved_checked.json');assert np.array_equal(n['node_ids'],m['node_ids']) and np.array_equal(v['node_ids'],m['node_ids']) and np.array_equal(s['cell_ids'],m['cell_ids']) and np.array_equal(s['cell_part_ids'],m['cell_part_ids']);x=m['xyz'];u=n['U_mm'];cells=m['cells'];pids=m['cell_part_ids'];fiber=s['mean_stress_local_MPa'][:,0];vm=s['mean_GP_von_Mises_MPa'];parts={p['part_id']:p for p in review['parts'] if p['structural']};kinds={name:i for i,name in enumerate(['carbon','aluminum','titanium','steel','adhesive'])};body=[];display_errors={}
 assert read(c/'preflight/review.json')['original_mesh_screen_passed'],'This viewer explicitly labels original mesh screens as passed'
 for pid,p in parts.items():
  take=(sf['part_ids']==pid)&sf['is_exterior'];faces=sf['faces'][take];owners=sf['owner_cells'][take]
  if not len(faces):continue
  row=dict(part_id=pid,name=p['name'],material=p['material'],kind=kinds[p['material']],xyz=enc(x[faces].reshape(-1,3),'<f4'),tri6_faces=len(faces),quantized_int16_scales={});errors={}
  for key,arr in [('u',u[faces].reshape(-1,3)),('fiber',np.repeat(fiber[owners],6)),('vm',np.repeat(vm[owners],6))]:row[key],row['quantized_int16_scales'][key],errors[key]=quant(arr)
  body.append(row);display_errors[str(pid)]=errors
 save(out/'surface.json.gz',dict(parts=body));del body
 # A clipped view uses real boundary triangles. No affine cut interpolation
 # is used on the isoparametric geometry; section filling is explicitly absent.
 save(out/'section.json.gz',dict(parts=[]))
 faces=sf['faces'][sf['is_exterior']];owners=sf['owner_cells'][sf['is_exterior']];params=r['geometry_parameters_mm'];rz,_=analytic_segments(params,(4,2500,16,2500,4));loads,tags,closure=native_boundary_loads(x,cells,faces,owners,params,rz);assert closure['passed'];load_errors={}
 for arr,key in zip(loads,['pressure_load_N','upper_cap_load_N','lower_cap_load_N']):
  stored=v[key];err=float(np.linalg.norm(arr-stored)/max(np.linalg.norm(stored),1e-30));assert err<1e-10;load_errors[key]=err
 del loads
 assert v['locked_components'].sum()==6;constraints=[]
 for k,node in enumerate(n['gauge_original_nodes']):
  node=int(node);assert np.allclose(x[node],n['gauge_coordinates_mm'][k],atol=1e-10);axes=np.flatnonzero(v['locked_components'][node]).tolist();constraints.append(dict(id='ABC'[k],type='nodal_translation_gauge',original_node_index=node,gmsh_node_id=int(m['node_ids'][node]),position_mm=x[node].tolist(),U_peak_mm=u[node].tolist(),locked_axes=axes,reaction_peak_N=v['reaction_N'][node].tolist(),values_mm=[0.]*len(axes)))
 assert [p['locked_axes'] for p in constraints]==[[0,1,2],[1,2],[2]]
 basis,deriv,w=triangle_rule();center=(-x[faces[:,:3]].sum(1)+4*x[faces[:,3:]].sum(1))/9;motion=(-u[faces[:,:3]].sum(1)+4*u[faces[:,3:]].sum(1))/9;groups=[];caps=[];sub=np.array([[0,3,5],[3,1,4],[5,4,2],[3,4,5]])
 for tag,name,color in [(1,'pressure',[.96,.34,.08,.19]),(2,'upper_cap',[.85,.12,.32,.65]),(3,'lower_cap',[.85,.12,.32,.65])]:
  ids=np.flatnonzero(tags==tag);nodes,inverse=np.unique(faces[ids],return_inverse=True);groups.append(dict(name=name,kind=5,xyz=enc(x[nodes],'<f4'),u=enc(u[nodes],'<f4'),tri=enc(inverse.reshape(-1,6)[:,sub].reshape(-1,3),'<u4'),color=color))
  if tag>1:
   coords=x[faces[ids]];j=np.einsum('eni,qnj->eqij',coords,deriv);weights=np.linalg.norm(np.cross(j[:,:,:,0],j[:,:,:,1]),axis=2)*w;area=float(weights.sum());pos=np.sum(np.einsum('qn,fnd->fqd',basis,coords)*weights[:,:,None],axis=(0,1))/area;disp=np.sum(np.einsum('qn,fnd->fqd',basis,u[faces[ids]])*weights[:,:,None],axis=(0,1))/area;key='upper' if tag==2 else 'lower';assert np.isclose(area,r['cap_areas_mm2']['top' if tag==2 else 'bottom'],rtol=1e-10);caps.append(dict(id=key,position_mm=pos.tolist(),U_peak_mm=disp.tolist(),direction=[0.,0.,1. if tag==2 else -1.],force_peak_N=r['cap_forces_N'][key],annular_area_mm2=area,facet_count=len(ids)))
 candidates=np.flatnonzero(tags==1);selected=set();end=params['straight_length_mm']/2+params['dome_depth_mm'];positive=rz[rz[:,1]>=0]
 for z in np.linspace(-.92*end,.92*end,7):
  near=candidates[abs(center[candidates,2]-z)<12]
  if not len(near):near=candidates
  radius=np.interp(abs(z),positive[:,1],positive[:,0])
  for phi in [0,np.pi/2,np.pi,3*np.pi/2]:selected.add(int(near[np.argmin(np.linalg.norm(center[near]-[radius*np.cos(phi),radius*np.sin(phi),z],axis=1))]))
 grad=np.array([[-1.,-1.],[1.,0.],[0.,1.]]);lam=np.ones(3)/3;d=np.empty((6,2));d[:3]=(4*lam[:,None]-1)*grad
 for k,(i,j) in enumerate([(0,1),(1,2),(2,0)],3):d[k]=4*(lam[i]*grad[j]+lam[j]*grad[i])
 arrows=[]
 for f in sorted(selected):
  tangent=x[faces[f]].T@d;normal=np.cross(tangent[:,0],tangent[:,1]);opp=next(int(i) for i in cells[owners[f],:4] if i not in set(faces[f,:3]));normal*= -1 if normal@(x[opp]-center[f])>0 else 1;normal/=np.linalg.norm(normal);arrows.append(dict(position_mm=center[f].tolist(),U_peak_mm=motion[f].tolist(),direction=(-normal).tolist(),length_mm=min(24.,float(np.linalg.norm(center[f,:2])*.65))))
 metadata=dict(schema='copv_bc_load_visualization_v1',mesh_sha256=r['mesh_sha256'],contract_sha256=r['contract_sha256'],peak_pressure_MPa=r['pressure_MPa'],constraints=constraints,pressure=dict(type='surface_traction',formula='t=-p*n_solid',facet_count=int((tags==1).sum()),arrows=arrows,external_pressure_MPa=0.),caps=caps,wall_supports=[],hinges=[],connections='Perfectly bonded shared structural interfaces',pressure_plus_caps_resultant_peak_N=v['consistent_load_N'].sum(0).tolist(),gauge_reaction_sum_peak_N=np.sum([p['reaction_peak_N'] for p in constraints],axis=0).tolist(),notes=['Actual saved native six-DOF masks/reactions; no entire end clamp.','Pressure directions evaluated from actual Tri6 tangent Jacobians at selected facet centers.','Highlight and cap areas use actual curved Tri6 geometry; opening-edge forces and native vectors match executed solve.','Arrow length is schematic. Every glyph follows saved actual nodal/P2 facet motion.','Curved display uses four triangles per original Tri6; clip only, no affine section fill.'])
 (out/'bc-loads.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');save(out/'bc-load-surfaces.json.gz',dict(parts=groups));carbon=np.isin(pids,[i for i,p in parts.items() if p['material']=='carbon']);metal=np.isin(pids,[i for i,p in parts.items() if p['material'] in ['aluminum','titanium','steel']]);manifest=dict(schema='c04_curved_native_GPU_display_v1',parts=[dict(part_id=i,name=p['name'],material=p['material']) for i,p in parts.items()],ranges=dict(U_mm=r['max_U_mm'],fiber_MPa=float(abs(fiber[carbon]).max()),metal_VM_MPa=float(vm[metal].max())),states=read(run/'load_unload_states.json'),summary=r,display_method='Actual isoparametric Tri6 nodes, four display triangles each; original-cell GP volume mean stress, no cross-material averaging. Clipped boundaries only, no affine section filling.',display_quantization_max_abs_error_at_1MPa=display_errors,state_method='Exact validated linear-elastic scaling; no image/frame blending',native_load_recompute_relative_errors=load_errors,source_run=str(run),production_eligible=False,seconds=time.perf_counter()-start)
 (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');(out/'export_review.json').write_text(json.dumps(dict(status='exported_pending_browser_QA',mesh_sha256=r['mesh_sha256'],actual_native_loads_checked=True,locked_components=6,pressure_arrows=len(arrows),curved_basis=True,affine_section_fill=False,seconds=time.perf_counter()-start,files={p.name:p.stat().st_size for p in out.iterdir()}),indent=2)+'\n');print(json.dumps(dict(output=str(out),seconds=time.perf_counter()-start,locked_components=6,pressure_arrows=len(arrows))),flush=True)
if __name__=='__main__':main()
