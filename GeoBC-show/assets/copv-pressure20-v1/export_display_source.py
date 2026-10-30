from pathlib import Path
import numpy as np,json,gzip,hashlib,shutil
ROOT=Path(__file__).resolve().parent
WORK=ROOT.parent/'experiments/copv_axisym_simplified_v1'
SRC=WORK/'pressure_sweep_20_v1'; OUT=ROOT/'site/GeoBC-show'; DATA=OUT/'assets/copv-pressure20-v1'
DATA.mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def pack(path,arrays):
    raw=bytearray(); desc={}
    for name,a in arrays.items():
        a=np.ascontiguousarray(a);desc[name]={'offset':len(raw),'shape':list(a.shape),'dtype':a.dtype.str};raw.extend(a.tobytes())
    path.write_bytes(gzip.compress(bytes(raw),compresslevel=7,mtime=0))
    return {'file':path.name,'sha256':sha(path),'bytes':path.stat().st_size,'arrays':desc}
proof=read(SRC/'FINAL_VERIFICATION.json');assert proof['all_20_actual_V100_cases_verified']
mesh=np.load(SRC/'cloud/homogeneous_E_input.npz');xy=mesh['xy'];cells=mesh['cells'];ij=mesh['node_ij'];nrad=int(ij[:,0].max()+1);nsta=int(ij[:,1].max()+1)
grid=np.empty((nrad,nsta),dtype=np.uint32);grid[ij[:,0],ij[:,1]]=np.arange(len(ij))
assert nrad==97 and nsta==2881
bycell=np.empty((48,1440),dtype=np.uint32);cij=mesh['original_cell_ij'];bycell[cij[:,-2],cij[:,-1]]=np.arange(len(cells))
# Native quadratic boundary edges split at their actual midpoint, never decimated.
boundary=[]
for radial in [0,nrad-1]:
    for s in range(nsta-1):boundary.append([grid[radial,s],grid[radial,s+1],bycell[0 if radial==0 else 47,s//2]])
for station in [0,nsta-1]:
    for r in range(nrad-1):boundary.append([grid[r,station],grid[r+1,station],bycell[r//2,0 if station==0 else 1439]])
geometry=pack(DATA/'geometry.bin.gz',dict(xy=xy.astype('<f4'),cells=cells.astype('<u4'),boundary=np.array(boundary,dtype='<u4')))
geometry.update(nodes=len(xy),cells=len(cells),native_boundary_segments=len(boundary),coordinate_quantization_max_mm=float(abs(xy.astype('f4').astype('f8')-xy).max()))
cases=[];load_nodes=[]
for s in np.linspace(31,nsta-32,13).astype(int):load_nodes.append(int(grid[0,s//2*2+1]))
for s in [0,nsta-1]:
    for r in [23,47,71]:load_nodes.append(int(grid[r,s]))
for row in proof['cases']:
    cid=row['case'];folder=SRC/'cloud/results'/cid;d=np.load(folder/'fields.npz');h=read(folder/'download_hashes.json');assert sha(folder/'fields.npz')==h['fields.npz']['sha256']
    u=d['displacement_rtz'];stress=d['stress_rtz_voigt'];w=d['quadrature_volume_weight']
    vm=np.sqrt(.5*((stress[...,0]-stress[...,1])**2+(stress[...,1]-stress[...,2])**2+(stress[...,2]-stress[...,0])**2)+3*(stress[...,3:]**2).sum(-1))
    mean=(vm*w).sum(1)/w.sum(1)
    item=pack(DATA/(cid+'.bin.gz'),dict(U=u.astype('<f4'),VM=mean.astype('<f4')))
    item.update(row);item.update(raw_field_sha256=sha(folder/'fields.npz'),U_quantization_max_mm=float(abs(u.astype('f4').astype('f8')-u).max()),VM_quantization_max_MPa=float(abs(mean.astype('f4').astype('f8')-mean).max()))
    cases.append(item)
    if cid=='p01':
        forces=d['nodal_load_rtz'][load_nodes];assert np.all(np.linalg.norm(forces,axis=1)>0)
        loads=[dict(node=n,direction_rtz=(f/np.linalg.norm(f)).tolist(),source_force_N_at_3_5MPa=f.tolist(),kind='internal_pressure' if k<13 else 'closed_end_pressure_resultant') for k,(n,f) in enumerate(zip(load_nodes,forces))]
manifest=dict(schema='copv-axisym-display-v1',geometry=geometry,cases=cases,loads=loads,display=dict(coordinates='Original Q2 nodes; native midpoints retained; rotation about z is display tessellation, not 3D FE elements',U='Actual saved displacement of each V100 solve; float32 display copy only',stress='von Mises computed at the 9 native GP, then volume-weighted mean within each original cell. No cross-cell smoothing. Not an orthotropic failure criterion.',triangulation='Each Quad9 is split into 8 triangles through its original center and edge midpoints, for display only',load_arrows='Actual original nodal_load_rtz directions; pressure and its closed-end resultants; orange, schematic lengths, x-ray depth for visibility',common_ranges=dict(U=[0,max(x['max_U_mm'] for x in cases)],VM=[0,max(x['max_GP_VM_MPa'] for x in cases)])),scope=proof['physical_scope'],current_geometry_certified_for_70MPa=False,production_labels_released=0,input_sha256=sha(SRC/'cloud/homogeneous_E_input.npz'),source_sha256=proof['source_sha256'])
dump(DATA/'manifest.json',manifest)
for name in ['FINAL_VERIFICATION.json','case_summary.csv','pressure_response.png','pressure_response.pdf']:
    shutil.copy2(SRC/name,DATA/name)
for name in ['FINAL_VERIFICATION.json','material_contract.json']:
    shutil.copy2(WORK/'homogeneous_orthotropic_v1'/name,DATA/('homogeneous_'+name))
for name in ['mesh_convergence','same_position_stress_comparison','model_and_stress','hotspot_profile_and_difference']:
    for ext in ['png','pdf']:shutil.copy2(WORK/'homogeneous_orthotropic_v1/figures'/(name+'.'+ext),DATA/(name+'.'+ext))
dump(ROOT/'export_verification.json',dict(passed=True,cases=len(cases),native_nodes=len(xy),native_cells=len(cells),all_raw_SHA_verified=True,geometry=geometry,maximum_display_U_quantization_mm=max(x['U_quantization_max_mm'] for x in cases),maximum_display_VM_quantization_MPa=max(x['VM_quantization_max_MPa'] for x in cases),assets_bytes=sum(p.stat().st_size for p in DATA.iterdir())))
print(json.dumps(read(ROOT/'export_verification.json'),indent=2))
