"""Reproducible C01 geometry delivery; no simulation labels or trained models."""
from pathlib import Path
import argparse, base64, csv, gzip, hashlib, json, time
import numpy as np
import gmsh
from geobc.copv_geometry import generate_copv_geometry
from geobc.copv_dataset import geometry_eligibility, physical_parameters, independent_parameters
from geobc.copv_structured_mesh import analytic_segments

REPO=Path(__file__).resolve().parents[1]
WS=REPO.parents[1]

def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def surface(rz,nphi=96):
    phi=np.arange(nphi)*2*np.pi/nphi
    xyz=np.stack([rz[:,0,None]*np.cos(phi),rz[:,0,None]*np.sin(phi),np.broadcast_to(rz[:,1,None],(len(rz),nphi))],axis=-1).reshape(-1,3)
    i,j=np.indices((len(rz)-1,nphi));a=i*nphi+j;b=i*nphi+(j+1)%nphi;c=b+nphi;d=a+nphi
    tri=np.stack([a,b,c,a,c,d],axis=-1).reshape(-1,3)
    return xyz,tri

def encoded(a,dtype):return base64.b64encode(np.asarray(a,dtype=dtype).tobytes()).decode('ascii')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    start=time.perf_counter();out=args.output
    config=json.loads((REPO/'configs/copv_dataset_candidate.json').read_text())['geometry']
    ranges=config['candidate_independent_ranges'];base=config['baseline_parameters_mm'];gates=config['eligibility_gates_before_any_solve']
    designs=[('baseline','基准',base)]
    for i,(key,bounds) in enumerate(ranges.items()):
        for side,value in zip(('min','max'),bounds):
            independent=independent_parameters(base);independent[i]=value
            designs.append((f'{key}_{side}',key+' '+side,physical_parameters(independent)))
    records=[];viewer=[]
    for name,label,p in designs:
        tick=time.perf_counter();dest=out/name
        eligibility=geometry_eligibility(p,gates)
        if not eligibility['passed']:raise RuntimeError((name,eligibility))
        report=generate_copv_geometry(dest,p,dome_samples=101)
        gmsh.initialize()
        try:
            gmsh.option.setNumber('General.Terminal',0)
            gmsh.model.occ.importShapes(str(dest/'geometry.step'));gmsh.model.occ.synchronize()
            volumes=sorted(gmsh.model.occ.getMass(3,t) for _,t in gmsh.model.getEntities(3))
            expected=sorted(report['checks']['positive_material_volumes_mm3'].values())
            volume_error=float(np.max(np.abs(np.array(volumes)/expected-1))) if len(volumes)==2 else float('inf')
            if len(volumes)!=2 or volume_error>1e-6:raise RuntimeError('STEP round trip failed')
        finally:gmsh.finalize()
        report['c01_review']={'eligibility':eligibility,'step_roundtrip_volume_count':len(volumes),'step_roundtrip_max_relative_volume_error':volume_error}
        write(dest/'geometry_manifest.json',report)
        rz,normal=analytic_segments(p,(2,32,24,32,2));parts=[]
        for offset in (0,p['liner_thickness_mm'],p['liner_thickness_mm']+p['overwrap_thickness_mm']):
            xyz,tri=surface(rz+offset*normal);parts.append({'xyz':encoded(xyz,'<f4'),'tri':encoded(tri,'<u4')})
        row={'case_id':name,'label':label,'parameters_mm':p,'checks':report['checks'],'review':report['c01_review'],'seconds':time.perf_counter()-tick}
        records.append(row)
        asset={'record':row,'surfaces':parts,'meridian':{k:v.tolist() for k,v in [('inner',rz),('interface',rz+p['liner_thickness_mm']*normal),('outer',rz+(p['liner_thickness_mm']+p['overwrap_thickness_mm'])*normal)]}}
        asset_path=out/'viewer'/f'{name}.json.gz';asset_path.parent.mkdir(parents=True,exist_ok=True)
        asset_path.write_bytes(gzip.compress(json.dumps(asset,separators=(',',':'),ensure_ascii=False).encode(),mtime=0))
        viewer.append({'id':name,'label':label,'asset':f'../assets/copv-c01/{name}.json.gz'})
        print(json.dumps({'case':name,'passed':True,'seconds':row['seconds']}),flush=True)
    # Quad9 meridian source for independent FEniCSx geometry import smoke.
    counts=(2,16,16,16,2);rz,n=analytic_segments(base,counts);offsets=np.r_[np.linspace(0,base['liner_thickness_mm'],3),base['liner_thickness_mm']+np.linspace(0,base['overwrap_thickness_mm'],5)[1:]]
    shape=(len(offsets),len(rz));ij=np.indices(shape).reshape(2,-1).T;k,s=ij.T;xy=rz[s]+offsets[k,None]*n[s]
    ref=np.array([[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[1,0],[0,1],[-1,0],[0,0]])
    grid=np.arange(len(xy)).reshape(shape);e=np.indices((3,sum(counts))).reshape(2,-1).T;idx=2*e[:,None,:]+ref[None,:,:]+1;cells=grid[idx[:,:,0],idx[:,:,1]]
    np.savez_compressed(out/'fenicsx_geometry.npz',xy=xy,cells=cells,node_ij=ij,cell_ij=e,gmsh_reference=ref,parameters=np.array(list(base.values())))
    summary={'schema':'copv_c01_geometry_v1','cases':records,'case_count':len(records),'all_passed':True,'parameter_ranges':ranges,'range_rule':'R in mm; all other independent parameters are ratios to R. Not every Cartesian-product combination is certified; eligibility required.','units':'mm, mm^3','profile':'C1 cubic Bezier, not isotensoid; constant analytic normal thickness; two open necks; no bosses/contact/winding manufacturability','viewer_geometry':'Analytic profile revolution tessellation, not FE mesh; quarter cutaway is display-only. CAD interpolation deviation recorded per case.','seconds':time.perf_counter()-start,'viewer_cases':viewer}
    write(out/'summary.json',summary)
    with (out/'parameters-and-quality.csv').open('w',newline='',encoding='utf-8-sig') as f:
        fields=['case_id',*base,'liner_volume_mm3','wrap_volume_mm3','offset_curvature_ratio','CAD_profile_error_mm','STEP_relative_volume_error','seconds'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in records:
            c=r['checks'];w.writerow(dict(case_id=r['case_id'],**r['parameters_mm'],liner_volume_mm3=c['positive_material_volumes_mm3']['liner_volume'],wrap_volume_mm3=c['positive_material_volumes_mm3']['overwrap_volume'],offset_curvature_ratio=c['max_offset_curvature_ratio'],CAD_profile_error_mm=max(c['midpoint_CAD_profile_error_mm']),STEP_relative_volume_error=r['review']['step_roundtrip_max_relative_volume_error'],seconds=r['seconds']))
    print(json.dumps({'completed':len(records),'seconds':summary['seconds']}),flush=True)

if __name__=='__main__':main()
