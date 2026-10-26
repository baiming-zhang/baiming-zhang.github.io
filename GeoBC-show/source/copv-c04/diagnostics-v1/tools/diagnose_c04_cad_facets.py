"""Compare original mesh facets with the same CAD faces; sampled diagnostic."""
from pathlib import Path
import argparse,json,time
import numpy as np
import gmsh
W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    combined=read(W/'data/copv/c04_research_v1/combined_review.json');results=[]
    for suffix in ['0002','0018']:
        tick=time.perf_counter();case=next(r for r in combined['cases'] if r['parent_case_id']=='copv_c04_'+suffix);root=W/case['artifact_root'];m=np.load(root/'mesh/mesh.npz');xyz=m['xyz'];nodeids=m['node_ids'];order=np.argsort(nodeids);ordered=nodeids[order];partreport={p['part_id']:p for p in case['quality']['parts']};meshreview=read(root/'mesh/mesh_review.json');partbyentity={r['tag']:r['part_id'] for r in meshreview['volume_records']}
        gmsh.initialize();gmsh.option.setNumber('General.Terminal',0)
        try:
            gmsh.open(str(root/'mesh/mesh.msh'));owners={};flux={p:0. for p in partreport};patches=[]
            sf=np.load(root/'mesh/surface/surface.npz');faces=sf['faces'][:,:3];faceparts=sf['part_ids'];nativefaces=nodeids[faces];v=xyz[faces];center=v.mean(1);area=.5*np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]);areas=np.linalg.norm(area,axis=1)
            for p in partreport:flux[p]=float(np.einsum('ij,ij->i',center[faceparts==p],area[faceparts==p]).sum()/3)
            covered=np.zeros(len(faces),bool)
            for _,face in gmsh.model.getEntities(2):
                member,_,_=gmsh.model.mesh.getNodes(2,face,includeBoundary=True);member=np.asarray(member)
                if len(member)<3:continue
                mask=np.isin(nativefaces[:,0],member)&np.isin(nativefaces[:,1],member)&np.isin(nativefaces[:,2],member);covered|=mask
                for p in np.unique(faceparts[mask]):
                    selection=np.flatnonzero(mask&(faceparts==p));selected=selection[np.linspace(0,len(selection)-1,min(240,len(selection)),dtype=int)]
                    patches.append(dict(face=face,count=len(selection),center=center[selected],area=area[selected],corner=v[selected,0],total_area=float(areas[selection].sum()),owners=[(int(p),1)]))
            if not covered.all():raise RuntimeError('CAD entity classification did not cover all surface facets: '+str((~covered).sum()))
            gmsh.clear();gmsh.model.add('same_CAD');gmsh.model.occ.importShapes(str(W/'data/copv/c04_pilot_v1/cases'/case['parent_case_id']/'partition_v2/structural_partition.brep'));gmsh.model.occ.synchronize();rows=[];approx={p:0. for p in partreport}
            for patch in patches:
                face=patch['face'];center=patch['center'];closest,_=gmsh.model.getClosestPoint(2,face,center.ravel());closest=np.asarray(closest).reshape(-1,3);corner_closest,_=gmsh.model.getClosestPoint(2,face,patch['corner'].ravel());corner_error=np.linalg.norm(np.asarray(corner_closest).reshape(-1,3)-patch['corner'],axis=1)
                area=patch['area'];areas=np.linalg.norm(area,axis=1);normal=area/areas[:,None];offset=np.einsum('ij,ij->i',closest-center,normal);mean=float(np.sum(offset*areas)/areas.sum());contribution=mean*patch['total_area']
                for p,sign in patch['owners']:approx[p]+=sign*contribution
                rows.append(dict(face=face,CAD_type=gmsh.model.getType(2,face),owners=patch['owners'],triangles=patch['count'],sampled=len(center),area_mm2=patch['total_area'],max_corner_projection_error_mm=float(corner_error.max()),max_centroid_CAD_distance_mm=float(np.linalg.norm(closest-center,axis=1).max()),signed_normal_offset_area_weighted_mm=mean,estimated_volume_correction_mm3=contribution,centroid_z_range_mm=[float(center[:,2].min()),float(center[:,2].max())],centroid_r_range_mm=[float(np.linalg.norm(center[:,:2],axis=1).min()),float(np.linalg.norm(center[:,:2],axis=1).max())]))
            parts=[]
            for p,report in partreport.items():parts.append(dict(part_id=p,name=report['name'],mesh_volume_mm3=report['mesh_volume_mm3'],CAD_volume_mm3=report['CAD_volume_mm3'],CAD_minus_mesh_mm3=report['CAD_volume_mm3']-report['mesh_volume_mm3'],boundary_flux_volume_mm3=flux[p],boundary_flux_relative_difference=abs(flux[p]/report['mesh_volume_mm3']-1),sampled_first_order_volume_correction_mm3=approx[p]))
            result=dict(case_id=case['parent_case_id'],seconds=time.perf_counter()-tick,scope='Facet-to-same-CAD-face projection at up to 240 deterministic samples/face. First-order area-weighted normal-offset estimate is diagnostic, not an exact volume integration or certified Hausdorff distance.',parts=parts,faces=rows)
            (out/(case['parent_case_id']+'.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');results.append(result);print(json.dumps(dict(case_id=result['case_id'],seconds=result['seconds'],faces=len(rows),max_corner_error=max(x['max_corner_projection_error_mm'] for x in rows))),flush=True)
        finally:gmsh.finalize()
    (out/'summary.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
