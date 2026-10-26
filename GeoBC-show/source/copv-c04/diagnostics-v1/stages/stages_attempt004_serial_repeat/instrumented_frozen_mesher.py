"""Mesh an already repaired C02 partition without repeating CAD operations."""
import argparse,json,time
from pathlib import Path
import numpy as np
import gmsh
from scipy.optimize import linear_sum_assignment
from scipy import sparse
from scipy.sparse.csgraph import connected_components


def capture_stage(out, label, pid, started):
    import hashlib
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.asarray(node_tags); xyz = np.asarray(coords).reshape(-1,3)
    order = np.argsort(node_tags); ordered = node_tags[order]
    ideal = np.array([[1.,.5,.5],[0.,np.sqrt(3)/2,1/(2*np.sqrt(3))],[0.,0.,np.sqrt(2/3)]])
    invideal = np.linalg.inv(ideal)
    parts=[]; defects=[]; geometry_hash=hashlib.sha256(); topology_hash=hashlib.sha256()
    for _, entity in gmsh.model.getEntities(3):
        types, ids, conns = gmsh.model.mesh.getElements(3, entity)
        for typ, element_ids, conn in zip(types, ids, conns):
            if typ not in (4,11): raise RuntimeError('Unexpected tetra type '+str(typ))
            width=4 if typ==4 else 10
            native=np.asarray(conn).reshape(-1,width)
            indices=order[np.searchsorted(ordered,native)]
            v=xyz[indices[:,:4]]; j=np.swapaxes(v[:,1:]-v[:,:1],1,2)
            det=np.linalg.det(j); a=j@invideal
            cofactor=np.stack([np.cross(a[:,:,1],a[:,:,2]),np.cross(a[:,:,2],a[:,:,0]),np.cross(a[:,:,0],a[:,:,1])],axis=2)
            denominator=np.sqrt(np.sum(a*a,axis=(1,2))*np.sum(cofactor*cofactor,axis=(1,2)))
            sicn=np.divide(3*np.linalg.det(a),denominator,out=np.zeros_like(det),where=denominator>0)
            sort=np.argsort(element_ids); geometry_hash.update(v[sort].tobytes());topology_hash.update(native[sort,:4].tobytes())
            api_j=np.asarray(gmsh.model.mesh.getElementQualities(element_ids,'minDetJac'))
            bad=np.flatnonzero((det<=0)|(api_j<=0)); worst=np.unique(np.r_[np.argsort(sicn)[:20],np.argsort(api_j)[:20]])
            selected=np.unique(np.r_[bad,worst]); center=v[selected].mean(1)
            edge=np.stack([np.linalg.norm(v[:,b]-v[:,a],axis=1) for a,b in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]],axis=1)
            areas=np.stack([.5*np.linalg.norm(np.cross(v[:,b]-v[:,a],v[:,c]-v[:,a]),axis=1) for a,b,c in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]],axis=1)
            altitude=abs(det)/2/areas.max(1)
            midpoint_error=0.
            if width==10:
                for mid,(aa,bb) in enumerate([(0,1),(1,2),(0,2),(0,3),(2,3),(1,3)],4):
                    midpoint_error=max(midpoint_error,float(np.linalg.norm(xyz[indices[:,mid]]-.5*(v[:,aa]+v[:,bb]),axis=1).max()))
            part=dict(entity=int(entity),part_id=int(pid[entity]),cells=len(det),negative=int((det<0).sum()),zero=int((det==0).sum()),min_corner_det_mm3=float(det.min()),api_minDetJac_mm3=float(api_j.min()),api_nonpositive=int((api_j<=0).sum()),api_vs_corner_max_abs_difference=float(abs(api_j-det).max()),min_affine_SICN=float(sicn.min()),midpoint_max_deviation_mm=midpoint_error,mesh_volume_signed_mm3=float(det.sum()/6))
            parts.append(part)
            defects.extend([dict(part_id=int(pid[entity]),cell_id=int(element_ids[k]),det_mm3=float(det[k]),api_minDetJac_mm3=float(api_j[k]),sicn=float(sicn[k]),centroid_mm=v[k].mean(0).tolist(),min_altitude_mm=float(altitude[k]),min_edge_mm=float(edge[k].min()),max_edge_mm=float(edge[k].max()),corner_node_ids=native[k,:4].tolist(),corner_xyz_mm=v[k].tolist()) for k in selected])
    stage=dict(stage=label,elapsed_seconds=time.perf_counter()-started,parts=parts,nonpositive_cells=sum(p['negative']+p['zero'] for p in parts),corner_geometry_sha256=geometry_hash.hexdigest(),corner_topology_sha256=topology_hash.hexdigest(),defects_and_worst20_per_part=defects)
    (out/(label+'.json')).write_text(json.dumps(stage,indent=2)+'\n',encoding='utf-8')
    print('C04_DIAGNOSTIC_STAGE '+json.dumps({k:stage[k] for k in ['stage','elapsed_seconds','nonpositive_cells','corner_geometry_sha256','corner_topology_sha256']}),flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partition',type=Path,required=True);ap.add_argument('--review',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--h',type=float,default=6.);ap.add_argument('--curvature',type=int,default=64);ap.add_argument('--algorithm',type=int,default=1);ap.add_argument('--straight-geometry',action='store_true',help='P2 displacement on straight-sided TET10; do not claim curved CAD interpolation');ap.add_argument('--optimize-netgen',action='store_true');args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();report=json.loads(args.review.read_text(encoding='utf-8'))
    gmsh.initialize();gmsh.option.setNumber("General.NumThreads",1);gmsh.option.setNumber("Mesh.RandomSeed",1)
    try:
        gmsh.model.occ.importShapes(str(args.partition));occ=gmsh.model.occ;occ.synchronize();vols=gmsh.model.getEntities(3);records=report['volume_records'];assert len(vols)==len(records)
        mass=np.array([occ.getMass(3,t) for _,t in vols]);center=np.array([occ.getCenterOfMass(3,t) for _,t in vols]);expected=np.array([r['volume_mm3'] for r in records]);target=np.array([r['center_mm'] for r in records]);cost=(np.log(mass[:,None]/expected[None,:])*1000)**2+np.sum((center[:,None,:]-target[None,:,:])**2,axis=-1);rows,cols=linear_sum_assignment(cost);pid={vols[int(i)][1]:records[int(j)]['part_id'] for i,j in zip(rows,cols)}
        used={abs(t) for d,t in gmsh.model.getBoundary(vols,combined=False,oriented=True) if d==2};orphan=[(2,t) for d,t in gmsh.model.getEntities(2) if t not in used]
        if orphan:occ.remove(orphan,recursive=True);occ.synchronize()
        for p in sorted(set(pid.values())):
            tags=[t for t,v in pid.items() if v==p];gmsh.model.addPhysicalGroup(3,tags,p);gmsh.model.setPhysicalName(3,p,'part_'+str(p))
        for key,val in [('Mesh.MeshSizeMin',.5),('Mesh.MeshSizeMax',args.h),('Mesh.MeshSizeFromCurvature',args.curvature),('Mesh.MeshSizeExtendFromBoundary',0),('Mesh.Algorithm',6),('Mesh.Algorithm3D',args.algorithm),('Mesh.MaxNumThreads3D',1),('Mesh.MaxNumThreads2D',1)]:gmsh.option.setNumber(key,val)
        gmsh.model.mesh.generate(3)
        capture_stage(out,"01_after_generate",pid,start)
        gmsh.model.mesh.optimize('Relocate3D',niter=3)
        capture_stage(out,"02_after_Relocate3D",pid,start)
        if args.optimize_netgen or __import__('os').environ.get('COPV_C02_NETGEN')=='1':gmsh.model.mesh.optimize('Netgen');report['additional_optimizer']='Netgen'
        linear_done=time.perf_counter();gmsh.option.setNumber('Mesh.SecondOrderLinear',int(args.straight_geometry));gmsh.model.mesh.setOrder(2)
        capture_stage(out,"03_after_straight_order2",pid,start)
        return
        tags,coord,_=gmsh.model.mesh.getNodes();xyz=np.array(coord).reshape(-1,3);lookup={int(t):i for i,t in enumerate(tags)};cells=[];ids=[];pids=[];quality=[]
        for d,t in vols:
            types,elementids,conn=gmsh.model.mesh.getElements(3,t)
            for typ,e,c in zip(types,elementids,conn):
                if typ!=11:raise RuntimeError('Expected TET10')
                cells.append(np.fromiter((lookup[int(v)] for v in c),dtype=np.int64).reshape(-1,10));ids.extend(e);pids.extend([pid[t]]*len(e));quality.extend(gmsh.model.mesh.getElementQualities(e,'minDetJac'))
        cells=np.concatenate(cells);quality=np.array(quality)
        if (not np.isfinite(quality).all() or quality.min()<=0) and not args.straight_geometry:
            gmsh.model.mesh.optimize('HighOrderFastCurving');quality=np.array(gmsh.model.mesh.getElementQualities(ids,'minDetJac'))
            tags,coord,_=gmsh.model.mesh.getNodes();xyz=np.array(coord).reshape(-1,3)
        if not np.isfinite(quality).all() or quality.min()<=0:raise RuntimeError(f'Nonpositive curved Jacobians: {quality.min()}')
        used=np.unique(cells);remap=np.full(len(xyz),-1,np.int64);remap[used]=np.arange(len(used));cells=remap[cells];xyz=xyz[used];tags=np.array(tags)[used]
        corners=np.unique(cells[:,:4]);inv=np.full(len(xyz),-1,int);inv[corners]=np.arange(len(corners));edge=np.concatenate([cells[:,[0,1]],cells[:,[0,2]],cells[:,[0,3]],cells[:,[1,2]],cells[:,[1,3]],cells[:,[2,3]]]);graph=sparse.coo_matrix((np.ones(len(edge)),(inv[edge[:,0]],inv[edge[:,1]])),shape=(len(corners),len(corners)));components,_=connected_components(graph,directed=False)
        if components!=1:raise RuntimeError(f'Structural mesh has {components} disconnected components')
        np.savez_compressed(out/'mesh.npz',xyz=xyz,cells=cells,node_ids=tags,cell_ids=np.array(ids),cell_part_ids=np.array(pids),min_det_jac=quality)
        gmsh.write(str(out/'mesh.msh'));report.update(status='prepared',nodes=len(xyz),cells=len(cells),geometry_mapping='straight-sided TET10 coordinate storage (affine geometry); displacement degree selected by solver' if args.straight_geometry else 'curved TET10 / quadratic displacement',min_det_jac_mm3=float(quality.min()),connected_components=int(components),h_mm=args.h,curvature_points_per_circle=args.curvature,linear_meshing_seconds=linear_done-start,meshing_total_seconds=time.perf_counter()-start)
        report['C04_scope']='quarantined_research_only';report['mesher_sha256']=__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest()
        (out/'mesh_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:report[k] for k in ['status','nodes','cells','min_det_jac_mm3','connected_components','meshing_total_seconds']}),flush=True)
    finally:gmsh.finalize()
if __name__=='__main__':main()
