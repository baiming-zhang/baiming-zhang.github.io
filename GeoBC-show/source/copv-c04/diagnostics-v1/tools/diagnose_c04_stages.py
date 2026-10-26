"""Instrument the frozen C04 mesher; never modify or resume an original attempt."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time
from datetime import datetime
from zoneinfo import ZoneInfo

W = Path(__file__).resolve().parents[1]
FROZEN = W / 'data/copv/c04_research_retry_v1/frozen_sources/mesh_copv_c04_research.py'

CAPTURE = r'''
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

'''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--threads',type=int,default=8);ap.add_argument('--seed',type=int);args=ap.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    text=FROZEN.read_text(encoding='utf-8')
    text=text.replace("('Mesh.MaxNumThreads3D',8)","('Mesh.MaxNumThreads3D',"+str(args.threads)+")").replace("('Mesh.MaxNumThreads2D',8)","('Mesh.MaxNumThreads2D',"+str(args.threads)+")")
    settings='gmsh.option.setNumber("General.NumThreads",'+str(args.threads)+')'
    if args.seed is not None:settings+=';gmsh.option.setNumber("Mesh.RandomSeed",'+str(args.seed)+')'
    text=text.replace('    gmsh.initialize()','    gmsh.initialize();'+settings,1)
    text=text.replace('def main():',CAPTURE+'\ndef main():',1)
    text=text.replace('        gmsh.model.mesh.generate(3)','        gmsh.model.mesh.generate(3)\n        capture_stage(out,"01_after_generate",pid,start)',1)
    text=text.replace("        gmsh.model.mesh.optimize('Relocate3D',niter=3)","        gmsh.model.mesh.optimize('Relocate3D',niter=3)\n        capture_stage(out,\"02_after_Relocate3D\",pid,start)",1)
    text=text.replace('gmsh.model.mesh.setOrder(2)','gmsh.model.mesh.setOrder(2)\n        capture_stage(out,"03_after_straight_order2",pid,start)',1)
    # A diagnostic attempt stops after its stage records; do not export another FE label.
    text=text.replace('        tags,coord,_=gmsh.model.mesh.getNodes();xyz=', '        return\n        tags,coord,_=gmsh.model.mesh.getNodes();xyz=',1)
    script=output/'instrumented_frozen_mesher.py';script.write_text(text,encoding='utf-8')
    env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join(str(W/x) for x in ['.rod_deps','.local_deps']);env.pop('COPV_C02_NETGEN',None)
    manifest=dict(created=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),threads=args.threads,seed=args.seed,scope='Diagnostic only; no original outputs changed; no FE label; no CAD edit',source=str(FROZEN.relative_to(W)),source_sha256=hashlib.sha256(FROZEN.read_bytes()).hexdigest(),instrumented_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),cases=[])
    for suffix in ['0000','0011']:
        case='copv_c04_'+suffix;partition=W/'data/copv/c04_pilot_v1/cases'/case/'partition_v2';dest=output/case;dest.mkdir()
        cmd=[sys.executable,'-B',str(script),'--partition',str(partition/'structural_partition.brep'),'--review',str(partition/'partition_review.json'),'--output',str(dest),'--h','6','--curvature','64','--algorithm','1','--straight-geometry']
        tick=time.perf_counter()
        with (dest/'execution.log').open('w',encoding='utf-8') as log:
            run=subprocess.run(cmd,env=env,cwd=W,stdout=log,stderr=subprocess.STDOUT,timeout=360,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        record=dict(case_id=case,exit_code=run.returncode,seconds=time.perf_counter()-tick,command=cmd,partition_sha256=hashlib.sha256((partition/'structural_partition.brep').read_bytes()).hexdigest())
        manifest['cases'].append(record);(output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8');print(json.dumps(record),flush=True)

if __name__=='__main__':main()
