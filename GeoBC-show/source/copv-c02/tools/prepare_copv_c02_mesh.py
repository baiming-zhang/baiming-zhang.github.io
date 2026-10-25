"""Prepare a conformal structural mesh from the accepted C01 CAD baseline.

Elastomer/PEEK seals are display-only in this ideal-bonded elastic baseline;
seal compression, friction, leakage and failure are outside this contract.
"""
import argparse,base64,gzip,hashlib,json,time
from pathlib import Path
import numpy as np
import gmsh
from scipy.optimize import linear_sum_assignment
from geobc.copv_structured_mesh import analytic_segments

def write(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cad',type=Path,required=True);ap.add_argument('--model',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--h',type=float,default=14.);ap.add_argument('--inspect-only',action='store_true');args=ap.parse_args();start=time.perf_counter();out=args.output;out.mkdir(parents=True,exist_ok=True)
    source=json.loads(gzip.decompress(args.model.read_bytes()));parts=source['parts']
    gmsh.initialize()
    try:
        gmsh.model.add('COPV_C02_structural_ideal_bonded');occ=gmsh.model.occ;occ.importShapes(str(args.cad));occ.synchronize();vols=gmsh.model.getEntities(3)
        if len(vols)!=18:raise RuntimeError('Accepted CAD no longer has 18 solids')
        mass=np.array([occ.getMass(3,t) for _,t in vols]);com=np.array([occ.getCenterOfMass(3,t) for _,t in vols]);expected=np.array([p['volume_mm3'] for p in parts]);positions=[np.frombuffer(base64.b64decode(p['xyz']),'<f4').reshape(-1,3) for p in parts];centers=np.array([(x.min(0)+x.max(0))/2 for x in positions])
        cost=(np.log(mass[:,None]/expected[None,:])*1000)**2+np.linalg.norm(com[:,None,:]-centers[None,:,:],axis=-1)**2/100
        rows,cols=linear_sum_assignment(cost);assignment={int(c):vols[int(r)] for r,c in zip(rows,cols)}
        roles=[];remove=[]
        for i,part in enumerate(parts):
            if abs(mass[vols.index(assignment[i])]/expected[i]-1)>1e-4:raise RuntimeError('CAD identity/volume mismatch')
            entry={k:v for k,v in part.items() if k not in ('xyz','normal','tri')};entry.update(part_id=i+1,source_volume_tag=assignment[i][1])
            if part['material'] in ('elastomer','peek'):remove.append(assignment[i]);entry['structural']=False
            else:entry['structural']=True
            roles.append(entry)
        occ.remove(remove,recursive=True);occ.synchronize()
        # Fit the nominal collar inner interface to the exact spline surface.
        # C01's sampled polygon is fine for display but oscillates across this
        # surface. Retain its outer envelope and boss trimming, replacing only
        # the mating interface by an exact boolean against the outer envelope.
        params=source['parameters_mm'];rz,n=analytic_segments(params,(2,50,16,50,2));offset=params['liner_thickness_mm']+.15+params['overwrap_thickness_mm'];profile=rz+offset*n
        cuts=np.cumsum([0,4,100,32,100,4]);cache={};curves=[]
        def point(q):
            k=tuple(np.round(q,11))
            if k not in cache:cache[k]=occ.addPoint(float(q[0]),0,float(q[1]))
            return cache[k]
        for i in range(5):
            a=profile[cuts[i]:cuts[i+1]+1];pts=[point(q) for q in a];curves.append(occ.addSpline(pts) if i in (1,3) else occ.addLine(pts[0],pts[-1]))
        at=occ.addPoint(0,0,float(profile[-1,1]));ab=occ.addPoint(0,0,float(profile[0,1]));loop=occ.addCurveLoop(curves+[occ.addLine(point(profile[-1]),at),occ.addLine(at,ab),occ.addLine(ab,point(profile[0]))]);face=occ.addPlaneSurface([loop]);enclosure=[dt for dt in occ.revolve([(2,face)],0,0,0,0,0,1,2*np.pi) if dt[0]==3]
        bosses=[assignment[i] for i,p in enumerate(parts) if p['group']=='boss'];repairs=[]
        for i,part in enumerate(parts):
            if part['group']!='reinforcement':continue
            sign=1 if np.frombuffer(base64.b64decode(part['xyz']),'<f4').reshape(-1,3)[:,2].mean()>0 else -1
            shoulder=params['straight_length_mm']/2+params['dome_depth_mm'];cutz=shoulder-8
            ids=np.flatnonzero((sign*rz[:,1]>=params['straight_length_mm']/2+.45*params['dome_depth_mm'])&(sign*rz[:,1]<=cutz+3));a=profile[ids];q=np.clip((sign*rz[ids,1]-(params['straight_length_mm']/2+.45*params['dome_depth_mm']))/(.38*params['dome_depth_mm']),0,1);extra=1.4*np.sin(q*np.pi/2)**2+.03;outer=a+extra[:,None]*n[ids]
            poly=np.vstack([[0,outer[0,1]],outer,[0,outer[-1,1]]]);pt=[occ.addPoint(float(r),0,float(z)) for r,z in poly];ln=[occ.addLine(pt[j],pt[(j+1)%len(pt)]) for j in range(len(pt))];f=occ.addPlaneSurface([occ.addCurveLoop(ln)]);raw=[dt for dt in occ.revolve([(2,f)],0,0,0,0,0,1,2*np.pi) if dt[0]==3]
            fixed=[dt for dt in occ.cut(raw,enclosure+bosses,removeTool=False)[0] if dt[0]==3]
            if len(fixed)!=1:raise RuntimeError('Collar interface fit did not retain one solid')
            before=occ.getMass(*assignment[i]);after=occ.getMass(*fixed[0]);occ.remove([assignment[i]],recursive=True);assignment[i]=fixed[0]
            repairs.append(dict(part_id=i+1,old_volume_mm3=before,new_volume_mm3=after,relative_change=abs(after/before-1),rule='retain C01 outer collar envelope; fit inner interface to exact outer-ply spline and original boss pockets'))
        occ.remove(enclosure,recursive=True);occ.synchronize();active=[assignment[i] for i,p in enumerate(parts) if p['material'] not in ('elastomer','peek')]
        result,mapping=occ.fragment(active[:1],active[1:]);occ.synchronize();valid=[dt for dt in result if dt[0]==3]
        tagged={};ambiguous=[]
        for entry,mapped in zip([r for r in roles if r['structural']],mapping):
            for d,t in mapped:
                if d!=3:continue
                if t in tagged and tagged[t]!=entry['part_id']:ambiguous.append((t,tagged[t],entry['part_id'],occ.getMass(3,t)))
                else:tagged[t]=entry['part_id']
        # C01 polar collars used sampled straight segments on a spline surface.
        # Resolve only their tiny same-carbon overlap by preserving the outer
        # structural ply and trimming the reinforcement. Other overlaps fail.
        policy=[]
        byid={r['part_id']:r for r in roles}
        for tag,first,second,volume in ambiguous:
            a,b=byid[first],byid[second]
            if a['material']!='carbon' or b['material']!='carbon' or {a['group'],b['group']}!={'ply','reinforcement'}:
                raise RuntimeError('Unsupported material overlap; no automatic ownership choice')
            owner=first if a['group']=='ply' else second;tagged[tag]=owner
            policy.append(dict(volume_tag=tag,volume_mm3=volume,kept_part_id=owner,rule='preserve outer ply; trim sampled polar collar overlap'))
        occ.synchronize();valid=gmsh.model.getEntities(3)
        volume_records=[dict(tag=t,part_id=tagged.get(t),volume_mm3=occ.getMass(3,t),center_mm=occ.getCenterOfMass(3,t)) for _,t in valid]
        report={'status':'CAD_partition_inspection','source_CAD_sha256':hashlib.sha256(args.cad.read_bytes()).hexdigest(),'source_model_sha256':hashlib.sha256(args.model.read_bytes()).hexdigest(),'parts':roles,'collar_interface_repairs':repairs,'partition_volume_count':len(valid),'volume_records':volume_records,'original_ambiguous_volume_ownership':ambiguous,'resolved_overlap_policy':policy,'ambiguous_volume_ownership':[],'gmsh_version':gmsh.__version__,'h_mm':args.h,'seconds':time.perf_counter()-start}
        write(out/'partition_review.json',report)
        gmsh.write(str(out/'structural_partition.brep'))
        print(json.dumps({'partition_volumes':len(valid),'resolved_overlap_mm3':sum(a['volume_mm3'] for a in policy),'minimum_volume_mm3':min(r['volume_mm3'] for r in volume_records),'seconds':report['seconds']}),flush=True)
        if args.inspect_only:return
        # Any real overlap needs explicit resolution, rather than silently assigning material.
        if report['ambiguous_volume_ownership']:raise RuntimeError('Material ownership overlaps; inspect partition_review.json')
        for entry in roles:
            if not entry['structural']:continue
            tags=[t for t,pid in tagged.items() if pid==entry['part_id']]
            if not tags:raise RuntimeError('Lost structural part '+entry['name'])
            gmsh.model.addPhysicalGroup(3,tags,entry['part_id']);gmsh.model.setPhysicalName(3,entry['part_id'],f'part_{entry["part_id"]}_{entry["material"]}')
        # All external structural boundary facets; wetted/cap tagging follows exact
        # original-node connectivity and the documented fluid cavity definition.
        for key,val in [('Mesh.MeshSizeMin',.7),('Mesh.MeshSizeMax',args.h),('Mesh.MeshSizeFromCurvature',24),('Mesh.MeshSizeExtendFromBoundary',0),('Mesh.Algorithm',6),('Mesh.Algorithm3D',10),('Mesh.ElementOrder',1),('Mesh.MaxNumThreads3D',8),('Mesh.MaxNumThreads2D',8)]:gmsh.option.setNumber(key,val)
        gmsh.model.mesh.generate(3);gmsh.model.mesh.setOrder(2)
        nodes,coords,_=gmsh.model.mesh.getNodes();xyz=np.array(coords).reshape(-1,3);lookup={int(t):i for i,t in enumerate(nodes)};cells=[];ids=[];pids=[];quality=[]
        for _,tag in valid:
            typ,etags,conn=gmsh.model.mesh.getElements(3,tag)
            for t,e,c in zip(typ,etags,conn):
                if t!=11:raise RuntimeError('Expected native Gmsh TET10')
                cells.append(np.array([lookup[int(v)] for v in c]).reshape(-1,10));ids.extend(e);pids.extend([tagged[tag]]*len(e));quality.extend(gmsh.model.mesh.getElementQualities(e,'minDetJac'))
        cells=np.concatenate(cells);quality=np.array(quality)
        if not np.isfinite(quality).all() or quality.min()<=0:raise RuntimeError('Nonpositive curved TET10 Jacobian')
        np.savez_compressed(out/'mesh.npz',xyz=xyz,cells=cells,node_ids=nodes,cell_ids=np.array(ids),cell_part_ids=np.array(pids),min_det_jac=quality)
        gmsh.write(str(out/'mesh.msh'));report.update(status='mesh_prepared',nodes=len(xyz),cells=len(cells),min_det_jac_mm3=float(quality.min()),seconds=time.perf_counter()-start);write(out/'mesh_review.json',report);print(json.dumps({k:report[k] for k in ['status','nodes','cells','min_det_jac_mm3','seconds']}),flush=True)
    finally:gmsh.finalize()

if __name__=='__main__':main()
