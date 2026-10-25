"""Mesh an already repaired C02 partition without repeating CAD operations."""
import argparse,json,time
from pathlib import Path
import numpy as np
import gmsh
from scipy.optimize import linear_sum_assignment
from scipy import sparse
from scipy.sparse.csgraph import connected_components

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partition',type=Path,required=True);ap.add_argument('--review',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--h',type=float,default=6.);ap.add_argument('--curvature',type=int,default=64);ap.add_argument('--algorithm',type=int,default=1);ap.add_argument('--straight-geometry',action='store_true',help='P2 displacement on straight-sided TET10; do not claim curved CAD interpolation');ap.add_argument('--optimize-netgen',action='store_true');ap.add_argument('--min-size',type=float,default=.15);ap.add_argument('--bond-h',type=float,default=2.);ap.add_argument('--collar-h',type=float,default=2.);ap.add_argument('--collar-root-h',type=float,default=.2);args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();report=json.loads(args.review.read_text(encoding='utf-8'))
    gmsh.initialize()
    try:
        gmsh.model.occ.importShapes(str(args.partition));occ=gmsh.model.occ;occ.synchronize();vols=gmsh.model.getEntities(3);records=report['volume_records'];assert len(vols)==len(records)
        mass=np.array([occ.getMass(3,t) for _,t in vols]);center=np.array([occ.getCenterOfMass(3,t) for _,t in vols]);expected=np.array([r['volume_mm3'] for r in records]);target=np.array([r['center_mm'] for r in records]);cost=(np.log(mass[:,None]/expected[None,:])*1000)**2+np.sum((center[:,None,:]-target[None,:,:])**2,axis=-1);rows,cols=linear_sum_assignment(cost);pid={vols[int(i)][1]:records[int(j)]['part_id'] for i,j in zip(rows,cols)}
        used={abs(t) for d,t in gmsh.model.getBoundary(vols,combined=False,oriented=True) if d==2};orphan=[(2,t) for d,t in gmsh.model.getEntities(2) if t not in used]
        if orphan:occ.remove(orphan,recursive=True);occ.synchronize()
        for p in sorted(set(pid.values())):
            tags=[t for t,v in pid.items() if v==p];gmsh.model.addPhysicalGroup(3,tags,p);gmsh.model.setPhysicalName(3,p,'part_'+str(p))
        for key,val in [('Mesh.MeshSizeMin',args.min_size),('Mesh.MeshSizeMax',args.h),('Mesh.MeshSizeFromCurvature',args.curvature),('Mesh.MeshSizeExtendFromBoundary',0),('Mesh.Algorithm',6),('Mesh.Algorithm3D',args.algorithm),('Mesh.MaxNumThreads3D',8),('Mesh.MaxNumThreads2D',8)]:gmsh.option.setNumber(key,val)
        params=report['parameters_mm'];fields=[]
        def restrict(base,volumes):
            f=gmsh.model.mesh.field.add('Restrict');gmsh.model.mesh.field.setNumber(f,'InField',base);gmsh.model.mesh.field.setNumbers(f,'VolumesList',volumes);gmsh.model.mesh.field.setNumber(f,'IncludeBoundary',1);fields.append(f)
        for role,h in [(2,args.bond_h),(7,args.collar_h),(8,args.collar_h),(9,1.5),(10,1.5),(14,1.5),(18,1.5)]:
            f=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(f,'F',str(h));restrict(f,[t for t,p in pid.items() if p==role])
        for sign,role in [(1,7),(-1,8)]:
            volume=next(t for t,p in pid.items() if p==role);bb=gmsh.model.getBoundingBox(3,volume);startz=bb[2] if sign>0 else -bb[5]
            for width,h in [(1.2*params['inner_radius_mm']/100,args.collar_root_h),(5*params['inner_radius_mm']/100,.45)]:
                z0,z1=sorted([sign*(startz-.1),sign*(startz+width)]);f=gmsh.model.mesh.field.add('Box')
                for key,value in [('VIn',h),('VOut',args.h),('XMin',-2*params['inner_radius_mm']),('XMax',2*params['inner_radius_mm']),('YMin',-2*params['inner_radius_mm']),('YMax',2*params['inner_radius_mm']),('ZMin',z0),('ZMax',z1)]:gmsh.model.mesh.field.setNumber(f,key,value)
                restrict(f,[volume])
        f=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(f,'FieldsList',fields);gmsh.model.mesh.field.setAsBackgroundMesh(f)
        report['C04_sizing_recipe']=dict(bond_h_mm=args.bond_h,collar_h_mm=args.collar_h,collar_root_h_mm=args.collar_root_h,min_size_mm=args.min_size,documentation='https://gmsh.info/doc/texinfo/#Gmsh-mesh-size-fields')
        # Resolve small turning/groove arcs without curvature-driven refinement
        # over their entire adjacent toroidal surface.
        for _,curve in gmsh.model.getEntities(1):
            length=occ.getMass(1,curve)
            if gmsh.model.getType(1,curve)=='Circle' and .5<length<4.:gmsh.model.mesh.setTransfiniteCurve(curve,12)
        gmsh.model.mesh.generate(3)
        if args.optimize_netgen or __import__('os').environ.get('COPV_C02_NETGEN')=='1':gmsh.model.mesh.optimize('Netgen');report['additional_optimizer']='Netgen'
        linear_done=time.perf_counter();gmsh.option.setNumber('Mesh.SecondOrderLinear',int(args.straight_geometry));gmsh.model.mesh.setOrder(2)
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
        (out/'mesh_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:report[k] for k in ['status','nodes','cells','min_det_jac_mm3','connected_components','meshing_total_seconds']}),flush=True)
    finally:gmsh.finalize()
if __name__=='__main__':main()
