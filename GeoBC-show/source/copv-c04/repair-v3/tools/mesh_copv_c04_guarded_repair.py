"""Versioned C04 mesh repair: deterministic CAD meshing and guarded relocation."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import gmsh
from scipy.optimize import linear_sum_assignment
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from geobc.copv_structured_mesh import analytic_segments

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partition',type=Path,required=True);ap.add_argument('--review',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--h',type=float,default=5.);ap.add_argument('--shell-h',type=float,default=4.);ap.add_argument('--bond-h',type=float,default=3.);ap.add_argument('--root-h',type=float,default=.75);ap.add_argument('--circle-segments',type=int,default=128);ap.add_argument('--recipe',choices=['matched','baseline'],default='matched');ap.add_argument('--seed',type=int,default=1);ap.add_argument('--max-cells',type=int,default=2000000);ap.add_argument('--skip-relocate',action='store_true');ap.add_argument('--adaptive-cloud',type=Path);ap.add_argument('--linear-checkpoint',type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();report=json.loads(args.review.read_text(encoding='utf-8'));stages=[]
    def write(name,value):
        (out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    def state(name,**kw):
        r=dict(stage=name,seconds=time.perf_counter()-start,**kw);stages.append(r);write('stages.json',stages);print(json.dumps(r),flush=True)
    def snapshot():
        tags,coords,_=gmsh.model.mesh.getNodes();return np.asarray(tags),np.asarray(coords).reshape(-1,3)
    def measure():
        eids=np.asarray(gmsh.model.mesh.getElements(3)[1][0]);j=np.asarray(gmsh.model.mesh.getElementQualities(eids,'minDetJac'));q=np.asarray(gmsh.model.mesh.getElementQualities(eids,'minSICN'));return dict(cells=len(eids),min_J=float(j.min()),nonpositive=int((j<=0).sum()),min_SICN=float(q.min()),below_005=int((q<.05).sum()))
    # Select the volume algorithm through a recorded environment override.
    algorithm=int(__import__('os').environ.get('COPV_C04_VOLUME_ALGORITHM','1'))
    if algorithm not in [1,10]:raise ValueError('Only declared Delaunay/HXT comparison is supported')
    gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1);gmsh.option.setNumber('General.AbortOnError',0);gmsh.logger.start()
    try:
        occ=gmsh.model.occ;occ.importShapes(str(args.partition));occ.synchronize();vols=gmsh.model.getEntities(3);records=report['volume_records'];assert len(vols)==len(records)
        mass=np.array([occ.getMass(3,t) for _,t in vols]);center=np.array([occ.getCenterOfMass(3,t) for _,t in vols]);expected=np.array([r['volume_mm3'] for r in records]);target=np.array([r['center_mm'] for r in records]);rows,cols=linear_sum_assignment((1000*np.log(mass[:,None]/expected[None,:]))**2+np.sum((center[:,None,:]-target[None,:,:])**2,axis=-1));pid={vols[int(i)][1]:records[int(j)]['part_id'] for i,j in zip(rows,cols)}
        used={abs(t) for d,t in gmsh.model.getBoundary(vols,combined=False,oriented=True) if d==2};orphan=[(2,t) for d,t in gmsh.model.getEntities(2) if t not in used]
        if orphan:occ.remove(orphan,recursive=True);occ.synchronize()
        for p in sorted(set(pid.values())):
            tags=[t for t,v in pid.items() if v==p];gmsh.model.addPhysicalGroup(3,tags,p);gmsh.model.setPhysicalName(3,p,'part_'+str(p))
        options={'Mesh.MeshSizeMin':.15,'Mesh.MeshSizeMax':args.h,'Mesh.MeshSizeFromCurvature':0,'Mesh.MeshSizeFromPoints':0,'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.Algorithm':6,'Mesh.Algorithm3D':algorithm,'Mesh.MaxNumThreads3D':1,'Mesh.MaxNumThreads2D':1,'Mesh.RandomSeed':args.seed,'Mesh.SecondOrderLinear':1,'Mesh.Binary':1,'Mesh.SaveAll':1}
        for key,val in options.items():gmsh.option.setNumber(key,val)
        params=report['parameters_mm'];roots=[];root_band=float(__import__('os').environ.get('COPV_C04_ROOT_BAND_MM','.75'));root_growth=float(__import__('os').environ.get('COPV_C04_ROOT_GROWTH','.6'));matched_root_band=__import__('os').environ.get('COPV_C04_MATCH_ROOT_CIRCLES')=='1';frontal_parts=set(map(int,__import__('os').environ.get('COPV_C04_FRONTAL_PARTS','7,8').split(',')));adjacent_band=__import__('os').environ.get('COPV_C04_ADJACENT_ROOT_BAND')=='1'
        if not 0<root_growth<=1.8:raise ValueError('Root size growth outside declared range')
        for sign,role in [(1,7),(-1,8)]:
            role_vols=[t for t,p in pid.items() if p==role];bounds=[gmsh.model.getBoundingBox(3,t) for t in role_vols];z0=min(b[2] for b in bounds) if sign>0 else -max(b[5] for b in bounds)
            # Exact collar boundary circle at its initial axial plane.
            curves={abs(c) for volume in role_vols for d,f in gmsh.model.getBoundary([(3,volume)],combined=False) for d,c in gmsh.model.getBoundary([(2,abs(f))],combined=False)}
            candidates=[]
            for c in curves:
                box=gmsh.model.getBoundingBox(1,c)
                if gmsh.model.getType(1,c)=='Circle' and abs(abs((box[2]+box[5])/2)-z0)<1e-4:candidates.append(occ.getMass(1,c)/(2*np.pi))
            radius=max(candidates) if candidates else np.sqrt(max(0,params['inner_radius_mm']**2-(z0-params['straight_length_mm']/2)**2))
            roots.append((z0,radius))
        fields=[]
        def restrict(base,volumes):
            f=gmsh.model.mesh.field.add('Restrict');gmsh.model.mesh.field.setNumber(f,'InField',base);gmsh.model.mesh.field.setNumbers(f,'VolumesList',volumes);gmsh.model.mesh.field.setNumber(f,'IncludeBoundary',1);fields.append(f)
        for roles,h in [([1,2,3,4,5,6,7,8],args.shell_h),([2],args.bond_h)]:
            f=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(f,'F',str(h));restrict(f,[t for t,p in pid.items() if p in roles])
        torus_faces=[t for _,t in gmsh.model.getEntities(2) if gmsh.model.getType(2,t)=='Torus']
        if torus_faces:
            f=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(f,'F','.25');g=gmsh.model.mesh.field.add('Restrict');gmsh.model.mesh.field.setNumber(g,'InField',f);gmsh.model.mesh.field.setNumbers(g,'SurfacesList',torus_faces);gmsh.model.mesh.field.setNumber(g,'IncludeBoundary',1);fields.append(g)
        port_faces=sorted({abs(face) for volume,part in pid.items() if part in [9,10,14,18] for dim,face in gmsh.model.getBoundary([(3,volume)],combined=False) if dim==2})
        port_face_h=float(__import__('os').environ.get('COPV_C04_PORT_FACE_H','1.0'))
        if not .5<=port_face_h<=2.:raise ValueError('Smooth port-face size outside declared range')
        if port_faces:
            f=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(f,'F',str(port_face_h));g=gmsh.model.mesh.field.add('Restrict');gmsh.model.mesh.field.setNumber(g,'InField',f);gmsh.model.mesh.field.setNumbers(g,'SurfacesList',port_faces);gmsh.model.mesh.field.setNumber(g,'IncludeBoundary',1);fields.append(g)
        for z0,radius in roots:
            f=gmsh.model.mesh.field.add('MathEval');expression=f'{args.root_h}+Min({args.h-args.root_h},1.8*Sqrt((Sqrt(x*x+y*y)-{radius})^2+(Abs(z)-{z0})^2))';gmsh.model.mesh.field.setString(f,'F',expression);restrict(f,[t for t,p in pid.items() if p in [1,2,3,4,5,6,7,8]])
            f=gmsh.model.mesh.field.add('MathEval');expression=f'{args.root_h}+Min({args.h-args.root_h},Max(0,Abs(Abs(z)-{z0})-{root_band})*{root_growth})';gmsh.model.mesh.field.setString(f,'F',expression);restrict(f,[t for t,p in pid.items() if p in [7,8]])
            if adjacent_band:
                adjacent_h=float(__import__('os').environ.get('COPV_C04_ADJACENT_ROOT_H',str(args.root_h)))
                f=gmsh.model.mesh.field.add('MathEval');expression=f'{adjacent_h}+Min({max(0,args.h-adjacent_h)},Max(0,Abs(Abs(z)-{z0})-{root_band})*{root_growth})';gmsh.model.mesh.field.setString(f,'F',expression);restrict(f,[t for t,p in pid.items() if p==6])
        # Resolve meridian curvature against the 0.15mm bondline, without applying
        # tiny port-fillet curvature to the entire vessel surface.
        chord_target=float(__import__('os').environ.get('COPV_C04_MAIN_CHORD_MM','.015'))
        if not .005<=chord_target<=.06:raise ValueError('Main chord sizing heuristic outside declared probe range')
        rz,normal=analytic_segments(params,(4,1500,16,1500,4));arc=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(rz,axis=0),axis=1))];unique=np.r_[True,np.diff(arc)>1e-10];rz=rz[unique];normal=normal[unique];arc=arc[unique];kappa=np.linalg.norm(np.gradient(normal,arc,axis=0),axis=1);meridian_sizes=[]
        for bucket in np.unique(np.floor(abs(rz[:,1]))):
            mask=(np.floor(abs(rz[:,1]))==bucket)&(abs(rz[:,1])>=params['straight_length_mm']/2-.5)
            if not mask.any():continue
            k=float(kappa[mask].max());size=min(args.shell_h,np.sqrt(8*chord_target/max(k,1e-12)));size=max(.6,size)
            if size<args.shell_h:meridian_sizes.append((float(np.mean(abs(rz[mask,1]))),float(size)))
        if meridian_sizes:
            expression=str(args.shell_h)
            for z0,size in meridian_sizes:expression=f'Min({expression},{size}+.25*Abs(Abs(z)-{z0}))'
            f=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(f,'F',expression);restrict(f,[t for t,p in pid.items() if p in [1,2,3,4,5,6,7,8]])
        thin_chord=__import__('os').environ.get('COPV_C04_THIN_CHORD_MM');thin_sizes=[]
        if thin_chord is not None:
            thin_chord=float(thin_chord)
            if not .005<=thin_chord<=chord_target:raise ValueError('Thin-layer chord sizing must be at least as fine as main sizing')
            for bucket in np.unique(np.floor(abs(rz[:,1]))):
                mask=(np.floor(abs(rz[:,1]))==bucket)&(abs(rz[:,1])>=params['straight_length_mm']/2-.5)
                if not mask.any():continue
                size=max(.6,min(args.shell_h,np.sqrt(8*thin_chord/max(float(kappa[mask].max()),1e-12))))
                if size<args.shell_h:thin_sizes.append((float(np.mean(abs(rz[mask,1]))),float(size)))
            expression=str(args.shell_h)
            for z0,size in thin_sizes:expression=f'Min({expression},{size}+.25*Abs(Abs(z)-{z0}))'
            f=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(f,'F',expression);restrict(f,[t for t,p in pid.items() if p in [1,2,3]])
        if args.adaptive_cloud:
            cloud=np.load(args.adaptive_cloud);points=np.asarray(cloud['xyz']);tags=[occ.addPoint(*map(float,q)) for q in points];occ.synchronize();f=gmsh.model.mesh.field.add('Distance');gmsh.model.mesh.field.setNumbers(f,'PointsList',tags);g=gmsh.model.mesh.field.add('Threshold')
            for key,value in [('InField',f),('SizeMin',1.),('SizeMax',args.h),('DistMin',.5),('DistMax',2.)]:gmsh.model.mesh.field.setNumber(g,key,value)
            fields.append(g)
        f=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(f,'FieldsList',fields);gmsh.model.mesh.field.setAsBackgroundMesh(f)
        circles=[]
        thin_curves={abs(curve) for volume,part in pid.items() if part in [1,2,3] for dim,face in gmsh.model.getBoundary([(3,volume)],combined=False) for dim,curve in gmsh.model.getBoundary([(2,abs(face))],combined=False)}
        for _,curve in gmsh.model.getEntities(1):
            if gmsh.model.getType(1,curve)!='Circle':continue
            length=occ.getMass(1,curve);box=gmsh.model.getBoundingBox(1,curve);z=(box[2]+box[5])/2
            if abs(box[5]-box[2])<1e-5 and length>4:
                radius=length/(2*np.pi);near=min(np.hypot(abs(z)-zz,radius-rr) for zz,rr in roots);local=args.root_h+min(args.h-args.root_h,1.8*near)
                if matched_root_band:
                    # A transfinite curve overrides volume sizing. Match the axial
                    # collar refinement on every circular material boundary too.
                    local=min(local,args.root_h+root_growth*max(0,min(abs(abs(z)-zz) for zz,rr in roots)-root_band))
                if radius>params['boss_bore_radius_mm']+3:local=min([local]+[size+.25*abs(abs(z)-z0) for z0,size in meridian_sizes])
                if curve in thin_curves and radius>params['boss_bore_radius_mm']+3:local=min([local]+[size+.25*abs(abs(z)-z0) for z0,size in thin_sizes])
                segments=max(args.circle_segments,int(np.ceil(length/min(args.h,local))));segments=8*int(np.ceil(segments/8));gmsh.model.mesh.setTransfiniteCurve(curve,segments+1);circles.append(dict(curve=curve,z=z,radius=radius,segments=segments))
            else:
                lo,hi=gmsh.model.getParametrizationBounds(1,curve);curvature=float(abs(gmsh.model.getCurvature(1,curve,[(float(lo[0])+float(hi[0]))/2])[0]));segments=max(4,int(np.ceil(length*curvature/(2*np.pi)*48)));gmsh.model.mesh.setTransfiniteCurve(curve,segments+1)
        torus_boundary_controls=[];torus_curves={abs(c) for face in torus_faces for dim,c in gmsh.model.getBoundary([(2,face)],combined=False)}
        for curve in sorted(torus_curves):
            length=occ.getMass(1,curve);lo,hi=gmsh.model.getParametrizationBounds(1,curve);k=float(np.max(abs(gmsh.model.getCurvature(1,curve,np.linspace(float(lo[0]),float(hi[0]),33).tolist()))));segments=max(4,int(np.ceil(length/.25)),int(np.ceil(length*k/(2*np.pi)*48)));found=next((c for c in circles if c['curve']==curve),None)
            if found is not None:segments=max(segments,found['segments']);found['segments']=segments
            gmsh.model.mesh.setTransfiniteCurve(curve,segments+1);torus_boundary_controls.append(dict(curve=curve,type=gmsh.model.getType(1,curve),length_mm=length,max_curvature_per_mm=k,segments=segments,target_max_edge_mm=.25))
        groups={}
        for circle in circles:groups.setdefault(round(abs(circle['z']),5),[]).append(circle)
        for group in groups.values():
            count=max(c['segments'] for c in group)
            for circle in group:circle['segments']=count;gmsh.model.mesh.setTransfiniteCurve(circle['curve'],count+1)
        report['repair_recipe']=dict(options=options,roots=roots,circles=circles,torus_surface_only_size_mm=.25,torus_faces=torus_faces,torus_boundary_controls=torus_boundary_controls,collar_frontal=__import__('os').environ.get('COPV_C04_COLLAR_FRONTAL')=='1',root_band_mm=root_band,matched_root_circle_sizing=matched_root_band,linear_optimizer=__import__('os').environ.get('COPV_C04_LINEAR_OPTIMIZER','none'),meridian_curvature_size_controls=meridian_sizes,nominal_meridian_chord_target_mm=.015,shell_h=args.shell_h,bond_h=args.bond_h,root_h=args.root_h,scope='Same BREP/materials. Explicit torus-boundary segmentation removes under-resolved transverse fillet chords; not a fully swept/multiblock layer topology.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),partition_sha256=hashlib.sha256(args.partition.read_bytes()).hexdigest())
        report['repair_recipe'].update(frontal_part_ids=sorted(frontal_parts),root_size_growth_per_mm=root_growth,adjacent_outer_ply_root_band=adjacent_band,adjacent_outer_ply_root_h_mm=adjacent_h if adjacent_band else None,nominal_meridian_chord_target_mm=chord_target,thin_layer_chord_target_mm=thin_chord,thin_layer_curvature_sizes=thin_sizes,smooth_port_face_h_mm=port_face_h,chord_scope='Linear boundary sizing heuristic only. Actual curved CAD-distance, Jacobian, SICN, material-volume and GPU resource criteria remain unchanged; torus .25mm boundary/surface controls retained.')
        if args.recipe=='baseline':
            gmsh.model.mesh.removeConstraints();gmsh.model.mesh.field.remove(f);uniform=gmsh.model.mesh.field.add('MathEval');gmsh.model.mesh.field.setString(uniform,'F',str(args.h));gmsh.model.mesh.field.setAsBackgroundMesh(uniform)
            for key,value in [('Mesh.MeshSizeMin',.5),('Mesh.MeshSizeFromCurvature',64),('Mesh.MeshSizeFromPoints',1)]:gmsh.option.setNumber(key,value);options[key]=value
            report['repair_recipe'].update(recipe='baseline_without_local_sizing',circles=[],scope='Original meshing settings with fixed seed/thread and guarded optimization; control candidate only')
        if args.linear_checkpoint:
            report['repair_recipe'].update(linear_checkpoint=str(args.linear_checkpoint.resolve()),linear_checkpoint_sha256=hashlib.sha256(args.linear_checkpoint.read_bytes()).hexdigest(),source='Resume the unchanged linear checkpoint in a NEW attempt; do not remesh or rewrite its parent.')
        surface_checkpoint=__import__('os').environ.get('COPV_C04_SURFACE_CHECKPOINT');surface_cap=int(__import__('os').environ.get('COPV_C04_MAX_SURFACE_NODES','450000'));report['repair_recipe']['max_surface_nodes']=surface_cap
        if surface_checkpoint:report['repair_recipe'].update(surface_checkpoint=surface_checkpoint,surface_checkpoint_sha256=hashlib.sha256(Path(surface_checkpoint).read_bytes()).hexdigest())
        write('recipe.json',report['repair_recipe'])
        if args.linear_checkpoint:
            gmsh.merge(str(args.linear_checkpoint));state('imported_retained_linear_checkpoint',**measure())
        else:
            if surface_checkpoint:gmsh.merge(surface_checkpoint);state('surface_checkpoint_imported',nodes=len(gmsh.model.mesh.getNodes()[0]))
            else:gmsh.model.mesh.generate(2);state('surface_generated',nodes=len(gmsh.model.mesh.getNodes()[0]))
            if len(gmsh.model.mesh.getNodes()[0])>surface_cap:raise RuntimeError('Surface node resource guard: retain attempt, do not start volume generation')
            if not surface_checkpoint:gmsh.write(str(out/'surface_checkpoint.msh'))
            if __import__('os').environ.get('COPV_C04_COLLAR_FRONTAL')=='1':
                gmsh.option.setNumber('Mesh.MeshOnlyVisible',1);gmsh.model.setVisibility(vols,0);gmsh.model.setVisibility([(3,t) for t,p in pid.items() if p not in frontal_parts],1);gmsh.model.mesh.generate(3);state('selected_Delaunay_generated',part_ids=sorted(set(pid.values())-frontal_parts));gmsh.model.setVisibility(vols,0);gmsh.model.setVisibility([(3,t) for t,p in pid.items() if p in frontal_parts],1);gmsh.option.setNumber('Mesh.Algorithm3D',4);gmsh.model.mesh.generate(3);gmsh.option.setNumber('Mesh.MeshOnlyVisible',0);gmsh.model.setVisibility(vols,1);assert all(sum(len(e) for e in gmsh.model.mesh.getElements(3,t)[1])>0 for t in pid),'Hybrid recipe left an empty material region';state('selected_Frontal_generated',part_ids=sorted(frontal_parts))
            else:gmsh.model.mesh.generate(3)
        log=gmsh.logger.get();write('gmsh_messages.json',log)
        if any(line.startswith('Error:') for line in log):raise RuntimeError('Gmsh reported a meshing error; inspect gmsh_messages.json and failed_mesh.msh')
        state('linear_generated',**measure())
        # Gmsh order conversion canonicalizes native tetra orientation; inspect it explicitly.
        gmsh.model.mesh.setOrder(2);state('straight_order2_canonicalized',**measure());gmsh.model.mesh.setOrder(1);baseline=measure();state('linear_checkpoint',**baseline)
        if baseline['nonpositive']:raise RuntimeError('Canonical linear mesh has nonpositive Jacobians')
        if baseline['cells']>args.max_cells:raise RuntimeError('Cell resource guard: '+str(baseline['cells']))
        if __import__('os').environ.get('COPV_C04_LINEAR_OPTIMIZER')=='Netgen':
            gmsh.write(str(out/'linear_before_Netgen.msh'));gmsh.model.mesh.optimize('Netgen',niter=3);optimized=measure();state('Netgen_linear_proposal',**optimized)
            if optimized['nonpositive'] or optimized['min_SICN']<baseline['min_SICN'] or optimized['below_005']>baseline['below_005']:raise RuntimeError('Netgen proposal failed positive-J/non-degradation guard; parent checkpoint retained')
            if optimized['cells']>args.max_cells:raise RuntimeError('Netgen proposal exceeded declared cell resource guard')
            baseline=optimized
        if not args.skip_relocate:
            tags,xyz=snapshot();gmsh.model.mesh.optimize('Relocate3D',niter=3);candidate=measure();newtags,newxyz=snapshot();assert np.array_equal(tags,newtags),'Relocate changed node IDs';moved=np.linalg.norm(newxyz-xyz,axis=1)>0;accepted=candidate['nonpositive']==0 and candidate['min_SICN']>=baseline['min_SICN'] and candidate['below_005']<=baseline['below_005'];np.savez_compressed(out/'relocation_proposal.npz',node_ids=tags[moved],before_xyz=xyz[moved],proposed_xyz=newxyz[moved]);state('relocation_proposal',accepted=accepted,moved_nodes=int(moved.sum()),**candidate)
            if not accepted:
                for node,position in zip(tags[moved],xyz[moved]):gmsh.model.mesh.setNode(int(node),position.tolist(),[])
                rt,rx=snapshot();assert np.array_equal(tags,rt) and np.array_equal(xyz,rx),'Node rollback did not restore the checkpoint';restored=measure();assert restored==baseline;state('relocation_rolled_back',**restored)
        gmsh.model.mesh.setOrder(2);final=measure();state('final_straight_TET10',**final)
        if final['nonpositive']:raise RuntimeError('Final straight TET10 has nonpositive Jacobians; no FE export')
        tags,xyz=snapshot();order=np.argsort(tags);sortedtags=tags[order];cells=[];ids=[];pids=[];qualities=[]
        for _,vol in vols:
            types,ee,cc=gmsh.model.mesh.getElements(3,vol)
            for typ,e,c in zip(types,ee,cc):
                assert typ==11;conn=np.asarray(c).reshape(-1,10);indices=order[np.searchsorted(sortedtags,conn)];assert np.array_equal(tags[indices],conn);cells.append(indices);ids.extend(e);pids.extend([pid[vol]]*len(e));qualities.extend(gmsh.model.mesh.getElementQualities(e,'minDetJac'))
        cells=np.concatenate(cells);used=np.unique(cells);remap=np.full(len(xyz),-1,np.int64);remap[used]=np.arange(len(used));cells=remap[cells];xyz=xyz[used];tags=tags[used];corners=np.unique(cells[:,:4]);inv=np.full(len(xyz),-1,np.int64);inv[corners]=np.arange(len(corners));edges=np.concatenate([cells[:,[i,j]] for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]);graph=sparse.coo_matrix((np.ones(len(edges)),(inv[edges[:,0]],inv[edges[:,1]])),shape=(len(corners),len(corners)));components,_=connected_components(graph,directed=False);assert components==1
        np.savez_compressed(out/'mesh.npz',xyz=xyz,cells=cells,node_ids=tags,cell_ids=np.asarray(ids),cell_part_ids=np.asarray(pids),min_det_jac=np.asarray(qualities));gmsh.write(str(out/'mesh.msh'));report.update(status='prepared',nodes=len(xyz),cells=len(cells),geometry_mapping='straight-sided TET10 coordinate storage (affine geometry); displacement degree selected by solver',min_det_jac_mm3=final['min_J'],connected_components=int(components),h_mm=args.h,meshing_total_seconds=time.perf_counter()-start,repair_stages=stages,production_eligible=False);write('mesh_review.json',report);state('exported_for_independent_quality_review',nodes=len(xyz),cells=len(cells))
    except BaseException as error:
        write('gmsh_messages.json',gmsh.logger.get());write('failure.json',dict(error=str(error),type=type(error).__name__,seconds=time.perf_counter()-start,stages=stages,error_entities=gmsh.model.mesh.getLastEntityError(),error_nodes=gmsh.model.mesh.getLastNodeError().tolist()));gmsh.write(str(out/'failed_mesh.msh'));raise
    finally:gmsh.finalize()
if __name__=='__main__':main()
