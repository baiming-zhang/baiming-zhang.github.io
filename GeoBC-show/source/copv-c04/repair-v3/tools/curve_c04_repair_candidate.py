"""Curve one versioned native mesh onto its unchanged CAD and integrate P2 volume."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,hashlib,json,os,shutil,time
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[_thread_key]='1'
import numpy as np
import gmsh
W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def optimize_sicn_patch(out):
    """Optimize only volume-owned nodes whose complete incident star is included."""
    from scipy.optimize import minimize
    import c04_shape_objective as objective
    validation=read(W/'data/copv/c04_repair_v3/sicn_objective_validation_v1/review.json');helper=Path(objective.__file__);assert validation['passed'] and hashlib.sha256(helper.read_bytes()).hexdigest()==validation['helper_sha256'];shutil.copyfile(helper,out/'c04_shape_objective.py')
    tags,coordinates,_=gmsh.model.mesh.getNodes();tags=np.asarray(tags);xyz=np.asarray(coordinates).reshape(-1,3);types,ids,conn=gmsh.model.mesh.getElements(3);assert list(types)==[11];ids=np.asarray(ids[0]);native=np.asarray(conn[0]).reshape(-1,10)
    if tags.max()<8*len(tags):
        inverse=np.full(int(tags.max())+1,-1,np.int64);inverse[tags]=np.arange(len(tags));cells=inverse[native]
        def dense(node_ids):return inverse[np.asarray(node_ids,dtype=np.int64)]
    else:
        order=np.argsort(tags);sorted_tags=tags[order]
        def dense(node_ids):return order[np.searchsorted(sorted_tags,np.asarray(node_ids))]
        cells=dense(native)
    assert np.all(cells>=0);quality=np.asarray(gmsh.model.mesh.getElementQualities(ids,'minSICN'));seed=np.flatnonzero(quality<.05);assert len(seed)>0,'No shape failures to optimize';target=float(os.environ.get('COPV_C04_SICN_OPT_TARGET','.065'));max_patch=int(os.environ.get('COPV_C04_SICN_MAX_PATCH_CELLS','50000'));rings=int(os.environ.get('COPV_C04_SICN_PATCH_RINGS','1'));mask=np.zeros(len(cells),bool);mask[seed]=True
    for _ in range(rings):
        hit=np.zeros(len(tags),bool);hit[np.unique(cells[mask])]=True;mask=np.any(hit[cells],axis=1)
    patch=np.flatnonzero(mask);assert len(patch)<=max_patch,'Local patch exceeds declared CPU guard';total_incidence=np.bincount(cells.ravel(),minlength=len(tags));patch_incidence=np.bincount(cells[patch].ravel(),minlength=len(tags));volume_owned=np.zeros(len(tags),bool)
    for _,vol in gmsh.model.getEntities(3):volume_owned[dense(gmsh.model.mesh.getNodes(3,vol,includeBoundary=False)[0])]=True
    free_global=np.flatnonzero(volume_owned&(patch_incidence>0)&(patch_incidence==total_incidence));assert len(free_global)>0,'No volume-owned free nodes in complete incident stars';used=np.unique(cells[patch]);remap=np.full(len(tags),-1,np.int64);remap[used]=np.arange(len(used));local_cells=remap[cells[patch]];free=remap[free_global];base=xyz[used].copy();assert len(free)<=20000,'Too many local coordinate variables'
    gp=np.asarray(gmsh.model.mesh.getIntegrationPoints(11,'Gauss4')[0]).reshape(-1,3);lattice=np.array([[i/4,j/4,k/4] for i in range(5) for j in range(5-i) for k in range(5-i-j)]);points=np.unique(np.round(np.vstack([gp,lattice]),14),axis=0);dN=objective.reference_derivatives(points);corners=base[local_cells[:,:4]];ref=np.linalg.det(np.swapaxes(corners[:,1:]-corners[:,:1],1,2)@objective.IDEAL_INV);assert ref.min()>0;length=np.min(np.stack([np.linalg.norm(corners[:,b]-corners[:,a],axis=1) for a,b in objective.EDGES],axis=1),axis=1);scale=np.full(len(base),np.inf);np.minimum.at(scale,local_cells.ravel(),np.repeat(length,10));scale=scale[free];assert scale.min()>0;iteration=0;calls=0;last_loss=None;start=time.perf_counter();movement_weight=float(os.environ.get('COPV_C04_SICN_MOVE_WEIGHT','1e-5'))
    def fun(delta):
        nonlocal calls,last_loss
        calls+=1;loss,g=objective.objective_gradient(delta,base,local_cells,free,dN,target,ref,scale,movement_weight=movement_weight);assert np.isfinite(loss) and np.isfinite(g).all();last_loss=loss;return loss,g
    native_gate_met=False;native_check_every=int(os.environ.get('COPV_C04_SICN_NATIVE_CHECK_EVERY','0'))
    def callback(delta):
        nonlocal iteration
        nonlocal native_gate_met
        iteration+=1
        if iteration%5==0:write(out/'SICN_optimizer_progress.json',dict(iteration=iteration,evaluations=calls,last_objective=last_loss,seconds=time.perf_counter()-start));print(json.dumps(dict(SICN_iteration=iteration,objective=last_loss,seconds=time.perf_counter()-start)),flush=True)
        if native_check_every and iteration%native_check_every==0:
            positions=base[free]+np.asarray(delta).reshape(-1,3)
            for node,position in zip(tags[free_global],positions):gmsh.model.mesh.setNode(int(node),position.tolist(),[])
            patch_q=np.asarray(gmsh.model.mesh.getElementQualities(ids[patch],'minSICN'));patch_j=np.asarray(gmsh.model.mesh.getElementQualities(ids[patch],'minDetJac'));native_gate_met=bool(np.isfinite(patch_q).all() and np.isfinite(patch_j).all() and patch_q.min()>=.0505 and patch_j.min()>0);write(out/'SICN_native_patch_check.json',dict(iteration=iteration,min_SICN=float(patch_q.min()),min_J=float(patch_j.min()),native_patch_gate_met=native_gate_met,scope='Outside patch unchanged; mandatory final full-mesh check still follows.'))
            if native_gate_met:raise StopIteration
    bound=np.repeat(.5*scale,3);initial=np.zeros(3*len(free));initial_loss=fun(initial)[0];write(out/'SICN_patch_scope.json',dict(full_cells=len(cells),seed_cells=len(seed),patch_cells=len(patch),free_volume_nodes=len(free),points_per_cell=len(points),target_sampled_SICN=target,original_acceptance_min_SICN=.05,all_CAD_and_incomplete_star_nodes_fixed=True,helper_sha256=validation['helper_sha256'],scope='Sampled quality objective only; full native Gmsh SICN/Jacobian and original CAD/volume/load/resource gates decide acceptance.'))
    solver_options={'maxiter':int(os.environ.get('COPV_C04_SICN_MAX_ITER','120')),'ftol':float(os.environ.get('COPV_C04_SICN_FTOL','1e-13')),'gtol':float(os.environ.get('COPV_C04_SICN_GTOL','1e-9')),'maxls':20,'maxcor':10};write(out/'SICN_optimizer_settings.json',dict(options=solver_options,movement_weight=movement_weight,target_sampled_SICN=target,original_acceptance_min_SICN=.05))
    result=minimize(fun,initial,jac=True,method='L-BFGS-B',bounds=list(zip(-bound,bound)),callback=callback,options=solver_options);assert np.isfinite(result.x).all();positions=base[free]+result.x.reshape(-1,3);np.savez_compressed(out/'SICN_patch_proposal.npz',node_ids=tags[free_global],before_xyz=xyz[free_global],after_xyz=positions,node_step_scales_mm=scale,seed_cell_ids=ids[seed],patch_cell_ids=ids[patch])
    for node,position in zip(tags[free_global],positions):gmsh.model.mesh.setNode(int(node),position.tolist(),[])
    rt,rx,_=gmsh.model.mesh.getNodes();assert np.array_equal(tags,rt);keep=np.ones(len(tags),bool);keep[free_global]=False;assert np.array_equal(xyz[keep],np.asarray(rx).reshape(-1,3)[keep]),'A fixed node changed';types,ri,rc=gmsh.model.mesh.getElements(3);assert list(types)==[11] and np.array_equal(ids,ri[0]) and np.array_equal(native.ravel(),rc[0]),'Native topology changed';write(out/'SICN_optimizer_review.json',dict(optimizer_success=bool(result.success),message=str(result.message),native_patch_gate_met=native_gate_met,iterations=int(result.nit),evaluations=calls,initial_objective=float(initial_loss),final_objective=float(result.fun),full_native_topology_preserved=True,fixed_nodes_unchanged=True,seconds=time.perf_counter()-start))

def optimize_mesh_roi(method,selected,partition,full_snapshot,out):
    """Separate model enforces scope; never rely on the unimplemented dimTags API."""
    full_name=gmsh.model.getCurrent();tags,coords,_=gmsh.model.mesh.getNodes();tags=np.asarray(tags);xyz=np.asarray(coords).reshape(-1,3);sort=np.argsort(tags);sorted_tags=tags[sort];expected=[];allowed=[]
    for _,vol in selected:
        types,ids,_=gmsh.model.mesh.getElements(3,vol);assert list(types)==[11];expected.extend(ids[0]);allowed.extend(gmsh.model.mesh.getNodes(3,vol,includeBoundary=False)[0])
    expected=np.sort(np.asarray(expected));allowed=np.asarray(allowed);original_types,original_ids,original_conn=gmsh.model.mesh.getElements(3);assert list(original_types)==[11];original_ids=np.asarray(original_ids[0]);original_conn=np.asarray(original_conn[0]);roi_name='c04_ROI_optimizer';applied=False
    try:
        gmsh.model.add(roi_name);gmsh.model.occ.importShapes(str(partition));gmsh.model.occ.synchronize();gmsh.merge(str(full_snapshot));selected_tags={t for _,t in selected};remove=[(3,t) for _,t in gmsh.model.getEntities(3) if t not in selected_tags];gmsh.model.occ.remove(remove,recursive=False);gmsh.model.occ.synchronize();assert set(gmsh.model.getEntities(3))==set(selected),'CAD ROI entity set differs';types,ids,_=gmsh.model.mesh.getElements(3);assert list(types)==[11] and np.array_equal(np.sort(ids[0]),expected),'ROI contains missing or unexpected volume cells'
        write(out/'ROI_scope.json',dict(full_cells=len(original_ids),ROI_cells=len(expected),selected_entities=selected,method='Separate model with unselected CAD volumes removed, keeping shared boundary geometry; no unsupported dimTags filtering. Partial ROI is never an FE mesh.'))
        print(json.dumps(dict(stage='ROI_scope_enforced',full_cells=len(original_ids),ROI_cells=len(expected),selected_entities=selected)),flush=True);gmsh.model.mesh.optimize(method)
        changed_tags,changed_coords,_=gmsh.model.mesh.getNodes();changed_tags=np.asarray(changed_tags);changed_xyz=np.asarray(changed_coords).reshape(-1,3);index=np.searchsorted(sorted_tags,changed_tags);assert np.all(index<len(sorted_tags)) and np.array_equal(sorted_tags[index],changed_tags),'Optimizer created or renumbered nodes';base=xyz[sort[index]];delta=np.linalg.norm(changed_xyz-base,axis=1);permitted=np.isin(changed_tags,allowed);outside=float(delta[~permitted].max()) if (~permitted).any() else 0.;assert outside<=1e-12,'Fixed boundary or non-selected node moved';moved=permitted&(delta>1e-13);update_tags=changed_tags[moved];updates=changed_xyz[moved];np.savez_compressed(out/'ROI_node_proposal.npz',node_ids=update_tags,before_xyz=base[moved],after_xyz=updates)
        gmsh.model.remove();gmsh.model.setCurrent(full_name);rt,rx,_=gmsh.model.mesh.getNodes();assert np.array_equal(tags,rt) and np.array_equal(xyz,np.asarray(rx).reshape(-1,3)),'Original full model changed during ROI optimization';types,ids,conn=gmsh.model.mesh.getElements(3);assert list(types)==[11] and np.array_equal(original_ids,ids[0]) and np.array_equal(original_conn,conn[0]),'Full native topology changed'
        for node,position in zip(update_tags,updates):gmsh.model.mesh.setNode(int(node),position.tolist(),[])
        applied=True;write(out/'ROI_restore.json',dict(full_native_topology_preserved=True,full_cells=len(original_ids),ROI_cells=len(expected),modified_volume_nodes=len(update_tags),max_nonselected_displacement_mm=outside,boundaries_fixed=True,proposal_applied=True))
    finally:
        if gmsh.model.getCurrent()==roi_name:gmsh.model.remove()
        gmsh.model.setCurrent(full_name)
    assert applied
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--safe-local',action='store_true');ap.add_argument('--high-order',action='store_true');ap.add_argument('--shape-optimize',action='store_true');ap.add_argument('--direct-only',action='store_true');ap.add_argument('--high-order-method',choices=['HighOrder','HighOrderElastic'],default='HighOrder');args=ap.parse_args();parent=args.parent.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'curver_executed.py');state=read(parent/'state.json');review=read(parent/'mesh/mesh_review.json');partition=Path(state['partition_file']) if 'partition_file' in state else Path(state['command'][state['command'].index('--partition')+1]);start=time.perf_counter();stages=[]
    def measure(label):
        types,ids,_=gmsh.model.mesh.getElements(3);assert list(types)==[11];ee=ids[0];j=np.asarray(gmsh.model.mesh.getElementQualities(ee,'minDetJac'));q=np.asarray(gmsh.model.mesh.getElementQualities(ee,'minSICN'));row=dict(stage=label,seconds=time.perf_counter()-start,cells=len(ee),min_J=float(j.min()),nonpositive=int((j<=0).sum()),min_SICN=float(q.min()),below_005=int((q<.05).sum()));stages.append(row);write(out/'stages.json',stages);print(json.dumps(row),flush=True);return row
    gmsh.initialize();gmsh.model.add('c04_full');gmsh.option.setNumber('General.NumThreads',1);gmsh.option.setNumber('Mesh.Binary',1);gmsh.option.setNumber('Mesh.SaveAll',1)
    try:
        gmsh.model.occ.importShapes(str(partition));gmsh.model.occ.synchronize();types_before={tag:gmsh.model.getType(2,tag) for _,tag in gmsh.model.getEntities(2)};gmsh.merge(str(parent/'mesh/mesh.msh'));assert all(gmsh.model.getType(2,tag)==typ for tag,typ in types_before.items()),'MSH merge replaced CAD surfaces';measure('parent_affine_imported');parent_ids=gmsh.model.mesh.getElements(3)[1][0];parent_q=np.asarray(gmsh.model.mesh.getElementQualities(parent_ids,'minSICN'));parent_sort=np.argsort(parent_ids)
        gmsh.option.setNumber('Mesh.SecondOrderLinear',0)
        if os.environ.get('COPV_C04_PROJECTED_CHECKPOINT'):
            checkpoint=Path(os.environ['COPV_C04_PROJECTED_CHECKPOINT']);gmsh.model.mesh.clear();gmsh.merge(str(checkpoint));write(out/'projected_checkpoint.json',dict(path=str(checkpoint),sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),reason='Reuse retained full CAD projection in a NEW attempt; no repeat order conversion.'))
        else:gmsh.model.mesh.setOrder(1);gmsh.model.mesh.setOrder(2)
        projected=measure('CAD_projected_order2')
        if args.shape_optimize:
            untangle=os.environ.get('COPV_C04_UNTANGLE_FULL_CAD')=='1'
            assert not args.high_order and not args.safe_local and (projected['nonpositive']==0 or untangle)
            if untangle:
                proof=read(W/'data/copv/c04_repair_v3/sicn_negative_J_validation_v1/review.json');assert proof['passed'] and proof['helper_sha256']==hashlib.sha256((W/'notes/c04_shape_objective.py').read_bytes()).hexdigest()
                assert not os.environ.get('COPV_C04_PROJECTED_CHECKPOINT'),'Untangling must begin from a fresh full CAD projection, never a retreated checkpoint'
                gmsh.write(str(out/'full_CAD_before_untangling.msh'));write(out/'full_CAD_untangling_policy.json',dict(starting_native_nonpositive=projected['nonpositive'],starting_native_low_SICN=projected['below_005'],all_CAD_nodes_fixed=True,only_complete_star_volume_nodes_may_move=True,negative_J_derivative_validation_sha256=hashlib.sha256((W/'data/copv/c04_repair_v3/sicn_negative_J_validation_v1/review.json').read_bytes()).hexdigest(),scope='Intermediate invalid mesh repair only. No FE output or acceptance until full native positiveJ/SICN and independent CAD/volume/load/resource screens pass. No guarded boundary retreat.'))
            before_shape=projected.copy();optimize_sicn_patch(out);projected=measure('after_sampled_SICN_patch_optimization')
            if projected['nonpositive'] or projected['min_SICN']<before_shape['min_SICN']-1e-12 or projected['below_005']>before_shape['below_005']:
                gmsh.write(str(out/'rejected_curved_mesh.msh'));write(out/'result.json',dict(status='SICN_patch_proposal_rejected_no_FE',reason='Full native positive-J/shape non-degradation check failed',stages=stages));return
        if args.high_order:
            before_high_order=projected.copy();mesh_roi=os.environ.get('COPV_C04_HIGH_ORDER_MESH_ROI')=='1'
            gmsh.write(str(out/'CAD_projected_before_optimization.msh'));options={'Mesh.Renumber':0,'Mesh.HighOrderFixBoundaryNodes':0,'Mesh.HighOrderDistCAD':1,'Mesh.HighOrderNumLayers':2,'Mesh.HighOrderIterMax':50,'Mesh.HighOrderPassMax':5,'Mesh.HighOrderThresholdMin':.1,'Mesh.HighOrderThresholdMax':2.}
            if args.high_order_method=='HighOrderElastic':options['Mesh.HighOrderFixBoundaryNodes']=1
            if 'COPV_C04_HIGH_ORDER_FIX_BOUNDARY' in os.environ:options['Mesh.HighOrderFixBoundaryNodes']=int(os.environ['COPV_C04_HIGH_ORDER_FIX_BOUNDARY'])
            if 'COPV_C04_HIGH_ORDER_INTERNAL_MIN' in os.environ:options['Mesh.HighOrderThresholdMin']=float(os.environ['COPV_C04_HIGH_ORDER_INTERNAL_MIN'])
            if 'COPV_C04_HIGH_ORDER_DIST_CAD' in os.environ:options['Mesh.HighOrderDistCAD']=int(os.environ['COPV_C04_HIGH_ORDER_DIST_CAD'])
            requested_parts=[int(v) for v in os.environ.get('COPV_C04_HIGH_ORDER_PARTS','').split(',') if v];selected_entities=[]
            if requested_parts:
                present={int(p) for _,p in gmsh.model.getPhysicalGroups(3)};assert set(requested_parts)<=present,'Unknown requested material part';selected_entities=sorted({(3,int(v)) for p in requested_parts for v in gmsh.model.getEntitiesForPhysicalGroup(3,p)})
                if not mesh_roi:raise RuntimeError('Installed Gmsh does not interface optimizer dimTags filtering; use a separately enforced mesh ROI in a NEW attempt')
            write(out/'high_order_options.json',dict(method=args.high_order_method,options=options,selected_material_parts=requested_parts,selected_entities=selected_entities,documentation='https://gmsh.info/doc/texinfo/#Mesh-options',scope='Explicit versioned candidate. Optimizer internal parameters are not SICN acceptance thresholds. Actual CAD/volume/quality checks still required; no release claim.'))
            for key,value in options.items():gmsh.option.setNumber(key,value)
            if mesh_roi:
                assert selected_entities and projected['nonpositive']==0 and options['Mesh.HighOrderFixBoundaryNodes']==1,'ROI requires a valid initial map and fixed boundary';optimize_mesh_roi(args.high_order_method,selected_entities,partition,out/'CAD_projected_before_optimization.msh',out)
            else:gmsh.model.mesh.optimize(args.high_order_method)
            projected=measure('after_'+args.high_order_method+'_optimization')
            if mesh_roi and (projected['nonpositive'] or projected['min_SICN']<before_high_order['min_SICN']-1e-12 or projected['below_005']>before_high_order['below_005']):
                gmsh.write(str(out/'rejected_curved_mesh.msh'));write(out/'result.json',dict(status='ROI_proposal_rejected_no_FE',reason='Full-mesh positive-J/shape non-degradation guard failed',stages=stages));return
        if args.safe_local:
            tags0,coordinates,_=gmsh.model.mesh.getNodes();tags0=np.asarray(tags0);target=np.asarray(coordinates).reshape(-1,3);sort0=np.argsort(tags0);sorted0=tags0[sort0];types0,ids0,cc0=gmsh.model.mesh.getElements(3);assert list(types0)==[11];ids0=ids0[0];native0=np.asarray(cc0[0]).reshape(-1,10);conn0=sort0[np.searchsorted(sorted0,native0)];base=target.copy();assigned=np.zeros(len(base),bool)
            for mid,(a,b) in enumerate([(0,1),(1,2),(0,2),(0,3),(2,3),(1,3)],4):
                indices=conn0[:,mid];values=.5*(target[conn0[:,a]]+target[conn0[:,b]]);unique,first=np.unique(indices,return_index=True);base[unique]=values[first];assigned[unique]=True
            delta=target-base;alpha=np.ones(len(base));baseline_j=np.linalg.det(np.swapaxes(base[conn0[:,1:4]]-base[conn0[:,:1]],1,2));assert baseline_j.min()>0;parent_indices=parent_sort[np.searchsorted(np.asarray(parent_ids)[parent_sort],ids0)];assert np.array_equal(np.asarray(parent_ids)[parent_indices],ids0);quality_floor=.98*np.minimum(.05,parent_q[parent_indices])
            det_only=os.environ.get('COPV_C04_CURVE_DET_ONLY')=='1'
            write(out/'curving_guard_policy.json',dict(positive_J_guard=True,per_cell_parent_SICN_floor_enabled=not det_only,original_production_min_SICN=.05,production_eligible=False,scope='Local validity candidate only. Disabling the parent-relative floor does not relax the original production threshold. Independent CAD distance and volume checks required.'))
            if det_only:quality_floor=np.full(len(ids0),-np.inf)
            for iteration in range(36):
                j=np.asarray(gmsh.model.mesh.getElementQualities(ids0,'minDetJac'));q=np.asarray(gmsh.model.mesh.getElementQualities(ids0,'minSICN'));bad=(~np.isfinite(j))|(j<.05*baseline_j)|(~np.isfinite(q))|(q<quality_floor)
                if not bad.any():break
                affected=np.unique(conn0[bad,4:].ravel());alpha[affected]=alpha[affected]*.5 if iteration<30 else 0.;positions=base[affected]+alpha[affected,None]*delta[affected]
                for node,position in zip(tags0[affected],positions):gmsh.model.mesh.setNode(int(node),position.tolist(),[])
                write(out/'safe_curving_progress.json',dict(iteration=iteration,unsafe_cells=int(bad.sum()),affected_midnodes=len(affected),seconds=time.perf_counter()-start))
                print(json.dumps(dict(safe_curving_iteration=iteration,unsafe_cells=int(bad.sum()),affected_midnodes=len(affected))),flush=True)
            final_j=np.asarray(gmsh.model.mesh.getElementQualities(ids0,'minDetJac'));final_q=np.asarray(gmsh.model.mesh.getElementQualities(ids0,'minSICN'));assert np.isfinite(final_j).all() and final_j.min()>0,'Guarded local curving could not certify positive J';assert np.isfinite(final_q).all() and np.all(final_q>=quality_floor-1e-12),'Guarded curving did not preserve its per-cell shape-quality floor'
            moved=np.linalg.norm(delta,axis=1)>1e-12;np.savez_compressed(out/'projection_fractions.npz',node_ids=tags0[moved],fraction=alpha[moved],affine_xyz=base[moved],CAD_target_xyz=target[moved]);write(out/'projection_summary.json',dict(method='Shared-midnode local line search from affine coordinates toward CAD projection; corners fixed; keep all element minDetJac >=5% affine baseline where possible, and strictly positive globally.',CAD_target_nodes=int(moved.sum()),fully_projected=int((moved&(alpha==1)).sum()),partially_projected=int((moved&(alpha>0)&(alpha<1)).sum()),affine_retained=int((moved&(alpha==0)).sum()),max_target_distance_mm=float(np.linalg.norm((1-alpha[moved,None])*delta[moved],axis=1).max()),geometry_scope='Quadratic approximation; partially projected nodes may remain off CAD. Must pass independent per-region volume and boundary/load checks before any FE field.'))
            projected=measure('after_guarded_local_curving')
        elif projected['nonpositive'] and not args.high_order and not args.shape_optimize and not args.direct_only:
            gmsh.model.mesh.optimize('HighOrderFastCurving');projected=measure('after_HighOrderFastCurving')
        if projected['nonpositive']:
            gmsh.write(str(out/'rejected_curved_mesh.msh'));write(out/'result.json',dict(status='nonpositive_curved_J_retained_no_FE',parent=str(parent.relative_to(W)),stages=stages));return
        tags,coords,_=gmsh.model.mesh.getNodes();tags=np.asarray(tags);xyz=np.asarray(coords).reshape(-1,3);order=np.argsort(tags);sortedtags=tags[order];pid={int(t):int(p) for _,p in gmsh.model.getPhysicalGroups(3) for t in gmsh.model.getEntitiesForPhysicalGroup(3,p)};conn=[];ids=[];pids=[];minj=[]
        for _,vol in gmsh.model.getEntities(3):
            typ,ee,cc=gmsh.model.mesh.getElements(3,vol)
            for t,e,c in zip(typ,ee,cc):
                assert t==11;native=np.asarray(c).reshape(-1,10);idx=order[np.searchsorted(sortedtags,native)];assert np.array_equal(tags[idx],native);conn.append(idx);ids.extend(e);pids.extend([pid[vol]]*len(e));minj.extend(gmsh.model.mesh.getElementQualities(e,'minDetJac'))
        cells=np.concatenate(conn);pids=np.asarray(pids);used=np.unique(cells);remap=np.full(len(tags),-1,np.int64);remap[used]=np.arange(len(used));cells=remap[cells];xyz=xyz[used];tags=tags[used];mesh=out/'mesh';mesh.mkdir();np.savez_compressed(mesh/'mesh.npz',xyz=xyz,cells=cells,node_ids=tags,cell_ids=np.asarray(ids),cell_part_ids=pids,min_det_jac=np.asarray(minj));gmsh.write(str(mesh/'mesh.msh'))
        layout=read(W/'data/copv/c04_research_v1/analysis/GP_layout.json');points=np.asarray(layout['points']);weights=np.asarray(layout['reference_weights']);lam=np.c_[1-points.sum(1),points];grad=np.array([[-1.,-1.,-1.],[1,0,0],[0,1,0],[0,0,1]]);dN=np.zeros((len(points),10,3));dN[:,:4]=(4*lam[:,:,None]-1)*grad[None,:,:]
        for index,(a,b) in enumerate([(0,1),(1,2),(0,2),(0,3),(2,3),(1,3)],4):dN[:,index]=4*(lam[:,a,None]*grad[b]+lam[:,b,None]*grad[a])
        volumes=np.empty(len(cells));min_gp=float('inf')
        for first in range(0,len(cells),20000):
            nodes=xyz[cells[first:first+20000]];jac=np.einsum('eni,qnj->eqij',nodes,dN);det=np.linalg.det(jac);min_gp=min(min_gp,float(det.min()));volumes[first:first+len(nodes)]=det@weights
        assert min_gp>0 and np.isfinite(volumes).all()
        expected={}
        for record in review['volume_records']:expected[record['part_id']]=expected.get(record['part_id'],0.)+record['volume_mm3']
        if 'quality_reference_part_volumes_mm3' in review:expected={int(k):v for k,v in review['quality_reference_part_volumes_mm3'].items()}
        parts=[dict(part_id=p,CAD_volume_mm3=expected[p],mesh_curved_volume_mm3=float(volumes[pids==p].sum()),relative_error=abs(float(volumes[pids==p].sum())/expected[p]-1)) for p in sorted(expected)];np.savez_compressed(out/'curved_cell_volumes.npz',cell_ids=np.asarray(ids),volume_mm3=volumes)
        review.update(geometry_mapping='curved quadratic TET10 isoparametric CAD geometry; requires curved Tri6 loads and GP-Jacobian volume integration',nodes=len(xyz),cells=len(cells),min_det_jac_mm3=projected['min_J'],meshing_total_seconds=time.perf_counter()-start,parent_attempt=str(parent.relative_to(W)),curving_stages=stages,production_eligible=False);write(mesh/'mesh_review.json',review);write(out/'result.json',dict(status='curved_positive_J_mesh_preflight_only',case_id=state['case_id'],parent_attempt=str(parent.relative_to(W)),stages=stages,min_GP_J=min_gp,volume_method='Original Gmsh TET10 coordinate basis, degree4 reference tetra quadrature; detJ varies by GP',parts=parts,max_part_volume_error=max(p['relative_error'] for p in parts),new_FE_labels=0,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.perf_counter()-start));print(json.dumps(read(out/'result.json')),flush=True)
    except BaseException as error:write(out/'failure.json',dict(error=str(error),seconds=time.perf_counter()-start,stages=stages));raise
    finally:gmsh.finalize()
if __name__=='__main__':main()
