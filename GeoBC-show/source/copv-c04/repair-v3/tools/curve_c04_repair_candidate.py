"""Curve one versioned native mesh onto its unchanged CAD and integrate P2 volume."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,hashlib,json,os,shutil,time
import numpy as np
import gmsh
W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--safe-local',action='store_true');ap.add_argument('--high-order',action='store_true');ap.add_argument('--direct-only',action='store_true');ap.add_argument('--high-order-method',choices=['HighOrder','HighOrderElastic'],default='HighOrder');args=ap.parse_args();parent=args.parent.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'curver_executed.py');state=read(parent/'state.json');review=read(parent/'mesh/mesh_review.json');partition=Path(state['partition_file']) if 'partition_file' in state else Path(state['command'][state['command'].index('--partition')+1]);start=time.perf_counter();stages=[]
    def measure(label):
        types,ids,_=gmsh.model.mesh.getElements(3);assert list(types)==[11];ee=ids[0];j=np.asarray(gmsh.model.mesh.getElementQualities(ee,'minDetJac'));q=np.asarray(gmsh.model.mesh.getElementQualities(ee,'minSICN'));row=dict(stage=label,seconds=time.perf_counter()-start,cells=len(ee),min_J=float(j.min()),nonpositive=int((j<=0).sum()),min_SICN=float(q.min()),below_005=int((q<.05).sum()));stages.append(row);write(out/'stages.json',stages);print(json.dumps(row),flush=True);return row
    gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1);gmsh.option.setNumber('Mesh.Binary',1);gmsh.option.setNumber('Mesh.SaveAll',1)
    try:
        gmsh.model.occ.importShapes(str(partition));gmsh.model.occ.synchronize();types_before={tag:gmsh.model.getType(2,tag) for _,tag in gmsh.model.getEntities(2)};gmsh.merge(str(parent/'mesh/mesh.msh'));assert all(gmsh.model.getType(2,tag)==typ for tag,typ in types_before.items()),'MSH merge replaced CAD surfaces';measure('parent_affine_imported');parent_ids=gmsh.model.mesh.getElements(3)[1][0];parent_q=np.asarray(gmsh.model.mesh.getElementQualities(parent_ids,'minSICN'));parent_sort=np.argsort(parent_ids)
        gmsh.option.setNumber('Mesh.SecondOrderLinear',0)
        if os.environ.get('COPV_C04_PROJECTED_CHECKPOINT'):
            checkpoint=Path(os.environ['COPV_C04_PROJECTED_CHECKPOINT']);gmsh.model.mesh.clear();gmsh.merge(str(checkpoint));write(out/'projected_checkpoint.json',dict(path=str(checkpoint),sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),reason='Reuse retained full CAD projection in a NEW attempt; no repeat order conversion.'))
        else:gmsh.model.mesh.setOrder(1);gmsh.model.mesh.setOrder(2)
        projected=measure('CAD_projected_order2')
        if args.high_order:
            gmsh.write(str(out/'CAD_projected_before_optimization.msh'));options={'Mesh.HighOrderFixBoundaryNodes':0,'Mesh.HighOrderDistCAD':1,'Mesh.HighOrderNumLayers':2,'Mesh.HighOrderIterMax':50,'Mesh.HighOrderPassMax':5,'Mesh.HighOrderThresholdMin':.1,'Mesh.HighOrderThresholdMax':2.}
            if args.high_order_method=='HighOrderElastic':options['Mesh.HighOrderFixBoundaryNodes']=1
            write(out/'high_order_options.json',dict(method=args.high_order_method,options=options,documentation='https://gmsh.info/doc/texinfo/#Mesh-options',scope='Explicit versioned candidate. Actual boundary-to-CAD check still required; no release claim.'))
            for key,value in options.items():gmsh.option.setNumber(key,value)
            gmsh.model.mesh.optimize(args.high_order_method);projected=measure('after_'+args.high_order_method+'_optimization')
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
        elif projected['nonpositive'] and not args.high_order and not args.direct_only:
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
