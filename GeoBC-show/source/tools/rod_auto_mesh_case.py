"""One-command, bounded automatic curved rod meshing with explicit quality gates."""
import argparse,hashlib,json,os,re,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
from tools.run_rod_v2_dataset import resources

REPO=Path(__file__).resolve().parents[1];WS=REPO.parents[1]
ROOT=REPO/'results/rod_auto_mesh_v1'
POLICIES=[dict(surface_algorithm=6,algorithm3d=1,size_factor=1.,curvature=8),dict(surface_algorithm=2,algorithm3d=1,size_factor=.85,curvature=8),dict(surface_algorithm=6,algorithm3d=1,size_factor=.70,curvature=12)]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def guard():
    ram,disk=resources()
    if ram<8*1024**3 or disk<20*1024**3:raise RuntimeError('8GiB RAM /20GiB disk resource reserve')
def inputs(case):
    assert re.fullmatch(r'[A-Za-z0-9_-]{1,80}',case),'Safe geometry name required'
    request=ROOT/'requests'/(case+'.json')
    if request.exists():
        request=read(request);record=dict(case_id=case,parameters=request['parameters']);meta=dict(source_step=request['step'],source_step_sha256=request['STEP_sha256'],h_mm=request['h_mm']);step=Path(meta['source_step']);assert sha(step)==meta['source_step_sha256'];return record,meta,step
    record=next(r for r in read(WS/'data/rod_mixed_1000_v1/selection_manifest.json')['new_geometry_records'] if r['case_id']==case)
    meta=read(WS/'data/rod_v2_sixload_v1/cases'/case/'h3/mesh.json');step=Path(meta['source_step']);assert sha(step)==meta['source_step_sha256']
    return record,meta,step
def worker(case,attempt):
    import gmsh
    from tools.repair_curved_edge_nodes import repair_edges
    from geobc.rod_v2_probes import shape
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    guard();record,meta,step=inputs(case);p=record['parameters'];policy=POLICIES[attempt];out=ROOT/'cases'/case/f'attempt_{attempt+1}';out.mkdir(parents=True,exist_ok=True);start=time.time()
    def status(phase,**kw):write(out/'status.json',dict(phase=phase,pid=os.getpid(),elapsed_seconds=time.time()-start,updated_unix=time.time(),no_FE_started=True,**kw))
    def qualities(typ):
        ids,conn=gmsh.model.mesh.getElementsByType(typ)
        return np.asarray(ids),np.asarray(conn),np.asarray(gmsh.model.mesh.getElementQualities(ids,'minDetJac')),np.asarray(gmsh.model.mesh.getElementQualities(ids,'minSICN'))
    gmsh.initialize()
    try:
        status('CAD_preflight');gmsh.option.setNumber('General.NumThreads',2);gmsh.model.add(case);gmsh.model.occ.importShapes(str(step));gmsh.model.occ.synchronize();volumes=gmsh.model.getEntities(3);assert len(volumes)==1,'CAD must be one solid'
        cad_volume=gmsh.model.occ.getMass(3,volumes[0][1]);assert cad_volume>0
        factor=policy['size_factor'];hmax=meta['h_mm']*factor;hmin=min(.45,p['tf']/6,p['re']*1.5)*factor;local=min(p['tf']/2,hmax)*factor
        assert hmin>0 and local>hmin
        for key,value in {'Mesh.Algorithm':policy['surface_algorithm'],'Mesh.Algorithm3D':policy['algorithm3d'],'Mesh.MeshSizeMin':hmin,'Mesh.MeshSizeMax':hmax,'Mesh.MeshSizeFromCurvature':policy['curvature'],'Mesh.MeshSizeFromPoints':0,'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.MinimumCircleNodes':24,'Mesh.RandomSeed':2026,'Mesh.Smoothing':10,'Mesh.OptimizeThreshold':.3}.items():gmsh.option.setNumber(key,value)
        low=(p['rs']+p['ds'])*.75;high=p['L']-(p['rb']+p['db'])*.65;transition=15.*factor
        def size(dim,tag,x,y,z,lc):
            distance=max(low-x,0,x-high);t=min(distance/transition,1);smooth=t*t*(3-2*t);return min(lc,local+(hmax-local)*smooth)
        gmsh.model.mesh.setSizeCallback(size)
        write(out/'protocol.json',dict(case_id=case,source_STEP=str(step),source_STEP_sha256=sha(step),policy=policy,global_h_mm=hmax,local_h_mm=local,min_h_mm=hmin,transition_mm=transition,positive_Jacobian_required=True,min_SICN_required=.05,max_relative_CAD_volume_error=.001,max_nodes=800000,no_FE_started=True))
        status('generating_volume');gmsh.model.mesh.generate(3);guard();gmsh.model.mesh.optimize('Netgen')
        linear_nodes=len(gmsh.model.mesh.getNodes()[0]);assert linear_nodes*8<1000000,'Estimated quadratic mesh exceeds memory/node cap'
        status('curving_and_automatic_boundary_repair');gmsh.option.setNumber('Mesh.SecondOrderLinear',0);gmsh.model.mesh.setOrder(2);edge_review=repair_edges(types=(9,11));write(out/'automatic_edge_repair.json',edge_review)
        gmsh.option.setNumber('Mesh.HighOrderFixBoundaryNodes',1);gmsh.option.setNumber('Mesh.HighOrderThresholdMin',.1);gmsh.option.setNumber('Mesh.HighOrderThresholdMax',2.);gmsh.option.setNumber('Mesh.HighOrderNumLayers',6 if attempt<2 else 10);gmsh.option.setNumber('Mesh.HighOrderIterMax',100 if attempt<2 else 250)
        status('optimizing_curved_volume');gmsh.model.mesh.optimize('HighOrderElastic');gmsh.model.mesh.optimize('HighOrder');guard()
        status('validating_all_elements');surface_ids,faces,sj,sq=qualities(9);cell_ids,cells,vj,vq=qualities(11)
        ids,xyz,_=gmsh.model.mesh.getNodes();order=np.argsort(ids);ids=np.asarray(ids)[order];xyz=np.asarray(xyz).reshape(-1,3)[order];faces=np.searchsorted(ids,faces.reshape(-1,6));cells=np.searchsorted(ids,cells.reshape(-1,10));assert len(ids)<=800000
        # Check all six nodes of shared faces, not only their corners.
        allfaces=np.concatenate([np.sort(cells[:,idx],axis=1) for idx in [(0,1,2,4,5,6),(0,1,3,4,9,7),(0,2,3,6,8,7),(1,2,3,5,8,9)]]);unique,counts=np.unique(allfaces,axis=0,return_counts=True);boundary=unique[counts==1];stored=np.sort(faces,axis=1)
        boundary_match=np.array_equal(boundary[np.lexsort(boundary.T[::-1])],stored[np.lexsort(stored.T[::-1])]);manifold=bool(np.all((counts==1)|(counts==2)))
        edges=np.concatenate([cells[:,idx] for idx in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]);corners=np.unique(cells[:,:4]);a=np.searchsorted(corners,edges[:,0]);b=np.searchsorted(corners,edges[:,1]);graph=coo_matrix((np.ones(len(a)),(a,b)),shape=(len(corners),len(corners)));components=int(connected_components(graph,directed=False,return_labels=False))
        bary=np.vstack([np.full(4,.25),np.full((4,4),1/6)+np.eye(4)/3]);_,D,_=shape(bary[:,1:]);weights=np.array([-2/15]+[.075]*4);integral=0.
        for first in range(0,len(cells),10000):
            J=np.einsum('qak,eaj->eqjk',D,xyz[cells[first:first+10000]]);integral+=float((np.linalg.det(J)*weights).sum())
        vol_error=abs(integral/cad_volume-1);finite=bool(np.isfinite(xyz).all() and np.isfinite(sj).all() and np.isfinite(vj).all())
        passed=bool(finite and sj.min()>0 and vj.min()>0 and sq.min()>=.05 and vq.min()>=.05 and boundary_match and manifold and components==1 and vol_error<=.001)
        report=dict(case_id=case,attempt=attempt+1,policy=policy,nodes=len(ids),Tri6_faces=len(faces),C3D10_elements=len(cells),surface_nonpositive_Jacobians=int(np.sum(sj<=0)),volume_nonpositive_Jacobians=int(np.sum(vj<=0)),min_surface_Jacobian=float(sj.min()),min_volume_Jacobian=float(vj.min()),min_surface_SICN=float(sq.min()),min_volume_SICN=float(vq.min()),volume_SICN_quantiles=np.quantile(vq,[.01,.05,.5]).tolist(),boundary_six_node_match=boundary_match,manifold=manifold,connected_components=components,CAD_volume_mm3=cad_volume,mesh_volume_mm3=integral,relative_volume_error=vol_error,source_STEP_sha256=sha(step),all_finite=finite,passed_all_declared_mesh_gates=passed,automatically_repaired_curve_nodes=len(edge_review['optimized_nodes']),initial_bad_surface_elements=edge_review['initial_bad_surface_elements'],no_FE_started=True,seconds=time.time()-start)
        np.savez_compressed(out/'mesh.npz',xyz=xyz,node_ids=ids,cells=cells,cell_ids=cell_ids,faces=faces,surface_ids=surface_ids,min_sicn=vq)
        report['mesh_sha256']=sha(out/'mesh.npz');gmsh.option.setNumber('Mesh.SaveParametric',1);gmsh.write(str(out/'mesh.msh'));write(out/'validation.json',report);status('passed' if passed else 'failed_mesh_gates',review=report);print(json.dumps(report))
    except Exception as exc:
        status('failed',error=repr(exc),traceback=traceback.format_exc());raise
    finally:gmsh.finalize()
def case_job(case):
    inputs(case);base=ROOT/'cases'/case;base.mkdir(parents=True,exist_ok=True);errors=[];start=time.time()
    for attempt in range(len(POLICIES)):
        guard();folder=base/f'attempt_{attempt+1}';folder.mkdir(exist_ok=True)
        path=folder/'validation.json'
        if not path.exists():
            with (folder/'driver.log').open('w',encoding='utf-8') as log:
                proc=subprocess.Popen([sys.executable,'-B','-m','tools.rod_auto_mesh_case','--case',case,'--attempt',str(attempt)],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));write(base/'active_process.json',dict(pid=proc.pid,attempt=attempt+1,started_unix=time.time()))
                try:
                    while proc.poll() is None:
                        guard()
                        if time.time()-read(base/'active_process.json')['started_unix']>900:raise TimeoutError('900 second bounded mesh attempt')
                        time.sleep(3)
                except Exception:
                    subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],stdout=log,stderr=log);proc.wait();raise
        if path.exists():
            report=read(path)
            if report['passed_all_declared_mesh_gates']:
                acceptance=dict(case_id=case,accepted_automatic_mesh=True,attempts=attempt+1,validation=str(path),validation_sha256=sha(path),mesh=str(folder/'mesh.npz'),mesh_sha256=report['mesh_sha256'],nodes=report['nodes'],C3D10_elements=report['C3D10_elements'],min_volume_SICN=report['min_volume_SICN'],elapsed_seconds=time.time()-start,no_FE_started=True,requires_visual_review_before_solver_adoption=True);write(base/'acceptance.json',acceptance);return acceptance
            errors.append(dict(attempt=attempt+1,validation=report))
        else:errors.append(dict(attempt=attempt+1,error='Mesh process failed; log and status preserved'))
    result=dict(case_id=case,accepted_automatic_mesh=False,attempts=len(POLICIES),errors=errors,elapsed_seconds=time.time()-start,no_FE_started=True,quarantined_not_sent_to_solver=True);write(base/'failure.json',result);return result
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--case',required=True);a.add_argument('--attempt',type=int,choices=range(3));a.add_argument('--step',type=Path);a.add_argument('--parameters',type=Path);a.add_argument('--h',type=float);args=a.parse_args()
    if any(x is not None for x in [args.step,args.parameters,args.h]):
        assert args.step is not None and args.parameters is not None and args.h is not None,'External STEP needs --parameters and --h'
        assert re.fullmatch(r'[A-Za-z0-9_-]{1,80}',args.case);assert args.h>0
        parameters=read(args.parameters);parameters=parameters.get('parameters',parameters)
        required=['L','rb','rs','db','ds','Tb','Ts','H','B','tw','tf','kb','ks','lb','ls','rB','rS','rp','re','e'];assert all(k in parameters and np.isfinite(parameters[k]) for k in required)
        assert all(parameters[k]>0 for k in required if k!='e') and parameters['H']>2*parameters['tf'] and parameters['B']>parameters['tw']
        request=dict(step=str(args.step.resolve()),STEP_sha256=sha(args.step),parameters=parameters,h_mm=args.h)
        requests=ROOT/'requests';requests.mkdir(parents=True,exist_ok=True);path=requests/(args.case+'.json')
        if path.exists():assert read(path)==request,'Name already bound to different input'
        else:write(path,request)
    if args.attempt is None:print(json.dumps(case_job(args.case)))
    else:worker(args.case,args.attempt)
