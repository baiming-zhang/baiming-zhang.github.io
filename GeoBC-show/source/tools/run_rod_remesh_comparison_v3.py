"""Explicitly authorized 0006 remesh comparison, same approved physical recipe."""
import argparse,hashlib,json,msvcrt,os,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
from geobc.rod_io import write_json
from tools.run_rod_v2_dataset import resources

REPO=Path(__file__).resolve().parents[1];WS=REPO.parents[1];ROOT=WS/'data/rod_remesh_comparison_0006_v3';CASE='rodv2_0006'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def status(phase,**kw):write_json(ROOT/'worker_status.json',dict(phase=phase,pid=os.getpid(),updated_unix=time.time(),**kw))
def prepare_mesh(protocol):
    import gmsh
    src=Path(protocol['new_mesh']);m=dict(np.load(src));source=read(WS/'data/rod_v2_sixload_v1/cases'/CASE/'h3/mesh.json');p=protocol['parameters']
    assert sha(src)==protocol['new_mesh_sha256'] and source['source_step_sha256']==protocol['STEP_sha256']
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.NumThreads',1);gmsh.option.setNumber('General.Terminal',0);gmsh.model.add('physical_face_adapter');gmsh.model.occ.importShapes(source['source_step']);gmsh.model.occ.synchronize();gmsh.merge(str(src.with_suffix('.msh')))
        groups=np.full(len(m['faces']),3,int);sid={int(x):i for i,x in enumerate(m['surface_ids'])};found={1:[],2:[]}
        for dim,tag in gmsh.model.getEntities(2):
            if gmsh.model.getType(2,tag)!='Cylinder':continue
            area=gmsh.model.occ.getMass(2,tag);group=1 if abs(area-2*np.pi*p['rb']*(p['Tb']-2*p['re']))<1e-3 else 2 if abs(area-2*np.pi*p['rs']*(p['Ts']-2*p['re']))<1e-3 else None
            if group:
                found[group].append(tag);ids,_=gmsh.model.mesh.getElementsByType(9,tag)
                for ident in ids:groups[sid[int(ident)]]=group
        assert len(found[1])==len(found[2])==1 and np.any(groups==1) and np.any(groups==2)
    finally:gmsh.finalize()
    m['face_group_ids']=groups;np.savez_compressed(ROOT/'mesh_native.npz',**m)
    members=np.unique(m['faces'][groups==1]);normal=np.zeros_like(m['xyz']);normal[members,:2]=m['xyz'][members,:2]-[p['L'],0];normal[members]/=np.linalg.norm(normal[members],axis=1,keepdims=True)
    oldbc=np.load(WS/'data/rod_v2_sixload_v1/cases'/CASE/'h3/constraints.npz');oldm=np.load(WS/'data/rod_v2_sixload_v1/cases'/CASE/'h3/mesh.npz');anchor=oldm['xyz'][int(oldbc['axial_gauge_original_index'])];idx=int(members[np.argmin(np.linalg.norm(m['xyz'][members]-anchor,axis=1))])
    bc=dict(node_ids=m['node_ids'],radial_normal=normal.astype('float32'),axial_gauge_original_index=idx);np.savez_compressed(ROOT/'constraints.npz',**bc)
    write_json(ROOT/'mesh_adapter_review.json',dict(source_mesh_sha256=sha(src),native_mesh_sha256=sha(ROOT/'mesh_native.npz'),unchanged_xyz_and_connectivity=True,big_bore_faces=int(np.sum(groups==1)),small_bore_faces=int(np.sum(groups==2)),big_bore_supported_nodes=len(members),physical_CAD_surface_tags=found,point_anchor_used_only_by_removed_legacy_deck_boundary=True))
    return m,bc
def worker():
    from geobc.rod_distributed_gauge import coefficients,replace_boundary,equilibrium
    from geobc.rod_revised_bending_load import assemble
    from geobc.rod_v3_ccx import write_deck,run_ccx,export_result
    from tools.run_rod_mixed1000_distributed import histories
    from geobc.rod_v2_probes import locate,evaluate
    protocol=read(ROOT/'protocol.json')
    for name,digest in protocol['source_sha256'].items():assert sha(REPO/name)==digest
    status('preparing_native_mesh');m,bc=prepare_mesh(protocol);p=protocol['parameters'];material=read(WS/'data/rod_v3_physics/material_reference.json');assert sha(WS/'data/rod_v3_physics/material_reference.json')==protocol['material_sha256']
    c=coefficients(m,bc,p['L']);status('assembling_same_physical_load');force,load=assemble(m,p);old=Path(protocol['old_target_fields']).parent;oldload=read(old/'load.json')
    for key in ['theta_deg','omega_rod_rad_s','alpha_rod_rad_s2','net_pressure_MPa','a_small_local_mm_s2']:assert np.array_equal(np.asarray(load['snapshot'][key]),np.asarray(oldload['snapshot'][key]))
    assert load['width_skew']==oldload['width_skew']==.3
    folder=ROOT/'solution';folder.mkdir(exist_ok=True);np.savez_compressed(folder/'input_force.npz',node_ids=m['node_ids'],F_N=force);write_json(folder/'load.json',load)
    deck=folder/'rod.inp'
    if not deck.exists():
        write_deck(deck,m,bc,force,material,plastic=True,nlgeom=True,unload=True,increment=.1);text=deck.read_text().replace('*CLOAD','*CONTROLS,PARAMETERS=FIELD\n1e-5,1e-4,,,1e-5,1e-5,1e-4,1e-8\n*CLOAD').replace('FREQUENCY=999999','FREQUENCY=1');deck.write_text(text,encoding='ascii');replace_boundary(deck,m,c)
        write_json(folder/'provenance.json',dict(study_protocol_sha256=sha(ROOT/'protocol.json'),mesh_sha256=sha(ROOT/'mesh_native.npz'),source_mesh=str(ROOT/'mesh_native.npz'),force_sha256=sha(folder/'input_force.npz'),deck_sha256=sha(deck),same_physical_recipe_as_V4_accepted_0006_bending=True,material_sha256=protocol['material_sha256'],boundary='area_integral_tangent_axial_plus_original_reference_radial'))
    else:assert sha(deck)==read(folder/'provenance.json')['deck_sha256']
    log=folder/'driver.log';finished=log.exists() and 'Job finished' in log.read_text(errors='replace') and '*ERROR' not in log.read_text(errors='replace')
    status('reusing_completed_solver_output' if finished else 'solving_FE',nodes=len(m['node_ids']),C3D10_elements=len(m['cells']),matrix_threads=protocol['matrix_threads'])
    seconds=None if finished else run_ccx(WS/'.runtime/calculix_2_23/ccx_static.exe',deck,timeout=protocol['timeout_seconds'],threads=protocol['matrix_threads'])
    if not (folder/'summary.json').exists():export_result(deck,m,force,material,dict(case_id='combined_bending',geometry_case_id=CASE,recipe='same_V4_bending15_distributed_remesh_comparison',plastic=True,unload=True,load=load),seconds)
    status('checking_native_fields_and_all_states');fields=np.load(folder/'fields.npz');assert np.array_equal(fields['node_ids'],m['node_ids']) and fields['U_peak_mm'].dtype==np.float64 and all(np.isfinite(x).all() for x in fields.values())
    assert np.array_equal(fields['force_peak_N'],force);eq=equilibrium(folder,m,bc,c,p);hist=histories(folder,material['plastic_strain_data_limit']);summary=read(folder/'summary.json')
    status('comparing_exact_physical_probes');reference=Path(protocol['old_reference_fields']).parent;refmesh=Path(read(WS/'data/rod_mixed_1000_v4_bending15/pilot'/CASE/'combined_bending/reference/provenance.json')['source_mesh']);points=np.load(refmesh.parent/'probes.npz')['points'];newel,newN,newerr=locate(m['xyz'],m['cells'],points);newU=evaluate(m['cells'],fields['U_peak_mm'],newel,newN);comparisons={}
    for label,fieldpath in [('old_target',Path(protocol['old_target_fields'])),('old_reference',Path(protocol['old_reference_fields']))]:
        data=np.load(fieldpath);meshpath=WS/'data/rod_v2_sixload_v1/cases'/CASE/'h3/mesh.npz' if label=='old_target' else refmesh;mesh=np.load(meshpath);el,N,err=locate(mesh['xyz'],mesh['cells'],points);oldU=evaluate(mesh['cells'],data['U_peak_mm'],el,N);diff=newU-oldU;oldsummary=read(fieldpath.parent/'summary.json');oldeq=read(fieldpath.parent/'equilibrium_review.json')
        changes=dict(U_probe_relative_L2=float(np.linalg.norm(diff)/np.linalg.norm(oldU)),U_probe_absolute_RMSE_mm=float(np.sqrt(np.mean(diff**2))),U_component_RMSE_mm=np.sqrt(np.mean(diff**2,axis=0)).tolist(),max_U_relative=abs(summary['max_U_mm']/oldsummary['max_U_mm']-1),force_dot_U_relative=abs(eq['endpoint_force_dot_U_Nmm']/oldeq['endpoint_force_dot_U_Nmm']-1),stress_ROI_relative={k:abs(eq['GP_reference_volume_weighted_stress_ROI_rms_MPa'][k]/v-1) for k,v in oldeq['GP_reference_volume_weighted_stress_ROI_rms_MPa'].items()},raw_GP_peak_VM_MPa=oldsummary['max_GP_von_Mises_MPa'],old_FE_seconds=oldsummary['seconds'])
        comparisons[label]=changes;np.savez_compressed(ROOT/(label+'_probe_comparison.npz'),points=points,old_U_mm=oldU,new_U_mm=newU,difference_mm=diff)
    report=dict(case_id=CASE,scenario='combined_bending',net_pressure_MPa=1.5,width_skew=.3,theta_deg=60,rpm=4500,new_nodes=len(m['node_ids']),new_elements=len(m['cells']),new_FE_seconds=summary['seconds'],new_raw_GP_peak_VM_MPa=summary['max_GP_von_Mises_MPa'],new_max_U_mm=summary['max_U_mm'],new_max_PEEQ=summary['max_PEEQ'],all_states=hist,equilibrium=eq,comparisons=comparisons,physical_probes=2048,no_rigid_alignment=True,new_mapping_max_error_mm=float(newerr),increment_review_required=summary['max_PEEQ']>0,new_fields_sha256=sha(folder/'fields.npz'),old_results_preserved=True,no_automatic_dataset_replacement=True,formal_training_started=False)
    write_json(ROOT/'comparison_review.json',report);status('complete_pending_root_review')
    try:
        subprocess.run([sys.executable,'-B','-m','tools.render_rod_remesh_comparison_v3'],cwd=REPO,check=True)
    except Exception as exc:write_json(ROOT/'plot_warning.json',dict(error=repr(exc),results_preserved=True))
def launch():
    ROOT.mkdir(exist_ok=True);start=time.time()
    lock=(ROOT/'run.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    # Serialize the large FE comparison after the geometry-only validation matrix.
    while True:
        matrix=read(REPO/'results/rod_auto_mesh_v1/matrix_status.json');ram,disk=resources()
        matrix_done=matrix['phase'].startswith('complete')
        if matrix_done and ram>=32*1024**3 and disk>=20*1024**3:break
        write_json(ROOT/'run_status.json',dict(phase='waiting_mesh_matrix_and_memory',matrix_phase=matrix['phase'],matrix_passed=matrix['passed'],matrix_failed=matrix['failed'],available_RAM_bytes=ram,elapsed_seconds=time.time()-start,updated_unix=time.time(),no_new_FE_started=True))
        time.sleep(10)
    # This comparison is explicitly authorized independently of the mesh matrix.
    with (ROOT/'worker.log').open('a',encoding='utf-8') as log:
        proc=subprocess.Popen([sys.executable,'-B','-m','tools.run_rod_remesh_comparison_v3','--worker'],cwd=REPO,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        write_json(ROOT/'active_process.json',dict(parent_pid=os.getpid(),worker_pid=proc.pid,started_unix=start))
        try:
            while proc.poll() is None:
                ram,disk=resources()
                if ram<16*1024**3 or disk<20*1024**3:raise RuntimeError('Comparison reserve guard: 16GiB available RAM/20GiB disk, protecting existing jobs')
                time.sleep(2)
            if proc.returncode:raise RuntimeError(f'Comparison worker exit{proc.returncode}; original logs preserved')
            write_json(ROOT/'run_status.json',dict(phase='complete_pending_root_review',elapsed_seconds=time.time()-start,updated_unix=time.time()))
        except Exception as exc:
            if proc.poll() is None:subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],stdout=log,stderr=log);proc.wait()
            write_json(ROOT/'run_status.json',dict(phase='stopped_or_failed',error=repr(exc),elapsed_seconds=time.time()-start,updated_unix=time.time(),original_outputs_preserved=True));raise
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--worker',action='store_true');a.add_argument('--launch',action='store_true');args=a.parse_args()
    if args.worker:worker()
    elif args.launch:launch()
    else:raise SystemExit('Specify --launch')
