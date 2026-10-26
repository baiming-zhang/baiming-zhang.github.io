"""One immutable C04 mesh attempt, then original gates with summed subvolumes."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,ast,hashlib,json,os,shutil,shlex,subprocess,sys,time
import numpy as np
import gmsh
W=Path(__file__).resolve().parents[1];P=W/'publication/GeoBC-Neural-Operator';BASE=W/'data/copv/c04_repair_v3';O=BASE
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
    global O
    ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--attempt',required=True);ap.add_argument('--partition',type=Path,required=True);ap.add_argument('--h',type=float,default=5.);ap.add_argument('--shell-h',type=float,default=4.);ap.add_argument('--bond-h',type=float,default=3.);ap.add_argument('--root-h',type=float,default=.75);ap.add_argument('--circle-segments',type=int,default=128);ap.add_argument('--recipe',choices=['matched','baseline'],default='matched');ap.add_argument('--max-cells',type=int,default=2000000);ap.add_argument('--backend',choices=['windows','linux'],default='windows');ap.add_argument('--linear-checkpoint',type=Path);args=ap.parse_args();start=time.perf_counter();root=BASE/args.case/args.attempt;root.mkdir(parents=True,exist_ok=False);O=root
    source=P/'tools/mesh_copv_c04_guarded_repair.py';frozen=root/'mesher_executed.py';shutil.copyfile(source,frozen);partition=args.partition.resolve();assert read(partition/'partition_review.json')['case_id']==args.case;cmd=[sys.executable,'-B',str(frozen),'--partition',str(partition/'structural_partition.brep'),'--review',str(partition/'partition_review.json'),'--output',str(root/'mesh'),'--h',str(args.h),'--shell-h',str(args.shell_h),'--bond-h',str(args.bond_h),'--root-h',str(args.root_h),'--circle-segments',str(args.circle_segments)]
    cmd.extend(['--recipe',args.recipe,'--max-cells',str(args.max_cells)])
    if args.linear_checkpoint:cmd.extend(['--linear-checkpoint',str(args.linear_checkpoint.resolve())])
    if args.backend=='linux':
        def linux(value):
            s=str(value).replace('\\','/');return '/mnt/'+s[0].lower()+s[2:] if len(s)>2 and s[1]==':' else s
        linuxcmd=[linux(W/'.runtime/fenicsx_cuda/env/bin/python')]+[linux(s) for s in cmd[1:]]
        setup='#!/bin/bash\nset -e\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\nexport COPV_C04_VOLUME_ALGORITHM='+os.environ.get('COPV_C04_VOLUME_ALGORITHM','1')+'\nexport PYTHONPATH='+shlex.quote(linux(W/'.runtime/c04_linux_mesher/deps')+':'+linux(P))+'\nexport LD_LIBRARY_PATH='+shlex.quote(linux(W/'.runtime/c04_linux_mesher/native-libs/usr/lib/x86_64-linux-gnu'))+':"${LD_LIBRARY_PATH:-}"\nexec '+shlex.join(linuxcmd)+'\n'
        (root/'run_mesh_linux.sh').write_bytes(setup.encode());cmd=['wsl','-d','cold-surge-ubuntu-24.04','--','bash',linux(root/'run_mesh_linux.sh')]
    if args.backend=='linux':
        shellfile=root/'run_mesh_linux.sh';shelltext=shellfile.read_text();shelltext=shelltext.replace('\nexec ', '\nexport COPV_C04_LINEAR_OPTIMIZER='+shlex.quote(os.environ.get('COPV_C04_LINEAR_OPTIMIZER','none'))+'\nexport COPV_C04_COLLAR_FRONTAL='+shlex.quote(os.environ.get('COPV_C04_COLLAR_FRONTAL','0'))+'\nexec timeout --signal=TERM --kill-after=10s 850s ');shellfile.write_bytes(shelltext.encode())
        extra=''
        for key in ['COPV_C04_SURFACE_CHECKPOINT','COPV_C04_MAX_SURFACE_NODES','COPV_C04_ROOT_BAND_MM','COPV_C04_MATCH_ROOT_CIRCLES','COPV_C04_FRONTAL_PARTS','COPV_C04_ADJACENT_ROOT_BAND']:
            if key in os.environ:extra+='export '+key+'='+shlex.quote(linux(os.environ[key]) if key.endswith('CHECKPOINT') else os.environ[key])+'\n'
        shellfile.write_bytes(shelltext.replace('\nexec ','\n'+extra+'exec ').encode())
    state=dict(case_id=args.case,attempt=args.attempt,pid=os.getpid(),phase='meshing',started=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),command=cmd,partition_file=str(partition/'structural_partition.brep'),backend=args.backend,max_cells=args.max_cells,linear_optimizer=os.environ.get('COPV_C04_LINEAR_OPTIMIZER','none'),source_sha256=hashlib.sha256(frozen.read_bytes()).hexdigest(),partition_sha256=hashlib.sha256((partition/'structural_partition.brep').read_bytes()).hexdigest(),new_FE_labels=0)
    def save():
        state['elapsed_seconds']=time.perf_counter()-start;state['updated']=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat();write(root/'state.json',state);progress=read(BASE/'progress.json');progress.update(updated=state['updated'],phase=state['phase'],current_attempt=str(root.relative_to(W)).replace('\\','/'),current_pid=os.getpid());progress['attempts']=[r for r in progress['attempts'] if not(r['case_id']==args.case and r['attempt']==args.attempt)]+[dict(case_id=args.case,attempt=args.attempt,status=state['phase'],state_file=str((root/'state.json').relative_to(W)).replace('\\','/'))];write(BASE/'progress.json',progress)
    save();env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(P),str(W/'.rod_deps'),str(W/'.local_deps')]);env.pop('COPV_C02_NETGEN',None)
    try:
        with (root/'mesh_execution.log').open('w',encoding='utf-8') as log:run=subprocess.run(cmd,cwd=W,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        if run.returncode:state.update(phase='mesh_failed_retained',exit_code=run.returncode);save();print(json.dumps(state),flush=True);return
        state['phase']='independent_surface_and_quality_review';save();params=W/'data/copv/c04_pilot_v1/cases'/args.case/'parameters.json';surface=root/'mesh/surface'
        with (root/'surface_execution.log').open('w',encoding='utf-8') as log:run=subprocess.run([sys.executable,'-B','-m','tools.prepare_copv_c04_surface','--mesh',str(root/'mesh/mesh.npz'),'--output',str(surface),'--parameters',str(params)],cwd=P,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        if run.returncode:raise RuntimeError('Independent surface topology/load closure failed')
        tree=ast.parse((W/'notes/export_c03_cases.py').read_text(encoding='utf-8'));node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='quality');code=ast.get_source_segment((W/'notes/export_c03_cases.py').read_text(encoding='utf-8'),node);code=code.replace("expected={r['part_id']:r['volume_mm3'] for r in review['volume_records']}","expected={}\n    for record in review['volume_records']: expected[record['part_id']]=expected.get(record['part_id'],0.)+record['volume_mm3']")
        code=code.replace("for record in review['volume_records']: expected[record['part_id']]=expected.get(record['part_id'],0.)+record['volume_mm3']", "for record in review['volume_records']: expected[record['part_id']]=expected.get(record['part_id'],0.)+record['volume_mm3']\n    if 'quality_reference_part_volumes_mm3' in review: expected={int(k):v for k,v in review['quality_reference_part_volumes_mm3'].items()}")
        (root/'quality_executed_function.py').write_text(code,encoding='utf-8');exec(compile(code,str(root/'quality_executed_function.py'),'exec'),globals());q=globals()['quality'](root/'mesh',surface,'quality');mr=read(root/'mesh/mesh_review.json');closure=read(surface/'load_closure_review.json');q['production_gate'].update(connected_components=mr['connected_components'],independent_load_closure_passed=closure['passed']);q['production_gate']['passed']=bool(q['production_gate']['passed'] and mr['connected_components']==1 and closure['passed'] and q['total_relative_volume_error']<=.001);q['production_gate']['status']='passed_mesh_only_no_FE_yet' if q['production_gate']['passed'] else 'quality_rejected_no_FE';q['mesh_sha256']=hashlib.sha256((root/'mesh/mesh.npz').read_bytes()).hexdigest();write(root/'quality/quality.json',q);state.update(phase='mesh_quality_passed_no_FE_yet' if q['production_gate']['passed'] else 'mesh_quality_rejected_retained',quality={k:q[k] for k in ['cells','nodes','min_volume_sicn','min_surface_sicn','below_005','max_part_relative_volume_error','total_relative_volume_error','production_gate']});save();print(json.dumps(state),flush=True)
    except BaseException as error:state.update(phase='review_failed_retained',error=str(error));save();raise
if __name__=='__main__':main()
