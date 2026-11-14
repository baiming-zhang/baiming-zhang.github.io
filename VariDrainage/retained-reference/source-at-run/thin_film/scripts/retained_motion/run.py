"""Candidate-only single-step controls, frozen references and unchanged gates."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,time
import numpy as np
from compare import load
HERE=Path(__file__).resolve().parent;R=HERE.parents[2];O=R/'results/retained_motion_20261113';RAW=R/'results/local-runs/retained_motion_20261113/native';LOG=R/'results/logs/retained_motion_20261113'
XMAKE=R.parents[1]/'.tools/xmake/xmake/xmake.exe';SDK=R.parents[1]/'.tools/msvc'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=True,allow_nan=False)+'\n',encoding='utf-8')
def inventory(root):return [dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(root.iterdir()) if p.is_file()]
ENV={k:v for k,v in os.environ.items() if not k.startswith('VF_')};ENV['VF_THREADS']='2';ENV['NO_COLOR']='1'
def main():
    assert not O.exists() and not RAW.parent.exists() and not LOG.exists()
    g=R/'results/retained_faces_20261113/summary.json';c=R/'results/retained_checks_20261113/summary.json';assert read(g)['fixed_geometry_matrix_admitted'] and read(c)['all_fixture_checks_passed'] and read(c)['new_unit_checks']==24
    original=read(R/'results/uniform5_pilot_20261113/protocol.json');assert all(sha(R/f)==h for f,h in original['source_hashes'].items()) and sha(Path(original['executable']))==original['binary_sha256']
    refs=read(R/'results/local_assembly_20261113/status.json')['cases'];assert len(refs)==3 and all(r['passed'] for r in refs)
    for rec in refs:
        for f in rec['files']:assert sha(R/f['path'])==f['sha256']
    check=load();O.mkdir();RAW.parent.mkdir(parents=True);LOG.mkdir()
    files=[Path(__file__),HERE/'PLAN.md',HERE/'compare.py',HERE/'xmake.lua',HERE/'audit_case.py',HERE.parents[1]/'bench/tf_retained_face_bench.cpp',R/'thin_film/include/vf/thinfilm/VoronoiRetainedFacesCandidate.h',HERE.parent/'local_assembly/equivalence.py',HERE.parent/'global_space/audit_ready.py',HERE.parent/'mechanical_work/analyze.py',HERE.parent/'curved_edge/analyze.py',HERE.parent/'prepared_motion/analyze.py']
    cases=[]
    for name,dx,dt,N in [('guard_u5_dt00625',5,.0625,2112),('guard_u5_dt003125',5,.03125,2112),('u10_dt00625',10,.0625,544)]:cases.append(dict(name=name,mode='uniform',dx_um=dx,dt_us=dt,T_us=dt,N=N,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=1,refine_start_um=200,observer=True))
    save(O/'protocol.json',dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources={str(p.relative_to(R)):sha(p) for p in files},original_sources=original['source_hashes'],original_binary_sha256=original['binary_sha256'],geometry_summary_sha256=sha(g),fixture_summary_sha256=sha(c),references=refs,cases=cases,gates=dict(force_ratio=1.000001,balance=1e-10,energy=1e-7,preparation=1e-8),scope='Actual candidate single-step controls only, not time/space/fullC1 qualification'))
    state=dict(state='building',cases=[],skipped=[],full_C1_complete=False);save(O/'status.json',state);attempts=[]
    for name,cmd,cap in [('configure',[str(XMAKE),'f','-y','-P',str(HERE),'-p','windows','-a','x64','-m','release','--sdk='+str(SDK)],300),('build',[str(XMAKE),'build','-y','-P',str(HERE),'-j2','tf_retained_face_bench'],600)]:
        log=LOG/(name+'.log');start=time.monotonic()
        with log.open('w',encoding='utf-8') as stream:
            try:rc=subprocess.run(cmd,cwd=R,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=cap).returncode
            except subprocess.TimeoutExpired:rc='timeout'
        attempts.append(dict(name=name,command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log)));save(O/'build_attempts.json',attempts);assert rc==0,'Preserve build failure, no retry'
    exe=RAW.parent/'bin/tf_retained_face_bench.exe';assert exe.exists();save(O/'build.json',dict(binary=str(exe.relative_to(R)),sha256=sha(exe)))
    ready=True;state['state']='running_matrix'
    for case in cases:
        name=case['name']
        if not ready:state['skipped'].append(dict(name=name,reason='Earlier candidate native/equivalence failure'));continue
        dest=RAW/name;assert not dest.exists();cmd=[str(XMAKE),'run','-P',str(HERE),'tf_retained_face_bench','uniform',str(case['dx_um']),str(case['dt_us']),str(case['T_us']),'500','0',str(dest),'6','variational','1','200','1'];log=LOG/(name+'.log');start=time.monotonic()
        with log.open('w',encoding='utf-8') as stream:
            try:rc=subprocess.run(cmd,cwd=R,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=240).returncode
            except subprocess.TimeoutExpired:rc='timeout'
        rec=dict(name=name,parameters=case,command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log),completed=False,native_gate=False)
        if (dest/'summary.json').exists():
            m=read(dest/'summary.json');rec.update(summary=m,preparation=read(dest/'preparation.json'),mesh=read(dest/'mesh.json'));rec['completed']=rc==0 and m['steps']==1 and m['requested_steps']==1 and m['failures']==0 and m['N']==case['N']
            if rec['completed']:
                mechanics=np.atleast_1d(np.genfromtxt(dest/'mechanics.csv',delimiter=',',names=True));ratio=float(max(mechanics['g1_N']/mechanics['target_N']));rec['maximum_force_target_ratio']=ratio;rec['native_gate']=len(mechanics)==1 and ratio<=1.000001 and max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and m['max_budget_excess']<=1e-7 and rec['preparation']['pressure_residual']<=1e-8
        rec['files']=inventory(dest)
        if rec['native_gate']:
            ref=next(r for r in refs if r['name']==name);root=next(R/f['path'] for f in ref['files'] if Path(f['path']).name=='states.csv').parent;eq=check(root,dest,case['dx_um']*1e-6);save(O/('equivalence_'+name+'.json'),eq);rec['first_step_equivalence_passed']=eq['passed']
        rec['passed']=bool(rec['native_gate'] and rec.get('first_step_equivalence_passed',False));ready=rec['passed'];state['cases'].append(rec);save(O/'status.json',state);print(name,'native',rec['native_gate'],'first-step',rec.get('first_step_equivalence_passed'),round(rec['wall_s'],3),flush=True)
    assert all(sha(R/f)==h for f,h in original['source_hashes'].items()) and sha(Path(original['executable']))==original['binary_sha256'];state.update(state='completed',original_native_capture_unchanged=True);save(O/'status.json',state)
if __name__=='__main__':main()
