"""Seven frozen inputs, original and candidate, one build and no native replay."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,time
import numpy as np
HERE=Path(__file__).resolve().parent;R=HERE.parents[2]
O=R/'results/retained_faces_20261113';LOCAL=R/'results/local-runs/retained_faces_20261113';LOG=R/'results/logs/retained_faces_20261113'
XMAKE=R.parents[1]/'.tools/xmake/xmake/xmake.exe';SDK=R.parents[1]/'.tools/msvc'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=True,allow_nan=False)+'\n',encoding='utf-8')
def source_check():
    original=R/'thin_film/include/vf/thinfilm/VoronoiDrainage.h';s=original.read_text(encoding='utf-8');a=s.index('inline VoronoiDrainageGeometry Make_Voronoi_Drainage(');b=s.index('// Region a->b',a);f=s[a:b]
    expected=f.replace('Make_Voronoi_Drainage(', 'Make_Voronoi_Drainage_RetainedFaces(',1).replace('if(length<=eps)continue; // only machine-scale degenerate segments','if(!(length>0.))continue; // candidate: retain every positive edge after original cleanup')
    candidate=R/'thin_film/include/vf/thinfilm/VoronoiRetainedFacesCandidate.h'
    assert candidate.read_text(encoding='utf-8')=='// Opt-in diagnostic candidate; original geometry/defaults unchanged.\n#pragma once\n#include "vf/thinfilm/VoronoiDrainage.h"\nnamespace vf::thinfilm {\n'+expected+'}\n'
    return dict(original_sha256=sha(original),candidate_sha256=sha(candidate),exact_two_text_changes_verified=True)
def main():
    assert not O.exists() and not LOCAL.exists() and not LOG.exists()
    native=read(R/'results/uniform5_pilot_20261113/protocol.json');assert all(sha(R/f)==h for f,h in native['source_hashes'].items());assert sha(Path(native['executable']))==native['binary_sha256']
    change=source_check();O.mkdir();LOCAL.mkdir();LOG.mkdir();inputs=O/'inputs';inputs.mkdir();cases=[]
    for name,delta in [('unperturbed',0),('short_face',8e-9),('resolved_face',8e-7)]:
        xy=np.array([(x,y) for y in [-1.,0.,1.] for x in [-1.,0.,1.]]);xy[8,1]-=delta;p=inputs/(name+'.txt')
        with p.open('w',encoding='utf-8') as stream:
            stream.write('-32 -16 32 16 9\n')
            for i,(x,y) in enumerate(xy):stream.write(f'{i} {x:.17g} {y:.17g} -1\n')
        cases.append(dict(name=name,N=9,kind='manufactured',input=str(p.relative_to(R)),input_sha256=sha(p)))
    status=read(R/'results/uniform5_pilot_20261113/status.json');assert status['state']=='completed'
    refs=read(R/'results/global_space_20261113/protocol.json')['reference_cases']+native['reference_cases']+status['cases']
    for name,step,N in [('u20_dt003125',320,144),('u10_dt003125',320,544),('u5_dt00625',160,2112),('u5_dt003125',40,2112)]:
        rec=next(c for c in refs if c['name']==name);f=next(f for f in rec['files'] if Path(f['path']).name=='states.csv');path=R/f['path'];assert sha(path)==f['sha256']
        data=np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True));a=data[(data['step']==step)&(data['stage']==2)];assert len(a)==N and np.array_equal(a['id'],np.arange(N));p=inputs/(name+'.txt')
        with p.open('w',encoding='utf-8') as stream:
            stream.write(f'0 0 {320e-6/1e-3:.17g} {160e-6/1e-3:.17g} {N}\n')
            for row in a:stream.write(f'{int(row["id"])} {row["x_m"]/1e-3:.17g} {row["y_m"]/1e-3:.17g} {row["A_m2"]/1e-6:.17g}\n')
        cases.append(dict(name=name,N=N,kind='saved_actual_geometry',step=step,time_us=step*rec['summary']['dt_s']*1e6,source=f,input=str(p.relative_to(R)),input_sha256=sha(p),trajectory_replayed=False,failed_case_prefix_only=not rec['native_gate']))
    files=[Path(__file__),HERE/'PLAN.md',HERE/'xmake.lua',HERE.parents[1]/'bench/tf_retained_face_probe.cpp',R/'thin_film/include/vf/thinfilm/VoronoiRetainedFacesCandidate.h']
    save(O/'protocol.json',dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources={str(p.relative_to(R)):sha(p) for p in files},change=change,original_sources=native['source_hashes'],original_binary_sha256=native['binary_sha256'],cases=cases,area_import_gate=1e-9,normal_closure_gate=1e-9,area_closure_gate=1e-10,orthogonality_certificate_gate=1,new_native_trajectory_steps=0,full_C1_complete=False,scope='Geometry-only candidate admission; no derivatives, dynamics or actual78 replay'))
    attempts=[];env={k:v for k,v in os.environ.items() if not k.startswith('VF_')};env['VF_THREADS']='2';env['NO_COLOR']='1'
    def execute(name,cmd,cap):
        start=time.monotonic();log=LOG/(name+'.log')
        with log.open('w',encoding='utf-8') as stream:
            try:rc=subprocess.run(cmd,cwd=R,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=cap).returncode
            except subprocess.TimeoutExpired:rc='timeout'
        attempts.append(dict(name=name,command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log)));save(O/'attempts.json',attempts)
        assert rc==0, name+' failed; preserve record, no automatic retry'
    save(O/'status.json',dict(state='building',cases=[],full_C1_complete=False))
    execute('configure',[str(XMAKE),'f','-y','-P',str(HERE),'-p','windows','-a','x64','-m','release','--sdk='+str(SDK)],300)
    execute('build',[str(XMAKE),'build','-y','-P',str(HERE),'-j2','tf_retained_face_probe'],300)
    exe=LOCAL/'bin/tf_retained_face_probe.exe';assert exe.exists();save(O/'build.json',dict(binary=str(exe.relative_to(R)),sha256=sha(exe)))
    done=[]
    for c in cases:
        assert sha(R/c['input'])==c['input_sha256'];execute(c['name'],[str(XMAKE),'run','-P',str(HERE),'tf_retained_face_probe',str(R/c['input']),str(LOCAL/'native'/c['name'])],60);done.append(c['name']);save(O/'status.json',dict(state='running_matrix',cases=done,full_C1_complete=False))
    assert all(sha(R/f)==h for f,h in native['source_hashes'].items()) and sha(Path(native['executable']))==native['binary_sha256']
    save(O/'raw_manifest.json',{str(p.relative_to(R)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in (LOCAL/'native').rglob('*') if p.is_file()})
    save(O/'status.json',dict(state='completed',cases=done,new_native_trajectory_steps=0,original_native_capture_unchanged=True,full_C1_complete=False));print('Completed seven paired geometry observations; no trajectories replayed.')
if __name__=='__main__':main()
