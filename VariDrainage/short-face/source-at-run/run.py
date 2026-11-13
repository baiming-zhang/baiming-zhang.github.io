"""One isolated build and fixed geometry observation matrix; no retries."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,time

HERE=Path(__file__).resolve().parent;R=HERE.parents[2]
O=R/'results/short_face_20261113';RAW=R/'results/local-runs/short_face_20261113/native';LOG=R/'results/logs/short_face_20261113'
XMAKE=R.parents[1]/'.tools/xmake/xmake/xmake.exe';SDK=R.parents[1]/'.tools/msvc'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def main():
    assert not O.exists() and not RAW.exists() and not LOG.exists()
    prior=R/'results/uniform5_pilot_20261113/protocol.json';native=json.loads(prior.read_text())
    assert all(sha(R/f)==h for f,h in native['source_hashes'].items())
    assert sha(Path(native['executable']))==native['binary_sha256']
    O.mkdir();LOG.mkdir();RAW.parent.mkdir(parents=True,exist_ok=True)
    files=[Path(__file__),HERE/'PLAN.md',HERE/'xmake.lua',HERE.parents[1]/'bench/tf_short_face_probe.cpp']
    env={k:v for k,v in os.environ.items() if not k.startswith('VF_')};env['VF_THREADS']='2';env['NO_COLOR']='1'
    save(O/'protocol.json',dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources={str(f.relative_to(R)):sha(f) for f in files},original_native_protocol_sha256=sha(prior),original_sources=native['source_hashes'],original_binary_sha256=native['binary_sha256'],xmake=str(XMAKE),xmake_sha256=sha(XMAKE),sdk=str(SDK),cases=[dict(name=n,delta=d) for n,d in [('unperturbed',0),('short_face',8e-9),('resolved_face',8e-7)]],tolerance=1e-10,normal_closure_gate=1e-9,new_native_trajectory_steps=0,environment={'VF_THREADS':'2','caller_VF':'excluded'},scope='Fixed manufactured counterexample, not actual failure78 replay or solver repair'))
    save(O/'status.json',dict(state='building',full_C1_complete=False))
    commands=[('configure',[str(XMAKE),'f','-y','-P',str(HERE),'-p','windows','-a','x64','-m','release','--sdk='+str(SDK)],300),('build',[str(XMAKE),'build','-y','-P',str(HERE),'-j2','tf_short_face_probe'],300),('observe',[str(XMAKE),'run','-P',str(HERE),'tf_short_face_probe',str(RAW)],30)]
    attempts=[]
    for name,cmd,cap in commands:
        start=time.monotonic();log=LOG/(name+'.log')
        with log.open('w',encoding='utf-8') as stream:
            try:rc=subprocess.run(cmd,cwd=R,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=cap).returncode
            except subprocess.TimeoutExpired:rc='timeout'
        record=dict(name=name,command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log));attempts.append(record)
        save(O/'attempts.json',attempts)
        if rc!=0:
            save(O/'status.json',dict(state='failed',failed_stage=name,full_C1_complete=False));raise RuntimeError(f'{name} failed; retain log, no automatic retry')
        if name=='build':
            exe=R/'build/windows/x64/release/tf_short_face_probe.exe';assert exe.exists()
            save(O/'build.json',dict(binary=str(exe.relative_to(R)),sha256=sha(exe)))
    assert all(sha(R/f)==h for f,h in native['source_hashes'].items()) and sha(Path(native['executable']))==native['binary_sha256']
    save(O/'raw_manifest.json',{str(f.relative_to(R)):dict(bytes=f.stat().st_size,sha256=sha(f)) for f in RAW.iterdir() if f.is_file()})
    save(O/'status.json',dict(state='completed',diagnostic_observations=3,new_native_trajectory_steps=0,original_native_capture_unchanged=True,full_C1_complete=False,scope='Probe executed, expected counterexample pattern not yet assessed'))
    print('Completed fixed geometry observations; original long-run sources/binary unchanged.')

if __name__=='__main__':main()
