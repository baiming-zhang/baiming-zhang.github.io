"""Execute the originally planned observations after a path-only driver fix."""
import json,subprocess,time,os
from pathlib import Path
from run import R,O,RAW,LOG,HERE,XMAKE,sha,save

assert not RAW.exists() and not (LOG/'observe.log').exists()
p=json.loads((O/'protocol.json').read_text());patch=json.loads((O/'driver_path_patch.json').read_text());attempts=json.loads((O/'attempts.json').read_text())
assert [a['name'] for a in attempts]==['configure','build'] and all(a['returncode']==0 for a in attempts)
assert sha(HERE/'run.py')==patch['corrected_driver_sha256']
for f,h in p['sources'].items():assert sha(O/'source-at-build/run.py' if f.endswith('short_face\\run.py') else R/f)==h
assert all(sha(R/f)==h for f,h in p['original_sources'].items())
old=json.loads((R/'results/uniform5_pilot_20261113/protocol.json').read_text());assert sha(Path(old['executable']))==p['original_binary_sha256']
exe=Path(patch['observed_binary']);assert exe.exists()
save(O/'build.json',dict(binary=str(exe.relative_to(R)),sha256=sha(exe),path_patch_sha256=sha(O/'driver_path_patch.json')))
env={k:v for k,v in os.environ.items() if not k.startswith('VF_')};env['VF_THREADS']='2';env['NO_COLOR']='1'
cmd=[str(XMAKE),'run','-P',str(HERE),'tf_short_face_probe',str(RAW)];log=LOG/'observe.log';start=time.monotonic()
with log.open('w',encoding='utf-8') as stream:
    try:rc=subprocess.run(cmd,cwd=R,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=30).returncode
    except subprocess.TimeoutExpired:rc='timeout'
attempts.append(dict(name='observe',command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log),resume_source_sha256=sha(Path(__file__))))
save(O/'attempts.json',attempts)
assert rc==0,'Preserve native failure; no automatic retry'
assert all(sha(R/f)==h for f,h in p['original_sources'].items()) and sha(Path(old['executable']))==p['original_binary_sha256']
save(O/'raw_manifest.json',{str(f.relative_to(R)):dict(bytes=f.stat().st_size,sha256=sha(f)) for f in RAW.iterdir() if f.is_file()})
save(O/'status.json',dict(state='completed',diagnostic_observations=3,new_native_trajectory_steps=0,original_native_capture_unchanged=True,driver_path_patch=True,full_C1_complete=False,scope='Probe executed once after path-only driver fix; expected pattern not yet assessed'))
print((RAW/'observations.csv').read_text())
