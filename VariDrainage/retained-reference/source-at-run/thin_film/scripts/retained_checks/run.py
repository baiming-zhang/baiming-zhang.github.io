"""Run original curved fixture bodies under an opt-in geometry candidate once."""
from pathlib import Path
import hashlib,json,os,re,subprocess,time,datetime
HERE=Path(__file__).resolve().parent;R=HERE.parents[2]
O=R/'results/retained_checks_20261113';LOCAL=R/'results/local-runs/retained_checks_20261113';LOG=R/'results/logs/retained_checks_20261113'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
def main():
    assert not O.exists() and not LOCAL.exists() and not LOG.exists()
    g=R/'results/retained_faces_20261113/summary.json';assert read(g)['fixed_geometry_matrix_admitted']
    native=read(R/'results/uniform5_pilot_20261113/protocol.json');assert all(sha(R/f)==h for f,h in native['source_hashes'].items()) and sha(Path(native['executable']))==native['binary_sha256']
    O.mkdir();LOCAL.mkdir();LOG.mkdir()
    files=[Path(__file__),HERE/'PLAN.md',HERE/'xmake.lua',R/'thin_film/include/vf/thinfilm/VoronoiRetainedFacesCandidate.h']
    names=['curved','preparation','observer']
    for name in names:files.append(R/('thin_film/tests/test_retained_'+name+'.cpp'));files.append(R/('thin_film/tests/test_local_assembly_'+name+'.cpp'))
    files += [R/'thin_film/tests'/name for name in ['test_curved_edge.cpp','test_edge_preparation.cpp','test_curved_mechanical_observer.cpp']]
    save(O/'protocol.json',dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources={str(p.relative_to(R)):sha(p) for p in files},geometry_summary_sha256=sha(g),original_sources=native['source_hashes'],original_binary_sha256=native['binary_sha256'],tests=names,VF_HAVE_CHOLMOD=True,VF_THREADS=2,scope='Existing fixtures under new geometry selection, not full C1 or actual78 replay'))
    env={k:v for k,v in os.environ.items() if not k.startswith('VF_')};env['VF_THREADS']='2';env['NO_COLOR']='1';attempts=[];xmake=R.parents[1]/'.tools/xmake/xmake/xmake.exe';sdk=R.parents[1]/'.tools/msvc'
    commands=[('configure',[str(xmake),'f','-y','-P',str(HERE),'-p','windows','-a','x64','-m','release','--sdk='+str(sdk)],300),('build',[str(xmake),'build','-y','-P',str(HERE),'-j2'],600)]+[(name,[str(xmake),'run','-P',str(HERE),'test_retained_'+name],120) for name in names]
    save(O/'status.json',dict(state='running',attempts=[],full_C1_complete=False))
    for name,cmd,cap in commands:
        log=LOG/(name+'.log');start=time.monotonic()
        with log.open('w',encoding='utf-8') as stream:
            try:rc=subprocess.run(cmd,cwd=R,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=cap).returncode
            except subprocess.TimeoutExpired:rc='timeout'
        a=dict(name=name,command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log));attempts.append(a)
        if name in names:
            counts=re.findall(r'(\d+) passed, (\d+) failed',log.read_text(encoding='utf-8'));a['counts']=counts;a['passed']=rc==0 and bool(counts) and all(int(f)==0 for _,f in counts);exe=LOCAL/'bin'/('test_retained_'+name+'.exe');a['binary_sha256']=sha(exe)
        save(O/'attempts.json',attempts)
        if name in ['configure','build']:assert rc==0,'Build failure retained, no retry'
    assert all(sha(R/f)==h for f,h in native['source_hashes'].items()) and sha(Path(native['executable']))==native['binary_sha256']
    tests=[a for a in attempts if a['name'] in names];save(O/'summary.json',dict(new_unit_checks=sum(int(n) for a in tests for n,_ in a['counts']),all_fixture_checks_passed=all(a['passed'] for a in tests),tests=tests,new_native_trajectory_replays=0,original_capture_unchanged=True,full_C1_complete=False,scope='Existing curved fixtures only; near-topology and actual dynamic acceptance still required'))
    save(O/'status.json',dict(state='completed',full_C1_complete=False));print((O/'summary.json').read_text(encoding='utf-8'))
if __name__=='__main__':main()
