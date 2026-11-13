"""Frozen bounded native moving-interface matrix. Never overwrites an attempt."""
from pathlib import Path
import sys,json,time,subprocess,datetime
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'moving_voronoi'))
import run as base
R=base.R;O=R/'results/moving_interface_20261113';RAW=O/'data';LOG=R/'results/logs/moving_interface_20261113'
save,sha=base.save,base.sha
def main():
    start=time.monotonic();O.mkdir(exist_ok=True);RAW.mkdir(exist_ok=True);LOG.mkdir(exist_ok=True);base.LOG=LOG
    cases=[]
    def add(name,mode='wall_local',fine=100,dt=1.25,amp=.005,radius=0,cut=800):
        cases.append(dict(name=name,mode=mode,fine_um=fine,dt_us=dt,T_us=10,amplitude=amp,radius_um=radius,cut_um=cut))
    add('equilibrium',dt=2.5,amp=0)
    for tag,dt in [('t2p5',2.5),('t1p25',1.25),('t0p625',.625)]:add('closed_l100_'+tag,dt=dt)
    add('closed_u100','wall_uniform');add('closed_l50',fine=50,dt=.625);add('closed_u50','wall_uniform',50,.625)
    add('open_l100',radius=500);add('observer_disabled',cut=None);add('alternative_cut',cut=775)
    assert len(cases)==10 and not (O/'protocol.json').exists()
    paths=[Path(__file__),Path(__file__).with_name('README.md'),R/'thin_film/bench/tf_moving_voronoi_bench.cpp',R/'thin_film/include/vf/thinfilm/MovingInterfaceLedger.h']
    save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),cases=cases,source_hashes={str(p.relative_to(R)):sha(p) for p in paths},binary_sha256=sha(base.BIN/'tf_moving_voronoi_bench.exe'),environment={'VF_THREADS':'2','caller_VF':'excluded','native_dll_paths':'local xmake package parents'},started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),whole_budget_s=1800,per_case_timeout_s=180))
    state=dict(state='running',tests=[],cases=[]);save(O/'status.json',state)
    for name,count in [('test_moving_interface',10),('test_moving_voronoi',22)]:
        rec=base.execute(name,[str(base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
        assert rec['returncode']==0 and f'{count} passed, 0 failed' in (R/rec['log']).read_text(),rec
    for c in cases:
        dest=RAW/c['name'];assert not dest.exists();dest.mkdir()
        command=[str(base.BIN/'tf_moving_voronoi_bench.exe'),c['mode']]+[str(c[k]) for k in ['fine_um','dt_us','T_us','amplitude','radius_um']]+[str(dest)]
        if c['cut_um'] is not None:command.append(str(c['cut_um']))
        left=1800-(time.monotonic()-start)
        if left<1:state['skipped']='whole budget';break
        rec=base.execute(c['name'],command,min(180,left));rec['parameters']=c
        if (dest/'summary.json').exists():rec['summary']=json.loads((dest/'summary.json').read_text())
        rec['completed']=rec['returncode']==0 and rec.get('summary',{}).get('steps')==round(10/c['dt_us'])
        rec['files']=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()]
        state['cases'].append(rec);save(O/'status.json',state);print(c['name'],rec['returncode'],f"{rec['wall_s']:.2f}s",flush=True)
    state.update(state='completed',wall_s=time.monotonic()-start);save(O/'status.json',state)
if __name__=='__main__':main()
