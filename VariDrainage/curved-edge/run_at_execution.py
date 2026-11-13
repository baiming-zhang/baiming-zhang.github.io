"""Frozen, bounded full-state curved-edge SI matrix; preserve failures."""
from pathlib import Path
import importlib.util,json,time,subprocess,datetime
spec=importlib.util.spec_from_file_location('curved_native_base',Path(__file__).resolve().parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/curved_edge_20261113';RAW=O/'data';LOG=R/'results/logs/curved_edge_20261113';sha,save=base.sha,base.save
def main():
    O.mkdir(exist_ok=True);RAW.mkdir(exist_ok=True);LOG.mkdir(exist_ok=True);base.LOG=LOG;assert not (O/'protocol.json').exists();start=time.monotonic()
    cases=[]
    def add(name,mode='local',dx=40,dt=.0625,T=.5,radius=500,amp=.005):cases.append(dict(name=name,mode=mode,dx_um=dx,dt_us=dt,T_us=T,radius_um=radius,amplitude=amp))
    add('guard',dt=.125,T=.25,amp=0)
    for tag,dt in [('t125',.125),('t0625',.0625),('t03125',.03125)]:add('l40_'+tag,dt=dt)
    add('u40','uniform');add('l20',dx=20,dt=.03125);add('u20','uniform',20,.03125);add('r1000',radius=1000)
    paths=[Path(__file__),Path(__file__).with_name('README.md')]+[R/'thin_film'/p for p in ['include/vf/thinfilm/CurvedEdge.h','include/vf/thinfilm/CurvedEdgeMotion.h','include/vf/thinfilm/MovingVoronoiFilm.h','include/vf/thinfilm/Units.h','bench/tf_curved_edge_bench.cpp']]
    save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in paths},binaries={n:sha(base.BIN/(n+'.exe')) for n in ['test_curved_edge','tf_curved_edge_bench']},cases=cases,whole_budget_s=1800,per_case_timeout_s=180,environment={'VF_THREADS':'2','caller_VF':'excluded','DLL_parents':'local xmake packages'},started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    status=dict(state='running',tests=[],cases=[]);save(O/'status.json',status)
    rec=base.execute('test_curved_edge',[str(base.BIN/'test_curved_edge.exe')],120);status['tests'].append(rec);save(O/'status.json',status);assert rec['returncode']==0 and '11 passed, 0 failed' in (R/rec['log']).read_text()
    for c in cases:
        if time.monotonic()-start>1800:status['skipped']='whole batch budget';break
        dest=RAW/c['name'];assert not dest.exists();dest.mkdir();cmd=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest)]
        rec=base.execute(c['name'],cmd,min(180,1800-(time.monotonic()-start)));rec['parameters']=c
        if (dest/'summary.json').exists():rec['summary']=json.loads((dest/'summary.json').read_text())
        d=rec.get('summary',{});rec['completed']=rec['returncode']==0 and d.get('steps')==round(c['T_us']/c['dt_us']);rec['native_gate']=rec['completed'] and max(d[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and d['max_budget_excess']<=1e-7
        rec['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()];status['cases'].append(rec);save(O/'status.json',status)
        print(c['name'],rec['returncode'],rec['native_gate'],f"{rec['wall_s']:.3f}s",flush=True)
        if c['name']=='guard' and not rec['native_gate']:status['skipped']='failed startup guard';break
    status.update(state='completed',wall_s=time.monotonic()-start);save(O/'status.json',status)
if __name__=='__main__':main()
