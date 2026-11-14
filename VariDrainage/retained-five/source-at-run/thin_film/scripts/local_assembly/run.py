"""Candidate controls before actual dynamics; original baseline immutable."""
from pathlib import Path
import importlib.util,json,subprocess,re,csv,math
from equivalence import check
spec=importlib.util.spec_from_file_location('base',Path(__file__).parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/local_assembly_20261113';RAW=R/'results/local-runs/local_assembly_20261113';LOG=R/'results/logs/local_assembly_20261113';save,sha=base.save,base.sha
def inventory(root):return [dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(root.iterdir()) if p.is_file()]
def main():
    for p in [O,RAW,LOG]:p.mkdir(parents=True,exist_ok=True)
    assert not (O/'protocol.json').exists();base.LOG=LOG
    baseline=json.loads((R/'results/drainage_cost_20261113/status.json').read_text());assert baseline['state']=='completed' and all(c['passed'] for c in baseline['cases'])
    old=json.loads((R/'results/global_five_20261113/status.json').read_text());old2=json.loads((R/'results/global_space_20261113/status.json').read_text());refs=old['cases']+[next(c for c in old2['cases'] if c['name']=='u10_dt00625')]
    for c in refs:
        assert c['native_gate'] and c['completed']
        for f in c['files']:assert sha(R/f['path'])==f['sha256']
    expected=json.loads((R/'results/drainage_cost_20261113/protocol.json').read_text())['original_binary_sha256'];original=base.BIN/'tf_curved_edge_bench.exe';assert sha(original)==expected
    targets=['test_local_curved_hessian','test_local_assembly_curved','test_local_assembly_preparation','test_local_assembly_observer','tf_local_assembly_probe','tf_local_assembly_bench'];sources=list((R/'thin_film/include/vf/thinfilm').glob('*.h'))+[p for p in (R/'thin_film/bench').glob('*.cpp')]+[p for p in (R/'thin_film/tests').glob('test_local*.cpp')]+[Path(__file__),Path(__file__).with_name('equivalence.py'),Path(__file__).parents[1]/'drainage_cost/LOCAL_ASSEMBLY.md',R/'thin_film/xmake.lua']
    save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),sources={str(p.relative_to(R)):sha(p) for p in sources},binaries={n:sha(base.BIN/(n+'.exe')) for n in targets},original_binary_sha256=expected,reference_cases=refs,operator_inputs=baseline['cases'],whole_native_budget_s=3600,operator_budget_s=600,actual_budget_s=1800,environment={'VF_THREADS':'2','caller_VF':'excluded'},scope='Opt-in local energy Hessian assembly, original fields and dense H1 backend unchanged; no default adoption/full C1 claim'))
    state=dict(state='running_tests',tests=[],operators=[],cases=[],skipped=[]);save(O/'status.json',state);used=0
    def run(name,cmd,cap):
        nonlocal used
        assert used<3600;rec=base.execute(name,cmd,min(cap,3600-used));used+=rec['wall_s'];return rec
    for name,count in zip(targets[:4],[32,14,4,6]):
        rec=run(name,[str(base.BIN/(name+'.exe'))],120);counts=re.findall(r'(\d+) passed, (\d+) failed',(R/rec['log']).read_text(encoding='utf-8'));rec['counts']=counts;rec['passed']=rec['returncode']==0 and bool(counts) and int(counts[-1][0])>=count and counts[-1][1]=='0';state['tests'].append(rec);save(O/'status.json',state);print(name,rec['passed'],counts,flush=True)
    ready=all(t['passed'] for t in state['tests']);operator_used=0
    if ready:
        state['state']='running_operators';save(O/'status.json',state)
        for c in baseline['cases']:
            if operator_used>=600:state['skipped'].append(dict(name='operator_'+c['name'],reason='operator budget'));continue
            dest=RAW/('operator_'+c['name']);assert not dest.exists();dest.mkdir();src=R/'results/local-runs/drainage_cost_20261113'/c['name']/'input.csv';assert sha(src)==next(f['sha256'] for f in c['files'] if Path(f['path']).name=='input.csv');inp=dest/'input.csv';inp.write_bytes(src.read_bytes());out=dest/'comparison.json'
            rec=run('operator_'+c['name'],[str(base.BIN/'tf_local_assembly_probe.exe'),str(inp),str(out)],min(240,600-operator_used));operator_used+=rec['wall_s'];rec['passed']=False
            if out.exists():
                rec['comparison']=json.loads(out.read_text());d=rec['comparison'];rec['passed']=rec['returncode']==0 and d['first_fields_exact'] and max(d[k] for k in ['matrix_relative_error','symmetry_relative_error','diagonal_relative_error','action_relative_error'])<=2e-13
            rec['files']=inventory(dest);state['operators'].append(rec);save(O/'status.json',state);print(rec['name'],rec['passed'],round(rec['wall_s'],3),flush=True)
    ready &= len(state['operators'])==4 and all(c['passed'] for c in state['operators']);actual_used=0
    state['state']='running_actual' if ready else 'stopped_controls';save(O/'status.json',state)
    for ref in refs:
        if not ready or actual_used>=1800:state['skipped'].append(dict(name=ref['name'],reason='prior control/equivalence failed or actual budget'));continue
        name=ref['name'];dest=RAW/name;assert not dest.exists();dest.mkdir();cmd=ref['command'].copy();cmd[0]=str(base.BIN/'tf_local_assembly_bench.exe');cmd[7]=str(dest);rec=run(name,cmd,min(900 if name=='u10_dt00625' else 300,1800-actual_used));actual_used+=rec['wall_s'];rec.update(parameters=ref['parameters'],completed=False,native_gate=False,passed=False)
        if (dest/'summary.json').exists():
            m=json.loads((dest/'summary.json').read_text());p=json.loads((dest/'preparation.json').read_text());rec.update(summary=m,preparation=p,mesh=json.loads((dest/'mesh.json').read_text()));rec['completed']=rec['returncode']==0 and m['steps']==ref['summary']['steps'] and m['failures']==0 and m['N']==ref['summary']['N']
            if rec['completed']:
                stats=list(csv.DictReader((dest/'mechanics.csv').open()));ratio=max(float(r['g1_N'])/float(r['target_N']) for r in stats);rec['maximum_force_target_ratio']=ratio;rec['native_gate']=len(stats)==m['steps'] and ratio<=1.000001 and max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and m['max_budget_excess']<=1e-7 and p['pressure_residual']<=1e-8
        rec['files']=inventory(dest);reference_root=R/'results/local-runs'/('global_space_20261113' if name=='u10_dt00625' else 'global_five_20261113')/name;expected_files={Path(f['path']).name:f['sha256'] for f in ref['files']};actual_files={Path(f['path']).name:f['sha256'] for f in rec['files']};rec['same_native_file_set']=set(actual_files)==set(expected_files);rec['byte_matches']={n:actual_files.get(n)==digest for n,digest in expected_files.items()}
        if rec['completed'] and rec['same_native_file_set']:
            eq=check(reference_root,dest,ref['summary']['dx_m']);save(O/('equivalence_'+name+'.json'),eq);rec['equivalence_passed']=eq['passed'];rec['passed']=rec['native_gate'] and eq['passed']
        state['cases'].append(rec);save(O/'status.json',state);print(name,'native',rec['native_gate'],'equiv',rec.get('equivalence_passed'),round(rec['wall_s'],3),flush=True);ready &= rec['passed']
    state.update(state='completed',compute_wall_s=used,original_binary_unchanged=sha(original)==expected);save(O/'status.json',state);assert state['original_binary_unchanged']
if __name__=='__main__':main()
