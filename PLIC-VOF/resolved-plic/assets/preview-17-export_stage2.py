"""Reuse verified stage1 exporter with isolated stage2 plan,paths and identities."""
from pathlib import Path
import sys,os
sys.dont_write_bytecode=True
import hashlib,importlib.util,json,py_compile,shutil,time,traceback,types
R=Path(__file__).resolve().parent;PARENT=R.parent;W=PARENT.parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2))
    for i in range(40):
        try:temp.replace(path);return
        except PermissionError:
            if i==39:raise
            time.sleep(.05)
def setup_aliases():
    pairs=[(PARENT/'data/validation.npz',R/'data/validation.npz'),(PARENT/'compact_models.py',R/'compact_models.py')]
    source=PARENT/'native_validation'
    pairs.extend((source/name,R/'native_validation'/name) for name in ['summary.json','libaphros_production.so','production_build.json','production_parity.json'])
    pairs.extend((p,R/'native_validation'/p.name) for p in source.glob('*_predictions.npz'))
    entries=[]
    for src,dest in pairs:
        assert src.exists(),src;dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():assert sha(dest)==sha(src),('Immutable stage2 alias differs',dest)
        elif src.suffix=='.py':shutil.copy2(src,dest)
        else:
            try:os.link(src,dest)
            except OSError:shutil.copy2(src,dest)
        entries.append(dict(source=str(src.relative_to(W)),alias=str(dest.relative_to(W)),sha256=sha(src),method='hard_link' if os.path.samefile(src,dest) else 'identical_copy',use='read-only input;never modified by stage2 postprocessing'))
    target=R/'input_aliases.json'
    metadata=dict(stage='stage2_courant_objective',workspace=str(W),parent_stage_inputs_read_only=True,entries=entries)
    if not target.exists() or json.loads(target.read_text())!=metadata:save(target,metadata)
    return metadata
def plan_records():
    path=R/'train_stage2.py'
    spec=importlib.util.spec_from_file_location('_stage2_plan_readonly',path);module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    rows=module.records();assert len(rows)==8 and len({r['group'] for r in rows})==8
    assert sum(r['operator']=='flux' for r in rows)==6 and sum(r['operator']=='curvature' for r in rows)==2
    return rows
def load_exporter():
    base=PARENT/'export_resolved.py';text=base.read_text()
    # Identity-only suffix replacement;all math/export/parity code is reused.
    assert text.count('_resolved_r1')==2
    text=text.replace('_resolved_r1','_stage2_r1')
    module=types.ModuleType('_stage2_export_base');module.__file__=str(base);sys.modules[module.__name__]=module
    exec(compile(text,str(base),'exec'),module.__dict__)
    module.R=R;module.W=W;module.OUT=R/'inference';module.write=save;module.planned_records=plan_records
    def prepare_only():
        torch=module.torch;assert not torch.cuda.is_initialized();torch.set_num_threads(1)
        if not (R/'train_stage2.py').exists():
            module.state(state='prepared_waiting_trainer_implementation',cuda_used=False);return
        sources=[Path(__file__),R/'train_stage2.py',R/'compact_models.py',base,W/'operator-ablation-20261122/03-data-stencil/models.py']
        for path in sources:py_compile.compile(str(path),cfile=str(module.OUT/'cpu_engineering'/(path.stem+'.pyc')),doraise=True)
        checks=[];seen=set()
        for record in plan_records():
            config=record['config'];key=json.dumps(config,sort_keys=True)
            if key in seen:continue
            seen.add(key);torch.manual_seed(9173);model=module.build_cpu(config);k=config['stencil_size']
            args=[]
            for b in [17,79]:
                x=torch.rand(b,k,k);c=torch.rand(b);args.append((x,c) if record['operator']=='flux' else (x,))
            with torch.inference_mode():
                compiled=torch.jit.trace(model,args[0],check_inputs=[args[1]],strict=True)
                error=max(float((model(*a)-compiled(*a)).abs().max()) for a in args)
            assert error<=module.ABS_TOL
            checks.append(dict(config=config,cpu_trace_max_abs=error,batches=[17,79],passed=True))
        assert not torch.cuda.is_initialized()
        save(module.OUT/'cpu_compile_check.json',dict(passed=True,cuda_used=False,architectures_checked=len(checks),checks=checks,stage='stage2'))
        module.state(state='prepared_waiting_stage2_selection',cuda_used=False,architectures_checked=len(checks))
        print(json.dumps(dict(state='prepared_waiting_stage2_selection',architectures_checked=len(checks),cuda_used=False)))
    module.prepare_only=prepare_only
    save(R/'export_wrapper_metadata.json',dict(stage='stage2_courant_objective',wrapper_path=str(Path(__file__).relative_to(W)),wrapper_sha256=sha(__file__),base_script=str(base.relative_to(W)),base_script_sha256=sha(base),
        output_root=str(R.relative_to(W)),workspace=str(W),plan_source='train_stage2.py:records()',selection='selection_locked.json',
        identity_only_transform='_resolved_r1 model-id suffix replaced by _stage2_r1',numerical_export_code_unchanged=True,aliases='input_aliases.json'))
    return module
def main():
    setup_aliases();module=load_exporter();module.main()
if __name__=='__main__':
    try:main()
    except Exception:
        save(R/'export_status.json',dict(state='failed',pid=os.getpid(),traceback=traceback.format_exc()));raise
