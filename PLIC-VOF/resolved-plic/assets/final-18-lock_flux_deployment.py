"""Freeze the single promising flux deployment before independent data exist."""
from pathlib import Path
import datetime, hashlib, json

R=Path(__file__).resolve().parent;W=R.parent
GROUP='flux_compact_s5_w64_b16__sixteen_courants'
SEEDS=[4101,4102,4103]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    out=R/'deployment_selection_locked.json'
    if out.exists():
        print('Existing immutable deployment lock retained');return
    assert not (R/'flux_blind_test/datasets/manifest.json').exists()
    assert not (R/'flux_blind_test/large_grid/manifest.json').exists()
    adapter=R/'stage3_runtime_symmetry/runtime_deployment_adapter.py'
    core=R/'stage3_runtime_symmetry/runtime_cuda_graph.py'
    assert adapter.exists(),'Wait for reviewed frozen graph adapter'
    manifest=R/'stage2_courant_objective/inference/model_manifest.json'
    all_models=json.loads(manifest.read_text())['models']
    selected=sorted([m for m in all_models if m['group']==GROUP],key=lambda m:m['seed'])
    assert [m['seed'] for m in selected]==SEEDS and all(m['status']=='ready' for m in selected)
    complete=json.loads((R/'stage3_runtime_symmetry/runtime_completion.json').read_text())
    assert complete['state']=='complete' and complete['models']==6
    completion_models={m['model_id']:m for m in complete['identity']['frozen_models']}
    frozen=[]
    for m in selected:
        p=W/m['torchscript_path'];assert sha(p)==m['torchscript_sha256']==completion_models[m['model_id']]['torchscript_sha256']
        frozen.append(dict(model_id=m['model_id'],seed=m['seed'],torchscript_path=m['torchscript_path'],torchscript_sha256=sha(p),checkpoint_path=m['checkpoint_path'],checkpoint_sha256=sha(W/m['checkpoint_path'])))
    ids=[m['model_id'] for m in selected]
    obj=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=-7))).isoformat(),operators=['flux'],validation_only=True,
        selection='Frozen5x5 sixteen-C candidate from validation accuracy and provisional development graph timing;all three seeds retained. No final-test selection.',
        model_manifest_path=str(manifest.relative_to(W)),model_manifest_sha256=sha(manifest),model_ids=ids,model_runtimes={mid:'cuda_graph' for mid in ids},
        runtime_adapter_path=str(adapter.relative_to(W)),runtime_adapter_sha256=sha(adapter),runtime_core_path=str(core.relative_to(W)),runtime_core_sha256=sha(core),
        chosen_groups={'flux':dict(group=GROUP,best_accuracy=GROUP,fastest_reaching_native_error=GROUP,most_accurate_within_native_time=GROUP,scout_only=True)},
        seeds=SEEDS,all_seeds_retained=True,models=frozen,test_protocol_path=str((R/'flux_final_protocol.json').relative_to(W)),test_protocol_sha256=sha(R/'flux_final_protocol.json'),
        formal_speed_evidence=False,final_test_generated=False,curvature='Not advanced: retained negative validation result,not claimed to match particles or HF.',
        deployment_policy='FP32,TF32off,unchanged weights,CUDAgraph for every seed;no padding;capture shape equals actual query count,never selected by test error/speed.')
    temp=out.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2));temp.replace(out)
    print(json.dumps(dict(state='locked',models=len(ids),operators=obj['operators'],new_test_generated=False)))

if __name__=='__main__':main()
