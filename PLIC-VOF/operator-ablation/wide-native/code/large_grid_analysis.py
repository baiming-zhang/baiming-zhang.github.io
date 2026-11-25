"""Aggregate actual-grid GPU/CPU records; retain seeds and unreachable targets."""
from pathlib import Path
import argparse,csv,json,collections
import numpy as np
from large_grid_operator_fixtures import ROOT,R,save

def load(path):return json.loads(path.read_text())
def csvout(path,rows):
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def run(require_complete=False):
    out=R/'analysis';out.mkdir(exist_ok=True)
    cases=load(R/'fixture_manifest.json')['cases']
    cpu=[load(p) for p in sorted((R/'cpu_timings').glob('*.json'))]
    gpu=[load(p) for p in sorted((R/'gpu_timings').glob('*.json'))]
    complete=len(cpu)==96 and len(gpu)==992
    status=dict(state='complete' if complete else 'waiting_for_measurements',cpu_configurations=len(cpu),expected_cpu=96,gpu_configurations=len(gpu),expected_gpu=992)
    save(out/'status.json',status)
    if not complete:
        pre=[load(p) for p in sorted((R/'cpu_preflight').glob('*.json'))]
        rows=[]
        for rec in pre:
            m=rec['metrics'];key='physical_kappa' if rec['operator']=='curvature' else 'flux_density'
            rows.append(dict(scene=rec['case']['scene'],N=rec['case']['N'],operator=rec['operator'],backend=rec['backend'],interface_cells=m['samples'],coverage=m['coverage'],full_task_mae=m[key]['mae'],available_only_mae=m[key]['available_mae']))
        csvout(out/'native_preflight_accuracy.csv',rows)
        if require_complete:raise RuntimeError('Formal measurements incomplete; no complete deployment conclusion may be drawn')
        print(json.dumps(status));return
    assert all(r['formal'] for r in cpu+gpu)
    by_case=collections.defaultdict(list)
    for rec in cpu+gpu:by_case[(rec['case']['scene'],rec['case']['N'],rec['operator'])].append(rec)
    records=[];coverage=[]
    for (scene,n,op),entries in by_case.items():
        case=entries[0]['case'];source=np.load(ROOT/case['source']);metric='physical_kappa' if op=='curvature' else 'flux_density'
        truth=source['kappa_truth'] if op=='curvature' else source['flux_density_truth'];fraction=source['center_fraction'] if op=='flux' else None
        predictions={r['prediction']:np.load(ROOT/r['prediction'])['prediction'] for r in entries}
        common=np.ones(len(truth),dtype=bool)
        for pred in predictions.values():common&=np.isfinite(pred)
        coverage.append(dict(scene=scene,N=n,operator=op,original_cells=len(truth),all_methods_finite_cells=int(common.sum()),removed=int((~common).sum())))
        for rec in entries:
            isgpu=rec['device']=='GPU';m=rec['model'] if isgpu else {};pred=predictions[rec['prediction']]
            values=pred*n if op=='curvature' else fraction*pred/.5
            available=float(abs(values[common]-truth[common]).mean()) if common.any() else None
            records.append(dict(scene=scene,N=n,total_grid_cells=n*n,interface_cells=len(truth),operator=op,device=rec['device'],
                implementation=m['model_id'] if isgpu else rec['backend'],mode=rec.get('mode','cpu'),workers=rec.get('workers',''),
                policy=m.get('data_strategy',m.get('data_policy','historical')) if isgpu else '',stencil=m.get('stencil_size','') if isgpu else (3 if op=='flux' else 7),
                seed=m.get('seed','') if isgpu else '',historical=m.get('historical_reference',False) if isgpu else False,
                time_median_ms=rec['wall_ms']['median'],time_p05_ms=rec['wall_ms']['p05'],time_p95_ms=rec['wall_ms']['p95'],
                full_task_coverage=rec['metrics']['coverage'],full_task_mae=rec['metrics'][metric]['mae'],common_finite_mae=available,
                common_finite_cells=int(common.sum()),metric=metric,cpu_worker_ratio=rec.get('worker_cpu_wall_ratio_median',''),
                source_prediction=rec['prediction']))
    csvout(out/'all_measured_records.csv',records);csvout(out/'coverage.csv',coverage)
    grouped=collections.defaultdict(list)
    for rec in records:
        if rec['device']=='GPU':
            key=(rec['scene'],rec['N'],rec['operator'],rec['mode'],rec['policy'],rec['stencil'],rec['historical'])
            grouped[key].append(rec)
    groups=[]
    for (scene,n,op,mode,policy,width,historical),members in grouped.items():
        expected=1 if historical else 3;assert len(members)==expected,(scene,n,op,mode,policy,width,len(members))
        full=[v['full_task_mae'] for v in members];conditional=[v['common_finite_mae'] for v in members];times=[v['time_median_ms'] for v in members]
        groups.append(dict(scene=scene,N=n,operator=op,mode=mode,policy=policy,stencil=width,historical=historical,seeds=len(members),
            time_mean_ms=float(np.mean(times)),time_seed_sd_ms=float(np.std(times,ddof=1)) if len(times)>1 else None,
            time_min_seed_ms=min(times),time_max_seed_ms=max(times),coverage_min=min(v['full_task_coverage'] for v in members),
            full_task_mae=float(np.mean(full)) if all(v is not None for v in full) else None,
            full_task_seed_sd=float(np.std(full,ddof=1)) if len(full)>1 and all(v is not None for v in full) else None,
            common_finite_mae=float(np.mean(conditional)) if all(v is not None for v in conditional) else None,
            member_error_max=max(full) if all(v is not None for v in full) else None,
            member_ids=';'.join(v['implementation'] for v in members)))
    csvout(out/'gpu_three_seed_groups.csv',groups)
    comparisons=[];targets=[]
    for case in cases:
        for op in ['flux','curvature']:
            refs=[r for r in records if r['scene']==case['scene'] and r['N']==case['N'] and r['operator']==op and r['device']=='CPU']
            assert len(refs)==6
            for backend in ['aphros','simliquid']:
                candidates=[r for r in refs if r['implementation']==backend];best=min(candidates,key=lambda r:r['time_median_ms'])
                for mode in ['resident','cpu_gpu_cpu']:
                    gs=[g for g in groups if g['scene']==case['scene'] and g['N']==case['N'] and g['operator']==op and g['mode']==mode]
                    for g in gs:
                        comparisons.append(dict(scene=case['scene'],N=case['N'],interface_cells=case['active_interface_cells'],operator=op,mode=mode,
                            gpu_policy=g['policy'],gpu_stencil=g['stencil'],historical=g['historical'],seeds=g['seeds'],
                            cpu_backend=backend,cpu_fastest_tested_workers=best['workers'],cpu_ms=best['time_median_ms'],gpu_ms=g['time_mean_ms'],
                            cpu_to_gpu_time_ratio=best['time_median_ms']/g['time_mean_ms'],cpu_full_mae=best['full_task_mae'],gpu_full_mae=g['full_task_mae'],
                            cpu_coverage=best['full_task_coverage'],gpu_coverage=g['coverage_min']))
                    qualified=[g for g in gs if g['coverage_min']==1 and g['full_task_mae'] is not None]
                    ontime=[g for g in qualified if g['time_mean_ms']<=best['time_median_ms']]
                    target_error=best['full_task_mae'];onerror=[g for g in qualified if target_error is not None and g['full_task_mae']<=target_error]
                    same_time=min(ontime,key=lambda g:g['full_task_mae']) if ontime else None
                    same_error=min(onerror,key=lambda g:g['time_mean_ms']) if onerror else None
                    targets.append(dict(scene=case['scene'],N=case['N'],operator=op,mode=mode,cpu_backend=backend,cpu_workers=best['workers'],
                        cpu_ms=best['time_median_ms'],cpu_error=target_error,cpu_coverage=best['full_task_coverage'],
                        equal_time_reachable=same_time is not None,equal_time_gpu=same_time,
                        equal_error_reachable=same_error is not None,equal_error_gpu=same_error,
                        all_gpu_seeds_reach_error=(same_error['member_error_max']<=target_error) if same_error else False))
    csvout(out/'gpu_vs_fastest_tested_cpu.csv',comparisons)
    save(out/'equal_time_equal_error.json',dict(no_interpolation=True,all_seeds_retained=True,records=targets))
    save(out/'summary.json',dict(state='complete',scope='actual local operators on complete static VOF fields;not DNS/PDE time',
        precision='GPU float32 input/native ONNX arithmetic;CPU numerical kernels float64',
        primary_mode='cpu_gpu_cpu includes H2D and D2H;resident mode separate',
        fastest_cpu='Minimum measured time over1/8/24 fixed workers per backend;parallel variants are explicit adapters',
        operating_system='CPU thread/wall indicators and device telemetry retained;do not infer precision-normalized algorithmic speedup',
        cases=cases,coverage=coverage,measured_record_count=len(records),gpu_group_count=len(groups),target_count=len(targets)))
    print(json.dumps(status))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--require-complete',action='store_true');run(parser.parse_args().require_complete)
