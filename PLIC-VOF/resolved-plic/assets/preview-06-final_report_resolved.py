"""Final flux/negative-curvature report. CPU rendering, no test generation.

--prepare-only only compiles and publishes the input contract. No HTML or final
claim is generated until frozen independent accuracy AND new-grid formal time
are complete. Interim reports and every previous render remain untouched.
"""
from pathlib import Path
import argparse,csv,datetime,fcntl,hashlib,json,math,py_compile,shutil,statistics,sys,traceback
import report_resolved as H

P=Path(__file__).resolve().parent;W=P.parent;S2=P/'stage2_courant_objective';S3=P/'stage3_runtime_symmetry';B=P/'flux_blind_test';F=S3/'runtime_formal_results'
STATUS=P/'final_report_status.json';SEEDS=[4101,4102,4103];GROUP='flux_compact_s5_w64_b16__sixteen_courants'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def status(name,**extra):H.write(STATUS,dict(state=name,updated_utc=now(),gpu_used=False,**extra))
def done(d):return d.get('state')=='complete' or d.get('complete') is True
def csvrows(path):return H.csv_rows(path)
def f(row,key):return H.number(row.get(key))
def avg(values):return H.mean_sd(values)
def allscope(row):return row.get('shape')=='all' and row.get('radius_bin')=='all'
def nid(row):return row['model_id'].startswith('native_')
def field(record,*names):
    for name in names:
        if name in record:return record[name]
    raise KeyError('Required explicit field missing: '+str(names))
def pathof(path):
    path=str(path).replace('\\','/')
    if len(path)>2 and path[1:3]==':/':path='/mnt/'+path[0].lower()+path[2:]
    value=Path(path);return value if value.is_absolute() else W/value

def contract():
    return dict(stage2=dict(summary='stage2_courant_objective/validation_pareto/scout/summary.json',training='stage2_courant_objective/training_protocol.json',scope='development only;paired source/objective constant'),
        stage3=dict(curvature='stage3_runtime_symmetry/curvature_summary.json',flux='stage3_runtime_symmetry/runtime_summary.json',scope='finite negative curvature diagnostic and seen-fixture deployment discovery;not final speed'),
        frozen=dict(lock='deployment_selection_locked.json',protocol='flux_final_protocol.json',three_seeds=SEEDS,runtime='cuda_graph',no_post_test_tuning=True),
        accuracy=dict(completion='flux_blind_test/completion.json',metrics='flux_blind_test/metrics.csv',primary='exact_cell_transport_mae',
            required_subsets=['planned_full','resolved_primary','resolved_allN_common','resolved_common_finite_conditional','resolved_allN_common_finite_conditional'],
            physical_primary='local+conservative whole-structure radius>=4;coverage retained',convergence='fixed resolved_allN_common,not varying perN population',curvature_tested=False),
        formal=dict(completion='stage3_runtime_symmetry/runtime_formal_results/completion.json',summary='stage3_runtime_symmetry/runtime_formal_results/summary.json',
            cases='new seed400 single_drop/many_drops,N512..8192,actualquerycount',required=['cases','per_seed_error_time','native_points','groups','equal_time_equal_error'],
            scope='same frozen artifacts/adapter;90s idle and clean telemetry;all3seeds and fastest actual native1/4/8 retained',capture='capture_seconds includes12 eager warmups+graphcapture per actualshape;excluded from hotcalls;new querycount needsrecapture'),
        output='One versioned HTML with embedded PNG/PDF/source/tables;no oldreport overwrite;root visual review required')

def load_inputs(require_formal=True):
    paths=[P/'deployment_selection_locked.json',P/'flux_final_protocol.json',B/'completion.json',B/'metrics.csv',B/'resolution_coverage.json',B/'datasets/manifest.json',
        S2/'training_protocol.json',S2/'validation_pareto/scout/summary.json',S3/'curvature_summary.json',S3/'runtime_summary.json']
    if require_formal:paths.extend([F/'completion.json',F/'summary.json',B/'large_grid/manifest.json'])
    missing=[str(p.relative_to(P)) for p in paths if not p.exists()]
    if missing:return None,missing
    lock=read(paths[0]);assert lock['operators']==['flux'] and len(lock['model_ids'])==3 and sorted(lock['seeds'])==SEEDS
    assert set(lock['model_runtimes'].values())=={'cuda_graph'}
    blind=read(B/'completion.json');formal=None;large=None
    if not done(blind):missing.append('Independent flux accuracy incomplete')
    if require_formal:
        formal=read(F/'summary.json');formal_done=read(F/'completion.json');large=read(B/'large_grid/manifest.json')
        assert formal['identity']==formal_done['identity']
        assert formal['identity']['deployment_lock_sha256']==sha(paths[0])
        assert formal['identity']['new_grid_manifest_sha256']==sha(B/'large_grid/manifest.json')
        if not done(formal_done) or not done(formal):missing.append('New-grid formal time incomplete')
        if not done(large):missing.append('New independent true-grid fields incomplete')
        if not formal.get('formal_speed_evidence',formal.get('speed_evidence_eligible',False)):missing.append('Formal time has not passed its independent idle/contamination gate')
    if missing:return None,missing
    assert blind['curvature_tested'] is False and blind['identity']['deployment_lock_sha256']==sha(paths[0])
    assert {m['model_id'] for m in blind['identity']['models']}==set(lock['model_ids'])
    for m in lock['models']:
        assert sha(pathof(m['torchscript_path']))==m['torchscript_sha256']
        assert sha(pathof(m['checkpoint_path']))==m['checkpoint_sha256']
    assert sha(pathof(lock['runtime_adapter_path']))==lock['runtime_adapter_sha256']
    assert sha(pathof(lock['runtime_core_path']))==lock['runtime_core_sha256']
    rows=csvrows(B/'metrics.csv');required={'dataset','model_id','seed','N','shape','radius_bin','subset','samples','failures','coverage','resolution_coverage','exact_cell_transport_mae','flux_volume_mae','ratio_mae'}
    assert rows and required.issubset(rows[0])
    assert {r['model_id'] for r in rows if not nid(r)}==set(lock['model_ids'])
    assert all(int(float(r['N']))==-1 for r in rows if r['dataset']=='circle')
    s2=read(S2/'validation_pareto/scout/summary.json');assert s2['validation_only'] and not s2['speed_evidence_eligible']
    curv=read(S3/'curvature_summary.json');dev=read(S3/'runtime_summary.json')
    assert curv['validation_only'] and not curv['blind_test_generated'] and not curv['new_training']
    assert dev['development_grid_only'] and not dev['formal_speed_evidence']
    if require_formal:
        assert all(k in formal for k in ('cases','per_seed_error_time','native_points','groups','equal_time_equal_error'))
        for case in formal['cases']:
            assert case['scene'] in ('single_drop','many_drops') and int(case['N']) in (512,1024,2048,4096,8192)
            assert int(field(case,'actual_queries','batch','active_interface_cells'))>0 and case['genuine_field'] and case['minimum_radius_cells']>=4
        assert len(formal['cases'])==10
    return dict(lock=lock,blind=blind,rows=rows,formal=formal,large=large,s2=s2,curv=curv,dev=dev,required=paths),[]

def aggregate_rows(bundle,dataset,N,shape='all',subset='resolved_primary',radius='all'):
    rows=[r for r in bundle['rows'] if (r['dataset'],int(float(r['N'])),r['shape'],r['subset'],r['radius_bin'])==(dataset,N,shape,subset,radius)]
    nn=sorted([r for r in rows if not nid(r)],key=lambda r:int(r['seed']));native=[r for r in rows if nid(r)]
    assert len(nn)==3 and [int(r['seed']) for r in nn]==SEEDS and len(native)==2,(dataset,N,shape,subset,radius)
    e=[f(r,'exact_cell_transport_mae') for r in nn];mean,sd=avg(e)
    return dict(nn=nn,native=native,mean=mean,sd=sd,values=e,samples=int(nn[0]['samples']),resolution_coverage=f(nn[0],'resolution_coverage'),
        coverage=min(f(r,'coverage') or 0 for r in rows),failures=sum(int(r['failures']) for r in rows))
def aphros(record):return next(r for r in record['native'] if 'aphros' in r['model_id'])

def stage2_pairs(s2):
    protocol=read(S2/'training_protocol.json');result=[]
    chosen=[g for g in s2['groups'] if g['batch']==65536 and g['mode']=='cpu_gpu_cpu']
    for source in sorted({g['source_group'] for g in protocol['groups'] if g['operator']=='flux'}):
        a=next(g for g in protocol['groups'] if g['source_group']==source and g['condition']=='single_courant');b=next(g for g in protocol['groups'] if g['source_group']==source and g['condition']=='sixteen_courants')
        assert a['sources']==b['sources'] and a['config']==b['config']
        x=next(g for g in chosen if g['group']==a['group']);y=next(g for g in chosen if g['group']==b['group'])
        paired=[]
        for seed in SEEDS:
            aa=next(r['error'] for r in s2['per_seed_error_time'] if (r['group'],r['seed'],r['batch'],r['mode'])==(a['group'],seed,65536,'cpu_gpu_cpu'))
            bb=next(r['error'] for r in s2['per_seed_error_time'] if (r['group'],r['seed'],r['batch'],r['mode'])==(b['group'],seed,65536,'cpu_gpu_cpu'))
            paired.append(100*(1-bb/aa))
        result.append(dict(config=a['config'],single=x,multi=y,reduction=100*(1-y['error_mean']/x['error_mean']),paired=paired))
    return result

def formal_native_time(row):
    if 'time_ms' in row:return row['time_ms']
    if 'median_ms' in row:return row['median_ms']
    return row['wall_ms']['median']
def formal_seed_time(row):return field(row,'median_ms','time_ms') if 'wall_ms' not in row else row['wall_ms']['median']
def formal_error(row):
    if 'exact_cell_transport_mae' in row:return row['exact_cell_transport_mae']
    if isinstance(row.get('metrics'),dict) and 'exact_cell_transport_mae' in row['metrics']:return row['metrics']['exact_cell_transport_mae']
    # Accept the explicit primary error only when its definition is supplied.
    metric=str(row.get('metric',row.get('primary_metric','')))
    if metric in ('exact_cell_transport_mae','exact_normalized_cell_transport_mae'):return field(row,'error')
    raise KeyError('Formal primary exact cell transport error must be explicit')
def formal_groups(bundle):
    summary=bundle['formal'];out=[]
    for case in summary['cases']:
        for mode in ('resident','cpu_gpu_cpu'):
            points=[r for r in summary['per_seed_error_time'] if (r['scene'],int(r['N']),r['mode'])==(case['scene'],int(case['N']),mode)]
            assert len(points)==3 and sorted(int(r['seed']) for r in points)==SEEDS
            assert all(r.get('formal',False) for r in points)
            native=[r for r in summary['native_points'] if (r['scene'],int(r['N']))==(case['scene'],int(case['N']))]
            assert len(native)==6 and all(r.get('formal',False) for r in native)
            fastest=min(native,key=formal_native_time);ts=[formal_seed_time(r) for r in points];es=[formal_error(r) for r in points]
            tm,td=avg(ts);em,ed=avg(es);ne=formal_error(fastest);nt=formal_native_time(fastest)
            finite=all(e is not None for e in es) and ne is not None
            success=finite and max(es)<=ne and max(ts)<=nt and (max(es)<ne or max(ts)<nt)
            captures=[r.get('capture_seconds') for r in points]
            cold=[r.get('cold_load_and_capture_seconds') for r in points]
            assert all(v is not None and v>=0 for v in captures),'Recorded warm capture cost required'
            assert all(v is not None and v>=0 for v in cold),'Recorded full adapter constructor cost required'
            amortization=[math.ceil(1000*c/(nt-t)) if success and nt>t else None for c,t in zip(captures,ts)]
            cold_amortization=[math.ceil(1000*c/(nt-t)) if success and nt>t else None for c,t in zip(cold,ts)]
            out.append(dict(case=case,scene=case['scene'],N=int(case['N']),mode=mode,points=points,native=native,fastest=fastest,time_mean=tm,time_sd=td,error_mean=em,error_sd=ed,native_time=nt,native_error=ne,
                all_seeds_pareto=success,capture_seconds=captures,cold_constructor_seconds=cold,break_even_reuses=amortization,cold_break_even_reuses=cold_amortization))
    return out

def plot_shapes(bundle,dest):
    plt=H.plot_style();from matplotlib.lines import Line2D
    fig,axes=plt.subplots(2,3,figsize=(14.2,8));fig.subplots_adjust(left=.08,right=.79,bottom=.1,top=.985,wspace=.4,hspace=.35)
    levels=sorted({int(float(r['N'])) for r in bundle['rows'] if r['dataset']=='physical'})
    for col,shape in enumerate(('circle','ellipse','weak_sine')):
        records=[aggregate_rows(bundle,'physical',n,shape) for n in levels]
        for j,seed in enumerate(SEEDS):
            values=[r['values'][j] if r['values'][j] is not None else float('nan') for r in records]
            axes[0,col].plot(levels,values,'o-',color='#9ecae1',alpha=.55,linewidth=.6,markersize=2.5)
        means=[r['mean'] if r['mean'] is not None else float('nan') for r in records]
        axes[0,col].plot(levels,means,'o-',color='#08519c',linewidth=1.3,markersize=4)
        for backend,color,marker in [('aphros','#263e53','v'),('simliquid','#6c8ca3','<')]:
            vals=[f(next(n for n in r['native'] if backend in n['model_id']),'exact_cell_transport_mae') for r in records]
            axes[0,col].plot(levels,[v if v is not None else float('nan') for v in vals],marker=marker,linewidth=.9,color=color,markerfacecolor='none',markersize=4)
        ratios=[r['mean']/f(aphros(r),'exact_cell_transport_mae') if r['mean'] is not None and f(aphros(r),'exact_cell_transport_mae') else float('nan') for r in records]
        axes[1,col].plot(levels,ratios,'o-',color='#2171b5',linewidth=1.2,markersize=4);axes[1,col].axhline(1,color='#6f8ca4',linestyle='--',linewidth=.8)
        for row in (0,1):
            ax=axes[row,col];ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks(levels);ax.set_xticklabels([str(n) for n in levels],rotation=45,ha='right');ax.set_xlabel('Grid N');H.box(ax)
        axes[0,col].set_ylabel('Exact cell-transport MAE');axes[1,col].set_ylabel('NN / Aphros PLIC MAE')
    handles=[Line2D([],[],color='#08519c',marker='o',label='NN, mean of 3 seeds'),Line2D([],[],color='#9ecae1',marker='o',linewidth=.6,label='Individual NN seed'),Line2D([],[],color='#263e53',marker='v',markerfacecolor='none',label='Aphros PLIC'),Line2D([],[],color='#6c8ca3',marker='<',markerfacecolor='none',label='SimLiquid PLIC'),Line2D([],[],color='#2171b5',marker='o',label='NN / PLIC ratio'),Line2D([],[],color='#6f8ca4',linestyle='--',label='Equal error')]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.805,.52),frameon=False,fontsize=8)
    files=dest/'blind_by_shape.png',dest/'blind_by_shape.pdf';fig.savefig(files[0],bbox_inches='tight',pad_inches=.08);fig.savefig(files[1],bbox_inches='tight',pad_inches=.08);plt.close(fig);return files

def plot_common(bundle,dest):
    plt=H.plot_style();fig,ax=plt.subplots(figsize=(10,4.8));fig.subplots_adjust(left=.1,right=.73,bottom=.18,top=.98)
    levels=sorted({int(float(r['N'])) for r in bundle['rows'] if r['dataset']=='physical'});records=[aggregate_rows(bundle,'physical',n,subset='resolved_allN_common') for n in levels]
    assert len({r['samples'] for r in records})==1,'Convergence plot must keep the fixed crossN support'
    for j in range(3):ax.plot(levels,[r['values'][j] for r in records],'o-',color='#b8d5eb',linewidth=.6,markersize=3)
    ax.errorbar(levels,[r['mean'] for r in records],yerr=[r['sd'] for r in records],fmt='o-',color='#08519c',capsize=2,label='NN, mean +/- seed SD')
    for backend,color,marker in [('aphros','#263e53','v'),('simliquid','#6c8ca3','<')]:
        ax.plot(levels,[f(next(n for n in r['native'] if backend in n['model_id']),'exact_cell_transport_mae') for r in records],marker=marker,linewidth=.9,color=color,markerfacecolor='none',label=('Aphros' if backend=='aphros' else 'SimLiquid')+' PLIC')
    ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks(levels);ax.set_xticklabels([str(n) for n in levels],rotation=40,ha='right');ax.set_xlabel('Grid N, fixed cross-N support');ax.set_ylabel('Exact cell-transport MAE');H.box(ax);ax.legend(loc='center left',bbox_to_anchor=(1.025,.5),frameon=False,fontsize=8)
    files=dest/'fixed_support_refinement.png',dest/'fixed_support_refinement.pdf';fig.savefig(files[0],bbox_inches='tight',pad_inches=.08);fig.savefig(files[1],bbox_inches='tight',pad_inches=.08);plt.close(fig);return files

def plot_formal_ratios(formal,dest):
    plt=H.plot_style();from matplotlib.lines import Line2D
    fig,axes=plt.subplots(2,2,figsize=(12.2,8));fig.subplots_adjust(left=.09,right=.76,bottom=.1,top=.985,wspace=.32,hspace=.3)
    for row,scene in enumerate(('single_drop','many_drops')):
        for mode,color,marker in [('resident','#8cb9db','s'),('cpu_gpu_cpu','#08519c','o')]:
            records=sorted([r for r in formal if (r['scene'],r['mode'])==(scene,mode)],key=lambda r:r['N']);ns=[r['N'] for r in records]
            for column in (0,1):
                values=[]
                for r in records:
                    values.append(r['native_time']/r['time_mean'] if column==0 else (r['error_mean']/r['native_error'] if r['error_mean'] is not None and r['native_error'] else float('nan')))
                axes[row,column].plot(ns,values,marker=marker,color=color,linewidth=1,markersize=4)
            for r in records:
                # Three individual frozen-seed ratios are retained visually.
                for point in r['points']:
                    axes[row,0].plot(r['N'],r['native_time']/formal_seed_time(point),marker=marker,color=color,alpha=.4,markersize=2.4)
                    if formal_error(point) is not None and r['native_error']:
                        axes[row,1].plot(r['N'],formal_error(point)/r['native_error'],marker=marker,color=color,alpha=.4,markersize=2.4)
        for column in (0,1):
            ax=axes[row,column];ax.set_xscale('log',base=2);ax.set_yscale('log');ax.axhline(1,color='#648098',linewidth=.8,linestyle='--');ax.set_xticks([512,1024,2048,4096,8192]);ax.set_xticklabels(['512','1024','2048','4096','8192']);ax.set_xlabel('Grid N');ax.set_ylabel('Fastest native / NN hot time' if column==0 else 'NN / fastest native exact MAE');H.box(ax)
    handles=[Line2D([],[],color='#8cb9db',marker='s',label='GPU-resident'),Line2D([],[],color='#08519c',marker='o',label='CPU-GPU-CPU'),Line2D([],[],color='#648098',linestyle='--',label='Equal time / error')]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.78,.5),frameon=False,fontsize=8)
    files=dest/'independent_formal_deployment.png',dest/'independent_formal_deployment.pdf';fig.savefig(files[0],bbox_inches='tight',pad_inches=.08);fig.savefig(files[1],bbox_inches='tight',pad_inches=.08);plt.close(fig);return files

def reproducible_sources():
    return [Path(__file__),P/'report_resolved.py',P/'data_generation.py',P/'new_large_grid.py',P/'lock_flux_deployment.py',P/'start_formal_waiter.py',
        P/'flux_final_protocol.json',P/'flux_blind_protocol.json',P/'flux_blind_evaluation.py',P/'blind_evaluation.py',
        S2/'train_stage2.py',S2/'export_stage2.py',S2/'validation_stage2.py',S2/'compact_models.py',
        S3/'runtime_formal_waiter.py',S3/'runtime_deployment_adapter.py',S3/'runtime_cuda_graph.py',S3/'runtime_formal_plan.json',S3/'runtime_protocol.json',
        S3/'curvature_models.py',S3/'curvature_evaluate.py',S3/'curvature_protocol.json']

def common_support_table(records,levels):
    rows=[]
    for n,r in zip(levels,records):
        sim=next(v for v in r['native'] if 'simliquid' in v['model_id'])
        rows.append([n,r['samples'],H.pm(r['mean'],r['sd']),H.fmt(f(aphros(r),'exact_cell_transport_mae')),H.fmt(f(sim,'exact_cell_transport_mae'))])
    return H.detail('固定跨N共同支持的聚合误差简表',H.table(['N','固定共同查询','NN精确MAE ± 标准差','Aphros PLIC精确MAE','SimLiquid PLIC精确MAE'],rows))

def render_accuracy_preview(bundle,dest,reuse_figures=None):
    circle=aggregate_rows(bundle,'circle',-1);reference=f(aphros(circle),'exact_cell_transport_mae');gain=100*(1-circle['mean']/reference)
    high={shape:aggregate_rows(bundle,'physical',8192,shape) for shape in ('circle','ellipse','weak_sine')};overall=aggregate_rows(bundle,'physical',8192)
    # Explicit formula kept separate from display rounding.
    shape_gain={shape:100*(1-r['mean']/f(aphros(r),'exact_cell_transport_mae')) for shape,r in high.items()}
    overall_gain=100*(1-overall['mean']/f(aphros(overall),'exact_cell_transport_mae'))
    levels=sorted({int(float(r['N'])) for r in bundle['rows'] if r['dataset']=='physical'});common=[aggregate_rows(bundle,'physical',n,subset='resolved_allN_common') for n in levels]
    curv=bundle['curv'];raw=next(g for g in curv['groups'] if g['method']=='raw');sym=next(g for g in curv['groups'] if g['method']=='d4_phase_odd');native=curv['native_q_mae_target']
    symgain=100*(1-sym['q_mae_mean']/raw['q_mae_mean']);gap=100*(sym['q_mae_mean']/native-1)
    parts=['<header><p>独立准确度已完成；正式推理计时等待GPU空闲</p><h1>独立准确度预览：改善有范围，收敛仍未解决</h1><ol class="lead">',
        f'<li><b>新圆精度改善，但不是所有几何都赢。</b>32,768个新圆查询中，NN精确运输体积MAE为{H.pm(circle["mean"],circle["sd"])}，原生PLIC为{H.fmt(reference)}，低{gain:.2f}%。物理N8192时圆和椭圆误差分别为PLIC的{high["circle"]["mean"]/f(aphros(high["circle"]),"exact_cell_transport_mae"):.0f}倍和{high["ellipse"]["mean"]/f(aphros(high["ellipse"]),"exact_cell_transport_mae"):.0f}倍；弱正弦低{shape_gain["weak_sine"]:.1f}%。三形状总体仅低{overall_gain:.1f}%，均值掩盖形状差异。</li>',
        '<li><b>正式推理时间仍待测。</b>新独立真实网格的正式计时正在等待GPU空闲，此预览不填入旧开发场或共享负载试测的时间，不宣称精度—时间Pareto已胜出。固定查询数CUDA Graph的热调用与新形状首次捕获成本会在正式结果中分别列出。</li>',
        f'<li><b>收敛与曲率目标未解决。</b>固定跨N共同支持每层{common[0]["samples"]}个查询，NN误差由N{levels[0]}的{H.fmt(common[0]["mean"])}到N{levels[-1]}的{H.fmt(common[-1]["mean"])}，仍有明显地板。曲率D4＋相补平均只改善{symgain:.2f}%，仍比particles高{gap:.1f}%，因此曲率未进入盲测。</li></ol></header>',
        '<p>冻结5×5通量网络的三个随机种子全部保留。模型、运行方式、精度与CUDA Graph适配器均在新数据前锁定，没有根据盲测换权重、调整算术或重新挑种子。输入仅为原始VOF体积分数与Courant数；解析几何只供离线真值，无冻结平面SVD层或几何回退。</p>']
    shape_rows=[]
    for shape,r in high.items():shape_rows.append([{'circle':'圆','ellipse':'椭圆','weak_sine':'弱正弦'}[shape],H.pm(r['mean'],r['sd']),H.fmt(f(aphros(r),'exact_cell_transport_mae')),H.fmt(r['mean']/f(aphros(r),'exact_cell_transport_mae')),str(r['samples']),H.fmt(100*r['resolution_coverage'])+'%'])
    parts.append('<h2>按形状拆开，而不只报总体平均</h2><p>主误差为原始float32中心体积分数×预测通量比与精确扫掠面积之差的绝对值均值；单位为单元面积归一化，没有乘h²。以下N8192的原始查询均已解析，负结果没有剔除。</p>'+H.table(['物理形状','NN精确MAE ± 标准差','Aphros PLIC精确MAE','NN/PLIC误差倍数','解析查询','解析覆盖'],shape_rows))
    if reuse_figures:
        saved=read(reuse_figures/'inputs_used.json')['files'];current={str(p.relative_to(P)):sha(p) for p in bundle['required']}
        assert all(current.get(v['path'])==v['sha256'] for v in saved),'Only unchanged scientific inputs permit figure reuse'
        figures=[]
        for name in ('blind_by_shape','fixed_support_refinement'):
            pair=[]
            for ext in ('png','pdf'):
                source=reuse_figures/f'{name}.{ext}';target=dest/source.name;shutil.copyfile(source,target);assert sha(source)==sha(target);pair.append(target)
            figures.append(tuple(pair))
    else:figures=[plot_shapes(bundle,dest),plot_common(bundle,dest)]
    parts.append(H.figure_html(figures[0],'左至右为圆、椭圆、弱正弦。上排为精确运输体积误差，下排为NN／Aphros PLIC误差比，对数纵轴；比值高于1表示网络更差。每N解析支持不同，这一图不能单独证明收敛；所有三个神经随机种子、两套原生PLIC均保留。'))
    parts.append(H.figure_html(figures[1],f'固定所有N共同解析的{common[0]["samples"]}个查询，不随网格换几何或样本。图中不使用h²缩放，不宣称已获得收敛阶。'))
    parts.append(common_support_table(common,levels))
    circles=[]
    for radius,label in [('all','全部4≤R/h≤512'),('4_to_8','4≤R/h<8'),('8_to_32','8≤R/h<32'),('32_to_128','32≤R/h<128'),('128_to_512','128≤R/h≤512')]:
        r=aggregate_rows(bundle,'circle',-1,radius=radius);ne=f(aphros(r),'exact_cell_transport_mae')
        circles.append([H.esc(label),str(r['samples']),H.pm(r['mean'],r['sd']),H.fmt(ne),('是' if all(v is not None and v<=ne for v in r['values']) else '否')])
    parts.append(H.detail('独立圆分层，未把新圆称为物理细化',H.table(['半径范围','独立查询','NN精确MAE ± 标准差','Aphros PLIC精确MAE','三种子均更准'],circles)))
    cov=[]
    for n in levels:
        r=aggregate_rows(bundle,'physical',n);cov.append([n,r['nn'][0]['original_samples'],r['samples'],H.fmt(100*r['resolution_coverage'])+'%',aggregate_rows(bundle,'physical',n,subset='resolved_allN_common')['samples'],r['failures']])
    parts.append('<h2>解析覆盖与固定支持</h2><p>保留所有原始物理查询。主解析支持同时要求局部曲率半径和整条结构的保守最小曲率半径≥4个单元；低N未解析查询只作描述，不支持主结论。共同支持曲线固定原始样本，不用随N扩大的解析集合制造收敛。</p>'+H.table(['N','全部原始查询','主解析查询','主解析覆盖','跨N共同解析查询','预测失败合计'],cov))
    parts.append('<h2>正式部署时间：等待GPU空闲</h2><p>这里尚无合格正式点，不填开发速度，不计算正式速度胜出或复用盈亏点。正式计时将使用新种子400的完整真网格、实际查询数、两套原生PLIC的1／4／8逻辑CPU线程，以及相同冻结网络的三个随机种子。驻留与包含H2D／D2H的往返边界分别报告。</p><p>CUDA Graph热调用适用于已捕获的固定查询数，新批量需重新捕获。capture_seconds实际包含12次准备调用和捕获；完整适配器构造另记录。正式数据齐后才列实测捕获成本及符合精度条件的有限复用盈亏点。</p>')
    pair_rows=[];pairs=stage2_pairs(bundle['s2'])
    for pair in pairs:
        c=pair['config'];pair_rows.append([f"{c['stencil_size']}×{c['stencil_size']}、宽{c['width']}、{c['num_bins']}bin",H.pm(pair['single']['error_mean'],pair['single']['error_seed_sd']),H.pm(pair['multi']['error_mean'],pair['multi']['error_seed_sd']),f"{pair['reduction']:+.1f}%",' / '.join(f'{v:+.1f}%' for v in pair['paired'])])
    parts.append('<h2>第二阶段的同源因果配对</h2><p>下表仅在相同源权重、架构、训练数据、运输体积L1目标与新增更新预算下比较单C和16C。不能把阶段1到阶段2的“改目标＋加C”总变化都归因于多C；三个配对随机种子按4101／4102／4103列出。</p>'+H.table(['源架构','单C验证MAE ± 标准差','16C验证MAE ± 标准差','组均MAE下降','逐随机种子下降'],pair_rows))
    curv_rows=[];names={'raw':'原网络','d4':'D4空间平均','phase_odd':'相补奇对称平均','d4_phase_odd':'D4＋相补平均'}
    for g in curv['groups']:curv_rows.append([names[g['method']],H.pm(g['q_mae_mean'],g['q_mae_std']),f"{100*(1-g['q_mae_mean']/raw['q_mae_mean']):+.2f}%",f"{100*(g['q_mae_mean']/native-1):+.1f}%",g['views']['logical_source_predictions_per_sample'],('通过' if g['all_guard_pass'] else '未通过')])
    parts.append('<h2>曲率有限负结果，未进入盲测</h2><p>这些仅对同一纯神经网络做空间／相补变换平均，没有引入几何曲率或回退。误差仍未达到真实particles，因此停止扩展，未用视图数推算速度，也未声称解决物理曲率收敛。</p>'+H.table(['纯神经诊断','验证q MAE ± 标准差','相对原网络改善','比particles误差高','每样本原网络预测数','原分布保留条件'],curv_rows))
    detail_rows=[];detail_identities=[]
    for r in bundle['rows']:
        if r['subset']!='resolved_primary' or r['radius_bin']!='all' or r['shape']=='all':continue
        method='神经网络' if not nid(r) else ('Aphros PLIC' if 'aphros' in r['model_id'] else 'SimLiquid PLIC')
        shape={'all':'全部','circle':'圆','ellipse':'椭圆','weak_sine':'弱正弦'}.get(r['shape'],r['shape']);support='逐N解析' if r['subset']=='resolved_primary' else '跨N共同解析'
        detail_rows.append([('新圆' if r['dataset']=='circle' else '独立物理'),r['N'],shape,support,method,r['seed'] if not nid(r) else '确定性',r['samples'],r['failures'],H.fmt(100*(f(r,'coverage') or 0))+'%',H.fmt(f(r,'exact_cell_transport_mae')),H.fmt(f(r,'ratio_mae'))])
        detail_identities.append((r['dataset'],r['N'],r['shape'],r['model_id']))
    assert len(detail_identities)==len(set(detail_identities))
    parts.append(H.detail('实际形状逐N解析结果：各随机种子与两套原生方法',H.table(['集','N','形状','预定义支持','方法','随机种子','查询','失败','有效覆盖','精确运输MAE','通量比MAE'],detail_rows,True)))
    parts.append('<h2>可复核证据与代码</h2><p>PDF、完整分层CSV与代码内嵌于本HTML。文件名仅在下载附件中展示；正文不展示机器标识或校验串。共同有限支持误差是条件统计，CSV同时保留原始失败率与移除数量。</p>'+H.attachment(B/'metrics.csv','独立精度完整分层CSV'))
    for source in reproducible_sources():
        if source.exists():parts.append(H.attachment(source))
    counts=read(P/'training_counts.json');datasets=read(B/'datasets/manifest.json')['files']
    split_rows=[['原解析子集·训练',f"{counts['original_resolved_train']:,}局部几何",'同源训练／对照'],['原解析子集·验证',f"{counts['original_resolved_validation']:,}局部几何",'保留条件与权重选择'],['新增圆·训练',f"{counts['added_train']:,}母圆",f"{counts['added_train']*16:,}监督对＝母圆×16C（相关多标签）"],['新增圆·验证',f"{counts['added_validation']:,}母圆",f"{counts['added_validation']:,}主要单C验证查询"]]
    for d in datasets:
        if d['dataset']=='circle':split_rows.append(['冻结后独立新圆','32,768母圆',f"{d['samples']:,}查询"])
        else:split_rows.append(['冻结后独立物理集','48母形状',f"{d['samples']:,}查询＝48×9偏移×4变体×9网格（跨N相关）"])
    split_rows.append(['独立真场，正式时间待测','2个固定物理场景','2场景×5分辨率＝10网格实例，非10个独立随机场景'])
    split_rows.append(['曲率最终盲测','—','未运行，验证目标未达'])
    parts.append('<h2>数据划分</h2>'+H.table(['数据部分','独立母几何／场景数','查询／标签数与用途'],split_rows))
    css='body{font-family:Arial,"Microsoft YaHei",sans-serif;color:#19374d;max-width:1450px;margin:36px auto;padding:0 30px;line-height:1.65}h1{font-size:30px}h2{border-top:1px solid #cbddeb;padding-top:18px;margin-top:34px}.lead{background:#f0f6fb;border-left:4px solid #2171b5;padding:18px 22px 18px 44px;font-size:18px}.lead li{padding:8px 0}img{width:100%;height:auto}figure{margin:28px 0}figcaption{font-size:13px;color:#4d6e87}table{border-collapse:collapse;font-size:13px;width:100%;font-variant-numeric:tabular-nums}th,td{padding:9px 11px;border-bottom:1px solid #d6e5ef;white-space:nowrap;text-align:left}th{background:#eaf2f8;position:sticky;top:0}.tablebox{overflow:auto;max-height:650px;border:1px solid #cbddea;margin:12px 0 22px}tbody tr:nth-child(even){background:#f8fbfe}details{padding:12px 15px;border:1px solid #cadae8;margin:18px 0}summary{cursor:pointer}.download{display:inline-block;border:1px solid #cbddea;padding:5px 9px;margin:5px;color:#08519c;text-decoration:none}.filter{padding:8px;width:350px;max-width:85%;border:1px solid #acc5d9}'
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>独立准确度预览：正式计时等待中</title><style>'+css+'</style>'+''.join(parts)+'<script>function filterTable(el){let q=el.value.toLowerCase();el.nextElementSibling.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}</script></html>'
    target=P/'accuracy_preview.html' if not (P/'accuracy_preview.html').exists() else dest/'accuracy_preview.html';tmp=dest/'preview.tmp.html';tmp.write_text(H.naturalize_html(doc),encoding='utf-8');tmp.replace(target)
    H.write(dest/'accuracy_numerical_checks.json',dict(passed=True,formal_time_used=False,human_detail_rows=len(detail_rows),human_detail_unique_keys=len(set(detail_identities)),scientific_csv_rows=len(bundle['rows']),reused_unchanged_figures=bool(reuse_figures),circle=dict(mean=circle['mean'],sd=circle['sd'],native=reference,reduction_percent=gain),physical8192={k:dict(mean=v['mean'],native=f(aphros(v),'exact_cell_transport_mae'),ratio=v['mean']/f(aphros(v),'exact_cell_transport_mae'),reduction_percent=shape_gain[k]) for k,v in high.items()},physical_average_gain_percent=overall_gain,common=dict(samples_per_N=common[0]['samples'],N_min=levels[0],N_max=levels[-1],error_first=common[0]['mean'],error_last=common[-1]['mean']),curvature=dict(improvement_percent=symgain,above_native_percent=gap)))
    return target,figures

def render(bundle,dest):
    rows=bundle['rows'];circle=aggregate_rows(bundle,'circle',-1);nr=f(aphros(circle),'exact_cell_transport_mae');gain=100*(1-circle['mean']/nr)
    high={shape:aggregate_rows(bundle,'physical',8192,shape) for shape in ('circle','ellipse','weak_sine')};high_all=aggregate_rows(bundle,'physical',8192)
    sine_gain=100*(1-high['weak_sine']['mean']/f(aphros(high['weak_sine']),'exact_cell_transport_mae'));overall_gain=100*(1-high_all['mean']/f(aphros(high_all),'exact_cell_transport_mae'))
    levels=sorted({int(float(r['N'])) for r in rows if r['dataset']=='physical'});common=[aggregate_rows(bundle,'physical',n,subset='resolved_allN_common') for n in levels]
    formal=formal_groups(bundle);roundtrips=[r for r in formal if r['mode']=='cpu_gpu_cpu'];formal_wins=sum(r['all_seeds_pareto'] for r in roundtrips)
    curv=bundle['curv'];raw=next(g for g in curv['groups'] if g['method']=='raw');sym=next(g for g in curv['groups'] if g['method']=='d4_phase_odd');native=curv['native_q_mae_target'];symgain=100*(1-sym['q_mae_mean']/raw['q_mae_mean']);gap=100*(sym['q_mae_mean']/native-1)
    parts=['<header><p>R/h≥4 · 冻结5×5通量网络 · 三个随机种子全部保留</p><h1>局部神经算子：精度、部署时间与收敛范围</h1><ol class="lead">',
        f'<li><b>精度改善有明确范围。</b>新圆盲测32,768例，NN精确运输体积MAE {H.pm(circle["mean"],circle["sd"])}，PLIC为{H.fmt(nr)}，低{gain:.2f}%。但N8192物理圆、椭圆误差分别约为PLIC的{high["circle"]["mean"]/f(aphros(high["circle"]),"exact_cell_transport_mae"):.1f}倍和{high["ellipse"]["mean"]/f(aphros(high["ellipse"]),"exact_cell_transport_mae"):.0f}倍；弱正弦则低{sine_gain:.1f}%。三形状总体仅低{overall_gain:.1f}%，平均掩盖形状差异，不能泛称全部胜出。</li>',
        f'<li><b>正式部署比较按实际真场逐案成立。</b>新独立真场的10个CPU→GPU→CPU工作负载中，{formal_wins}/10个同时满足三个随机种子的精度与热调用时间不劣于最快实际原生参考。以下逐案保留未达标结果；这是固定查询数已捕获图的热推理，不是任意新形状的首次调用。</li>',
        f'<li><b>收敛与曲率仍未解决。</b>固定跨网格共同支持每N {common[0]["samples"]}例，通量MAE从N{levels[0]}的{H.fmt(common[0]["mean"])}到N{levels[-1]}的{H.fmt(common[-1]["mean"])}，仍有明显误差地板。曲率纯网络群对称平均只改善{symgain:.2f}%，仍比particles高{gap:.1f}%；曲率未进入最终盲测。</li></ol></header>',
        '<p>神经输入仅为VOF体积分数与Courant数；曲率诊断仅用VOF。没有冻结平面SVD层、外部几何特征、几何修正或回退。解析几何仅离线提供真值。最终冻结模型、三个随机种子、FP32精度、CUDA Graph与适配器在新测试前锁定，未根据测试换权重或推理方法。</p>']
    # Main negatives remain visible before aggregate development averages.
    shape_rows=[]
    for shape,record in high.items():
        ne=f(aphros(record),'exact_cell_transport_mae');shape_rows.append([{'circle':'圆','ellipse':'椭圆','weak_sine':'弱正弦'}[shape],H.pm(record['mean'],record['sd']),H.fmt(ne),H.fmt(record['mean']/ne),str(record['samples']),H.fmt(100*record['resolution_coverage'])+'%'])
    parts.append('<h2>按形状拆开的独立精度</h2><p>主误差为 |原始float32中心体积分数×预测通量比−精确扫掠面积| 的均值，不乘h²。下表N8192已100%解析，负结果不能用三形状总体平均掩盖。</p>'+H.table(['形状','NN MAE均值 ± 标准差','Aphros PLIC MAE','NN/PLIC误差倍数','解析样本','解析覆盖'],shape_rows))
    circle_bins=[]
    for radius,label in [('4_to_8','4≤R/h<8'),('8_to_32','8≤R/h<32'),('32_to_128','32≤R/h<128'),('128_to_512','128≤R/h≤512')]:
        a=aggregate_rows(bundle,'circle',-1,radius=radius);native_error=f(aphros(a),'exact_cell_transport_mae')
        passed=all(e is not None and e<=native_error for e in a['values']) if native_error is not None else False
        circle_bins.append([H.esc(label),str(a['samples']),H.pm(a['mean'],a['sd']),H.fmt(native_error),('是' if passed else '否')])
    parts.append(H.detail('独立圆盲测按解析半径分层',H.table(['解析半径','独立查询','NN精确MAE ± 标准差','Aphros PLIC精确MAE','三个随机种子均达精度'],circle_bins)))
    shapes_plot=plot_shapes(bundle,dest);common_plot=plot_common(bundle,dest)
    parts.append(H.figure_html(shapes_plot,'左至右为圆、椭圆、弱正弦；上排精确运输体积误差，下排NN／Aphros PLIC误差比，对数纵轴。该图每N使用其解析支持，覆盖率不同，不能单独证明收敛；两套原生PLIC及三个神经种子全部保留。'))
    parts.append(H.figure_html(common_plot,f'仅使用所有N共同解析且固定的{common[0]["samples"]}个原始查询。没有随N换样本，也没有用h²缩放制造下降；不拟合或宣称已获得收敛阶。'))
    parts.append(common_support_table(common,levels))
    cov=[]
    for n in levels:
        record=aggregate_rows(bundle,'physical',n);cov.append([str(n),str(record['samples']),str(record['nn'][0]['original_samples']),H.fmt(100*record['resolution_coverage'])+'%',str(aggregate_rows(bundle,'physical',n,subset='resolved_allN_common')['samples']),str(record['failures'])])
    parts.append(H.detail('逐N解析覆盖与固定共同支持',H.table(['N','主解析查询','全部原始查询','主解析覆盖','共同解析查询','预测失败合计'],cov)))
    # Formal measurements stay separate from old-fixture discovery.
    parts.append('<h2>新真场上的正式部署时间</h2><p>每场为新种子生成的真实完整VOF场，包含全部实际界面查询，不复制验证样本。CPU保留两套PLIC各1／4／8线程，逐案取实际最快者；GPU三个随机种子全保留。驻留与含H2D／D2H的往返分开测量，均同步到完成输出；不含压力／推进／可视化。</p>')
    time_plot=plot_formal_ratios(formal,dest)
    parts.append(H.figure_html(time_plot,'上排单圆滴，下排多圆滴；左为实际最快原生CPU／NN热调用时间，超过1表示网络更快；右为NN／同一原生参考的精确误差，小于1表示网络更准。两种部署边界与全部三个随机种子保留；所有查询来自新独立完整网格，图中仅比热调用，冷构造成本另表。'))
    timing_rows=[];capture_rows=[]
    for r in formal:
        c=r['case'];fast=r['fastest'];timing_rows.append([H.esc(r['scene']),str(r['N']),str(field(c,'actual_queries','batch','active_interface_cells')),H.boundary(r['mode']),H.pm(r['time_mean'],r['time_sd']),H.fmt(r['native_time']),H.esc(fast['backend'])+' / '+str(fast['workers'])+'线程',H.pm(r['error_mean'],r['error_sd']),H.fmt(r['native_error']),('全部达标' if r['all_seeds_pareto'] else '未达，保留')])
        for point,capture,cold,uses,cold_uses in zip(r['points'],r['capture_seconds'],r['cold_constructor_seconds'],r['break_even_reuses'],r['cold_break_even_reuses']):capture_rows.append([H.esc(r['scene']),str(r['N']),H.boundary(r['mode']),str(point['seed']),H.fmt(capture),H.fmt(cold),str(uses) if uses is not None else '未同时达标',str(cold_uses) if cold_uses is not None else '未同时达标'])
    parts.append(H.table(['真实场景','N','实际查询','边界','NN热调用ms ± 标准差','最快CPU ms','实际原生参考','NN精确MAE ± 标准差','原生精确MAE','三种子均精度／时间达标'],timing_rows,True))
    native_detail=[];seed_detail=[]
    for p in bundle['formal']['native_points']:
        native_detail.append([p['scene'],p['N'],p['actual_queries'],p['backend'],p['workers'],formal_error(p),formal_native_time(p),field(p,'p05_ms') if 'p05_ms' in p else p['wall_ms']['p05'],field(p,'p95_ms') if 'p95_ms' in p else p['wall_ms']['p95'],p['metrics']['coverage']])
    for p in bundle['formal']['per_seed_error_time']:
        seed_detail.append([p['scene'],p['N'],p['actual_queries'],p['seed'],p['mode'],formal_error(p),formal_seed_time(p),field(p,'p05_ms') if 'p05_ms' in p else p['wall_ms']['p05'],field(p,'p95_ms') if 'p95_ms' in p else p['wall_ms']['p95'],p['metrics']['coverage']])
    native_headers=['场景','N','实际查询','原生实现','CPU线程','精确MAE','热调用ms中位数','P05 ms','P95 ms','覆盖率'];seed_headers=['场景','N','实际查询','随机种子','边界','精确MAE','热调用ms中位数','P05 ms','P95 ms','覆盖率']
    def formatted_detail(data):
        return [[H.esc(H.boundary(v) if v in ('resident','cpu_gpu_cpu') else v) if isinstance(v,str) else H.fmt(v) for v in row] for row in data]
    parts.append(H.detail('全部60个原生正式点，不仅最快参考',H.table(native_headers,formatted_detail(native_detail),True)))
    parts.append(H.detail('全部60个神经正式点，保留每个随机种子',H.table(seed_headers,formatted_detail(seed_detail),True)))
    parts.append('<p>神经组时间均值与标准差来自三个随机种子的各31次调用中位数；表内P05／P95来自各自真实调用分布。它们与随机种子间标准差分开，不当作同一个不确定度。</p>')
    for filename,headers,data in [('formal_native_points.csv',native_headers,native_detail),('formal_neural_seed_points.csv',seed_headers,seed_detail)]:
        with (dest/filename).open('w',encoding='utf-8-sig',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(headers);writer.writerows(data)
        parts.append(H.attachment(dest/filename))
    parts.append('<h3>固定形状的热推理与实测捕获成本</h3><p>CUDA Graph的热时延仅适用于已捕获的实际查询数。不同批量必须重新捕获；相同形状的新输入仍需按对应边界更新输入。capture_seconds按实现实际记录12次普通预热与图捕获，不计入热调用；它未包含此前的模型加载与缓冲分配。另列cold_load_and_capture_seconds，覆盖这次完整适配器构造，包括加载、缓冲和捕获；它也不是从原始磁盘数据开始的完整PDE首次步。仅当同一案正式精度与时间全部达标时，分别列ceil(对应初始化秒数×1000／(CPU毫秒−NN热调用毫秒))，只摊销所标明的实测初始化段。</p>'+H.detail('逐场、逐随机种子初始化成本与有限复用盈亏点',H.table(['场景','N','边界','随机种子','12次预热＋捕获秒','完整适配器构造秒','预热捕获段复用次数','完整构造段复用次数'],capture_rows,True)))
    # Causal training effects are within-stage pairs only.
    pairs=stage2_pairs(bundle['s2']);pair_rows=[]
    for pair in pairs:
        c=pair['config'];pair_rows.append([f"{c['stencil_size']}×{c['stencil_size']}、宽{c['width']}、{c['num_bins']}bin",H.pm(pair['single']['error_mean'],pair['single']['error_seed_sd']),H.pm(pair['multi']['error_mean'],pair['multi']['error_seed_sd']),f"{pair['reduction']:+.1f}%",' / '.join(f'{v:+.1f}%' for v in pair['paired'])])
    parts.append('<h2>数据、目标与部署优化分别解释</h2><p>下表只比较第二阶段同源权重、架构、数据、目标与新增预算相同的单C／16C；两者均为运输体积L1。阶段1到阶段2的总变化不能全归因于多C；它还更换了训练目标。</p>'+H.table(['源架构','单C验证MAE ± 标准差','16C验证MAE ± 标准差','组均MAE下降','三个配对随机种子下降'],pair_rows))
    curv_rows=[]
    names={'raw':'原网络','d4':'D4空间平均','phase_odd':'相补奇对称平均','d4_phase_odd':'D4＋相补平均'}
    for g in curv['groups']:curv_rows.append([names[g['method']],H.pm(g['q_mae_mean'],g['q_mae_std']),f"{100*(1-g['q_mae_mean']/raw['q_mae_mean']):+.2f}%",f"{100*(g['q_mae_mean']/native-1):+.1f}%",str(g['views']['logical_source_predictions_per_sample']),('通过' if g['all_guard_pass'] else '未通过')])
    parts.append('<h3>曲率有限负结果：停止扩展，未做盲测</h3><p>这些平均只对同一纯神经网络的空间／相补输入作对称变换，不加入几何曲率或回退。全部三随机种子与原分布保留条件均检查；误差仍高于真实particles，因此结束该有限诊断，没有按视图数估算或宣称等精度速度。</p>'+H.table(['纯神经诊断','验证q MAE ± 标准差','相对原网络改善','比particles误差高','每样本原网络预测数','原分布保留条件'],curv_rows))
    dev_rows=[]
    for g in bundle['dev']['groups']:
        if g['group']!=GROUP or g['mode']!='cpu_gpu_cpu':continue
        refs=[n for n in bundle['dev']['native'] if n['workload_id']==g['workload_id']];fast=min(refs,key=lambda n:n['wall_ms']['median'])
        dev_rows.append([H.esc(g['workload_id'].replace('real_grid_N4096_','开发N4096 · ')),('CUDA Graph' if g['runtime']=='cuda_graph' else '普通TorchScript'),H.pm(g['mean_ms'],g['seed_sd_ms']),H.fmt(fast['wall_ms']['median']),H.pm(g['mean_volume_mae'],g['seed_sd_volume_mae']),H.fmt(fast['metrics']['error'])])
    parts.append('<h3>部署发现来自已看开发场，而非最终胜出证据</h3><p>此表是旧开发场与共享负载下的发现过程，只改变同权重推理调度。它解释为何锁定5×5候选，不代替上面的新独立真场正式结果。</p>'+H.table(['已看开发工作负载','同权重推理方式','往返ms ± 标准差','最快CPU ms','开发MAE ± 标准差','原生开发MAE'],dev_rows))
    allrows=[]
    for r in rows:
        if r['subset']!='resolved_primary' or r['radius_bin']!='all' or r['shape']=='all':continue
        method='NN' if not nid(r) else ('Aphros PLIC' if 'aphros' in r['model_id'] else 'SimLiquid PLIC')
        allrows.append([H.esc(r['dataset']),r['N'],H.esc(r['shape']),H.esc(r['subset']),method,r['seed'] if not nid(r) else '确定性',r['samples'],r['failures'],H.fmt(100*(f(r,'coverage') or 0))+'%',H.fmt(f(r,'exact_cell_transport_mae')),H.fmt(f(r,'ratio_mae')),H.fmt(f(r,'flux_volume_mae'))])
    parts.append(H.detail('全部最终独立精度：形状、支持集与每个随机种子',H.table(['集','N','形状','预定义支持','方法','随机种子','查询','失败','有效覆盖','精确运输MAE','ratio MAE','选择期口径MAE'],allrows,True)))
    parts.append('<h2>内嵌证据与代码</h2><p>正文不展示机器标识或校验串；下列PDF、完整分层CSV与代码均内嵌，不依赖外部站点。全原始查询、未解析描述性结果与共同有限支持条件误差在CSV中分别保留。</p>')
    parts.append(H.attachment(B/'metrics.csv','独立精度完整分层CSV'))
    for source in reproducible_sources():
        if source.exists():parts.append(H.attachment(source))
    splits=read(P/'training_counts.json');datasets=read(B/'datasets/manifest.json')['files']
    split_rows=[['原解析子集·训练',f"{splits['original_resolved_train']:,}局部几何",'训练与同源对照'],['原解析子集·验证',f"{splits['original_resolved_validation']:,}局部几何",'保留条件与权重选择'],['新增圆·训练',f"{splits['added_train']:,}母圆",f"{splits['added_train']*16:,}监督对＝母圆×16C（相关多标签）"],['新增圆·验证',f"{splits['added_validation']:,}母圆",f"{splits['added_validation']:,}主要单C验证查询"]]
    for d in datasets:
        if d['dataset']=='circle':split_rows.append(['冻结后独立新圆','32,768母圆',f"{d['samples']:,}查询"])
        else:split_rows.append(['冻结后独立物理集','48母形状',f"{d['samples']:,}查询＝48×9偏移×4变体×9网格（跨N相关）"])
    split_rows.append(['冻结后独立完整真场','2个固定物理场景','2场景×5分辨率＝10网格实例；实际界面查询逐场见表'])
    parts.append('<h2>数据划分</h2>'+H.table(['数据部分','独立母几何／场景数','查询／标签数与用途'],split_rows))
    css='body{font-family:Arial,"Microsoft YaHei",sans-serif;color:#19374d;max-width:1450px;margin:36px auto;padding:0 30px;line-height:1.65}h1{font-size:30px}h2{border-top:1px solid #cbddeb;padding-top:18px;margin-top:34px}.lead{background:#f0f6fb;border-left:4px solid #2171b5;padding:18px 22px 18px 44px;font-size:18px}.lead li{padding:8px 0}img{width:100%;height:auto}figure{margin:28px 0}figcaption{font-size:13px;color:#4d6e87}table{border-collapse:collapse;font-size:13px;width:100%;font-variant-numeric:tabular-nums}th,td{padding:9px 11px;border-bottom:1px solid #d6e5ef;white-space:nowrap;text-align:left}th{background:#eaf2f8;position:sticky;top:0}.tablebox{overflow:auto;max-height:650px;border:1px solid #cbddea;margin:12px 0 22px}tbody tr:nth-child(even){background:#f8fbfe}details{padding:12px 15px;border:1px solid #cadae8;margin:18px 0}summary{cursor:pointer}.download{display:inline-block;border:1px solid #cbddea;padding:5px 9px;margin:5px;color:#08519c;text-decoration:none}.filter{padding:8px;width:350px;max-width:85%;border:1px solid #acc5d9}'
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>纯神经局部算子：最终冻结评估</title><style>'+css+'</style>'+''.join(parts)+'<script>function filterTable(el){let q=el.value.toLowerCase();el.nextElementSibling.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}</script></html>'
    target=P/'final_report.html' if not (P/'final_report.html').exists() else dest/'final_report.html';temp=dest/'final.tmp.html';temp.write_text(H.naturalize_html(doc),encoding='utf-8');temp.replace(target)
    H.write(dest/'numerical_checks.json',dict(passed=True,circle= dict(mean=circle['mean'],sd=circle['sd'],native=nr,reduction_percent=gain),physical8192={k:dict(mean=v['mean'],native=f(aphros(v),'exact_cell_transport_mae'),ratio=v['mean']/f(aphros(v),'exact_cell_transport_mae')) for k,v in high.items()},curvature=dict(symmetry_gain_percent=symgain,remaining_above_native_percent=gap),formal_roundtrip_all_seed_joint_pass=formal_wins,formal_cases=10))
    return target,[shapes_plot,common_plot,time_plot]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--check-inputs',action='store_true');ap.add_argument('--accuracy-preview',action='store_true',help='Explicit preview:independent accuracy complete,formal timing still waiting;does not weaken default finalgate');ap.add_argument('--reuse-figures-from',default=None,help='Preview only:reuse PNG/PDF bytes when previous scientific inputs match exactly');args=ap.parse_args()
    folder=P/'figures';folder.mkdir(exist_ok=True);guard=(folder/'final_report.lock').open('a+')
    try:fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print(json.dumps(dict(state='already_running_no_duplicate')));return
    if args.prepare_only:
        py_compile.compile(str(Path(__file__)),cfile=str(folder/'final_report_compile.pyc'),doraise=True);H.write(folder/'FINAL_REPORT_INPUTS.json',contract())
        status('prepared_waiting_independent_formal_data',compiled=True,report_generated=False,final_conclusions_generated=False,input_contract='figures/FINAL_REPORT_INPUTS.json');print(json.dumps(dict(state='prepared',report_generated=False,gpu_used=False)));return
    bundle,missing=load_inputs(require_formal=not args.accuracy_preview)
    if missing:status('waiting_for_independent_accuracy_and_formal_time',missing=missing,report_generated=False);print(json.dumps(dict(state='waiting',missing=missing)));return
    if args.check_inputs:status('inputs_ready_no_report',report_generated=False);print(json.dumps(dict(state='inputs_ready',report_generated=False)));return
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');dest=folder/(('accuracy_preview_run_' if args.accuracy_preview else 'final_run_')+stamp);dest.mkdir(exist_ok=False)
    reuse=None
    if args.reuse_figures_from:
        assert args.accuracy_preview
        value=Path(args.reuse_figures_from);reuse=(value if value.is_absolute() else P/value).resolve();reuse.relative_to(P)
    target,figures=render_accuracy_preview(bundle,dest,reuse) if args.accuracy_preview else render(bundle,dest);H.write(dest/'inputs_used.json',dict(files=[dict(path=str(p.relative_to(P)),sha256=sha(p)) for p in bundle['required']]))
    status('accuracy_preview_generated_unreviewed' if args.accuracy_preview else 'generated_unreviewed',preview=args.accuracy_preview,report_generated=True,final_report_generated=not args.accuracy_preview,formal_time_used=not args.accuracy_preview,report_path=str(target.relative_to(P)),figures=[str(p.relative_to(P)) for pair in figures for p in pair],visual_review_required=True,curvature_blind_tested=False,old_interim_report_preserved=True)
    print(json.dumps(dict(state='generated_unreviewed',report=str(target),visual_review_required=True)))
if __name__=='__main__':
    try:main()
    except BaseException as error:status('failed',error=repr(error),traceback=traceback.format_exc());raise
