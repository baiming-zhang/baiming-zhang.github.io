"""Build one self-contained scientific HTML only after final blind accuracy.

CPU-only renderer. --prepare-only compiles this source and records readiness;
it does not read outcomes, create figures/conclusions, or initialize CUDA.
Existing reports are retained; later renders use figures/run_*/report.html.
"""
from pathlib import Path
import argparse,base64,csv,datetime,fcntl,hashlib,html,io,json,math,py_compile,re,statistics,sys,traceback

P=Path(__file__).resolve().parent;W=P.parent;R=P;FIG=P/'figures';STATUS=P/'report_status.json'
B=P/'blind_test';DEPLOYMENT_LOCK=P/'deployment_selection_locked.json'
GROUP_META={};GROUP_LABELS={};GROUP_ENGLISH={};GROUP_STYLES={}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8');temp.replace(path)
def state(name,**extra):write(STATUS,dict(state=name,updated_utc=now(),**extra))
def esc(x):return html.escape(str(x))
def number(x):
    try:value=float(x)
    except (ValueError,TypeError):return None
    return value if math.isfinite(value) else None
def fmt(x,digits=4):
    value=number(x)
    if value is None:return '—'
    if value==0:return '0'
    if abs(value)<1e-5 or abs(value)>=1e5:return f'{value:.{digits-1}e}'
    return f'{value:.{digits}g}'
def mean_sd(values):
    if any(v is None for v in values) or not values:return None,None
    return statistics.mean(values),statistics.stdev(values) if len(values)>1 else 0.
def pm(mean,sd):return fmt(mean)+' ± '+fmt(sd,3) if mean is not None else '未完整覆盖'
def human_group(group):
    if group in GROUP_LABELS:return GROUP_LABELS[group]
    if 'compact' in group:
        found=re.search(r's(\d+)_w(\d+)_b(\d+)',group)
        return f'紧凑MLP {found[1]}×{found[1]} · w{found[2]} · {found[3]} bins' if found else group
    return '原MLP 7×7 · 加圆数据' if 'resolved_data_augmented' in group else '原MLP 7×7 · 原数据续训'
def english_group(group):
    if group in GROUP_ENGLISH:return GROUP_ENGLISH[group]
    if 'compact' in group:
        m=re.search(r's(\d+)_w(\d+)_b(\d+)',group);return f'NN {m[1]}x{m[1]}, w{m[2]}, {m[3]} bins' if m else group
    return 'NN 7x7 + circle data' if 'resolved_data_augmented' in group else 'NN 7x7, original data'

def short_legend(group):
    if R==P:return english_group(group)
    m=GROUP_META[group];k=m['config']['stencil_size'];condition=str(m.get('condition','')).lower()
    if m['operator']=='flux':suffix='16C' if 'sixteen' in condition or 'multi' in condition else '1C'
    else:suffix='qL1' if 'l1' in condition else 'qMSE'
    return f'NN {k}×{k}, {suffix}'

def architecture_note():
    descriptions=[]
    for op in ('flux','curvature'):
        configs=sorted({(m['config']['stencil_size'],m['config']['width'],m['config'].get('num_bins',64)) for m in GROUP_META.values() if m['operator']==op})
        values=[f'{k}×{k}宽{width}'+(f'、{bins}个bin' if op=='flux' else '') for k,width,bins in configs]
        descriptions.append(('通量：' if op=='flux' else '曲率：')+'；'.join(values))
    return '。'.join(descriptions)+'。1C／16C表示训练监督的Courant点数，推理接口仍为单个Courant。'

def naturalize_html(doc):
    # Change visible prose only: never mutate embedded base64/code/PDF attributes.
    words=[('onnx_cuda','ONNX GPU'),('torchscript_cuda','TorchScript GPU'),('native_cpu_f64','原生CPU双精度'),('非正式共享GPU scout','共享负载下的非正式GPU试测'),('共享GPU的非正式scout','共享负载下的非正式GPU试测'),('stage2','第二阶段'),('Stage2','第二阶段'),('stage1','第一阶段'),('Stage1','第一阶段'),('阶段1','第一阶段'),('阶段2','第二阶段'),('本stage','本阶段'),('raw VOF','原始VOF体积分数'),('rawVOF','原始VOF体积分数'),('flux另含Courant','通量另含Courant数'),('scout','共享负载试测'),('三seed','三个随机种子'),('3seed','三个随机种子'),('seed','随机种子'),('runtime','推理后端'),('Batch','批量')]
    pieces=re.split('(<[^>]*>)',doc);skip=False
    for i,piece in enumerate(pieces):
        if piece.startswith('<'):
            if re.match(r'<(?:style|script)\b',piece,re.I):skip=True
            elif re.match(r'</(?:style|script)\s*>',piece,re.I):skip=False
        elif not skip:
            if re.fullmatch(r'[A-Za-z0-9_.-]+\.(?:py|csv|pdf|json|txt)',piece):continue
            for source,target in words:piece=piece.replace(source,target)
            pieces[i]=piece
    return ''.join(pieces)

def initialize_group_labels(models):
    GROUP_META.clear();GROUP_LABELS.clear();GROUP_ENGLISH.clear();GROUP_STYLES.clear()
    palette=['#b8d5eb','#6baed6','#4292c6','#2171b5','#1764a3','#0b4a88','#08366c','#527fa8']
    groups=sorted({m['group'] for m in models.values()})
    for index,group in enumerate(groups):
        m=next(m for m in models.values() if m['group']==group);c=m['config'];op=m['operator'];k=c['stencil_size'];width=c['width'];bins=c.get('num_bins',64)
        GROUP_META[group]=m
        descriptor=' '.join(str(m.get(key,'')) for key in ('condition','objective','loss','group')).lower()
        if R!=P:
            if op=='curvature':
                objective='qL1' if any(v in descriptor for v in ('ql1','q_l1','q-l1')) else ('qMSE' if any(v in descriptor for v in ('qmse','q_mse','q-mse')) else f'方案{index+1}')
                label=f'纯MLP {k}×{k} · w{width} · {objective}';en=f'NN {k}x{k}, w{width}, {objective}'
            else:
                multi=any(v in descriptor for v in ('multi','16c','cdf','courant16','sixteen','16_courant'))
                single=any(v in descriptor for v in ('single','1c','one_c','courant1'))
                objective='16C' if multi else ('单C' if single else f'方案{index+1}')
                label=f'纯MLP {k}×{k} · w{width} · {bins} bins · {objective}';en=f'NN {k}x{k}, w{width}, b{bins}, '+('16C' if multi else ('1C' if single else f'candidate{index+1}'))
        else:
            label=human_group(group);en=english_group(group)
        GROUP_LABELS[group]=m.get('display_name',label);GROUP_ENGLISH[group]=en
        GROUP_STYLES[group]=group_style(group) if R==P else (palette[index%len(palette)],{3:'^',5:'s',7:'o'}.get(k,'D'))
def boundary(mode):return 'GPU驻留' if mode=='resident' else 'CPU→GPU→CPU'
def native_label(backend,op):return ('Aphros' if backend=='aphros' else 'SimLiquid')+(' 生产PLIC通量' if op=='flux' else ' 原生particles曲率')
def table(headers,rows,search=False):
    top='<input class="filter" placeholder="筛选表格文字" oninput="filterTable(this)">' if search else ''
    return top+'<div class="tablebox"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def detail(title,body):return '<details><summary>'+esc(title)+'</summary>'+body+'</details>'
def attachment(path,label=None):
    path=Path(path);mime='application/pdf' if path.suffix=='.pdf' else ('text/csv' if path.suffix=='.csv' else 'text/plain')
    return f'<a class="download" download="{esc(path.name)}" href="data:{mime};base64,{base64.b64encode(path.read_bytes()).decode()}">{esc(label or path.name)}</a>'

def csv_rows(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:return list(csv.DictReader(stream))
def full_scope(row):return str(row.get('shape','all')).lower() in ('all','aggregate','') and str(row.get('radius_bin','all')).lower() in ('all','aggregate','')
def row_n(row):return int(float(row['N']))
def is_native(row,models):return row['model_id'] not in models
def row_group(row,models):return models[row['model_id']]['group'] if row['model_id'] in models else row['model_id']
def primary(row):
    if row['operator']=='flux':return number(row.get('flux_volume_mae'))
    return number(row.get('kappa_mae' if row['dataset']=='physical' else 'q_mae'))
def row_backend(row):
    backend=row.get('backend','').lower()
    if backend in ('aphros','simliquid'):return backend
    ident=row['model_id'].lower()
    if 'aphros' in ident:return 'aphros'
    if 'simliquid' in ident:return 'simliquid'
    raise ValueError('Unrecognized native identity '+row['model_id'])

def load_inputs():
    required=[B/'completion.json',B/'metrics.csv',B/'summary.json',R/'validation_pareto/scout/summary.json',DEPLOYMENT_LOCK,R/'inference/model_manifest.json']
    missing=[str(p.relative_to(P)) for p in required if not p.exists()]
    if missing:return None,missing
    if read(required[0]).get('state')!='complete':return None,['blind_test/completion.json: notcomplete']
    scout=read(required[3]);assert scout['state']=='complete' and scout['validation_only'] and not scout['independent_final_test']
    assert not scout['speed_evidence_eligible'],'Scout must never be labeled formal'
    lock=read(required[4]);manifest=read(required[5]);models={m['model_id']:m for m in manifest['models']}
    manifest_hash=hashlib.sha256(required[5].read_bytes()).hexdigest()
    assert scout['identity']['model_manifest_sha256']==manifest_hash,'Scout belongs to different exported weights/stage'
    assert read(required[2])['identity']['export_manifest_sha256']==manifest_hash,'Blind belongs to different exported weights/stage'
    assert len(models)==24 and set(models)==set(lock['model_ids'])
    assert all(len([m for m in models.values() if m['group']==g])==3 for g in {m['group'] for m in models.values()})
    rows=csv_rows(required[1]);assert rows
    needed={'model_id','operator','dataset','N','samples','valid_samples','failures','coverage','q_mae','kappa_mae','ratio_mae','flux_volume_mae'}
    assert needed.issubset(rows[0]),('Blind CSV schema missing',sorted(needed-set(rows[0])))
    rows=[r for r in rows if full_scope(r)];assert rows
    for dataset in ('circle','physical'):
        present={r['model_id'] for r in rows if r['dataset']==dataset}
        assert set(models).issubset(present),(dataset,'Missingseedmodels')
    assert all(row_n(r)==-1 for r in rows if r['dataset']=='circle'),'Normalized circles cannot be relabeled physical refinement'
    physical_levels=sorted({row_n(r) for r in rows if r['dataset']=='physical'});assert len(physical_levels)>=3
    timing=scout;timing_kind='scout'
    formal=R/'validation_pareto/formal/summary.json'
    if formal.exists():
        f=read(formal)
        if f.get('state')=='complete' and f.get('speed_evidence_eligible'):
            assert f['identity']['model_manifest_sha256']==scout['identity']['model_manifest_sha256']
            timing=f;timing_kind='formal'
    initialize_group_labels(models)
    return dict(scout=scout,timing=timing,timing_kind=timing_kind,lock=lock,models=models,rows=rows,blind=read(required[2]),required=required,levels=physical_levels),[]

def blind_aggregates(bundle):
    models=bundle['models'];rows=bundle['rows'];out=[]
    for dataset in ('circle','physical'):
        levels=sorted({row_n(r) for r in rows if r['dataset']==dataset})
        for op in ('flux','curvature'):
            for n in levels:
                scope=[r for r in rows if r['dataset']==dataset and r['operator']==op and row_n(r)==n]
                for group in sorted({row_group(r,models) for r in scope}):
                    selected=[r for r in scope if row_group(r,models)==group];nn=selected[0]['model_id'] in models
                    if nn:assert sorted(int(float(r['seed'])) for r in selected)==[4101,4102,4103],(dataset,n,group)
                    else:assert len(selected)==1,(dataset,n,group,'Duplicate aggregate native rows')
                    values=[primary(r) for r in selected];mu,sd=mean_sd(values);coverage=min(number(r['coverage']) or 0. for r in selected)
                    out.append(dict(dataset=dataset,operator=op,N=n,group=group,is_nn=nn,backend=None if nn else row_backend(selected[0]),
                        error=mu,sd=sd,coverage=coverage,failures=sum(int(float(r['failures'])) for r in selected),rows=selected,
                        label=human_group(group) if nn else native_label(row_backend(selected[0]),op)))
    return out

def locked_comparisons(bundle,aggregates):
    result=[];lock=bundle['lock']
    for op,choices in lock['chosen_groups'].items():
        backend=choices['native_reference']['backend']
        for role in ('best_accuracy','fastest_reaching_native_error','most_accurate_within_native_time'):
            group=choices.get(role)
            if group is None:continue
            for dataset in ('circle','physical'):
                levels=sorted({r['N'] for r in aggregates if r['dataset']==dataset})
                for n in levels:
                    nn=next(r for r in aggregates if (r['dataset'],r['operator'],r['N'],r['group'])==(dataset,op,n,group))
                    native=next(r for r in aggregates if not r['is_nn'] and (r['dataset'],r['operator'],r['N'],r['backend'])==(dataset,op,n,backend))
                    errors=[primary(r) for r in nn['rows']]
                    comparable=nn['coverage']==1 and native['coverage']==1 and native['error'] is not None and all(e is not None for e in errors)
                    allseeds=comparable and all(e<=native['error'] for e in errors)
                    result.append(dict(operator=op,role=role,group=group,dataset=dataset,N=n,backend=backend,comparable=comparable,all_seeds_reach_native=allseeds,nn=nn,native=native))
    return result

def headline(bundle,comparisons):
    formal=bundle['timing_kind']=='formal';lines=[]
    for op in ('flux','curvature'):
        candidates=[r for r in comparisons if r['operator']==op and r['role']=='best_accuracy']
        passed=sum(r['all_seeds_reach_native'] for r in candidates);total=sum(r['comparable'] for r in candidates);missing=len(candidates)-total
        label='通量' if op=='flux' else '曲率'
        lines.append(f'验证期锁定的最低误差{label}组，在新圆盲测及{len(bundle["levels"])}个物理网格层级中，有{passed}/{total}个可比较单元的三种子均达到对应CPU原生误差；另{missing}个单元存在覆盖问题而无法完整比较。')
    # Assess only candidates frozen before the blind data. This never changes
    # weights/runtimes or chooses a new deployment winner on final-test scores.
    precision_pass=any(c['all_seeds_reach_native'] for c in comparisons)
    joint_pass=False
    for c in comparisons:
        if not c['all_seeds_reach_native']:continue
        for g in bundle['timing']['groups']:
            if (g['group'],g['operator'])!=(c['group'],c['operator']) or not g['full_coverage']:continue
            refs=[r for r in bundle['timing']['native_points'] if (r['backend'],r['operator'],r['batch'])==(c['backend'],c['operator'],g['batch']) and r['coverage']==1]
            if not refs:continue
            ref=min(refs,key=lambda r:r['time_ms'])
            strict=g['time_max_seed_ms']<ref['time_ms'] or max(primary(r) for r in c['nn']['rows'])<c['native']['error']
            if g['time_max_seed_ms']<=ref['time_ms'] and strict:joint_pass=True
    if not any(c['comparable'] for c in comparisons):
        opening='完整误差比较存在覆盖问题，当前不能确认胜出；失败与共同支持条件误差分别保留。'
    elif not precision_pass:
        opening='本轮没有证明胜过CPU原生方法：预锁候选尚未在盲测比较单元中实现三种子全部达到对应原生精度。'
    elif not joint_pass:
        opening='本轮出现局部盲测精度达标，但尚未观察到对应候选的精度与部署时间同时不劣于原生参考，未证明Pareto胜出。'
    elif not formal:
        opening='限定场景中观察到盲测精度达标和更快的scout部署时间，但共享GPU计时不能作为正式Pareto胜出的证据。'
    else:
        opening='相对明确列出的CPU原生参考，已观察到预锁候选在限定场景的正式精度—时间Pareto胜出；成立范围限定于对应数据、Batch与部署边界。'
    if not formal:opening+=' 本报告的所有部署时间均为非正式共享GPU scout。'
    return opening+' '+' '.join(lines)

def plot_style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':.75,'axes.edgecolor':'#253e55','axes.labelcolor':'#253e55','xtick.color':'#253e55','ytick.color':'#253e55','pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':240})
    return plt
def group_style(group):
    if group in GROUP_STYLES:return GROUP_STYLES[group]
    if 'compact' not in group:return ('#08519c','o') if 'resolved_data_augmented' in group else ('#9ecae1','o')
    m=re.search(r's(\d+)_w(\d+)',group);return {('3','32'):('#c6dbef','^'),('3','64'):('#6baed6','s'),('5','32'):('#4292c6','D'),('5','64'):('#2171b5','P')}[(m[1],m[2])]
def box(ax):
    for spine in ax.spines.values():spine.set_visible(True)
    ax.tick_params(direction='in',top=True,right=True);ax.grid(True,which='major',linewidth=.35,color='#d8e5f0');ax.set_axisbelow(True)
def plot_pareto(summary,batch,dest):
    plt=plot_style();from matplotlib.lines import Line2D
    fig,axes=plt.subplots(2,2,figsize=(12.8,8.4));fig.subplots_adjust(left=.09,right=.72,bottom=.09,top=.985,wspace=.37,hspace=.31)
    legends=[{},{}]
    for row,op in enumerate(('flux','curvature')):
        for col,mode in enumerate(('resident','cpu_gpu_cpu')):
            ax=axes[row,col];ys=[]
            for g in summary['groups']:
                if (g['operator'],g['batch'],g['mode'])!=(op,batch,mode):continue
                if g['error_mean'] is None or not g['full_coverage']:continue
                color,marker=group_style(g['group']);mean=g['error_mean'];sd=g['error_seed_sd'];ys.append(mean)
                ax.errorbar(g['time_mean_ms'],mean,xerr=g['time_seed_sd_ms'],yerr=min(sd,mean*.99) if mean>0 else sd,fmt=marker,color=color,ecolor=color,capsize=2,markersize=5.2,linewidth=.9,zorder=3)
                for seed in summary['per_seed_error_time']:
                    if (seed['group'],seed['batch'],seed['mode'])==(g['group'],batch,mode) and seed['error'] is not None:
                        ax.plot(seed['median_ms'],seed['error'],marker=marker,markersize=2.5,color=color,alpha=.55,linestyle='none',zorder=2)
                label=short_legend(g['group']);legends[row][label]=Line2D([],[],color=color,marker=marker,linestyle='none',markersize=5,label=label)
            for native in summary['native_points']:
                if (native['operator'],native['batch'])!=(op,batch) or native['error'] is None:continue
                color='#253e55' if native['backend']=='aphros' else '#64869e';marker={1:'v',4:'<',8:'>'}[native['workers']]
                ax.scatter(native['time_ms'],native['error'],marker=marker,s=36,facecolors='none',edgecolors=color,linewidths=1.1,zorder=4);ys.append(native['error'])
                method='PLIC' if op=='flux' else 'particles';brand='Aphros' if native['backend']=='aphros' else 'SimLiquid'
                threads='thread' if native['workers']==1 else 'threads'
                label=f"{brand} {method}, {native['workers']} CPU {threads}";legends[row][label]=Line2D([],[],color=color,marker=marker,markerfacecolor='none',linestyle='none',markersize=6,label=label)
            ax.set_xscale('log')
            from matplotlib.ticker import LogLocator,FuncFormatter,NullFormatter
            ax.xaxis.set_major_locator(LogLocator(base=10,subs=(1,2,5)));ax.xaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'));ax.xaxis.set_minor_formatter(NullFormatter())
            if ys and min(ys)>0:ax.set_yscale('log')
            ax.set_xlabel('GPU-resident call time (ms)' if mode=='resident' else 'CPU-GPU-CPU call time (ms)')
            ax.set_ylabel('Transported-volume MAE' if op=='flux' else 'Curvature q MAE');box(ax)
    for row,legend in enumerate(legends):
        nn=[(label,handle) for label,handle in legend.items() if label.startswith('NN ')]
        native=sorted([(label,handle) for label,handle in legend.items() if not label.startswith('NN ')],key=lambda x:x[0])
        fig.legend(handles=[handle for _,handle in nn+native],loc='center left',bbox_to_anchor=(.735,.755 if row==0 else .285),frameon=False,fontsize=8,handletextpad=.6,labelspacing=.7)
    png=dest/f'validation_pareto_B{batch}.png';pdf=png.with_suffix('.pdf');fig.savefig(png,bbox_inches='tight',pad_inches=.08);fig.savefig(pdf,bbox_inches='tight',pad_inches=.08);plt.close(fig)
    return png,pdf

def plot_refinement(aggregates,dest):
    plt=plot_style();from matplotlib.lines import Line2D
    fig,axes=plt.subplots(1,2,figsize=(12.8,4.7));fig.subplots_adjust(left=.085,right=.73,bottom=.17,top=.98,wspace=.34);legend={}
    for ax,op in zip(axes,('flux','curvature')):
        rows=[r for r in aggregates if r['dataset']=='physical' and r['operator']==op];ys=[]
        for group in sorted({r['group'] for r in rows}):
            series=sorted([r for r in rows if r['group']==group],key=lambda r:r['N']);nn=series[0]['is_nn']
            color,marker=group_style(group) if nn else (('#253e55','v') if series[0]['backend']=='aphros' else ('#64869e','<'))
            label=english_group(group) if nn else f"{series[0]['backend'].capitalize()} native"
            x=[r['N'] for r in series];y=[r['error'] if r['error'] is not None else float('nan') for r in series];ys.extend(v for v in y if math.isfinite(v))
            ax.plot(x,y,marker=marker,markersize=4,linewidth=1,color=color,markerfacecolor=color if nn else 'none')
            if nn:
                for r in series:
                    if r['error'] is not None:ax.errorbar(r['N'],r['error'],yerr=min(r['sd'],r['error']*.99) if r['error']>0 else r['sd'],fmt='none',color=color,linewidth=.7,capsize=1.5)
            legend[label]=Line2D([],[],color=color,marker=marker,linestyle='-',markersize=4,label=label)
        ax.set_xscale('log',base=2)
        levels=sorted({r['N'] for r in rows});ax.set_xticks(levels);ax.set_xticklabels([str(n) for n in levels],rotation=40,ha='right')
        if ys and min(ys)>0:ax.set_yscale('log')
        ax.set_xlabel('Grid resolution N');ax.set_ylabel('Cell-normalized transported-volume MAE' if op=='flux' else 'Physical curvature MAE (1/length)');box(ax)
    fig.legend(handles=list(legend.values()),loc='center left',bbox_to_anchor=(.745,.54),frameon=False,fontsize=8,labelspacing=.85)
    png=dest/'blind_physical_refinement.png';pdf=png.with_suffix('.pdf');fig.savefig(png,bbox_inches='tight',pad_inches=.08);fig.savefig(pdf,bbox_inches='tight',pad_inches=.08);plt.close(fig)
    return png,pdf

def figure_html(files,caption):
    png,pdf=files;encoded=base64.b64encode(png.read_bytes()).decode()
    return '<figure><img alt="'+esc(caption)+'" src="data:image/png;base64,'+encoded+'"><figcaption>'+esc(caption)+' '+attachment(pdf,'下载独立PDF')+'</figcaption></figure>'

def save_html(doc,dest):
    target=P/'report.html' if not (P/'report.html').exists() else dest/'report.html'
    temporary=dest/'rendered.tmp.html';temporary.write_text(naturalize_html(doc),encoding='utf-8');temporary.replace(target)
    return target

def load_validation_failure():
    sources=[R/'autofinish_status.json',R/'optimization_status.json']
    present=[p for p in sources if p.exists()]
    if not any(read(p).get('state')=='validation_goal_not_met' for p in present):
        return None,['Explicit validation_goal_not_met status required;validation-only export refused']
    summary_path=R/'validation_pareto/scout/summary.json';manifest_path=R/'inference/model_manifest.json'
    missing=[str(p.relative_to(P)) for p in (summary_path,manifest_path) if not p.exists()]
    if missing:return None,missing
    summary=read(summary_path);manifest=read(manifest_path);assert summary['state']=='complete' and summary['validation_only'] and not summary['independent_final_test'] and not summary['speed_evidence_eligible']
    assert summary['identity']['model_manifest_sha256']==hashlib.sha256(manifest_path.read_bytes()).hexdigest(),'Negative-result scout belongs to different exported weights/stage'
    models={m['model_id']:m for m in manifest['models']};assert len(models)==24
    assert all(len([m for m in models.values() if m['group']==g])==3 for g in {m['group'] for m in models.values()})
    registration=B/'registration.json'
    if registration.exists():
        identity=read(registration).get('identity',{})
        assert identity.get('export_manifest_sha256')!=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),'Thisstage already registeredblind;do not label it untested'
    initialize_group_labels(models)
    return dict(scout=summary,timing=summary,timing_kind='scout',models=models,required=[summary_path,manifest_path,*present],failure_records=[read(p) for p in present]),[]

def render_validation_failure(bundle,dest):
    summary=bundle['scout'];models=bundle['models'];figures=[plot_pareto(summary,b,dest) for b in (4096,65536)]
    parts=['<header><p>R/h ≥ 4 · stage2验证探索</p><h1>验证探索未达目标</h1><p class="lead">本stage的GPU纯网络尚未达到预注册的CPU原生精度—时间目标，保留全部24模型与各seed负结果。没有因此解锁或生成本stage最终盲测；以下准确度来自参与过选择的验证集，部署时间为非正式共享GPU scout，不能作为正式速度胜出证据。</p></header>',
        '<p>通量比较对象是真正生产PLIC重构与通量；曲率比较对象为原生particles。主线仅输入raw VOF（flux另含Courant），没有冻结平面SVD层或几何回退。阶段1只作开发背景；离线CDF理论可表示性不计为实际部署胜利。</p>']
    # Recalculate the summary numbers from each retained seed before reporting.
    chosen=[g for g in summary['groups'] if g['batch']==65536 and g['mode']=='cpu_gpu_cpu']
    checks=[]
    for g in chosen:
        seeds=[s for s in summary['per_seed_error_time'] if (s['group'],s['batch'],s['mode'])==(g['group'],65536,'cpu_gpu_cpu')]
        assert sorted(s['seed'] for s in seeds)==[4101,4102,4103]
        em,es=mean_sd([s['error'] for s in seeds]);tm,ts=mean_sd([s['median_ms'] for s in seeds])
        for recomputed,reported in ((em,g['error_mean']),(es,g['error_seed_sd']),(tm,g['time_mean_ms']),(ts,g['time_seed_sd_ms'])):
            assert recomputed is not None and math.isclose(recomputed,reported,rel_tol=1e-11,abs_tol=1e-14)
        checks.append(dict(group=g['group'],error_mean=em,error_seed_sd=es,time_mean_ms=tm,time_seed_sd_ms=ts,seeds=[s['seed'] for s in seeds]))
    selected={op:min([g for g in chosen if g['operator']==op and g['full_coverage'] and g['error_mean'] is not None],key=lambda g:g['error_mean']) for op in ('flux','curvature')}
    references={op:min([n for n in summary['native_points'] if (n['operator'],n['batch'])==(op,65536) and n['coverage']==1],key=lambda n:n['time_ms']) for op in ('flux','curvature')}
    f=selected['flux'];nf=references['flux'];q=selected['curvature'];nq=references['curvature']
    volume_reduction=100*(1-f['error_mean']/nf['error']);flux_slowdown=100*(f['time_mean_ms']/nf['time_ms']-1);q_increase=100*(q['error_mean']/nq['error']-1)
    fk=GROUP_META[f['group']]['config']['stencil_size']
    factual=f'通量方面，{fk}×{fk}、16点Courant监督网络的验证运输体积MAE为{fmt(f["error_mean"])}，相对实际原生PLIC的{fmt(nf["error"])}低{volume_reduction:.1f}%；但在Batch65,536的CPU→GPU→CPU边界，调用时间{f["time_mean_ms"]:.3f} ms，原生{nf["workers"]}线程为{nf["time_ms"]:.3f} ms，仍慢{flux_slowdown:.1f}%。曲率方面，qMSE网络的q MAE为{fmt(q["error_mean"])}，比particles的{fmt(nq["error"])}高{q_increase:.1f}%；GPU调用虽快，同精度目标尚未通过。'
    parts.insert(1,'<div class="lead">'+esc(factual)+'</div>')
    training=read(R/'training_protocol.json');pairs=[]
    for source in sorted({r['source_group'] for r in training['groups'] if r['operator']=='flux'}):
        a=next(r for r in training['groups'] if r['source_group']==source and r['condition']=='single_courant')
        b=next(r for r in training['groups'] if r['source_group']==source and r['condition']=='sixteen_courants')
        assert a['sources']==b['sources'] and a['config']==b['config'],'Courant causal comparison must use identical source/config'
        single=next(g for g in chosen if g['group']==a['group']);multi=next(g for g in chosen if g['group']==b['group'])
        changes=[]
        for seed in (4101,4102,4103):
            ae=next(s['error'] for s in summary['per_seed_error_time'] if (s['group'],s['seed'],s['batch'],s['mode'])==(a['group'],seed,65536,'cpu_gpu_cpu'))
            be=next(s['error'] for s in summary['per_seed_error_time'] if (s['group'],s['seed'],s['batch'],s['mode'])==(b['group'],seed,65536,'cpu_gpu_cpu'))
            changes.append(100*(1-be/ae))
        reduction=100*(1-multi['error_mean']/single['error_mean']);config=a['config']
        pairs.append(dict(stencil=config['stencil_size'],width=config['width'],bins=config['num_bins'],single_group=a['group'],multi_group=b['group'],single_mean=single['error_mean'],single_sd=single['error_seed_sd'],multi_mean=multi['error_mean'],multi_sd=multi['error_seed_sd'],reduction_percent=reduction,paired_seed_reduction_percent=changes))
    pair_rows=[[f"{p['stencil']}×{p['stencil']} · w{p['width']} · {p['bins']} bins",pm(p['single_mean'],p['single_sd']),pm(p['multi_mean'],p['multi_sd']),f"{p['reduction_percent']:+.1f}%",' / '.join(f'{v:+.1f}%' for v in p['paired_seed_reduction_percent'])] for p in pairs]
    parts.append('<h2>单C／16C的同源配对比较</h2><p>这里只比较stage2内部相同源权重、架构、数据与新增更新预算下的单C和16C，两臂均采用运输体积L1目标。不是把阶段1到stage2中“改目标＋加C”的总变化都归因于多C。正号表示误差下降，负号表示变差；三seed按4101／4102／4103对齐。</p>'+table(['源架构','单C MAE均值 ± SD','16C MAE均值 ± SD','组均MAE下降','逐seed MAE下降'],pair_rows))
    write(dest/'numerical_checks.json',dict(passed=True,scope='stage2 validation only;no blind or formal inference',group_recomputations=checks,reference_comparison=dict(flux_error_reduction_percent=volume_reduction,flux_roundtrip_slowdown_percent=flux_slowdown,curvature_error_increase_percent=q_increase),courant_pairs=pairs))
    goals=next((record['goals'] for record in bundle['failure_records'] if record.get('goals')),{});goal_rows=[]
    for op,goal in goals.items():
        comparison=goal.get('comparison') or {};native=comparison.get('native') or {}
        goal_rows.append([esc('通量' if op=='flux' else '曲率'),('有局部scout候选；尚非正式' if goal.get('validation_scout_candidate_dominance') else '验证目标未达'),str(comparison.get('batch','—')),esc(boundary(comparison['mode'])) if comparison.get('mode') else '—',fmt(native.get('error')),fmt(native.get('time_ms'))])
    if goal_rows:parts.append('<h2>两个目标分别判断</h2><p>整体门槛未满足时，仍保留单项出现的正面scout信号；它不替代正式时间或最终盲测。</p>'+table(['算子','本轮验证判断','Batch','边界','原生MAE','原生调用ms'],goal_rows))
    group_rows=[]
    for g in summary['groups']:
        group_rows.append([esc(g['operator']),esc(human_group(g['group'])),str(g['batch']),esc(boundary(g['mode'])),pm(g['error_mean'],g['error_seed_sd']),pm(g['time_mean_ms'],g['time_seed_sd_ms']),('完整' if g['full_coverage'] else '失败保留'),esc('/'.join(sorted(set(g['runtimes']))))])
    parts.append('<h2>全部验证候选</h2><p>每组三seed全部保留；表内误差与时间同时列出，不把NN相对旧NN的改善当成胜过原生。</p>'+table(['算子','候选','Batch','边界','MAE均值 ± SD','调用ms均值 ± SD','覆盖','runtime'],group_rows,True))
    for batch,files in zip((4096,65536),figures):parts.append(figure_html(files,f'验证集Batch{batch:,}；上通量／下曲率，左GPU驻留／右CPU→GPU→CPU。小点保留各seed，误差条为三seed样本标准差；CPU为实际1／4／8逻辑CPU线程原生点。全部部署时间为非正式scout。'+architecture_note()))
    iso=[]
    for row in summary['equal_time_equal_error']:
        if row['reference_scope']!='fastest_tested_native':continue
        et=row.get('equal_time_group');ee=row.get('equal_error_group');ref=row['native']
        iso.append([esc(row['operator']),str(row['batch']),esc(boundary(row['mode'])),esc(ref['backend'])+f" {ref['workers']}线程",fmt(ref['time_ms']),fmt(ref['error']),esc(human_group(et['group'])) if et else '不可达',esc(human_group(ee['group'])) if ee else '不可达',('是' if row['all_seeds_fit_time'] else '否'),('是' if row['all_seeds_reach_error'] else '否')])
    parts.append('<h2>实际原生参考下的等时／等误差</h2><p>离散实测点，无插值；不存在满足条件的组时明确写“不可达”。组均值与三seed全部达到条件分开列。</p>'+table(['算子','Batch','边界','CPU参考','参考ms','参考MAE','等时最低误差组','等误差最快组','三seed时间达标','三seed误差达标'],iso,True))
    native_rows=[[esc(native_label(n['backend'],n['operator'])),str(n['workers']),str(n['batch']),fmt(n['error']),fmt(n['time_ms']),fmt(n.get('call_p05'))+'—'+fmt(n.get('call_p95')),fmt(100*n['coverage'])+'%'] for n in summary['native_points']]
    parts.append(detail('全部CPU原生点',table(['实际方法','线程','Batch','MAE','调用ms中位数','P05—P95 ms','覆盖'],native_rows,True)))
    seed_rows=[[esc(s['operator']),esc(human_group(s['group'])),str(s['seed']),str(s['batch']),esc(boundary(s['mode'])),fmt(s['error']),fmt(s['median_ms']),fmt(s['p05_ms'])+'—'+fmt(s['p95_ms']),fmt(100*s['coverage'])+'%',esc(s['runtime'])] for s in summary['per_seed_error_time']]
    parts.append(detail('全部24模型各seed，不删除负结果',table(['算子','候选','seed','Batch','边界','MAE','调用ms中位数','P05—P95 ms','覆盖','runtime'],seed_rows,True)))
    parts.append('<h2>最终盲测与物理收敛</h2><p>本stage未解锁最终盲测，未生成最终circle／独立物理集结果，也没有独立收敛图。验证集不被改名为最终测试，不用h²缩放产生的误差下降替代物理收敛证据。</p>')
    paths=[R/name for name in ('train_stage2.py','compact_models.py','export_stage2.py','validation_stage2.py','autofinish_stage2.py')]+[P/name for name in ('data_generation.py','compact_models.py','train_resolved.py','export_resolved.py','native_validation.py','validation_pareto.py')]+[Path(__file__)]
    parts.append('<h2>可复核代码与验证明细</h2><div>'+''.join(attachment(path,path.name) for path in dict.fromkeys(paths) if path.exists())+'</div>')
    for path in (R/'validation_pareto/scout/per_seed_error_time.csv',R/'validation_pareto/scout/group_error_time.csv',R/'validation_pareto/scout/native_error_time.csv'):
        if path.exists():parts.append(attachment(path))
    counts=read(R/'training_counts.json') if (R/'training_counts.json').exists() else read(P/'training_counts.json')
    parts.append('<h2>数据划分</h2>'+table(['数据部分','数量','用途'],[['原resolved train',str(counts['original_resolved_train']),'训练／对照'],['原resolved validation',str(counts['original_resolved_validation']),'保留条件与选择'],['新增圆train',str(counts['added_train']),'训练'],['新增圆validation',str(counts['added_validation']),'验证探索；不是盲测'],['最终盲测','未生成','不为填充报告而解锁']]))
    css='body{font-family:Arial,"Microsoft YaHei",sans-serif;color:#18384e;max-width:1450px;margin:36px auto;padding:0 28px;line-height:1.65}h1{font-size:30px}.lead{font-size:19px;border-left:4px solid #2171b5;padding:15px;background:#f0f6fb}h2{font-size:22px;margin-top:32px;border-top:1px solid #cfdfeb;padding-top:18px}img{width:100%;height:auto}figure{margin:26px 0}figcaption{font-size:13px;color:#4b6f88}table{border-collapse:collapse;font-size:13px;width:100%;font-variant-numeric:tabular-nums}th,td{padding:8px 10px;border-bottom:1px solid #dae6ef;text-align:left;white-space:nowrap}th{background:#eaf2f9;position:sticky;top:0}tbody tr:nth-child(even){background:#f8fbfe}.tablebox{overflow:auto;max-height:620px;border:1px solid #cbddea;margin:10px 0 20px}details{border:1px solid #cadbe8;padding:12px;margin:15px 0}summary{cursor:pointer}.download{display:inline-block;border:1px solid #cadae8;color:#08519c;padding:5px 9px;margin:5px;text-decoration:none;font-size:13px}.filter{padding:8px;border:1px solid #adc7db;width:350px;max-width:85%}'
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>第二阶段验证探索：目标未达</title><style>'+css+'</style>'+''.join(parts)+'<script>function filterTable(el){let q=el.value.toLowerCase();el.nextElementSibling.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}</script></html>'
    return save_html(doc,dest),figures

def render(bundle,dest):
    aggregates=blind_aggregates(bundle);comparisons=locked_comparisons(bundle,aggregates);scout=bundle['scout'];timing=bundle['timing'];models=bundle['models']
    figures=[plot_pareto(timing,b,dest) for b in (4096,65536)]+[plot_refinement(aggregates,dest)]
    formal=bundle['timing_kind']=='formal';timeword='正式部署测量' if formal else '非正式共享GPU scout'
    training=read(R/'training_protocol.json')
    counts=read(R/'training_counts.json') if (R/'training_counts.json').exists() else read(P/'training_counts.json');parts=[]
    parts.append('<header><p class="eyebrow">R/h ≥ 4 · 单次局部算子 · 全部三种子</p><h1>GPU纯神经算子与CPU原生PLIC／particles</h1><p class="lead">'+esc(headline(bundle,comparisons))+'</p></header>')
    parts.append('<div class="scope"><b>本报告的证据边界</b><p>通量基线是真正的生产PLIC重构与扫掠通量；曲率基线是原生particles，不把PLIC、HF或未测方法混为一谈。主网络仅输入原始VOF，通量另输入Courant；没有冻结平面SVD层、几何先验或回退。</p><p>验证集参与训练、权重／架构选择；新盲测只在部署选择锁定后生成。GPU FP32与CPU float64（相同raw FP32输入提升为double）属于部署实现比较。'+esc(timeword)+'不能自动等同于跨硬件算法优越性。</p></div>')
    # All locked decisions are reported; never choose a new winner on blind data.
    rows=[]
    roles={'best_accuracy':'验证最低误差','fastest_reaching_native_error':'验证达原生误差的最快组','most_accurate_within_native_time':'验证原生时间内最低误差'}
    for c in comparisons:
        rows.append([esc('通量' if c['operator']=='flux' else '曲率'),esc(roles[c['role']]),esc(human_group(c['group'])),esc('新圆' if c['dataset']=='circle' else '独立物理集'),str(c['N']) if c['N']>0 else 'cell-normalized',
            pm(c['nn']['error'],c['nn']['sd']),fmt(c['native']['error']),esc(c['backend']),('不可比较' if not c['comparable'] else ('是' if c['all_seeds_reach_native'] else '否'))])
    parts.append('<h2>锁定部署选择在新盲测上的结果</h2><p>下面逐项审计验证期预先锁定的组，未根据盲测重新选择权重、种子或runtime。“三种子均达标”按各seed与对应实际native比较；物理曲率用κ误差，新圆用q误差。</p>'+table(['算子','预锁角色','神经组','盲测集','N','NN均值 ± SD','原生误差','原生实现','三seed均达标'],rows,True))
    parts.append('<h2>部署精度—时间：'+esc(timeword)+'</h2><p>误差来自16,384个独立验证查询，未拿重复计时样本扩充准确度样本数。每次调用完成后同步；CPU原生计算含normal／plane／particle处理和持久在线程池dispatch。GPU驻留边界不含传输；roundtrip包含H2D、一次算子及D2H。crop、dtype staging与预分配不计入两端，模型runtime列保留ONNX／TorchScript混合情况。B=65,536重复查询只构成工作负载，不是独立准确度样本，也不是完整PDE时间步。</p>')
    for b,files in zip((4096,65536),figures[:2]):parts.append(figure_html(files,f'Batch {b:,}；上排通量，下排曲率；左GPU驻留，右CPU→GPU→CPU。点为组均值，误差条为三seed样本SD，小点为各seed；原生点保留1／4／8线程。图内无标题，所有轴保留完整边框。'+('正式计时。' if formal else '共享GPU非正式scout，不作为正式速度胜出证据。')))
    grows=[]
    for g in timing['groups']:
        grows.append([esc(g['operator']),esc(human_group(g['group'])),str(g['batch']),esc(boundary(g['mode'])),pm(g['error_mean'],g['error_seed_sd']),pm(g['time_mean_ms'],g['time_seed_sd_ms']),esc('/'.join(sorted(set(g['runtimes'])))),('完整' if g['full_coverage'] else '失败保留'),('正式' if g['formal'] else '非正式')])
    parts.append(detail('全部验证组均值、误差与时间（不省略负结果）',table(['算子','组','Batch','边界','MAE均值 ± SD','调用时间ms均值 ± SD','runtime','覆盖','时间性质'],grows,True)))
    nrows=[[esc(native_label(n['backend'],n['operator'])),str(n['workers']),str(n['batch']),fmt(n['error']),fmt(n['time_ms']),fmt(n.get('call_p05'))+'—'+fmt(n.get('call_p95')),fmt(100*n['coverage'])+'%',('正式' if n['formal'] else '非正式')] for n in timing['native_points']]
    parts.append(detail('全部CPU原生点：1／4／8线程',table(['真正原生方法','线程','Batch','MAE','调用ms中位数','P05—P95 ms','覆盖','时间性质'],nrows,True)))
    iso=[]
    for row in timing['equal_time_equal_error']:
        if row['reference_scope']!='fastest_tested_native':continue
        et=row.get('equal_time_group');ee=row.get('equal_error_group');ref=row['native']
        iso.append([esc(row['operator']),str(row['batch']),esc(boundary(row['mode'])),esc(ref['backend'])+f" {ref['workers']}线程",fmt(ref['time_ms']),fmt(ref['error']),esc(human_group(et['group'])) if et else '不可达',esc(human_group(ee['group'])) if ee else '不可达',('是' if row['all_seeds_fit_time'] else '否'),('是' if row['all_seeds_reach_error'] else '否')])
    parts.append('<h3>等时／等误差的离散比较</h3><p>参考是本轮实际测过的最快原生点，不进行曲线插值。“存在组”使用组均值；末两列另给三种子全部满足的更严格判断。这里只是验证选择证据。</p>'+table(['算子','Batch','边界','实际CPU参考','参考ms','参考误差','等时可达的最低误差组','等误差可达的最快组','三seed均满足时间','三seed均满足误差'],iso,True))
    seedrows=[]
    for s in timing['per_seed_error_time']:
        seedrows.append([esc(s['operator']),esc(human_group(s['group'])),str(s['seed']),str(s['batch']),esc(boundary(s['mode'])),fmt(s['error']),fmt(s['median_ms']),fmt(s['p05_ms'])+'—'+fmt(s['p95_ms']),fmt(100*s['coverage'])+'%',esc(s['runtime'])])
    parts.append(detail('全部24模型、各seed和两个部署边界',table(['算子','组','seed','Batch','边界','MAE','调用ms中位数','P05—P95 ms','覆盖','runtime'],seedrows,True)))
    parts.append('<h2>独立物理几何下的细化行为</h2><p>圆的cell-normalized盲测N=-1仅衡量解析尺度精度，不画成网格收敛。以下仅画独立physical集的真实N层级；曲率使用physical κ MAE，不把q=κh随h减小当物理收敛。通量用不带h²的cell-normalized transported-volume误差；h²物理运输面积误差仅在明细中作辅助。覆盖失败保留为缺口，不用available-only结果替代主误差。</p>'+figure_html(figures[2],'独立物理集网格细化；左为cell-normalized运输体积误差，右为physical κ误差。三seed组均值与SD；不拟合或宣称未验证的收敛阶。'))
    blindrows=[]
    for row in bundle['rows']:
        nn=row['model_id'] in models;op=row['operator'];label=human_group(models[row['model_id']]['group']) if nn else native_label(row_backend(row),op)
        blindrows.append([esc(row['dataset']),str(row_n(row)),esc(op),esc(label),esc(row.get('seed','—')) if nn else '确定性',str(row['samples']),str(row['valid_samples']),str(row['failures']),fmt(100*(number(row['coverage']) or 0))+'%',fmt(primary(row)),fmt(row.get('ratio_mae')),fmt(row.get('q_mae')),fmt(row.get('exact_cell_transport_mae')),fmt(row.get('physical_transport_area_mae'))])
    parts.append(detail('完整盲测各seed／各N与覆盖率',table(['集','N','算子','方法','seed','样本','有效','失败','覆盖','主MAE','ratio MAE','q MAE','exact cell transport MAE','physical area MAE（辅助）'],blindrows,True)))
    common=B/'common_finite_metrics.csv'
    if common.exists():parts.append(detail('共同有限支持的条件误差（辅助，不替代失败率）','<p>主表仍使用全范围失败与覆盖率；下面文件只描述所有方法共同有限样本上的条件误差。</p>'+attachment(common)))
    endpoints=B/'endpoint_scopes.json'
    if endpoints.exists():parts.append(detail('C=0、C=1与低Courant分层',attachment(endpoints,'下载端点分层完整JSON')))
    if R==P:
        description='原MLP从既有enhanced7权重继续；compact从头初始化，总历史训练量不同，架构结果作为部署候选比较。原数据续训／加数据对照共享源权重。旧池过滤 |κh|≤0.25；加数据组混入完整圆。曲率损失含相对误差尺度，通量含ratio与运输体积项。'
    else:
        description='主线为stage2：三种通量源架构分别比较单Courant与16点Courant监督，使用运输体积L1目标；两组曲率分别采用qMSE与qL1。每个源内的配对组共享起点与新增更新预算；不同源的历史训练量因来源而异，架构结果作为成品部署比较。阶段1仅为开发对照，不与本stage合并统计或混合归因；CDF理论可表示性仅为离线潜力诊断，不算部署胜利。'
    parts.append('<h2>模型、训练与范围</h2><p>本stage新增'+str(training.get('steps',training.get('updates',12000)))+'次统一更新，batch='+str(training.get('batch','见协议'))+'。'+description+'第一次float32 coverage mean舍入导致的训练尝试保留于开发阶段attempts，未计入有效结果。</p><p>所有有效模型的权重选择只用验证。主线是纯MLP，没有冻结平面SVD层或外部几何算子；通量保留正bin单调累计及容量限制。完整协议与代码嵌在下方下载区。</p>')
    modelrows=[[esc(m['operator']),esc(human_group(m['group'])),str(m['seed']),str(m['stencil_size']),str(m['parameter_count']),str(m['selected_step']),('通过' if m.get('selection',{}).get('eligible') else '未通过或未标注，保留'),esc(bundle['lock']['model_runtimes'][m['model_id']])] for m in models.values()]
    parts.append(detail('全部有效模型与验证选择',table(['算子','组','seed','stencil','参数','选中新增step','原分布保留条件','锁定runtime'],modelrows,True)))
    sources=['data_generation.py','train_resolved.py','compact_models.py','export_resolved.py','native_validation.py','validation_pareto.py','select_deployment.py','blind_evaluation.py','report_resolved.py']
    codepaths=list(dict.fromkeys([R/name if (R/name).exists() else P/name for name in sources]+list(R.glob('*.py'))))
    parts.append('<h2>可复核文件</h2><p>下列代码、明细和PDF均作为数据嵌入本HTML，不依赖外部网页；正文不展示文件校验串或机器标识。</p><div class="attachments">'+''.join(attachment(path,path.name) for path in codepaths if path.exists())+attachment(B/'metrics.csv','盲测完整分层CSV')+'</div>')
    # The split table is intentionally the final scientific section.
    splitrows=[['原resolved池 · train',str(counts['original_resolved_train']),'训练／对照','原增强cubic，|κh|≤0.25'],['原resolved池 · validation',str(counts['original_resolved_validation']),'保留条件与权重选择','不作盲测'],['新增圆 · train',str(counts['added_train']),'训练','R/h=4—512，独立base'],['新增圆 · validation',str(counts['added_validation']),'多次验证／架构和部署选择','不是最终独立测试']]
    for dataset in ('circle','physical'):
        path=B/'datasets'/f'{dataset}.npz'
        if path.exists():
            import numpy as np
            with np.load(path) as z:queries=len(z['query_id']);bases=len(np.unique(z['geometry_id']))
            splitrows.append(['最终盲测 · '+dataset,f'{queries:,}查询 / {bases:,}基础几何','部署锁后新生成','新seed与参数；变体不跨split；全部24模型保留'])
    parts.append('<h2>数据划分</h2>'+table(['数据部分','数量','用途','范围与独立性'],splitrows)+'<p class="foot">circle是新种子与参数的同家族盲测；physical包含独立固定物理几何与多级网格，不把同家族新参数自动称为家族外泛化。计时重复不增加独立准确度样本。报告生成后仍需实际查看PNG／PDF，检查图例、边框、文字和重叠后交付。</p>')
    css='''body{margin:0;background:#fff;color:#152f46;font-family:Arial,"Microsoft YaHei",sans-serif;font-size:15px;line-height:1.65}main{max-width:1450px;margin:42px auto;padding:0 36px}h1{font-size:30px;font-weight:650;line-height:1.25}h2{font-size:22px;margin-top:38px;border-top:1px solid #c8d8e6;padding-top:20px}h3{font-size:18px}.eyebrow{color:#607e98;font-size:13px;letter-spacing:.08em}.lead{font-size:19px;max-width:1200px}.scope{border-left:4px solid #2171b5;background:#f1f6fb;padding:16px 24px;margin:24px 0}.scope p{margin:8px 0}figure{margin:26px 0}img{width:100%;height:auto;display:block}figcaption{font-size:13px;color:#405e76;margin-top:8px}table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}th,td{border-bottom:1px solid #dce6ef;padding:9px 11px;text-align:left;vertical-align:top;white-space:nowrap}th{position:sticky;top:0;background:#e9f1f8;color:#153e61}tbody tr:nth-child(even){background:#f8fbfd}.tablebox{overflow:auto;max-height:620px;border:1px solid #cad9e6;margin:12px 0 20px}details{margin:18px 0;border:1px solid #cfdeea;padding:12px 16px}summary{cursor:pointer;font-weight:600}.filter{width:360px;max-width:85%;padding:8px;border:1px solid #adc5d9;border-radius:2px;margin-top:8px}.download{display:inline-block;color:#08519c;padding:5px 9px;border:1px solid #c9dcea;margin:5px 6px 5px 0;text-decoration:none;font-size:13px}.foot{font-size:13px;color:#506f87}@media(max-width:700px){main{padding:0 16px;margin-top:24px}h1{font-size:24px}.lead{font-size:17px}.tablebox{max-height:450px}}'''
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Resolved neural/native operator comparison</title><style>'+css+'</style><main>'+''.join(parts)+'</main><script>function filterTable(el){let q=el.value.toLowerCase();el.nextElementSibling.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}</script></html>'
    return save_html(doc,dest),figures

def main():
    global R,B,DEPLOYMENT_LOCK
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare-only',action='store_true');parser.add_argument('--stage-root',default=None,help='Stage path relative to the resolved-project root;workspace remains its parent');parser.add_argument('--blind-root',default=None);parser.add_argument('--deployment-lock',default=None);parser.add_argument('--validation-only',action='store_true',help='Explicit negative-results exit;requires validation_goal_not_met status,never unlocks blind');args=parser.parse_args()
    if args.stage_root:
        candidate=Path(args.stage_root);R=(candidate if candidate.is_absolute() else P/candidate).resolve();R.relative_to(P)
    B=R/'blind_test' if (R/'blind_test/completion.json').exists() else P/'blind_test'
    DEPLOYMENT_LOCK=R/'deployment_selection_locked.json' if (R/'deployment_selection_locked.json').exists() else P/'deployment_selection_locked.json'
    if args.blind_root:
        path=Path(args.blind_root);B=(path if path.is_absolute() else P/path).resolve();B.relative_to(P)
    if args.deployment_lock:
        path=Path(args.deployment_lock);DEPLOYMENT_LOCK=(path if path.is_absolute() else P/path).resolve();DEPLOYMENT_LOCK.relative_to(P)
    FIG.mkdir(parents=True,exist_ok=True);guard=(FIG/'report.lock').open('a+')
    try:fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print(json.dumps(dict(state='already_rendering_no_duplicate')),flush=True);return
    if args.prepare_only:
        FIG.mkdir(parents=True,exist_ok=True);py_compile.compile(str(Path(__file__)),cfile=str(FIG/'report_compile.pyc'),doraise=True)
        state('ready_waiting_blind',stage_root=str(R.relative_to(P)) or '.',blind_root=str(B.relative_to(P)),deployment_lock=str(DEPLOYMENT_LOCK.relative_to(P)),workspace_root_inferred_from_project=True,cpu_compile_passed=True,gpu_used=False,report_generated=False,validation_failure_exit='--validation-only requires validation_goal_not_met;no blind generation or test claim',required_inputs=['blind_test/completion.json complete','blind_test/metrics.csv','blind_test/summary.json','validation_pareto/scout/summary.json','deployment_selection_locked.json'],visual_review_required_after_render=True)
        print(json.dumps(dict(state='ready_waiting_blind',gpu_used=False,conclusions_generated=False)),flush=True);return
    bundle,missing=load_validation_failure() if args.validation_only else load_inputs()
    if missing:state('waiting_for_blind_and_locked_results',missing=missing,gpu_used=False,report_generated=False);print(json.dumps(dict(state='waiting',missing=missing)));return
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');dest=FIG/('run_'+stamp);dest.mkdir(parents=True,exist_ok=False)
    state('rendering_cpu',gpu_used=False,run=str(dest.relative_to(P)))
    target,figures=render_validation_failure(bundle,dest) if args.validation_only else render(bundle,dest)
    write(dest/'inputs_used.json',dict(files=[dict(path=str(p.relative_to(P)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in bundle['required']],timing_kind=bundle['timing_kind'],report_path=str(target.relative_to(P))))
    state('generated_unreviewed',gpu_used=False,report_generated=True,report_scope='validation_exploration_goal_not_met' if args.validation_only else 'blind_accuracy_and_deployment_time',final_blind_test_claim=not args.validation_only,report_path=str(target.relative_to(P)),figures=[str(p.relative_to(P)) for pair in figures for p in pair],visual_review_required=True,
        note='Root must open the real PNG/PDF before delivery;oldreport retained')
    print(json.dumps(dict(state='generated_unreviewed',report=str(target),visual_review_required=True)),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException as error:
        state('failed',error=repr(error),traceback=traceback.format_exc(),gpu_used=False);raise
