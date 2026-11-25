"""Insert completed real-grid GPU/CPU results into the existing readable report."""
from pathlib import Path
import csv,json,html,re,shutil
from large_grid_operator_fixtures import ROOT,R,save

SOURCE=ROOT/'project-web-operator-ablation/PLIC-VOF/operator-ablation/wide-native'

def rows(name):
    with (R/'analysis'/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def num(r,key):return float(r[key]) if r.get(key) not in ('',None) else None
def fmt(v):return '未覆盖全部样本' if v is None else f'{float(v):.4g}'
def historical(g):return str(g['historical']).lower()=='true'
def name(g):
    if historical(g):return '此前训练的3×3联合网络'
    p='基准数据' if g['policy']=='control' else '调整后的数据'
    k=int(g['stencil']);return f'{k}×{k}网络 · {p}'
def spread(g,value,sd):
    text=fmt(num(g,value));s=num(g,sd)
    return text+(f' ± {fmt(s)}' if s is not None else '')
def table(headers,body):
    return '<div class="table-scroll"><table class="reader-table"><thead><tr>'+''.join('<th scope="col">'+h+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in body)+'</tbody></table></div>'
def figure(stem,caption):
    base='figures/'+stem
    return f'<figure><a href="{base}.png" target="_blank" rel="noopener"><img src="{base}.png" alt="{html.escape(caption)}" loading="lazy"></a><figcaption>{caption}</figcaption><div class="figure-links"><a href="{base}.png">原尺寸PNG</a><a href="{base}.pdf">矢量PDF</a></div></figure>'

def main():
    summary=json.loads((R/'analysis/summary.json').read_text())
    assert summary['state']=='complete' and summary['measured_record_count']==1088
    groups=rows('gpu_three_seed_groups.csv');records=rows('all_measured_records.csv')
    targets=json.loads((R/'analysis/equal_time_equal_error.json').read_text())['records']
    body=['<!-- large-grid-results:start -->','<section id="gpu-deployment" class="hero latest-overview"><div class="eyebrow">最新补充 · 大网格 · GPU神经网络对CPU传统数值算法</div><h1>换成GPU、大网格以后，结果如何？</h1>',
        '<p class="lead">这里直接比较GPU上的纯神经网络和CPU上的SimLiquid、APHROS局部算法。主结果包含数据传到GPU、推理和结果传回CPU的时间。</p>']
    conclusions=[];highlights=[]
    for op,cn in [('flux','通量'),('curvature','曲率')]:
        cpu=[r for r in records if r['scene']=='many_drops' and int(r['N'])==4096 and r['operator']==op and r['device']=='CPU' and num(r,'full_task_coverage')==1]
        if not cpu:
            conclusions.append(f'{cn}：传统算法存在未覆盖样本，不能用部分样本误差替代完整任务精度。');continue
        best=min(cpu,key=lambda r:num(r,'time_median_ms'))
        target=next(t for t in targets if t['scene']=='many_drops' and t['N']==4096 and t['operator']==op and t['mode']=='cpu_gpu_cpu' and t['cpu_backend']==best['implementation'])
        title='APHROS' if best['implementation']=='aphros' else 'SimLiquid'
        base=f'{cn}：在4096²多液滴场景，最快的已测CPU配置是{title}（{best["workers"]}线程），耗时{fmt(num(best,"time_median_ms"))}毫秒，误差{fmt(num(best,"full_task_mae"))}。'
        equal=target['equal_error_gpu']
        if equal:
            ratio=target['cpu_ms']/equal['time_mean_ms']
            timing=f'快{ratio:.2f}倍' if ratio>=1 else f'耗时为CPU的{1/ratio:.2f}倍'
            seed='三个种子都达到该误差目标' if target['all_gpu_seeds_reach_error'] else '三个种子的平均误差达标，但并非每个种子都达标'
            base+=f'按平均误差达到同一目标的最快神经组为{name(equal)}，{fmt(equal["time_mean_ms"])}毫秒，{timing}；{seed}。'
        else:base+='这批神经候选的平均误差均未达到该CPU精度目标，因此没有可报告的“同精度加速比”。'
        same=target['equal_time_gpu']
        if same:
            improvement=100*(1-same['full_task_mae']/target['cpu_error']) if target['cpu_error'] else None
            change=f'误差低{improvement:.1f}%' if improvement is not None and improvement>=0 else (f'误差高{-improvement:.1f}%' if improvement is not None else '')
            base+=f'在相同时间预算内，误差最低的已测神经组为{name(same)}，误差{fmt(same["full_task_mae"])}，{change}。'
        else:base+='没有神经组进入该CPU的相同时间预算。'
        conclusions.append(base)
        if equal:
            comparison=f'快{ratio:.1f}倍' if ratio>=1 else f'耗时为CPU的{1/ratio:.2f}倍'
            highlights.append(f'{cn}在达到CPU误差目标时，GPU{fmt(equal["time_mean_ms"])}毫秒，CPU{fmt(target["cpu_ms"])}毫秒，{comparison}。')
        else:highlights.append(f'{cn}：这批神经组尚未达到CPU的误差目标。')
    body.extend('<p class="reader-take">'+text+'</p>' for text in conclusions)
    body.append('<p class="overview-context">以上重点数字来自4096×4096、256个液滴的场景，实际处理127090个界面格。CPU在每份传统实现的1、8、24线程中比较实测时间；8、24线程是持久线程池适配器，不是整模拟器计时。新网络用三个种子，旧联合网络单列。GPU使用单精度输入及原模型算术，传统数值算法使用双精度。</p>')
    body.append(figure('actual_vof_fields',(R/'figures/actual_vof_fields-caption.txt').read_text(encoding='utf-8')))
    body.append(figure('large_grid_accuracy_time',(R/'figures/large_grid_accuracy_time-caption.txt').read_text(encoding='utf-8')))
    body.append(figure('large_grid_time_scaling',(R/'figures/large_grid_time_scaling-caption.txt').read_text(encoding='utf-8')))
    body.append('<h2>同一个场景，把各个方法放在一起看</h2><p>下表为4096²多液滴场景。“±”表示三个训练种子之间的标准差；单个旧网络与传统算法不填写±0。表中保留全部窗口和两种数据策略。</p>')
    for op,cn in [('flux','通量密度'),('curvature','物理曲率')]:
        body.append('<h3>'+cn+'</h3>')
        cpu=[r for r in records if r['scene']=='many_drops' and int(r['N'])==4096 and r['operator']==op and r['device']=='CPU']
        cpu.sort(key=lambda r:(r['implementation'],int(r['workers'])))
        native=[[('APHROS' if r['implementation']=='aphros' else 'SimLiquid')+f' · {r["workers"]}线程',fmt(num(r,'full_task_mae')),fmt(num(r,'time_median_ms')),fmt(num(r,'time_p05_ms'))+'–'+fmt(num(r,'time_p95_ms'))] for r in cpu]
        body.append(table(['CPU传统数值算法','平均绝对误差','时间中位数［毫秒］','5%–95%测量范围［毫秒］'],native))
        for mode,label in [('cpu_gpu_cpu','GPU：传入、推理、传回'),('resident','GPU：数据驻留设备')]:
            gs=[g for g in groups if g['scene']=='many_drops' and int(g['N'])==4096 and g['operator']==op and g['mode']==mode]
            gs.sort(key=lambda g:(not historical(g),g['policy'],int(g['stencil'])))
            if mode=='resident':body.append('<details><summary>另看数据已在GPU上的时间</summary>')
            body.append(table([label,'平均绝对误差','时间［毫秒］'],[[name(g),spread(g,'full_task_mae','full_task_seed_sd'),spread(g,'time_mean_ms','time_seed_sd_ms')] for g in gs]))
            if mode=='resident':body.append('</details>')
    case_rows=[[f'{c["N"]}×{c["N"]}',('单液滴' if c['scene']=='single_drop' else '256个液滴'),f'{c["total_grid_cells"]:,}',f'{c["active_interface_cells"]:,}'] for c in summary['cases']]
    body.append('<details><summary>全部网格与实际工作量</summary>'+table(['网格','场景','网格总格数','参与计算的界面格'],case_rows)+'</details>')
    body.append('<p>本轮所有方法都覆盖全部测试界面格，没有因失败输出删掉样本。CPU通量计时波动较大，因此表中同时给出5%–95%的测量范围；速度比是本机中位数的比较，不是与精度和硬件无关的算法常数。</p><p>只计单步局部算子：CPU从已准备的体积分数邻域到CPU输出；GPU主结果从已准备的锁页CPU缓冲区，经传输和推理，再回传CPU输出。邻域提取、模型加载、内存分配、压力求解和流体推进均在计时范围之外。每组先预热10次，再测量93次，表中单模型时间取中位数。本轮大网格场景均为圆形界面，其他界面类型的精度见后面的独立几何实验。</p><p>这轮固定已有模型，不重新训练。大窗口仍使用13万基础训练样本和12000步预算；这些结果不能证明训练数据已经足够，也不能作为大窗口充分训练后的性能上限。</p>')
    export_names=['all_measured_records.csv','coverage.csv','gpu_three_seed_groups.csv','gpu_vs_fastest_tested_cpu.csv']
    links=[]
    labels=['全部模型、种子与CPU线程结果','各网格完整任务覆盖率','全部神经组均值与波动','对最快已测CPU配置的比较']
    for filename,label in zip(export_names,labels):
        target=SOURCE/'tables'/('large-grid-'+filename);target.parent.mkdir(exist_ok=True);shutil.copyfile(R/'analysis'/filename,target)
        links.append(f'<a class="button" href="tables/{target.name}">{label} · CSV</a>')
    code=SOURCE/'code';code.mkdir(exist_ok=True)
    scripts=['large_grid_operator_fixtures.py','large_grid_gpu_cpu.py','large_grid_gpu_preflight.py','large_grid_analysis.py','large_grid_figures.py','large_grid_field_figure.py','large_grid_reader.py']
    for filename in scripts:shutil.copyfile(ROOT/'agent_experiments'/filename,code/filename)
    for filename in ['protocol.json','interpretation_notes.json']:shutil.copyfile(R/filename,code/filename)
    code_links=''.join(f'<a class="button" href="code/{filename}">{html.escape(filename)}</a>' for filename in scripts)
    body.append('<div class="download-row">'+''.join(links)+'</div><details><summary>计算、分析与绘图代码</summary><p>代码沿用本项目已选网络权重和原生C++算子；生成测试场不需要已有模拟输出。文件可单独保存，图表和CSV已直接内置于离线HTML。</p><div class="download-row">'+code_links+'</div></details></section><!-- large-grid-results:end -->')
    text=(SOURCE/'index.html').read_text(encoding='utf-8')
    text=re.sub(r'<!-- large-grid-results:start -->.*?<!-- large-grid-results:end -->','',text,flags=re.S)
    text=text.replace('<main class="shell">','<main class="shell">'+''.join(body),1)
    # Earlier evidence remains available, with its CPU-only scope explicit.
    text=text.replace('<h1>你关心的三个问题，结论是什么？</h1>','<h2>此前的收敛、窗口与同CPU实验</h2>',1)
    text=text.replace('扩大邻域和调整数据带来了局部精度收益；稳定的细网格收敛和相对传统数值方法的CPU推理优势，仍未实现。','此前的等预算实验显示部分局部精度改善，但未获得稳定的细网格收敛；下面的旧计时仅比较同一CPU，与上面的GPU部署结果分开阅读。')
    if '<a href="#gpu-deployment">' not in text:text=text.replace('<nav class="side" aria-label="报告目录">','<nav class="side" aria-label="报告目录"><a href="#gpu-deployment">大网格GPU对CPU</a>',1)
    text=text.replace('<title>9×9、11×11与原生数值算子 · PLIC-VOF</title>','<title>大网格GPU对CPU、窗口与收敛实验 · PLIC-VOF</title>')
    (SOURCE/'index.html').write_text(text,encoding='utf-8',newline='\n')
    for stem in ['actual_vof_fields','large_grid_accuracy_time','large_grid_time_scaling']:
        for ext in ['.png','.pdf']:shutil.copyfile(R/'figures'/(stem+ext),SOURCE/'figures'/(stem+ext))
    (R/'analysis/reader-conclusions.txt').write_text('\n\n'.join(conclusions),encoding='utf-8')
    home=ROOT/'project-web-single-html-download/PLIC-VOF/index.html'
    if home.exists():
        source=home.read_text(encoding='utf-8')
        overview='<!-- OPERATOR_ABLATION_REPORT --><section id="operator-ablation" class="latest-overview" style="border-top:3px solid #155da5"><div class="eyebrow" style="color:#155da5">最新结果 · 大网格GPU对CPU与纯神经算子消融</div><h2>曲率用GPU有明显速度收益；通量仍是CPU更快。</h2><div class="overview">'
        overview+='<div class="overview-row"><div><h3>网格加密，精度会一直提高吗？</h3><p>目前没有稳定收敛。此前独立几何测试中，通量停在约0.000914；3×3联合网络的物理曲率误差由N128的0.0722增至N512的0.1879。</p></div></div>'
        overview+='<div class="overview-row"><div><h3>GPU对CPU，精度与时间怎样？</h3><p>'+html.escape(' '.join(highlights))+'</p></div></div>'
        overview+='<div class="overview-row"><div><h3>扩大窗口、改数据有用吗？</h3><p>有局部收益：此前N128测试中，基准数据7×7扩大至11×11，曲率误差低31.9%；固定11×11调整数据再低26.4%。但所有窗口仍使用同样约13万训练样本，尚未测出数据量是否足够。</p></div></div></div><p class="overview-context">GPU重点数字来自4096²、256个圆滴的单步算子测试，实际处理127090个界面格，包含数据传入和结果传回。新网络保留三个种子；计时不含完整流体模拟。</p><div class="downloads"><a class="download" href="operator-ablation/wide-native/#gpu-deployment">查看完整图表与结果</a><a class="download" href="operator-ablation/wide-native/download.html">下载单文件离线HTML</a></div></section><!-- /OPERATOR_ABLATION_REPORT -->'
        source=re.sub(r'<!-- OPERATOR_ABLATION_REPORT -->.*?<!-- /OPERATOR_ABLATION_REPORT -->',overview,source,flags=re.S)
        from humanize_operator_names import document
        home.write_text(document(source),encoding='utf-8',newline='\n')
    save(R/'reader_status.json',dict(state='prepared_visual_check_pending',conclusions=conclusions,cpu_configurations=96,gpu_configurations=992,scope='GPU vs CPU single operator, real large static fields'))
    print('\n\n'.join(conclusions))

if __name__=='__main__':
    import sys
    sys.stdout.reconfigure(encoding='utf-8')
    main()
