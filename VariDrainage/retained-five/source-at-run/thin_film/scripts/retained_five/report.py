"""Terminal-only actual-data report; failures and missing comparisons stay visible."""
from pathlib import Path
import argparse,html,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from run import O,RAW,REF,read,save,sha
import archive
import compare_terminal as compare
BLUE='#3C5488';RED='#E64B35'

def eligible():
    state=archive.eligible();a=read(O/'archive_verified.json')
    assert a['all_restored_byte_exact'] and a['manifest_sha256']==sha(O/'download_manifest.json')
    assert not (O/'index.html').exists(),'Preserve the existing terminal report'
    return state,a

def main():
    state,arc=eligible();rows=[]
    if (O/'comparisons.json').exists():
        compare.eligible();rows=read(O/'comparisons.json');assert len(rows)==3
    coverage=read(O/'audit_coverage.json') if (O/'audit_coverage.json').exists() else {}
    native_pass=len(state['cases'])==2 and all(c['native_gate'] for c in state['cases'])
    summary=dict(native_pair_passed=native_pass,new_native_cases=len(state['cases']),skipped=state['skipped'],accepted_steps=sum(c.get('summary',{}).get('steps',0) for c in state['cases']),native_failures=sum(c.get('summary',{}).get('failures',int(c['returncode']!=0)) for c in state['cases']),native_batch_wall_s=state['compute_wall_s'],scalar_saved_audits=coverage.get('scalar_saved_steps',0),mechanical_saved_audits=coverage.get('mechanical_saved_steps',0),comparisons_available=bool(rows),comparisons=rows,original50_200us_complete=False,continuous_physics_validated=False,full_C1_complete=False)
    summary['missing_native_summaries']=[c['name'] for c in state['cases'] if 'summary' not in c]
    summary['accepted_steps_scope']='Counts proven by available native terminal summaries only; missing summary prevents a full accepted-step claim'
    save(O/'summary.json',summary)
    plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    comparison_figure=''
    if rows:
        fig,axes=plt.subplots(1,2,figsize=(10.8,4.5),layout='constrained')
        for ax,kind,gate,title in [(axes[0],'time',2,'Time-step halving at dx=5 um'),(axes[1],'space',10,'Same-dt spatial comparison:10 -> 5 um')]:
            group=[r for r in rows if r['kind']==kind];y=np.arange(len(group));w=.32
            for offset,key,color,label in [(-w/2,'dh',BLUE,'Thickness change'),(w/2,'displacement',RED,'Displacement')]:
                bars=ax.barh(y+offset,[100*r[key] for r in group],w,color=color,label=label);ax.bar_label(bars,fmt='%.3f',padding=3,fontsize=9)
            ax.axvline(gate,color='#555f69',ls='--',label=f'Original {gate}% gate');ax.set_yticks(y,['dx5 um'] if kind=='time' else ['dt=.0625 us','dt=.03125 us']);ax.invert_yaxis();ax.set(xlabel='Original moving-overlap difference (%)',title=title);ax.margins(x=.3);ax.grid(axis='x',alpha=.15)
        axes[0].legend(frameon=False,fontsize=9,loc='best');fig.suptitle('Same candidate at10us; all comparison failures retained',fontsize=12)
        for ext in ['png','pdf']:fig.savefig(O/('time-space.'+ext),dpi=200)
        plt.close(fig)
        comparison_figure='<figure><a href="time-space.png"><img src="time-space.png" alt="原范数下5微米时间配对及10到5微米空间比较"></a><figcaption>蓝色：膜厚变化；红色：位移；虚线：原验收标准。失败行完整保留。<a href="time-space.pdf">矢量PDF</a></figcaption></figure>'
    plots=[];snapshots=[]
    for rec in state['cases']:
        root=RAW/rec['name'];f=root/'states.csv'
        if not f.exists() or f.stat().st_size<200:continue
        states=np.atleast_1d(np.genfromtxt(f,delimiter=',',names=True));saved=states[states['stage']==2]
        if not len(saved):continue
        initial=np.atleast_1d(np.genfromtxt(root/'initial_states.csv',delimiter=',',names=True));assert len(initial)==2112
        complete=[int(step) for step in sorted(set(saved['step']),reverse=True) if len(saved[saved['step']==step])==2112]
        if not complete:continue
        step=complete[0];final=saved[saved['step']==step]
        assert np.array_equal(final['id'],initial['id'])
        dt=rec['parameters']['dt_us'];film=final['bath']==0;dh=(final['h_m']-initial['h_m'])*1e9
        snapshots.append(dict(case=rec['name'],saved_step=step,time_us=step*dt,last_accepted_step=rec.get('summary',{}).get('steps'),separate_failed_state_file_exists=(root/'failed_states.csv').exists(),scope='Latest complete accepted saved frame only, never labeled as failure geometry'))
        plots.append((rec,step,final,film,dh))
    particles=''
    if plots:
        limit=max(float(max(abs(dh[film]))) for _,_,_,film,dh in plots);assert limit>0
        cmap=LinearSegmentedColormap.from_list('nature_white',[BLUE,'#f5f5f3',RED]);fig,axes=plt.subplots(1,len(plots),figsize=(5.4*len(plots),4.2),layout='constrained',squeeze=False)
        for ax,(rec,step,final,film,dh) in zip(axes[0],plots):
            sc=ax.scatter(final['x_m'][film]*1e6,final['y_m'][film]*1e6,c=dh[film],cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit),s=8,linewidths=0);ax.scatter(final['x_m'][~film]*1e6,final['y_m'][~film]*1e6,c='#a0a7ae',s=8,linewidths=0);ax.axvline(240,color=BLUE,ls='--',lw=1);ax.set(xlabel='x(um)',ylabel='y(um)',aspect='equal',title=f"dt={rec['parameters']['dt_us']}us; saved step{step}, t={step*rec['parameters']['dt_us']:.5g}us")
        fig.colorbar(sc,ax=list(axes[0]),label='Thickness change from own initial state(nm)',shrink=.75);fig.suptitle('Actual saved nodes; N=2112; no interpolation or motion amplification',fontsize=12)
        for ext in ['png','pdf']:fig.savefig(O/('actual-particles.'+ext),dpi=200)
        plt.close(fig)
        particles='<figure><a href="actual-particles.png"><img src="actual-particles.png" alt="实际终态算例最后完整保存帧"></a><figcaption>逐点画真实粒子，共用膜厚变化色标；灰色为浴液。每图标注保存步和时刻，时刻不同不能作为空间收敛证据。失败算例的保存帧不等于失败现场。<a href="actual-particles.pdf">矢量PDF</a></figcaption></figure>'
    save(O/'snapshot_scope.json',snapshots)
    title='计算完成，继续看时间和空间能否过关。' if native_pass else '本批未完成，通过与失败分别保留。'
    failures=''.join('<li><b>'+html.escape(c['name'])+'</b>：'+html.escape(json.dumps(dict(returncode=c['returncode'],timed_out=c['timed_out'],native_gate=c['native_gate'],failure=c.get('failure')),ensure_ascii=False))+'</li>' for c in state['cases'])
    skipped=''.join('<li>'+html.escape(c['name']+'：'+c['reason'])+'</li>' for c in state['skipped'])
    table=''.join(f'<tr><td>{"5微米时间配对" if r["kind"]=="time" else "10→5微米／"+r["candidate"]}</td><td>{100*r["dh"]:.3f}%</td><td>{100*r["displacement"]:.3f}%</td><td>{"通过" if r["passed"] else "失败"}</td></tr>' for r in rows)
    manifest=read(O/'download_manifest.json');downloads=''.join(f'<li><a href="{part["path"]}">原始数据分卷{i+1}</a></li>' for i,part in enumerate(manifest['parts']))
    comparison_text='<table><tr><th>比较</th><th>膜厚变化差异</th><th>位移差异</th><th>原标准</th></tr>'+table+'</table>' if rows else '<p class="note">配对比较不可用：尚无满足完整原生计算、独立审计和原范数输入要求的结果。不能把缺失写成零误差或通过。</p>'
    css='*{box-sizing:border-box}body{margin:0;background:#fff;color:#213345;font:17px/1.9 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1100px;margin:auto;padding:30px 24px}h1{font-size:clamp(30px,5vw,46px);line-height:1.4}section{padding:24px 0;border-bottom:1px solid #d9dfe2}a{color:#3c5488}img{width:100%;height:auto}figure{margin:24px 0}figcaption{font-size:14px;color:#687584}.note{padding:20px;border-left:3px solid #e64b35;background:#faf3f1}table{width:100%;border-collapse:collapse;font-size:14px}td,th{text-align:left;padding:10px 5px;border-bottom:1px solid #d9dfe2}li,code{overflow-wrap:anywhere}'
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VariDrainage · 候选5微米配对</title><style>{css}</style></head><body><main><a href="https://baiming-zhang.github.io/VariDrainage/">项目首页</a><section><h1>{title}</h1><p>同一保留短边候选，2112粒子，两个计划时间步，目标时长10微秒。原程序、材料和验收标准未改变，旧失败没有重跑或覆盖。</p><p>原生终态摘要中已确认接受{summary['accepted_steps']}步，运行{summary['new_native_cases']}个原生算例，总耗时{summary['native_batch_wall_s']/3600:.2f}小时。独立保存步审计覆盖：守恒{summary['scalar_saved_audits']}项、机械功{summary['mechanical_saved_audits']}项；原生通过与独立审计分别记录。缺少终态摘要的算例：{html.escape(str(summary["missing_native_summaries"]))}；缺失时不能确认其完整接受步数。</p><ul>{failures}{skipped}</ul></section><section><h2>时间与空间：都需要证据</h2><p>时间步减半的原标准为2%，粒子加密的原标准为10%，膜厚变化和位移两项都要满足。差异是离散解之间的比较，不是相对实验真值的误差。</p>{comparison_figure}{comparison_text}<p><a href="https://baiming-zhang.github.io/VariDrainage/retained-reference/">此前同候选20／10微米结果</a>中的空间失败仍保留，不重复计算或隐藏。</p></section><section><h2>实际保存的粒子状态</h2>{particles if particles else '<p>本批没有可用于画图的完整阶段2保存帧，不生成替代粒子图。</p>'}<p><a href="snapshot_scope.json">各图实际步数与时刻</a>。内部240微米切面是账本切面，尚非已验证的物理接触线。</p></section><section><h2>接下来仍要回答什么</h2><p>先依据原生失败或时间／空间差异判断下一步，不能通过放宽阈值换取成功。10微秒只是首检查点：50／200微秒、连续物理通量和牵引力、独立交换系数与薄膜元厚度比、整体条纹及主实验／保留工况验证仍未完成。</p></section><section><details><summary>完整证据、下载与复现</summary><p><a href="summary.json">精确摘要</a> · <a href="status_at_native_completion.json">原生终态</a> · <a href="protocol.json">原生冻结协议</a> · <a href="postprocess_protocol.json">独立审计协议</a> · <a href="comparison_protocol.json">原范数比较协议</a> · <a href="archive_report_protocol.json">归档与报告协议</a>。</p><p>{manifest['file_count']}份原始输出、关闭日志、已有审计与源码，共{manifest['raw_bytes']/1e6:.3f}MB，已从{len(manifest['parts'])}个小于5MB的ZIP实际恢复，逐文件字节核验通过。<a href="download_manifest.json">字节清单</a> · <a href="archive_verified.json">恢复核验</a> · <a href="restore_raw.py">恢复脚本</a>。</p><ul>{downloads}</ul><p>恢复：<code>python restore_raw.py REPORT_DIRECTORY NEW_OUTPUT_DIRECTORY</code>，拒绝覆盖现有目录。</p></details></section></main></body></html>'''
    save(O/'status_at_native_completion.json',state);(O/'index.html').write_text(page,encoding='utf-8')
    save(O/'report_source.json',dict(source_sha256=sha(Path(__file__)),summary_sha256=sha(O/'summary.json'),archive_manifest_sha256=sha(O/'download_manifest.json'),full_C1_complete=False))
    print('Rendered actual terminal report; visual inspection still required before publication.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--run',action='store_true');g.add_argument('--readiness',action='store_true');a=parser.parse_args()
    if a.run:main()
    else:
        try:eligible()
        except (AssertionError,FileNotFoundError) as e:print(json.dumps(dict(ready=False,reason=str(e),report_exists=(O/'index.html').exists(),new_native_executions=0)))
        else:print(json.dumps(dict(ready=True,new_native_executions=0)))
