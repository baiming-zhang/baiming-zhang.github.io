"""Report a preserved terminal failure, never treat its prefix as a passed case."""
from pathlib import Path
import json, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from run import R,O,RAW,save,sha
import archive_case

NAME='u5_dt003125';D=O/'terminal-evidence'/NAME
BLUE='#3C5488';RED='#E64B35'

def main():
    state,rec=archive_case.terminal(NAME)
    assert state['state']=='completed' and not rec['native_gate'] and rec['returncode']!=0
    assert not (D/'failure_summary.json').exists()
    archive=json.loads((D/'archive_verified.json').read_text())
    assert archive['all_restored_byte_exact'] and archive['manifest_sha256']==sha(D/'download_manifest.json')
    assert json.loads((D/'case_record.json').read_text())['cases'][0]==rec
    read=lambda name:np.atleast_1d(np.genfromtxt(RAW/NAME/name,delimiter=',',names=True))
    initial=read('initial_states.csv');states=read('states.csv');hist=read('history.csv');force=read('mechanics.csv');m=rec['summary'];failure=rec['failure']
    assert len(hist)==len(force)==m['steps']==failure['last_accepted_step']==77
    assert m['requested_steps']==320 and failure['step']==78 and m['failures']==1
    last=int(max(states['step']));snapshot=states[(states['step']==last)&(states['stage']==2)]
    assert last==40 and len(snapshot)==len(initial)==m['N'] and np.array_equal(initial['id'],snapshot['id'])
    film=snapshot['bath']==0;dh=(snapshot['h_m']-initial['h_m'])*1e9
    maximum_force=float(max(force['g1_N']/force['target_N']))
    summary=dict(case=NAME,native_gate=False,native_returncode=rec['returncode'],accepted_steps=m['steps'],requested_steps=m['requested_steps'],failed_step=failure['step'],accepted_time_us=m['T_s']*1e6,last_saved_step=last,last_saved_time_us=last*m['dt_s']*1e6,failure_message=failure['message'],native_driver_wall_s=rec['wall_s'],maximum_accepted_VS_residual=max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res']),maximum_accepted_force_target_ratio=maximum_force,maximum_accepted_energy_budget_excess=m['max_budget_excess'],failed_state_exported=False,independent_failed_case_audit_executed=False,time_pair_available=False,space_comparison_qualified=False,full_C1_complete=False,scope='Failed native trajectory; accepted-prefix native diagnostics only; saved snapshot at1.25us is NOT the failure state')
    save(D/'failure_summary.json',summary)
    save(O/'batch_outcome.json',dict(native_terminal_cases=2,passed_native_cases=1,failed_native_cases=1,accepted_steps=sum(c['summary']['steps'] for c in state['cases']),time_pair_available=False,space_comparison_qualified=False,comparison_executed=False,reason='Fine-dt native case failed at attempted step78 after77accepted; did not reach10us',original_checkpoints_us=[10,50,200],full_C1_complete=False))
    plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    cmap=LinearSegmentedColormap.from_list('nature_white',[BLUE,'#f5f5f3',RED]);fig,ax=plt.subplots(figsize=(9.5,4.4),layout='constrained');limit=float(max(abs(dh[film])))
    sc=ax.scatter(snapshot['x_m'][film]*1e6,snapshot['y_m'][film]*1e6,c=dh[film],s=9,cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit),linewidths=0)
    ax.scatter(snapshot['x_m'][~film]*1e6,snapshot['y_m'][~film]*1e6,c='#a0a7ae',s=9,linewidths=0);ax.axvline(240,color=BLUE,ls='--',lw=1);ax.set(xlabel='x (um)',ylabel='y (um)',aspect='equal',title='Last SAVED state: step40, t=1.25 us; NOT failure state at step78')
    fig.colorbar(sc,ax=ax,label='Thickness change from own initial state (nm)',shrink=.8)
    for ext in ['png','pdf']:fig.savefig(D/('last-saved.'+ext),dpi=200)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained');balance=np.maximum.reduce([hist[k] for k in ['V_res','S_res','regionV_res','regionS_res']]);ratio=force['g1_N']/force['target_N']
    axes[0].semilogy(hist['step'],np.maximum(balance,1e-18),color=BLUE);axes[0].axhline(1e-10,color=RED,ls='--',label='Original ledger gate');axes[0].set(ylabel='Native V/S residual')
    axes[1].semilogy(force['step'],np.maximum(ratio,1e-18),color=BLUE);axes[1].axhline(1.000001,color=RED,ls='--',label='Original force gate');axes[1].set(ylabel='Force residual / own target')
    for ax in axes:ax.set_xlabel('Accepted native step');ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.15)
    fig.suptitle('77 accepted steps only; geometry failed on attempted step78',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(D/('accepted-prefix.'+ext),dpi=200)
    plt.close(fig)
    manifest=json.loads((D/'download_manifest.json').read_text());links=''.join(f'<li><a href="{p["path"]}">原始数据分卷{i+1}</a></li>' for i,p in enumerate(manifest['parts']))
    css='*{box-sizing:border-box}body{margin:0;background:#fff;color:#213345;font:17px/1.9 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1060px;margin:auto;padding:30px 24px}h1{font-size:clamp(30px,5vw,48px);line-height:1.4}h2{font-size:26px}a{color:#3c5488}section{padding:26px 0;border-bottom:1px solid #d9dfe2}.note{padding:20px;border-left:3px solid #e64b35;background:#faf3f1}img{width:100%;height:auto}figure{margin:25px 0}figcaption{font-size:14px;color:#687584}code,pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:14px}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}.card{padding:20px;background:#f4f7fa;border-top:3px solid #e64b35}.card b{display:block;font-size:28px}@media(max-width:650px){.cards{grid-template-columns:1fr}body{font-size:16px}}'
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VariDrainage · 第二组失败现场</title><style>{css}</style></head><body><main><a href="https://baiming-zhang.github.io/VariDrainage/">项目首页</a><section><h1>第二组在第78步停止。<br>本轮尚不能判断收敛。</h1><p>5微米网格、2112节点、时间步0.03125微秒。原定320步到10微秒，实际接受77步到2.40625微秒，随后原生几何检查报错退出。</p><div class="note"><b>这是失败报告。</b>首组160步完成不变，但两组没有在10微秒配齐，因此没有执行本轮时间／空间比较。没有放宽门槛、重启作业或把接受步前缀当成合格终态。</div></section><section><h2>哪里出了问题？</h2><p>求解器在为运动项计算二阶响应时，会观察轻微扰动后的粒子几何。这次“正向扰动”触发了Voronoi面几何一致性检查。日志指出法向闭合残差超出原几何门槛；这定位了报错环节，尚未证明根因或修复方案。</p><pre>{failure['message']}</pre><p>已接受步的液量／表活最大账本残差为{summary['maximum_accepted_VS_residual']:.3e}，能量预算超额最大为0；这些数值不能消除第78步的几何失败，也不能证明连续物理模型正确。</p></section><section><h2>真实保存到了哪里？</h2><div class="cards"><div class="card"><b>77 / 320</b>原生接受步</div><div class="card"><b>2.40625 μs</b>接受步推进终点</div><div class="card"><b>1.25 μs</b>最后完整保存的粒子帧</div></div><figure><img src="last-saved.png" alt="实际第40步1.25微秒的保存粒子，非失败状态"><figcaption>实际第40步节点，膜厚变化相对各自初态，未插值或放大位移；灰色为浴液节点。<b>没有导出失败扰动几何，也没有导出第77步完整粒子帧。</b><a href="last-saved.pdf">矢量PDF</a></figcaption></figure><figure><img src="accepted-prefix.png" alt="仅77个接受步的原生账本和力残差"><figcaption>仅绘制已接受步；红线是原门槛，对数显示下限10⁻¹⁸仅用于显示零值。第78步不属于通过记录。未对本失败算例执行独立保存步审计。<a href="accepted-prefix.pdf">矢量PDF</a></figcaption></figure></section><section><h2>下一步怎么做？</h2><p>先查阅稳健Voronoi构造和数值几何资料，再用有界诊断定位触发检查的几何。任何修复都需保持原几何、守恒和动力学门槛，并在独立记录中验证；不会直接改阈值后重跑长作业。原10／50／200微秒、物理通量／牵引力和完整C1目标保留。</p><p><a href="https://baiming-zhang.github.io/VariDrainage/uniform5-first/">首组完成结果</a> · <a href="https://baiming-zhang.github.io/VariDrainage/uniform5-flux/">首组局部切面诊断</a>。先前局部一致性改善不替代此次失败的时间配对。</p></section><section><details><summary>失败现场、完整数据与复现依据</summary><p><a href="failure_summary.json">失败摘要</a> · <a href="case_record.json">精确原生终态记录</a> · <a href="native_protocol.json">原始参数、源码与二进制</a> · <a href="logs/index.html">关闭后的完整日志</a>。本组耗时{rec['wall_s']/3600:.3f}小时。</p><p>{manifest['file_count']}份文件／日志、{manifest['raw_bytes']/1e6:.6f}MB，已从{len(manifest['parts'])}卷ZIP实际恢复，逐文件字节核验一致。<a href="download_manifest.json">字节清单</a> · <a href="archive_verified.json">恢复核验</a> · <a href="restore_raw.py">恢复脚本</a>。</p><ul>{links}</ul><p>恢复：<code>python restore_raw.py REPORT_DIRECTORY NEW_OUTPUT_DIRECTORY</code>。拒绝覆盖现有目录。完整原生终态与失败保留；处理文本可能经Git换行规范化，ZIP内为原始字节。</p></details></section></main></body></html>'''
    (D/'index.html').write_text(page,encoding='utf-8')
    shutil.copy2(Path(__file__),D/'failure_report.py')
    save(D/'report_source.json',dict(source_sha256=sha(Path(__file__)),archived_manifest_sha256=sha(D/'download_manifest.json'),native_status_sha256=sha(O/'status.json'),figure_sha256={f.name:sha(f) for f in D.iterdir() if f.suffix in ['.png','.pdf']},scope=summary['scope']))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
