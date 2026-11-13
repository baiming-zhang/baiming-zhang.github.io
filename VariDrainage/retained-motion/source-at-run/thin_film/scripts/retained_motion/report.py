"""Actual coupled one-step controls and original independent audit evidence."""
from pathlib import Path
import json,shutil,zipfile,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from run import R,O,RAW,LOG,HERE,read,save,sha
BLUE='#3C5488';RED='#E64B35'

def main():
    assert not (O/'index.html').exists()
    p=read(O/'protocol.json');state=read(O/'status.json');coverage=read(O/'audit_coverage.json')
    assert state['state']=='completed' and coverage['actual_single_step_controls_qualified'] and len(state['cases'])==3 and not state['skipped']
    assert all(sha(R/f)==h for f,h in p['sources'].items()) and all(sha(R/f)==h for f,h in p['original_sources'].items())
    original=read(R/'results/uniform5_pilot_20261113/protocol.json');assert sha(Path(original['executable']))==p['original_binary_sha256']
    summaries=[];plots=[]
    for rec in state['cases']:
        assert rec['passed'] and all(sha(R/f['path'])==f['sha256'] for f in rec['files'])
        d=O/'audits/cases'/rec['name'];seal=read(d/'verified.json');assert all(sha(d/n)==h for n,h in seal['outputs'].items())
        m=rec['summary'];eq=read(O/('equivalence_'+rec['name']+'.json'));errors=[v['normalized_max'] for c in eq['checks'].values() for v in c.get('errors',{}).values()];assert eq['passed']
        fa=read(d/'saved_formula_audit.json')['audits'][rec['name']];ma=read(d/'mechanical_audit.json')[rec['name']];assert fa['passed'] and len(ma)==1 and ma[0]['passed']
        root=RAW/rec['name'];initial=np.atleast_1d(np.genfromtxt(root/'initial_states.csv',delimiter=',',names=True));states=np.atleast_1d(np.genfromtxt(root/'states.csv',delimiter=',',names=True));final=states[(states['step']==1)&(states['stage']==2)];assert len(final)==len(initial)==m['N'];dh=(final['h_m']-initial['h_m'])*1e9;film=final['bath']==0;plots.append((final,dh,film,rec))
        summaries.append(dict(name=rec['name'],N=m['N'],dx_um=m['dx_m']*1e6,dt_us=m['dt_s']*1e6,native_gate=True,first_step_numerical_equivalence=True,maximum_compared_state_face_error=max(errors),scalar_saved_steps=len(fa['audited_saved_steps']),mechanical_saved_steps=len(ma),maximum_VS_residual=max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res']),force_target_ratio=rec['maximum_force_target_ratio'],maximum_energy_budget_excess=m['max_budget_excess'],native_wall_s=rec['wall_s'],minimum_film_h_m=float(min(final['h_m'][film])),maximum_abs_dh_m=float(max(abs(dh[film]))*1e-9)))
    summary=dict(cases=summaries,new_candidate_accepted_steps=3,independent_scalar_saved_steps=3,independent_mechanical_saved_steps=3,actual_single_step_controls_qualified=True,reused_candidate_fixture_checks=24,new_fixture_checks_executed_here=0,original_binary_sources_unchanged=True,old_reference_trajectories_replayed=0,actual_failed78_reproduced=False,time_space_qualified=False,full_C1_complete=False,scope='Three actual candidate single-step controls; no full trajectory equivalence or time/space/physical acceptance')
    save(O/'summary.json',summary)
    plt.rcParams.update({'font.size':10,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    cmap=LinearSegmentedColormap.from_list('nature_white',[BLUE,'#f5f5f3',RED]);limit=max(float(max(abs(dh[film]))) for _,dh,film,_ in plots)
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for ax,(final,dh,film,rec) in zip(axes,plots):
        sc=ax.scatter(final['x_m'][film]*1e6,final['y_m'][film]*1e6,c=dh[film],cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit),s=6,linewidths=0);ax.scatter(final['x_m'][~film]*1e6,final['y_m'][~film]*1e6,c='#a0a7ae',s=7,linewidths=0);ax.axvline(240,color=BLUE,ls='--',lw=.8);ax.set(xlabel='x (um)',ylabel='y (um)',aspect='equal',title=f'dx={rec["parameters"]["dx_um"]} um, t={rec["parameters"]["dt_us"]} us')
    fig.colorbar(sc,ax=axes,label='Thickness change from own initial state (nm)',shrink=.7);fig.suptitle('Actual candidate one-step saved nodes; different physical end times',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(O/('actual-particles.'+ext),dpi=200)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10.7,4.2),layout='constrained');labels=['5um / .0625us','5um / .03125us','10um / .0625us']
    for ax,key,gate,label in [(axes[0],'maximum_VS_residual',1e-10,'Native V/S residual'),(axes[1],'force_target_ratio',1.000001,'Force residual / own target')]:
        ax.semilogy(range(3),np.maximum([c[key] for c in summaries],1e-18),'o',color=BLUE);ax.axhline(gate,color=RED,ls='--',label='Original gate');ax.set_xticks(range(3),labels,rotation=15);ax.set_ylabel(label);ax.legend(frameon=False);ax.grid(alpha=.15)
    fig.suptitle('Three accepted steps; independent SI/GEOS audits reported separately',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(O/('native-checks.'+ext),dpi=200)
    plt.close(fig)
    source=O/'source-at-run';source.mkdir();files={**p['sources'],str(HERE/'audit_ready.py'):sha(HERE/'audit_ready.py'),str(Path(__file__)):sha(Path(__file__))}
    for name,h in files.items():
        f=Path(name) if Path(name).is_absolute() else R/name;assert sha(f)==h;dest=source/f.relative_to(R);dest.parent.mkdir(exist_ok=True,parents=True);shutil.copy2(f,dest)
    save(O/'source_index.json',dict(base_commit_before_uncommitted_controls=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),sources={str((Path(k) if Path(k).is_absolute() else R/k).relative_to(R)):v for k,v in files.items()},scope='Exact new sources snapshotted; old original sources/binary unchanged'))
    (O/'logs').mkdir()
    for f in LOG.iterdir():shutil.copy2(f,O/'logs'/f.name)
    members={};archive=O/'downloads/raw.zip';archive.parent.mkdir()
    candidates=[*(('native/'+f.relative_to(RAW).as_posix(),f) for f in RAW.rglob('*') if f.is_file()),*(('logs/'+f.name,f) for f in LOG.iterdir()),*(('audits/'+f.relative_to(O/'audits').as_posix(),f) for f in (O/'audits').rglob('*') if f.is_file())]
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,f in candidates:z.writestr(name,f.read_bytes());members[name]=dict(original_path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=sha(f))
    assert archive.stat().st_size<5_000_000
    target=RAW.parent/'restored';assert not target.exists();target.mkdir()
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==set(members)
        for name,item in members.items():
            f=target/name;f.parent.mkdir(exist_ok=True,parents=True);f.write_bytes(z.read(name));assert sha(f)==item['sha256']
    save(O/'download_manifest.json',dict(archive='downloads/raw.zip',archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),members=members))
    save(O/'archive_verified.json',dict(files=len(members),raw_bytes=sum(v['bytes'] for v in members.values()),all_restored_byte_exact=True,manifest_sha256=sha(O/'download_manifest.json'),restore_directory=str(target.relative_to(R))))
    rows=''.join(f'<tr><td>{c["dx_um"]:.0f}</td><td>{c["dt_us"]:.5f}</td><td>{c["N"]}</td><td>{c["maximum_compared_state_face_error"]:.1e}</td><td>{c["native_wall_s"]:.2f}s</td></tr>' for c in summaries)
    css='*{box-sizing:border-box}body{margin:0;background:#fff;color:#213345;font:17px/1.9 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1100px;margin:auto;padding:30px 24px}h1{font-size:clamp(30px,5vw,46px);line-height:1.4}h2{font-size:26px}section{padding:28px 0;border-bottom:1px solid #d9dfe2}a{color:#3c5488}.note{padding:20px;border-left:3px solid #e64b35;background:#faf3f1}img{width:100%;height:auto}figure{margin:24px 0}figcaption{font-size:14px;color:#687584}table{width:100%;border-collapse:collapse;font-size:14px}td,th{text-align:left;padding:10px 5px;border-bottom:1px solid #d9dfe2}'
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VariDrainage · 候选实际单步对照</title><style>{css}</style></head><body><main><a href="https://baiming-zhang.github.io/VariDrainage/">项目首页</a><section><h1>从几何反例，<br>走到实际算例的单步检查。</h1><p>保留短边候选在三组真实参数算例中各推进一步；原生门槛、第一步输出对照、独立液量／表活审计及机械功审计均通过。</p><div class="note"><b>这是单步控制，不是完整动力学验收。</b>没有重放原第78步失败，没有完成新的10微秒时间配对，也没有验收空间收敛或连续物理通量。</div></section><section><h2>真实粒子发生了什么变化？</h2><figure><img src="actual-particles.png" alt="三组实际候选单步保存节点与膜厚变化，终点物理时间不同"><figcaption>逐点绘制实际保存节点；颜色为相对各自初始状态的膜厚变化，三图共用色标，无插值或位移放大。灰色是浴液节点，虚线是内部切面；三组终点时间不同，不能把图像差异当作收敛误差。<a href="actual-particles.pdf">矢量PDF</a></figcaption></figure></section><section><h2>与旧数据的第一步相同吗？</h2><table><thead><tr><th>网格／μm</th><th>时间步／μs</th><th>节点</th><th>比较列最大归一化差异</th><th>原生耗时</th></tr></thead><tbody>{rows}</tbody></table><p>三组比较列的差异均为0。使用原有数值对照标准：离散结构精确一致，速度和面通量容差10⁻⁸，其余连续列10⁻¹⁰，位置用网格间距归一化。5微米参考来自旧单步算例，10微米参考只取旧160步轨迹保存的第一步。完整轨迹摘要不作对照，候选原生摘要另按原门槛检查。</p><p>本轮只执行候选；旧参考未重新计算。24项先前候选检查通过哈希复用，本页没有再次执行它们。</p></section><section><h2>账本和机械功是否对得上？</h2><figure><img src="native-checks.png" alt="三组单步原生收支和力残差与原门槛比较"><figcaption>蓝色是实际原生单步记录，红色是原门槛。独立SI／GEOS审计另外对每组保存步核对液量／表活与机械功，3组标量和3组机械审计均通过。<a href="native-checks.pdf">矢量PDF</a></figcaption></figure><p>这些是离散公式和单步一致性证据，不能证明所有拓扑变化附近的导数正确，也不能替代连续牵引力、实验或实际时间／空间收敛。原默认求解器、材料、容差和归档二进制均未修改。</p></section><section><h2>下一步：同一候选实现下的参考配对</h2><p>先计算候选的20／10微米参考算例并检查时间配对，再按结果决定5微米计算。不能把旧实现的粗时间步和候选实现的细时间步拼成合格配对。原10／50／200微秒、交换系数、薄膜元厚度比、整体条纹和保留工况实验仍待完成。</p><p><a href="https://baiming-zhang.github.io/VariDrainage/retained-faces/">此前几何与24项检查</a> · <a href="https://baiming-zhang.github.io/VariDrainage/uniform5-second/">原失败现场</a>。</p></section><section><details><summary>原始数据、独立审计和完整复现</summary><p><a href="summary.json">本页精确数值</a> · <a href="protocol.json">执行前协议与参考哈希</a> · <a href="status.json">原生记录</a> · <a href="audit_coverage.json">独立审计范围</a> · <a href="audit_attempts.json">审计命令</a>。</p><p><a href="downloads/raw.zip">完整原始输出、关闭日志及审计字节ZIP</a> · <a href="download_manifest.json">逐文件SHA清单</a> · <a href="archive_verified.json">实际恢复核验</a>。全部{len(members)}份文件已实际恢复，字节一致；每个原始文件和ZIP均小于5MB。<a href="source_index.json">源码SHA与基点</a> · <a href="source-at-run/thin_film/scripts/retained_motion/PLAN.md">固定计划</a> · <a href="source-at-run/thin_film/scripts/retained_motion/run.py">候选运行脚本</a> · <a href="source-at-run/thin_film/scripts/retained_motion/compare.py">第一步选择适配</a>。处理文本可能经Git换行规范化，ZIP保留原始字节。</p></details></section></main></body></html>'''
    assert all(v['bytes']<5_000_000 for v in members.values())
    (O/'index.html').write_text(page,encoding='utf-8');save(O/'report_source.json',dict(source_sha256=sha(Path(__file__)),summary_sha256=sha(O/'summary.json'),coverage_sha256=sha(O/'audit_coverage.json'),scope=summary['scope']))
    print('Reported3actual controls,3scalar+3mechanical independent audits; rawfiles',len(members),'ZIP bytes',archive.stat().st_size)
if __name__=='__main__':main()
