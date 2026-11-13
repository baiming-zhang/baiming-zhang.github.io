"""Frozen physical-cut recovery on a completed case; never consume live frames."""
from pathlib import Path
import argparse, datetime, json, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import R, O, save, sha
import audit_terminal as audits
import compare_terminal as comparisons

HERE=Path(__file__).resolve().parent
D=R/'results/uniform5_cut_20261113'
BLUE='#3C5488'; RED='#E64B35'

def sources():
    paths=[Path(__file__),HERE/'CUT.md',HERE/'audit_terminal.py',HERE/'compare_terminal.py',HERE.parent/'prepared_motion/analyze.py',HERE.parent/'curved_edge/analyze.py']
    return {str(p.relative_to(R)):sha(p) for p in paths}

def prepare():
    assert not D.exists(), 'Preserve earlier attempts'
    audits.preflight()
    rec=audits.eligible('u5_dt00625')
    audit=O/'audits/cases/u5_dt00625'
    wrapper=json.loads((audit/'postprocess_verified.json').read_text())
    assert wrapper['passed'] and wrapper['sources']==audits.inventory()
    assert wrapper['inherited_seal_sha256']==comparisons.verify_audit(audit)
    protocol=json.loads((O/'protocol.json').read_text())
    refs=protocol['reference_cases']
    assert {r['name'] for r in refs}=={'u10_dt00625','u10_dt003125'}
    for r in refs:comparisons.verify_raw(r)
    D.mkdir()
    save(D/'protocol.json',dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=sources(),case=rec,references=refs,raw={f['path']:f['sha256'] for r in [rec]+refs for f in r['files']},audit_seal_sha256=sha(audit/'verified.json'),old_cut_sha256=sha(R/'results/global_space_20261113/physical_cut.json'),native_protocol_sha256=sha(O/'protocol.json'),timeout_s=120,scope='One terminal saved step; same-dt resolution diagnostic, not physical validation or convergence acceptance'))

def cut(prep, rec):
    initial,states,m=comparisons.read_case(prep,rec)
    step=round(10e-6/m['dt_s'])
    a,b,c=[states[(states['step']==step)&(states['stage']==stage)] for stage in range(3)]
    assert all(len(s)==m['N'] for s in [a,b,c])
    root=next(R/f['path'] for f in rec['files'] if Path(f['path']).name=='faces.csv')
    faces=prep.old.read(root);faces=faces[(faces['step']==step)&(faces['stage']==2)]
    mean=float(initial['S_mol'].sum()/(m['Gamma_inf_mol_m2']*prep.old.domain.area))
    q=prep.physical_cut(a,b,c,faces,m,mean)
    for key,suffix in [('relative','_m3_s'),('material','_m3_s'),('surface','_mol_s')]:
        q[key+'_single_step_relative_difference']=abs(q['recovered_'+key+suffix]-q['ledger_'+key+suffix])/max(1e-60,abs(q['ledger_'+key+suffix]))
    return dict(time_us=10,scope='One saved step, nonconforming recovery, not physical normal-flux validation or temporal L2.',**q)

def run():
    p=json.loads((D/'protocol.json').read_text())
    assert not (D/'physical_cut.json').exists()
    assert p['sources']==sources() and p['native_protocol_sha256']==sha(O/'protocol.json')
    assert all(sha(R/f)==h for f,h in p['raw'].items())
    assert p['old_cut_sha256']==sha(R/'results/global_space_20261113/physical_cut.json')
    assert p['audit_seal_sha256']==comparisons.verify_audit(O/'audits/cases/u5_dt00625')
    assert audits.eligible('u5_dt00625')==p['case']
    prep=comparisons.load();old=json.loads((R/'results/global_space_20261113/physical_cut.json').read_text())
    controls=[]
    for rec in p['references']:
        actual=cut(prep,rec);expected=old[rec['name']]
        errors={k:abs(v-expected[k])/max(abs(expected[k]),1e-60) for k,v in actual.items() if isinstance(v,(float,int))}
        assert max(errors.values())<1e-12, errors
        controls.append(dict(case=rec['name'],max_relative_reproduction_difference=max(errors.values())))
    result=cut(prep,p['case']);rows={**old,'u5_dt00625':result}
    save(D/'physical_cut.json',rows)
    save(D/'summary.json',dict(new_terminal_observations=1,new_native_runs=0,reference_reproductions=controls,case='u5_dt00625',same_dt_cases=['u20_dt00625','u10_dt00625','u5_dt00625'],new_case=result,full_C1_complete=False,physical_flux_validated=False,time_space_qualified=False,scope=p['scope']))
    plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    names=['u20_dt00625','u10_dt00625','u5_dt00625'];x=np.arange(3);width=.36
    for ax,key,suffix,label in zip(axes,['relative','material','surface'],['_m3_s','_m3_s','_mol_s'],['Relative liquid flux (nL/s)','Material liquid flux (nL/s)','Material surfactant flux (pmol/s)']):
        ax.bar(x-width/2,[rows[n]['ledger_'+key+suffix]*1e12 for n in names],width,color=BLUE,label='Crossing ledger')
        ax.bar(x+width/2,[rows[n]['recovered_'+key+suffix]*1e12 for n in names],width,color=RED,label='Recovered field')
        ax.axhline(0,color='#889099',lw=.7);ax.set_xticks(x,['20','10','5']);ax.set(xlabel='Particle spacing (um)',ylabel=label);ax.grid(axis='y',alpha=.15)
    axes[0].legend(frameon=False,fontsize=9)
    fig.suptitle('Same dt=0.0625 us; one saved crossing at 10 us; no continuum truth',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(D/('cut-flux.'+ext),dpi=200)
    plt.close(fig)
    cells=''.join('<tr><td>'+dx+'</td>'+''.join(f'<td>{100*rows[n][k+"_single_step_relative_difference"]:.3f}%</td>' for k in ['relative','material','surface'])+'</tr>' for dx,n in zip(['20','10','5'],names))
    css='*{box-sizing:border-box}body{margin:0;background:#fff;color:#213345;font:17px/1.9 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1080px;margin:auto;padding:30px 24px}h1{font-size:clamp(30px,5vw,46px);line-height:1.4}h2{font-size:25px;margin-top:40px}a{color:#3c5488}img{width:100%;height:auto}.note{border-left:3px solid #e64b35;padding:18px;background:#faf3f1}table{width:100%;border-collapse:collapse;font-size:15px}td,th{padding:10px 5px;border-bottom:1px solid #d9dfe2;text-align:left}figcaption{font-size:14px;color:#687584}figure{margin:24px 0}'
    page=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VariDrainage · 切面通量诊断</title><style>{css}</style><main><a href="../index.html">VariDrainage 项目首页</a><h1>网格更细后，<br>局部流量的两种算法接近了吗？</h1><p>同一条 x=240 微米切面、同一时间步 0.0625 微秒、同一时刻 10 微秒。新增 5 微米网格终态观测，复用已有 20／10 微米结果。</p><h2>一边数穿过去多少，一边从局部状态重建</h2><p>蓝色是运动单元穿过切面的离散账本；红色是用保存的膜厚、粒子速度与压力梯度重建的流量。“相对液流”表示液体相对粒子的流动，“随粒子输运”表示运动本身带过切面的液量或表面活性剂。正负号表示方向，不能只比较绝对大小。</p><figure><img src="cut-flux.png" alt="相同时间步的20、10、5微米网格带符号通量对照"><figcaption>白底红蓝；实际保存步，无重跑轨迹、拟合或新增物理参数。液量单位 nL/s，表活单位 pmol/s。<a href="cut-flux.pdf">矢量 PDF</a></figcaption></figure><h2>两种算法相差多少？</h2><table><thead><tr><th>网格／微米</th><th>相对液流</th><th>随粒子液量</th><th>随粒子表活</th></tr></thead><tbody>{cells}</tbody></table><p>差异 = |重建值 − 账本值| / |账本值|。这衡量两种离散观察的差异，不是相对真实实验的误差；分母接近零时也不宜单看百分比。</p><div class="note"><b>结论边界：</b>这是一个保存步的非相容场重建。两种观察涉及不同的时间层和运动单元穿越，接近或不接近都不能单独证明连续物理通量、牵引力或动态空间收敛。第二组更小时间步尚未纳入；没有改变原来的时间／空间验收门槛。</div><h2>下一步由完整配对决定</h2><p>先完成 5 微米网格的时间配对，再用原始重叠区域范数比较膜厚变化与位移。如果空间检查仍失败，就定位离散或物理通量问题；不会用单个切面结果替代验收，也不会仅因差异减小就放大计算规模。原 50／200 微秒及交换系数、薄膜元厚度比、整体条纹和保留工况实验仍待完成。</p><p><a href="../uniform5-first/">首组实际粒子、守恒与机械功报告</a> · <a href="../global-space/">已有空间失败和原始参考数据</a></p><details><summary>计算依据与复现</summary><p>两个已有 10 微米切面按原公式重算并核对；新增原生步数为 0，只处理已完成首组。<a href="summary.json">结果摘要</a> · <a href="physical_cut.json">全部带符号原始数值</a> · <a href="protocol.json">处理前冻结输入与源码哈希</a> · <a href="CUT.md">诊断计划</a> · <a href="cut_report.py">后处理代码</a>。</p><p>5 微米完整原始数据和独立审计在首组报告；旧参考数据保留在已有空间报告。完整 C1 与物理实验尚未完成。</p></details></main></html>'''
    (D/'index.html').write_text(page,encoding='utf-8')
    for src in [Path(__file__),HERE/'CUT.md']:shutil.copy2(src,D/src.name)
    assert all(sha(R/f)==h for f,h in p['raw'].items())
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',action='store_true');g.add_argument('--run',action='store_true');a=parser.parse_args();prepare() if a.prepare else run()
