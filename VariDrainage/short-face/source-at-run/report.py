"""Describe actual fixed observations, archive raw bytes, never claim a repair."""
from pathlib import Path
import csv,json,shutil,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import R,O,RAW,LOG,HERE,sha,save
BLUE='#3C5488';RED='#E64B35'

def main():
    assert not (O/'summary.json').exists()
    protocol=json.loads((O/'protocol.json').read_text(encoding='utf-8'));status=json.loads((O/'status.json').read_text(encoding='utf-8'));assert status['state']=='completed'
    assert all(sha(R/f)==h for f,h in protocol['original_sources'].items())
    assert sha(O/'source-at-build/run.py')==protocol['sources'][str(HERE.relative_to(R)/'run.py')]
    manifest=json.loads((O/'raw_manifest.json').read_text(encoding='utf-8'));assert all(sha(R/f)==v['sha256'] for f,v in manifest.items())
    records=list(csv.DictReader((RAW/'observations.csv').open(encoding='utf-8')));assert [r['case'] for r in records]==[c['name'] for c in protocol['cases']]
    for row in records:
        for k in list(row):
            if k not in ['case','message']:row[k]=float(row[k])
    base,short,resolved=records
    pattern=bool(base['original_completed']==resolved['original_completed']==1 and short['original_completed']==0 and short['original_normal_closure']>1e-9 and short['central_filtered_closure']>1e-9 and short['central_full_closure']<1e-12 and short['short_edge_count']==1)
    summary=dict(expected_fixed_counterexample_pattern_observed=pattern,observations=records,new_native_trajectory_steps=0,native_observation_wall_s=json.loads((O/'attempts.json').read_text(encoding='utf-8'))[-1]['wall_s'],actual_failed_step78_reproduced=False,solver_repair_adopted=False,dynamic_convergence_qualified=False,full_C1_complete=False,scope='A manufactured ORIGINAL-geometry counterexample; polygon full-edge sum is a telescoping diagnostic, not a repaired operator or physical validation')
    save(O/'summary.json',summary)
    native=O/'native';native.mkdir();logs=O/'logs';logs.mkdir()
    for f in RAW.iterdir():shutil.copy2(f,native/f.name)
    for f in LOG.iterdir():shutil.copy2(f,logs/f.name)
    source=O/'source-at-run';source.mkdir()
    for f in [HERE.parents[1]/'bench/tf_short_face_probe.cpp',HERE/'xmake.lua',HERE/'PLAN.md',HERE/'run.py',HERE/'resume_after_build.py',Path(__file__)]:shutil.copy2(f,source/f.name)
    members={};archive=O/'downloads/raw.zip';archive.parent.mkdir()
    candidates=[*(('native/'+f.name,f) for f in RAW.iterdir()),*(('logs/'+f.name,f) for f in LOG.iterdir()),*(('source-at-build/'+f.name,f) for f in (O/'source-at-build').iterdir())]
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,f in candidates:z.writestr(name,f.read_bytes());members[name]=dict(original_path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=sha(f))
    target=R/'results/local-runs/short_face_20261113/restored';assert not target.exists();target.mkdir()
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==set(members)
        for name,rec in members.items():
            p=target/name;p.parent.mkdir(exist_ok=True,parents=True);p.write_bytes(z.read(name));assert sha(p)==rec['sha256']
    save(O/'download_manifest.json',dict(archive='downloads/raw.zip',archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),members=members))
    save(O/'archive_verified.json',dict(files=len(members),all_restored_byte_exact=True,restore_directory=str(target.relative_to(R)),manifest_sha256=sha(O/'download_manifest.json')))
    plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    sites=list(csv.DictReader((RAW/'sites.csv').open(encoding='utf-8')));edges=list(csv.DictReader((RAW/'central_edges.csv').open(encoding='utf-8')))
    fig,axes=plt.subplots(1,2,figsize=(10,4.3),layout='constrained');s=[r for r in sites if r['case']=='short_face'];axes[0].scatter([float(r['x']) for r in s],[float(r['y']) for r in s],s=50,color=BLUE)
    e=[r for r in edges if r['case']=='short_face'];delta=short['delta']
    for row in e:
        x=np.array([float(row['ax']),float(row['bx'])]);y=np.array([float(row['ay']),float(row['by'])]);omitted=int(row['omitted_by_original_cutoff'])==1
        axes[0].plot(x,y,color=RED if omitted else BLUE,lw=2)
        axes[1].plot((x-.5)/delta,(y-.5)/delta,color=RED if omitted else BLUE,lw=3 if omitted else 1.5,marker='o',ms=5)
    axes[0].set(xlabel='x / site spacing',ylabel='y / site spacing',title='Actual nine manufactured sites',aspect='equal',xlim=(-1.3,1.3),ylim=(-1.3,1.3))
    axes[1].axhline(0,color='#999999',ls='--',lw=.7);axes[1].axvline(0,color='#999999',ls='--',lw=.7);axes[1].set(xlabel='(x - 0.5) / delta',ylabel='(y - 0.5) / delta',title='Actual corner coordinates: explicit geometric zoom',aspect='equal',xlim=(-.8,.15),ylim=(-.8,.15))
    for ext in ['png','pdf']:fig.savefig(O/('geometry.'+ext),dpi=200)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.4,4.5),layout='constrained');x=np.arange(3)
    ax.semilogy(x,np.maximum([r['original_normal_closure'] for r in records],1e-18),'o-',color=RED,label='Original native assembled closure')
    ax.semilogy(x,np.maximum([r['central_full_closure'] for r in records],1e-18),'s--',color=BLUE,label='Central polygon full-edge identity only')
    ax.axhline(1e-9,color='#555f69',ls=':',label='Unchanged original closure gate');ax.set_xticks(x,['Unperturbed','Short edge','Resolved edge']);ax.set(ylabel='Normalized normal closure',title='One original rejection; neither gate nor native code changed');ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=9)
    for ext in ['png','pdf']:fig.savefig(O/('closure.'+ext),dpi=200)
    plt.close(fig)
    css='*{box-sizing:border-box}body{margin:0;background:white;color:#213345;font:17px/1.9 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1060px;margin:auto;padding:30px 24px}h1{font-size:clamp(30px,5vw,46px);line-height:1.4}h2{font-size:26px}section{padding:28px 0;border-bottom:1px solid #d9dfe2}a{color:#3c5488}.note{padding:20px;background:#faf3f1;border-left:3px solid #e64b35}img{width:100%;height:auto}figure{margin:24px 0}figcaption{font-size:14px;color:#687584}code{overflow-wrap:anywhere}'
    verdict='出现了预先指定的短边反例' if pattern else '没有出现完整预期模式，保留负结果'
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VariDrainage · 极短边几何诊断</title><style>{css}</style></head><body><main><a href="https://baiming-zhang.github.io/VariDrainage/">项目首页</a><section><h1>只用9个点，<br>检验一次几何失败候选。</h1><p>{verdict}。原生几何函数没有修改、原闭合门槛没有放宽，三组观测只执行一次，没有推进任何模拟时间步。</p><div class="note"><b>这是新构造的小反例，不是第78步的重放。</b>原失败扰动几何没有导出；本页证明原规则存在一个问题场景，尚未证明它就是长算例的确切根因，也没有采用修复。</div></section><section><h2>丢掉一条很短的边，会发生什么？</h2><p>对一个闭合多边形，把每条边的“长度×外法向”相加应互相抵消。原实现会略过短于某个界限的边；这可能破坏后续面组装的闭合。我们固定9个点，只改变右上角一点的位置，包含未扰动、极小扰动、较大扰动三个预定条件。</p><figure><img src="geometry.png" alt="九点实际构造和明确放大的局部角点，红色表示被省略的边"><figcaption>红色是原规则省略的实际短边。右图按已声明的delta归一化放大几何角点，不是实际膜运动或粒子位移放大。单位为点间距，无实验长度映射。<a href="geometry.pdf">矢量PDF</a></figcaption></figure><p>极小扰动组中，短边长度为{short['short_edge_total']:.6e}个点间距，小于原截止值{short['cutoff']:.6e}。省略后中央单元闭合残差{short['central_filtered_closure']:.6e}，与原生报错数值一致；完整多边形边向量和为{short['central_full_closure']:.3e}。</p></section><section><h2>两组控制通过，一组被原门槛拒绝</h2><figure><img src="closure.png" alt="三组真实原生闭合残差，极短边组超过原门槛"><figcaption>红色为实际原生组装的闭合残差；蓝色仅是中央多边形全边向量和恒等式。零值以10⁻¹⁸作为绘图下限，不修改原数据。蓝线不是合格的新求解器或物理验证。<a href="closure.pdf">矢量PDF</a></figcaption></figure><p>未扰动和较大扰动组的原生构造通过；极短边组被原闭合门槛10⁻⁹拒绝。原生观测耗时{summary['native_observation_wall_s']:.3f}秒。构建后曾修正驱动输出路径和UTF-8读取，全部原记录保留；原生观测没有重试。</p></section><section><h2>下一步检验候选修复，不直接重跑长作业</h2><p>候选修复必须保持几何闭合、面积、正交、守恒和动力学的原门槛；还要处理真实运动几何、导数和实际时间／空间比较。完整边向量和接近零，只解释这个构造，不代表修复通过。原10／50／200微秒及完整C1仍未完成。</p><p>已查阅<a href="https://www.cs.cmu.edu/~quake/robust.html">Shewchuk稳健谓词资料</a>和<a href="https://doc.cgal.org/latest/Kernel_23/index.html">CGAL几何内核文档</a>：分类判断与新坐标构造需要分别考虑；这些资料不能替代本项目的反例或动力学证据。</p><p><a href="https://baiming-zhang.github.io/VariDrainage/uniform5-second/">原第二组失败现场</a> · <a href="https://baiming-zhang.github.io/VariDrainage/uniform5-first/">首组真实终态</a>。</p></section><section><details><summary>原始数据、构建记录与完整来源</summary><p><a href="summary.json">完整观测与科学范围</a> · <a href="protocol.json">观测前固定方案</a> · <a href="attempts.json">实际命令和耗时</a> · <a href="native/observations.csv">原生观测CSV</a> · <a href="native/central_edges.csv">完整中央边坐标</a> · <a href="native/sites.csv">所有点坐标</a>。</p><p><a href="downloads/raw.zip">原始文件／关闭日志及原驱动源码ZIP</a> · <a href="download_manifest.json">逐文件SHA清单</a> · <a href="archive_verified.json">实际恢复核验</a>。全部{len(members)}个文件已实际恢复、字节一致。<a href="source-at-run/PLAN.md">固定计划</a> · <a href="source-at-run/tf_short_face_probe.cpp">原生探针源码</a> · <a href="driver_path_patch.json">驱动路径修正</a> · <a href="metadata_encoding_patch.json">UTF-8修正记录</a>。ZIP保留原始字节，处理后的文本可能经Git换行规范化。</p></details></section></main></body></html>'''
    (O/'index.html').write_text(page,encoding='utf-8')
    save(O/'report_source.json',dict(source_sha256=sha(Path(__file__)),raw_manifest_sha256=sha(O/'raw_manifest.json'),summary_sha256=sha(O/'summary.json'),scope=summary['scope']))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
