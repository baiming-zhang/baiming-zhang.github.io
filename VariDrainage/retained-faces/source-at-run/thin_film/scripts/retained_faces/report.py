"""Publish geometry and existing-fixture candidate admission, including limits."""
from pathlib import Path
import json,re,shutil,zipfile,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import R,O,LOCAL,LOG,HERE,read,save,sha
BLUE='#3C5488';RED='#E64B35'

def main():
    assert not (O/'index.html').exists()
    summary=read(O/'summary.json');p=read(O/'protocol.json');cdir=R/'results/retained_checks_20261113';checks=read(cdir/'summary.json');cp=read(cdir/'protocol.json')
    assert summary['fixed_geometry_matrix_admitted'] and checks['all_fixture_checks_passed'] and checks['new_unit_checks']==24
    assert cp['geometry_summary_sha256']==sha(O/'summary.json')
    assert all(sha(R/f)==h for f,h in p['original_sources'].items())
    assert all(sha(R/f)==h for f,h in p['sources'].items()) and all(sha(R/f)==h for f,h in cp['sources'].items())
    original=read(R/'results/uniform5_pilot_20261113/protocol.json');assert sha(Path(original['executable']))==original['binary_sha256']
    for a in read(O/'attempts.json')+read(cdir/'attempts.json'):assert sha(R/a['log'])==a['log_sha256']
    raw=read(O/'raw_manifest.json');assert all(sha(R/f)==v['sha256'] for f,v in raw.items())
    dst=O/'checks';shutil.copytree(cdir,dst)
    (O/'native').mkdir();(O/'logs').mkdir()
    for case in p['cases']:shutil.copytree(LOCAL/'native'/case['name'],O/'native'/case['name'])
    for f in LOG.iterdir():shutil.copy2(f,O/'logs'/f.name)
    (dst/'logs').mkdir()
    for a in read(cdir/'attempts.json'):shutil.copy2(R/a['log'],dst/'logs'/Path(a['log']).name)
    sources={**p['sources'],**cp['sources'],str(Path(__file__).relative_to(R)):sha(Path(__file__)),str(HERE/'analyze.py'):sha(HERE/'analyze.py')}
    source_dir=O/'source-at-run';source_dir.mkdir()
    for name,h in sources.items():
        f=Path(name) if Path(name).is_absolute() else R/name;assert sha(f)==h
        relative=f.relative_to(R);target=source_dir/relative;target.parent.mkdir(exist_ok=True,parents=True);shutil.copy2(f,target)
    source=R/'thin_film/include/vf/thinfilm/VoronoiDrainage.h';shutil.copy2(source,source_dir/'VoronoiDrainage.original.h')
    save(O/'source_index.json',dict(base_commit_before_uncommitted_candidate=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),sources={str((Path(n) if Path(n).is_absolute() else R/n).relative_to(R)):h for n,h in sources.items()},original_factory_sha256=sha(source),scope='New candidate files snapshotted by SHA; original captured source/binary unchanged'))
    members={};archive=O/'downloads/raw.zip';archive.parent.mkdir()
    originals=[*(('inputs/'+f.name,f) for f in (O/'inputs').iterdir()),*(('native/'+f.relative_to(LOCAL/'native').as_posix(),f) for f in (LOCAL/'native').rglob('*') if f.is_file()),*(('logs/geometry/'+f.name,f) for f in LOG.iterdir()),*(('logs/checks/'+Path(a['log']).name,R/a['log']) for a in read(cdir/'attempts.json'))]
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for name,f in originals:z.writestr(name,f.read_bytes());members[name]=dict(original_path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=sha(f))
    assert archive.stat().st_size<5_000_000
    restored=LOCAL/'restored';assert not restored.exists();restored.mkdir()
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==set(members)
        for name,item in members.items():
            f=restored/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_bytes(z.read(name));assert sha(f)==item['sha256']
    save(O/'download_manifest.json',dict(archive='downloads/raw.zip',archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),members=members))
    save(O/'archive_verified.json',dict(files=len(members),raw_bytes=sum(v['bytes'] for v in members.values()),all_restored_byte_exact=True,manifest_sha256=sha(O/'download_manifest.json'),restore_directory=str(restored.relative_to(R))))
    plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    y=np.arange(7);orig=[];cand=[]
    for case in summary['cases']:
        a,b=case['observations'];v=a['normal_closure']
        if a['completed']==0:v=float(re.search(r'normal_closure=([^ ]+)',a['message']).group(1))
        orig.append(v);cand.append(b['normal_closure'])
    fig,ax=plt.subplots(figsize=(9.5,5.2),layout='constrained')
    ax.semilogx(np.maximum(orig,1e-18),y-.1,'o',color=BLUE,label='Original (short-face case rejected)');ax.semilogx(np.maximum(cand,1e-18),y+.1,'s',color=RED,label='Retained-face candidate')
    ax.axvline(1e-9,color='#555f69',ls='--',label='Unchanged original closure gate');ax.set_yticks(y,['9 sites: unperturbed','9 sites: short face','9 sites: resolved face','Actual dx20: 10 us','Actual dx10: 10 us','Actual dx5: 10 us','Failed-case saved dx5: 1.25 us']);ax.invert_yaxis();ax.set(xlabel='Normalized native normal closure',title='Same gate; geometry-only candidate observations');ax.legend(frameon=False,fontsize=9,loc='lower right');ax.grid(axis='x',alpha=.15)
    for ext in ['png','pdf']:fig.savefig(O/('closure-comparison.'+ext),dpi=200)
    plt.close(fig)
    css='*{box-sizing:border-box}body{margin:0;background:#fff;color:#213345;font:17px/1.9 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1060px;margin:auto;padding:30px 24px}h1{font-size:clamp(30px,5vw,46px);line-height:1.4}h2{font-size:26px}a{color:#3c5488}section{padding:28px 0;border-bottom:1px solid #d9dfe2}img{width:100%;height:auto}figure{margin:24px 0}figcaption{font-size:14px;color:#687584}.note{padding:20px;background:#faf3f1;border-left:3px solid #e64b35}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}.card{padding:20px;background:#f4f7fa;border-top:3px solid #3c5488}.card b{display:block;font-size:28px}code{overflow-wrap:anywhere}@media(max-width:650px){.cards{grid-template-columns:1fr}}'
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VariDrainage · 保留短边候选</title><style>{css}</style></head><body><main><a href="https://baiming-zhang.github.io/VariDrainage/">项目首页</a><section><h1>保留短边的候选，<br>通过几何与既有检查。</h1><p>候选只在原坐标去重之后，保留所有正长度边。面积、正交、面配对和闭合检查照常执行；原默认求解器与已归档二进制均未修改。</p><div class="note"><b>尚未验收完整动力学。</b>这不是实际第78步重放，也没有证明那次失败的确切根因。旧时间配对仍不可用，未启动新的长作业。</div></section><section><h2>固定反例现在怎样了？</h2><p>原函数在9点短边反例上被10⁻⁹闭合门槛拒绝；候选在同一输入、同一门槛下得到闭合残差 <b>5.876×10⁻¹⁷</b>。未扰动、较大扰动两个控制组也通过。原失败记录完整保留。</p><figure><img src="closure-comparison.png" alt="七个固定输入中原函数与保留短边候选的原生几何闭合残差"><figcaption>两种实现同进程比较；对数图以10⁻¹⁸显示零值，未修改原数据。曲线点表示原生面组装结果，不是多边形恒等式。<a href="closure-comparison.pdf">矢量PDF</a></figcaption></figure><p>另取四份实际保存几何：20／10微米网格的10微秒、5微米首组10微秒，以及失败组最后保存的1.25微秒。两种实现的单元顶点、面积、邻接和导通系数在这些快照上逐项相同；并未重跑这些轨迹。</p></section><section><h2>原有曲面与机械功检查仍成立吗？</h2><div class="cards"><div class="card"><b>14</b>曲面能量等既有检查通过</div><div class="card"><b>4</b>初态准备检查通过</div><div class="card"><b>6</b>机械功观察检查通过</div></div><p>共24项是本轮新执行，使用原测试内容、容差和CHOLMOD/H1配置，通过编译选择候选几何函数。检查中包含原有小规模耦合步；没有重放冻结的长期算例。通过这些固定工况，不代表所有拓扑变化附近的导数或实际时间／空间收敛都已成立。</p></section><section><h2>下一道关卡：真实算例的有界检查</h2><p>先检验实际算例的短步原生账本、机械功和结果对照，再按证据决定是否重新计算同一候选实现下的时间配对。不能直接用旧实现的粗时间步和新实现的细时间步拼成合格配对。任何超过2小时的后续计算仍须每半小时监视。</p><p>物理通量／牵引力、原10／50／200微秒、交换系数、薄膜元厚度比和整体实验仍待完成，C1长期目标继续。</p><p><a href="https://baiming-zhang.github.io/VariDrainage/short-face/">原9点反例</a> · <a href="https://baiming-zhang.github.io/VariDrainage/uniform5-second/">原第78步失败</a>。</p></section><section><details><summary>全部固定输入、原始输出和代码</summary><p><a href="summary.json">七组几何对照</a> · <a href="checks/summary.json">24项新执行检查</a> · <a href="protocol.json">几何观测前协议</a> · <a href="checks/protocol.json">检查前协议</a> · <a href="attempts.json">几何实际命令</a> · <a href="checks/attempts.json">检查实际命令</a>。</p><p><a href="downloads/raw.zip">完整原始输入／输出与关闭日志ZIP</a> · <a href="download_manifest.json">字节清单</a> · <a href="archive_verified.json">实际恢复核验</a>。{len(members)}份原始文件已恢复且逐文件字节一致。<a href="source_index.json">源码SHA与基点</a> · <a href="source-at-run/thin_film/scripts/retained_faces/PLAN.md">固定几何计划</a> · <a href="source-at-run/thin_film/scripts/retained_checks/PLAN.md">既有检查计划</a> · <a href="source-at-run/thin_film/include/vf/thinfilm/VoronoiRetainedFacesCandidate.h">候选源码</a> · <a href="source-at-run/VoronoiDrainage.original.h">原函数源码</a>。ZIP保留原始字节，处理后文本可能经Git换行规范化。</p></details></section></main></body></html>'''
    (O/'index.html').write_text(page,encoding='utf-8')
    save(O/'report_source.json',dict(source_sha256=sha(Path(__file__)),geometry_summary_sha256=sha(O/'summary.json'),checks_summary_sha256=sha(cdir/'summary.json'),full_C1_complete=False,scope='Candidate geometry and existing-fixture admission only; original long-run failure preserved'))
    print('Reported seven paired geometry inputs and24new fixture checks; no actual78 replay or full dynamics qualification.')
if __name__=='__main__':main()
