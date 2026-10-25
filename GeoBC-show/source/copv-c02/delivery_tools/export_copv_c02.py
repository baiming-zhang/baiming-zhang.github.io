"""Publish completed C02 evidence; never promote an unconverged run."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import ast, hashlib, html, json, re, shutil, zipfile

WS = Path(__file__).resolve().parents[1]
REPO = WS / 'publication/GeoBC-Neural-Operator'
SITE = WS / 'showcase-public-site/GeoBC-show'
ROOT = WS / 'data/copv/c02_structural_v1'

def read(p): return json.loads(p.read_text(encoding='utf-8'))
def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)

def source_dependencies(files):
    """Include local Python imports, without shipping unrelated project files."""
    pending = list(files); seen = set()
    while pending:
        relative = pending.pop()
        if relative in seen: continue
        path = REPO / relative
        if not path.is_file(): raise FileNotFoundError(path)
        seen.add(relative)
        if path.suffix != '.py': continue
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8-sig'))):
            modules = [node.module] if isinstance(node, ast.ImportFrom) else [x.name for x in node.names] if isinstance(node, ast.Import) else []
            if isinstance(node, ast.ImportFrom) and node.level:
                package = list(Path(relative).parent.parts)
                base = package[:len(package)-node.level+1]
                modules = ['.'.join(base + (node.module.split('.') if node.module else [x.name])) for x in node.names]
            for module in modules:
                if module and module.startswith(('geobc.', 'tools.')):
                    candidate = module.replace('.', '/') + '.py'
                    if (REPO / candidate).is_file(): pending.append(candidate)
    return sorted(seen)

def main():
    s = read(ROOT / 'fenicsx/summary.json')
    ref = read(ROOT / 'reference/comparison.json')
    patch_gpu = read(ROOT / 'patch/fenicsx_review.json')
    patch_ccx = read(ROOT / 'patch/calculix_review.json')
    closure = read(ROOT / 'surface/load_closure_review.json')
    assert s['status'] == 'completed' and s['peak_reason'] > 0
    assert s['true_relative_residual'] < 1e-7 and s['energy_work_relative_error'] < 1e-7
    assert s['half_load_relative_error'] < 1e-6 and s['zero_load_U_norm'] < 1e-12
    assert ref['passed'] and patch_gpu['passed'] and patch_ccx['passed'] and closure['passed']
    contract_path = REPO / 'configs/copv_c02_structural_v1.json'
    solver_path = REPO / 'tools/copv_fenicsx_structural.py'
    assert s['contract_sha256'] == digest(contract_path)
    assert s['solver_sha256'] == digest(ROOT / 'fenicsx/solver_executed_snapshot.py')
    assert s['mesh_sha256'] == digest(ROOT / 'mesh_h6_netgen/mesh.npz')
    asset = SITE / 'assets/copv-c02'
    for name in ['bc-loads.json','bc-load-surfaces.json.gz']:
        copy(ROOT / 'visualization' / name, asset / name)
    for name in ['manifest.json', 'surface.json.gz', 'section.json.gz']:
        assert (asset / name).is_file(), name
    for name in ['reference-01-U.jpg','reference-02-fiber-section.jpg','reference-03-metal-section.jpg','reference-04-bc-loads.jpg']:
        if (ROOT / name).is_file(): copy(ROOT / name, asset / name)

    evidence = {'schema': 'copv-c02-evidence-v1', 'detailed_3D': s,
        'independent_axisymmetric_specimen': ref,
        'analytical_rotated_TET10_patch_FEniCSx': patch_gpu,
        'analytical_rotated_TET10_patch_CalculiX': patch_ccx,
        'independent_pressure_and_caps_closure': closure,
        'comparison_scope': 'CalculiX checks a rotated straight TET10 patch; skfem checks a geometry-matched axisymmetric vessel. Neither is a second full detailed-assembly solve.',
        'stress_scope': 'Raw original-cell GP tensors retained locally; public display uses per-cell GP volume means without cross-material averaging. Composite von Mises is not a failure criterion.'}
    write(SITE / 'data/copv-c02-review.json', evidence)
    copy(contract_path, SITE / 'data/copv-c02-material-load-contract.json')
    copy(ROOT / 'mesh_h6_netgen/mesh_review.json', SITE / 'data/copv-c02-mesh-review.json')
    copy(ROOT / 'fenicsx/load_unload_states.json', SITE / 'data/copv-c02-load-states.json')

    bundle = SITE / 'source/copv-c02'
    seed = ['tools/copv_fenicsx_structural.py', 'tools/prepare_copv_c02_mesh.py',
        'tools/mesh_copv_c02_partition.py', 'tools/prepare_copv_c02_surface.py',
        'tools/prepare_copv_c02_reference.py', 'tools/copv_fenicsx_axisymmetric.py',
        'tools/prepare_copv_c02_patch.py', 'tools/copv_fenicsx_patch.py',
        'tools/export_copv_c02_viewer.py', 'tools/export_copv_c02_bc.py', 'configs/copv_c02_structural_v1.json',
        'configs/copv_material_candidates.json']
    for name in source_dependencies(seed): copy(REPO / name, bundle / name)
    for name in ['geobc', 'tools']: (bundle / name / '__init__.py').write_text('', encoding='utf-8')
    for name in ['mesh/structural_partition.brep', 'mesh/partition_review.json',
        'mesh_h6_netgen/mesh_review.json', 'surface/load_closure_review.json',
        'reference/comparison.json', 'reference/fenicsx/summary.json',
        'reference/skfem/solution.json', 'patch/input.npz', 'patch/patch.inp',
        'patch/patch.dat', 'patch/fenicsx_review.json', 'patch/calculix_review.json',
        'fenicsx/summary.json', 'fenicsx/peak_algebra.json', 'fenicsx/linearity.json',
        'fenicsx/load_unload_states.json', 'fenicsx/ksp_view.txt', 'fenicsx/petsc_options.json', 'fenicsx/residual_refinement.json', 'fenicsx/solver_executed_snapshot.py', 'visualization/bc-loads.json']:
        copy(ROOT / name, bundle / 'evidence' / name)
    for name in ['compare_copv_c02_reference.py', 'export_copv_c02.py']:
        copy(WS / 'notes' / name, bundle / 'delivery_tools' / name)

    readme = f'''# C02 detailed COPV elastic benchmark

Accepted C01 construction: 12 structural solids (Al liner, epoxy bondline, four orthotropic ply groups, two polar collars, two Ti bosses, two steel inserts). Four O-rings and two PEEK rings remain display-only. Shared structural interfaces are perfectly bonded. No contact opening, seal compression, plasticity, damage, burst or rated-pressure qualification is calculated.

Units: mm, N, MPa, N mm. Pressure is 1 MPa on the wetted inner wall and port shoulders. Separate closed-end resultants p*pi*r^2 are applied once to the two annular ends: 326.8513 N and 213.8246 N. A rank-six 3-2-1 nodal gauge removes free rigid motion; it is not a mounting clamp. The independent discrete-surface closure check records the small force remainder caused by analytical circular cap areas versus polygonized bores.

The accepted P1 solve has {s['nodes']:,} corner nodes and {s['dofs']:,} displacement DOF on {s['cells']:,} tetrahedra. The coordinate/output mesh stores {s['geometry_nodes']:,} TET10 points. Midside displacements are exact affine interpolation, not additional P2 unknowns. Geometry is straight. Full detailed P2 AMG attempts exceeded the current resource budget. The two sampled polar collars required an inner-interface fit to the exact outer-ply surface; source hashes and volume changes are recorded in partition_review.json. Mesh, P1/P2 and local peak-stress convergence belong to C03.

Material constants and assumptions are explicit in configs/copv_c02_structural_v1.json. The carbon lamina uses all nine 3-D constants, engineering shear and GP-dependent [fiber, transverse, normal] axes. Ply angles are [55,-55,-55,55] degrees; polar collars use +25/-25 degrees. These are research direction fields, not qualified winding paths. Ti E is inherited from R02, while nu=0.28 is the retained modeling assumption. Steel E=200 GPa is sourced, nu=0.30 is assumed. Adhesive E=3 GPa and nu=0.35 are both numerical assumptions.

Backend: DOLFINx {s['dolfinx']}, PETSc {'.'.join(map(str,s['petsc']))}. {s['backend']}. Matrix/vector types: {s['matrix_type']} / {s['vector_type']}. Exact hierarchy, coarse solve and options are retained in evidence/fenicsx/ksp_view.txt and petsc_options.json. Failed solver attempts are not accepted results.

Validation: independent skfem versus GPU FEniCSx axisymmetric vessel U relative L2={ref['U_relative_L2']:.8e}, volume-weighted raw GP tensor difference={ref['raw_GP_stress_volume_weighted_relative_Frobenius']:.8e}. A rotated straight-sided TET10 patch matches analytical fields in both CalculiX and GPU FEniCSx; CalculiX stress output is ASCII precision. These P2 verification specimens are distinct from the full detailed P1 solve. Neither check is a second full detailed-assembly CalculiX calculation.

Detailed solve: true residual={s['true_relative_residual']:.8e}, energy/work error={s['energy_work_relative_error']:.8e}, half-load independent solve error={s['half_load_relative_error']:.8e}. Full peak solve {s['peak_solve_seconds']:.2f} s; accepted-run total {s['total_seconds']:.2f} s (excludes earlier rejected CAD/mesh/solver attempts and independent reference preparations). Twenty detailed display states are exact linear superpositions, independently checked at half and zero load. The smaller axisymmetric reference performed twenty actual solves. No nonlinear unloading or hysteresis is claimed.

Raw original-element GP local/global stress, coordinates and weights are saved locally as float64 .npy files. Public stress colors use per-cell volume-weighted GP means, retain material/ply discontinuities and are not raw peak values. P1 displacement is evaluated for the FE cut section. No temporal frame blending is used. Von Mises is displayed only for metal parts; fiber stress is displayed for composites. Composite strength allowables/failure indices are not supplied.

Reproduce in a DOLFINx/PETSc CUDA environment with CPU HYPRE, NumPy, SciPy, Gmsh 4.15.2, scikit-fem and optional CalculiX 2.23:
```
python -m tools.mesh_copv_c02_partition --partition evidence/mesh/structural_partition.brep --review evidence/mesh/partition_review.json --output mesh_h6_netgen --h 6 --curvature 64 --algorithm 1 --straight-geometry --optimize-netgen
python -m tools.copv_fenicsx_structural --mesh mesh_h6_netgen/mesh.npz --review mesh_h6_netgen/mesh_review.json --contract configs/copv_c02_structural_v1.json --output fenicsx
```
Set PETSC_OPTIONS='-use_gpu_aware_mpi 0', OPENBLAS_NUM_THREADS=4 and OMP_NUM_THREADS=4 for this single-rank run. The download includes the repaired partition, source and small verification evidence; the volume mesh and large raw GP arrays are generated locally, not embedded in the website download. The exact executed solver is retained in evidence/fenicsx/solver_executed_snapshot.py and matches solver_sha256. The delivered tools/copv_fenicsx_structural.py only corrects a stale displacement-description string in report metadata; numerical code and saved fields are unchanged. CLI help on the reference and patch scripts gives their independent reproduction steps.

The live viewer has independent boundary-condition and load visibility buttons. A/B/C mark the actual saved lower-port nodes: A locks Ux/Uy/Uz, B locks Uy/Uz, C locks Uz in global coordinates. This is a numerical 3-2-1 gauge, not a fixed end face, wall or hinge. Orange highlights all 57,544 selected wetted facets; direction arrows are a 28-facet sample. Magenta highlights the actual annular end faces and displays the two equivalent cap resultants. Markers are shown through the model and follow the displayed displacement; arrow lengths are schematic. Loads and gauge reactions scale with the selected elastic state, and zero load has no force arrows. Visibility buttons do not change the solved boundary conditions. No FEM calculation was rerun for this display addition.

Re-export these locations after the mesh/solution commands above:
```
python -m tools.prepare_copv_c02_surface --mesh mesh_h6_netgen/mesh.npz --output surface
python -m tools.export_copv_c02_bc --root . --output visualization
```
The export verifies the mesh hash, six saved constraint components, selected facet counts and the integrated pressure/cap resultant against the accepted C02 evidence.
'''
    (bundle / 'README.md').write_text(readme, encoding='utf-8')
    hashes = {str(p.relative_to(bundle)).replace('\\','/'): digest(p) for p in bundle.rglob('*') if p.is_file() and p.name != 'SHA256.json'}
    write(bundle / 'SHA256.json', hashes)
    archive = SITE / 'data/copv-c02-source-evidence.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in bundle.rglob('*'):
            if p.is_file(): z.write(p, p.relative_to(bundle))

    summary = f"C02 已完成：详细装配的 12 个承力分区，P1 位移、{s['dofs']:,} 自由度；1 MPa 内压与闭口端合力，峰值位移 {s['max_U_mm']:.4f} mm。支持位移/应力剖面，以及边界条件和载荷的独立显示开关；源码与核对记录可复现。"
    review = f"参考 R02 的材料/载荷合同及原始场核对方法：独立轴对称模型位移差 {ref['U_relative_L2']:.3e}、原始 GP 应力差 {ref['raw_GP_stress_volume_weighted_relative_Frobenius']:.3e}；另有 CalculiX 三维旋转铺层补丁核对。详细模型真残差 {s['true_relative_residual']:.3e}。CPU 组装，CUDA FGMRES/GAMG，最粗层为小型 CPU SVD；理想粘结、直边几何及 P1 弹性基准的范围见详情，网格及 P1/P2 收敛留给 C03。"
    links = [{'label': a, 'url': b} for a,b in [
        ('三维结果与剖面', 'viewers/copv-c02.html'), ('核对指标与计算范围', 'data/copv-c02-review.json'),
        ('材料 / 载荷 / 约束合同', 'data/copv-c02-material-load-contract.json'),
        ('求解与验证源码 ZIP', 'data/copv-c02-source-evidence.zip'),
        ('复现与结果解释', 'source/copv-c02/README.md'), ('20 个弹性载卸载状态', 'data/copv-c02-load-states.json')]]
    stamp = datetime.now(ZoneInfo('Asia/Shanghai'))
    tasks = read(SITE / 'data/tasks.json')
    t = next(t for t in tasks['tasks'] if t['id'] == 'C02')
    t.update(status='已完成', issues='弹性基准、材料方向、载荷闭合及求解核对通过；密封接触/损伤/额定压力不属于本阶段，网格收敛由 C03 完成。', results=links, actual_start='2026-10-25', completed=stamp.strftime('%Y-%m-%d'))
    c01 = next(t for t in tasks['tasks'] if t['id'] == 'C01')
    c01['issues'] = '用户已确认 C01 完成。详细构造的承力分区、材料及求解结果见 C02；原 v1 参数样例与接入验证保留。'
    tasks['updated'] = stamp.strftime('%Y-%m-%d %H:%M 北京时间'); write(SITE / 'data/tasks.json', tasks)
    notes = read(SITE / 'data/task-notes.json')
    notes['notes']['C02'].update(summary=summary, review=review)
    notes['notes']['C01']['review'] = '用户已确认 C01 构造完成。纹理与螺纹示意的范围保留；详细构造承力分区的正式弹性计算及材料/载荷合同见 C02。'
    write(SITE / 'data/task-notes.json', notes)
    updates = read(SITE / 'data/updates.json')
    title = 'C02：详细 COPV 三维正交各向异性 GPU 求解与独立核对'
    if not any(u['task']=='C02' and u['title']==title for u in updates['updates']):
        updates['updates'].append(dict(task='C02', date=stamp.strftime('%Y-%m-%d'), title=title, summary=summary, review=review, links=links))
    bc_title='C02：新增实际边界条件与载荷显示开关'
    if not any(u['task']=='C02' and u['title']==bc_title for u in updates['updates']):
        updates['updates'].append(dict(task='C02',date=stamp.strftime('%Y-%m-%d'),title=bc_title,summary='三点 3–2–1 约束、完整受压面、压力方向采样与上下端口等效合力可独立显示，随状态更新数值。',review='位置取自已执行网格和保存的实际约束 DOF；墙面/铰链未设置。显示开关不改变求解工况，未重跑 FE。',links=[dict(label='打开三维结果',url='viewers/copv-c02.html'),dict(label='边界与载荷位置数据',url='assets/copv-c02/bc-loads.json')]))
    write(SITE / 'data/updates.json', updates)

    body = f'''<section id="copv-c02-method"><h2>边界条件 / 载荷可视化</h2><p>在查看器打开“边界条件”可显示下端 A/B/C 三点及各自锁定的全局位移分量；打开“载荷”可透视高亮真实受压面与环形端面，并显示压力方向及上下端合力。可配合剖面查看。箭头长度仅示意方向，载荷数值随状态变化；显示开关不改变求解工况。本例没有墙面或铰链支承。<a href="assets/copv-c02/bc-loads.json">实际节点、分量与载荷数据</a> · <a href="assets/copv-c02/reference-04-bc-loads.jpg">BC 与载荷参考图</a>。</p><h2>材料、载荷与求解范围</h2>
<p>计算 12 个承力实体：Al 6061-T6 内衬、环氧粘接层、四个碳纤维/环氧铺层分组、两个极区补强、Ti-6Al-4V boss 与 316L 接口套。六个密封/挡圈实体只在 C01 构造中展示。承力界面为共享网格的理想粘结。</p>
<p>单位 mm–N–MPa；内压 1 MPa 是代码基准。内湿壁与端口内肩承受压力，上下端口另各施加一次 pπr² 闭口端合力，分别 326.8513 N 和 213.8246 N。三点 3–2–1 约束只去除六个刚体模态。解析圆孔端面积与直边网格之间保留 0.181476 N 的轴向余量，约为 pπR² 的 5.777×10⁻⁶，由数值定位约束反力平衡；这不是实际安装反力。离散端面积闭合属于 C03 的检查项。</p>
<p>AS4/3501-6 使用完整九个三维弹性常数，按原始积分点的局部纤维/横向/法向基旋转；铺层角从内到外为 +55/−55/−55/+55°，上下极区为 +25/−25°。这是研究用方向场。碳材料与 Al 弹性参数引用 NASA 文献；Ti 的 E 沿用 R02，ν=0.28 为建模假设；钢 E=200 GPa 引用材料表，ν=0.30 为假设；胶层 E=3 GPa、ν=0.35 均为数值假设。完整数值与来源见材料合同。</p>
<table><thead><tr><th>材料</th><th>E [MPa]</th><th>ν</th><th>数值来源 / 假设</th></tr></thead><tbody><tr><td>Al 6061-T6</td><td>69000</td><td>0.33</td><td><a href="https://ntrs.nasa.gov/citations/20180003005">NASA 工程分析表 4</a></td></tr><tr><td>Ti-6Al-4V</td><td>115100</td><td>0.28</td><td><a href="data/material-reference.json">R02 来源记录</a>；ν 为沿用假设</td></tr><tr><td>316L</td><td>200000</td><td>0.30</td><td><a href="https://www.alfalaval.no/produkter/vaskebehandling/en-1-4404-aisi-316l/">材料弹性模量表</a>；ν 为假设</td></tr><tr><td>环氧粘接层</td><td>3000</td><td>0.35</td><td>两者均为数值基准假设</td></tr></tbody></table><p>AS4/3501-6：E₁/E₂/E₃ = 114452.971 / 10342.136 / 10342.136 MPa；G₁₂/G₁₃/G₂₃ = 5998.439 / 5998.439 / 3792.117 MPa；ν₁₂/ν₁₃/ν₂₃ = 0.33 / 0.33 / 0.47。<a href="https://ntrs.nasa.gov/citations/19940028433">NASA-TM-109108 表 1</a>转引上游参考文献，本次未独立取得该上游报告；来源与转换精度保存在合同中。</p>
<p>{s['nodes']:,} 个角点、{s['cells']:,} 个四面体、{s['dofs']:,} 自由度。位移为 P1，几何为直边；坐标/输出网格保留 {s['geometry_nodes']:,} 个 TET10 点，中边点位移按仿射关系插值。完整 P2 装配的预条件层超出资源预算，本次不声称详细 P2 已通过。两处极区补强的内表面拟合修复及体积变化有单独记录；线性体网格经 Netgen 优化，修正了近退化单元。CPU FFCx 组装，CUDA FGMRES/GAMG 矩阵和向量运算，最粗层为小型 CPU SVD，完整层级见复现包。</p>
<h2>参考 R02 的核对方法</h2><table><thead><tr><th>对象</th><th>独立核对</th><th>结果</th></tr></thead><tbody>
<tr><td>旋转正交各向异性 TET10 补丁</td><td>解析场 / GPU FEniCSx / CalculiX</td><td>GPU 应力误差 {patch_gpu['stress_relative_L2']:.3e}；Cal 应力误差 {patch_ccx['stress_relative_L2']:.3e}（ASCII 输出）</td></tr>
<tr><td>同几何轴对称闭口多层容器</td><td>独立 skfem / GPU FEniCSx</td><td>U 差 {ref['U_relative_L2']:.3e}；原始 GP 张量差 {ref['raw_GP_stress_volume_weighted_relative_Frobenius']:.3e}</td></tr>
<tr><td>完整详细三维装配</td><td>真残差 / 力矩平衡 / 功能关系 / 独立半载及零载</td><td>残差 {s['true_relative_residual']:.3e}；功能误差 {s['energy_work_relative_error']:.3e}；半载位移误差 {s['half_load_relative_error']:.3e}</td></tr>
</tbody></table><p>上述 CalculiX 与 skfem 核对属于两个明确的验证试件，没有声称完成另一套详细装配全规模 CalculiX 计算。详细三维模型峰值解耗时 {s['peak_solve_seconds']:.1f} 秒，本次通过的运行共 {s['total_seconds']:.1f} 秒；不包含此前被拒绝的构造、网格和预条件器尝试。</p>
<h2>场与载卸载状态如何解释</h2><p>原始单元/积分点的全局及局部应力、位置和积分权重以 float64 保存在本地。网页应力色为各单元原始 GP 的体积加权均值，不跨材料或铺层平均；切面位移在 FE 域中按 P1 插值。金属显示 von Mises，复合层显示纤维方向应力。</p><p>详细模型的 20 个状态按线弹性解精确比例生成，另从零初值独立求解半载与零载验证；小型轴对称参考模型实际执行了 20 次求解。这些结果不表示塑性回线、密封压缩、接触开闭、损伤或爆破强度。局部应力峰值及网格/几何/P1-P2 收敛由 C03 处理。</p><p>少量参考视图：<a href="assets/copv-c02/reference-01-U.jpg">位移全景</a> · <a href="assets/copv-c02/reference-02-fiber-section.jpg">纤维应力剖面</a> · <a href="assets/copv-c02/reference-03-metal-section.jpg">金属应力剖面</a>。</p></section>'''
    page = (SITE / 'C02.html').read_text(encoding='utf-8')
    stage = '<section class="delivery-stage" aria-label="可视化交付物"><iframe src="viewers/copv-c02.html?v=c02-bc-1" title="COPV 三维位移、边界条件与载荷" loading="eager" style="width:100%;border:0"></iframe></section>'
    page = re.sub(r'<section class="delivery-stage".*?</section>', lambda _:stage, page, count=1, flags=re.S)
    for key,value in [('public-summary',summary),('public-review',review),('task-issues',t['issues'])]:
        page = re.sub(f'(<p id="{key}">).*?(</p>)', lambda m:m[1]+html.escape(value)+m[2], page, count=1)
    page = page.replace('data-meta="status">未开始', 'data-meta="status">已完成')
    page = re.sub(r'(<div class="link-row" id="task-results">).*?(</div>)',lambda m:m[1]+''.join(f'<a href="{l["url"]}">{l["label"]} ↗</a>' for l in links)+m[2],page,count=1)
    if 'id="copv-c02-method"' in page: page = re.sub(r'<section id="copv-c02-method">.*?</section>',lambda _:body,page,count=1,flags=re.S)
    else: page = page.replace('<section id="daily">',body+'<section id="daily">')
    (SITE / 'C02.html').write_text(page,encoding='utf-8')
    c01page = (SITE / 'C01.html').read_text(encoding='utf-8')
    c01review = notes['notes']['C01']['review']
    c01page = re.sub(r'(<p id="public-review">).*?(</p>)',lambda m:m[1]+c01review+m[2],c01page,count=1)
    c01page = re.sub(r'(<p id="task-issues">).*?(</p>)',lambda m:m[1]+c01['issues']+m[2],c01page,count=1)
    c01page = c01page.replace('未对这套 boss / 密封 / 铺层装配计算新 FE 结果。正式材料、铺层方向、接触、端部合力和收敛仍由 C02/C03 验收。','详细装配的 12 个承力分区现已完成 C02 弹性计算，见 <a href="C02.html">材料、载荷与三维结果</a>；密封接触及局部峰值收敛范围按 C02/C03 说明解释。')
    (SITE / 'C01.html').write_text(c01page,encoding='utf-8')
    print(json.dumps(dict(published_evidence=True,source_files=len(hashes),zip_bytes=archive.stat().st_size,max_U_mm=s['max_U_mm']),ensure_ascii=False))

if __name__ == '__main__': main()
