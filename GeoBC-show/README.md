# GeoBC-show — daily research directory

19 task pages, one stable URL per task. Published under the existing personal GitHub Pages site.

- `data/tasks.json`: canonical task status, dates, delivery and review criteria, result links.
- `data/updates.json`: append daily updates; never overwrite old review history.
- `assets/`: curated figures and compressed real mesh payloads.
- `R01.html` … `M03.html`: stable per-task pages.
- `guide.html`: daily update instructions.

The plan metadata uses the user-provided table. Existing scientific results are clearly separated from incomplete acceptance. Large raw solve files and resumable optimizer/RNG checkpoints remain local; selected inference weights are linked inside the task details. No external analytics or fonts.

## 任务页展示规范

- 首屏只放一个主要图片或交互模型，以及几句简短的验收／review评语。
- 禁止在任务页首屏再次嵌套该任务的完整页面或数字看板。
- 多 case 对比可以使用专用嵌入式查看器，但查看器不得复制任务页标题、评论区或详情结构；网格显示应作为独立开关。
- `error` 必须写明定义和基准。没有真值场时只显示求解残差／失败状态，不能把不同几何之间的场差冒充误差。
- 样本数、误差、耗时、参数表、质控矩阵、下载链接和完整记录统一放在页尾的“查看详情”折叠区内。
- 图片必须来自可追溯的真实结果；变形放大、应力平均等显示处理必须在图注中说明，不能暗示计算精度提高。
- 失败和隔离结果保留在详情中，不用只展示成功结果的方式掩盖范围限制。
- 所有展示应力或位移的三维视图都应提供实际受力箭头，参照 C04 的三维实体箭头格式。载荷独立开关，默认可见；压力为橙色，惯性体载荷为紫红色。方向来自原施加载荷，不能用位移方向代替力方向；箭头长度仅示意，位置随当前显示构形移动。

## R05–R08 展示维护

- 四页首屏使用 `viewers/rod-review-embed.html`：原节点表面与真实场的 WebGL 展示，不嵌套任务页或数值看板。
- R07/R08 真值与预测共享色阶；误差单独定义、单独色阶，在参考构形显示。默认展示该训练子集中联合相对误差最大的样本，允许切换工况和独立开启网格。
- 简短评语明确标注“AI审阅 · Codex”；完整审阅依据在 `assets/rod-review-visuals-v1/ai-review.json` 和页尾默认关闭的“查看详情”内。
- 数据导出源：研究代码仓库的 `tools/build_rod_review_visuals.py`；布局整理源：`tools/organize_rod_visual_pages.py`。训练发布器在有这些可视化资产时复用该布局，避免再次把曲线和指标放回正文。
- 载荷箭头由 `tools/export_rod_review_loads.py` 从冻结载荷／原节点与单元记录导出，使用 C04 相同的 `FESpatialGlyphs` 三维几何。真值、预测用相同外载，分别跟随各自保存位移；误差图在参考构形显示载荷。原图和求解数据不改写。

## R06 工况场景说明

- R06 选择器切换时，同步显示该工况的可能场景、受力原因和可能响应；仅在 R06 显示，不改 R07/R08 的训练场对比。
- 逐例来源为 `assets/r06-case-context-v1.json`，由研究仓库 `tools/export_rod_case_context.py` 从冻结接口和原载荷说明生成；记录原载荷、几何和说明文件哈希。
- 燃气净压与销孔承压须区分；纯轴向对照不因元数据存在 engine 字段就声称施加了惯性。厚度偏压的正负侧按原 κ 值说明。
- 可能响应是条件性机理解释，不是已观测损伤、失稳、疲劳或实测循环的结论。原标签、输入和训练结果不变。

## M01 连杆阶段汇总

- `assets/m01-checkpoint-v1/` 保存数据分布、代表场、模型误差、同例对比、1000条逐例CSV及四页离线PDF。网页正文只显示一张可切换的汇总图与简短AI审阅，明细和下载折叠到底部。
- 统计源为最终选定索引和权威accepted记录，图表源为large_v2已保存端点；`tools/build_m01_summary.py`提取统计，`tools/publish_m01_checkpoint.py`排版网页与PDF。没有重新求解FE或训练模型。
- `viewers/m01-render.html`和`assets/m01-render.js`只作真实原节点场的静态插图导出。代表图共用0–1000MPa节点VM色阶；模型同行共用数据范围，保留负预测，均按原位移×1并附实际载荷方向箭头。
- 代表图按物理机制和来源有目的选取；模型难例沿用TS3 TRAIN20既有联合误差最大规则，不按MGN表现挑选。所有原失败、检查点和源场继续保留。

- M01当前为“进行中”：已有报告只是数据与本地小规模拟合的阶段交付，V100上全部1000个场的完整训练尚待完成，不能因报告交付而结项。
