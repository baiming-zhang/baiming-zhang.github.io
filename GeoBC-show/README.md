# GeoBC-show — daily research directory

19 task pages, one stable URL per task. Published under the existing personal GitHub Pages site.

- `data/tasks.json`: canonical task status, dates, delivery and review criteria, result links.
- `data/updates.json`: append daily updates; never overwrite old review history.
- `assets/`: curated figures and compressed real mesh payloads.
- `R01.html` … `M03.html`: stable per-task pages.
- `guide.html`: daily update instructions.

The plan metadata uses the user-provided table. Existing scientific results are clearly separated from incomplete acceptance. Large raw solve files and model weights are not included. No external analytics or fonts.

## 任务页展示规范

- 首屏只放一个主要图片或交互模型，以及几句简短的验收／review评语。
- 禁止在任务页首屏再次嵌套该任务的完整页面或数字看板。
- 多 case 对比可以使用专用嵌入式查看器，但查看器不得复制任务页标题、评论区或详情结构；网格显示应作为独立开关。
- `error` 必须写明定义和基准。没有真值场时只显示求解残差／失败状态，不能把不同几何之间的场差冒充误差。
- 样本数、误差、耗时、参数表、质控矩阵、下载链接和完整记录统一放在页尾的“查看详情”折叠区内。
- 图片必须来自可追溯的真实结果；变形放大、应力平均等显示处理必须在图注中说明，不能暗示计算精度提高。
- 失败和隔离结果保留在详情中，不用只展示成功结果的方式掩盖范围限制。
