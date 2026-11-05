# 被动输运质量账本与近纯残留

12 项诊断场与 54 次正式计时已完成。完整账本来自诊断运行，不拿开启逐扫掠日志的耗时作为性能结论。面平流、压缩、端点清理、扫掠间变化和更新舍入分别保存。余额闭合只证明账目核对，不能代替实际面积漂移指标。

| 算例 | 面积相对漂移 | 平流变化 | 压缩变化 | 清理变化 | 扫掠间变化 | 余额 |
|---|---:|---:|---:|---:|---:|---:|
| translation_32_geometric_ledger | -3.267e-17 | -4.938e-20 | 0 | -4.662e-21 | 0 | 0 |
| translation_32_neural_cuda_ledger | -3.267e-17 | 7.902e-20 | 0 | 7.264e-18 | -3.388e-21 | 0 |
| translation_64_geometric_ledger | -9.802e-17 | 4.124e-21 | 0 | -1.095e-19 | -3.388e-21 | 0 |
| translation_64_neural_cuda_ledger | -1.634e-16 | -1.367e-20 | 0 | 0 | 3.388e-21 | 0 |
| translation_128_geometric_ledger | 3.267e-17 | 1.89e-20 | 0 | -1.298e-19 | 1.016e-20 | 0 |
| translation_128_neural_cuda_ledger | 1.023e-13 | 1.262e-20 | 0 | 5.443e-15 | 0 | 0 |
| vortex_return_32_geometric_ledger | -2.475e-12 | -7.58e-19 | 1.7e-18 | -1.749e-13 | 0 | 3.123e-23 |
| vortex_return_32_neural_cuda_ledger | -2.331e-10 | 4.354e-19 | -9.963e-19 | -1.648e-11 | 8.132e-20 | 6.872e-22 |
| vortex_return_64_geometric_ledger | -6.837e-13 | 1.444e-19 | -1.344e-18 | -4.832e-14 | 0 | 1.113e-22 |
| vortex_return_64_neural_cuda_ledger | -1.706e-10 | -3.653e-20 | -1.729e-18 | -1.206e-11 | -5.421e-20 | -1.104e-22 |
| vortex_return_128_geometric_ledger | 2.778e-13 | -6.976e-20 | -4.374e-19 | 1.964e-14 | 6.776e-21 | 3.118e-23 |
| vortex_return_128_neural_cuda_ledger | -3.111e-10 | -2.269e-19 | 7.971e-21 | -2.199e-11 | -1.016e-19 | 1.329e-22 |

![不同门限的残留格](figures/mixed_tail.png)

门限 ε 决定统计口径，图中每一条曲线对应同一份终态。近纯门控阈值 0.002 与端点 cleanup 阈值 1e-10 不相同。神经输出变小、被门控跳过、实际清理为零三个概念应分开。

耦合动力学账本另见 ../20261104_coupled/metrics.json；本报告补齐被动平移与反转涡旋分支。实际性能与精确解误差见 ../../E07/20261104_passive/README.md。
