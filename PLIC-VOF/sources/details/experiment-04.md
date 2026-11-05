# 同一模拟器中的通量—曲率对照

冻结重训模型；共 42 项预先声明的测量。所有方法使用相同初始化、材料平滑、边界、时间步、压力后端与几何重建；统一关闭 sharpening。速度不在本轮并发测量中作结论。

| 算例 | 网格 | 神经通量 | 神经曲率 | 干预 | 主指标 | 最大速度 |
|---|---|---|---|---|---:|
| static | 64×64 | False | False | none | 压差相对误差 1.587% | 0.000134353 |
| static | 64×64 | False | True | none | 压差相对误差 5.121% | 0.000430261 |
| static | 64×64 | True | False | none | 压差相对误差 1.531% | 0.000140136 |
| static | 64×64 | True | True | none | 压差相对误差 5.237% | 0.000454092 |
| static | 128×128 | False | False | none | 压差相对误差 2.174% | 0.000361988 |
| static | 128×128 | False | True | none | 压差相对误差 6.370% | 0.000228033 |
| static | 128×128 | True | False | none | 压差相对误差 2.143% | 0.000344065 |
| static | 128×128 | True | True | none | 压差相对误差 5.096% | 0.000225548 |
| wave | 64×64 | False | False | none | 幅值 RMS 0.000716238 | 0.0140115 |
| wave | 64×64 | False | True | none | 幅值 RMS 0.000325811 | 0.0175461 |
| wave | 64×64 | True | False | none | 幅值 RMS 0.000755173 | 0.0150551 |
| wave | 64×64 | True | True | none | 幅值 RMS 0.000168852 | 0.0219071 |
| wave | 128×128 | False | False | none | 幅值 RMS 0.000415126 | 0.0176353 |
| wave | 128×128 | False | True | none | 幅值 RMS 0.000184847 | 0.0199278 |
| wave | 128×128 | True | False | none | 幅值 RMS 0.000386351 | 0.0198073 |
| wave | 128×128 | True | True | none | 幅值 RMS 0.000764151 | 0.078585 |
| bubble | 128×32 | False | False | none | 末帧质心 x=1.08191（数值观测） | 0.476931 |
| bubble | 128×32 | False | True | none | 末帧质心 x=1.08295（数值观测） | 0.47014 |
| bubble | 128×32 | True | False | none | 末帧质心 x=1.07571（数值观测） | 0.351793 |
| bubble | 128×32 | True | True | none | 末帧质心 x=1.1019（数值观测） | 0.476202 |
| bubble | 256×64 | False | False | none | 末帧质心 x=1.08271（数值观测） | 0.478259 |
| bubble | 256×64 | False | True | none | 末帧质心 x=1.08783（数值观测） | 0.475265 |
| bubble | 256×64 | True | False | none | 末帧质心 x=1.07341（数值观测） | 0.329088 |
| bubble | 256×64 | True | True | none | 末帧质心 x=1.10794（数值观测） | 0.513279 |
| static | 64×64 | False | False | constant_native_mask | 压差相对误差 0.000% | 3.47085e-13 |
| static | 64×64 | False | True | constant_native_mask | 压差相对误差 0.000% | 3.47085e-13 |
| static | 64×64 | True | False | constant_native_mask | 压差相对误差 0.000% | 3.47085e-13 |
| static | 64×64 | True | True | constant_native_mask | 压差相对误差 0.000% | 3.47085e-13 |
| wave | 64×64 | False | False | zero | 幅值 RMS 0.00647049 | 0 |
| wave | 64×64 | False | True | zero | 幅值 RMS 0.00647049 | 0 |

![静滴](figures/static_2x2.png)

![毛细波](figures/wave_2x2.png)

![气泡](figures/bubble_2x2.png)

固定场交叉保存每种输入上两种曲率的输出，初始干净圆有解析曲率；失真后的界面只比较方法差异，不再把理想圆当其真值。解析常曲率干预保留原掩码、哨兵和面力选择，区别于给解析路径额外外推。

波面直接记录每个时间步的逐列高度和复 Fourier 系数。幅值参考来自归档 Prosperetti 数值表，仅提供正幅值；已对单元列平均的 sinc 衰减作明确换算，不评分未知相位误差。非基频能量并非自动等同于数值噪声。

气泡使用冻结源码的旋转半域、Hysing 参数设置。表中质心和速度是本机计算结果，尚未引入独立高精度参考，不能称为精确解误差。

完整指标、所有退步、质量账本和原始场见 summary.json、metrics.csv 与每个配置目录。
