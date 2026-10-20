# RPI Wall Boiling：平面长方形通道测试

本算例参考 `papers/phaseChanges/energies-17-04225.pdf`，验证压力基
Eulerian–Eulerian 主循环中的 RPI 壁面沸腾、相间作用力、重力和壁面热源耦合。
流体工况沿用 DEBORA DEB7，但几何改为宽 `0.5 m`、长 `5 m` 的二维平面
长方形通道，用于延长 Eulerian–Eulerian/RPI 耦合回归。全部输入采用
Pa、m、kg、s、K；本算例不是圆管，也不使用轴对称体积权重。

## 几何与工况

- R12 饱和压力：`1.46e6 Pa`
- 质量通量：`2024 kg/(m2 s)`
- 壁面热流：`76260 W/m2`
- 入口液温：`317.36 K`
- 饱和温度：`331.2491346726 K`
- 平面通道宽度：`0.5 m`
- 总轴向长度：`5.0 m`
- 未加热入口段：`0 <= y < 1.0 m`
- 加热段：`1.0 <= y <= 4.5 m`
- 未加热出口段：`4.5 < y <= 5.0 m`
- 网格：`30 x 300 x 1` 横向—流向单层结构网格
- 均匀网格尺寸：`dx=dy=1.6667e-2 m`
- 几何：`geometryModel cartesian`
- 左壁：无滑移、绝热；右壁：无滑移，仅 `1.0 <= y <= 4.5 m` 加热

网格采用横向:流向 `1:10` 的单元数比例，物理尺寸同样为 `1:10`，因此单元近似
正方形。三个流向区段共用一个 `heatedWall` patch，
`constant/phaseChange` 的 `coordinate y` 与 `range (1.0 4.5)` 只把 RPI 热源
映射到中间加热段；入口和出口壁面保持绝热。

## 物性与闭式模型

`constant/phaseProperties` 分开定义 liquid、vapor 物性和相间模型：

```text
phaseProperties
├── phases (liquid vapor)
├── liquid/vapor
│   ├── rho(T), mu(T), Cp(T), k(T)
│   ├── 有效温度范围
│   └── 气相 constant diameterModel
└── interphaseModels.phasePairs.liquid_vapor
    ├── drag/lift/virtualMass
    ├── wallLubrication/turbulentDispersion
    ├── heatTransfer
    └── surfaceTension
```

R12 温变物性由 CoolProp 8.0.0 的 HEOS 后端在 `1.46 MPa` 下取样后拟合；液相
适用 `300–340 K`，气相适用 `315–400 K`。多项式超出显式有效范围或给出非正
物性时直接终止，不外推、不截断。`Cp(T)` 同时用于积分显焓，求解后的
`phaseEnthalpy` 通过数值反演恢复温度。

RPI 使用论文中的 Hibiki–Ishii 成核密度、Unal 脱离直径和 Cole 脱离频率。
`q_conv`、`q_quench`、`q_evap` 按论文 Table 1 的闭式关系分配；所需摩擦速度和
温度壁函数量在 `constant/phaseChange` 中显式给出。
`wallTemperatureModel heatFluxBalance` 将壁面热流作为边界约束，逐壁面单元求解
局部 `Tw`，保证 `q_conv+q_quench+q_evap=qWall`；字典中的 `334.25 K` 是求根
参考值，不再与恒热流同时作为独立约束。相间动量和热传递采用
Schiller–Naumann、constant lift、constant virtual mass、Antal wall lubrication、
GDB turbulent dispersion 和 Ranz–Marshall。

初始壁面状态的相关式结果为：

- `d_w = 3.72902171577e-4 m`
- `f = 181.848113225 1/s`
- `n'' = 9.57022291743e4 1/m2`
- `m'' = 3.02881147791e-2 kg/(m2 s)`
- `q_conv = 50365.4976216 W/m2`
- `q_quench = 21957 W/m2`
- `q_evap = 3937.45492128 W/m2`
- `q_total = 76259.9525429 W/m2`

## 运行与输出

~~~sh
./build/sonicSolver test/rpiWallBoilingCase
~~~

求解结果写入 `result/`。VTK 除每相
`alpha/rho/U/T/phaseMass/phaseEnthalpy` 外，还写出每相
`massSource`、`momentumSourceX/Y/Z`、`energySource`，以及统一壁面沸腾
ledger 的 `wallBoilingMdot`、三部分热源、气泡脱离直径和相间机械加热。
液相湍流方程还输出 `k.liquid`、`epsilon.liquid` 和 `mut.liquid`。

当前回归推进 `20` 步到 `2e-4 s`，`maxDeltaT=1e-5 s`。`2e-5 s` 和
`5e-5 s` 在首步都会被温变物性闭合拒绝，说明此工况的时间步主要受相焓/壁面源
尺度限制，不是常规对流 CFL。`CFL` 控制各相输运，`phaseSourceCFL` 限制质量、
相焓和强阻力源时间尺度，`maxDeltaT` 是用户硬上限。

本次 20 步回归中压力残差保持在约 `1.0e-9`，`sum(alpha_k)` 误差为零，
HYPRE 压力矩阵结构只建立一次并复用 20 次。两条湍流方程各自建立一次矩阵
结构，随后只更新系数和 RHS；每步两条方程共 2 次 Krylov 迭代，最终相对残差
约 `3.2e-12`。最终 `alpha.vapor` 范围为
`9.173e-4` 至 `1.0139e-3`；加热段 `y=2 m` 处，近加热壁气相体积分数约
`1.01299e-3`，近左壁约 `0.99929e-3`，已经出现由壁面汽化驱动的横向相分布差异。
最终 `k.liquid` 为 `4.99966e-3–5.02638e-3 m2/s2`，
`epsilon.liquid` 为 `1.69978e-3–1.71140e-3 m2/s3`，全部有限且高于字典下界。

## 适用边界

本算例已包含论文的流向分段和温变 R12 物性，但长方形几何不再是 DEBORA
圆管定量复现。此外仍没有 population-balance/size-group 方程和实验原始径向
剖面。

本算例通过 `constant/turbulenceProperties` 的 `phases (liquid)` 开启液相
RAS/k-epsilon；两条守恒输运方程进入同一个 PIMPLE outer corrector，对流使用
液相 canonical mass flux，`mu_t` 同时耦合进液相动量和焓扩散。HYPRE 设置来自
`system/solverProperties.linearSolvers.turbulence`。

同一方程系统也接受 RAS/k-omega SST 和 LES/Smagorinsky。算例已提供 `0/omega`
以及 `walls (leftWall heatedWall)`，可直接用于 SST 的显式壁距构造；默认仍保留
k-epsilon，保证下述 RPI 回归结果具有固定基线。

当前壁面 `k/epsilon` 仍是显式 zero-gradient 回归边界，没有标准壁函数、低雷诺数
近壁处理或气泡诱导湍流源；GDB turbulent dispersion 也尚未改为直接消费本次求解
得到的 `mu_t`。因此该结果证明“湍流方程已守恒装配并参与耦合流程”，但仍不能当作
完整的 DEBORA 湍流定量复现。
