# A6：纯反应、多扩散、交叉扩散及向量方程

日期：2026-10-23。本轮扩展选定的显式数值能力，不重建 Math/WHAT，不宣称任意 PDE、隐式或分布式新能力完成。

## 1. 起点、范围与复现标识

起点分支 `codex/core-single-path-closure-20261022`，HEAD `761712ef1d4f3f140178c08c24954702ed15bb39`。工作树包含尚未提交的两轮 A5，另有六份用户旧报告修改；未 reset 或覆盖。A5 依赖与本轮分开提交，依赖提交为 `5d65dbdf21b17a8e6c64075fecf1ad2a2e9b8644`；六份旧报告不纳入提交。

- 本轮 before：`test/explicitA6/baseline/sonicSolver`，SHA256 `26da6e5bde9ddd94aa03d0148cfb1d163fb76690dacfa8a2f0f6fab3e2fb9ae4`。
- final：`build/sonicSolver` 与不可变 `test/explicitA6/candidate/sonicSolver`，SHA256 `329703d6ff050ca93ed364d65dc3d8986e28cb52f5915a49bd3fd7e6861b40f1`。
- 起点状态、diff、逐源码内容：`test/explicitA6/evidence/start-status.txt`、`start.patch`、`start-sources.json`、`start.json`。
- 对照的是实际未提交 A5 的 executable，不仅是 HEAD。新增能力原先 Unsupported，没有把首次可运行结果重新冻结成“旧算法基线”。

## 2. 实现与责任边界

调用链保持：typed contribution/现有 native AST → EquationRegistry/STATE/HOW/WHICH → ScalarMethod → 原 FormulaCompiler 与具体 Central2 provider → compiled physicalResources → 同一 bindScalarResources → 原同步 temporal group/PlanExecutor → 实际 kernel。

| 层/对象 | before → after |
|---|---|
| Math/WHAT | 同一签名和 AST → 不变，不新增数学分析器或展示方程 |
| 普通方法能力 | 一个 target diffusion、一个分量 → 零/多个常系数 diffusion、交叉 Stage symbol、标量/三分量向量 |
| 编译扩散布局 | 单一 diffusivity → diffusions 各自保存 occurrence/input/coefficient/components；kernel 独立冻结系数/slot |
| physical storage | 原 ScalarResources 数组 → 相同 owner，向量由三个唯一 ScalarField 分量组成一个 physical view |
| stage/workspace | solver old/stage/rhs/RK → 同一 owner，component-major 全分量索引 |
| 边界 | 每个 scalar 都要求 → 仅实际 stencil 输入要求自己的边界；向量 fixedValue 显式三分量、zeroGradient 分量独立 |
| HOW/clock | 原 group barrier/Commit → 不变，只有 StateBundle clock |
| 输出 | 原 registeredStateFields/writer → 复用已有多分量输出，不复制 physical state |
| 依赖/MPI | 原依赖和 COPY/SUM → 不变；没有新增 MPI 实现或并行能力 |

类别是“扩展所选方法怎样计算既有数学”，包括真实新增的反应/交叉扩散/向量数值能力；不是纯外观重构。旧 flow、transport/shared/SST 数学、floors、CFL、stage coefficients/abscissae、boundary/halo 顺序及 canonical COPY/GlobalDof SUM 未修改。

普通常系数允许有限的 signed 值，由作者决定数学；没有自动修改扩散矩阵、限制权重、假扩散、隐式降阶或替换算法。原 conservative/shared/SST 的 rho-weighted coefficient 检查仍保留。专用融合流体和湍流 kernel 未拆。

**公共 NumericalCompiler、PlanExecutor、中央模块选择链均无 A6 修改或新增名称分支。** 控制层只扩展现有 native field metadata 的 vector layout；binder 只分配/核验/注册，不求 RHS、不刷新 stage、不推进时间。

## 3. 删除项、API 与代码量

- 删除普通路径“必须有一个 diffusion”“只能 diffusion 自己”“只能 scalar target”的限制；transport/shared/SST 真实限制保留。
- 删除 CompiledScalarEquation 的单一 diffusivity 字段；不留兼容 getter/转发 wrapper。
- 删除 target-only scalarCentral2Kernel、旧 central2ScalarDiffusion(double) 和无人消费的 FormulaValues.neighbor；现役测试直接迁至 slot-aware provider。
- ScalarBoundaryCondition 新增真实消费的 componentValue；原 scalar 回调仍是其自身实现，没有新旧并行执行入口。
- 同一 symbol/storage 仍唯一 owner；重复 occurrence 不按 AST 去重。所有资源验证成功后才发布 borrowed views，失败不留下悬空注册。
- 原 negative-D “必须拒绝”断言不再适用于普通 signed 常系数，改为真实的非有限值/变系数能力负例；原数值参考与容差没有替换。

相对本轮 start-sources.json：**11 份已有 C++ 源码/头，+159/-102，净 +57 行；生产源码新增文件 0**。不含文档、CMake、测试。记录见 `evidence/source-metrics.json`。另增加一个集中 C++ 数值验收目标及测试/对照脚本；没有新增生产类、Manager、Registry、solver family 或第二份 FormulaCompiler。A5 的三个新增生产文件属于依赖提交，不重复统计为 A6。

## 4. 内部精度算例：实际设置

以下均在 `test/explicitA6/cases/<名称>/`；源输入是 `test_equations.cpp` 的 typed contribution、初值函数和具体 binder，生成 explain、metrics、clock；数值输出均为各 case/result。全部 serial、single patch，没有物体/IBM/ILW、EOS、流体 Q、湍流或物性 closure。scalar 无量纲，向量量纲 L/T；D 为 L²/T，反应率 T⁻¹。坐标/时间分别按 m/s 解释，数值均使用 SI 一致量纲。无 CFL 选步，实际固定 dt；无 MPI，没有假 halo 交换。

### 4.1 纯反应（8 例）

共同网格：5×1×1 物理节点，x=[0,1]，其余坐标为 0，NG=0。不是退化网格 CFD，方程无空间导数，不计算几何 Jacobian/邻域；无 BC、无 halo 需求。配置终止 t=1，实际均到 t≈1；初始 t=0/step=0。每例输出初始/最终两帧。

|case 名称|WHAT / IC|recipe、dt / 实际 steps|独立参考与最终 max error|
|---|---|---|---|
|linearReaction、linearReverse|a'=0.7b，b'=0.3a；a0=1,b0=2|RK4，.025 / 40|2×2 矩阵指数，2.52070365e-10；反转 authored HOW 后数组完全相同|
|logisticEuler0.100000|c'=2c(1-c)，c0=.2|Euler，.1 / 10|1/(1+4exp(-2t))，8.86252749e-3|
|logisticEuler0.050000|同上|Euler，.05 / 20|4.23017910e-3|
|logisticEuler0.025000|同上|Euler，.025 / 40|2.06463981e-3|
|logisticRK40.100000|同上|RK4，.1 / 10|7.84146788e-7|
|logisticRK40.050000|同上|RK4，.05 / 20|4.98318720e-8|
|logisticRK40.025000|同上|RK4，.025 / 40|3.13935533e-9|

时间减半，Euler 误差比约 2.10/2.05，RK4 约 15.74/15.87。不以“有限/不 crash”替代精度。双变量反转只测试共同 stage 的等价 RHS 读取，未将有数学依赖的顺序迭代随意重排。

### 4.2 空间算例（10 例）

共同为 x=[0,1] 的节点 Cartesian 网格，NY=NZ=1、NG=1，x 两端 fixedValue(0)，无其他活动轴/边界。Central2，RK4，dt=1e-4，配置 t=.05、实际 .05/500 steps，输出初始/最终两帧。物理初值和 referenced Stage 输入均直接写入同一 binder 创建的数组，无 partner copy。

|case|方程、初值、网格|独立参考/结果|
|---|---|---|
|multi17、multi33、multi65|a_t=.03 a_xx+.07 a_xx；a0=sin(πx)；17/33/65 节点|离散特征值 λh=4 sin²(πh/2)/h² 的解析半离散指数，误差 <2e-13|
|duplicateDiffusion|.05 a_xx+.05 a_xx；33 节点，同初值|独立总系数 .1 的指数参考，6.66134e-16；证明同表达式没有去重|
|crossDiffusion、crossReverse|D=[[.07,.02],[.02,.04]]；[a,b]0=[1,2]sin(πx)；33 节点|解析 2×2 exp(-λhDt)，误差约 1.11e-15；反向 HOW 数组完全相同|
|signedCross|D=[[.07,-.02],[-.02,.04]]，其余同上|独立矩阵指数，1.11e-15；适定矩阵的负交叉项未被修改|
|vector17、vector33、vector65|V_t=.03 V_xx-.3V；V0=[sin(πx),-2sin(2πx),.5sin(3πx)]；17/33/65 节点|逐分量连续解 Ac sin(mcπx)exp[-(.03(mcπ)²+.3)t]；max error=1.64494428e-3、4.14226459e-4、1.03744058e-4，空间误差比约 3.97/3.99|

每个 cross 方程读取对方的空间导数；半离散参考不调用生产 stencil/RK，连续向量参考另行验证空间收敛。矩阵不是从数值结果拟合。非局部向量 dot/cross 等没有在此被实现。

## 5. 普通 native IO 六例

路径 `test/explicitA6/native/<名称>`，完整输入由 `native.py` 生成，严格沿 solvers/models/fields/mesh 单目录注册。执行 explain 后实际 run，最终读取 VTS 数值，不只检查配置可通过。共同 serial、单块，models registry 为空，无 EOS/Q/IBM；输出初始和最终。无修改正式案例。

网格 x=[0,1]、y=z=[0,.1]，无 scale，单位 m；reaction/logistic 为 5×3×3 节点（4×2×2 cells），其余 33×3×3 节点（32×2×2 cells）。NG=1 是现有 SFM loader 的要求；反应方程仍声明 halo=0、无 boundaries 注册文件。Loader 的体积/NG约束没有为 A6 放宽；零 NG 由第4节内部执行证明。

空间例初值沿 x 由显式 node mesh sets 给出，y/z 同值。除 vectorFixed 外，x 两端各 target fixedValue(0)，其余四面 zeroGradient；vectorFixed 的 x fixedValue=[1,-2,.5]。每个交叉输入自己的初值/边界均绑定。time recipe 全 classicalRK4。反应 dt=.025、endTime=1、40 steps、writeInterval=40；空间 dt=1e-4、endTime=.05、500 steps、writeInterval=500；实际终点与配置一致。

|case|实际 WHAT/STATE 与初值|独立最终误差|
|---|---|---|
|reaction|第4.1节双变量线性反应，a=1,b=2，标量|2.52070365e-10|
|logistic|a'=2a(1-a)，a0=.2，标量|3.13935533e-9|
|multi|.03 a_xx+.07 a_xx，a0=sin(πx)|6.66134e-16（半离散）|
|cross|第4.2节正交叉矩阵，[1,2]sin(πx)|1.11022e-15（半离散矩阵指数）|
|vector|Vexperiment_t=.03ΔVexperiment-.3Vexperiment，三种幅值/频率同第4.2节|2.44249e-15（各分量半离散指数）|
|vectorFixed|Vexperiment_t=.03ΔVexperiment；全域 [1,-2,.5]；显式 all 初始化 set、非零分量壁值|0；独立常量解，验证 boundary componentValue 与 output 三分量|

向量输出实际 `NumberOfComponents=3`，不冒充 U。纯反应没有零/极小 D。命令和 binary SHA、最终时间、输出帧及误差在 `evidence/native.json`；各目录 explain.log/run.log/result 保留。

## 6. 原功能 before/after 保护：12 对

共同输入来源与差异在脚本明确复制：`module_cases.py` 从保留的 A5 stageSST 配置构造八个独立副本；`protect.py` 复制四份 A5 保护输入。仅 before/after binary 不同，同配置的两个副本不混写结果。以下均 serial、单 patch，无 IBM/ILW；所有 output/result、配置、命令、二进制哈希和日志在 test/explicitA6 下。

### 6.1 八例模块保护

路径 `before/<name>` 和 `after/<name>`。共同域 [0,100]×[0,20]×[0,.1] m，100×20×1 cells（101×21×2 节点），NG=4。PerfectGas γ=1.4、R=287.05，μ=1e-5 Pa·s、Pr=.72；rho=1.2，U=[.3,0,0]，p=100000，k=.01、epsilon=.001、omega=1 的文件配置按启用模型消费；闲置字段不宣称被积分。

原 nativeRight/nativeSides 的 rho/U/p zeroGradient，frontBack empty；nativeLeft 对这些 flow 量没有显式 BC，保留现有默认/继承行为，不能称 fixed inlet。k/omega 的 nativeLeft/Right/Sides zeroGradient，epsilon 仅 Right/Sides 显式 zeroGradient、frontBack empty，Left 未显式指定。WHAT 为原 mass/momentum/energy 与启用时原 once-step k/epsilon 或 k/omega source/diffusion，未新增 RAS convection。WENO5/Rusanov、Central2。

全部 startTime=0、endTime=.001、maxDeltaT=.0001、CFL=.5；实际 10steps/t=.001，每步输出，共11帧。

|name|具体启用/覆盖差异|flow recipe|
|---|---|---|
|absent|turbulence 未注册，gravity false|Euler|
|laminar|明确 laminar，gravity false|Euler|
|gravity|laminar + g=[.01,0,0]；原动量及机械功贡献|Euler|
|kepsilon|RAS kEpsilon、legacyOnceStep、无 gravity|Euler|
|kepsilonGravity|上例 + g=[.01,0,0]|Euler|
|disabled|kEpsilon configured but enabled:false；无专属回调/资源|Euler|
|override|kEpsilon once-step + 显式已有合法 WHICH 覆盖；无 gravity|RK4|
|sstOnce|kOmegaSST legacyOnceStep；无 gravity|RK4|

开关 gravity 改变物理输入，不要求开/关相等；本轮只要求每项相同输入 before/after 一致。A5 的 source 预算证据仍为历史，不冒称本轮重做预算。

### 6.2 四例执行保护

|path 相对 protection/{before,after}|具体设置与执行终点|
|---|---|
|flowPair|来源 resourceA5/composition/flowPair；[0,1]×[0,.25]² m，33×5×5 节点，NG=3；PG γ1.4/R287.05/μ0；rho1,U=[.1,0,0],p1，flow 六面 zeroGradient。dye_t=.03Δdye+.7salt，salt_t=.08Δsalt+.7dye；dye 左半1/其余0，salt 右半2/其余.5；dye 左1右0，其余zeroGradient，salt全zeroGradient。RK4，WENO5/StegerWarming，仅标量Central2；dt=.0005、end=.005，实际10steps/11帧|
|scalarPair|与 flowPair 同网格、同 scalar 方程/IC/BC/recipe/终止；仅移除流体三方程及其物理资源，不初始化Q/EOS。实际10steps/.005/11帧|
|stageSST|来源 resourceA5/after/stageSST；[0,100]×[0,20]×[0,.1]，101×21×2节点，NG=4；PG1.4/R287.05、μ1e-5/Pr.72，rho1.2,U=.3,p1e5,k.01,omega1，epsilon0闲置；flow BC同6.1，k/omega 左也zeroGradient。原 stageConservative SST + shared mass flux，RK4/WENO7/Rusanov/Central2。start0/endTime=.001/CFL.5，maxDeltaT未在该runtime显式设置；endStep=1，实际 t=7.316037672350e-4/1step/2帧；这是首步保护，不是跑完endTime|
|eulerianSST|来源 resourceA5/after/eulerianSST；[0,1]×[0,.2]×[0,.01]，21×5×2节点（20×4×1 cells），NG=3；water/air相质量加权momentum/energy/phase与共享pressure，SST各相k=.01/omega1。water:rho998.2,alpha.9,T300,U0,μ.001002,Cp=Cv4182,kcond.598；air:rho1.1766,alpha.1,T300,U0,μ1.846e-5,Cp1007,Cv718,R287.05,kcond.0262。p101325；Left p/rho/alpha/T/k/omega固定上述值，U0；RightzeroGradient；top/bottom U0，其余zeroGradient；frontBackempty。MRF/all,center[.5,.1,0],axis z,omega1；SchillerNaumann(d=.001)、virtualMass.5、RanzMarshall，lift0，无wall lubrication/dispersion。Euler+PIMPLE outer3/pressure2/nonorth1，relax U.7/p.3，原Upwind/Central2；HYPRE flexGMRES/BoomerAMG max100/rel1e-8/abs1e-12/krylov30。start0/end=.0002/maxDt=.0001/CFL.5/endStep1，实际3.624299498372e-5/1step/2帧。不是长期多相湍流验收|

全部12对 **result 数值文件字节一致、Plan/clock trace一致、物理体积梯形积分 mass/momentum/energy、rho/p extrema 和可用场范数一致**。实际指标见 `evidence/regression.json`；trace不比较路径或存储地址。未为此新增完整 IBM CFD 测试，相关数值/绑定没有迁移；集中 viscous/IBM 契约保护列入 CTest，不能代替各 IBM 方法完整 no-slip 验收。

## 7. 测试、负例及明确限制

完整 configured build 成功，最后 build-confirm.log 为 no work to do；compile/link failure=0。本轮最终 **45/45 定向 CTest**，包括新增精度目标、A5 contribution/resources、A/B 数学覆盖、standalone/multi/coupled/flow scalar、stage coefficient/SST、once-step RAS、Eulerian、pressure、producer/resource lifetime 和相关 IBM 契约。不是全仓任意所有测试通过。新 C++ 精度目标另在默认 CTest 工作目录运行1/1，以保留 cases/ 下的最终18份数据。

7个 native负例：partner缺Right BC、vector初值只有两分量、未注册internalField、系数量纲错误、普通variable D、SSPRK3无实现、MPI后端未实现，分别在编译/真实binding阶段拒绝，无物理step提交。另有集中 rank/签名、未声明STATE、向量cross缺numerical provider、绑定原子性负例。负例日志是正确拒绝证据，不算完整PDE或物理验收。

架构检查通过：709源文件、2条既有allowlisted methods→solver依赖、11组同名头；没有新反向依赖。无新A/B类fallback，不修改原合法floor/正则化；普通signed coefficient是显式能力扩展，不是数值修补。

**未覆盖/Unsupported**：新的 MPI/multipatch scalar/vector；普通变系数/各向异性/曲线网格；三阶以上或一般tensor PDE、dot/cross执行；向量 known-source callback；非单位质量/隐式/Jacobian/DAE/Subcycle；所有phase/face/constraint资源自动装配、任意外部公式、production KKT/新IBM数学。forwardEuler/classicalRK4之外的recipe保持拒绝。

**历史问题保留，未在A6重跑/解决**：原压力serial/MPI shape差 `1.0851072147488594e-6 > 1e-6`；Level-set历史isolated-kernel参考失配。不改阈值、不替换参考、不用本轮45/45概括为它们通过。本轮不扩展并行能力，以native拒绝及集中原契约测试保护，不以serial数值代替MPI验收。

## 8. 命令、证据与后续

从仓库根目录（所有TMPDIR/测试cwd在test/内）：

```sh
DEVELOPER_DIR=/Library/Developer/CommandLineTools TMPDIR="$PWD/test/explicitA6/evidence" cmake --build build --parallel 4
python3 -B tools/check_architecture.py
python3 -B test/explicitA6/protect.py before
python3 -B test/explicitA6/module_cases.py before
python3 -B test/explicitA6/run_tests.py
python3 -B test/explicitA6/native.py
python3 -B test/explicitA6/negative.py
python3 -B test/explicitA6/module_cases.py after
python3 -B test/explicitA6/protect.py after
python3 -B test/explicitA6/verify.py
DEVELOPER_DIR=/Library/Developer/CommandLineTools TMPDIR="$PWD/test/explicitA6/evidence" ctest --test-dir build -R '^explicitEquationCapabilities$' --output-on-failure
```

本地 before/after 脚本使用冻结binary及保留的A5输入；新克隆无法重建未上传的历史二进制/原始结果。新增 C++ 目标和 native.py 自行构造全部新能力算例，默认使用当前 build/sonicSolver，也可用第一个位置参数指定冻结 executable，不依赖历史数值场。结果、日志、build、binary、大型源码快照不提交；必要源码、生成脚本及文档提交到工作分支，不合并 main/master，不创建PR。

最终证据：build-final/build-confirm、architecture-final、ctest-final、equations-final、native.json、negative.json、regression.json、source-metrics.json，以及各case输入/clock/explain/run/result。浮点数据对照用实际数值，没有用自动色标云图夸大微差。

后续分别推进：明确需求的 variable/aniso coefficients → 非单位质量与具体隐式/块实现 → 对应真实资源与并行绑定。不能把“Typed WHAT 可表达”当作这些能力已有实现；不在本轮继续设计新全局框架。
