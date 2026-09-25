# Phase 26 — 常密度压力约束数值算子

## 1. Scope

**Implemented**：串行、单块、层流、常密度 `Momentum + div(U)=0` 的 PISO 数值执行。入口仍为 `ExecutableEquationSystem → CompiledSolvePlan → PlanExecutor → OpRegistry`。Phase25B 的保守变量压力修正、密度显式推进、Eulerian 与 IBM 路径保留。

## 2. Mathematical equation change

常密度系统把 `rho=rho0` 作为模型常量；主未知量为 `U`，压力 `p` 为约束乘子。动量预测采用 `ddt(U) + div(momentumFlux) + diffusion(nu,U) = 0` 中实际存在的项；`nu=mu/rho0`。空间对流选择新增且显式命名的 `primitiveUpwind1`，扩散选择 `central2Explicit`。该新路径的数值精度不等于既有 WENO 路径。

## 3. Continuity transport vs incompressibility constraint

原常密度 preset 中伪装成方程的 `E_CONTINUITY` 已移除。Raw system 只含物理 `E_MOMENTUM` 与 `C_INCOMPRESSIBILITY: div(U)=0`；`rho` 无 runtime transport storage，不推进 `rhoE`。

## 4. Derived algorithmic equations

压力约束 transformer 从 raw Momentum 定义导出 `E_MOMENTUM_PREDICTOR`，并生成压力修正方程、`pPrime` workspace 和 correction operations；用户不输入 PressurePoisson。Compiled resource 使用 `U` 的真实 storage key/offset，常密度 workspace 为 `HbyA`，保守变量路径继续使用原有 residual。

## 5. Shared SIMPLE/PISO/PIMPLE operations

三种 preset 的数学 operation ID 共用 `momentum.assemble/solve`、`pressure.assemble/solve`、`pressure.update.prepare`、`velocity.correct`、`flux.correct` 等。数值算子不读取 preset 名称或持有 corrector loop。

## 6. Momentum predictor implementation

`PressureOperators` 每个公开方法对应一个 Op。`momentum.assemble` 基于旧 `U` 和已编译的对流/扩散 term 计算显式 `HbyA` 与 `rAU=dt/rho0`；`momentum.solve` 只做一次 `U*=HbyA-rAU∇p_old`，后续 PISO 修正不重复完整 `dt` 动量推进。workspace 由单 patch 执行期持有。

## 7. Compiled Momentum term consumption

provider 检查并消费 `BoundTerm`：`div(momentumFlux) → primitiveUpwind1`，可选 `diffusion(nu,U) → central2Explicit`。Capability resolution 还验证 raw Momentum 恰有 `ddt(U)`、该 divergence、至多一个 diffusion、零右端；额外项、额外物理方程或约束不被静默忽略，而是 **Unsupported**。当前尚无常密度 source-term numerical provider。

## 8. Pressure equation discretization

离散压力修正来自 corrected face-flux 的连续性缺陷和每面 `rAU` pressure response：内部相邻行形成正对角、负邻接系数的稀疏系统；固定压力边界提供 Dirichlet correction 条件。该矩阵无 `1/(rho c² dt²)` 项，不从旧可压缩 corrector 复制数学公式。

## 9. Rhie–Chow / face flux design

`CompiledNumericalSystem::pressureFaceCoupling` 对具有常量 rho 与压力乘子的系统选择 `RhieChow`；`flux.correct` 独立绑定 `flow.rhie-chow`。独立的 `SF_rhieChow.h` 计算面预测通量、压力梯度差响应及修正。`PressureOperators` 持有 `predictedFlux`、`faceResponse` 与跨 corrector 保留的 `correctedFlux`；下一次修正直接读取已修正面通量。单元测试覆盖均匀 p、均匀 U、线性 p 与 checkerboard 面差。

## 10. Velocity correction

`velocity.correct` 按 `U ← U-rAU∇pPrime` 更新已实现的语义 U view，随后应用 velocity BC。它不读 `rhoU/rho`，也不假设五分量保守布局。

## 11. Flux correction

`flux.correct` 对每个物理面只执行一次 `F ← F-faceResponse·ΔpPrime`，保留此 `F` 作为下一 pressure corrector 的输入。面索引采用 workspace 中唯一的 `(axis, lowerCell)` 标识；没有面通量平均操作。

## 12. Pressure reference

全 Neumann 压力系统在配置的 `referenceCell` 钉住 correction 的零空间，更新后整体平移压力到 `referencePressure`。固定压力边界按 Dirichlet 处理，不额外强加无关的 gauge。腔体测试初值 `p=101300`、参考值 `101325`，两步后参考点恰为 `101325`。

## 13. Linear backend binding

压力算子只构造 `LinearAlgebra::SparseSystem` 并调用 `SolverSession`。当前测试配置选 HYPRE；算子文件不包含 HYPRE 或 MPI 头。`explain` 的 operation/provider 与线性 backend requirement 均可见。

## 14. State realization binding

`CompiledStateRealization` 传递 `U[3@0,velocity]`、`p[1@0,pressure]`、`rho[1@0,rhoConst=1]`。`Field` 的三分量数组承载 U，独立 `ScalarField` 承载 p，rho 仅为 derived constant。旧网格初始化仍经过既有 Field 初始化步骤，随后在 application composition 时改绑三分量 U；这不是第二份运行时物理状态。

## 15. Turbulence decoupling

压力数值算子不 include、不查询 turbulence 类型。当前 turbulence 引入额外方程/term 的组合在 provider resolution 为 **Unsupported**；未来需相应 term/provider 和 plan 支持，不能静默跳过。

## 16. MRF/source decoupling

算子无 MRF、gravity 或其它模型分支。当前 Momentum source 因缺少常密度 source numerical binding 而 **Unsupported**。加入模型时应由 composition 生成 term，再给 term 实现 provider。

## 17. IBM decoupling

Ghost/ILW 及约束 IBM 的既有 production 路径未改。常密度 pressure provider 不直接调用 IBM；本阶段常密度+IBM 为 **Unsupported**，现有保守变量 Ghost 20 步对照保持一致。

## 18. MPI boundary

常密度压力 provider 仅对 serial single-block 解析为 Runnable；多块/MPI 为 **Unsupported**。新算子无 raw MPI；后续必须经现有 ExecutionRuntime 和 distributed linear algebra 明确验证 halo、canonical COPY、row ownership 与 residual SUM 后再开放。

## 19. PISO plan execution

PISO 的两次压力修正完全由 `CompiledSolvePlan` 的 pressure-corrector Loop 表达。`SingleFluidStepper` 仅把 operation ID 绑定到一次性 numerical leaf；`PlanExecutor` 执行 loop。`explain` 同时显示 raw equations、derived operations、`RhieChow`、`flow.pressure-operators`/`flow.rhie-chow`、plan 与 `runnable`。

## 20. SIMPLE status

**Unsupported**：现有显式 momentum predictor 是真实物理时间推进，没有 fixed-point predictor、relaxation/convergence 数值语义。不得将其在 outer loop 内重复使用。

## 21. PIMPLE status

**Unsupported**：同上；虽然 pressure-corrector kernel 可复用，outer fixed-point 及 convergence semantics 未实现。

## 22. Cavity validation

新增 `test/pressure/constantDensityPiso`。80×80、移动顶盖、两步 forward Euler、每步两次 PISO correction；压力与速度有限，顶盖 Ux=1、底壁 Ux=0，rho 始终为 1。最终压力 checkerboard parity/RMS 约 `4.6e-5`，低于测试阈值 `0.05`。这是短时稳定性/耦合测试，不宣称腔体 benchmark 精度。

## 23. Poiseuille validation

新增 `test/pressure/constantDensityPoiseuille` 及双网格回归。通道 `-0.5≤y≤0.5`，`mu=0.1`、`rho=1`、左右压差 `0.1`；解析式 `Ux(y)=0.5(0.25-y²)`。同一物理时间约 `9.766` 的中心截面结果：

| 网格 | L2(Ux−Uexact) | Linf | 最大 Ux | 压差 |
|---|---:|---:|---:|---:|
| 16×16 | 7.4524e-4 | 9.8684e-4 | 0.126003 | 0.100000 |
| 32×32 | 3.5327e-4 | 4.7575e-4 | 0.125480 | 0.100000 |

解析中心速度为 `0.125`。误差随网格加密下降，未为解析曲线调参。两级 run 在 `t=10` 末的最大离散连续性缺陷分别约 `2.24e-9`、`1.20e-9`。

## 24. Divergence/correction diagnostics

腔体 step 1：correction #1 的 maxDiv `2 → 3.39853e-8`，#2 保持 `3.39853e-8`；step 2：#1 `0.64824 → 4.63737e-8`，#2 保持 `4.63737e-8`。第二次修正的已收敛线性系统可返回零迭代；测试只要求有效首修正降低缺陷。

## 25. Existing regression invariance

Phase25B 源码 commit `4972cc7` 在独立 `/tmp` 构建并作为对照。相同 case 下，4-rank Sod 20 步的 time/dt 与 state-closure 日志逐项相同；serial Ghost IBM 20 步的对应日志逐项相同；Eulerian 两步诊断逐项相同。`pressureConstraintPiso` 最终 VTS SHA-256 仍为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`，与 frozen legacy result 字节一致。WENO7/Rusanov Sod `t=0.2` 冻结回归通过。架构守卫未出现新依赖违反。

## 26. Remaining numerical limitations

**Implemented**：本阶段限定的显式一阶 primitive upwind + 可选二阶中心扩散、serial/one-block 正交网格 PISO、现有 fixedValue/zeroGradient/empty U/p 边界。**Interface-only**：pressure face-coupling selection 可扩展其它方法，但当前只有 Rhie–Chow。**Legacy**：保守变量压力修正仍由原 corrector 执行，保持冻结数值行为。**Unsupported**：SIMPLE/PIMPLE fixed-point、常密度 IBM/模型源项/额外方程、常密度 MPI/多块、非正交网格及修正、其它 primitive convection/face-coupling 方法。非正交网格在绑定时明确报错。Field 的底层网格初始化仍沿用旧流程，但运行时 U/p/rho authority 已分离。

## 27. Recommended next phase

Phase 27 可增加 SIMPLE/PIMPLE 所需的 relaxation、convergence 与 outer fixed-point plan semantics，复用本阶段 momentum、pressure、Rhie–Chow、velocity/flux correction kernels；开放 MPI 前需单独完成 distributed pressure matrix 和 face ownership 验证。

## 验证与所有权摘要

变更层级是 composition、transformation、numerical compilation、operation provider 与 state realization；性质是“改变怎么解方程”，同时把常密度连续性表达为约束。State authority 由旧保守 Q 假定改为 U/p/rhoConst 语义绑定；workspace 由新算子持有，plan 仍拥有时间与 correction 顺序。公开 timestep API 未新增；依赖方向和既有 MPI COPY/SUM 语义未改变。新常密度路径的数值语义是新增方法，旧 production 数值路径没有改变。
