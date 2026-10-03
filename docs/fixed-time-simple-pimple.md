# Phase 26B / 27 — 固定物理时间的 SIMPLE/PIMPLE 压力约束执行

## 1. Scope

**Implemented**：收窄 Phase26 常密度压力算子的配置依赖；为 serial、single-block、正交网格的 `rhoConst + Momentum + div(U)=0` 增加 SIMPLE/PIMPLE 固定物理时间迭代。没有新增 solver family。**Unsupported**：本报告第 31 节列出的其它组合。

## 2. Phase 26 starting architecture

原 PISO 已由 `RawEquationSystem → pressureConstraint transformation → ExecutableEquationSystem → CompiledNumericalSystem → CompiledSolvePlan → OpRegistry → PlanExecutor` 执行。Phase26 的 `PressureOperators` 接收整个 `SolverConfig`，SIMPLE/PIMPLE Plan 使用具名缺失 leaf 报告 Unsupported，原因是原 `momentum.solve` 不能在 outer loop 中重复以当前 U 为时间基推进完整 `dt`。

## 3. PressureOperators raw-config cleanup

`PressureOperators` 现在只接收不可变 `PressureOperatorConfig` 值：`mu`、两种 relaxation factor、pressure reference、pressure linear config、U/p typed boundary；数值 recipe 和 CFL 仍读取已有 `CompiledNumericalSystem`，`rhoConst` 来自 `StateRealization`。它不能访问完整 `SolverConfig`、turbulence、IBM、MRF 或 phase 配置。绑定时验证 Momentum 的每个已编译 term；一个 `primitiveUpwind1` divergence 和至多一个 `central2Explicit` diffusion 有 provider，Source 或其它 term 明确失败。

## 4. Pressure namespace/ownership cleanup

仅迁移 Phase26 新算子到 `src/solver/algorithm/pressure/SF_pressureOperators.*`、`SF::Pressure`。Rhie–Chow 仍在 `solver/discretization/pressure`，namespace 改为 `SF::Pressure::RhieChow`。Legacy `PressureBased::Corrector`、KKT 与 Eulerian 实现均未移动。

## 5. Provider binder cleanup

`SF_pressureProviderBinding.cpp` 根据已解析的 `pressure.prepare` binding 一次性创建常密度新算子或保守变量 legacy corrector；未知 provider fail-fast。`SingleFluidStepper` 不再比较 concrete provider 名称，也不选择 coupling preset。Operation callbacks 仍由 Stepper 注入 `OpRegistry`，PlanExecutor 执行编译的 Plan。

## 6. Pressure reference semantics

`validatePressureCorrectionConfig` 对 multiplier pressure 只要求 `referencePressure` 有限，`0` 合法。全 Neumann 仍在 pressure correction 后做均匀 gauge 平移；SIMPLE/PIMPLE 压力 under-relaxation 后再次恢复 gauge，均匀平移不改 `grad(p)`。零 gauge 的两步 cavity 运行及参考单元断言通过。

## 7. Pressure boundary contract cleanup

Typed `BoundaryConfig::energyFromPressure` 改为 `BoundaryConfig::pressure`；native `p` field 的兼容解码在 application 层进入该 typed 成员。边界数值类型、施加顺序和 legacy corrector 接收的值未变。历史输入字段名仍由既有 compatibility decoder 处理。

## 8. Fixed-time iteration model

`pressure.prepare` 计算一次 `dt`；固定点 Plan 随后以 `pressure.step.begin` 冻结 `U_n`。每个 `pressure.iteration.begin` 保存当前 `U_k,p_k`，共用的 Momentum/pressure/face-flux operations 更新 working U/p，`relaxation.apply`、`convergence.evaluate`、`iteration.end` 完成一次 outer。循环结束后只有一个 `pressure.step.commit` 增加 time 和 step。

## 9. Physical state vs iteration state

`StateBundle` 保持唯一 physical-state/clock authority。已有 U/p Field storage 在 step 内承载 working iterate；`PressureOperators` 只拥有 `U_n` 快照、上一 outer 的 U/p、标量诊断和既有压力/面通量 workspace。没有复制 `StateBundle` 或新建第二份全场物理状态 authority。

## 10. Temporal term semantics

首次 outer 的 U 与 `U_n` 相同，保留 Phase26 PISO 的原始浮点求值顺序。后续 outer 的空间通量、扩散和边界读取 `U_k`，但 `ddt(U)` 的时间基始终是冻结的 `U_n`：`HbyA = U_n + dt · RHS_spatial(U_k)`。单元测试用 `(U-U_n)/dt+aU=b` 验证第二次 outer 的候选为 `2.5`，而错误地把首次结果作为时间基会得到 `3.5`。

## 11. Iteration workspace ownership

`baseVelocity_`、`iterationVelocity_`、`iterationPressure_` 仅在 fixed-time Plan 首次 `step.begin` 按 patch layout 分配；普通 PISO 不分配这些数组。`oldFlux` 保留 Phase26 的局部生命周期，以保护冻结 PISO 求值与输出哈希。每个 outer 只覆盖已有 workspace。

## 12. SIMPLE formulation

Raw system 仍为 `rho=rhoConst`、物理 `E_MOMENTUM` 和 `C_INCOMPRESSIBILITY`。SIMPLE 注册的是同一个 pressure-constraint transformation；导出的 `E_MOMENTUM_PREDICTOR`、`E_PRESSURE`、U/p correction 与 PISO 同源。

## 13. SIMPLE SolvePlan

`pressure.prepare → pressure.step.begin → Loop outer { iteration.begin → momentum.assemble/solve → Loop pressureCorrectors { pressure boundary/assemble/solve → velocity/pressure/flux correct → correction.commit } → relaxation.apply → convergence.evaluate → iteration.end } → pressure.step.commit`。目前 `Loop` 是配置的固定最大次数，没有 hidden while loop 或提前退出。

## 14. SIMPLE relaxation mathematics

第一版为明确的 solution relaxation：`U_new=U_k+alpha_U(U_candidate-U_k)`、`p_new=p_k+alpha_p(p_candidate-p_k)`，只对 interior U/p 执行，再恢复 typed BC 与 Neumann gauge。`alpha=1` 直接返回 candidate，保证 PIMPLE outer=1 与 PISO 的浮点结果等价。`momentumRelaxation` 和 `pressureRelaxation` 是 solve policy 配置，不是 Equation term，也不进入 provider identity。已用纯算术单元测试与 `alpha=0.8` 的 Poiseuille 运行验证。

## 15. SIMPLE convergence

`pressure.convergence.evaluate` 记录 U 与 p 的最大 outer 差、按上次迭代幅值归一的 `maxDelta`、以及 corrected face-flux 的最大连续性缺陷；非有限值 fail-fast。固定次数是本阶段的执行语义，diagnostic 不触发提前停止。Poiseuille 513 个物理步骤中，每步第二次 outer 的 `maxDelta` 均小于第一次。

## 16. PIMPLE formulation

PIMPLE 贡献相同的 pressure-constraint formulation 和 executable operations；没有 PIMPLE 方程副本，也没有 `PimpleStepper`。

## 17. PIMPLE SolvePlan

PIMPLE 的 outer 是固定时间的 iteration loop；其内 pressure-corrector loop 可以多次执行相同 pressure/velocity/flux operations。当前 non-orthogonal loop 仍在 Plan 中，但非零 non-orthogonal 数学对新常密度 provider 为 **Unsupported**，不能冒充 Runnable。

## 18. PISO unchanged path

PISO Plan、`momentum.assemble` 的原浮点求值顺序、Rhie–Chow、压力矩阵、速度/面通量修正、pressure reference、CFL 和 HYPRE 设置均未改。独立 Phase26 构建与当前版本的常密度 PISO 3 个 VTS 文件哈希逐一相同；cavity 与双网格 Poiseuille 原测试通过。

## 19. Shared numerical operators

SIMPLE/PISO/PIMPLE 均引用同一组 `momentum.assemble/solve`、`pressure.assemble/solve`、`pressure.update.prepare`、`velocity.correct`、`flux.correct`、`pressure.correction.commit` 和 `pressure.step.commit`。固定点只多 `step.begin`、`iteration.begin/end`、`relaxation.apply`、`convergence.evaluate`。没有复制压力矩阵或 Rhie–Chow 数学。

## 20. Rhie–Chow independence

`SF_rhieChow.h` 不读取 preset、outer count 或 coupling config；`NumericalCompiler` 仍绑定 `RhieChow` face coupling，`flux.correct` 继续解析为 `flow.rhie-chow`。

## 21. Linear-backend independence

压力 operation 继续只使用 `SolverSession`/`SparseSystem`；HYPRE 是编译出的 linear backend requirement，Plan 与 provider 不 include HYPRE 或 `mpi.h`。

## 22. Turbulence/MRF/IBM decoupling

新算子不按模型名 dispatch。尚无数值 provider 的 Momentum source/额外方程以及 constant-density IBM 继续 **Unsupported**；既有 Ghost/ILW/KKT production 数学不变。

## 23. MPI boundary

新常密度 pressure 路径仍只解析 serial single-block。分布式压力矩阵、halo、face owner 和 row ownership 未在本阶段实现；4-rank Sod 只验证现有保守变量 MPI 路径，COPY/SUM 语义未变。

## 24. Provider resolution

Resolver 依据 Plan 所需 OpId、typed capability、state realization、编译 term 与 serial/patch 条件判定 Runnable。常密度 SIMPLE/PISO/PIMPLE 在受支持配置下分别 Runnable；保守变量 SIMPLE/PIMPLE 因缺 fixed-time provider 为 **Unsupported**。架构测试逐项确认三种 preset 的 pressure assemble/solve、velocity correct、flux correct 绑定同一 provider，并确认常密度 MPI 与非正交修正仍为 **Unsupported**。`explain` 能同时显示 raw、derived equation/operation、三层 Plan、provider 与运行状态。

## 25. SIMPLE numerical validation

以 Phase26 Poiseuille 为源，仅改变 preset、outer=2、pressureCorrectors=1、两项 relaxation=0.8。16×16、512 步时，解析速度 L2 误差 `0.00822508807`，左右压差 `0.1`，最后 corrected face-flux 最大连续性缺陷 `2.40635e-10`；全部 U/p 有限，513 步的第二次 outer normalized delta 均下降。

## 26. PIMPLE numerical validation

两步 cavity 设置 outer=2、pressureCorrectors=2；trace 每步有一次 `step.begin`、两次 `iteration.begin`、四次 pressure correction、一次 `step.commit`，通量连续性保持有限并降低。物理时间序列与 PISO 均为 `0.000390625, 0.00078125`，每步 `dt=0.000390625`。

## 27. PIMPLE→PISO degeneration test

PIMPLE outer=1、pressureCorrectors=2、relaxation=1 与相同 PISO 两步 cavity 的 U/p/rho 数组逐元素在 `1e-12` 内一致；压力 reference 相同。fixed-time 的第一次 outer 使用与 Phase26 PISO 相同的求值顺序。

## 28. Fixed-time temporal-base test

`test_fixedTimePressure` 验证两次 outer 的时间基恒为 `U_n` 并验证 solution relaxation；`fixedTimePressureRegression` 检查 Plan trace 中每个物理步只有一次 base capture、一次 commit，outer count 不改变 time/dt。

## 29. Phase 26 PISO regression

`constantDensityPressureRegression` 与 `constantDensityPoiseuilleRegression` 通过。Phase26 HEAD 独立构建的常密度 PISO 对照：初始与两次 commit 的所有 VTS SHA-256 均相同；逐条 step/closure 诊断也相同。`referencePressure=0` 的两步 cavity 另行运行，参考单元最终恰为 0。

## 30. Production regression invariance

当前完整构建成功，`ctest --test-dir build --output-on-failure --parallel 1` 为 **19/19**；architecture checker 通过。Phase26 HEAD 独立编译、相同 case/配置对照：4-rank Sod 20 步、serial Ghost IBM 20 步、Eulerian 2 步、legacy `pressureConstraintPiso` 与新常密度 PISO 的 VTS 哈希均一致；Sod/Ghost/PISO 的逐步诊断也一致。legacy pressureConstraintPiso 最终 SHA-256 保持 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`。WENO7/Rusanov frozen Sod CTest 通过。MPI 是在允许 PRTE socket 的 host 环境执行。

## 31. Remaining Unsupported capabilities

**Unsupported**：constant-density MPI/multi-block、IBM、turbulence/MRF/gravity/custom source、额外物理方程、非正交压力修正、新 convection/face-coupling recipe、隐式时间方法。Conservative SIMPLE/PIMPLE 仍缺 fixed-time predictor；不得用 repeated full-`dt` legacy operation 伪装。**Legacy**：保守变量 pressure corrector 与 Eulerian/IBM 专用数值 kernel 仍在现有路径，但 Plan/operation authority 未回退。

AppleDouble 清理：工作树中 `.git` 以外的 `._*` 已删除；`.gitignore` 原有 `._*` 和 `.DS_Store` 规则。`.git/` 内约 2610 个 AppleDouble 元数据文件仍存在：对该目录的递归删除被自动审批拒绝，本轮没有绕过或触碰 Git 元数据。

## 32. Recommended next phase

可独立选择 **Phase 28A — Distributed Pressure Operators**，或 **Phase 28B — Equation-Term Providers for Pressure-Constrained Momentum**。不要再改 SIMPLE/PISO/PIMPLE 的主架构，也不要把 non-orthogonal、IBM 或模型源项暗中降级。
