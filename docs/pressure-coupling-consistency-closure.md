# Phase 27B — 压力耦合数值一致性闭环

日期：2026-09-26。状态：**Implemented**（限定于 serial、single-block、正交网格、常密度 pressure system）。本文记录实际执行能力；`Unsupported` 不代表接口缺失。

## 1. Scope

本轮收敛松弛后 U/p/face flux 的同一迭代态、PISO/SIMPLE/PIMPLE preset 计数语义、通用 Loop 收敛退出，以及 constant-density pressure provider 的编译配置绑定。属于“改变怎么解方程”；没有增加物理 term、equation 或新的 solver family。未进入 Phase 28。

## 2. Starting state

Phase 26B/27 已有 `RawEquationSystem → pressure formulation → ExecutableEquationSystem → operations → PlanFragment → CompiledSolvePlan → OpRegistry → PlanExecutor`，并有共享压力数值实现。此前 fixed-time 路径的 `convergence.evaluate` 只做诊断，outer 循环总是执行配置次数；PISO 的旧测试还把 outer=2 当作有效拓扑。旧 production pressure 与 conservative execution 仍需原样保持。

## 3. Relaxed U/p/flux inconsistency

压力修正给出候选 `U*、p*、F*` 后，单独松弛 U/p 会留下 `U_relaxed、p_relaxed、F*`。F* 对应的是候选态，不能作为本次 outer 迭代结束的 face-flux authority。此问题只涉及 fixed-time SIMPLE/PIMPLE；PISO 不执行 solution relaxation。

## 4. Final iteration-state contract

`pressure.iteration.begin` 捕获 `U_k、p_k、F_k`；压力修正产生候选；`pressure.relaxation.apply` 写入松弛的 U/p；`pressure.flux.consistency.restore` 恢复同一 U/p 的 face flux；`pressure.convergence.evaluate` 读取最终态。物理状态仍由 `StateBundle` 拥有，`iterationVelocity_、iterationPressure_、iterationFlux_` 只是 stepper workspace，没有第二份 physical-state authority。

## 5. Flux consistency implementation

`PressureOperators::reconstructFaceFlux` 复用既有 `RhieChow::predict`、当前 `rAU_`、几何与 pressure gradient。松弛系数任一小于 1 时，恢复 OpId 根据松弛后的 U/p 重建 `correctedFlux_`；均为 1 时沿用 pressure correction 的 F*，保留 PIMPLE outer=1 与 PISO 的原有离散路径。`beginFixedTimeStep` 也为第一轮建立与 U_n/p_n 对应的 F_k。没有引入另一种 Rhie–Chow 公式或通量 fallback。

## 6. Continuity diagnostic authority

候选 `div(F*)` 只保存为 `candidateDiv` 调查值。fixed-time 的 `maxDiv`、`Pressure continuity maxDiv(before/after)` 与 convergence 判断均来自松弛后最终 `correctedFlux_` 的 `max |div F|`。PISO 的既有 continuity 消息与输出顺序保持不变。

## 7. PISO preset semantics

PISO 要求 `outerCorrectors == 1`；Plan 是一次 momentum predictor 加固定 N 次 pressure correction，不生成 outer fixed-point wrapper。`pressureCorrectors >= 1`，PISO 不使用 outer convergence signal。`explain` 不把 `outerCorrectors` 展示为 PISO 执行循环。

## 8. SIMPLE preset semantics

SIMPLE 要求 `pressureCorrectors == 1`。每个最多 N 次 outer 迭代包含一次 predictor、一次 pressure correction、一次 relaxation、一次 flux consistency restore、一次 convergence evaluation；`outerCorrectors` 是最大迭代数。

## 9. PIMPLE preset semantics

PIMPLE 的 outer 最多 N 次，每个 outer 固定 M 次 pressure correction；只允许 outer 根据 convergence signal 提前退出，inner correction 不提前退出。SIMPLE/PIMPLE 使用相同 OpId 和数值 provider，没有独立 numerical kernel。

## 10. Preset validation

typed config 与 transformation match 均拒绝 `PISO outerCorrectors != 1` 和 `SIMPLE pressureCorrectors != 1`，错误建议选用 PIMPLE；没有 clamp、自动升级或按 solver family 决策。`relativeTolerance` 必须有限且大于 0，`absoluteTolerance` 必须有限且非负。旧 architecture 测试的 PISO outer=2 假设已经修正，并新增 PIMPLE 双层拓扑验证。

## 11. Generic loop termination IR

通用 `Loop` 节点携带 `terminationSignal` 与 `minimumIterations`，fragment 到 compiled plan 逐项降低并验证范围。`PlanExecutor` 每次进入 loop 重置 signal，完整执行一次 loop body 后检查发布值；未发布即检查会 fail-fast。executor 不识别 pressure、preset 或 residual 公式。

## 12. Convergence signal lifecycle

`pressure.step.begin` 和 loop 入口重置 `pressure.outer.converged`；`pressure.iteration.begin` 清除 provider 内的 converged 状态；`pressure.convergence.evaluate` 每个 outer 发布本次结果。每次 `execute` 创建独立 signal state，不能从上一 physical timestep 继承。

## 13. Convergence formula

第一版定义为：`ΔU=max_cells,components |U_(k+1)-U_k|`，`Δp=max_cells |p_(k+1)-p_k|`；`scaleU=max(1,max|U_k|)`，`scaleP=max(1,max|p_k|)`；`maxDelta=max(ΔU/scaleU,Δp/scaleP)`；`maxDiv=max_cells |div(F_(k+1))|`。`converged = (maxDelta <= relativeTolerance) && (maxDiv <= absoluteTolerance)`。`fluxDelta=max_faces |F_(k+1)-F_k|` 是辅助诊断，不是第二套停止准则。非有限值 fail-fast。

## 14. Fixed-time temporal-base preservation

`pressure.step.begin` 冻结 `U_n`；后续 outer 的空间 RHS 使用 `U_k`，但候选仍由 `U_n + dt·RHS_spatial(U_k)` 构成。没有把每次 predictor 变成从 U_k 再走完整 dt 的多次物理推进。first-iterate 的求值顺序保持 Phase 26 PISO degeneration。

## 15. Early exit + step commit semantics

收敛检查位于 `iteration.end` 之前；完整提交当前 relaxed U/p/F iterate 后，generic Loop 可以提前结束。`pressure.step.commit` 位于 outer Loop 之外，每个 timestep 只执行一次。cavity 两步对照中 PISO、SIMPLE、PIMPLE 时钟均为 `(0.000390625, 0.000390625)`、`(0.00078125, 0.000390625)`。

## 16. Runtime binder cleanup

`bindPressureProvider()` 接收已解析的 operation binding、`CompiledNumericalSystem` 和仅供旧 conservative corrector 使用的窄 `LegacyPressureInputs`；不再接收整个 `SolverConfig`，不从 raw model/phase/turbulence 配置重新推断 provider。旧 conservative corrector 的 boundary/pressure/gamma 输入仍标记 **Legacy**，未在本轮改其数值实现。

## 17. Pressure operator compiled config

constant-density provider 读取 `CompiledPressureOperatorConfig`：粘度、U/p 松弛因子、相对/绝对收敛阈值、pressure reference、线性求解配置、typed U/p 边界。`rhoConst` 来自 compiled state realization 的派生常量，未复制成第二份密度 authority。缺 compiled config 时绑定失败。

## 18. Rhie–Chow independence

`RhieChow::predict/correct` 仍只处理 collocated face coupling。其实现不读取 SIMPLE/PISO/PIMPLE、outer 次数、松弛或收敛阈值；恢复操作只是复用当前 face predictor。正交网格以外仍为 **Unsupported**。

## 19. Linear backend independence

压力 provider 把 sparse matrix 和 `LinearSolverConfig` 交给已有 `SolverSession`；线性 backend 不知道 preset 或 outer loop。没有改变 pressure matrix、reference cell 或 pressure difference 语义。

## 20. Numerical tests

`cmake --build build --clean-first --parallel 4` 完成 277 个 Ninja 步骤，无编译/链接失败。`ctest --test-dir build --output-on-failure --parallel 1` 最终 **19/19 通过**；首轮唯一失败为旧 `formulationArchitecture` 测试的 PISO outer=2 假设，更新拓扑断言后复测及全套测试均通过。`tools/check_architecture.py` 通过，`git diff --check` 通过。

定向测试检查 typed invalid preset、PISO/SIMPLE/PIMPLE Plan topology、通用 mock Loop 在两次 physical step 分别提前退出、固定 U_n、一次 clock commit、trace 中 `relaxation → flux restore → convergence`、relaxed candidate/final continuity 差异，以及最终 continuity 消息与最终 flux diagnostic 一致。每次 fixed-time convergence decision 均输出 `converged=true/false`，early-exit 测试要求退出前报告 true。

## 21. SIMPLE Poiseuille result

既有几何/粘度/dt/pressure drop 不变。512-step SIMPLE 变体：`L2=0.00822508807`、`Linf=0.0112953895`、末步 `maxDiv=1.37627e-09`、pressure drop `0.1`、每步 2 次 outer、末次 `maxDelta=9.26491e-06`。每步第二次迭代的 maxDelta 均下降；末步 `converged=false`，原因是该原始变体最多只允许 2 次 outer，且最终 continuity 尚高于默认 `absoluteTolerance=1e-10`。它稳定完成并如实报告未达停止准则，没有伪称收敛。既有 PISO Poiseuille 16×16：`L2=0.000746526192`、`Linf=0.000988599105`；32×32：`L2=0.000354528614`、`Linf=0.00047749875`。SIMPLE 精度与 PISO 不同，本轮未调整格式以逼近 PISO。

## 22. PIMPLE cavity result

`outerCorrectors=4`、`pressureCorrectors=2`、U/p relaxation=0.8、`relativeTolerance=0.004`、`absoluteTolerance=5e-8`：两步各在第 2 次 outer 退出。每步一次 step.begin、两次 iteration.begin、每 outer 两次 pressure.solve、一次 step.commit。SIMPLE 的相同最大次数与阈值也在每步第 2 次 outer 退出。

## 23. PIMPLE→PISO degeneration

PIMPLE outer=1、pressureCorrectors=2、relaxation=1 与 PISO 的最终 `U、p、rho` 每个样本差值均不超过 `1e-12`，physical clock 相同。

## 24. PISO frozen regression

Phase26 独立可执行文件与当前文件对同一 constantDensityPiso case：6 条时钟/closure diagnostics 逐字相等，3 个 VTS 文件 SHA-256 逐一相等。PISO 数值路径未受本轮 fixed-time flux restore 影响。

## 25. Production regressions

Phase26 独立构建与当前构建在相同输入下比较：4-rank Sod 20 步 diagnostics `40/40`、VTS `4/4`；serial Ghost 20 步 `120/120`、VTS `1/1`；Eulerian 2 步 VTS `3/3`；legacy pressureConstraintPiso diagnostics `4/4`、VTS `2/2`。所有可比 diagnostics 和 VTS 哈希均相等。legacy PISO 最终 SHA-256 为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`。CTest 的 WENO7/Rusanov frozen Sod regression 通过。MPI 对照在允许 PRTE 打开本机 socket 的 host 环境完成；初次 sandbox 启动失败属于 runtime socket 权限，非测试断言失败。

## 26. Architecture guards

`tools/check_architecture.py` 检查 PISO 无 outer wrapper、fixed-time schedule 包含 flux restore、binder 不持有 raw `SolverConfig`、generic executor 具有中立 signal contract，且 numerical leaf 中没有新的 SIMPLE/PISO/PIMPLE family dispatch。没有新增 fallback、MPI 直接调用或第二个 pressure lifecycle。`src/ARCHITECTURE.md` 同步记录已实现的控制语义。

## 27. Remaining Unsupported capabilities

**Unsupported**：constant-density MPI、多 block pressure coupling、constant-density IBM、turbulence、MRF、gravity、custom Momentum source、额外方程、非正交 correction mathematics、implicit time、BackwardEuler、SDIRK、VOF 与 multiphase pressure coupling。能力检查继续显式拒绝，不静默忽略或替换数学。

## 28. Architecture freeze statement

在上述限定能力内，SIMPLE/PISO/PIMPLE 核心架构冻结：preset 贡献 pressure-constraint formulation 与 SolvePlan；`CompiledSolvePlan` 拥有顺序和循环；`PlanExecutor` 解释通用控制流；共享 provider 只执行单个 pressure operation。此次对 fixed-time SIMPLE/PIMPLE 的数值语义有明确改变：松弛后重建一致 face flux，且可提前收敛；PISO、旧 conservative、Eulerian、Ghost 和 MPI numerical baselines 保持原样。公共 API 增加的是 generic Loop signal 与独立 flux restore operation，不是新 solver family。

## 29. Recommended Phase 28

建议下一阶段单独处理 distributed pressure operators：U/p halo、canonical face-flux owner、pressure matrix row ownership、distributed HYPRE、global continuity reduction；再独立处理 pressure-constrained Momentum term providers。以上均属于后续工作，本轮未实现，也不要求重新定义三个 preset 的生命周期。
