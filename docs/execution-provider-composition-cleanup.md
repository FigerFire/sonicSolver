# Execution Provider Composition Cleanup Report

日期：2026-09-20

## 1. Scope

状态：Implemented

本轮只迁移 execution composition、provider requirement、runtime service requirement 和失效架构。未修改 Euler/SSPRK3/RK4、WENO/TENO、通量分裂、CFL、PISO/PIMPLE、IBM/KKT、phase/level-set/turbulence 方程或 MPI COPY/SUM 数学。

## 2. Before execution graph

状态：Legacy

```text
CaseConfig
  -> ResolvedSimulationSystem
  -> SF_executionBuilder
       if PressureVelocityCoupling -> Eulerian executor
       else                       -> conservative executor
  -> application checks E_LEVEL_SET / E_PHASE_MASS / E_LEGACY_ALPHA
  -> manually populated SolverServices
  -> OpRegistry
  -> PlanExecutor
```

审计分类如下：

| 判断 | 原类别 | 处理 |
|---|---|---|
| composition 中解析 physics、unknown、equation | A 数学 composition | 保留 |
| multi-patch、HYPRE、IBM topology 能力检查 | B capability validation | 改为读取 compiled requirements |
| equation ID 决定 level-set/mixture/thermodynamics 对象 | C provider selection | 移入 `providerRequirements` |
| `PressureVelocityCoupling` 选择 Eulerian executor | D solver-family routing | 删除 |

## 3. After execution graph

状态：Implemented

```text
CaseConfig
  -> SystemBuilder
  -> ExecutableEquationSystem + CompiledSolvePlan
  -> derived providerRequirements/runtimeServiceRequirements
  -> concrete state/workspace/provider binders
  -> OpRegistry
  -> PlanExecutor::validateBindings
  -> PlanExecutor
```

`SF_executionBuilder.cpp` 只匹配 `flow.conservative` 或 `flow.eulerian-pressure` provider group；它不再读取 equation ID、physics template、formulation enum 或 solve strategy。

## 4. Remaining runtime authority

状态：Implemented

- 数学对象：`ExecutableEquationSystem`。
- 控制流：`CompiledSolvePlan`。
- operation 集合：从 Plan 递归得到的 `RuntimeReport::requiredOperations`。
- 数值实现：`CompressibleAlgorithm`、`PressureStepper` 及 IBM/phase/interface concrete providers 注册的 callbacks。
- 并行语义：`IExecutionRuntime` / `IParallelCoordinator` 的 backend implementation。

## 5. ExecutionRequirements design

状态：Implemented

没有新增通用 Provider 基类。`ResolvedSimulationSystem` 增加两个简单 value-object 列表：

- `ProviderRequirement { id, responsibility }`
- `RuntimeServiceRequirement { id, reason }`

其余需求不复制：operation 仍在 `runtime.requiredOperations`，unknown/state、constraint、closure、compiled equation binding 仍在 `ExecutableEquationSystem`，workspace 仍在 `workspaceRequirements`。`deriveExecutionComposition()` 根据 executable equations、closures、Plan operation IDs、已有 capability requirements 与 IBM descriptor 推导 provider/service keys。

## 6. Provider composition

状态：Implemented

当前 provider keys 包括 `flow.conservative`、`flow.eulerian-pressure`、`thermodynamics.single-fluid`、`thermodynamics.homogeneous`、`equation.legacy-mixture`、`equation.level-set`、`equation.phase-change`、`closure.turbulence`、`equation.turbulence-transport`、`ibm.boundary`、`ibm.constraint` 和 `linear.hypre`。未请求的 key 不创建对应 model/workspace。

Validator 要求每个 resolved system 恰好有一个 flow provider group，并拒绝空、重复 provider/service requirement。`PlanExecutor::validateBindings()` 继续在首个 operation 前验证所有 Plan leaf。

## 7. single-fluid path

状态：Implemented

single-fluid explicit 与 PISO 都选择 `flow.conservative`。PISO 额外要求 `linear.hypre` 和 `mpi.initialized`，不会因为 `PressureVelocityCoupling` 被送入 Eulerian executor。SIMPLE/PIMPLE fixed-point Plan 仍可形成，但因 dedicated predictor provider 未实现而在执行前报告 Unsupported。

## 8. Eulerian path

状态：Implemented

只有 Plan 出现 `ee.*` operation group 时才要求 `flow.eulerian-pressure`。`PressureStepper` 不决定 equations 或三层循环，只持有真实 N-phase state、pressure/phase workspace，并把现有数值 helper 绑定为 `ee.*` callbacks。2 相改为 3 相的 architecture test 无需修改 runner、PlanExecutor 或 application branch。

## 9. homogeneous multiphase path

状态：Implemented

`E_PHASE_MASS` 产生 `thermodynamics.homogeneous` requirement。原 `HomogeneousPhaseChangeCoupling` 更名为 `HomogeneousPhaseChangeProvider`；它保留真实 phase-change RHS/commit 数值职责。原 Legacy mixture concrete types 更名为 `MixtureEquationProvider` 与 `MultiPatchMixtureEquationProvider`，其辅助状态推进、RHS、halo 和 commit 行为未变。

## 10. level-set path

状态：Implemented

`E_LEVEL_SET` 产生 `equation.level-set` requirement；移除该 equation 后不会创建 interface model/provider。`InterfaceEquationProvider` 与 multi-patch 版本拥有真实 phi transport、reinitialization、curvature publish、jump/transport port 和 halo contract，因此不是 forwarding wrapper。

## 11. turbulence/source path

状态：Implemented

closure metadata 产生 `closure.turbulence`；只有 `E_TURB_*` 存在时才产生 `equation.turbulence-transport`。DNS closure 与 transported RAS equations 不再共用一个模糊开关。MRF、gravity 等仍作为 equation source contribution，由现有 assembler 计算，没有新增 executor。

## 12. IBM path

状态：Implemented

Ghost/ILW descriptor 产生 `ibm.boundary`。`ibm.constraint.project` 或 `ibm.kkt.solve` 产生 `ibm.constraint`。application 会核对 compiled IBM requirement 与实际已初始化的 IBM resource，不匹配即 fail fast；不再以 nullable pointer 静默跳过。MPI-2 forcing 与 surface KKT 均继续走原 constraint/HYPRE implementation。

## 13. StateBundle/resource ownership

状态：Implemented

`StateBundle` 仍只拥有 persistent physical/distributed state 和 clock。`realizeState()` 从 executable unknown storage bindings 与 workspace requirements 验证每个 participating patch 的资源覆盖。RK、pressure、phase 和 KKT scratch 继续由其 concrete provider/stepper 生命周期持有，没有加入 universal workspace 或第二份 physical state。

## 14. SolverServices cleanup

状态：Implemented

删除从未被 production/test 设置的 `SolverServices::boundary`。density RHS 在没有 `BoundaryPipeline` 的直接测试路径使用自身稳定构造的 `Boundary::Applicator`。剩余字段均有 production reader：boundary pipeline、execution runtime、transport contribution、equation contribution、IBM 三端口和 observer；绑定仍只允许一次，timestep 中不切换。

## 15. Deleted interfaces

状态：Implemented

本轮删除四个位于 interface umbrella、全仓零调用的 concrete no-op 类型：`NoopImmersedBoundary`、`NoopBoundaryPipeline`、`LocalParallelCoordinator`、`NoopTransportModel`。它们不是 interface contract，也没有 production/test caller。

保留的 core interfaces 都有独立跨层语义：`INavierStokesStepper` 有 conservative 与 Eulerian 两个实现；boundary、execution runtime、parallel、transport、equation coupling、observer、IBM boundary/constraint/system ports 各自隔离不同 lifecycle 或 backend。

## 16. Deleted adapters/wrappers

状态：Implemented

删除旧命名 `LegacyMultiphaseEquationCoupling`、`MultiPatchLegacyEquationCoupling`、`registerLegacyState` 和 `InterfaceEquationCoupling`；具体实现以 provider/state registration 职责命名，不保留 alias 或 forwarding compatibility symbol。

`IBBoundaryAdapter`、`IBConstraintAdapter`、`MultiPatchIBAdapter` 仍由 `SF_singleFluid.cpp` / `SF_multiPatch.cpp` 调用。它们把 concrete `IBM::IB`/`CompositeIB` 分别绑定到 boundary、constraint/system、multi-patch boundary 三个窄 port，并保持不同 apply/constraint/KKT 生命周期，不能由一个现有 authority 直接替代。

## 17. Deleted factories

状态：Implemented

未发现已经失去多态数值职责的 execution factory，因此本轮删除数为 0。`makeInterfaceModel`、phase-system、EOS 和 turbulence factories 仍选择真正不同的 numerical model implementation，不属于 solver-family lifecycle routing。

## 18. Deleted fields/enums

状态：Implemented

删除 `SolverServices::boundary` 以及零调用查询 API `hasSolveBlock()`、`hasSolveStrategy()`、`hasRequirement()`。本轮没有删除 enum value；`PhysicsTemplateKind`、formulation 与 strategy enums 仍用于 mathematical composition、capability validation 或 explain provenance，不再用于 execution builder lifecycle selection。

## 19. Deleted files

状态：Implemented

本轮新增需求类型位于现有 resolved-system 文件，provider rename 位于现有 coupling 文件，因此没有新增或删除 source file，也没有 CMake entry 变化。此前迁移已删除的 flow adapter/workflow 文件继续由 architecture checker 的 obsolete-path guard 保护。

## 20. Remaining specialized executors

状态：Legacy

仍保留 `CompressibleAlgorithm`、`EulerianEulerian::PressureStepper` 与 topology-specific `runMultiPatch()`。它们不再选择数学系统或 control-flow；其职责限于 concrete state/workspace ownership、operation callback binding 和现有 numerical helper 调用。

## 21. Why each remaining executor still exists

状态：Legacy

- `CompressibleAlgorithm`：拥有 per-patch `FluxField`/`Residual`、explicit stage storage、generic PISO corrector/KKT workspace，并绑定 `flow.*`、`explicit.*`、`pressure.*`、`ibm.*` operations。
- `PressureStepper`：拥有 N-phase state、interphase/source/turbulence/pressure workspaces、HYPRE reuse state，并绑定 `ee.*` operations。
- `runMultiPatch()`：绑定 mesh partition、canonical face storage、distributed registry 和 multi-patch boundary/interface providers；这是 storage topology，不是 solver-family authority。

## 22. Architecture guards

状态：Implemented

`tools/check_architecture.py` 新增守卫：

- application execution 禁止 `System::hasEquation()` 和 `System::hasSolveStrategy()` provider routing；
- execution builder 禁止 `PressureVelocityCoupling`、`SolverAlgorithm::`、`PhysicsTemplateKind`；
- environment 禁止由 `solverConfig.numerics.solver` 推断 runtime service；
- 禁止重新加入四个 zero-caller no-op types 和已删除 legacy provider symbols；
- 要求 resolved system 存在 provider/runtime-service requirements 且由 builder 推导。

依赖检查保持 `allowlisted dependency edges = 10`，无新增 dependency violation。

## 23. Numerical invariants

状态：Implemented

stage 系数/时间、flux/reconstruction、CFL、pressure matrix/corrector、source clear、IBM projection/KKT、canonical face COPY、GlobalDof SUM、state/geometry owner COPY 均未修改。变更只发生在 requirement 推导、resource/provider 选择与死符号删除。

## 24. Regression results

状态：Implemented

| 检查 | 结果 |
|---|---|
| `cmake --build build --parallel 4` | 完整成功 |
| `cmake --build build-tests --parallel 4` | 完整成功 |
| architecture checker | 成功；10 条 allowlist |
| `pisoArchitecture` | 成功，含 provider、3-phase、missing-provider tests |
| `explicitStageMathematics` | 成功 |
| frozen PISO numerical regression | 成功；最终 VTS SHA-256 `11cba96d719f1b00ced3a01237ba489c636ffec13de3c0c516144af47ca88040`，byte-identical |
| WENO7 Euler | 3-step smoke 成功；完整既有 case 在 `t=0.1297871` 仍因 negative pressure fail-fast，未修改 tolerance/flux |
| Ghost IBM RK4 | 3-step smoke 成功 |
| Velocity Forcing | 3-step smoke 成功 |
| fractional DFM | 3-step smoke 成功 |
| MPI-2 Velocity Forcing | host 环境 3-step smoke 成功 |
| surface KKT | 原 case 3 steps 成功 |
| `numericalFluxContract` | 既有失败保持：constant stencil reconstructed Rusanov 与 physical flux 不一致 |

CTest 汇总：4 项中 3 项成功，唯一失败为已知 `numericalFluxContract`。

## 25. Remaining Legacy

状态：Legacy

- `SolverServices` 仍是向两个 steppers 注入稳定窄端口的聚合 value object；所有剩余字段均有实际 production reader。
- level-set、mixture、turbulence concrete resource construction 仍位于 conservative execution binder 内，但创建条件完全来自 provider requirements；它们不拥有 Plan 或 timestep loop。
- IBM adapters 仍承担 concrete model 到 solver port 的必要 impedance mapping。
- `CompressibleAlgorithm` 名称仍带历史 formulation 色彩，但其保留理由是 workspace/provider ownership，未在本轮仅为命名制造 churn。

## 26. Remaining Unsupported

状态：Unsupported

- single-fluid SIMPLE/PIMPLE dedicated fixed-point predictor provider；
- constant-density pressure numerical provider；
- generic pressure schedule 与 auxiliary level-set/turbulence 的联合 provider；
- homogeneous thermodynamics 的 multi-patch execution；
- multi-patch transported turbulence 与 legacy alpha canonical scalar flux；
- Eulerian variational IBM phase-wise fluid-port assembly；
- multi-patch variational forcing routing。

这些路径继续显式 fail fast，没有 fallback 到其他数学算法。

## 27. Global dead-symbol search

状态：Implemented

production/test source 搜索确认以下符号为 0：`NoopImmersedBoundary`、`NoopBoundaryPipeline`、`LocalParallelCoordinator`、`NoopTransportModel`、`LegacyMultiphaseEquationCoupling`、`MultiPatchLegacyEquationCoupling`、`registerLegacyState`、`InterfaceEquationCoupling`（这些字符串只保留在 architecture guard 中）。application execution 中 `System::hasEquation(` 和 `System::hasSolveStrategy(` 也为 0。

保留的 `Adapter`/`Manager` 名称逐项有 caller：三个 IBM adapters 由 single/multi-patch binder 使用；`Turbulence::Manager` 拥有 runtime-selectable turbulence model/state，并非空 facade。未发现可安全删除且不改变 numerical ownership 的 zero-caller production class。

## 28. Final authority table

状态：Implemented

| 问题 | 唯一 authority |
|---|---|
| 谁决定有哪些 equations？ | `ResolvedSimulationSystem / ExecutableEquationSystem` |
| 谁决定 timestep 做什么？ | `CompiledSolvePlan` |
| 谁决定 numerical operation 怎么算？ | 对应 concrete provider callback |
| 谁选择 provider？ | executable content + required operations 推导的 execution requirements |
| application 是否按 single-fluid/Eulerian 选择完整 lifecycle？ | 否；只绑定 compiled provider group 的资源 |
| 加入第三相是否需要改主 run loop？ | 原则上否；architecture test 已覆盖 2→3 phase |
| 加入 Level Set 是否需要新 solver？ | 否；增加 equation/state/provider requirement |
| 加入 MRF 是否需要新 solver？ | 否；增加 source contribution |
| 加入 IBM constraint 是否切换 solver family？ | 否；增加 constraint 与 `ibm.*` operation |
| missing provider 怎么处理？ | execution 前 fail fast |
| 旧 adapter 是否作为 fallback？ | 否 |

## 29. Recommended next task

状态：Interface-only

下一步应把 conservative binder 中 level-set/mixture/turbulence 的 concrete construction 分拆为同目录的小型 resource-binding functions，并让每个 function 接受单个 provider requirement；同时保持 `StateBundle` 与 provider workspace 边界。不要引入通用 provider interface、DI container 或新的 solver family。数值工作应单独处理已知 WENO7 长程失稳与 `numericalFluxContract`，不能混入本架构提交。
