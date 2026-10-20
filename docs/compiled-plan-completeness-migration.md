# Compiled Plan 完整性与 Built-in IBM Runtime 迁移报告

日期：2026-09-18

## 1. Scope

本轮只迁移 runtime authority，不修改控制方程、离散公式、压力修正、IBM forcing、KKT 矩阵、RK 系数或 MPI ownership 语义。目标主链为：

```text
RawEquationSystem
  -> TransformationPipeline
  -> ExecutableEquationSystem
  -> SolvePlanner
  -> CompiledSolvePlan
  -> Run::PlanExecutor
  -> Run::OpRegistry
  -> existing numerical helpers
```

本轮不迁移 density RK、BindingPass 和 AssemblyPlanRegistry。

## 2. Source audit

审计确认修改前仍有以下 authority 分裂：

- `src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp` 的 `advance()` 仍直接计算 `dt` 并提交 clock。
- `ee.pressure.solve` 同时求解并发布 `p'`，而 `ee.pressure.publish` 是空 callback。
- Eulerian turbulence 同时以 `S_TURBULENCE` policy 和 PressureStepper runtime condition 表示。
- `compileEulerianPimple()` 的特殊分支可以提前返回，从而不检查其他 active policy。
- `RuntimeReport.requiredOperations` 有独立于 Plan 的浅层遍历。
- `SolveStage` / `SolveStageKind` 仍残留在 core interface。
- single-fluid forcing IBM 的数学 descriptor 已存在，但实际 correction/KKT 调用没有由 IBM Plan Op 定位。
- Eulerian-Eulerian variational IBM 缺少 phase-wise fluid-port assembly，却存在被静默遗漏的风险。

## 3. Before architecture

```text
PressureStepper::advance
  -> stableTimeStep
  -> PlanExecutor(inner PIMPLE loops)
  -> state clock commit

CompressibleAlgorithm
  -> density/pressure lifecycle
  -> flowAlgorithm->correct
  -> built-in IBM formula dispatch
```

Plan 已能描述三层 pressure loop，但尚未描述完整 pressure timestep，也没有强制证明每个 active policy 已被消费。

## 4. After architecture

Pressure production path 现在为：

```text
ExecutableEquationSystem
  -> CompiledSolvePlan
  -> PlanExecutor
  -> OpRegistry
  -> PressureStepper numerical callbacks
  -> StateBundle / ExecutionRuntime / HYPRE
```

single-fluid IBM 的 execution requirement 现在由同一份 `CompiledSolvePlan` 给出：

```text
ConstraintProjection -> ibm.constraint.project
MonolithicKKT        -> ibm.kkt.solve
```

Density RK 的 stage loop 仍由 `DensityBasedTime` 执行，状态明确为 `Legacy`；本轮仅用 `executeOperation()` 执行已经编译出的 IBM leaf，没有迁移 RK lifecycle。

## 5. Dead authority removed

已删除：

- `SolveStage`
- `SolveStageKind`
- core interface 中基于 SolveStage 的 binding contract
- `LegacyEulerianExecutionAdapter` fallback
- runtime 对 Workflow-era stage 的依赖

`PressureStepper::stepImpl()` 在上一轮已经删除；本轮继续删除其在 `advance()` 中残留的 dt/clock schedule authority。

## 6. requiredOperations authority

唯一来源现在是：

```cpp
SolvePlanner::requiredOperations(result.solvePlan)
```

该函数递归访问所有 structured node，收集非空 `OpId`。`SystemBuilder` 不再维护浅层遍历或 operation 名称副本。

## 7. Pressure timestep lifecycle before/after

Before：

```text
PressureStepper::advance
  dt compute
  PlanExecutor
  time/step commit
```

After：

```text
PressureStepper::advance
  bind active state
  register numerical callbacks
  PlanExecutor::execute(compiledPlan)
  return snapshot
```

`stableTimeStep()` 数学与副作用未改，只作为 `ee.dt.compute` callback。`state.dt`、`state.time += dt`、`state.step++` 的原顺序未改，只作为 `ee.time.commit` callback，由 Plan 决定执行位置。

## 8. EE Plan structure

实际 compiled tree 为：

```text
EE.step
  ee.dt.compute
  ee.step.begin
  EE.outer x outerCorrectors
    ee.interphase.compute
    ee.sources.assemble
    ee.sources.validate
    ee.turbulence.prepare          [closure/equation 存在时]
    ee.momentum.diagonal
    ee.momentum.flux
    ee.faceFlux.canonical
    ee.continuity.assemble
    ee.boundary.prepare
    ee.momentum.solve
    ee.interphase.correct
    ee.boundary.afterMomentum
    ee.diagonal.sync
    ee.momentum.flux.after
    ee.faceFlux.canonical.after
    EE.pressure x pressureCorrectors
      EE.nonOrthogonal x (nonOrthogonalCorrectors + 1)
        ee.pressure.solve
        ee.pressure.publish
        ee.pressure.sync
        ee.phase.correct
        ee.faceFlux.correct
        ee.boundary.afterPressure
        ee.faceFlux.canonical.pressure
    ee.energy.solve
    ee.turbulence.solve            [transport equations 存在时]
    ee.boundary.final
    ee.outer.validate
  ee.step.commit
  ee.time.commit
```

该顺序来自迁移前真实 `PressureStepper`，没有按教科书顺序重排。

## 9. Pressure operation truthfulness fixes

`ee.pressure.solve` 现在只执行 pressure correction solve，并保存 `LinearAlgebra::SolveResult`。

`ee.pressure.publish` 读取 pending result，发布 pressure correction 和原有 diagnostics，然后清除 pending result。缺少 pending result 或重复发布均 fail-fast。

`ee.pressure.sync` 继续只负责原有同步。不存在空 callback 或把 publish 工作隐藏在 solve 中的情况。

## 10. Turbulence contribution lowering

Eulerian turbulence transport equations 现在扩展 `S_EE_PIMPLE`。Plan 是否包含 turbulence operation 由 `ExecutableEquationSystem` 中的 closure/equation 决定：

- closure 存在：生成 `ee.turbulence.prepare`；
- `E_TURB_*` transport equation 存在：生成 `ee.turbulence.solve`；
- equation 不存在：不生成 solve Op。

已删除 runtime 的 `turbulenceStageActive_` 二次判断。turbulence 模型仍提供现有数值 helper，不拥有第二个 global lifecycle。

## 11. Policy consumption rules

SolvePlanner 维护本地 consumed policy ID 集合。规划结束时，所有 active policy 必须满足：

```text
consumed by compiled plan
or
fail-fast with policy ID
```

Eulerian lowering只接受 `S_EE_PIMPLE` 以及被合并进该 schedule 的 transitional `S_TURBULENCE` metadata。任何其他 active policy 都会报告未消费，不再因 EE special branch 提前返回而丢失。

## 12. Unsupported contribution validation

`SystemBuilder` 在 runtime execution 前拒绝以下组合：

- Eulerian-Eulerian + Ghost IBM：缺少 phase-wise ghost-state boundary closure；
- Eulerian-Eulerian + Peskin/DFM/Velocity Forcing/KKT：缺少 phase-wise IBM fluid-port assembly；
- distributed monolithic single-fluid KKT：现有 distributed pressure/KKT contract 尚未完成。

诊断同时说明 IBM constraint、Eulerian multiphase system 和缺失 capability。不存在运行 Eulerian 而忽略 IBM、切换到 single-fluid IBM 或 Ghost fallback 的路径。

## 13. IBM mathematical classification

| 数学类别 | built-in method | runtime 表达 | Status |
| --- | --- | --- | --- |
| Boundary Closure | Ghost / ILW | boundary composition 与原 boundary order | Implemented |
| Sequential Constraint | Peskin original | `ibm.constraint.project` | Implemented |
| Sequential Constraint | Velocity Forcing FTS | `ibm.constraint.project` | Implemented |
| Sequential Constraint | Brinkman/BP | `ibm.constraint.project` | Implemented |
| Sequential Constraint | explicit/fractional DFM | `ibm.constraint.project` | Implemented |
| Coupled Constraint Block | implicit DLM | `ibm.kkt.solve` | Implemented |
| Coupled Constraint Block | augmented-Lagrangian KKT | `ibm.kkt.solve` | Implemented |
| Eulerian phase-wise variational IBM | all forcing/KKT variants | capability rejection | Unsupported |

## 14. Single-fluid IBM Plan lowering

`CompressibleAlgorithm::bindSolvePlan()` 从 recursively derived required operations 中绑定唯一 IBM OpId。执行 correction 时：

```text
CompiledSolvePlan leaf
  -> PlanExecutor::executeOperation
  -> OpRegistry callback
  -> existing IFlowAlgorithm::correct
  -> existing projectPredictedState/KKT helper
```

`executeOperation()` 要求 Plan 中恰好存在一个同名 leaf；不存在或重复均 fail-fast。它是 density StageLoop 迁移前的过渡入口，不创建第二份 schedule。

## 15. Ghost IBM behavior

Ghost/ILW 仍是 boundary/interface closure，没有被转换为 constraint source。真实顺序仍为：

```text
physical boundary -> halo -> ghost/ILW reconstruction -> spatial operator
```

本轮没有改 Ghost stencil、ILW order、低阶处理、geometry classification 或 boundary/halo 顺序。

## 16. Sequential IBM behavior

Plan 使用一个与现有数值责任相符的 leaf：

```text
ibm.constraint.project -> projectPredictedState()
```

内部 built-in switch 仅选择 Peskin、DFM、FTS 或 BP 的局部数值公式，不选择 timestep、PIMPLE 或 RK lifecycle。没有新增 IBM manager、registry 或 plugin framework。

## 17. KKT behavior

当前 KKT helper 紧密组合 `prepareMonolithicSystem()`、HYPRE solve 和 `acceptMonolithicSolution()`。为避免改变 matrix assembly 和 correction 时机，本轮以一个 truthful leaf 表示：

```text
ibm.kkt.solve
```

serial `dfmImplicitPrescribed` 与 `dfmAugmentedLagrangian` 实际通过该 Op 执行。KKT block、row ownership、preconditioner 和 HYPRE 参数未改。

## 18. Eulerian + IBM unsupported combinations

Architecture test 构造 Eulerian variational IBM contribution，确认 build 阶段返回明确 Unsupported reason。Ghost 和 variational forcing 分开诊断，没有把 Ghost 当作 forcing/KKT capability。

## 19. Runtime operation table

| OpId | callback responsibility | owner of order | Status |
| --- | --- | --- | --- |
| `ee.dt.compute` | 原 `stableTimeStep()` | CompiledSolvePlan | Implemented |
| `ee.step.begin` | 原 timestep begin helper | CompiledSolvePlan | Implemented |
| `ee.pressure.solve` | pressure solve，保留 pending result | CompiledSolvePlan | Implemented |
| `ee.pressure.publish` | publish `p'` 与 diagnostics | CompiledSolvePlan | Implemented |
| `ee.time.commit` | dt/time/step commit | CompiledSolvePlan | Implemented |
| `ibm.constraint.project` | existing sequential IBM projection | CompiledSolvePlan | Implemented |
| `ibm.kkt.solve` | existing coupled IBM/KKT path | CompiledSolvePlan | Implemented |
| density RK stage ops | `DensityBasedTime` internal stages | DensityBasedTime | Legacy |

## 20. Numerical invariants preserved

未修改：

- WENO/TENO/Rusanov/Steger-Warming 公式；
- RK/Euler coefficients 和 stage time；
- pressure matrix、predictor、velocity/flux correction；
- Peskin `J/J^T`、DFM、FTS、BP、KKT 公式；
- boundary/halo ordering；
- canonical face `COPY`；
- residual/source/load/GlobalDof `SUM`；
- state owner-to-replica `COPY`；
- constraint/multiplier ownership；
- force、torque、power reduction。

Serial runtime 在 MPI backend 可用但 world size 为 1 时继续使用 serial execution ownership；rank-1 reduction 是 identity。

## 21. Architecture tests

通过：

- `cmake --build build --parallel 4`：完整 build 成功，`sonicSolver` 与 GUI targets 均完成 compile/link；
- `python3 tools/check_architecture.py`：608 files，allowlisted dependency edges = 10；
- `test_pisoArchitecture`：Plan recursion、single-leaf execution、EE timestep ops、turbulence presence、unconsumed policy、Eulerian IBM Unsupported 均通过；
- guard 确认 PlanExecutor 不检查 Eulerian/IBM/turbulence/formulation/solver family；
- guard 确认 `PressureStepper::advance()` 不含 dt calculation 或 clock commit；
- guard 确认 System/Transformer 不直接访问 MPI backend。

当前 `build-phase1a-tests` 重新生成后包含 3 项测试：`pisoArchitecture` 和
`pisoNumericalRegression` 通过；`numericalFluxContract` 因第 22 节记录的既有
high-order constant-stencil 数值差异失败。该失败不是 compile/link、Plan、PISO
或 IBM assertion failure。

## 22. Serial regression

### Frozen PISO

`test/pressureConstraintPiso` 通过，final VTS SHA-256 与 frozen baseline 完全相同：

```text
11cba96d719f1b00ced3a01237ba489c636ffec13de3c0c516144af47ca88040
```

因此所有输出数组的 difference `L∞ = 0`、`L2 = 0`。关键结果：

```text
dt/time       = 1.659315e-06
iterations    = 2
p residual    = 4.38769e-11
maxDiv        = 2.01179 -> 2.01052
HYPRE         = 1 rebuild / 1 solve
```

### Eulerian pressure

新 Plan 实际执行一 timestep：

```text
dt/time              = 3.6242994983723128e-05
outer/pCorr/pIter     = 1 / 1 / 2
pressure residual    = 4.645606e-14
maxAlphaSumError     = 0
phaseMass            = 2.358562e+02
phaseEnthalpy        = 2.958757e+08
HYPRE                 = 1 rebuild / 1 solve
```

迁移前在同一工作树中该 case 在 timestep 前以 `Invalid drag/virtual-mass coefficient` 失败，因此没有可计算的 before/after field `L∞/L2`。本报告不把该 case 宣称为 frozen numerical equivalence；Plan trace 证明执行顺序与旧 source lifecycle 逐项一致。

### High-order contract

`numericalFluxContract` 仍报告独立的既有 numerical issue：constant stencil 经 characteristic reconstruction 后的 Rusanov flux 与 physical flux 不一致，例如 mass flux `36.0551` 对 `36`。本轮没有修改 WENO、characteristic matrix 或 Rusanov 数学，也没有调整 tolerance。

## 23. MPI regression

host MPI/PRTE 环境实际运行：

- 4-rank Sod：通过，`dt/time = 6.681531e-04`，TENO5 + Steger-Warming 路径保持；
- Velocity Forcing FTS MPI-2：通过一步，两个 rank 均执行 `ibm.constraint.project`；
- Velocity Forcing FTS MPI-4：通过一步，四个 rank 均执行 `ibm.constraint.project`；
- fractional DFM MPI-2：通过一步；
- fractional DFM MPI-4：通过一步；
- Peskin MPI-2：两个 rank 均进入 `ibm.constraint.project` 后停在现有 Peskin distributed helper；
- distributed monolithic KKT：capability validation 明确 Unsupported。

MPI-4 Velocity Forcing 关键结果：

```text
surface dofs       = 12
max|Ju-Us|         = 1.96005e-12
SchurCG iterations = 9
dt/time            = 1.0e-04 (one-step diagnostic copy)
```

MPI-4 fractional DFM 关键结果：

```text
body dofs          = 138
max|Ju-Us|         = 0
dt/time            = 1.0e-04 (one-step diagnostic copy)
```

上述 MPI IBM 运行使用 case copy 将 endTime 缩短为一步，仅用于 runtime/order/ownership diagnostic；按 IBM validation contract，它们不是 5 s final-field acceptance。canonical case 的 numerics、partition 和源码没有被修改。

## 24. IBM regression

所有 repository `cylinderFlow*` serial built-in family 均完成一步 runtime diagnostic：

| Case | YAML algorithm | 关键结果 | Status |
| --- | --- | --- | --- |
| cylinderFlowGhost | Ghost | RK4 final min rho=1.22383, min p=101208 | Implemented |
| cylinderFlowPeskin | `peskinOriginal` | markers=12, pre-next-step max constraint=2.36643 | Implemented |
| cylinderFlowVelocityForcing | `velocityForcingFTS` | max constraint=1.9613e-12, SchurCG=9 | Implemented |
| cylinderFlowVelocityForcingBP | `velocityForcingBP` | body cells=138, coefficient=1000 | Implemented |
| cylinderFlowFictitiousDomain | `dfmFractionalStepPrescribed` | body dofs=138, max constraint=0 | Implemented |
| cylinderFlowDFMExplicitSelfPropelled | `dfmExplicitSelfPropelled` | body cells=138, lagged multiplier advanced | Implemented |
| cylinderFlowDFMFractionalStepSelfPropelled | `dfmFractionalStepSelfPropelled` | body dofs=138, max constraint=0 | Implemented |
| cylinderFlowDFMImplicitSelfPropelled | `dfmImplicitSelfPropelled` | KKT iter=3, residual=1.74224e-14 | Implemented |
| cylinderFlowDFMAugmentedLagrangian | `dfmAugmentedLagrangian` | KKT iter=3, residual=1.2545e-12 | Implemented |
| cylinderFlowFullyImplicitDLM | `dfmImplicitPrescribed` | KKT iter=2, residual=1.07973e-09 | Implemented |
| Peskin MPI-2/MPI-4 | `peskinOriginal` | distributed helper does not complete | Unsupported |
| monolithic KKT MPI | implicit/augmented | distributed contract unavailable | Unsupported |

Peskin 与 Ghost 的迁移前/后日志中 dt、stage extrema 和 method diagnostics 文本相同；这些 logged checkpoint 的 difference 为 0。未保存独立的迁移前完整 VTS，因此不虚构 IBM full-field `L∞/L2`。

## 25. First-difference analysis

本轮发现三类差异：

1. Brinkman/BP 起初要求 multiplier transformer，但该方法实际没有 multiplier constraint。修正 mathematical lowering 后，它生成 `ibm.constraint.project`；数值公式未改。
2. fully implicit DLM 起初因 pressure entry 只接受 predictor + pressure correction 而在执行前失败。binding 现在接受 compiled `MonolithicKKT`，实际进入原 KKT helper；数值矩阵未改。
3. Peskin MPI-2/4 的所有 rank 都到达同一个 `ibm.constraint.project`，随后停在 existing distributed Peskin helper。首差异位于 leaf 内部，不是 Plan order、loop count 或 callback omission。本轮按 numerical freeze 不修改该算法。

## 26. Deleted files/types

本轮删除死类型 `SolveStage` / `SolveStageKind`。此前已经失去 authority 的 Workflow source files 继续保持删除状态；数值 helper 没有因 wrapper 删除而删除。

## 27. Remaining Legacy

| Concern | Current authority | Status |
| --- | --- | --- |
| density RK lifecycle | `DensityBasedTime` / `CompressibleAlgorithm::stepDensity` | Legacy |
| density IBM leaf insertion | `PlanExecutor::executeOperation` inside legacy stage tail | Legacy |
| equation term lowering | `AssemblyPlanRegistry` / EquationId dispatch | Legacy |
| binding pass | transformer-time binding | Interface-only |
| distributed Peskin execution | existing distributed helper | Unsupported |
| distributed monolithic KKT | no complete execution contract | Unsupported |
| Eulerian variational IBM | no phase-wise fluid-port assembly | Unsupported |

## 28. Next task

下一任务是 `Density RK Plan Authority Migration`：把 `DensityBasedTime::ddtDispatch` 的 StageLoop 降到 `CompiledSolvePlan -> PlanExecutor -> OpRegistry`。之后再处理 BindingPass 和 `Term -> OperatorId`，本轮未实现这些内容。

## Final authority table

| Concern | Before | After | Status |
| --- | --- | --- | --- |
| global runtime flow | mixed algorithm wrappers | CompiledSolvePlan（density RK 例外单列） | Implemented |
| pressure timestep | PressureStepper + Plan | CompiledSolvePlan | Implemented |
| outer loop | Plan | Plan | Implemented |
| pressure loop | Plan | Plan | Implemented |
| nonOrth loop | Plan | Plan | Implemented |
| operation requirements | duplicated/manual | recursive Plan-derived | Implemented |
| pressure publish | hidden in solve/no-op node | explicit truthful Op | Implemented |
| turbulence scheduling | policy + hard-coded runtime | system/plan-derived | Implemented |
| active policy handling | possible silent ignore | consume or fail-fast | Implemented |
| Ghost IBM | boundary path | unchanged boundary contribution | Implemented |
| sequential single-fluid IBM | legacy correction path | Plan Op + existing helper | Implemented |
| monolithic single-fluid IBM | legacy KKT path | Plan Op + existing helper | Implemented |
| Eulerian variational IBM | silently at risk of omission | explicit capability rejection | Unsupported |
| state | StateBundle | StateBundle | Implemented |
| MPI | ExecutionRuntime | ExecutionRuntime | Implemented |
| density RK | DensityBasedTime | DensityBasedTime | Legacy |
| BindingPass | transformer binding | transformer binding | Interface-only |
| equation lowering | AssemblyPlanRegistry | AssemblyPlanRegistry | Legacy |
