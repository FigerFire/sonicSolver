# Flow Algorithm Adapter Removal Report

日期：2026-09-19

## 1. Scope

本轮只清理 single-fluid runtime 的 control authority。它属于“改变怎么解 equation”的执行架构迁移，不改变 equation、离散公式、状态存储或 MPI 语义。

状态：**Implemented**。

## 2. Before runtime graph

迁移前同时存在两条 authority：

```text
ResolvedSimulationSystem -> CompiledSolvePlan
CompressibleAlgorithm
  -> density/pressure booleans
  -> RuntimeStatus routing
  -> stepDensity / stepPressureGeneric / stepPressure
  -> IFlowAlgorithm::correct
  -> DensityBased::Algorithm / PressureBased::Algorithm
```

`stepPressure()` 自行控制边界、dt、时间推进、IBM correction、commit；`flowAlgorithm_->correct()` 又按 solver family 选择 IBM correction。这使 Plan 不能单独回答一个 timestep 执行哪些操作。

状态：**Legacy**（已删除）。

## 3. After runtime graph

```text
ResolvedSimulationSystem
  -> CompiledSolvePlan
  -> required OpId
  -> OpRegistry numerical callback
  -> PlanExecutor
  -> StateBundle
```

`CompressibleAlgorithm::advance()` 只绑定当前 Plan 所需的 explicit、pressure、IBM、commit 操作，先验证 operation coverage，再执行同一份 Plan。它不再根据 density/pressure 或 `RuntimeStatus` 选择 timestep lifecycle。

状态：**Implemented**。

## 4. Deleted interfaces/files

已删除：

- `core/interfaces/SF_flowAlgorithm.h`
- `solver/algorithm/SF_solverAlgorithm.h/.cpp`
- `solver/algorithm/densityBased/SF_correction.h/.cpp`
- `solver/algorithm/densityBased/SF_densityBased.h`
- `solver/algorithm/pressureBased/SF_adapter.h/.cpp`
- `solver/algorithm/pressureBased/SF_pressureBased.h`

同时删除 `IFlowAlgorithm`、`FlowAlgorithmContext`、`FlowAlgorithmResult`、`makeFlowAlgorithm()`、`DensityBased::Algorithm`、`PressureBased::Algorithm`、`flowAlgorithm_` 和 `CompressibleAlgorithm::stepPressure()`。

状态：**Implemented**。

## 5. Operation provider ownership

| Operation family | Numerical owner | Control owner | Status |
|---|---|---|---|
| `flow.step.*`, `explicit.stage.execute` | `Time::Explicit`、`DensityBasedRHS` 和现有 boundary/transport helpers | `CompiledSolvePlan` | Implemented |
| `pressure.*`, `momentum.*`, `velocity.correct`, `flux.correct` | `PressureBased::Corrector` 和现有 RHS helper | `CompiledSolvePlan` | Implemented |
| `ibm.constraint.project` | `ImmersedAlgorithm::projectConstraint()` + `IImmersedConstraint` | `CompiledSolvePlan` | Implemented |
| `ibm.kkt.solve` | `PressureBased::MonolithicKKT` + `IImmersedConstraint` | `CompiledSolvePlan` | Implemented |
| `flow.step.commit`, `time.commit` | `CompressibleAlgorithm` 的 state publication callbacks | `CompiledSolvePlan` | Implemented |

未增加新的 flow algorithm 总接口。

## 6. IBM constraint migration

新增小型 numerical helper `solver/algorithm/immersed/SF_constraintOps.*`。`ibm.constraint.project` 直接完成 provider identity、capability、distributed ownership 校验，绑定 execution runtime，并调用 `projectPredictedState()`。

`ibm.kkt.solve` 保持独立 OpId，直接调用现有 `MonolithicKKT`。它不会 fallback 到 projection。Ghost/ILW 继续留在 boundary closure，没有变成 equation correction。

IBM applicability 现在由 `requiresPredictedState`、`requiresPressureCorrection`、surface/body/transfer/distributed capabilities 校验，不再由 `SolverAlgorithm::DensityBased/PressureBased` 选择。

状态：**Implemented**。

## 7. Pressure execution migration

原 `stepPressureGeneric()` 的数值 callback 被保留并重构为 `bindPressureOps()`。PISO 的 assemble、solve、pressure update、velocity correction、flux correction 和 commit 顺序仍来自 compiled Plan。完整 hidden `stepPressure()` 已删除。

没有 lowering 为 operation 的旧 single-fluid pressure 配置现在标记 `Unsupported`，不会回退到旧 lifecycle。

状态：PISO 为 **Implemented**；未 lowering 的旧配置为 **Unsupported**。

## 8. CompressibleAlgorithm responsibility before/after

| Responsibility | Before | After |
|---|---|---|
| timestep control flow | solver-family booleans和多个 `step*()` | `CompiledSolvePlan` |
| flow/IBM correction dispatch | `IFlowAlgorithm` | distinct IBM OpId callbacks |
| pressure dispatch | `RuntimeStatus` + fallback | required OpId + provider coverage |
| RK stage loop | `CompiledSolvePlan`/`Time::Explicit` | 不变 |
| numerical callback binding | 分散在各 `step*()` | `bindExplicitOps()` / `bindPressureOps()` / IBM callbacks |
| state/workspace binding | `CompressibleAlgorithm` | 不变 |

状态：**Implemented**。

## 9. RuntimeStatus remaining role

`RuntimeStatus` 只用于 resolved-system 报告和 capability 诊断。它不再出现在 `CompressibleAlgorithm::advance()` 或 `PlanExecutor` 的 numerical routing 中。

`LegacyAdapterRequired` 已随 pressure schedule capability 收敛删除；状态现在只有 `Runnable` 与 `Unsupported`。

状态：**Implemented**。

## 10. Remaining formulation-dependent branches

`SolverAlgorithm`/formulation 判断仍存在于 system composition、compatibility validation、single-phase pressure-corrector capability check、Eulerian pressure stepper validation和报告中。它们没有出现在 `CompressibleAlgorithm::advance()` 与 `PlanExecutor` 的 lifecycle dispatch 中。

状态：composition/validation 为 **Implemented**；`PressureBased::Corrector` 与 Eulerian stepper 内的 formulation capability check 为 **Legacy**。

## 11. Remaining Legacy / Unsupported paths

- `CompressibleAlgorithm` 类名仍带旧物理命名，但不再代表 runtime family：**Legacy**。
- single-fluid SIMPLE/PIMPLE 等尚未 lowering 为完整 operations 的组合：**Unsupported**。
- distributed `ibm.kkt.solve` 缺少 pressure-block `ConstraintGlobalDof` provider：**Unsupported**。
- Eulerian multiphase IBM constraint execution：**Unsupported**。
- `PlanNodeKind::Subcycle` execution：**Interface-only**。
- Eulerian pressure stepper 已使用 Plan/OpRegistry，但仍是独立 physical-system executor：**Legacy**。

不存在以 legacy single-fluid timestep fallback 伪装 runnable 的路径。

## 12. Architecture guards

`tools/check_architecture.py` 现在阻止旧文件和以下 symbols 回归：

```text
IFlowAlgorithm
FlowAlgorithmContext
FlowAlgorithmResult
makeFlowAlgorithm
flowAlgorithm_
DensityBased::Algorithm
PressureBased::Algorithm
CompressibleAlgorithm::stepPressure
```

守卫还检查 `CompressibleAlgorithm::advance()` 不得用 density/pressure enum 选择 control flow，并要求执行前调用 `PlanExecutor::validateBindings()`。

`test_pisoArchitecture` 增加 operation coverage 测试：缺少 `ibm.constraint.project` 时在任何 operation 执行前 fail fast；绑定后按 Plan 正常执行。

状态：**Implemented**。

## 13. Numerical invariants

本轮没有修改 Euler/SSPRK3/RK4 系数、stage time、WENO/TENO、flux splitting、CFL、RHS、pressure/KKT matrix、IBM multiplier、Ghost/ILW、boundary/halo 顺序、canonical face COPY、GlobalDof SUM、state layout 或 workspace ownership。

显式 path 的顺序仍为 prepare → dt → begin → stage loop → IBM operation（若有）→ flow commit → time commit。PISO operation 顺序与迁移前 compiled Plan 相同。

状态：**Implemented**。

## 14. Regression results

| Check | Result | Evidence/status |
|---|---|---|
| full build | pass | 299 build steps，无 compile/link failure；Implemented |
| architecture checker | pass | allowlisted dependency edges = 10；Implemented |
| WENO7 Euler, one step | pass | `1ab9ccf4cfb3083459de82fa2312e9eaafc0dbe4be896cdf354d9aafd9d55d07`；Implemented |
| Ghost IBM RK4, one step | pass | `b9d70a4855121102dfb351841332986ae12080862f0facc2aa696d194f329b0e`；Implemented |
| Pressure PISO control | pass | `11cba96d719f1b00ced3a01237ba489c636ffec13de3c0c516144af47ca88040`；Implemented |
| SSPRK3 explicit stage mathematics | pass | targeted CTest；Implemented |
| Velocity Forcing serial | pass | final SHA-256 `d63bc3a5601d8356d4f0f10e4e02b961c224c4724b4b7e2c2fc3d86b68465ab8`；Implemented |
| fractional DFM serial | pass | final SHA-256 `3724cb99fe5c067fcfdc70cff77e62be5a091b7d80773496089c42bf6d75c6d8`；Implemented |
| Velocity Forcing MPI-2 | pass | closure `max|Ju-Us| = 1.9613e-12`；Implemented |
| surface KKT serial | pass | 72 iterations，relative residual `8.04775e-10`；Implemented |
| `numericalFluxContract` | known failure | 迁移前 Rusanov constant-stencil mismatch；Legacy |

## 15. Global symbol search results

对 `src/`、`test/` 和 architecture checker 的全局搜索未发现已删除的 flow-algorithm symbols。`CompressibleAlgorithm::advance()` 与 `PlanExecutor` 中未发现 `SolverAlgorithm::DensityBased/PressureBased` lifecycle branch，也未发现 `RuntimeStatus` numerical routing。

状态：**Implemented**。

## 16. Final authority table

| Question | Authority | Status |
|---|---|---|
| timestep 做哪些事情 | `CompiledSolvePlan` | Implemented |
| operation 如何计算 | 对应 numerical provider callback | Implemented |
| callback 如何查找 | `OpRegistry` | Implemented |
| structured control flow 如何执行 | `PlanExecutor` | Implemented |
| explicit stage 数学 | `Time::Explicit` | Implemented |
| pressure correction 数学 | `PressureBased::Corrector` | Implemented |
| sequential IBM constraint 数学 | `IImmersedConstraint` provider | Implemented |
| monolithic IBM 数学 | `MonolithicKKT` provider | Implemented |
| physical state/time | `StateBundle` | Implemented |
| unsupported operation 如何处理 | pre-execution coverage failure | Implemented |

## 17. Recommended next task

下一任务应为仍标记 **Unsupported** 的 single-fluid pressure policies 实现 dedicated numerical providers。随后可收敛 `PressureBased::Corrector` 和 Eulerian stepper 中只用于 capability validation 的 formulation checks。不要在该任务中迁移 density RK 数学、修改 IBM/KKT 数值公式或增加新的 flow algorithm family。
