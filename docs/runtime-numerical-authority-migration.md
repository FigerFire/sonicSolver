# Runtime Numerical Authority Migration

## 1. Scope

本轮只迁移 numerical runtime authority，不改变 governing equations、离散格式、RK 系数/时间、pressure correction、IBM/HYPRE、boundary/halo、canonical COPY 或 GlobalDof SUM。

## 2. Before dependency graph

```text
CompiledSolvePlan -> CompressibleAlgorithm / PressureStepper
  -> hidden step()/stepImpl()/ddtDispatch()
  -> DensityBasedRHS / Corrector / Eulerian equations
  -> StateBundle + ExecutionRuntime
```

## 3. After dependency graph

当前已实现的 minimal pressure route：

```text
CompiledSolvePlan -> Run::PlanExecutor -> Run::OpRegistry
  -> existing DensityBasedRHS / forwardEuler / PressureBased::Corrector
```

Eulerian pressure 和 density RK 仍为 Legacy 路线。

## 4. Pressure lifecycle before

`EulerianEulerian::PressureStepper::stepImpl()` 的真实顺序：

```text
reset diagnostics -> refresh row map -> boundary+sync -> validate
outerCorrector:
  interphase -> sources -> turbulence prepare/sync -> source validation
  momentum diagonal -> interpolated flux -> canonical phase flux
  continuity -> boundary+sync -> momentum predictors -> semi-implicit interphase
  boundary+sync -> diagonal sync -> interpolated flux -> canonical phase flux
  pressureCorrector:
    nonOrthogonal loop:
      pressure solve -> set p' -> p' sync -> phase correction
      canonical face flux correction -> boundary+sync -> canonical phase flux
  phase energy -> turbulence -> boundary+sync -> validate
commit time level -> turbulence commit -> diagnostics
```

`PressureStepper` 因此仍拥有 outer/correction/non-orthogonal loop authority。

## 5. Pressure lifecycle after

**Implemented**：single-field minimal PISO plan 表达 `prepare -> momentum assemble/solve -> pressure correction loop -> commit`，并由 `Run::PlanExecutor` 执行。

**Legacy**：Eulerian PIMPLE/PISO 的完整 outer、correction 和 non-orthogonal loops 尚未进入 plan；它们仍在 `PressureStepper::stepImpl()`。

## 6. Density lifecycle before

`CompressibleAlgorithm::stepDensity()` 计算 dt、开始 trace、准备 boundary 后调用 `DensityBasedTime::ddtDispatch()`。每个 stage 的 RHS 路径为：boundary/halo/IBM ghost -> equation-system RHS prepare -> candidate convection -> canonical F* COPY and +/- residual -> diffusion/sources -> coupling RHS -> local residual -> GlobalDof SUM；stage update 后 publish state、refresh/validation。integrator 后才运行 `IFlowAlgorithm::correct()`、commit、final boundary 与 clock commit。

## 7. Density lifecycle after

**Legacy**：RK stage lifecycle 仍由 `DensityBasedTime::ddtDispatch()` 驱动，未被错误拆分或重排。`SF_PLAN_TRACE` 等价 trace 尚未加入；已有 `HighOrderTrace` 保留 checkpoint。

## 8. Operation provider table

| OpId | Provider | Reads/writes | Status |
|---|---|---|---|
| `pressure.prepare` | `CompressibleAlgorithm` lambda | StateBundle, boundary/runtime | Implemented |
| `momentum.assemble` | `DensityBasedRHS::assembleAllPatches` | Q -> residual/flux | Implemented |
| `momentum.solve` | `Time::forwardEuler` | residual -> Q | Implemented |
| `pressure.*` / correction ops | `PressureBased::Corrector` | p', velocity, flux | Implemented |
| Eulerian pressure operations | `PressureStepper::stepImpl` | phase system/HYPRE | Legacy |
| density RK stage operations | `DensityBasedTime` | Q/RK storage/RHS | Legacy |

## 9. OpId table

开放 `OpId` 已覆盖 minimal pressure plan：`pressure.prepare`、`momentum.assemble`、`momentum.solve`、`pressure.boundary.prepare`、`pressure.assemble`、`pressure.solve`、`pressure.update.prepare`、`velocity.correct`、`flux.correct`、`pressure.correction.commit`、`pressure.step.commit`。

## 10. Which old wrappers lost authority

**Implemented**：Workflow、ModuleGraph、SolveStage lowering、backend enum 和 central operation enum 已失去并删除 runtime authority。

## 11. Deleted files/types

**Implemented**：见 `docs/runtime-authority-collapse.md`。

## 12. Remaining Legacy

`PressureStepper`、`CompressibleAlgorithm`、`IFlowAlgorithm`、`makeFlowAlgorithm`、`DensityBasedTime`、`AssemblyPlanRegistry`。

## 13. Binding status

**Interface-only**：`OperatorId` 已开放；`PressureConstraintTransformer` 仍生成 resource binding，BindingPass 尚未抽出。

## 14. AssemblyPlanRegistry status

**Legacy**：EquationId dispatch 仍存在，尚未将单条 equation 迁到 term -> OperatorId provider。

## 15. Numerical invariants preserved

源码迁移未修改 convection/diffusion/source formula、RK tableau、stage time、pressure/HYPRE formula、boundary/halo、canonical COPY、GlobalDof SUM。

## 16. Serial regression

**Unsupported**：本阶段仅完成源码审计，未运行 production serial numerical regression。

## 17. MPI regression

**Unsupported**：本阶段仅完成源码审计，未运行 1/2/4-rank regression。

## 18. First-difference analysis if any

**Unsupported**：没有执行新数值路径，故无差异。

## 19. Architecture tests

**Implemented**：现有 `pisoArchitecture` 测试通过，architecture guard 通过。

## 20. Unsupported paths

constant-density pressure numerical provider 继续 fail-fast；不回退到 PerfectGas 或其他 backend。

## 21. Next step

先为 `PressureStepper::stepImpl()` 提取不改公式的 operation callbacks，并在 `CompiledSolvePlan` 编译完整 outer/correction/non-orthogonal loop；每一个 provider 调用现有 helper。完成 pressure serial/MPI 回归后，才处理 density RK 和 BindingPass。

## Authority table

| Concern | Before authority | After authority | Status |
|---|---|---|---|
| control flow | algorithm + Workflow | CompiledSolvePlan（minimal pressure） | Implemented / Legacy |
| leaf execution | algorithm methods | OpRegistry provider（minimal pressure） | Implemented |
| pressure loops | PressureStepper | PressureStepper | Legacy |
| RK stage control | DensityBasedTime | DensityBasedTime | Legacy |
| numerical formula | kernels | unchanged kernels | Implemented |
| runtime state | StateBundle | StateBundle | Implemented |
| MPI | ExecutionRuntime | ExecutionRuntime | Implemented |
| resource binding | transformer | transformer | Interface-only |
| equation lowering | AssemblyPlanRegistry | AssemblyPlanRegistry | Legacy |
