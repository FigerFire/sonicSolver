# Single-Fluid Pressure Schedule Semantics Completion

日期：2026-09-19

## 1. Scope

本轮只修正 pressure execution policy、structured Plan lowering、provider capability report 和 fail-fast guard。没有实现新的 SIMPLE/PIMPLE predictor，也没有修改 pressure、momentum、flux 或 time-integration 数值公式。

## 2. Previous inconsistency

single-fluid `S_PRESSURE` 原来把：

```text
repeatCount       = pressureCorrectors
nestedRepeatCount = nonOrthogonalCorrectors + 1
```

而 Eulerian policy 使用：

```text
repeatCount       = outerCorrectors
nestedRepeatCount = pressureCorrectors
innerRepeatCount  = nonOrthogonalCorrectors + 1
```

因此相同字段表达不同控制层级，single-fluid PIMPLE 也丢失了 outer-corrector 信息。

## 3. Canonical policy semantics

所有 pressure-velocity policies 现在统一为：

```text
repeatCount       = outer correctors
nestedRepeatCount = pressure correctors
innerRepeatCount  = non-orthogonal passes
innerRepeatCount  = configured nonOrthogonalCorrectors + 1
```

该映射同时用于 built-in single-fluid、semantic constant-density composition 和 Eulerian-Eulerian policy。

状态：**Implemented**。

## 4. Typed pressure policy

| Algorithm | ExecutionPolicyKind | SolveStrategyKind | Runtime provider |
|---|---|---|---|
| PISO | `SegregatedPressureCorrection` | `PressureCorrection` | Implemented for the current supported schedule |
| SIMPLE | `PressureVelocityFixedPoint` | `PressureVelocityCoupling` | Unsupported |
| PIMPLE | `PressureVelocityFixedPoint` | `PressureVelocityCoupling` | Unsupported |

`strategyName` 只保留给 explain/provenance。单元测试会修改该字符串并验证 Plan lowering 不变。

## 5. PISO lowering

当前可执行 PISO Plan 明确显示：

```text
prepare
outerCorrectors = 1
  momentum predictor
  pressureCorrectors
    nonOrthogonalCorrectors = 1
      pressure assemble/solve
    velocity/pressure/flux correction
commit
```

当前 numerical provider 的能力边界保持为：Euler time、一个 outer pass、至少一个 pressure corrector、一个实际 non-orthogonal pass、无 auxiliary schedule。已有 `momentum.solve` 仍只执行一次 full-dt physical update。

状态：**Implemented**。

## 6. SIMPLE/PIMPLE structural lowering

SIMPLE/PIMPLE 现在会生成完整三层 control-flow Plan，所有 executable leaves 都具有明确的 `pressure.schedule.*` OpId。例如：

```text
PressureSchedule.outerCorrectors
  pressure.schedule.predictor.assemble
  pressure.schedule.predictor.solve
  PressureSchedule.pressureCorrectors
    PressureSchedule.nonOrthogonalCorrectors
      pressure.schedule.pressure.assemble
      pressure.schedule.pressure.solve
```

这些 OpId 描述未来 dedicated fixed-point provider 的接口。目前没有绑定实现，因此不会执行。

状态：Plan 为 **Interface-only**；numerical provider 为 **Unsupported**。

## 7. Full-dt predictor guard

fixed-point Plan 不包含现有 `momentum.solve`。该 operation 调用 `Time::forwardEuler(..., dt)`，会真实推进 conservative state；在 outer loop 中重复它会重复使用完整 dt。本轮通过 dedicated OpId、runtime capability check、architecture guard 和单元测试共同冻结这一约束。

状态：**Implemented**。

## 8. Transitional lowering removal

`compilePolicy()` 不再为 pressure policies 生成 `operation == ""` 的 Assemble/Solve/Correct leaves。Boundary closure 也不能作为无 OpId 的独立 executable leaf；它必须由 owning numerical schedule 消费。

`System::validate()` 现在拒绝任何没有 OpId 的 Plan leaf。

状态：**Implemented**。

## 9. Runtime capability reporting

`RuntimeReport` 现在包含：

```text
status
requiredOperations
missingOperations
reason
```

状态只剩 `Runnable` 和 `Unsupported`；未使用的 `LegacyAdapterRequired` 与 `requiredAdapters` 已删除。

对于 SIMPLE/PIMPLE，`explain` 会同时打印完整 Plan、全部 required operations、全部 missing providers，并明确说明 current `momentum.solve` 不能作为 outer fixed-point predictor。

状态：**Implemented**。

## 10. Provider availability decision

single-fluid pressure 的 `Runnable` 不再由“Plan 中存在 OpId”推导。只有 typed PISO policy 同时满足当前 equation、time 和 schedule capability 时才为 `Runnable`。其余 pressure schedules 保留 Plan，但标记 `Unsupported`。

Eulerian PIMPLE 的现有 operation provider 和 Plan 不受此 single-fluid capability 判断影响。

状态：**Implemented**。

## 11. Architecture guards

新增或扩展的守卫包括：

- pressure Plan lowering 禁止比较 `strategyName == SIMPLE/PISO/PIMPLE`；
- single-fluid policy 必须保持统一的三层 count mapping；
- unsupported fixed-point Plan 必须使用 dedicated predictor OpId；
- compiled Plan 的每个 executable leaf 必须有 OpId；
- `LegacyAdapterRequired` 和 `requiredAdapters` 禁止重新出现；
- 单元测试验证 SIMPLE/PIMPLE 不包含 `momentum.solve`。

状态：**Implemented**。

## 12. Numerical behavior

没有改变：

```text
Time::forwardEuler
PISO pressure assembly/solve/correction
momentum residual
pressure matrix
flux correction
boundary/halo order
canonical face COPY
GlobalDof SUM
MPI backend
```

PISO operation 顺序与迁移前相同，仅在 Plan 中显式增加 repetition=1 的 outer/non-orthogonal structural nodes。

状态：**Implemented**。

## 13. Validation

| Check | Result |
|---|---|
| full build | pass |
| architecture checker | pass; allowlist remains 10 |
| `pisoArchitecture` | pass |
| `explicitStageMathematics` | pass |
| `pisoNumericalRegression` | pass; frozen VTS hash unchanged |
| PISO explain | Runnable; missing providers = none |
| PIMPLE 2×3×2 explain | Unsupported; all `pressure.schedule.*` providers listed |
| PIMPLE production launch | fail-fast before timestep |
| `numericalFluxContract` | known pre-existing Rusanov constant-stencil failure |

## 14. Deferred numerical implementation

- dedicated SIMPLE/PIMPLE fixed-point momentum predictor：**Unsupported**。
- relaxation and convergence provider for fixed-point schedules：**Unsupported**。
- single-fluid non-orthogonal pressure provider beyond the current one-pass PISO contract：**Unsupported**。
- constant-density pressure numerical provider：**Unsupported**。

下一任务应先定义 dedicated fixed-point predictor 的 state/time semantics，再逐个实现 `pressure.schedule.*` providers。不得把 full-dt `momentum.solve` 包进 outer loop。

