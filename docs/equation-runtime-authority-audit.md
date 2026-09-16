# Equation Runtime Authority 审计

日期：2026-09-16。

## 1. 当前真实 runtime authority graph

```text
CaseConfig
  -> ResolvedSimulationSystem
      -> Equation::Definition / SolveBlock
      -> Workflow::Plan
          -> INavierStokesStepper
              -> CompressibleAlgorithm | Eulerian::PressureStepper
                  -> specialized numerical operators
                      -> ExecutionRuntime -> StateBundle
```

`ResolvedSimulationSystem` 是 unknown、方程、约束、solve block、capability 和
workspace requirement 的描述 authority。`StateRealizer` 在时间循环前将它绑定到
既有 `StateBundle` storage。`ExecutionRuntime` 是 halo、canonical face COPY 与
GlobalDof SUM 的唯一执行边界。

## 2. Density-based call graph

```text
runFlow
  -> bindSolveStages / prepare
  -> CompressibleAlgorithm::advance
  -> stepDensity
  -> DensityBasedTime::advance (Euler/SSPRK3/RK4 stage timing)
  -> DensityBasedRHS::assembleAllPatches
  -> boundary -> halo -> convection -> canonical F* COPY
     -> diffusion/source -> GlobalDof SUM
  -> stage update -> flow/IBM correction -> commit
```

现有 production WENO/TENO、flux split、central viscous、source、IBM 与 RK kernel
仍属于 specialized operator implementation；它们不应被通用解释器重写。

## 3. EquationDefinition 到 RHS 的数据流

density 的 `E_MASS/E_MOMENTUM/E_ENERGY` 由
`ResolvedSimulationSystem::equationDefinitions` 创建；`Equation::Compressible::System`
从同一对象建立 `AssemblyPlan`。plan 缺少对流、扩散或 source term 时，相应现有
specialized operator 不会被请求。Eulerian `PhaseEquationAssembler` 也从同一个
`AssemblyPlanRegistry` 读取 phase continuity/momentum/enthalpy 和 shared-pressure
plan，并在调用 kernel 前校验所需 term。

当前限制是 fused conservative convection/viscous kernel 仍以完整 conservative
component bundle执行。单独删去 energy diffusion 并不安全地等同于“保留 momentum
viscous 而屏蔽 energy component”；这需要先把 fused operator 的输入/输出 slice
作为 compiler-lowering 的明确 contract，不能靠在 stepper 内加入 equation-id 分支。

## 4. String semantic dispatch

本次迁移前，`Workflow::Plan` 用 `strategy.find("PIMPLE")`、`find("correction")`
等从显示字符串推断 stage。现在 `SolveBlock::strategyKind` 是唯一 runtime strategy
authority，workflow 只 switch 该 enum；显示字符串仅用于旧输入兼容和 explain。

仍保留的 string identity：equation/unknown/constraint ID 用于 registry、诊断和
explain。hot numerical path 尚未完全消除 symbol-name lookup，尤其是 fused
compressible operator 与 phase/turbulence specialized storage。

## 5. Duplicated authority

| 信息 | 当前 authority | 仍需收敛处 |
|---|---|---|
| equation mathematics | `Equation::System` | lifecycle coupling 的 auxiliary equation 仍由 specialized hook 执行 |
| solve strategy | typed `SolveBlock::strategyKind` | nested placement（RK stage/PIMPLE outer loop）尚未编码 |
| runtime order | DensityBasedTime、PressureStepper | 尚无能保持 nested placement 的 `ExecutionGraph` |
| state storage | `StateBundle` + StateRealizer binding | turbulence/IBM/auxiliary specialized storage 尚未统一 handle |
| MPI synchronization | `ExecutionRuntime` | 无重复 authority；必须保持 COPY/SUM 不变 |

## 6. 本轮最小修改文件

- `src/core/interfaces/SF_solveStrategy.h`
- `src/core/interfaces/SF_interfaces.h`
- `src/solver/system/SF_resolvedSimulationSystem.h`
- `src/solver/system/SF_systemBuilder.cpp`
- `src/solver/system/SF_systemValidator.cpp`
- `src/solver/algorithm/workflow/SF_workflow.cpp`
- `src/solver/algorithm/SF_compressible.cpp`
- `src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp`
- `src/app/application/run/SF_executionBuilder.cpp`

## 7. 下一安全迁移顺序

1. 将 `AssemblyPlan` lowering 扩展为带 resolved symbol/storage/operator/read-write
   handle 的 `CompiledEquation`。
2. 为 fused compressible convection/viscous kernel 写出 term-to-component-slice
   contract；默认配置保持现有全量 fused path。
3. 给 `SolveBlock` 增加已验证的 nested placement metadata，并以 RK/PIMPLE
   checkpoint 固定顺序。
4. 再建立真正遍历 compiled plan 的 density execution graph。

在第 3 步之前，禁止按现有顶层 block 列表创建通用 executor；那会改变 turbulence
和 auxiliary equation 的 stage/outer-iteration 位置。
