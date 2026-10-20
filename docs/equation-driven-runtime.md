# Equation-Driven Runtime Migration Report

日期：2026-09-16。

## 1. Before architecture

```text
ResolvedSimulationSystem -> Workflow stage text -> specialized stepper
```

workflow 曾从 strategy 字符串推断 predictor、correction 和 constraint，因而 solve
strategy 的 runtime authority 并不稳定。

## 2. Authority changes

现在每个 `SolveBlock` 同时保留 explain/compatibility 文本和
`SolveStrategyKind`。所有 runtime 分派均读取 enum：

```text
ResolvedSimulationSystem
  -> typed SolveBlock
  -> Workflow::Plan / SolveStage
  -> stepper composition and validation
```

`CompressibleAlgorithm` 以 explicit/predictor/pressure-correction strategy 选择
既有已验证的 lifecycle；`PressureStepper` 以 pressure-velocity-coupling strategy
识别 PIMPLE block；application runner 也不再用 `S_EE_PIMPLE` ID 选择 executor。

## 3. Compiled equation status

`AssemblyPlanRegistry` 已在启动阶段为每个 resolved equation 建表，并被 density
和 Eulerian assembler 消费；term vocabulary 已有 transient/divergence/gradient/
diffusion/explicit/implicit source/constraint/algebraic relation 以及 mode/time/
closure/boundary requirement metadata。

它尚未是完整 `CompiledEquation`：symbol 到 storage handle、fused operator 的
component write set 和 per-term synchronization requirement 尚待 lowering。不能把
当前 plan 误报为已经支持任意 scalar equation 的 generic executor。

## 4. Execution graph status

目标图保持为：

```text
Input
  ↓
Equation IR
  ↓
CompiledEquation
  ↓
ExecutionGraph
  ↓
ExecutionRuntime
  ↓
StateBundle
```

当前 density 的实际顺序已经是这个图的 production specialization：boundary/halo →
convection → canonical F* COPY → diffusion/source → GlobalDof SUM → RK stage update →
commit。它还由 `DensityBasedRHS` 和 `DensityBasedTime` 持有，而没有被一个空壳
`SystemExecutor` 包装。

## 5. Numerical and MPI behavior

本次仅改变 strategy authority。没有改动 WENO/TENO/Rusanov/central diffusion、RK
系数或 stage time、boundary/halo 顺序、IBM correction、canonical face COPY、
GlobalDof SUM、residual 符号或 storage layout。

`cmake --build build --parallel 4` 通过；`python3 tools/check_architecture.py` 通过，
allowlist 仍为 10。守卫新增禁止 workflow 使用 `strategy.find/rfind`。

`test/thermalCase/case.yaml`（density/Euler）完整推进至 `t=0.02`；其解释输出为
`S_FLUID ... kind=explicitTimeIntegration`。`test/eulerianEulerianCase/case.yaml`
（pressure/PIMPLE）完整推进至 `t=2e-4`；其解释输出为
`S_EE_PIMPLE ... kind=pressureVelocityCoupling`，末步 pressure residual 为
`3.413064e-16`。

Host CTest 的 5 个现存 executable 全部通过，包括三个 MPI/HYPRE 测试。第 6 个
`numericalFluxContract` 未运行，因为 `build-phase1a-tests/test_numericalFlux`
不存在；这不是测试断言失败。sandbox 下 MPI/PRTE socket 被拒绝，host 重跑后 MPI
测试通过。

## 6. Remaining specialized authority

- `CompressibleAlgorithm::stepDensity` 仍拥有 density lifecycle；
- `PressureStepper::stepImpl` 仍拥有 PIMPLE outer/non-orthogonal/turbulence 的嵌套顺序；
- homogeneous、level-set、phase-change 与 IBM auxiliary equation 仍通过 specialized
  coupling hooks；
- output/diagnostics 仍未按 unknown registry 自动注册。

## 7. Completion assessment

本 checkpoint 完成了 typed solve-strategy authority 和已有 AssemblyPlan 的 runtime
binding；尚未完成完整 equation-driven runtime。下一步应先编码已验证的 nested
placement，再让 compiler 显式 lower fused operator 的 read/write slices。这样能在不
改变数值基线的前提下迁移 density execution graph。
