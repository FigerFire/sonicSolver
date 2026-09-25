# Runtime Authority Audit

日期：2026-09-17。审计范围是当前 `src/` 的运行期组合、规划和执行链；本报告记录迁移前事实，不改变数值实现。

## 1. 当前调用与数据流

```text
CaseConfig / CaseIO
  -> Application::inspectCase
  -> System::build + validate
  -> ResolvedSimulationSystem
  -> Workflow::makePlan / ModuleGraph
  -> runFlow
  -> INavierStokesStepper::bindSolveStages
  -> CompressibleAlgorithm / PressureStepper
  -> numerical kernels, boundary pipeline, ExecutionRuntime
```

`CompiledSolvePlan` 同时包含结构化 `root`，并含有 `legacyBlocks` 和 `backend`。前者被 `GenericPlanExecutor` 消费，后者经 `Workflow::Plan` 降级成 `SolveStage` 后由 `INavierStokesStepper::bindSolveStages()` 消费。因此同一执行顺序存在两个可写、可解释的表示。

## 2. Authority table

| 层 | 当前 authority | 输入 | 输出 | 重复/冲突 | 调用者 | 迁移动作 |
|---|---|---|---|---|---|---|
| case input | `CaseConfig` 与解析期 composition | YAML 与旧 config | `BuildRequest` | EOS/algorithm 枚举和字符串来源并存 | `inspectCase` | 将 model/algorithm identifier 保持为开放 ID |
| composition | `System::build` | feature bag `BuildRequest` | raw/executable/plan | `templateOrigin`、formulation、EOS 同时影响系统形状 | application | 让方程贡献先决定 raw system，EOS 只贡献 closure/simplification |
| transformation | transformer registry | raw system | executable system | transformer 直接写 concrete `EquationOperatorBinding` 和 storage slice | system builder | 用公开 operator ID 和独立 binding pass |
| solve planning | `CompiledSolvePlan` | executable + policy | `root`、`legacyBlocks`、`backend` | plan 同时是 IR 与 backend selector | workflow、algorithms | plan 仅描述 control flow；runtime status 移到编译报告 |
| workflow lowering | `Workflow::Plan` | `legacyBlocks` | `SolveStage` | 复制计划顺序，丢失 loop/operation identity | `runFlow` | 删除 Workflow plan，stepper 直接绑定 compiled plan |
| module graph | `Workflow::ModuleGraph` | resolved system | names/modules | 只作并行解释、非 runtime authority | application | 删除；explain 从 resolved system 直接生成 |
| runtime dispatch | `ExecutionBackendKind` 与 `CompressibleAlgorithm::flowAlgorithm_` | plan/backend + formulation | density/pressure family path | backend 和 solver-family 均可选择 lifecycle | `run.cpp`, compressible algorithm | 以 operation binding / runtime capability 报告替代；legacy providers 隔离 |
| operation binding | `PlanOperation` enum | plan leaf | function callback | 增加 operation 要修改中心 enum 和所有 consumers | generic PISO execution | 改为 `OpId` 字符串和开放 `OpRegistry` |
| equation operation binding | `EquationOperatorBinding` enum | transformer | compiled equation | 新 operator 需要修改中心 enum | pressure transformer | 改为 `OperatorId`，执行期 registry 解析 |
| numerical execution | density RHS, pressure stepper, IBM/MPI kernels | bound state/workspace | updated physical state | 目前被 flow-family wrapper 间接选择 | legacy algorithms | 保留公式与顺序；只由 operation provider/adaptor 调用 |

## 3. 文件级问题

- `src/solver/algorithm/workflow/SF_workflow.{h,cpp}`：把 `CompiledSolvePlan::legacyBlocks` 降级为第二份 `SolveStage` 执行计划；违反“plan 是唯一运行期 authority”。
- `src/solver/algorithm/workflow/SF_moduleGraph.{h,cpp}`：根据 `formulation` 和 template 派生 application 模块图，形成与 resolved system 平行的解释入口。
- `src/core/interfaces/SF_interfaces.h`：`INavierStokesStepper::bindSolveStages()` 使 interface 依赖 Workflow lowering，而非 compiled plan。
- `src/core/interfaces/SF_flowAlgorithm.h`、`src/solver/algorithm/SF_solverAlgorithm.{h,cpp}`：`IFlowAlgorithm` 以 density/pressure family 选择 correction lifecycle。
- `src/solver/system/SF_resolvedSimulationSystem.h`：`PlanOperation`、`EquationOperatorBinding`、`ExecutionBackendKind` 把可扩展 operation/provider 固化在中心 enum；同时承载 raw、executable、plan、runtime routing 多种概念。
- `src/solver/algorithm/SF_planExecutor.{h,cpp}`：名字为 generic，但限定 `GenericPisoPlan` backend，且只接受 PISO central enum。
- `src/solver/system/SF_systemBuilder.cpp`：EOS 选择 `rhoConst`/`perfectGas` 时直接选择不同 equation composition，混淆 closure 与方程系统 authority。
- `src/app/application/run/SF_runFlow.cpp`：在 time loop 前把 Workflow stage bind 给 stepper，形成 `CompiledPlan -> Workflow -> Stepper` 的强耦合。

## 4. 数据 ownership

`StateBundle` 仍是物理状态、clock、equation binding 的唯一 authority；`PatchWorkspace` 仍归 solver execution。审计未发现本轮计划要求复制 `Field`、`FluxField`、`Residual` 或改写 canonical face COPY / GlobalDof SUM 的必要性。迁移只应替换 plan/interface binding，不可改变 boundary、halo、stage、canonical flux、residual 或 MPI 顺序。

## 5. 目标链

```text
Input::Case
  -> SystemComposer
  -> RawEquationSystem
  -> TransformationPipeline
  -> ExecutableEquationSystem
  -> CompiledSolvePlan
  -> Run::PlanExecutor + Run::OpRegistry
  -> numerical operation providers / explicitly marked legacy adapters
```

计划本身不保存 backend 选择；运行期可执行性、缺失 operation 与 legacy adapter 使用情况属于 compile/runtime report。未注册 operation 必须带 operation ID fail-fast。
