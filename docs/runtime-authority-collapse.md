# Runtime Authority Collapse Report

日期：2026-09-18。本轮只重构组合、计划和运行期 authority；数值公式、RK stage time、边界/halo 次序、canonical face COPY、GlobalDof SUM 与线性矩阵结构均未修改。

## 1. Before dependency graph

`CaseConfig -> System::build -> CompiledSolvePlan -> Workflow::Plan -> SolveStage -> INavierStokesStepper::bindSolveStages -> numerical stepper`。

## 2. Duplicate authorities found

`CompiledSolvePlan::legacyBlocks` 和 `Workflow::Plan::stages` 都表达同一执行序列；`ExecutionBackendKind` 与 flow-family factory 都可路由 runtime；`PlanOperation` 和 `EquationOperatorBinding` 均是封闭的中心 enum。

## 3. Deleted files/types

**Implemented**：删除 `solver/algorithm/workflow/SF_workflow.*`、`SF_moduleGraph.*` 和旧 `solver/algorithm/SF_planExecutor.*`；删除 `PlanOperation`、`EquationOperatorBinding`、`ExecutionBackendKind`。

## 4. Moved legacy files

**Implemented**：通用 executor 移到 `solver/run/SF_planExecutor.*`。`densityBased`、`pressureBased` 和 IBM numerical provider 仍为 **Legacy**，尚未移动，因为移动本身不应改变已验证的数值执行顺序。

## 5. Rename table

| Before | After | Status |
|---|---|---|
| `legacyBlocks` | `blocks` | Implemented |
| `PlanOperation` | `OpId` | Implemented |
| `EquationOperatorBinding` | `OperatorId` | Implemented |
| `GenericPlanExecutor` | `Run::PlanExecutor` | Implemented |
| `PlanOperationRegistry` | `Run::OpRegistry` | Implemented |
| `bindSolveStages` | `bindSolvePlan` | Implemented |

## 6. New namespaces

**Implemented**：`SF::Run` 只包含 `PlanExecutor` 和 `OpRegistry`，不判断 PISO、SIMPLE、Eulerian、IBM 或 turbulence。

## 7. New Program model

**Interface-only**：当前 `ResolvedSimulationSystem` 仍是过渡期 Program value object，包含 Raw、Executable、Plan 与 RuntimeReport。文件拆分和 `Sys::Program` 正式改名尚未完成。

## 8. New Input model

**Implemented**：`EquationCompositionConfig` 用字符串 provider ID 保存 EOS、caloric thermo、transport 与 algorithm，不再使用四组 central enum。

## 9. Registry model

**Implemented**：`CompositionProviderRegistry` 可注册 EOS、thermo、transport、algorithm ID。注册 dummy provider 的测试不修改 core enum。生产 composition 目前仅注册 built-in provider；插件注入到 system compiler 属于后续工作。

## 10. Transform pipeline

**Implemented**：Raw -> transformer -> Executable -> planner 保持单向；transformer 没有 MPI implementation。

## 11. Binding pass

**Interface-only**：`CompiledEquation` 已改用开放 `OperatorId`；resource/storage binding 仍由 `PressureConstraintTransformer` 生成，尚未抽成独立 binding pass。

## 12. Plan model

**Implemented**：`CompiledSolvePlan` 只保存 structured root 与 `blocks`，不保存 backend 或 legacy adapter。

## 13. PlanExecutor

**Implemented**：executor 递归执行 `Sequence`、`Loop` 和叶节点，叶节点只通过 `OpId` 找 provider；无 backend gate。

## 14. Operation registry

**Implemented**：`Run::OpRegistry` 为开放 ID 到 callback 的显式映射。未知 operation 抛出包含 ID 的错误。

## 15. Removed solver-family dispatch

**Implemented（plan/runtime boundary）**：plan 与 `Run::PlanExecutor` 不再依赖 solver-family enum。**Legacy**：`CompressibleAlgorithm` 和 old correction adapters 仍在执行层按现有配置组合；彻底移除 `IFlowAlgorithm` factory 是下一步，不应在本轮把其数值校正直接改写。

## 16. Native vs legacy IO

**Interface-only**：native model path 现已保存 semantic composition ID，但 CaseAdapter 仍把部分 native object 映射给既有 compatibility dictionaries。尚未达到“native Input::Case 完全不经过 legacy dictionary”。

## 17. Preserved numerical kernels

**Implemented**：DensityBasedRHS/Time、WENO/TENO/flux、pressure corrector、IBM、boundary pipeline、MPI runtime 与 HYPRE provider 未改公式或调用次序。

## 18. Architecture tests

**Implemented**：architecture guard 验证 plan direct binding、开放 `OpId`、没有 plan backend/legacy-block authority，并继续确认 allowlist 为 10。

## 19. Production tests still working

**Implemented**：`SF_application` 静态库与 `test_pisoArchitecture` 构建成功；`pisoArchitecture` CTest 通过。未在本轮运行 production numerical case，因此不将数值回归标记为已完成。

## 20. Interface-only algorithms

**Interface-only**：constant-density U/p composition 能完整解析 Raw、transform、Executable、binding requirement、Plan 和 RuntimeReport；没有 numerical provider 时不执行。

## 21. Unsupported algorithms

**Unsupported**：constant-density pressure execution仍 fail-fast，原因明确为缺少 U/p predictor/pressure-correction provider；不会回退到 PerfectGas 或 legacy pressure path。

## 22. Remaining legacy debt

**Legacy**：`ResolvedSimulationSystem` 巨型头、`BuildRequest` feature bag、`PhysicsTemplateKind` provenance、`IFlowAlgorithm`/`makeFlowAlgorithm`、AssemblyPlanRegistry equation-ID authority、native-to-compatibility dictionary materialization、specialized runner dispatch。它们没有被包装为新的 authority。

## 23. Next numerical implementation step

在不改数学的前提下，将已有 pressure/density correction code 提取为显式 operation provider，并让 `CompressibleAlgorithm` 只执行已绑定 Plan；之后拆出 binding pass，再迁移 native `Input::Case`。
