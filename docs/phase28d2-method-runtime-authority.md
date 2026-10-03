# Phase 28D.2 — Structured HOW 与 Method Object Production Authority

日期：2026-09-29。

## Before / After

| 方面 | Before | After |
|---|---|---|
| 方程级 HOW | 平面 `ExecutionProgram.steps` | `ProgramNode` root：Sequence、Step、Repeat、LoopUntil、Group、Commit；生产 composition 只提交 root |
| density 时间计划 | `SolvePlanner` 按 `TimeRecipe` 自行构造固定显式 stage | 已解析 `ITemporalMethod` 编译时间片段；Planner 只组合其他 policy leaf |
| density 方程目标 | Plan 的 legacy fallback 可从 `solvedUnknowns` 重建 | `Flow -> Q` 的 `ConservativeResidual` 输出与 FormulaGroup 顺序冻结；无片段时 fail-fast |
| fused backend 的公式选择 | Stepper 按 `E_MASS/E_MOMENTUM/E_ENERGY` 名称比对 | Stepper 比对方法编译产物与 Plan leaf；仅保留 `rho/rhoU/rhoE` 固定存储布局能力限制 |
| 旧 Equation 定义 | 可独立供应 term-provider 编译 | Flow 成员按 FormulaGroup 选择，Formula 与旧定义不一致在 numerical compilation 时失败；旧 term 表示仍供 provider 迁移使用 |
| 时间 backend 选择 | `Time::Explicit` 读取 recipe ID 决定算法 | TemporalMethod 冻结 `ExplicitStageBackend`、系数和 stage 数，旧 kernel 按冻结 backend 执行 |

## HOW hierarchy

`ExecutionProgram` 只描述方程级步骤与先后/重复/收敛；`ProgramStep` 的 primary `OutputRef` 仍是唯一目标来源。`compileExecutionProgram` 验证步骤、恰好一个方法绑定、workspace producer/consumer、时间 residual group，并产生 `CompiledFormulaStep`。通用 `compileMethodProgram` 可将 HOW root 降为 `SolvePlanNode`，嵌套 Repeat/LoopUntil 由现有 `PlanExecutor` 执行。`CompiledSolvePlan` 仍是较低层的 runtime topology；边界、halo、canonical face 和 RK stages 不进入用户级 HOW。

目前 pressure production 的 PISO/SIMPLE/PIMPLE 循环仍来自既有 PlanFragment。PressureUpdate、VelocityCorrection、FluxCorrection 在 structured HOW 中已有步骤与 Output 分类，但还未由通用 DirectEvaluation 方法代替已验证的 PressureOperators callback。此边界在 `explain` 标注为 legacy pressure numerical callback。通用 HOW 可以表达下一阶段的 predictor/corrector 嵌套；本阶段没有重写压力矩阵、Rhie–Chow 或 fixed-time 数值公式。

## Method object hierarchy 与编译职责

TemporalMethod（ForwardEuler、SSPRK3、ClassicalRK4）冻结 stage 系数、所选 backend 和带 prepare/dt/begin/stage/commit 的时间片段。EquationMethod（ConservativeResidual、DirectEvaluation、LinearEquation）验证数学 Formula 与 Output，给出读写、workspace、backend operation。Spatial/term provider 独立负责算子离散；编译后的方法步骤记录其 resolved provider ID。线性后端负责 `GlobalDofSystem` 的解算。`FormulaMode` 只留在旧 FormulaCompiler/PlanFragment 适配路径，不参与 migrated density 的 runtime dispatch。

已删除 `CompiledTimeRecipe::compile()` 对内置方法 registry 的回调入口；生产 NumericalCompiler 和 stage 数学测试均先解析具体 TemporalMethod，再由该对象编译时间数据。`CompiledTimeRecipe` 仅保存冻结结果。

Central compiler 的职责仅为 resolve、validate、调用方法对象 compile、compose；它不推断主未知量、不解释微积分、不按预设名选择另一套 solver。

通用 Formula 的单 patch `StateRealization`/`DistributedFieldView` 绑定已覆盖真实 Field scalar/vector 存储，直接赋值和 Cartesian Central2 线性组装/求解写回均经过测试。它不复制 Field，也不声称已支持任意并行 GlobalDof 编号、face scalar 通用 kernel 或用户 native Formula 输入。缺少边界 closure、存储或组件布局立即失败。

## Production authority 与准确边界

1. 改变 TemporalMethod binding 会改变生产的 stage 数、系数、后端和 `CompiledSolvePlan` 片段；FE/SSPRK3/RK4 仍复用原 `Time::Explicit` 数值 kernel，公式与执行顺序未改。
2. 改变 Flow FormulaGroup 会改变参与编译的 FormulaCalls 和 term-provider selection；当前 fused backend 仅接受现有 `rho/rhoU/rhoE` layout，不支持任意新公式的 fused 计算。
3. 旧 `Equation::Definition` 仍提供 term-provider 所需的过渡项记录，**不能独立改变**已迁移 Flow 的生产数学：与 Formula AST 不一致时启动编译失败。这是严格 migration guard，尚非完全删除旧表示。
4. `FormulaCall.mode` 仍用于通用 FormulaCompiler 和未迁移 PlanFragment；它不是 HOW 用户语义，也不决定 density production backend。

| Component | 状态 | Authority / 边界 |
|---|---|---|
| Formula / FormulaGroup | Implemented | WHAT；Flow 成员选择生产 term 编译 |
| ProgramStep / OutputRef / ProgramNode | Implemented | HOW；生产 density 使用 root |
| FE / SSPRK3 / RK4 | Implemented + legacy kernel | TemporalMethod 编译片段，`Time::Explicit` 保留已验证更新数学 |
| ConservativeResidual | Implemented | EquationMethod 绑定 fused conservative backend |
| DirectEvaluation / LinearEquation | 单 patch generic proof | 测试中以真实 Field 与 GlobalDof 写回；pressure 专用 provider 尚未迁移 |
| ConservativeRHS | Backend | WENO/TENO、flux、边界、canonical/SUM 数值实现，不选择 Formula |
| FormulaMode | Legacy/internal | 仅旧通用编译/PlanFragment 适配 |
| PressureOperators | Legacy backend | 生产压力数值公式未改 |
| `ExecutionProgram.steps` | Legacy input | 编译时归一为 root；不能与 root 同时提交 |

## 数值与并行不变量

未修改 WENO/TENO、通量分裂、Riemann 特征结构、RK 系数/顺序、边界/halo/IBM 顺序、canonical face 的 COPY、GlobalDof 的 SUM、残差符号或源项 clear 时机。单 patch generic Field 绑定使用已有非拥有 `DistributedFieldView`；不产生第二份 Q。新架构 guard 禁止 fused stepper 重新按固定 Formula ID 或 `FormulaMode` 选择执行。

## 验证

| 检查 | 结果 |
|---|---|
| `cmake --build build --parallel 4` | 通过 |
| `python3 tools/check_architecture.py` | 通过；2 条既有 allowlist，未新增依赖方向违规 |
| Host `ctest --test-dir build --output-on-failure --parallel 1` | 最终源码完整重跑 27/27 通过，包含 Field/structured HOW、MPI/HYPRE 与 pressure regression |
| 4-rank Sod，两步 | `dt=6.681531e-4, 5.575894e-4`；最终 `min(rho)=0.273944`、`min(p)=27143.7`，与 28D.1 相同 |
| serial Ghost IBM，两步 | `dt=4.392680e-4, 4.392260e-4`；最终 `min(rho)=1.22285`、`min(p)=101109`，与 28D.1 相同 |
| `git diff --check` | 通过 |

受限 sandbox 中完整 CTest 的 3 个 MPI 项曾因 PRTE 无法绑定 socket 而未启动；host 环境重跑 27/27 通过。这不是数值断言失败。

## Remaining implementation work

Production pressure 的 structured HOW 尚未取代 PlanFragment 调度；专用 PressureOperators 仍是数值 backend。通用 FormulaField binding 暂仅单 patch cell scalar/vector；face scalar、分布式编号和任意 Equation 的 production operator lowering 尚未闭环。Flow 的旧 Equation term 表示仍由 NumericalCompiler 消费，但受到 Formula equality guard。`Time::Explicit` 的三个已验证 kernel 仍执行实际更新；方法对象编译它们的 stage 参数与 backend，而非复制更新公式。Native 自定义公式、通用湍流/VOF/level-set/KKT 和新的 implicit RK 都未在本轮实现。

后续进展：Phase 28D.3 已把 constant-density pressure 的五个 equation-level ProgramStep 绑定到专用 EquationMethod 编译片段，PlanFragment 只保留宏观 loop/control skeleton；详见 [Phase 28D.3 报告](phase28d3-pressure-method-fragments.md)。本报告以上正文保留 Phase 28D.2 当时状态，不能将其 pressure DirectEvaluation 描述理解为最新 production authority。旧 Equation term-provider 元数据及 PressureOperators 算术仍未迁移。

SonicSolver 的 central compiler 不决定数学意图：WHAT 已给出数学，HOW 已给出步骤与 Output，HOW NUMERICALLY 已选择具体方法对象；compiler 只解析、验证、编译并组合这些已经明确的契约，数值数学由各方法对象负责，runtime 只执行冻结后的计划。
