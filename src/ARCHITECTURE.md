# SonicSolver 架构 — phase25B

**版本：** 2026-09-25 · phase25B
**状态：** 架构编译链已收口；部分 numerical provider 仍未实现。本文描述当前源码的实际 authority，不把计划接口等同于已经可运行的算法。

## 1. 系统以方程而非 solver family 为中心

一个 case 独立说明主未知量、方程/约束、模型贡献、耦合 preset、时间 recipe、空间离散和运行环境。内置 Navier–Stokes 方程、湍流、IBM 及用户声明通过同一 composition path 汇入 RawEquationSystem。Preset 可以给出普通贡献与控制流片段，但不能暗中选择第二套运行流程。

```text
CaseConfig + BuildRequest
        │
        ├─ Built-in presets
        ├─ Model SystemContribution
        └─ User additions / modifications
        ↓
RawEquationSystem                         WHAT: physical equations/constraints
        ↓ TransformationPipeline
ExecutableEquationSystem                  WHAT: derived equations + operations
        ↓ compileStateRealization
CompiledStateRealization                  STATE: unknown/storage roles
        ↓ NumericalCompiler
CompiledNumericalSystem                   HOW: term recipes, dt/time/phase policy
        ↓ PlanFragment + SolvePlanner
CompiledSolvePlan                         ORDER: sequence/loop/stage/OpId
        ↓ ProviderResolver
ResolvedOperationBinding[]               BINDING: OpId → provider/Unsupported
        ↓ RuntimeRequirements + validation
OpRegistry → PlanExecutor → committed StateBundle
```

`ResolvedSimulationSystem` 保存这些不同阶段的结果，而非一份可任意重复修改的平行 state。`Program` 只以引用暴露 executable、numerical、plan 和 runtime 四类编译产物。

## 2. 各层的唯一 authority

| 问题 | 当前 authority | 主要位置 |
|---|---|---|
| 解什么方程、有哪些约束 | `RawEquationSystem`；变换后为 `ExecutableEquationSystem` | `core/system/SF_equationIR.h` |
| 模型增加什么 | 中立 `SystemContribution` 值，经 `SystemCompositionBuilder` 验证合并 | `core/system/SF_systemContribution.h`、`solver/system/SF_equationContribution.cpp` |
| pressure/IBM 约束如何变换 | transformation descriptor 与 pipeline | `solver/system/SF_transformation.cpp` |
| 每项的离散和时间 recipe | `CompiledNumericalSystem` | `solver/system/SF_numericalCompiler.cpp` |
| 循环层级和执行顺序 | `CompiledSolvePlan` | `core/system/SF_planFragment.h`、`solver/system/SF_solvePlan.cpp` |
| 哪个实现执行每个 OpId | `ResolvedOperationBinding` | `solver/system/SF_providerResolver.cpp` |
| backend 与 MPI 能力是否可用 | `RuntimeRequirements`，由 `BuildCapabilities` 验证 | `solver/system/SF_runtimeRequirements.h`、`SF_systemValidator.cpp` |
| 物理状态、时钟与 solver workspace | `StateBundle`；workspace 归 execution lifetime | `core/state/`、`solver/algorithm/` |

`ExecutableOperation` 只声明 OpId、stage 和 typed `OperationCapability`。它不写 `flow.conservative` 之类的 concrete provider。ProviderResolver 以 Plan 中实际引用的 OpId、operation 需求和 state realization 解析现有实现；未声明或不支持的操作保留在报告中，状态为 `Unsupported`。`RuntimeReport` 是绑定结果的只读摘要，不能再发明 missing-operation 列表。Stepper 的 OpRegistry 仅保留已分配给它的 operation。

## 3. Composition 与 transformation

中立 IR 位于 `core/system`；模型可以贡献未知量、方程、term、constraint、closure、transformation request 和 execution policy，但不能 include/link `solver/system` 编译实现。`solver/system` 执行组合、校验和 lowering。Builtin、Model 与 User 使用同一 `Equation::Definition`/descriptor 形式；`add`/`extend` 已有 typed path，尚未完整实现的 `replace`/`disable` 必须显式失败，不能按同名覆盖。

Pressure constraint transformer 根据方程与 `C_INCOMPRESSIBILITY` 匹配，生成 predictor、pressure correction 与 velocity/flux correction 的方程、算子和 executable operations。Eulerian shared-pressure transformer 根据 `C_SHARED_PRESSURE` 声明相级 `ee.*` operations。Coupling preset 的 PlanFragment 给出顺序和三层循环：outer corrector → pressure corrector → non-orthogonal pass。`SF_solvePlan.cpp` 只降低结构化 fragment，不包含 PISO/PIMPLE 或 `ee.pressure.*` 的算法分支。Ghost/ILW 保持 boundary closure；约束型 IBM 才进入 constraint operation/耦合路径。

## 4. Numerical recipe、时间与配置

方程层只表达 `ddt/div/diffusion/source/constraint` 等数学 term。TimeRecipe 的当前内置选项为 `forwardEuler`、`SSPRK3`、`classicalRK4`；`SolvePlanner` 编译 stage topology，`Time::Explicit` 实现单 stage 数学。NumericalCompiler 为有真实消费者的 term 或 operation 绑定 recipe；选择了未被任何消费者使用的 diffusion recipe 会 fail-fast，不能因系统含压力约束而放宽检查。

Typed pressure 配置将 `PressureCouplingConfig`、`PhaseTransportConfig`、`PressureReference` 和 `LinearSolverSet` 分开。`PressureCouplingPreset` 是 preset 名称，不是顶层 Algorithm 类。native 输入中的 `algorithm: PISO` 仍可被解析为 typed preset。Eulerian phase convection 和 phase-source CFL 属于 phase transport 数值策略，编译为 `CompiledNumericalSystem::phaseTransport`；运行时读取编译值，当前仅 Upwind 可执行。

## 5. Production execution 与真实支持状态

`app/` 解析 case、建立 state/services 与编译结果，然后进入共同的 compiled-plan 执行 contract。`PlanExecutor` 只解释 Sequence、Loop、StageLoop 和 leaf OpId；现有 numerical helper 仍实现真实 kernel。Single-fluid explicit path 与 Eulerian shared-pressure path 使用各自适用的 provider，但 provider 选择已经在运行前完成。Pressure、IBM、turbulence 和 phase 源项仍按数学职责区别处理，不因为它们都叫 model 而合并数值实现。

| 路径 | phase25B 状态 |
|---|---|
| Conservative single-fluid explicit time，现有密度格式 | 已有 production provider |
| 现有 conservative PISO pressure-correction 能力范围 | 已有 provider；仅在匹配的 state、schedule 和 time recipe 下解析 |
| Eulerian shared-pressure PIMPLE 现有路径 | 已有相级 provider；plan 由 shared-pressure fragment 降低 |
| 内置 Ghost/ILW 及现有 IBM constraint 数值路径 | 保留现有实现，通过对应 boundary/operation contract 接入 |
| Constant-density single-fluid PISO | `Unsupported`：缺 dedicated predictor/pressure/correction provider |
| SIMPLE/PIMPLE single-fluid fixed-point predictor | `Unsupported`：不得重复执行 full-`dt` conservative momentum update |
| 其他未提供的时间、离散或 backend 组合 | capability validation 显式 fail-fast，不静默降级 |

`explain` 应显示 contributions、raw/executable equations、operations、recipe、Plan、operation bindings、运行要求及其状态。`Active` 的数学变换不等于对应 numerical provider 一定 Runnable；两者必须分开报告。

## 6. State、workspace 与并行数值不变量

每个 physical state、clock、equation binding 只有一个 owner。`StateBundle` 保存参与计算的 patch 状态及物理时钟；`PatchWorkspace`、FluxField、Residual、RK stage arrays、线性和约束 workspace 归求解执行期，不回塞进 `Field`，也不复制另一份 authoritative Q。

```text
state / geometry / labels          owner → COPY → replicas
shared numerical face              candidate → one canonical F* → COPY → +F*/−F*
residual / source / load / J row    local contributors → SUM → owner
```

MPI 实现留在 infrastructure/backend；equation、discretization 与 plan 只声明同步/归约要求。架构调整不得改变 stage time、boundary/halo 顺序、flux 或 residual 符号。Phase25B 四个 production 回归和限制见 [迁移报告](../docs/provider-resolution-shared-pressure-closure.md)。

## 7. 模块依赖与源码入口

```text
             core/system + core interfaces
                   ↑             ↑
                   │             │
                models      solver/system compilers
                                  ↓
                       algorithm / run execution
                                  ↓
             equation + discretization + boundary
                                  ↓
                       linearAlgebra / backend
```

这是语义方向；部分旧数值目标尚有链接债，不能把目录图误写成已完全无环。主要文件：

- `solver/system/SF_systemBuilder.cpp`：组合与编译阶段入口。
- `solver/system/SF_transformation.cpp`、`SF_pressureCoupling.cpp`：压力约束和 shared-pressure contributions。
- `solver/system/SF_solvePlan.cpp`、`SF_providerResolver.cpp`：通用计划 lowering 与唯一 provider 解析。
- `solver/run/SF_planExecutor.cpp`：结构化执行器。
- `solver/algorithm/SF_singleFluidStepper.cpp`、`eulerian/SF_eulerianStepper.cpp`：现有 provider 到 OpRegistry 的绑定。
- `app/application/`：case composition、state/services 装配和输出。

## 8. 验证和未完成工作

Phase25B 的干净 clang++/Ninja 构建、14/14 CTest 和 `tools/check_architecture.py` 均通过。4-rank Sod 与 serial Ghost 的逐步诊断记录一致；PISO 最终 VTS SHA-256 相同；Eulerian 两步诊断、Plan 树和 OpId 顺序一致。Sod/Ghost 迁移前未保存场文件哈希，因此报告只证明记录精度下诊断逐行相同。

仍需独立处理：constant-density single-fluid PISO numerical provider（Phase 26）、SIMPLE/PIMPLE dedicated fixed-point provider、Eulerian phase equation 的更通用 recipe lowering、`SF_equation ↔ SF_turbulence` 旧静态库链接关系，以及部分兼容输入。不要用 fallback、clamp、修改 CFL 或改变现有数值公式来让这些路径表面可运行。
