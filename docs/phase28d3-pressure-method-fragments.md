# Phase 28D.3 — Pressure EquationMethod Fragment Authority

日期：2026-09-29。

> 后续状态：Phase 29 已删除本报告所述常密度 pressure PlanFragment macro authority、single-fluid BoundTerm 和 backend AssemblyPlan 数学开关；本报告保留历史阶段含义。最终现状与验证见 [Phase 29 收尾报告](phase29-single-fluid-authority-closure.md)。


## 目标与层级

本轮属于 **HOW NUMERICALLY 的 equation-level method 编译**，不是新增物理项或方程，也没有更换压力离散。physical state 仍由原有 state/Field 持有；`PressureOperators` 仍持有 `HbyA`、`rAU`、pressure matrix、`pPrime` correction 与 face-flux workspace。`CompiledSolvePlan` 仍是运行顺序 authority，`PlanExecutor` 只执行 OpId。

| 层 | Before | After |
|---|---|---|
| WHAT | Transformation 生成 predictor、pressure 与 correction Formula；旧 `E_MOMENTUM` 同时供应 BoundTerm | Formula 不变；`E_MOMENTUM` 作为显式 source-math input 记录，BoundTerm 仍是 legacy 空间项源 |
| HOW | pressure PlanFragment 同时写 FormulaCall/target/mode 与 assemble/solve 叶子 | 五个 ProgramStep 显式写 Formula → Output；PlanFragment 的数学位置只引用 step subject，仍持有 outer/pressure/non-orthogonal loop 与 fixed-time controls |
| HOW NUMERICALLY | pressure 的三个 correction step 仅分类为 DirectEvaluation，生产操作顺序由 PlanFragment 手写 | 五个 Pressure EquationMethod 生成 `SolvePlanNode` 局部片段；fragment 是 operation 顺序 authority |
| Runtime | OpId → SingleFluidStepper callbacks → PressureOperators | 同一 callbacks 与数值后端；只改变 OpId 的编译来源 |

## 五个方法片段

| ProgramStep → Output | EquationMethod | Compiled fragment | 关键 workspace/数值要求 |
|---|---|---|---|
| `E_MOMENTUM_PREDICTOR → U` | `PressureMomentum` | `MomentumAssemble → MomentumSolve` | 显式输入 `E_MOMENTUM`；产生内部 `HbyA`、`rAU`、`faceResponse`、`correctedFlux`；velocity boundary、Rhie–Chow 与 canonical COPY 保持原回调 |
| `E_PRESSURE → pPrime` | `PressureCorrection` | `PressureBoundaryPrepare → PressureAssemble → PressureSolve` | 需要 `faceResponse`、`correctedFlux`；写 `pressureCorrection` workspace；保留 gauge、GlobalDof、线性后端与 boundary |
| `F_VELOCITY_CORRECTION → U` | `VelocityCorrection` | `VelocityCorrect` | 需要 `pressureCorrection`、`rAU`；保留 velocity boundary |
| `F_PRESSURE_UPDATE → p` | `PressureUpdate` | `PressureUpdatePrepare` | 需要 `pressureCorrection`；保留 gauge restore 与 pressure boundary |
| `F_FLUX_CORRECTION → phi` | `FluxCorrection` | `FluxCorrect` | 需要 `pressureCorrection`、`faceResponse`；产出 `pressureFaceFlux`；保留 Rhie–Chow、canonical COPY 与 continuity diagnostic |

`HbyA` 是动量方法内部构造的 numerical workspace，不是新的 physical field 或 restart state；它在 compiled method 的 workspace contract 中显式展示。`pPrime` 仍使用 transformation 声明的 transient `pressureCorrection` storage。`phi` 保持现有 face-flux workspace，未新增独立存储副本。

`PressureUpdate`、`VelocityCorrection`、`FluxCorrection` 保留专用方法，因为现有回调除 Formula assignment 外还执行 gauge、边界、Rhie–Chow、canonical 同步或诊断。`PressureCorrection` 还需要矩阵、GlobalDof 和线性求解。通用 `DirectEvaluation` / `LinearEquation` 尚不能表达这些完整 numerical semantics，因此没有拿它们替换生产压力路径。

## 权威边界与过渡状态

`CompiledFormulaStep.fragment` 复用 `SolvePlanNode`；非空 fragment 与单 `backendOperation` 互斥。MethodStepRef 必须恰好匹配一个已编译 ProgramStep 且其 fragment 非空；缺失/重复直接失败，不回退到 formulation stage。工作区依赖由编译步骤的 requires/provides 验证，方法的 numerical requirements 在 `explain` 展示。

**Implemented**：五个 ProgramStep、Output、EquationMethod、fragment、MethodStepRef 精确解析与 fail-fast 检查；constant-density PISO/SIMPLE/PIMPLE 生产计划的方程级叶子来自方法片段。`explain` 同时显示 Formula、HOW step、method fragment、workspace 与 legacy source-math dependency。

**Interface-only**：Formula AST 的一般 operator-to-provider binding；当前 pressure 方法只是把已有 callback 的操作拓扑编译出来，尚不能通用组装任意 pressure Formula 的矩阵。

**Legacy**：`PressureOperators` 仍执行原数值算术；`E_MOMENTUM` 的 `Equation::Definition`/BoundTerm 仍供应动量空间项；outer loop、relaxation、convergence、iteration snapshot 和 commit 仍由 pressure PlanFragment control skeleton 描述。`FormulaMode` 仍服务通用/过渡后端，但已从 pressure coupling source 移除。

**Unsupported**：任意 pressure Formula 直接由 generic FormulaCompiler 取代现有回调、任意新 pressure operator/provider、完全由 structured HOW 生成 PISO/SIMPLE/PIMPLE 宏观循环。这些路径没有被标记成 Runnable 或静默降级。

## 数值与并行不变量

本轮未改 `PressureOperators.cpp`、Rhie–Chow、矩阵系数、gauge、fixed-time candidate、boundary callback、HYPRE 或 MPI 实现。方法片段展开后 PISO 的叶子顺序仍是 `MomentumAssemble → MomentumSolve → PressureBoundaryPrepare → PressureAssemble → PressureSolve → VelocityCorrect → PressureUpdatePrepare → FluxCorrect → CorrectionCommit`；SIMPLE/PIMPLE 仍先冻结 time base、进入 outer iteration，再按原顺序执行这些叶子、relaxation、flux restore、convergence 和 iteration end。Loop count、termination signal 与 physical commit 位置未变。

State/geometry 仍 owner COPY；canonical face 仍 COPY；GlobalDof contributions 仍 SUM。没有增加 Field/Q 副本或第二个压力矩阵。压力边界、gauge 与 canonical 同步实际时机由未改动的生产 callback 保证。

## 验证

迁移前：`cmake --build build --parallel 4` 通过；host 完整 CTest 27/27 通过；architecture checker 通过且保留 2 条既有 allowlist。迁移后：`cmake --build build --parallel 4`、`python3 tools/check_architecture.py`、`git diff --check` 通过；方法单测与 PISO/formulation 架构测试通过。完整 numerical/serial/MPI 回归结果以本报告最终更新为准。

## 后续工作

下一步先把 pressure 的 **macro control skeleton** 迁入 structured HOW，再独立迁移 Formula AST operator occurrence 到 term-provider binding，最后才删除旧 `Equation::Definition` / BoundTerm 双表示。此前不得声称 pressure HOW 或 generic Formula execution 已完全闭环，也不得为了减少 legacy 数量更改已验证的压力算法。
