# Phase 28D — Formula WHAT / HOW / Numerical Binding 实施记录

> 2026-09-29 更新：本文保留 Phase 28D 当时的接口与数值边界记录。其中 `FormulaCall.mode` 作为 HOW 的说法已被 Phase 28D.1 取代：用户层 HOW 现在是 `ProgramStep: Formula/FormulaGroup -> Output`，时间/方程/空间方法在 HOW NUMERICALLY 选择。当前实现与剩余 legacy 见 [Phase 28D.1 报告](phase28d1-method-object-architecture.md) 和 [三契约说明](formula-execution-architecture.md)。下文历史状态不代表当前 authority。

日期：2026-09-28。此报告区分已投入 production 的执行路径与仅在通用接口测试中运行的路径；**本阶段尚未达到完整 Formula execution authority closure**。

## Before

`Equation::Definition` 的平面 term、`NumericalRecipeSet`、`ExecutionPolicy` 和各 operation provider 分别描述方程、离散及实际执行。显式 Mass/Momentum/Energy 由 fused conservative kernel 装配；常密度 PISO/SIMPLE/PIMPLE 由同一 pressure operation provider 执行。`PlanExecutor` 决定顺序，但叶子未保存可独立编译的 Formula 调用。因此 equation/term metadata 不能直接作为任意公式的通用执行输入。

## After：三个正交输入

| 契约 | 主要类型与路径 | 当前作用 |
|---|---|---|
| WHAT | `src/core/system/SF_formula.h` 的 `FormulaExpr`、`Formula`、`FormulaRegistry` | 存储纯数学等式 AST、ID、来源；不存 target/mode/order。支持 add/replace/disable 的显式语义。 |
| HOW | `src/core/system/SF_solveProgram.h` 的 `FormulaCall`，`SF_planFragment.h` / `SF_solvePlan.cpp` | Plan 叶子保存 Formula ID、目标与 Assign/Explicit/Implicit 模式。显式 stage 明列 rho/rhoU/rhoE；常密度 pressure 计划明列 predictor、pressure correction、velocity/pressure/flux update。 |
| NUMERICAL | `src/solver/system/SF_formulaCompiler.*`、`src/solver/discretization/SF_formulaOperator.h` | 按 Formula/算子位置解析 provider；全局默认 < Formula 默认 < 位置覆盖。缺失、多义、重复位置或指向不存在位置的绑定报错。 |
| Runtime | 现有 `Run::PlanExecutor` + `OpRegistry` | 保持已验证 production 数值路径；独立测试中的通用 FormulaCall 可通过同一 PlanExecutor 调用。 |

`SystemContribution::addFormula` 使内置与用户 C++ Formula 都能经 `SystemCompositionBuilder` 写入同一 RawEquationSystem FormulaRegistry；Formula ID 和 provenance 不参与 `FormulaCompiler` 后端分派。旧 `Equation::Definition` 在迁移期通过 `formulaFromEquation` 生成 AST。**Raw/Executable 中暂时同时保留平面 equationDefinitions 和 FormulaRegistry**：前者仍被 production numerical compiler/旧 kernel 使用，不能把这对表示说成已经收敛为单一 authority。Native 输入的自定义 Formula lowering 尚未闭环。

## Formula 编译语义

- `Assign`：只允许直接左侧符号 `X = expression`；对全部 owned cell/component 先计算快照，再写回，避免同一步读写串扰。它是无 explicit/implicit 前缀的代数赋值语义。
- `Explicit`：当前支持左侧恰有一个正号 `ddt(target)` 的方程，编译成 `rhs - lhsWithoutDdt` 的导数评价；stage 系数与 stage time 继续由现有 `CompiledTimeRecipe` 管理。算子数值由 provider 给出。
- `Implicit`：编译器只进行目标依赖检查与仿射 AST 组合。`div(nu*grad(X))` 的 Cartesian Central2 provider 返回每行的 stencil 系数及 Dirichlet 边界常数；生成 `GlobalDofSystem`，交由现有 linear algebra 求解并写回。非线性目标乘积、目标出现在分母、没有线性能力的 provider 均 fail-fast。**generic implicit = generic compilation/assembly contract, not arbitrary nonlinear formula solver.**
- `rAU`、`HbyA`、`faceResponse`、Rhie–Chow 属于 pressure numerical workspace/provider；Formula 可引用这些符号或数学 primitive，但通用 compiler 不自行发明其计算方式、边界处理或压力 gauge。

## Backend / capability

| Formula 形态 | Explicit/赋值 | Implicit | 实际后端与状态 |
|---|---|---|---|
| density Mass/Momentum/Energy | Implemented | Unsupported | Plan 有 FormulaCalls；production 仍由既有 fused conservative operation provider 执行，尚非通用 FormulaCompiler lowering。绑定时要求三项调用及其 AST 与 backend 实际读取的定义严格一致。 |
| 代数赋值 `C=C+2` | Implemented | 不适用 | 通用 AST expression evaluation，经 PlanExecutor 测试。 |
| scalar `-div(nu grad(Y_O2))=S` | 可评价 | Implemented | 通用 Central2 provider + GlobalDofSystem/HYPRE，3×3 Dirichlet 解测试通过。 |
| vector `-div(nu grad(V))=S` | 可评价 | Implemented | 同一 provider/component 循环、两分量线性解测试通过，无 VectorSolver 分支。 |
| 常密度 pressure correction | Legacy | Legacy | Plan 标注 `implicit::pPrime`；实际 matrix/Rhie–Chow/gauge 仍在 `PressureOperators`。 |
| WENO/TENO implicit | Unsupported | Unsupported | 没有隐式 stencil/linearization provider；通用编译器明确报缺少 linear assembly capability，不回退显式或 Rusanov。 |
| 混合对流+扩散 passive scalar | 独立 C++ 测试 Implemented | Interface-only | 测试用显式 upwind provider 与 Central2 provider 经同一 AST 组合为 `ddt(C)+div(phi,C)-div(nu grad C)=S`；production storage、通量 provider、边界与时间配方尚未接通。 |

## Pressure 与 legacy 边界

`src/solver/system/SF_transformation.cpp` 为现有常密度 pressure transformation 加入 predictor、pressure equation 与三项更新 Formula：`U=HbyA-rAU·grad(p)`、`p=p+pPrime`、`U=U-rAU·grad(pPrime)`、`phi=phi-pressureFlux(...)` 等。`PressureCoupling` 将它们引用到已有 PISO/SIMPLE/PIMPLE 结构化计划。其实际 `HbyA`、`rAU`、Rhie–Chow、pressure reference、边界及修正顺序没有改变，仍由 `src/solver/algorithm/pressure/SF_pressureOperators.cpp` 的 operation callbacks 实现。`explain` 新增 Formula 与 FormulaCall execution 状态，明确标注 operation-bound production provider，避免把 Plan 注解误当作通用 compiler 已接管数值执行。

保留为 Legacy 的还有 `src/solver/algorithm/SF_conservativeRHS.cpp`、Eulerian transport、turbulence、level-set 和 IBM 的已有 numerical helper。Ghost/ILW 继续是 boundary closure；本阶段不触及 DLM/KKT、MPI 或 stage 数学。

## 验证与数值不变量

- `cmake --build build --parallel 4`：完整编译与链接成功。
- `python3 tools/check_architecture.py`：通过；已知 allowlisted dependency edges 仍为 2。增加了 Formula WHAT 不能 include MPI/HYPRE/solver implementation，以及通用 compiler 不能按 PISO/SIMPLE/kEpsilon/Momentum 特判的守卫。
- `formulaExecutionContract`：任意命名 scalar、vector 扩散实际 linear solve；赋值快照；显式对流+扩散 scalar residual；implicit WENO fail-fast；数值绑定三级优先级及坏路径拒绝；内置/用户 Formula 进入同一 composition 和 compiler，得到相同 backend/provider/数值导数。
- 最终 host 完整 CTest：26/26 通过，包括 Sod、PISO、SIMPLE/PIMPLE fixed-time、Poiseuille、MPI pressure primitive/fail-stop/regression。受限 sandbox 中 PRTE socket bind 会阻止 `haloFailureContract` 等 MPI 测试启动；host 重跑全部通过。fused FormulaCall 严格绑定及坏路径检查完成后，已重新完整构建和运行这轮完整测试。
- serial Ghost IBM 两步：WENO5/Lax–Friedrichs/RK4/ILW 路径完成；`dt=4.392680e-4, 4.392260e-4`，最终 `min(rho)=1.22285`、`min(p)=101109`。
- 4-rank Sod 两步：TENO5/Steger–Warming 路径完成；`dt=6.681531e-4, 5.575894e-4`，最终 `min(rho)=0.273944`、`min(p)=27143.7`。

本阶段没有修改 WENO/TENO、flux splitting、RK coefficients/stage time、pressure matrix/Rhie–Chow、boundary/halo 顺序、canonical face COPY、GlobalDof SUM 或 HYPRE ownership。对新通用扩散 provider 的验证是独立小网格测试，**不是**旧生产格式的数值等价证明。

## 下一步与未完成标准

要完成整份 Phase 28D 的 authority closure，还需要：

1. 将 production conservative stage 的 fused backend 真正绑定到已编译的 FormulaCall group 和稳定 operator occurrence，而非只校验调用列表后继续使用原有方程入口；保持 canonical F* 与 frozen 数值结果。
2. 将 pressure 现有 numerical helper 分解成可编译的 Formula primitive/provider 与必要的 boundary/gauge operations；先证明首差异，再切换，不重写 pressure 数学。
3. 给通用 Formula 编译结果补 read/write/workspace、halo/canonical/reduction 能力，接入真实 symbol/storage view；当前 `FormulaValues` 只在独立测试中绑定。
4. 打通 native 用户 Formula 和 add/extend/replace/disable 到同一 Raw/Executable/Plan 链，迁走平面 `equationDefinitions` 的剩余 production authority 后才能删除双表示。
5. 加入完整 scalar convection+diffusion、实际 SIMPLE/PIMPLE Formula lowering、Eulerian/IBM/turbulence/level-set 的分阶段数值回归。当前不应把这些组合标为通用 Formula execution 已实现。
