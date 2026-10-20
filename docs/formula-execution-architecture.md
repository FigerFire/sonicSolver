# Formula / Execution / Numerical Method 架构契约

日期：2026-09-29。本文描述 Phase 28D.1 的目标契约及当前可运行边界。

## 三个契约

| 契约 | 问题 | 中立表示 | 不得包含 |
|---|---|---|---|
| WHAT | 数学上是什么 | `Formula` 等式 AST、`FormulaGroup` 引用集合、Unknown/Constraint/Closure | output、时间方法、执行模式、顺序 |
| HOW | 哪个步骤产生什么，何时执行 | `ProgramStep {subject, OutputRef}`、`ExecutionProgram`、循环/提交结构 | 显式/隐式标志、RK 系数、WENO/HYPRE 名称 |
| HOW NUMERICALLY | 如何实现步骤与算子 | `ITemporalMethod`、`IEquationMethod`、空间 provider 与后端绑定 | 新的物理方程 authority、全局 lifecycle |

`OutputRef` 显式给出 primary output 与其符号，必要时给出 numerical workspace。额外的 `rAU`、`HbyA`、corrected face flux 是方法/provider 的派生产物，不能成为另一份 physical state。编译器只解析引用、验证方法能力、调用方法对象并组合 plan。它不能从 `ddt` 或左侧符号猜测输出，不能从 Formula ID 或 EquationRole 猜显式/隐式，也不能自己发明离散 stencil。

`FormulaCall {formula,target,mode}` 是现有后端的过渡 IR。`mode` 由 Equation Method 降低得到，用户层不声明 `explicit::X` 或 `implicit::X`。显式/隐式是方法的数值属性：时间推进是否显式，与空间算子是否需要矩阵装配是两个不同维度；代数赋值也不是显式时间积分。

## 当前 lowering

```text
Flow = {E_MASS, E_MOMENTUM, E_ENERGY}       WHAT
Flow -> Q{rho,rhoU,rhoE}                    HOW
TemporalMethod = FE / SSPRK3 / RK4          HOW NUMERICALLY
EquationMethod = ConservativeResidual      HOW NUMERICALLY
term providers = configured convection / diffusion / source
                 ↓
compiled stage FormulaCalls + CompiledTimeRecipe
                 ↓
existing fused conservative RHS and stage backend
```

`Flow` 是普通公式组，不持有 Q 或 timestep。`ConservativeResidual` 验证 group/output arity 和每条公式左侧对对应输出的 `ddt`，再形成过渡 `FormulaCall`；没有按 E_MASS 或 EquationRole 分派。三个 temporal method 对象分别拥有经验证的 stage 系数。生产 backend 仍按固定 conservative slices 做严格能力检查，不支持的 group fail-fast；这不是任意 Formula 的生产执行能力。

压力约束 transformation 已产生 `F_PRESSURE_UPDATE`、`F_VELOCITY_CORRECTION`、`F_FLUX_CORRECTION`。它们当前映射为 `DirectEvaluation`、输出分别为 p、U、phi；phi 显式绑定 pressure-face-flux workspace。Plan 的实际压力更新仍调用原有 `PressureOperators` 回调；`DirectEvaluation` 的运行时替换尚未完成。Pressure matrix、rAU/HbyA、Rhie–Chow、boundary 与 gauge 没有被通用编译器重算。

## 验证和限制

- `ProgramStep` 必须声明 subject 与 primary Output；未知方法、缺失存储、FormulaGroup arity 不符、错误的 `ddt(output)` 均报错。
- 一个 WHAT/HOW 可以分别选择 ForwardEuler 或 ClassicalRK4；stage 数学不同，Formula 和 Output 不变。
- `LinearEquation` Equation Method 能把 `DiffuseC -> C` 降到线性装配调用；已有 FormulaCompiler/Central2/GlobalDofSystem 测试用这个方法输出实际求解任意命名 scalar。用户层不提供 `Implicit` 标志。
- `FormulaCompiler` 独立小测试仍覆盖任意 scalar/vector Central2 线性扩散和 GlobalDofSystem；不等于 production 压力或高阶通量已迁移。
- Runtime 仍由 `PlanExecutor`/`OpRegistry` 执行已绑定 provider。canonical face COPY、residual SUM、stage time 与 boundary/halo 顺序不得由三契约重定义。
