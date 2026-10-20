# Phase 28D.1 — 三契约细化与方法对象编译报告

日期：2026-09-29。

后续 Phase 28D.2 已将 production HOW 转为结构化 root，让 TemporalMethod 编译 density 时间片段，并增加单 patch Field/StateRealization 的通用公式绑定。本文保留 28D.1 当时的验证记录；当前 authority 与剩余 Legacy 以 [28D.2 报告](phase28d2-method-runtime-authority.md) 为准。

## 已实现

WHAT 的 `FormulaGroupRegistry` 可将多条已有 Formula 组成无执行语义的命名组。单流体预设把 Mass/Momentum/Energy 组成 `Flow`。Raw 到 Executable transformation 保留该组；模型和用户 `SystemContribution` 可走同一组注册入口。

HOW 新增 `ProgramStep`、`OutputRef`、`ExecutionProgram`。单流体 production composition 显式给出 `Flow -> Q{rho,rhoU,rhoE}`。Output 不由 Formula AST、Formula ID、EquationRole 或求解器家族推断。既有 `FormulaCall.mode` 只保留为编译后端过渡字段。

HOW NUMERICALLY 新增独立 `ITemporalMethod` 和 `IEquationMethod` 及轻量注册表。ForwardEuler、SSPRK3、ClassicalRK4 各自生成原有 stage 系数；`ConservativeResidual` 验证 formula-group/output 对应关系并生成 stage 调用，`DirectEvaluation` 验证赋值左侧与 Output，`LinearEquation` 声明线性装配、GlobalDof 与求解后端要求。独立 scalar 扩散测试已经让 `DiffuseC -> C` 经过 `LinearEquation` 再进入现有 FormulaCompiler/Central2/GlobalDofSystem。缺方法、缺 output、缺 storage、错配组长度、缺 `ddt(output)` 均立即失败。方法编译结果记录读/写和最低要求；production 后端数值仍由原有实现执行。

数值方法绑定必须恰好对应一个被执行的 ProgramStep；重复绑定和未被任何步骤消费的绑定会直接报错，不能静默忽略用户选择。

新架构测试位于 `test/test_methodObjects.cpp`，现有通用公式测试 `test/test_formulaExecution.cpp` 也接入了 `LinearEquation`。仓库默认忽略整个 `test/`，因此 `.gitignore` 只对这两个尚未跟踪的测试源文件开放例外；`docs/` 仍按现有约定保持本地忽略。

`SolvePlanner` 的 single-fluid explicit stage 现在消费 `CompiledExecutionProgram` 生成的 FormulaCalls；不再从 EquationDescriptor 的 solvedUnknowns 推断这条生产路径的调用。现有 fused kernel 仍检查固定三分量 conservative layout 与 Formula/旧 Equation 定义的一致性，保护数值实现边界。`Flow -> Q` 与时间方法分别绑定，同一 WHAT/HOW 可在 FE 与 RK4 间切换。

常密度压力 transformation 的 PressureUpdate、VelocityCorrection、FluxCorrection 已形成 `FormulaStep -> Output` + `DirectEvaluation` 分类；`phi` 的 workspace 由 pressure numerical provider 确认。其数值执行仍是原 pressure Plan 与 `PressureOperators` 回调，不宣称 DirectEvaluation 已替换压力 kernel。

## 数值语义和未完成

本轮未改 WENO/TENO、通量分裂、conservative RHS、RK stage 更新公式、压力矩阵、Rhie–Chow、边界/halo、canonical COPY、GlobalDof SUM、IBM。方法对象目前编译 stage tableau 与公式调用，不直接执行全部时间数学；`Time::Explicit` 内的 FE/SSPRK3/RK4 分支仍是 transition backend。`FormulaCall` 和旧 Equation::Definition 也尚未删除。

Interface-only：自定义 Temporal Method 注册后的 production stage backend 扩展、通用 PressureCorrection/LinearEquation Method、任意 Formula 的 production storage/provider lowering、native 用户 Formula 输入。Legacy：pressure update callbacks、固定 conservative storage slices、旧 Equation::Definition 与 Formula 双表示。Unsupported：未提供相应 method/backend 能力的组合，禁止 fallback。

## 验证

| 验证 | 结果 |
|---|---|
| `cmake --build build --parallel 4` | 完整编译、链接通过；最后加入 LinearEquation 后再次增量构建通过 |
| `python3 tools/check_architecture.py` | 通过；已知 allowlisted dependency edges = 2，无新增方向违规 |
| Host `ctest --test-dir build --output-on-failure --parallel 1` | 27/27 通过，含压力 PISO/SIMPLE/PIMPLE、Poiseuille、MPI pressure regression、原有 IBM/数值测试 |
| LinearEquation 补充后的定向 CTest | `formulaExecutionContract`、`explicitStageMathematics`、`methodObjectCompilation`、`formulationArchitecture` 4/4 通过 |
| serial Ghost IBM，两步 | WENO5/Lax–Friedrichs/RK4/ILW；`dt=4.392680e-4, 4.392260e-4`；最终 `min(rho)=1.22285`、`min(p)=101109`，与 Phase 28D 记录一致 |
| 4-rank Sod，两步 | TENO5/Steger–Warming/ForwardEuler；`dt=6.681531e-4, 5.575894e-4`；最终 `min(rho)=0.273944`、`min(p)=27143.7`，与 Phase 28D 记录一致 |
| `git diff --check` | 通过 |

首次 CTest 发现 `phi` workspace 只在 pressure numerical provider 中声明，以及压力 DirectEvaluation 步骤不属于时间积分方程块。这两处迁移边界已修复并重新测试。随后 `resolvedSystemFieldCount` 守卫发现向 `ResolvedSimulationSystem` 添加了三个聚合成员；现已将 source/compiled program 放进既有 `CompiledSolvePlan` typed product，恢复原来的 11 个顶层字段，守卫通过。

完整 CTest 在最后新增 `LinearEquation` 方法前已通过。此方法仅用于独立公式测试，后续重建与上述 4 项定向测试再次通过；production 数值路径没有在该补充后变更。

数值差异必须定位首个不同 checkpoint，不能修改容差或数值公式掩盖差异。本轮没有进行 Git 提交或上传。
