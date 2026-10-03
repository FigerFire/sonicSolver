# src 源码导航

`src/` 以 WHAT（EquationRegistry）、STATE（base variables / lazy compiled views）、HOW（ordered EquationCall / generic scopes）、WHICH（Numerics/providers）组织。当前实现和 legacy 边界见 [ARCHITECTURE.md](ARCHITECTURE.md) 与 [四模块迁移报告](../docs/four-module-state-migration.md)。

| 目录 | 主要职责 |
|---|---|
| `core/` | Field、StateBundle、强类型配置、运行接口及中立 `core/system` 值类型 |
| `models/` | 物理源项、湍流、IBM、相模型；向系统提供 contribution，不拥有全局 timestep |
| `solver/system/` | 三契约的 composition、transformation、state/term compilation、SolvePlanner、ProviderResolver 和验证；尚未通用化任意方程 |
| `solver/algorithm/` | 现有数值 provider、solver workspace、显式 stage 运算和 Eulerian 相执行；部分专用 numerical helper 仍在迁移 |
| `solver/equation/`、`discretization/`、`boundary/`、`linearAlgebra/` | 方程装配、空间离散、边界闭合与线性后端 |
| `infrastructure/` | IO、网格和 MPI/backend 实现 |
| `methods/` | 可复用数学/几何方法及仍在迁移的数值工具 |
| `app/` | case 解析、application composition、CLI/GUI、输出 |

当前重点入口：

- `core/system/SF_equationIR.h`：Raw/Executable equation system 与 executable operation 的中立描述。
- `core/system/SF_systemContribution.h`：模型和用户可共享的方程贡献值；native 自定义方程还不能完整降低。
- `core/config/types/SF_timeRecipe.h`、`SF_termRecipe.h`：当前时间与空间 recipe；RK2/隐式及逐 term 自由选择尚未实现。
- `core/system/SF_solveProgram.h`：纯 source HOW 与独立 CompiledTarget/plan；`SF_planFragment.h` 仅保留显式 LegacyPlanFragment。
- `solver/system/SF_systemBuilder.cpp`：collect modules 与 composition 入口；`SF_numericalCompiler.cpp::compileSystem` 负责 target、numerics 与 lowering。
- `solver/system/SF_numericalSelection.h`：source WHICH；`SF_builtinProviders.cpp`：domain capabilities、storage 与 numerical lifecycle；`SF_providerResolver.cpp`：builtin runtime provider matching。
- `solver/system/SF_solvePlan.cpp`：通用 plan lowering；压力与 Eulerian 的顺序来自 contribution/fragment。
- `solver/run/SF_planExecutor.cpp`：执行结构化计划，不选择物理方程。
- `solver/algorithm/SF_singleFluidStepper.cpp`、`eulerian/SF_eulerianStepper.cpp`：绑定已分配的 operation 并调用现有数值 kernel。

模型贡献依赖 `core/system`，不能为了取得中立 IR 反向依赖 `solver/system`。预设在组合阶段展开，运行时不能按预设名选择另一个方程执行入口。single-fluid 调用现有 numerical kernels；Eulerian、湍流和 IBM 专用 backend 保留明确的 migration boundary。并行状态和几何使用 owner→COPY，共享面只有一个 canonical flux，残差/源项/载荷使用 SUM；详见架构文档。
