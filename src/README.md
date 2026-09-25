# src 源码导航 — phase25B

`src/` 按职责分为六个域。完整的编译与运行契约见 [ARCHITECTURE.md](ARCHITECTURE.md)。

| 目录 | 主要职责 |
|---|---|
| `core/` | Field、StateBundle、强类型配置、运行接口及中立 `core/system` 值类型 |
| `models/` | 物理源项、湍流、IBM、相模型；向系统提供 contribution，不拥有全局 timestep |
| `solver/system/` | 方程组合、变换、state realization、numerical compiler、SolvePlanner、ProviderResolver 和验证 |
| `solver/algorithm/` | 现有数值 provider、solver workspace、显式单 stage 运算和 Eulerian 相执行 |
| `solver/equation/`、`discretization/`、`boundary/`、`linearAlgebra/` | 方程装配、空间离散、边界闭合与线性后端 |
| `infrastructure/` | IO、网格和 MPI/backend 实现 |
| `methods/` | 可复用数学/几何方法及仍在迁移的数值工具 |
| `app/` | case 解析、application composition、CLI/GUI、输出 |

当前重点入口：

- `core/system/SF_equationIR.h`：Raw/Executable equation system 与 executable operation 的中立描述。
- `core/system/SF_planFragment.h`、`SF_solveProgram.h`：控制流片段与编译计划值类型。
- `solver/system/SF_systemBuilder.cpp`：启动阶段编译入口。
- `solver/system/SF_providerResolver.cpp`：OpId 到 numerical provider 的唯一编译绑定。
- `solver/system/SF_solvePlan.cpp`：通用 plan lowering；压力与 Eulerian 的顺序来自 contribution/fragment。
- `solver/run/SF_planExecutor.cpp`：执行结构化计划，不选择物理方程。
- `solver/algorithm/SF_singleFluidStepper.cpp`、`eulerian/SF_eulerianStepper.cpp`：绑定已分配的 operation 并调用现有数值 kernel。

模型贡献依赖 `core/system`，不能为了取得中立 IR 反向依赖 `solver/system`。并行状态和几何使用 owner→COPY，共享面只有一个 canonical flux，残差/源项/载荷使用 SUM；详见架构文档。
