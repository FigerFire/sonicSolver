# SonicSolver

SonicSolver 是科研型结构网格有限差分 CFD 框架。case 组合 WHAT（EquationRegistry 数学方程）、STATE（基础变量注册表与物理存储/closure 契约）、HOW（EquationCall → typed Target、scope-local order 与通用 Loop）、WHICH（时间/空间/通量/线性 numerical providers）。Compiler 接收冻结的 WHAT/STATE/HOW/WHICH，验证 capability、绑定 storage 并 lower 到现有 kernels；它不选择 solver。single-fluid density 与常密度压力主路径已迁入此模型；Eulerian、湍流与 IBM 的 legacy backend 边界明确保留。详见 [架构文档](src/ARCHITECTURE.md) 与 [四模块迁移报告](docs/four-module-state-migration.md)。

| 契约 | 用户选择的含义 | 当前对应物 |
|---|---|---|
| WHAT — Equations | 稳定数学定义与 AST | `EquationRegistry`、`Equation`、`SystemContribution` |
| HOW — Execution | `equation -> target`，scope-local order、循环、提交 | `EquationCall`、`Target`、`ExecutionScope` |
| WHICH — Numerics | FE/SSPRK3/RK4、ConservativeResidual/DirectEvaluation/LinearEquation、WENO/TENO/Central 等 | `NumericalSelection`、`TermRecipe`、provider catalog |

`incompressible` 在目标设计中是展开常密度闭合与不可压缩约束的**预设名**，不是顶层 solver family 开关。不可压缩 laminar 不输运密度，也不增加湍流方程，但仍必须满足 `div(U)=0`；压力是约束乘子，SIMPLE/PISO/PIMPLE 是可选的约束求解策略。`kEpsilon` 则应增加 `k`、`epsilon` 方程和动量应力闭合，不应启动独立的 `KEpsilonSolver`。

目标用法允许“预设”“预设 + add/extend/replace/disable”和“完全自定义方程”最终得到等价的 compiled system。**这是设计验收标准，不是当前 native 输入已支持的语法。** 例如 `incompressible + navierStokes + kEpsilon` 配 `RK2 + SIMPLE` 和 `WENO5 + Central2` 目前不能作为可运行组合照抄；RK2、用户方程编译及该组合的 stage/数值 provider 尚未闭环。现有输入格式请以 `test/` 下的真实 case 为准，使用 `sonicSolver check`、`explain` 核实当前 capability。

## 构建与检查

需要 CMake、C++17 编译器、Python 3；MPI/HYPRE 路径取决于构建环境。下列配置只构建 CLI：

```sh
cmake -S . -B build -G Ninja -DBUILD_CLI=ON -DBUILD_GUI=OFF -DBUILD_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --parallel 1
python3 tools/check_architecture.py
```

仓库中的实际 case 可以这样检查和运行：

```sh
./build/sonicSolver check test/Sod/sodCase
./build/sonicSolver explain test/pressure/constantDensityPiso
./build/sonicSolver run --steps 20 test/IBM/cylinderFlowGhost
mpirun -np 4 ./build/sonicSolver run --steps 20 test/Sod/sodCase
```

`sonicSolver --help` 列出 `run`、`check`、`explain`、`doctor`、`models`、`recipes` 与 `init` 等入口。当前 native case 使用 `case.yaml` 及 `fields/`、`equations/`、`models/`、`solvers/`、`mesh/` 下的 typed 文件；结果写入 case 的 `result/`。查看真实配置示例可从 [常密度 PISO case](test/pressure/constantDensityPiso) 或 [Sod case](test/Sod/sodCase) 开始。

## 当前可运行范围与限制

已有的 density single-fluid 显式时间实现支持 `forwardEuler`、`SSPRK3`、`classicalRK4` 的适配组合；高阶 conservative convection 仍有融合 kernel。常密度单流体 PISO、SIMPLE/PIMPLE 在已验证的正交网格、单逻辑源块能力范围内复用同一组 pressure operators。Eulerian shared-pressure 与内置 IBM 保留经过验证的专用 numerical helpers，通过已编译 operation/Plan 运行。Ghost/ILW 是边界闭合，不是另一套时间循环。

当前并不存在任意用户方程的通用 production lowering；native `replace/disable` 尚未实现，空间 recipe 主要按 convection/diffusion 类别选择，湍流输运和部分相方程仍有专用执行实现。不能从 `explain` 能打印方程就推断任意方程均可运行；`check` 的 provider 与 capability 结果才是当前可执行性的依据。内置预设与完全等价的用户显式方程尚无编译产物和数值结果相同的验收测试。缺口和源码证据集中列在 [审计报告](docs/three-contract-architecture-audit.md)。

当前 authority 与删除项见 [WHAT/HOW/WHICH 迁移报告](docs/equation-execution-migration.md)。较早 Phase 29 报告是历史记录，不再定义 production authority。
