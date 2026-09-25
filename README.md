# SonicSolver — phase25B

SonicSolver 是科研型结构网格有限差分 CFD 框架。本仓库包含可构建的 CLI、求解器、模型、基础设施、测试案例和架构检查工具。当前版本以**方程组合 → 系统变换 → 数值与执行计划编译 → provider 解析 → 运行**为主线；case 不通过一个顶层“solver family”名称决定全部数学和时间流程。

## 构建与验证

需要 CMake、C++17 编译器和 Python 3；Ninja 可用于下列命令。MPI 与 HYPRE 由 CMake 探测，依赖它们的路径须在可用的 host 环境运行。以下配置只构建 CLI，避免 GUI 的额外依赖：

```sh
cmake -S . -B build -G Ninja \
  -DBUILD_CLI=ON -DBUILD_GUI=OFF -DBUILD_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
python3 tools/check_architecture.py
```

运行或检查仓库中的 case：

```sh
./build/sonicSolver check test/Sod/sodCase
./build/sonicSolver explain test/Sod/sodCase
./build/sonicSolver run --steps 20 test/IBM/cylinderFlowGhost
mpirun -np 4 ./build/sonicSolver run --steps 20 test/Sod/sodCase
```

`sonicSolver --help` 列出 `run`、`check`、`explain`、`doctor`、`models`、`recipes` 和 `init` 等入口。case 输入位于各案例的 `fields/`、`models/`、`solvers/` 和 `mesh/` 目录；运行输出写入案例的 `result/`，已被 Git 忽略。

## 当前架构

```text
内置 preset + model contribution + 用户声明
                    ↓
             RawEquationSystem
                    ↓ formulation / transformation
         ExecutableEquationSystem
                    ↓ state realization + term recipe binding
          CompiledNumericalSystem
                    ↓ PlanFragment + generic SolvePlanner
             CompiledSolvePlan
                    ↓ ProviderResolver
        ResolvedOperationBinding[]
                    ↓ RuntimeRequirements / OpRegistry
               PlanExecutor
                    ↓
              committed state
```

`ExecutableOperation` 声明 OpId 和所需能力，不直接选择具体 provider。`ProviderResolver` 根据已实现的能力及实际 state realization 冻结绑定；无法执行的 operation 明确报告 `Unsupported`。`sonicSolver explain` 展示方程、操作、绑定、计划与运行时要求。详见 [架构文档](src/ARCHITECTURE.md) 和 [Phase 25B 迁移报告](docs/provider-resolution-shared-pressure-closure.md)。

## 支持边界与回归

现有 density single-fluid 的显式时间路径、已实现的 conservative PISO 路径、Eulerian shared-pressure 路径及内置 IBM 数值实现继续由各自现有 kernel 执行。Ghost/ILW 是边界闭合；约束型 IBM 通过已编译 operation 接入执行计划。**Constant-density single-fluid PISO、SIMPLE/PIMPLE dedicated fixed-point predictor 尚未实现 numerical provider**，保留可解释的 `Unsupported`，不会静默退回其他算法。

Phase 25B 验证包括干净构建、14/14 CTest、架构检查，以及 4-rank Sod、serial Ghost、PISO VTS SHA-256 和 Eulerian 两步回归。数值比较的精度与限制见迁移报告；这次架构迁移未修改数值公式或 MPI owner/COPY/SUM/canonical 语义。

源码目录导航见 [src/README.md](src/README.md)。
