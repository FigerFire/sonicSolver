# SonicSolver

SonicSolver 是科研型结构网格有限差分 CFD 框架。case 组合 WHAT（EquationRegistry 数学方程）、STATE（基础变量注册表与物理存储/closure 契约）、HOW（EquationCall → typed Target、scope-local order 与通用 Loop）、WHICH（时间/空间/通量/线性 numerical providers）。Compiler 接收冻结的 WHAT/STATE/HOW/WHICH，验证 capability、绑定 storage 并 lower 到现有 kernels；它不选择 solver。single-fluid density、常密度压力、Eulerian shared-pressure 和 single-fluid/Eulerian RAS 主路径已迁入此模型；专用 numerical kernels 保留，IBM native contribution 与尚未实现的组合能力分别标注。详见 [架构文档](src/ARCHITECTURE.md) 与 [四模块迁移报告](docs/four-module-state-migration.md)。

目标架构正式定义为 **双数学入口、统一方程驱动**：Equation Input 直接贡献方程；规划中的 Variational Input 通过显式变分规则生成同一 WHAT，复用 STATE/HOW/WHICH、Compiler、Runtime Binding 与 PlanExecutor。变分前端尚未实现，安排在 turbulence、IBM、level-set、Eulerian–Eulerian 的选定能力与验证收口之后。详见 [架构与数学边界](src/ARCHITECTURE.md) 和 [后续实施路线](docs/dual-mathematical-input-roadmap.md)。

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

`sonicSolver help` 列出 `run`、`check`、`explain`、`doctor`、`models`、`recipes` 与 `init` 等入口。当前 native case 使用 `case.yaml` 及 `fields/`、`equations/`、`models/`、`solvers/`、`mesh/` 下的 typed 文件；结果写入 case 的 `result/`。查看真实配置示例可从 [常密度 PISO case](test/pressure/constantDensityPiso) 或 [Sod case](test/Sod/sodCase) 开始。

## 当前可运行范围与限制

已有的 density single-fluid 显式时间实现支持 `forwardEuler`、`SSPRK3`、`classicalRK4` 的适配组合；高阶 conservative convection 仍有融合 kernel。常密度单流体 PISO、SIMPLE/PIMPLE 在已验证的正交网格、单逻辑源块能力范围内复用同一组 pressure operators。Eulerian shared-pressure 与内置 IBM 保留经过验证的专用 numerical helpers，通过已编译 operation/Plan 运行。Ghost/ILW 是边界闭合，不是另一套时间循环。

当前并不存在任意用户方程的通用 production lowering；native `replace/disable` 尚未实现，空间 recipe 主要按 convection/diffusion 类别选择，湍流输运和部分相方程仍有专用执行实现。不能从 `explain` 能打印方程就推断任意方程均可运行；最终 provider/capability 绑定决定当前可执行性；`check` 仍有检查覆盖缺口，不能把 `configuration valid` 独立视为 Runnable。Eulerian 湍流已有同等 C++ user contribution 的编译/provider/order 一致性测试，以及原生接线与旧内核的数值对照；这不代表任意 YAML 方程输入已实现。缺口和源码证据集中列在 [审计报告](docs/three-contract-architecture-audit.md)。

当前 authority 与删除项见 [WHAT/HOW/WHICH 迁移报告](docs/equation-execution-migration.md)。较早 Phase 29 报告是历史记录，不再定义 production authority。

Eulerian turbulence 原生迁移已完成：kEpsilon/kOmegaSST 的 phase-mass 加权输运、显式 source/隐式 sink，以及 mu_t 代数闭合通过同一 WHAT/STATE/HOW/WHICH 编译路径执行。STATE aliases 原数组，selected phases 的 subset/反序映射按真实 phase slot 绑定；独立 flow.eulerian-turbulence provider 保留既有 prepare/solve/commit 时序，不由 flow provider 隐式插入。原始 production 入口因 identity/HOW 断点无法运行；本轮与冻结旧内核的隔离接线对照比较，两个模型、多个 PIMPLE passes 及 LES closure 的全部输出/trace 严格相同。完整 host CTest 38/38 通过；Eulerian turbulence MPI/multi-patch 明确 Unsupported。详见[中文报告](docs/eulerian-turbulence-native-migration.md)。

IBM 原生贡献已接入同一 WHAT/STATE/HOW/WHICH：Ghost/ILW 是显式 stage boundary contract；BP 是动量 penalty/机械功；Peskin 与 explicit DFM 保留 lagged force；FTS/fractional DLM 是 post-predictor projection。原数组与数值 kernels 保留，legacy equation/policy/strategy routing 和空 transformer 已删除。KKT 目前有 native structural block contract，但 canonical production cases 仍存在明确 recipe/predictor 能力缺口；除后述显式静止黏性 Ghost＋SST 受限组合外，其他 turbulence＋IBM、pressure＋IBM、Eulerian＋IBM 未自动开放。Peskin MPI 的合法局部零支撑现在允许进入全局 SUM，完整 diagonal 仍严格检查；物理 ILW 的切分法向、切向 halo 输入，以及 force/mask 输出 COPY 已修正。上一冻结版本 47/47 主机 CTest 通过，7 个串行案例与冻结输出逐字节一致，8 个 MPI 案例完整运行到 endTime=5；串行/MPI 全时域最大速度分量差异 6.75e-13，输出 replica 冲突为 0。物理体积指标及旧 partitionOfUnity 的 surface 一阶矩/力矩缺口见[并行验证报告](docs/ibm-parallel-consistency-validation.md)。迁移阶段的历史同配置对照见[中文 IBM 迁移报告](docs/ibm-native-contributions-migration.md)。

Surface transfer 现在可以显式选择 `surfaceOperator.normalization: linearReproducing`，在原 WendlandC2 支撑与 adjoint spreading 上同时重现常量与坐标一阶矩，使 marker/Eulerian force、torque、power 使用一致的离散交换。旧 `partitionOfUnity` 继续保留；新选择会改变数值方法，不会隐式启用。当前要求 full-rank 3D affine support，退化支撑 fail-fast，修正权重允许有符号。50/50 主机 CTest 通过；15 个旧选择 serial/MPI controls 的完整输出及 time/dt 日志序列保持冻结结果，6 个新选择 serial/MPI-2/MPI-4 案例完整运行到 endTime=5，并验证一阶矩、力矩/功率配对和实际角动量交换。完整时域场对照及能力边界见[一阶矩／力矩报告](docs/ibm-first-moment-report.md)；该 transfer 工作本身不开放 SST viscous wall 或 production KKT。

Level-set 现已原生表达 HJ-WENO 物理输运、冻结 phi0、伪时间 Loop、normal/curvature 与 CSF/jump 关系。伪时间不推进物理钟。新增 native `executionProgram` / `providerBindings` 允许显式配置已有 HOW/WHICH；未知或未实现契约 fail-fast，任意 YAML PDE 仍未支持。Mixture/homogeneous 已移除内置 flat producer，但 numerical provider 能力仍明确 Unsupported。串行、MPI-2/4 同配置旧内核输出保持；MPI level-set 仍有旧分区边界 phi/派生场 replica 差异，不能将保持结果称为并行物理一致性证明。范围、53 项测试和证据见[本轮报告](docs/native-remaining-migration-report.md)。

跨模块组合现在由 scope-local placement、provider capabilities 和版本化 STATE 契约验证。显式 HOW 数组顺序保持；默认只在既有 recipe 主干中按相对要求插入模块，歧义要求显式选择。验证区分 current、stage、old-time、outer 冻结系数与上一步 lagged state；U/p/T 惰性缓存保留，STATE 不执行隐藏刷新。IBM 的本次选择冻结一次，provider 独立核对数学支持范围。移除模块组合黑名单并不开放缺少 immersed turbulence wall/boundary 数学的 RAS＋IBM；后述新壁面选项有独立数学与验证。实现与验收范围见[跨模块组合报告](docs/generic-cross-module-composition.md)。


## 静止黏性 Ghost 与 SST 的受限组合

现在可显式选择 `IBM.wallClosure: stationaryNoSlipAdiabatic`。它采用约束二次拟合实现静止无滑移速度和绝热热边界，与默认 `eulerSlip` 分开；必须关闭 Euler ILW，并启用 Central2 黏性项。支撑不满秩或 rho/p 非物理时失败，不自动降阶或切回 slip。网格节点恰好落在壁面上时按边界点施加零速度。

在 single-fluid、串行、单 patch 的显式 flow 中，这个选项提供 SST 的 k=0、网格相关 omega 壁值、mu_t=0 和原几何 wall-distance view。SST 每物理步读取 boundary-ready 速度并更新一次，新组合的 mu_t 在 flow stages 中冻结。方程、fluid port、壁面三者必须匹配，不能由其他壁面的同名 capability 代替。

可运行输入见 `test/IBM/cylinderFlowViscousGhost`、`cylinderFlowSSTGhost` 和 `cylinderFlowSSTGhostAuthored`。默认/显式 HOW 使用同一 provider；2000 步输出与操作轨迹一致。解析边界 profile 和实际 Central2/Newtonian 壁面剪切有网格收敛验证。**现有 single-fluid SST 内核仍没有对流项；本轮没有证明充分发展的湍流近壁剖面、y+ 或长期壁面预测精度。** MPI、移动壁、其他 IBM 的 SST 壁面、pressure/Eulerian 组合仍未开放。实现和证据边界见 [本轮报告](docs/viscous-immersed-wall-sst-native-composition.md)。


命令统一使用无前导横杠的子命令，例如 `sonicSolver postProcessing [CASE]`、
`sonicSolver cleanResult [CASE]`、`sonicSolver initialOutput [CASE]`、
`sonicSolver steps N [CASE]` 和 `sonicSolver explainModel MODEL`。
命令参数继续使用 `--recipe`、`--with` 等选项；顶层命令只接受上述统一拼写。

终端环境采用可 source 的 `etc/bashrc`，同时配置构建目录的 PATH 和 Bash/Zsh 补全。
在 `~/.bashrc`（Bash）或 `~/.zshrc`（Zsh）中加入以下一行，路径换成实际 checkout：

```sh
source "/Volumes/SSD_LPF/sonicSolver/sonicSolver/etc/bashrc"
```

当前终端执行同一行立即生效；以后新终端自动加载。Bash 登录终端需在
`~/.bash_profile` 中加载 `~/.bashrc`。输入 `sonicSolver po` 后按 Tab 会补全为
`sonicSolver postProcessing`，case 参数也支持目录补全。重复 source 不会重复添加 PATH。
构建脚本会打印环境加载命令；它不修改用户的 shell 启动文件。

`postProcessing` 在 macOS 通过 LaunchServices 打开 ParaView，启动成功后立即返回终端，
ParaView 独立运行。其他平台等待 `paraview` 退出并返回其退出状态；启动失败均报错。
