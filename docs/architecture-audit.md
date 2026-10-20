# src 架构审计

审计日期：2026-09-11。依据：根目录 `AGENTS.md` 与 `ARCHITECTURE.md`（下文简称 A）。本次只生成本文，不修改源代码、API、CMake 或目录结构。

## 1. 范围、方法与结论

扫描 `src/` 全部 724 个文件，检查所有 `.h/.hpp` 同名项、源码 include、实现文件行数、legacy 关键词和状态入口；对关键算法、状态和服务路径人工追踪。使用已有 `build/compile_commands.json` 的 249 条编译记录，按当前文件目录和 `-I` 顺序递归核对了 1325 个项目内 include 位置；同时用全量词法扫描覆盖未进入该构建配置的文件。编译数据库可能反映已有构建配置，不代表所有平台和预处理宏组合；未运行预处理器或重新配置构建。

结论：目录主骨架已形成，但尚不符合完整架构约束。确认 24 条重点跨层 include（其中两条是外部线性代数后端实现与接口的归属冲突，不是 CFD 算法反向控制）。另有 algorithm 内实现离散公式、methods 内承载 CFD 算法、models 控制温度推进、重复状态入口及双套显式 lifecycle。发现 17 组同名头文件、13 个至少 800 行的 `.cpp` 和 3 个至少 800 行的头文件。

风险定义：P0 = 已确认的严重数值错误或不可运行；P1 = 可造成状态、算法或并行语义分歧的架构问题；P2 = 依赖、维护或潜在 include 歧义。本次没有足够运行证据认定 P0。下文“可能改变数值”描述后续处理的敏感性，本文本身不改变任何数值结果。

## 2. 实际依赖方向

实际主要路径：

```text
CompressibleAlgorithm -> Equation::Compressible::System
                     -> discretization dispatch -> methods/numerics 实现
MultiPatchAlgorithm  -> 同一个 equation System
                     + 自己实现 Euler/SSPRK3/RK4 更新公式
PressureBased::Corrector -> 自己组装压力离散/线性系统
Eulerian PressureStepper -> PhaseEquationAssembler -> discretization/linearAlgebra
models -> boundary/discretization 几何与微分函数
methods/numerics -> boundary/discretization（形成反向边）
infrastructure/io CaseConfig -> algorithm/time RunControl
```

### 2.1 重点禁止方向的完整直接 include 清单

每一行是一个 include 位置；目标给出实际解析路径。表中的规则、处理及数值影响适用于该行，而非仅列关键词命中。

- D1：A §3.3、§4.2、§11；CFD 算子不能由 methods 反向调用 solver。把已有 CFD 实现归回 discretization/boundary；真正 Field-free 的数学函数留在 methods。保持公式及调用顺序时不应改变结果。
- D2：A §3.4、§4.3、§11；models 调用具体边界/微分实现，closure 与离散耦合。将边界闭合和梯度求值放到现有 boundary/discretization，由 equation 提供模型所需值；不新增适配层。改变梯度、边界执行时机可能改变结果。
- D3：A §4.4、§11、§17；IO 值对象引入 algorithm 头。将现有纯 RunControl 值定义置于中立配置位置，Driver 留在 algorithm；不改变结果。
- D4：A §4.4、§11 的字面依赖规则冲突；backend 实现使用线性代数公共接口/块消元代码。A §3.5 又允许 HYPRE backend，因此不能误判成 backend 控制 equation。先明确后端实现 target 与公共代数契约归属，保留现有数学算法；单纯整理依赖不改变结果。

| include 位置 | 解析目标 / 当前依赖 | 规则与推荐处理 | 风险 | 可能改变数值 |
|---|---|---|---|---|

| `src/infrastructure/io/SF_caseConfig.h:9` | `src/solver/algorithm/time/SF_time.h` | D3 | P2 | 仅整理依赖不应改变 |

| `src/infrastructure/mpi/backend/SF_hypreBackend.cpp:4` | `src/solver/linearAlgebra/hypre/SF_hypre.h` | D4 | P2 | 仅整理依赖不应改变 |

| `src/infrastructure/mpi/backend/SF_schurPreconditioner.h:5` | `src/solver/linearAlgebra/SF_blockSchur.h` | D4 | P2 | 仅整理依赖不应改变 |

| `src/methods/numerics/scalar/SF_scalarTransport.h:10` | `src/solver/boundary/SF_boundaryGeometry.h` | D1 | P2 | 仅整理依赖不应改变 |

| `src/methods/numerics/scalar/SF_scalarTransport.h:11` | `src/solver/boundary/reconstruction/ILW/SF_boundaryClosure.h` | D1 | P2 | 仅整理依赖不应改变 |

| `src/methods/numerics/structured/SF_iteration.h:6` | `src/solver/discretization/structured/SF_dimension.h` | D1 | P2 | 仅整理依赖不应改变 |

| `src/methods/numerics/viscous/SF_viscous.h:18` | `src/solver/discretization/structured/SF_vectorCalculus.h` | D1 | P2 | 仅整理依赖不应改变 |

| `src/models/ibm/method/ghost/SF_ilwClosure.cpp:11` | `src/solver/discretization/structured/SF_dimension.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/ibm/method/ghost/SF_weightBuilder.h:10` | `src/solver/discretization/structured/SF_dimension.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/heat/SF_wallHeatSource.h:10` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/interfaceModel/levelSet/SF_hjWeno.cpp:10` | `src/solver/discretization/structured/SF_dimension.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/phaseChange/RPI/SF_wallMapping.cpp:10` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/phaseSystem/SF_phaseBoundary.cpp:6` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/phaseSystem/interphase/SF_lift.cpp:6` | `src/solver/discretization/structured/SF_vectorCalculus.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/phaseSystem/interphase/SF_turbulentDispersion.cpp:6` | `src/solver/discretization/structured/SF_vectorCalculus.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/physics/phaseSystem/interphase/SF_wallLubrication.cpp:6` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/RAS/SF_kOmegaSST.cpp:12` | `src/solver/discretization/structured/SF_vectorCalculus.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_equationModelOps.cpp:6` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_equationModelOps.cpp:7` | `src/solver/discretization/structured/SF_dimension.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_equationModelOps.cpp:8` | `src/solver/discretization/structured/SF_vectorCalculus.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_equationSystem.cpp:6` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_kOmegaSSTEquation.cpp:8` | `src/solver/discretization/structured/SF_vectorCalculus.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_turbulence.cpp:13` | `src/solver/boundary/SF_boundaryGeometry.h` | D2 | P2 | 可能，若改变边界/微分时机 |

| `src/models/turbulence/SF_turbulence.cpp:14` | `src/solver/discretization/structured/SF_dimension.h` | D2 | P2 | 可能，若改变边界/微分时机 |

未发现已解析的 `core -> solver/models/IBM concrete/infrastructure`、`linearAlgebra -> equation`、`equation -> algorithm`、`discretization -> algorithm/equation` 或 `boundary -> algorithm/equation` 直接 include。core 中 `Physics::EquationSet::Model` 是 `src/core/state/SF_equationSet.h` 定义的核心契约，不是 models 实现依赖。

全量词法扫描的额外候选已排除：`methods/math/SF_math.h:39` 和 `models/ibm/method/ghost/SF_ilwClosure.cpp:14` 的 `SF_polynomial.h` 当前解析到 `methods/math/discrete`，不是 boundary 中的同名扩展点；core 内裸 `SF_field.h` 当前解析到 core，不能算 core -> IO。

MPI 头仅出现在 `src/infrastructure/mpi/backend/SF_mpiBackend.cpp:17`、`SF_hypreBackend.cpp:11`、`SF_schurPreconditioner.h:9`，未发现 solver 五层直接包含 `mpi.h`。

### 2.2 include 之外的分层问题与旁路

| ID / 文件路径 | 当前依赖或调用关系 | 为什么不符合 A | 风险 | 推荐处理 | 可能改变数值 |
|---|---|---|---|---|---|
| L1 `src/methods/numerics/CMakeLists.txt`；`convection/SF_WENO3.cpp`、`SF_WENO5.cpp`、`SF_WENO7.cpp`、`SF_TENO5.cpp`；`time/SF_time.h`；`scalar/SF_scalarTransport.h`（均在 methods/numerics） | discretization dispatch 向 methods 下放 WENO、通量、标量 PDE、时间积分；这些实现读取 Field/residual | §3.3、§4.2：这不是 solver-independent 数学库；D1 是此归属问题造成的循环 | P1 | 后续按现有算法族归回 discretization；保持纯数学函数；不为移动新增 public class | 仅移动不应；整合 RK/边界可能 |
| L2 `src/solver/algorithm/SF_multiPatchTime.cpp:193`；`src/solver/algorithm/pressureBased/SF_corrector.cpp:32`；`src/solver/algorithm/pressureBased/SF_kkt.cpp` | algorithm 实现 RK 状态公式、邻点度规/间距/压力行组装及 KKT 耦合组装 | §3.1–3.3、§5.3：algorithm 应控制时序，equation 定义耦合，discretization 负责公式 | P1 | 将既有数值函数归到现有 equation/discretization，algorithm 保留调用顺序；内部实现单元即可 | 是，矩阵符号、行索引和 RK 顺序敏感 |
| L3 `src/solver/linearAlgebra/SF_blockSchur.h:15`；`src/solver/linearAlgebra/SF_linearAlgebra.h` | blockRoles 数值 1/2/3/4 驱动压力、速度、乘子、刚体块的不同近似；pressureSystem 暴露物理分块含义 | §3.5：虽然没有 include equation，却仍知道物理块语义 | P2 | equation 明确决定块分组和近似方式，既有代数类只消费块索引/矩阵；先确认是否能简化现有接口 | 是，若改变预条件近似、容差或消元顺序 |
| L4 `src/solver/discretization/source/SF_sourceTerm.h:20–22` | 直接包含 `models/physics/gravity/SF_gravity.h`、`heat/SF_wallHeatSource.h`、`mrf/SF_mrf.h`，Sp 分派具体物理源 | §4.3、§11：equation 应决定物理 source 如何进入方程；discretization 不应成为物理模型组合根 | P2 | equation 评估已有模型源，离散层消费贡献；不用另加 SourceManager | 是，若更改源项评估时间或符号 |
| L5 `src/methods/numerics/Flux/SF_face.h:14`；`src/solver/boundary/reconstruction/ILW/SF_ILW.h:12` | 包含 `models/ibm/SF_cellType.h` | §4.2、§5、§11：通用方法依赖 IBM；边界依赖具体 IBM 辅助实现而非中立分类数据 | P2 | 读取 core 已有分类值/最小谓词；具体 IBM 分类生成留 subsystem | 保持分类判定不应 |
| L6 `src/solver/algorithm/pressureBased/SF_kkt.cpp:9` | 包含 `models/ibm/constraint/variational/SF_constraintBlock.h` | §5.3、§11：algorithm 穿透 IBM 接口实施 constraint assembly | P1 | 让现有 equation/coupled assembly 消费约束贡献，algorithm 只发起求解 | 是，涉及 J/Jᵀ 和行布局 |
| L7 `src/solver/equation/eulerian/SF_equations.h:10`、`SF_workspace.h:6`；`src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.h:7,9`、`.cpp:6–9` | equation 包含具体 turbulence EquationSystem/PhaseSystem；algorithm 直接组合 PhaseBoundary、PhaseSource、gravity、MRF、heat、phaseChange | §4.3、§11、§17：部分物理组合仍落在 algorithm/equation 具体实现，未完全由 application 装配 | P2 | 优先复用现有 equation 与模型接口，把配置组合上移 application；不为每种模型新建接口 | 保持既有调用顺序不应 |
| L8 `src/models/physics/multiphase/SF_multiphase.cpp:336,426`；`src/app/application/run/SF_equationCoupling.cpp:161` | beginStep -> model.advance -> advanceTemperature -> Scalar::advanceForwardEuler | §4.3、§9：物理模型实际推进温度 PDE；主流场 RK 与 thermal Euler 分离 | P1 | 模型提供系数/源；温度方程及其显式时间更新进入既有 equation/积分器 | 是；现有分裂方式和时间精度会受影响 |
| L9 `src/models/turbulence/SF_equationSystem.cpp:67,146,279` | model 模块自己实施标量边界、初始化/prepare/commit；微分闭式依赖具体 stencil | §3.4、§4.3：closure 与边界/时间层存储职责混合；并非发现 model 调用 NS 主 stepper | P2 | closure 数值保留 models，标量边界/方程 workspace 归既有 solver 层 | 是，边界与 previousConserved 时机敏感 |
| L10 `src/solver/equation/compressible/SF_compressible.cpp:47,63,74`；`src/solver/system/SF_systemBuilder.cpp`、`SF_resolvedSimulationSystem.h` | definition.contains(TermKind) 只开关整族算子；具体配置决定 div/laplacian/Sp；未发现实际 AssemblyPlan 驱动链 | §10：尚未实现逐 term binding/assembly，描述信息不能当成执行权威 | P2 | 冻结纯描述扩张，先用已有 equation 主线落实 term 到实际算子绑定；不新增 metadata 层 | 当前只澄清契约不应；真正接管装配时可能 |

## 3. 同名头文件完整清单

A §12 与 AGENTS §7 要求明确 include 路径、避免重复公共 basename。下表全部为 P2；每行推荐先把调用点写为 `src` 相对明确路径，再评估是否确需统一名称/删除未使用入口。本次不重命名。仅路径消歧且保持解析目标不改变数值；错误切换到另一个实现可能影响编译或数值。

当前 1325 个已核对位置中，没有发现同一个 include 因编译单元不同而解析到不同目标。**同名不等于已经编错**。局部同目录 include 有确定优先级；不同目录均 PUBLIC 导出时，异目录调用仍依赖搜索顺序。

| basename | 所有文件路径 | 歧义判断（均 P2） |
|---|---|---|

| `SF_IO.h` | `src/app/gui/IO/SF_IO.h`<br>`src/infrastructure/io/SF_IO.h` | infrastructure/mesh/SF_mesh.h:19 使用裸名；GUI IO 与 infrastructure IO 均导出目录，当前目标为 infrastructure IO。 |

| `SF_application.h` | `src/app/application/SF_application.h`<br>`src/app/gui/Application/SF_application.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_compressible.h` | `src/solver/algorithm/SF_compressible.h`<br>`src/solver/equation/compressible/SF_compressible.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_empty.h` | `src/solver/boundary/reconstruction/ILW/SF_empty.h`<br>`src/solver/boundary/reconstruction/linear/SF_empty.h` | linear 与 ILW 同时 PUBLIC 导出；局部确定，外部裸名风险。 |

| `SF_equationCoupling.h` | `src/app/application/run/SF_equationCoupling.h`<br>`src/core/interfaces/SF_equationCoupling.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_eulerian.h` | `src/solver/algorithm/pressureBased/eulerian/SF_eulerian.h`<br>`src/solver/discretization/eulerian/SF_eulerian.h`<br>`src/solver/equation/eulerian/SF_eulerian.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_executionRuntime.h` | `src/core/interfaces/SF_executionRuntime.h`<br>`src/infrastructure/execution/SF_executionRuntime.h` | linearAlgebra 两个布局头仍用裸名；core/interfaces 和 infrastructure/execution 分别导出。当前解析 core，潜在搜索顺序风险。 |

| `SF_field.h` | `src/core/field/SF_field.h`<br>`src/infrastructure/io/field/SF_field.h` | 大量跨目录裸名当前解析 core；IO/field 子目录未由 SF_io PUBLIC 导出，同目录 IO 实现解析自身。当前未证实冲突，改 include roots 时有风险。 |

| `SF_fixedValue.h` | `src/solver/boundary/reconstruction/ILW/SF_fixedValue.h`<br>`src/solver/boundary/reconstruction/linear/SF_fixedValue.h` | linear 与 ILW 同时 PUBLIC 导出；各实现同目录引用确定，外部裸名会依赖顺序。 |

| `SF_model.h` | `src/app/application/model/SF_model.h`<br>`src/core/model/SF_model.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_output.h` | `src/app/application/model/SF_output.h`<br>`src/app/application/run/SF_output.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_polynomial.h` | `src/methods/math/discrete/SF_polynomial.h`<br>`src/solver/boundary/reconstruction/polynomial/SF_polynomial.h` | math、IBM、ILW 有裸名；当前指向 math/discrete。boundary/polynomial 未单独导出且只有扩展点声明，当前不会抢占。 |

| `SF_state.h` | `src/core/state/SF_state.h`<br>`src/models/physics/interfaceModel/levelSet/SF_state.h` | 当前使用同目录或限定路径，未发现实际误解析；同名公共入口增加异目录调用与未来 include roots 调整风险。 |

| `SF_structured.h` | `src/methods/numerics/structured/SF_structured.h`<br>`src/solver/discretization/structured/SF_structured.h` | 两个聚合头职责交错；目前显式路径消歧，潜在公共 basename 冲突。 |

| `SF_symmetry.h` | `src/solver/boundary/reconstruction/ILW/SF_symmetry.h`<br>`src/solver/boundary/reconstruction/linear/SF_symmetry.h` | linear 与 ILW 同时 PUBLIC 导出；局部确定，外部裸名风险。 |

| `SF_time.h` | `src/methods/numerics/time/SF_time.h`<br>`src/solver/algorithm/time/SF_time.h`<br>`src/solver/discretization/time/SF_time.h` | 三个不同职责；跨层调用检查为显式路径，未发现裸名误选；methods time 与 algorithm time 导出同名风险仍在。 |

| `SF_zeroGradient.h` | `src/solver/boundary/reconstruction/ILW/SF_zeroGradient.h`<br>`src/solver/boundary/reconstruction/linear/SF_zeroGradient.h` | linear 与 ILW 同时 PUBLIC 导出；局部确定，外部裸名风险。 |

## 4. authoritative state 与 Field 职责

| ID / 文件路径 | 当前状态/调用关系 | 为什么不符合 A / 证据边界 | 风险 | 推荐处理 | 可能改变数值 |
|---|---|---|---|---|---|
| S1 `src/core/interfaces/SF_interfaces.h:164`；`src/core/state/SF_stateBundle.h:20`；`src/solver/algorithm/SF_compressible.cpp:274`；`src/solver/algorithm/SF_multiPatch.cpp:267` | SolverState 同时有 bundle、field、fields；bundle 又有 canonical、patches。单块优先 canonical，多块优先 patches 再回退；注册通信字段是另一组指向内存的 views | §7、§9：可独立改写的主状态入口不止一个，且没有统一相等性校验；不是证明数据数组复制了五份 | P1 | 以已有 patches 作为唯一状态枚举，单块也是一个 patch；注册表只引用它；迁移调用后删除回退入口 | 是，若现有路径实际上选择了不同 Field |
| S2 同上；`src/solver/algorithm/SF_compressible.h:195`；`src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp:304` | stepper 的 physicalTime/stepCounter/dt、SolverState 和 StateBundle 分别存时间；仅首步从 state.time 初始化，以后反向回写 | §7–9：外部时间/重启状态与内部计数均可写，bundle 时间不是实际推进的唯一权威 | P1 | 一个状态对象持时间；StepResult 仅结果快照；明确 restart 初始化入口 | 是，影响 stage time、移动壁面与输出时刻 |
| S3 `src/app/application/run/SF_singleFluid.cpp:359–363`；`src/solver/algorithm/SF_compressible.cpp:306`；`src/core/state/SF_variableRegistry.h` | multiPhaseVariables 复制为 bundle.transported；coupling 保留原 registry，stepper 临时改用 bundle 副本；distributed 另由原 registry 注册 | §7、§8：变量描述/列表存在多个可变容器，后续 add/remove 可令积分和 halo 集合不同；数组本身仍通过指针共享 | P1 | 唯一 registry，积分与分布式注册引用同一描述来源；禁止运行中各自修改注册集合 | 是，可能改变被推进/同步变量集合 |
| S4 `src/core/field/SF_field.h:44,99`；`src/core/state/SF_stateBundle.h:23`；`src/solver/algorithm/SF_compressible.h:145` | Field.equationSet、bundle.equations、stepper.equationSet_ 三个可独立绑定指针；EOS CFL/通量用 stepper 指针，Field.thermodynamicState 用自身指针 | §7：同一闭合/布局选择存在多个入口，没有强制同一模型；共享 const 指针可避免内容写入，但无法保证选的是同一个对象 | P1 | 保留一种绑定权威，调用阶段从该处读取；初始化校验布局与闭合一致 | 是，EOS/声速/能量闭合敏感 |
| S5 `src/models/physics/multiphase/SF_multiphase.h:60–69`；`.cpp:516–542,660` | alpha 与 phaseMass 都有可写 accessor；初始化 alpha -> phaseMass，后续 mass -> alpha；温度可推进也可由 conserved 派生 | §7：衍生量与主量的写权限未分开；当前转换路径已有方向，不能直接判定两个数组均为每步积分主量 | P1 | 明确 mixture 的 phaseMass 主量与 alpha cache；初始化转换一次，运行期 cache 只读；区分 thermal 主量与 mixture 派生温度 | 是，禁止直接删温度/alpha 或改变热力学重构 |
| S6 `src/core/field/SF_field.h:50–63,85,115`；`src/core/field/SF_fieldComponents.h` | Field 同时拥有守恒存储/布局、EOS 缓存、网格度规、BC sets、GlobalDof ownership、IBM signedDistance/wall/image/normal/velocity、flux/residual；外部可直接写底层数组，缓存失效靠手工通知 | §4.1、§7.1：Field 仍是跨几何、物理、workspace、IBM 的 God Object；已有 storage struct 只是内部分类，尚未形成职责边界 | P1 | 优先利用已有 GeometryStorage、BoundaryMetadata、IBMGeometryStorage 分清 owner 与访问权限；IBM 专属数据由 subsystem 管理，flux/residual 作为 workspace；不为每个 accessor 加新接口 | 是，数组布局、生命周期、缓存失效和 owner 写权限必须保留 |
| S7 `src/models/ibm/SF_IBM.h:64–76`；`src/models/ibm/SF_IBM.cpp:21–38`；`src/models/ibm/method/SF_method.h:114–134` | IB 保存 config/selection/capabilities/descriptor/fluidPorts；forcing 又存 config/descriptor；selection.fluidPorts 被复制；Field 另持 ghost 几何 | §7、§8：方法描述重复派生并缓存；当前 setup 一次生成、查询只读，未证明实际漂移，因此低于流场重复入口风险 | P2 | 保留一个冻结方法选择；其余由其推导或引用。分别记录几何来源、patch 映射与 forcing workspace，不合并不同物理量 | 配置去重不应；几何/乘子迁移可能 |

以下内容经检查不应列为“重复 authoritative 流场”：

- Runtime 和 ParallelCoordinator 的 `StateBundle*` 是借用指针；`src/infrastructure/execution/SF_executionRuntime.cpp:13` 把同一 bundle 传下去，freshness 由该 bundle.distributed 维护，未发现 Runtime 自己复制流场数组。
- `src/models/ibm/method/SF_method.h` 的 multiplier 与 laggedMultiplier 是不同时间层；surfaceSystem 是约束空间状态，不能和 Eulerian 主守恒场混为同一概念。
- `src/models/ibm/method/SF_methodState.cpp` 对表面乘子执行 canonical COPY，对贡献执行 SUM；未发现可据此认定“平均主状态”的证据。
- `src/solver/algorithm/SF_multiPatchTime.cpp` 的 q0/k1…k4、`src/solver/equation/eulerian/SF_workspace.h` 的 previous time level 是合法算法 workspace。
- Eulerian PressureStepper 以 PhaseSystem 的相变量为主状态，canonical Field 主要用于 geometry；`SF_pressureStepper.cpp:299` 有指针一致性检查。它说明 SolverState.field 命名契约过载，但不是两个流场都被推进的直接证据。

## 5. legacy / compatibility API

| ID / 文件路径 | 当前接口/调用链 | 违反规则或保留理由 | 风险 | 推荐处理 | 可能改变数值 |
|---|---|---|---|---|---|
| C1 `src/solver/algorithm/SF_compressible.cpp:282–364`；`src/solver/algorithm/SF_multiPatch.cpp:281–336` | setter 保存默认服务；advance 将 Services 非空项临时覆盖，再在成功/异常路径恢复 old 指针；单块还替换 registry 与 EOS | §8 明确禁止 temporary override + restore；null 不能表达本次禁用已有服务；下层 runtime attach 状态不会随指针恢复而回滚 | P1 | 每个 stepper 生命周期一次稳定注入，advance 只接受状态；迁移调用后删除 setter/旧服务回退 | 是，若已有调用依赖混合默认值 |
| C2 `src/solver/algorithm/SF_compressible.h:70–88`；`src/solver/algorithm/SF_multiPatch.h`；`src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.h:24–28` | public step 与 advance 同存；advance 包装 step 并承担状态绑定、dt 上限与回写 | §8–9：两个公共入口具有不同前置条件；step 可绕过 bundle/runtime 契约 | P1 | 保留一个外部推进入口，内部 step 改为实现细节；先查完 src/test 调用点 | 仅收口不应；遗漏旧初始化会改变 |
| C3 `src/solver/algorithm/SF_compressible.h:159–176` | setMultiPhaseAuxFields 建 legacyTransportedVariables，和 setTransportedVariables 同存 | §7–8：明确的单辅助场 compatibility API；src 未找到调用该旧 setter 的位置 | P2 | 确认测试/外部使用后删除旧入口和 legacy registry，不继续转发 | 保持注册描述一致不应 |
| C4 `src/methods/numerics/time/SF_time.h:163,348`；`src/solver/discretization/time/SF_time.h:17` | auxField/auxRHS 与 VariableRegistry 两套 Euler/RK4/ddtDispatch overload；标量组合也分别实现；stage 名字符串驱动更新 | §8–9：旧接口保留了额外算法路径；旧标量路径与 registry 的校验不完全一致 | P1 | 统一已有 registry 积分路径，删 aux overload；用明确 stage 调度避免字符串协议，勿新增 public 策略层 | 是，校验与标量更新顺序不同 |
| C5 `src/solver/algorithm/SF_compressible.cpp:46`；`src/solver/algorithm/SF_multiPatch.cpp:85`；`src/solver/boundary/SF_pipeline.cpp` | boundaryPipeline 优先，否则各 stepper 手写 physical BC -> halo -> IBM -> halo；boundaryApplicator 与 override 同存 | §8–9：统一管线已有实现，旧管线仍独立维护 | P1 | 所有运行路径稳定注入已有 Pipeline，删除手写回退；保留明确 operator contract | 是，边界先后和重复通信敏感 |
| C6 `src/infrastructure/execution/SF_executionRuntime.h:26–31` | Runtime(parallel) 构造注入与 setParallelCoordinator 同存；attachState 已把 bundle 绑定到旧 coordinator 后仍可换服务 | §8：稳定依赖可在运行中改变，setter 未重新 attach 现有 state；src 无 setter 调用 | P2 | 保留构造注入，移除未使用 setter；不是所有 setter 都应机械删除 | 正常配置收口不应 |
| C7 `src/core/interfaces/SF_interfaces.h:193`；`src/solver/algorithm/SF_compressible.cpp:139`；`src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp:327` | dt reducer 与 Runtime.globalMinimum 并存；单块优先 runtime，Eulerian 优先 reducer | §8–9：同一归约服务有两个入口且优先级不同 | P1 | 使用统一 runtime 归约入口，测试替身实现同一契约 | 是，全局 dt 可能改变 |
| C8 `src/app/application/run/SF_services.cpp:14–69`；`src/models/ibm/SF_IBM.cpp:72–107` | IBConstraintAdapter -> IB -> forcing，多个方法只是转发；IBBoundaryAdapter 同样转发 ghost | §13：调用层次可简化；但不同 interface 端口的隔离是真实职责，不能只因 Adapter 名字认定全部过时 | P2 | 检查能否直接复用现有端口，删除纯重复转发；保留 ghost/constraint 语义区别，不再新增 Adapter | 仅转发消除不应 |
| C9 `src/app/application/model/compatibility/` 全目录；`src/app/application/model/SF_model.cpp`；`src/app/application/run/SF_equationCoupling.cpp:153–301`；`src/models/physics/multiphase/SF_multiphase.h` | 旧 case decoder 与 native ModelDescription 并存；另有 LegacyMultiphaseEquationCoupling 和 MultiPatch 对应实现继续驱动模型 | §8、§17：外部格式 decoder 位于 application 是合理边界，不能等同于内部兼容 lifecycle；真正需收口的是运行期双接口和 L8 的模型推进 | P2（格式）/P1（运行期） | 保留明确外部格式支持至弃用决定；统一内部状态/方程路径，避免 decoder 产生第二条 solver 流程 | 格式等价转换不应；运行期合并可能 |

未发现 CompressibleAlgorithm 构造函数注入某个 service、同时又 setter 注入同一 service 的证据：其构造参数是 config。它的问题是 setter + per-call Services；确实存在 constructor + setter 同一依赖的是 Runtime。ResultWriter 构造 config、setter 注入 parallel 是不同依赖，不能误报。

## 6. single-patch / multi-patch timestep 调用链

以下比较同一可压缩显式主线。Eulerian SIMPLE/PISO/PIMPLE 的迭代在数学上不同，不因它有自己的 outer/pressure corrector 就判为重复 RK。

| 阶段 | single：`src/solver/algorithm/SF_compressible.cpp` | multi：`src/solver/algorithm/SF_multiPatch.cpp` 与 `SF_multiPatchTime.cpp` |
|---|---|---|
| 入口/选状态 | advance:271 -> canonical/field -> step:130 | advance:264 -> patches/fields/canonical/field -> step:210 |
| 初始边界 | prepareBoundaryState:46（dt=0） | prepareBoundaryState:85（dt=0） |
| CFL | EOS RusanovEOS::deltaT 或理想气体 deltaT；runtime/reducer | 遍历 patch 理想气体 deltaT -> runtime globalMinimum |
| 模型 begin | equationSystem.beginStep；另有 transportModel BC/correct/BC | equationSystem.beginStep(fields)，模型工作委托 coupling |
| RK stage | ddtDispatch -> methods/numerics/time；RHS lambda:182 | 自己 stepEuler:251 / stepSSPRK3:279 / stepRK4:321 |
| RHS | 边界 -> prepareRHS -> convection -> canonical finalize -> diffusion/sources -> coupling RHS -> SUM residual | assembleRHSAll:133；先所有 patch convection，再统一 canonical finalize，再各 patch diffusion/source，最后统一 residual SUM |
| stage 提交 | postStage 发布 conservative/transported writeOwned -> validate | publishIntegratedState 去重变量名 -> validateStateClosureAll |
| 校正 | flowAlgorithm.correct；forcing 校正后额外 WriteOwned | flowAlgorithm.correct；没有同位置的 forcing-specific WriteOwned |
| 步尾 | commitStep -> 非 forcing 路径边界；forcing 路径推迟到下一步边界 -> time/step/observer | 无条件边界 -> commitStep -> time/step/observer |

**T1 — P1，重复 lifecycle。** 文件：上述三个实现文件及 `src/methods/numerics/time/SF_time.h`、`SF_RK4.cpp`、`SF_Euler.cpp`。重复内容包括 CFL 裁剪、beginStep、边界/halo/RHS、Euler/SSPRK3/RK4、stage 发布、状态验证、flow correction、commit、时间回写及 observer。违反 A §3.1、§9；推荐保留一条 patches 生命周期，单 patch 作为长度一容器，由 Runtime 处理执行差异。不能简单对每个 patch 单独跑完整 RK：跨 patch 的 stage 和 canonical flux 必须保持全局阶段屏障。后续整合**可能改变数值**，需验证 stage 节点 `{0,1,1/2}` / `{0,1/2,1/2,1}`、源项时间、COPY/SUM/face owner。

**T2 — P1，已存在顺序差异。** 文件：`SF_compressible.cpp:242–263` 与 `SF_multiPatch.cpp:238–258`。single 在 commit 后边界闭合，multi 在 commit 前闭合；single 对 forcing 校正有显式发布并跳过步尾边界，multi 没有对等处理。违反 A §6、§9 的共用 lifecycle/freshness 契约。现有 application 多块 IBM 绑定 ghost 端口，不能据此宣称当前多块 forcing 已产生错误；但 MultiPatchAlgorithm 暴露 constraint setter，接口能力与通用流程保证不一致。推荐先统一校正输出/commit 的数学契约，再合并流程。**可能改变数值及 MPI 通信序列**。

**T3 — P1，EOS 与变量服务差异。** 文件：`SF_compressible.cpp:137,98,195` 与 `SF_multiPatch.cpp:149,163,189,219`。single 支持 equationSet-aware CFL/通量/验证，multi 的 convection context 传空 thermodynamics 并使用 PhysicalState 验证。违反 A §7–9 中执行模式不应引入不同 state/lifecycle 的目标；也可能是 workflow 明确限制的能力边界，不能推断所有多块 EOS case 都已获支持。推荐对合法工作流先显式核实能力矩阵，再统一每 patch 的 EOS/transport 来源。**扩展或改变闭合路径会改变数值**。

Eulerian 路径定位：`src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp:296,347`：advance -> stableTimeStep -> step -> boundary/sync -> outer(interphase/source/turbulence -> continuity -> momentum -> pressure/nonOrthogonal corrections -> energy/turbulence) -> workspace/turbulence commit -> diagnostics。其 dt/clock/services 问题归 S2/C7；其相方程主线不应强行替换成显式 RK。

## 7. 大文件及实际阶段

行数按文本换行统计，阈值取 >=800；长度本身不是数值错误。下表均对应 A §14 的实现组织要求；推荐保留现有 public class/函数入口，按列出的阶段划分 internal implementation units。不建议新增 Manager/Adapter/Interface。纯拆文件不应改变数值，表中“敏感”说明搬动时必须保持的不变量。

| 文件路径 | 行数 | 实际包含阶段 / 当前调用关系 | 风险 | 推荐拆分边界与数值影响 |
|---|---:|---|---|---|
| `src/infrastructure/mesh/communication/SF_communicationPlan.cpp` | 2241 | buildHaloExchangePlan:2048 汇总：dual volume、donor 空间搜索、trilinear 反演、Lagrange tensor 权重、conformal ghost map、GlobalPoint/GlobalDof owner、canonical face plan | P1 | 搜索/插值、接口拓扑、DOF/face ownership、校验分 internal 单元；保持 owner 选择与索引，误改可能影响数值 |
| `src/models/ibm/method/ghost/SF_ilwClosure.cpp` | 2222 | 几何/曲率、采样扩展、2D/3D 模板、WENO fit/光滑度、预计算 plan、特征/壁面一阶导数、高阶导数、Taylor ghost 写出及诊断 | P1 | 以预处理 plan 与运行期 closure 为边界，纯数学复用 methods；ILW 阶数/权重/失败策略不可随拆分改变 |
| `src/infrastructure/mpi/SF_haloExchange.cpp` | 1510 | 邻居 payload、主量/标量 pack、owner COPY、tensor 插值、单 Field 与 block 数组交换、canonical face flux 装配 | P1 | packing、state COPY、canonical flux 分单元；保留 collective 次序与 COPY/SUM 区别 |
| `src/infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.cpp` | 1354 | SFM 解析/集合规范化、拓扑校验、source canonical metrics、接口报告、loadFiles 分解/组装、路径处理 | P2 | 解析、metrics、组合构建、校验分单元；几何/分区顺序变更可能影响数值 |
| `src/infrastructure/mesh/createMesh/SF_meshGen.cpp` | 954 | structured block 点生成、global point 去重、边界/接口报告、combined mesh 输出、曲边表与 TFI | P2 | 生成、拓扑编号、输出分单元；保持坐标/ID 不应改变数值 |
| `src/infrastructure/mesh/decompose/SF_meshDecompose.cpp` | 832 | extrusion 判定、logical cuts、source slicing、覆盖校验、composite partition、跨 source ghost 坐标 COPY | P2 | 切分规划、数据切片、覆盖验证分单元；owner/坐标映射变更可能影响数值 |
| `src/models/physics/multiphase/SF_multiphase.cpp` | 863 | 配置/初始化、thermal Euler、mixture properties、mass/alpha/T 闭合、守恒量初始化、RHS、commit/BC、源项/诊断 | P1 | 首先按 L8/S5 纠正职责，再拆初始化/闭合/贡献/诊断；时间分裂改变会改变数值 |
| `src/app/gui/VTKView/SF_widget.cpp` | 2053 | SFM/VTK/PVD 读取、geometry/array 复用、网格预览、时序监视、播放、切片、标量选择、相机/高亮 | P2 | 读取、场景/切片、播放/监视分内部实现；不改变 solver 数值，可能影响显示 |
| `src/app/gui/Controller/SF_controller.cpp` | 1529 | 配置解析/后端检测、UI 连接、项目编辑、mesh/solver 进程、输出 task、并行命令、preview case 拷贝/修改 | P2 | 项目操作、进程/输出、preview 分内部实现；纯拆分不改数值，误改生成配置会影响输入 |
| `src/app/gui/IO/SF_project.cpp` | 1300 | 文件模板、dictionary 编辑、project open/reload/save、模块/field 管理、网格/STL 导入、几何写出 | P2 | 模板/解析、项目读写、导入分单元；配置字节应保持，不应影响数值 |
| `src/app/gui/GUI/SF_mainWindowTree.cpp` | 1221 | tree rebuild:178、条目编辑提交、output task、enum choices、module selector/context menu/highlight | P2 | tree 构建、编辑动作、输出导航分单元；不改变 solver 数值 |
| `src/app/gui/IO/SF_configDocument.cpp` | 910 | 字符/注释解析、配置条目定位与修改、Foam dictionary/blockMesh 解析 | P2 | 解析与文本编辑分内部实现；§14；保持读写语义不应改变数值 |
| `src/app/gui/Geometry/SF_editor.cpp` | 933 | blockMesh 字符解析、顶点/弧边/面投影、绘图交互、模式切换/loadFile | P2 | 解析与绘图分内部实现；纯分离不改数值，坐标编辑行为应保持 |

还需纳入实现体审查的头文件：

| 文件路径 | 行数 | 阶段/职责与 A 冲突 | 风险 | 推荐处理 | 可能改变数值 |
|---|---:|---|---|---|---|
| `src/models/ibm/method/ghost/SF_weightBuilder.h` | 930 | 权重构建、几何/样本选择、ILW 预计算组织；大量具体实现暴露于头中，§14 | P2 | 保留 WeightBuilder 公共入口，将非模板实现放内部 cpp | 仅移动不应；样本与容错规则敏感 |
| `src/methods/numerics/scalar/SF_scalarTransport.h` | 802 | 标量 BC/ILW、bounds、对流扩散 RHS、Euler update；§3.3/4.2，见 D1/L1 | P1 | 归既有离散/边界职责，分 BC/RHS/更新实现 | 是，bounds、边界时机与积分敏感 |
| `src/core/config/SF_config.h` | 987 | 多模块配置、枚举转换/校验、理想气体约束；§16 巨型配置倾向 | P2 | 复用已有 configTypes 与子配置；调用方只接收必需子对象 | 保持默认值/校验不应 |

## 8. 建议处理顺序与验证边界

1. 先处理同名头/显式路径与 D1/D3；这些可以作为保持行为的独立整理。
2. 明确 S1–S4 的唯一状态、时间、registry、EOS 来源及 C1 的稳定服务契约。
3. 对 T1–T3 建立相同输入的单 patch / 多 patch / MPI 对照，再统一显式 lifecycle；保留每个 stage 的 canonical face 计算及 COPY/SUM 屏障。
4. 处理 L8 的模型温度推进和 S6 的 Field 职责；大文件按阶段分内部实现，不增加公共抽象。
5. equation IR 真正接管装配前暂停描述层扩张。

后续数值敏感改动应覆盖时间积分阶数、非定常边界 stage 时间、EOS/transport、重启、同物理问题不同分区的残差/全局守恒。IBM 要分别验证 ghost/ILW、forcing、DLM/KKT；不能用 GhostCell 成功代替约束正确性。本次仅架构静态审计，不据此声称上述数值回归已通过。

## 9. 本次变更纪律与验证

- 修改 solver lifecycle：无。
- 增加 public API / 删除 legacy API：均无。
- 新 state source / compatibility path / 跨模块依赖：均无。
- MPI ownership semantics / 数值结果：均未改变。
- 移动文件、重命名类、修改 dispatch：均无。
- fallback：未新增、删除或修改；未进行全库 A/B/C fallback 专项判定。标量 bounds 等行为只作为后续重构必须保留并单独审查的数值敏感项。
- 执行检查：全 src 文件/include/同名头/行数扫描；CMake 与已有 compile_commands 的 include 搜索路径核对；关键状态、RK、压力及 IBM 服务链人工阅读；完成后对 src 全部文件做 SHA-256 前后比较。
- 未运行构建或 solver tests：本次仅文档，运行测试不能替代静态架构问题的定位。
- Git 状态不可用：`.git` 指向 `/Users/lpf/.sonicSolver-git/sonicSolver.git`，`git status` 报 not a git repository；不据此推断工作区原先是否干净。用 src 文件清单和哈希确认本次未修改源代码。
