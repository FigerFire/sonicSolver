# 三契约方程驱动架构：历史差距审计（2026-09-28）

> 此文是迁移前审计。当前 production authority、已删除概念和剩余 legacy 以 [ARCHITECTURE.md](../src/ARCHITECTURE.md) 和 [迁移报告](equation-execution-migration.md) 为准。

**日期：** 2026-09-28。**范围：** 当前工作树的 `src/`、对应 case 输入和现有架构测试。工作树在本轮开始前已有 Phase 28A–28C 未提交修改；本报告以读取到的**当前文件**为准，不把 HEAD 当作运行中的源码。本轮仅更新文档，未改 C++、输入或数值行为。

## 结论与判定方法

目标是三个正交入口：WHAT = Equation System；HOW = Execution Recipe；HOW NUMERICALLY = 每个 term 的 Numerical Recipe。内置 preset 和等价用户声明必须进入同一 IR、编译与运行链。现在已经有 Raw/Executable IR、SystemContribution、压力约束 transformation、结构化 Plan、term/operation provider 解析；**不能据此认定任意用户方程已可运行**。当前最大差距是 native 自定义方程无法 lower、term compiler/operation provider 只覆盖少数方程模式，以及湍流/相方程仍有专用执行实现。

下文用 **P1** 表示阻挡“三契约 + preset/user 同一路径”的核心差距，**P2** 表示局部硬编码或迁移债。未发现本轮文档任务必须立即修复的 P0 数值缺陷。`Unsupported` 且 fail-fast 的组合是**能力缺口**，不能误写成静默数值错误；已验证的专用 kernel 是**Legacy 实现**，不能为了设计整齐而直接替换。建议均留待用户审核后再实施。

本轮列出 **30 项**有源码位置的差距：WHAT 14 项、HOW 6 项、HOW NUMERICALLY 8 项、跨契约验收/诊断 2 项；其中 P1 21 项、P2 9 项。这是按**独立责任缺口**合并的清单，同一问题的每个调用点没有重复计数。

### WHAT：方程与闭合

| ID | 优先级/性质 | 当前源码与关系 | 与目标的差距；建议审核的方向 | 数值风险 |
|---|---|---|---|---|
| W1 | P1 / 能力缺口 | `src/app/application/model/SF_nativeBinding.cpp:70-98,174-183` 的 `equationRegistry` 只接受 Momentum/Continuity/Energy；自定义方程名和表达式明确抛错，`BuiltinEquationSystem` 只接受 compressible/eulerianEulerian。 | 用户无法把等价于预设的方程从 native input lower 到 `Equation::Definition`。应让预设和用户输入产出同一 contribution，保留显式 unsupported 直到 compiler 可执行。 | 仅文档无变化；未来开放输入可能改变方程与结果，须做 A/B。 |
| W2 | P1 / 能力缺口 | `src/solver/system/SF_equationContribution.cpp:127-141` 的 `applyModification` 对 `replace/disable` 一律报未实现；上述 native binding 对 add/extend/replace/disable 也一律拒绝。 | “预设 + 修改”尚不能表达，`SystemModification` 目前主要是语义占位。需 typed payload、明确目标和冲突检查；不可按同名后写覆盖。 | 未来替换方程会改变数值，必须显式请求并冻结旧 baseline。 |
| W3 | P1 / Legacy 分叉 | `src/app/application/SF_inspection.cpp:46-66,137-153` 仍以多相类型和兼容 `pressureBase` 标签决定 `PhysicsTemplateKind` / preset request；`src/solver/system/SF_systemBuilder.cpp:263-300` 按 `singleFluidPreset`、`composition.declared`、Eulerian template、pressure request 分叉装入不同方程包。 | Case label 仍能决定 WHAT 的构造路线，等价显式声明没有同一个 lowering 入口。应把兼容标签限制在入口翻译，最终统一为声明式贡献。 | 高；迁移时可能误改方程包、EOS 或初值绑定。 |
| W4 | P1 / 能力缺口 | `src/solver/system/SF_systemBuilder.cpp:190-241` 的 semantic composition 强制 Continuity + Momentum，perfectGas 强制 Energy，并按 EOS 名做分支。 | Validator 仍围绕既有单流体模板，而非先验证任意数学系统，再检查所选实现能力。应区分数学合法性和现有 provider 能力；不应按“不常见”拒绝。 | 高；不能把旧拒绝简单改为可运行。 |
| W5 | P1 / 不对称 | `src/solver/system/SF_presets.cpp:100-133`、`SF_singleFluidPreset.cpp:36-72` 生成预设方程及固定 storageKey/offset；`SF_stateRealization.cpp:53-69` 仍识别 `rho`、`rhoU`、`p` 和 `phaseMass*` 名字。 | 虽使用中立 descriptor，storage/realization 与若干预设 ID 紧耦合；用户写数学等价但命名/绑定不同的系统不一定得到相同执行产物。需定义可比较的语义 role、storage 绑定和规范化规则。 | 高；改 layout 或字段映射可能改变结果。 |
| W6 | P1 / 不对称 | `src/solver/system/SF_transformation.cpp:27-48,67-96` 以 `C_INCOMPRESSIBILITY`、`E_MOMENTUM`、`p`、`U` 精确 ID 匹配并读取预设定义；`SF_pressureCoupling.cpp:26-54` 同样识别固定约束/方程 ID。 | 等价的用户方程若 ID 不同，无法触发同一 pressure transformation。需要声明式 role/constraint capability 的匹配，同时保留 required storage/operator validation。 | 高；压力约束和 correction 顺序敏感。 |
| W7 | P1 / 双执行来源 | `src/models/turbulence/SF_turbulenceSystemContribution.cpp:40-97` 声明 `k/epsilon/omega` 方程；`src/models/turbulence/SF_turbulence.cpp:24-47,275-298` 的 Manager 按模型名创建并独立 `correct`；`src/solver/algorithm/SF_singleFluidStepper.cpp:231-270,536` 在固定位置调用 transportModel。 | 同一个湍流概念一边进入方程 IR，一边由模型专用校正接口执行；不能证明用户显式写的等价 `k/epsilon` 方程会走相同数值链。应先完成 generic equation/term/provider 绑定，再逐步迁移经过验证的 kernel。 | 很高；湍流应力、输运时序和边界不可机械合并。 |
| W8 | P1 / 双执行来源 | `src/solver/algorithm/eulerian/SF_eulerianStepper.cpp:133-140,502-559` 初始化 `Turbulence::EquationSystem`，在固定 OpId 上调用 `prepare`、`solveTurbulence`、`commit`；`src/solver/system/SF_transformation.cpp:276-288` 通过 `E_TURB_` 前缀声明 operation。 | Plan 反映了阶段，但执行仍按内置湍流对象/命名约定绑定；用户显式等价方程没有相同 provider 路径。需把数学 term 和 schedule 的能力声明与模型名/前缀解耦。 | 很高；Eulerian phase 数值与 MPI 顺序须冻结。 |
| W9 | P2 / 正交性 | `src/solver/system/SF_presets.cpp:203-207` 的 Eulerian 方程模板同时记录 `preset.coupling.PIMPLE`；`SF_systemBuilder.cpp:284-286` 从模板选择它。 | WHAT preset 仍携带特定 coupling 身份；应把默认 HOW recipe 作为独立默认值显式展开，而非“Eulerian 必然 PIMPLE”。现有有效 case 不能突然失去默认行为。 | 中；改默认会改变 loop 次数与结果。 |
| W10 | P2 / 扩展入口 | `src/app/application/model/SF_model.cpp:36-55` 的外部对象注册回调直接修改 `CaseConfig`，不是返回统一的 Equation/Execution/Numerical contribution。 | 扩展模块入口尚不能保证加入 WHAT 后与预设同编译。可以在 composition root 做中立值 lowering；不要新建 runtime ServiceLocator。 | 中；取决于扩展当前行为。 |
| W11 | P1 / 跨轴耦合 | `src/app/application/SF_inspection.cpp:127-150` 把 turbulence enabled/legacy mixture 与 `viscousEnabled` 合成 `fluidDiffusion`；`src/solver/system/SF_systemBuilder.cpp:253-256` 也从 numerics flag 推导 preset 的 diffusion；`SF_singleFluidPreset.cpp:53-59` 据此决定 Momentum/Energy 方程中是否存在 diffusion term。 | 目前 WHAT 的物理 term presence 会受 HOW NUMERICALLY 的配置位和模型开关影响。长期应由物理/闭合预设声明应有的应力/热通量项，再由 Numerical Recipe 选择其离散；不能通过换 stencil 让方程项消失。 | 很高；未来分离两个开关时必须保护 inviscid/viscous 既有结果。 |
| W12 | P1 / 方程表达不完整 | `src/solver/system/SF_presets.cpp:117-132` 的常密度 Raw Momentum 只有 `ddt(U)+div(momentumFlux)[+diffusion]=0`，同时声明 `div(U)=0`；`src/solver/algorithm/pressure/SF_pressureOperators.cpp:414,594-626` 的实际 predictor 使用压力梯度。 | Raw 方程没有显式表达该 pressure-gradient/constraint coupling，导致 WHAT 的展示与数值执行不完全同源。应定义压力项在 raw/derived equation 中的明确数学归属，再保持现有 pressure operator 的公式与次序；不能只给 explain 补一行文字。 | 很高；压力梯度与 Rhie–Chow/矩阵不可在文档重构中更改。 |
| W13 | P1 / 专用方程执行 | `src/models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.cpp:18-41` 声明 `E_LEVEL_SET`；`src/app/application/execution/SF_singleFluid.cpp:282-300,515-520` 又按 legacy mixture/interface/homogeneous 类型构造 `IEquationSystemCoupling` provider；`src/solver/algorithm/SF_singleFluidStepper.cpp:532-536` 在固定位置调用 `beginStep`。 | Level-set/mixture 方程进入了 WHAT，但执行仍依赖模型专用状态注册与服务入口；用户显式等价方程无法仅凭 IR 获得相同 execution。后续应将现有数值 helper 绑定为通用 equation/term operations，先保护原顺序。 | 很高；interface/phase transport、守恒量与 stage 绑定敏感。 |
| W14 | P2 / 约束 preset 专用 | `src/models/ibm/SF_ibmSystemContribution.cpp:76-151` 按具体 `strategy` 字符串决定 ConstraintProjection/MonolithicKKT 和固定 `ibm.*` OpId，monolithic 情况还插入 `E_MOMENTUM/E_PRESSURE` 名字。 | 当前 IBM 是真实的 constraint contribution，但“用户手写等价约束”尚不能自动 lower 到同一 transformation/operation。保留 Ghost boundary 与 DLM/KKT 数学区别，再让约束 capability 而非 preset 名字决定适用性。 | 很高；J/Jᵀ、乘子 ownership 和 MPI SUM 不可改变。 |

### HOW：时间、耦合和方程调度

| ID | 优先级/性质 | 当前源码与关系 | 与目标的差距；建议审核的方向 | 数值风险 |
|---|---|---|---|---|
| H1 | P1 / 能力缺口 | `src/core/config/types/SF_timeRecipe.h:11-37`、`src/app/application/model/parsing/SF_numericsParser.h:115-125` 只注册 forwardEuler、SSPRK3、classicalRK4；`TimeTopology` 只有 ExplicitStages。 | 用户举例的 RK2 以及隐式/半隐式时间并无 recipe/provider。下一步应先确定 stage/constraint contract，再逐个实现，不把 RK2 名字映射到现有方法。 | 很高；stage 系数、时间与约束处理直接影响结果。 |
| H2 | P1 / 组合能力 | `src/solver/system/SF_providerResolver.cpp:33-137` 对 conservative pressure 和 constant-density pressure 使用 stageCount==1、单一 pressure policy、固定方程/约束集合等能力门槛。 | `RK2 + SIMPLE + kEpsilon` 不是当前可运行的正交组合；这属于诚实的 Unsupported。需为所选阶段显式提供 predictor、约束修正与湍流方程 provider，不能重复 full-dt predictor 伪装外迭代。 | 很高。 |
| H3 | P1 / 调度缺口 | `src/solver/system/SF_solvePlan.cpp:142-190,226-310` 的无 fragment 显式路径把全部 transient equation 收入一个 fused stage OpId；有 fragment 时取首个 active fragment 并要求 policy 已被消费。 | 目前不能从任意用户声明的 E1→E2、block(E1,E2)、subcycle(E3) 自动编译通用日程。应以可组合的方程级 plan contribution 扩展，保留已有 RK/pressure 顺序。 | 高；顺序与物理时间提交敏感。 |
| H4 | P1 / 专用 provider | `src/solver/system/SF_providerCatalog.cpp:55-94` 按 `flow.conservative`、`flow.pressure-operators`、`flow.eulerian-pressure` 三组状态特征解析 OpId；`src/app/application/execution/SF_flowLoop.cpp:100-106` 与 `SF_eulerian.cpp:64-68` 构造两种 stepper。 | 统一 PlanExecutor 已有，但 operation 的生产实现仍由专用 provider/stepper 提供；预设等价用户系统不能自动复用任意 OpId。长期应按 operation/term capability 装配，而非恢复 solver-family lifecycle。 | 高；专用 kernel 有已验证基线，不能直接删除。 |
| H5 | P2 / 模型调度 | `src/solver/algorithm/SF_singleFluidStepper.cpp:231-270` 的 transport correction 在显式 step begin 调用；`src/solver/algorithm/eulerian/SF_eulerianStepper.cpp:502-559` 的湍流准备/求解/commit 虽绑定 OpId，仍由内置方程对象执行。 | k/epsilon 与 Momentum 的相对顺序未完全由可替换的 equation schedule 表达。需将模型必要的 boundary/halo 与 commit 声明为 operation 依赖，再保持原顺序迁移。 | 很高；不应把闭合刷新提前或延后。 |
| H6 | P2 / 契约混层 | `src/core/config/types/SF_termRecipe.h:133-137` 把 `TimeRecipe time` 放在 `NumericalRecipeSet`；`src/solver/system/SF_systemBuilder.cpp:474-477` 再以 resolved time recipe 覆写这一字段。 | 当前唯一实际时间选择仍可追踪，但类型上把 HOW 的时间方法包装进 HOW NUMERICALLY 的空间 recipe 集，未来易形成第二个 stage authority。宜在保持旧输入兼容的前提下明确单一 Execution Recipe 输入与 compiled time output。 | 中；迁移时 stage time/系数需数值对照。 |

### HOW NUMERICALLY：逐项离散与 provider

| ID | 优先级/性质 | 当前源码与关系 | 与目标的差距；建议审核的方向 | 数值风险 |
|---|---|---|---|---|
| N1 | P1 / 配置粒度 | `src/core/config/types/SF_termRecipe.h:133-137` 的 `NumericalRecipeSet` 仅有全局 convection/diffusion；`src/app/application/model/SF_nativeBinding.cpp:51-69` 的 `terms` 只允许这两个键。 | 无法为 `Momentum.div`、`k.div`、`epsilon.diffusion` 分别选 recipe，也不能独立配置 `grad(p)`。应增加方程/term 定位与缺省覆盖优先级，编译成逐项绑定。 | 高；默认优先级和浮点顺序必须冻结。 |
| N2 | P1 / 编译覆盖 | `src/solver/system/SF_numericalCompiler.cpp:16-27,69-71,124-183` 只为已知 Mass/Momentum/Energy 角色和精确 solvedUnknown 编译，其他方程被跳过；term 类型只绑定 divergence/diffusion/source。 | 用户新增方程或 turbulence/phase equation 出现在 IR，并不意味着其 term 获得数值 provider。需以每条 equation 的声明和 provider 能力编译，缺项 fail-fast；不能让“能 explain”冒充“能 execute”。 | 高。 |
| N3 | P1 / 融合范围 | `src/solver/algorithm/SF_conservativeRHS.cpp:44-45,146-167` 通过 `requireUniqueRecipe` 取整组 convection/diffusion，再交给 fused compressible assembly。 | 现有高性能融合 kernel 可保留，但尚不能按多个方程/term 的不同 recipe 与 component write-slice 选择贡献。应先建立编译态 slice/一致性契约，再决定可融合的组合。 | 很高；改 flux、source 次序或 residual 符号会改变结果。 |
| N4 | P2 / 模型清单残留 | `src/core/config/types/SF_termRecipe.h:45,61-63,229-240` 仍有封闭 `SourceKind {Gravity,MRF,WallHeat}` 与源项 recipe switch；`src/app/application/model/parsing/SF_modelConfigParser.h:14-25` 按固定名字解析。 | Phase 28C 已将实际 gravity/MRF/wallHeat term provider 放到模型侧，但用户入口/recipe 类型仍知道具体模型清单。应让 model-owned descriptor 提供 ID/参数编译，兼容解析保留在边界。 | 中；迁移注册顺序和 legacy 输入需验证。 |
| N5 | P1 / provider 适用性 | `src/solver/system/SF_termProviderCatalog.cpp:91-121` 的 primitive convection 只匹配 `Momentum(U)` + Upwind1，diffusion 只匹配 Central2/4；`src/solver/system/SF_numericalCompiler.cpp:74-84` 对压力 Momentum 需要特定 predictor binding。 | “WENO5 + Central2”不能被解释成当前不可压缩 `Momentum(U)` 的通用离散选择。应针对数学 term、状态表示、EOS/特征结构和边界能力匹配，保持不支持时 fail-fast。 | 很高；不能把 high-order 改成 Rusanov/Upwind 假装支持。 |
| N6 | P2 / Eulerian 专用 | `src/core/config/types/SF_pressureConfigTypes.h:25-29` 的 phase convection 只有 Upwind；`src/solver/algorithm/eulerian/SF_eulerianStepper.cpp:31-52,98-124` 自有固定 SourceKind→phase source factory 清单。 | 相方程与源项尚未消费与单流体同一个开放的逐项数值 provider 契约。这是 Legacy，不表示现有 Eulerian 数值错误。 | 高；相级源项/CFL/MPI 顺序不可顺手替换。 |
| N7 | P2 / 时间与空间混层 | `src/core/config/types/SF_termRecipe.h:48-49,84-108` 的空间 `TermRecipe` 固定携带 `TemporalRole::ExplicitResidual`。 | 这对当前显式 term 是真实能力描述，但长期不能把 WENO/Central 的身份定义为“只能属于 explicit solver”；隐式/半隐式装配需要独立的时间级和 operator capability。 | 高；不可仅改 enum 宣称隐式可运行。 |
| N8 | P1 / raw 配置残留 | `src/solver/algorithm/SF_conservativeRHS.cpp:112-170` 在运行期仍接整份 `SolverConfig`，把 `ibmBoundary`、`ilwOrder`、viscosity、Pr、gamma 等原始数值值传给 fused assembly；`SF_singleFluidStepper.cpp:529-557` 仍将该 config 带入每个 stage。 | Phase 28C 已编译模型源项参数，但 conservative numerical path 尚未完全只消费编译后的 term/boundary/thermo binding。应先逐项证明编译态 authority，再迁移参数读取；保留已验证的 kernel。 | 很高；任何值或读取时机改变都可能影响高阶数值。 |

### 跨契约验收与文档

| ID | 优先级/性质 | 当前源码与关系 | 与目标的差距；建议审核的方向 | 数值风险 |
|---|---|---|---|---|
| X1 | P1 / 验收缺失 | `test/test_pisoArchitecture.cpp:42-63,718-746` 有 test-only custom Momentum source provider；现有测试未构造“内置 preset vs 等价用户方程”两条完整 case 并比较 executable、numerical、plan、bindings 与结果。 | 当前证明的是**加一个项无需改中央编译器**，不是**任意等价方程同编译同执行**。先选可运行的常密度 PISO 或 conservative Euler 作为最小等价 A/B，再逐步扩展到 turbulence。比较忽略 provenance，保留项顺序、stage 与数值系数。 | 测试本身无变化；它将暴露高风险差异。 |
| X2 | P2 / 文档和诊断 | 更新前 `README.md:1,50-56` 仍称 phase25B 且把常密度 SIMPLE/PIMPLE 写为未实现；`src/app/cli/` 的 `--help` 首行仍称 “compressible FDM solver”；`sonicSolver recipes` 当前只列压缩流/IBM 等模板。 | 对外叙述仍以老版本或 solver-family recipe 为中心，易把目标能力与现状混淆。README/架构文档已在本轮修正；CLI 文案和 recipe catalog 留待源码修改阶段。 | 无数值风险。 |

## 已符合目标、应当保护的实现

- `src/core/system/SF_equationIR.h` 将 Raw 与 Executable 方程分开；`src/core/system/SF_systemContribution.h` 可用中立值贡献未知量、方程、term、约束和 policy。
- `src/solver/system/SF_equationContribution.cpp` 对已支持的 add/extend 有显式顺序和冲突检查；未实现的 replace/disable **报错而非静默覆盖**。
- `src/solver/system/SF_transformation.cpp` 把 pressure constraint 的 algorithmic equation 与原物理方程分开；`SF_pressureCoupling.cpp` 的 plan fragment 表达已实现的 predictor/corrector 顺序。
- `src/solver/system/SF_termProviderCatalog.cpp` 与 Phase 28C 的模型源项 descriptor 已有编译期选择、唯一匹配和冻结 callback；`src/solver/run/SF_planExecutor.cpp` 执行结构化 Plan，不根据 preset 名称选择整个 lifecycle。
- `src/solver/system/SF_providerResolver.cpp` 对缺数值能力的组合设置 Unsupported；这些 fail-fast 检查在实现新 provider 前应保留。

## 建议冻结的验收规则

1. **同一语义、同一 IR：** 选择一个当前真实可运行的方程包，分别由内置预设和用户显式声明产生；逐项比较 unknown role/storage、方程/term 顺序、约束、closure、transformation、execution recipe、numerical binding、Plan 和 provider ID。只剔除 provenance 与展示文字；不可忽略数值系数、读写范围或同步要求。
2. **同一执行：** 在同一网格、初值、边界与 MPI 配置下比较 `dt`、stage time、关键残差、pressure matrix、final field 和数值诊断。发生差异先定位第一个不同 checkpoint，不调整 tolerance 让结构重构通过。
3. **诚实能力：** 用户可声明目标组合，但缺 storage、term、stage、constraint 或 backend provider 时必须在编译/检查阶段报告精确原因；不能自动换 SIMPLE/PISO、WENO/Upwind、EOS 或时间 recipe。
4. **阶段顺序：** 下一阶段先闭环受支持的通用方程编译和执行，再用 k-ε 作更难的验收；`RK2 + SIMPLE + kEpsilon + WENO5` 是长期组合目标，不是下一次提交必须伪造的 green test。

## 本轮验证与边界

这是文档审计，无 C++/CMake/input 修改，未运行数值回归，也未声称修复上表缺口。审计使用 `rg`、源码与 case/test 阅读，并实际检查了 `sonicSolver --help`、`sonicSolver recipes` 的当前输出。之前 Phase 28C 的 host CTest 25/25 与冻结 A/B 结果记录在 [Phase 28C 报告](phase28c-model-provider-ownership.md)；它们不能证明本报告提出的 preset/user 同一性。
