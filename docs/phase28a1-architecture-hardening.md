# Phase 28A.1 — Architecture Hardening Closure

日期：2026-09-27。状态：**Implemented（本阶段范围）**。本报告只描述已落地源码和实际执行结果；未实现的数学能力标为 **Unsupported**，遗留耦合标为 **Legacy**。

## 1. 范围与冻结基线

本轮只处理时间 recipe authority、operation provider 解析、pressure solve-policy 重复语义、equation/physics/turbulence target 环及通信错误契约。未更改 pressure/IBM/WENO/RK 公式、stage 顺序、boundary/halo 顺序、canonical COPY、GlobalDof SUM 或压力参考行。Phase 28A 的单源块、每 rank 至少一个参与 patch/pressure row、正交结构网格支持边界保持原样。

修改前已冻结 HEAD `537889b05234c34217aa1a8eff7ab7844b69b4b9`、全部可追踪/未追踪源码及输入、dirty patch、可执行文件 SHA-256 `1d758ff0539cb4a2f9c298505b98d63f78f146af837de831203dee800efa1ce4`，位置 `/private/tmp/sonic-phase28a1-baseline-20260927`。修改前 architecture checker 为 2 条既有 allowlist，完整构建及 22/22 CTest 通过。未清理或覆盖前一阶段的工作树修改。

## 2. Time authority：before → after

| 属性 | Before | After |
|---|---|---|
| recipe id/family/order/stage count | `FDM::TimeRecipe` | `FDM::TimeRecipe` 为选择；`CompiledTimeRecipe` 校验并冻结执行值 |
| SSPRK3 stage abscissa、base/increment 权重 | `Time::Explicit::executeStage` 中三个静态数组 | `SF_compiledTimeRecipe.h`，Plan 与执行器消费同一编译值 |
| RK4 stage abscissa、中间更新及最终权重 | `Time::Explicit` 中字面量 | `CompiledTimeRecipe`；执行器按原有表达式顺序读取 |
| stage 次数 | raw recipe 与 provider 各自判断 | compiled recipe 校验，provider 在 step 前 `requireProviderStages` |
| stage time | execution 内按常数生成 | execution 按 compiled `c_i` 生成 |
| physical time commit | Plan 的 `time.commit` | 不变 |

**Implemented：** `NumericalCompiler` 将 built-in `forwardEuler`、`SSPRK3`、`classicalRK4` 编译为 `CompiledTimeRecipe`。`SolvePlanner` 读取编译值生成 StageLoop；`Time::Explicit` 读取其系数。`sonicSolver explain` 打印 stage abscissa、SSP 权重、RK4 最终权重。RK4 仍使用同一 field/scalar 表达式结构和运算次序；没有通用 RK 执行框架。

**Legacy：** `ResolvedSimulationSystem::timeRecipe` 和 `NumericalRecipeSet::time` 仍保留选择快照，供配置/provenance 与旧调用点读取，但不再是实际 stage 系数 authority。固定时间压力 outer predictor 仍是 `U_n + dt R(U_k)`；不能把每次 outer 误称作一次新的物理 Forward Euler 时间步。

**Unsupported：** DIRK、IMEX、BDF 和多 stage pressure-coupling。SSPRK3/RK4 与 PISO/SIMPLE/PIMPLE 组合不报告 Runnable，不静默改用 Forward Euler。

## 3. Operation provider：before → after

**Before：** `SF_providerResolver.cpp` 内按 conservative → constant-density pressure → Eulerian → IBM 的有序 `if/else` 选择实现；第一匹配具有隐式优先权。

**After — Implemented：** `SF_providerCatalog.h/.cpp` 登记 `flow.conservative`、`flow.pressure-operators`、`flow.rhie-chow`、`flow.eulerian-pressure`、`ibm.constraint` 的 capability 和匹配函数。Resolver 只从 `ExecutableEquationSystem`、`CompiledStateRealization`、`CompiledNumericalSystem`、编译后的 pressure schedule 与 `BuildCapabilities` 构造匹配上下文，然后逐个 OpId 查询 catalog。匹配 0 个为 `Unsupported`，1 个为 `Resolved`，超过 1 个为 `Invalid`，错误文字按排序后的 ID 列出全部候选，注册顺序不影响结果。`RuntimeReport` 保持与逐 OpId 绑定状态一致；explain 显示 `Invalid`。旧 provider ID 与实际 callbacks 保持不变。分布式压力仍使用 `flow.pressure-operators`/`flow.rhie-chow`，没有伪造 MPI 专用 provider 身份。

**界限：** 本轮 catalog 只解析 executable operations；gravity/MRF/turbulence 等 term provider 尚未迁入，这是 Phase 28B 范围。`BuildCapabilities` 可供 matcher 使用，HYPRE/MPI 的生产可用性仍由原有 `RuntimeRequirements`/validator 联合检查；不能把单个 catalog 匹配等同于 callback、拓扑和 backend 全部可用。

## 4. Pressure solve-policy authority

**Before：** provider 对 `ExecutionPolicyKind` 与 legacy `SolveStrategyKind` 作双重匹配。**After — Implemented：** runtime selection 只读取 `ExecutionPolicyKind`；`SolvePlanner` 在降低前校验 legacy strategy metadata 是否与 typed pressure policy 一致，冲突立即抛出明确错误。`strategyName` 只服务 explain/provenance。`SolveStrategyKind` 尚保留在 `SolveBlock`/旧配置中，状态为 **Legacy compatibility metadata**，不是第二个执行分派 authority。TimeRecipe、linear-solver policy 和 pressure preset 仍是互相独立的概念。

## 5. CMake 图和具体 adapter

修改前 target 图存在一个 SCC：`SF_equation → SF_discretization → SF_physics → SF_equation`，另有 `SF_equation ↔ SF_turbulence` 与 `SF_equation → SF_mesh → SF_physics → SF_equation`。修改后：

```text
SF_equation -> SF_discretization, core targets
SF_physics -> SF_equation
SF_turbulence -> SF_equation, SF_physics
SF_eulerianEquationAdapter -> SF_equation, SF_physics,
                              SF_turbulence, SF_mesh, SF_linearAlgebra
```

**Implemented：** 通用 `SF_equation` 不再链接模型、湍流或 mesh target；`SF_discretization` 的 INTERFACE target 不再无故传递 `SF_physics`。具体 multiphase/level-set coupling 实现从 `solver/equation/coupling` 移至 `models/physics/equationRuntime`，湍流耦合实现移至 `models/turbulence`。现有 Eulerian 相方程数值实现仍保留在 `solver/equation/eulerian` 源码目录，但作为明确过渡的 `SF_eulerianEquationAdapter` target 编译；这避免让 generic `SF_equation` 获得反向模型依赖，也没有把相方程数值循环交给物理模型 target。CMake 全 target 图 SCC 从 1 组降为 0；干净构建无 equation/physics/turbulence duplicate-library warning。尚有数值库的重复链接 warning，属于独立的静态链接元数据问题，本轮没有改动通量库数学或链接顺序。

`core/system` 未新增 models include/target 依赖。`SystemBuilder` 仍只组织 composition、transformation、numerics/plan/provider 编译，不包含新系数表或 MPI 错误处理。`BuildRequest` 仍含 level-set/turbulence 的具体 contribution spec include；若要消除需改造 model contribution 输入契约，故标 **Legacy**。`SystemContribution` 持有 `const RawEquationSystem&`，仍存在生命周期耦合；本轮按审计要求只记录，不额外制造 context/adapter。

## 6. MPI/communication failure contract

**Implemented：** scalar halo 在打包前检查缺失/重复 block view、值数量及非有限值，并以 backend `allRanksAgree` 使单 rank 输入错误在双方形成同一失败决定；canonical face 同样预检 owner、block/rank、face 索引、方向和本地 canonical flux。新增 2-rank fixture 验证 rank-local scalar 长度错误、重复 block 映射与 canonical-face 引用错误：报出 rank/context，所有 rank 结束，没有挂起。`CommunicationPlan` 的非法插值 stencil 不再深层 `std::exit`。GlobalDof COPY owner/replica 和若干打包/偏移错误改为带原始上下文的异常，已有 application catch 先输出 rank，MPI backend 在异常展开时 fail-stop。

`SF_haloExchange.cpp` 的直接 `std::exit(1)` 从 41 减到 25，`SF_communicationPlan.cpp` 从 1 减到 0。剩余 25 个涉及其他 halo、张量插值、state/flux 后处理与不同通信阶段，标 **Legacy**；没有机械全局替换，因为需要逐条证明 collective-safe 或 fail-stop。preflight 增加控制通信但不改 payload、owner、COPY/SUM 或 flux 算术。原有 MPI Peskin singular mass diagonal 仍是既存问题，本轮未用数值 fallback 掩盖。

## 7. Architecture guards 与测试

`tools/check_architecture.py` 新增 CMake target link cycle 检查，禁止把具体 provider family selector 放回 `SF_providerResolver.cpp`，禁止恢复显式执行器的局部 SSP 系数数组。允许的 include debt 仍为 2 条，无新增 `core → solver` 或 `infrastructure → solver` 边。

测试新增：`providerCatalogResolution` 覆盖单一/零/多候选、serial/distributed 和缺少 HYPRE capability；`haloFailureContract` 覆盖 2-rank 通信输入失败。扩充 `explicitStageMathematics` 覆盖全部 built-in recipe 系数与 provider-stage mismatch；扩充 `pisoArchitecture` 覆盖 typed/legacy policy 冲突和多 stage pressure 组合 Unsupported。

| 验证 | 修改前 | 修改后 |
|---|---:|---:|
| 完整普通构建 | 通过 | 通过 |
| 干净 Debug 源码、CLI、GUI 与测试目标构建 | — | 通过（309 个基础目标 + 28 个测试目标构建步骤） |
| CTest（含 host MPI/HYPRE） | 22/22 | 24/24，普通与干净 Debug 树均通过 |
| architecture checker | 2 条既有 allowlist | 2 条既有 allowlist，0 新违规 |
| CMake SCC | 1 组 | 0 组 |

独立冻结二进制对照（相同 case、步数与执行配置）显示：4-rank Sod、single WENO7、serial Ghost IBM、Eulerian、legacy pressureConstraintPiso、constant-density PISO 的 diagnostics 与 VTS SHA-256 全部逐字节一致。legacy PISO 最终 SHA-256 仍为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`。2/4-rank fractional DLM 和 velocity forcing 的 20-step 输出哈希也与修改前一致。Phase 28A 分布式 pressure 脚本在完整 CTest 中通过：cavity/channel 的 2×1、1×2、2×2 分解及 SIMPLE/PIMPLE 覆盖均沿用原阈值。一次增量组合测试曾出现不可复现的 `explicitStageMathematics` Bus error；单项重跑、重复组合及之后普通和干净 Debug 的完整 CTest 均通过，未观察到数值首差异。

## 8. Phase 28B 准备状态与遗留项

**Implemented：** 当前 operation provider 可以在不修改 central family selector 的情况下登记并被确定性匹配；时间 stage 数学由编译值驱动；pressure typed policy 是运行时唯一控制语义；原有 equation/physics/turbulence target 环已消除；已覆盖的 rank-local 通信输入错误能可见地、有限时失败。

**Interface-only：** catalog matcher 有 build/distributed 上下文，但 term-level contribution/provider lowering 尚未实现；不能据此宣称 gravity/MRF/turbulence 可自由组合执行。

**Legacy：** 固定时间 predictor 的时间解释、`BuildRequest` 具体 model spec、`SystemContribution` 对 raw 的引用、Eulerian 专用数值 adapter、25 处 halo 直接退出、二维 inactive-z 压力 gauge 与单 reference row、数值库重复链接 warning。

**Unsupported：** DIRK/IMEX/BDF、多 stage pressure coupling、空 pressure-row rank、任意多源块 pressure 拓扑、现有 MPI Peskin 质量对角奇异路径。后续应独立处理 pressure nullspace/连通分量与 term provider composition；本轮不启动 Phase 28B。
