# Phase 25B — Provider Resolution 与 Shared-Pressure 收口

## 1. Scope

**Implemented**：本轮只收敛运行前架构编译链。数值公式、CFL、线性求解容差、RK stage time、IBM 数学及 MPI owner/COPY/SUM/canonical 规则未改。生产回归使用相同 case、步数和执行配置。

| AGENTS.md 检查项 | 本轮结论 |
|---|---|
| 层与数学作用 | `core/system` 中立 IR；`solver/system` 编译 HOW。属于“改怎么解方程”的绑定/调度收口，不添加物理 term 或 equation。 |
| ownership before → after | operation 的 concrete provider 字段、builder 缺失清单 → 唯一 `ResolvedOperationBinding`；Eulerian plan 的 specialized compiler → shared-pressure fragment。physical state 与 solver workspace ownership 不变。 |
| execution order before → after | Eulerian 34 行 Plan 树及 28 个 OpId 的次序相同；density/IBM stage、boundary、halo 顺序不变。 |
| public API before → after | 删除 `compileEulerianPimple()` 与 `strictCoreTermCoverage` 参数；`PressureAlgorithm` typed 名称改为 `PressureCouplingPreset`，native `algorithm:` 输入仍可解析。 |
| dependency before → after | models 依赖 solver/system 编译实现 → models 只依赖 core/system 值类型；编译器仍留在 solver/system。 |
| MPI 与数值语义 | canonical COPY、GlobalDof SUM、stage time、pressure correction、IBM 算法均未变；验证见第 19–23 节。 |
| deferred | constant-density PISO、SIMPLE/PIMPLE fixed-point provider 及第 24 节的旧链接债。 |

## 2. Phase 25 starting state

**Legacy**：`ExecutableOperation` 直接携带 `provider`；`SF_systemBuilder.cpp` 拼装 missing-operation 报告；Eulerian shared-pressure 由 `compileEulerianPimple()` 在通用 planner 内生成；`strictCoreTermCoverage` 用一个含约束与否的布尔量放宽 recipe 检查；模型贡献头文件依赖 `solver/system`；`SolverPropertiesConfig` 混装耦合、相输运、参考压力及线性求解配置。Eulerian 相输运 CFL 还从 raw pressure workflow 读取。

## 3. Operation/provider separation

**Implemented**：`src/core/system/SF_equationIR.h` 的 `ExecutableOperation` 只保存 OpId、stage、typed `OperationCapability`、provenance 和可选 recipe 消费声明。压力 formulation 在 `src/solver/system/SF_transformation.cpp` 声明 predictor、pressure、velocity/flux correction 能力，不再指定 `flow.conservative`。一个 OpId 可因 state realization 不同解析为不同结果。

## 4. ResolvedOperationBinding design

**Implemented**：`src/solver/system/SF_runtimeRequirements.h` 以 `ResolvedOperationBinding {operation, provider, status, reason}` 保存编译后的单一绑定；状态为 `Resolved` 或 `Unsupported`。绑定位于 RuntimeRequirements，执行器只读取，不重新判断 preset 名称。

## 5. Provider resolution algorithm

**Implemented**：`src/solver/system/SF_providerResolver.cpp` 从 CompiledSolvePlan 递归取得 OpId，查 ExecutableOperation，并根据 typed requirements、`CompiledStateRealization`、时间 recipe 和已有 numerical capability 选择现存 conservative、Eulerian shared-pressure 或 IBM constraint provider。未声明 operation、缺少匹配 provider、或 fixed-point predictor 不可用，均保留该 OpId 并给出原因。MPI/HYPRE 的主机可用性仍由 BuildCapabilities 编译到独立 RuntimeRequirement，并在 validation 阶段检查；它不伪装成另一种数值 provider。

## 6. Runtime requirement derivation

**Implemented**：`src/solver/system/SF_systemBuilder.cpp` 先完成 plan 和 operation bindings，再从已解析绑定派生 provider requirement；`RuntimeReport` 由 `reportOperationBindings()` 生成。`src/solver/system/SF_systemValidator.cpp` 校验每个 plan OpId 恰有一条绑定、已解析 provider 存在于 requirements，且 missing list 与 Unsupported bindings 精确一致。Coupling preset 的数学适用状态与 numerical provider 的 Runnable/Unsupported 状态分开报告。

## 7. Eulerian shared-pressure migration

**Implemented**：shared-pressure transformer 声明 Eulerian `ee.*` executable operations；`src/solver/system/SF_pressureCoupling.cpp` 的 `sharedPressurePlanFragment()` 提供 outer、pressure-corrector、non-orthogonal 三层循环及完整既有 operation 顺序。通用 planner 只降低 fragment 值。可选 turbulence operations 仅在对应方程存在时声明和调度。

## 8. compileEulerianPimple removal

**Implemented**：删除 specialized `compileEulerianPimple()`。迁移前后 `COMPILED SOLVE PLAN` 的 34 行，以及 `REQUIRED OPERATIONS` 的 30 行逐行一致；`EE.begin`、`EE.commit` 的 node identity 也保持不变。

## 9. SolvePlanner neutrality

**Implemented**：`src/solver/system/SF_solvePlan.cpp` 不再按 SIMPLE/PISO/PIMPLE 或 `ee.pressure.*`、`ee.momentum.*` 分支。它处理 Sequence、Loop、StageLoop、fragment leaf 与 generic policy leaf；没有可用 fragment 或 explicit stage operation 时，不会把一个活跃 policy 默默当成可运行的显式阶段。显式时间 recipe 的现有 fused-stage 入口仍属于既定 compiled time contract。

## 10. Recipe consumer validation

**Implemented**：`src/solver/system/SF_numericalCompiler.cpp` 记录 `CompiledRecipeBinding`，区分 equation term 和 executable operation 消费者。选择了 diffusion recipe 却没有真实消费者会 fail-fast；约束存在本身不再豁免检查。自动默认 diffusion recipe 仅在当前 core fluid equation 的 diffusion term 确实由该 compiler 绑定时产生。Eulerian 专用 phase diffusion 数学仍由现有 phase executor 实现，不冒称它消费 core diffusion recipe。

## 11. strictCoreTermCoverage removal/replacement

**Implemented**：删除 `strictCoreTermCoverage` 参数和含约束即整体放宽的分支。测试分别覆盖 equation-term 消费、operation 消费、未消费 recipe 报错、以及存在约束时仍报错。`OperationRecipeRole` 是显式声明机制；本轮没有为了让 validation 通过而给未使用 recipe 虚构生产消费者。

## 12. core/system IR migration

**Implemented**：中立值类型移至 `src/core/system/`：`SF_expression.h`、`SF_equationIR.h`、`SF_solveProgram.h`、`SF_planFragment.h`、`SF_transformationTypes.h`、`SF_systemContribution.h`。`SystemBuilder`、TransformationPipeline、NumericalCompiler、SolvePlanner、ProviderResolver 和 validator 保留在 `solver/system`。核心 IR 不依赖 solver 实现；旧同名 forwarding header 已删除，避免新增 include basename 歧义。

## 13. models/system dependency graph

**Implemented**：物理源项、level-set、turbulence 和 IBM 模型输出中立 `SystemContribution` 值，`SystemCompositionBuilder::applyContribution()` 统一验证并合成。`src/models/**` 不再 include/link `solver/system`；`SF_turbulence → SF_solverSystem → SF_turbulence` 的直接静态库循环已消除。**Legacy**：`SF_equation ↔ SF_turbulence` 的独立旧链接关系仍会使 Darwin linker 提示重复静态库；本轮未移动其数值实现。

## 14. PressureConfig decomposition

**Implemented**：`PressureAlgorithm` 改为 `PressureCouplingPreset`；typed config 将 `PressureCouplingConfig`、`PhaseTransportConfig`、`PressureReference` 和 `LinearSolverSet` 各保存一次，集中归入 `PressureCorrectionConfig`。native 输入仍可写 `algorithm: PISO`，解析后只产生 typed `preset`。原有数值参数默认值和容差不变。

## 15. Eulerian raw config cleanup

**Implemented**：`phaseConvection` 是相输运数值格式选择，`phaseSourceCfl` 是相源项 timestep 限制；二者从 pressure coupling 移至 phase transport config，并编译进 `CompiledNumericalSystem::phaseTransport`。Eulerian timestep 使用编译值；当前仅有 Upwind provider，其他格式显式失败。raw config 仍作为初始化及已有方程/线性求解参数载体，但不再作为这两项 timestep 选择的第二 authority。

## 16. Architecture guards

**Implemented**：`tools/check_architecture.py` 继续禁止旧 operation/provider 列表及旧压力 fallback，新增 `compileEulerianPimple`、planner 算法级字符串分支、`ExecutableOperation` concrete provider、`strictCoreTermCoverage`、模型 include/link `solver/system` 的守卫。当前检查通过：622 个源文件；现存 allowlisted 依赖为 2 条，未新增违反方向的边。

## 17. Constant-density PISO status

**Unsupported**：constant-density `U/p` 与 pressure-multiplier realization 可以生成有效方程、operations 和 PISO plan，但 `pressure.solve` 等没有匹配的 single-fluid constant-density numerical provider。resolver 明确给出缺失 OpId 和 realization 原因；未借用 conservative full-`dt` momentum update 充当 predictor。

## 18. Unsupported providers

**Unsupported**：SIMPLE/PIMPLE dedicated fixed-point predictor、constant-density PISO、以及未声明或未实现的自定义 operation 继续 fail-fast。`RuntimeReport` 列出的 missing operations 来自 unresolved bindings，而非手写字符串清单。**Legacy**：Eulerian phase kernels 和 IBM constraint kernels 仍是现有专用数值实现，但现在通过已编译绑定接入 OpRegistry；本轮没有改写其数学。

## 19. Sod numerical invariance

**Implemented**：host 上 `mpirun -np 4 /tmp/sonic-phase25b/sonicSolver run --steps 20 test/Sod/sodCase` 成功。迁移前/后 40 条 closure 与 timestep 诊断逐行相同，含每步 `dt`、`time`、`min(rho/p/T/c)`；最后 `time=9.818765e-03`、`dt=4.572193e-04`。未调整高阶格式及后续已知负压问题。

## 20. cylinderFlow numerical invariance

**Implemented**：`cylinderFlowGhost` 20 步成功；迁移前/后 120 条 RK4 stage/final closure 与 timestep 诊断逐行相同。最后 `time=8.778858e-03`、`dt=4.386629e-04`。Ghost/ILW 仍是 boundary closure。

## 21. pressureConstraintPiso SHA regression

**Implemented**：`tools/check_piso_regression.py` 通过。最终 VTS 的 SHA-256 为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`，与迁移前及脚本冻结结果完全相同。

## 22. Eulerian regression

**Implemented**：`eulerianEulerianCase` 2 步成功。迁移前/后 2 条完整 step summary 逐行相同；第一步 `dt=3.624299e-05`，第二步 `dt=3.624215e-05`，`phaseMass=2.358562e+02`，pressure residual、HYPRE rebuild/solve 计数一致。Plan 树与 OpId 列表逐行相同。架构测试另外冻结了 28 个 OpId 的顺序及三层循环计数。

## 23. Clean build + CTest

**Implemented**：按指定 clang++/Ninja 参数在 `/tmp/sonic-phase25b` 完成干净构建（271 个 build steps；最后增量重建通过）。完整 CTest 为 **14/14 通过**；`python3 tools/check_architecture.py` 通过。构建无 compile/link failure；Darwin linker 的旧重复静态库警告见第 13 节。

## 24. Remaining architecture debt

**Legacy**：`SF_equation ↔ SF_turbulence` 链接关系、Eulerian phase 方程的专用 numerical lowering、部分压力配置的兼容输入及现存 numerical backend 尚未泛化。Sod/Ghost 的迁移前数据只保存了文本诊断，故本报告证明的是**记录精度下逐行完全相同**，不宣称未保存前置哈希的场文件逐字节相同。**Interface-only**：operation recipe 消费声明已可表达 model/operation 消费，但本轮没有制造虚假的生产绑定。

## 25. Recommendation for Phase 26

**Unsupported**：下一阶段应独立实现 **Constant-Density Single-Fluid PISO Provider**：由真实的 fixed-time predictor、pressure linear system、velocity/flux correction kernels 满足本轮生成的 OpId/capability contract，并建立数值基线。不要在该阶段顺手替换现有 conservative、Eulerian 或 IBM 数值实现。
