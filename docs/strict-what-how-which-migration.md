# Strict WHAT / HOW / WHICH / Compiler migration report

STATE 拆分及 builtin 参数目录的后续实现见 [四模块迁移报告](four-module-state-migration.md)；本文件保留此前 WHAT/HOW/WHICH 阶段的记录。

2026-10-01。继续现有 Equation Registry + Ordered Execution 架构，收口职责，没有建立新 solver family 或另一套 SimulationPlan。本报告描述本轮实际生产实现；工作区原有大量未提交修改已保留，不将全部 workspace diff 归因于本轮。

## 1. Before

审计发现的跨层 authority：

| 旧位置 | 越界 |
| --- | --- |
| applyPressureExecution | 同时变换 HOW 并插入 PressureMomentum/PressureCorrection 等 numerical bindings。 |
| ExecutionProgram.temporal | source HOW 保存时间方法 WHICH。 |
| Target | source target 保存 symbols/display、workspace、CompiledResourceBinding。 |
| EquationCall.requiresWorkspace | source HOW 保存具体 numerical workspace requirement。 |
| ExecutionScope.before/after | source HOW 保存 PressurePrepare/IterationBegin/CorrectionCommit/TimeCommit 等 runtime IDs。 |
| PressureConstraintTransformer | 同时注册数学关系、声明 numerical operations 和 compiled equation resources。 |
| NumericalCompiler | 按 constantDensity/pressureMultiplier 选择 RhieChow，按 equation ID 搜索 diffusion consumer，按 output==U 选择源项/工作区约定。 |
| execution compiler | 直接填入 ExplicitStageExecute；调用前修改 source target 来填资源。 |
| SystemBuilder | 用 StateRealization flags 声明时间/压力操作，决定 structured pressure HOW 与 legacy fragments。 |
| numerical precedence | 遇到低优先级重复时可能在看到 occurrence override 前抛错，结果依赖注册顺序。 |
| explain | WHICH section 混入 compiled storage 和 numerical micro-topology。 |

旧 ExecutionPolicy/PlanFragment 仍有真实 legacy 用户，没有将重命名冒充删除。density 的 physical Commit 也补为显式 source HOW，避免 temporal provider 独立决定提交位置。

## 2. Final authority graph

```text
MODULES
  ├── WHAT  ── EquationRegistry ──────┐
  ├── HOW   ── ExecutionProgram ─────┼── const inputs / composition complete
  └── WHICH ── NumericalSelection ───┘            │
                                                ▼
                                             Compiler
                                                │
                                   capability / Storage Binding
                                                │
                                         CompiledSolvePlan
                                                │
                                           PlanExecutor
                                                │
                                            providers / kernels
```

WHAT defines equations. HOW invokes equations. WHICH selects numerical implementations. Compiler binds and lowers; it does not make solver decisions. **The compiler is not a fourth authority.**

EquationRegistry 与 source execution/selection 可以独立创建和检查。Raw system 是 composition 前的历史视图；最终 executable registry 是编译所读的数学定义。compiled operations/resources 是输出 metadata，不会反向注册新数学方程。

## 3. Pressure HOW / WHICH split

实际接口：

```cpp
void applyPressureExecution(ExecutionProgram&, const CouplingPresetRequest&);
std::vector<NumericalBinding> pressureNumerics(const CouplingPresetRequest&);
```

WHAT 通道：PressureConstraintTransformer 注册 pSimple/correctU/correctP/correctFluxp 与 fixed-time mathematical relations。

HOW 通道：applyPressureExecution 找已有 Physical U 的默认 momentum occurrence，把 target 改为 Working U，赋予 predictor occurrence，插入 generic loops、correction calls 和最终 Commit，保留其他模块 entries 与数学定义。它不接收 numerical collection，不引用 numerical provider IDs，也没有 OpIds。

WHICH 通道：pressureNumerics 独立贡献 momentum@predictor→PressureMomentum、pSimple→PressureCorrection、correctU→VelocityCorrection、correctP→PressureUpdate、correctFluxp→FluxCorrection；fixed-time preset 再贡献三个 relation providers。它不接收或修改 ExecutionProgram。

composeContributions 是模块 composition 协调点，分别调用上述独立通道；native applicability 来自已贡献的默认 equation occurrence，不来自 StateRealization。legacy numerical declarations 在该边界显式标注，native 不进入它们。

PISO 的 source tree 为 momentum→U*、pressure Loop（包含 nonOrthogonal Loop 与 correction calls）、Commit；SIMPLE 使用带 termination signal 的 outer Loop 和一次 correction；PIMPLE 在 outer 内再包含 pressure Loop。loop counts 只有这一份 source HOW authority。

## 4. Temporal ownership

source ExecutionProgram.temporal / TemporalMethodBinding 已删除。NumericalSelection.recipes.time 保存完整、不可变的 builtin recipe。

selectNumerics 在 composition 阶段保存 configured time choice。compiler 查找已选 temporal implementation 并生成 stage topology、coefficients/workspace realization；StageLoop 是 lowering 结果，runtime 读取 frozen CompiledTimeRecipe。源程序没有 RK4/Euler 专用 node。

density source HOW 明确有最终 Commit。temporal fragment 不再凭空插入 state/clock Commit；compiler 保留 source Commit 的位置，ConservativeResidual provider 实现其 state/clock operations。当前 fused capability要求 Commit 位于 root 最后；nested Commit 和缺少 provider implementation 的 Commit 明确 fail fast，不能 silently drop。

压力 providers 声明现有 ForwardEuler 与各自可实现的 target-kind capability。pressure + SSPRK3/RK4 编译时返回具名 Unsupported provider/occurrence/temporal error，不改变 target，不替换 recipe，也不转入 legacy path。

## 5. Target ownership

| 层 | 实际类型和职责 |
| --- | --- |
| Source HOW | Target {symbol, TargetKind}；Physical/Working/Correction/Workspace 仅是语义。 |
| Compiler/provider output | CompiledTarget {symbol, kind, workspace, resources}；storage、offset、components、access/synchronization 在这里绑定。 |
| Runtime | SolvePlanNode/ExecutionContext 携带 CompiledTarget；kernel 消费 typed kind 和 frozen decisions。 |

`targetFromSyntax("U*")` 返回 U/Working，`targetFromSyntax("p'")` 返回 p/Correction。压力数学 AST 的 correction unknown 仍是 pPrime；PressureCorrection provider 将 source p/Correction 绑定到 pressureCorrection workspace。runtime 不解析星号或撇号。

physical target 由已声明的 storage-bound UnknownDescriptor backing；Working/Correction/Workspace 由 provider 实现。generic DirectEvaluation/LinearEquation workspace 用 occurrence 建立独立 identity；pressure provider 使用已有 numerical workspace。源码不保存其名字。

编译器检查 provider 不得改变 source equation、occurrence 或 target semantics。已有 PressureOperators 的 workingVelocity storage/publication 保留，未新增 Field、clock 或 physical state registry。formula cell runtime adapter 改为接收 CompiledTarget，并明确拒绝未实现的非 Physical materialization，避免误写 physical storage。

## 6. Pressure transformer

native PressureConstraintTransformer 现在仅：

- 匹配已声明的 momentum/U/pressure constraint。
- 注册 mathematical pPrime symbol 与压力/correction/fixed-time relations。
- 保存 mathematical generated-operator/provenance view。

它不再要求具体 U storage backing、不设置 correction storageKey、不创建 CompiledResourceBinding、不声明 runtime operations、不注册 compiled equation resources。

PressureMomentum/PressureCorrection 等 numerical implementations 位于 SF_builtinProviders.cpp。它们实现 target storage、capability、local numerical fragment 和 operation declaration。所有既有 pressure arithmetic 留在 PressureOperators/RhieChow/HYPRE 下。

每个旧 hook 的分类和去向：

| 旧 source hook | 当前所属层 |
| --- | --- |
| PressurePrepare | pressure provider root workspace lifecycle。 |
| PressureStepBegin | fixed-time provider root snapshot lifecycle。 |
| PressureIterationBegin/End | provider 在拥有对应 predictor 的 generic termination Loop 内附着 compiled lifecycle。 |
| PressureCorrectionCommit | FluxCorrection compiled fragment 的 publication/bookkeeping。 |
| PressureStepCommit / TimeCommit | source Commit 的 provider implementation。 |
| relaxIterate / restoreFlux / checkConvergence | 仍是 WHAT relations + source EquationCall + WHICH numerical implementation。 |

没有 source before/after，也没有 public pressure lifecycle node。fake callback tests验证了迁移后的顺序、重复次数、终止和 Commit；这不证明 CFD numerical equivalence。

variable-density pressure 的旧 compiled declaration 移到显式 Legacy::lowerPressure。shared-pressure transformer 更名 LegacySharedPressureTransformer，其旧混合职责作为未完成迁移边界保留，没有声称 Eulerian 已 native 化。

## 7. Compiler neutrality

SF_methodObjects.cpp 的 generic compilation 现在：排序 source 副本、resolve occurrence/binding、调用选中 provider、验证 target/temporal/workspace capability、lower generic scope，并保留 provenance。fusion backend 来自 provider declaration，检查所有 fused calls 相容，不再写死 ExplicitStageExecute。

SF_numericalCompiler.cpp 的 spatial binding 只消费 provider 声明的 spatialTerms、source kernel requirement 和 residual workspace。移除了：

- momentum/continuity/energy ID 的 permanent compiler scope rule。
- constantDensity/pressureMultiplier→RhieChow inference。
- output==U→primitive source / pressure workspace 分派。
- compiler 自动选择 diffusion default 和 pressure numerical configuration。
- generic compiler 导入 LegacyEquationRole/TermKind 来猜旧 mathematical inputs。

NumericalSelection 在 composition 前端已保存 recipes、pressure numerical inputs、term providers 与 operator overrides；legacy spatial inputs 由显式 Legacy adapter提供。compiler 接收 const HOW/WHICH，providers 接收 const WHAT/call/binding。编译排序在独立 source copy 上进行；source AST、order、target 和 selection 不变。

binding precedence 是 global < equation < occurrence。只有最高匹配优先级的重复报错；低优先级 duplicate defaults 不再使更具体 override 随注册顺序失效。不存在 unknown provider / unsupported capability 的 fallback。

允许保留的 domain recognition：

| 文件 | 原因 |
| --- | --- |
| SF_builtinProviders.cpp | numerical provider 验证 U/pPrime 数学输出、pressure storage 和既有 kernel capability；这是 provider matching。 |
| SF_builtinNumerics.cpp | composition numerical defaults，包含 pressure provider 的 RhieChow/HYPRE configuration；不修改 WHAT/HOW。 |
| SF_providerResolver.cpp / SF_providerCatalog.cpp | builtin runtime backend capability adapter，核对所需 compiled operation、字段 layout、recipe 和 host backend；不创建方程、target、循环或改变 selected temporal method。variable/shared 分支仍为 legacy capability debt。 |
| SF_executionComposition.cpp / SF_pressureCoupling.cpp | module composition/matching 与显式 legacy adapters；压力 domain 属于模块，不藏在 generic compiler。 |
| SF_solvePlan.cpp | generic lowering 与明确 LegacyExecutionPolicy/LegacyPlanFragment adapter；旧 explicit inference 仍在 compileLegacyExplicitStep，native 不使用它。 |

两个 generic compiler 文件没有 momentum/continuity/energy/U/rhoU/rhoE/p 的 hard-coded authority规则，也没有 PISO/SIMPLE/PIMPLE/RhieChow 或 StateRealization dispatch。静态检查覆盖这些高价值边界。

## 8. SystemBuilder

build 的 native职责为 collect module WHAT/HOW/WHICH → mathematical composition → module coupling composition → scope ordering →保存 NumericalSelection → invoke compiler → capability/runtime report。

删除了 Builder 中用 conservativeTransportedMass/pressureMultiplier/phaseTransportedState 决定 native pressure topology和时间/压力 operation declaration 的逻辑。StateRealization 仍编译 physical storage/role view；它不传给 native HOW applicability。

Builder 仍有 explicit template/model composition、mathematical applicability validation、runtime services 与 legacy reporting。这些不构成 native solver lifecycle dispatch。thermodynamic runtime requirement 读取已编译 provider requirement，不按 constantDensity 重新选择 HOW。

## 9. Source HOW vs compiled IR

Source ExecutionProgram：root + explicitly pending legacyEntries。source scopes只有 generic kind、equation call、symbol/kind target、order、origin、loop counts/minimum/termination signal。source没有 temporal method、workspace requirement、compiled resource、raw OpId、before/after或 consumedPolicies。

CompiledSolvePlan：operation tree、typed CompiledTarget、numerical fragments/lifecycle、temporal StageLoop、compiled equation occurrences、recipe/workspace capability与 metadata views。compiled plan不是第二份 authored HOW。

explain 先打印独立 WHAT / EQUATIONS、HOW / EXECUTION、WHICH / NUMERICS source selections，再打印 COMPILED OCCURRENCES / STORAGE / PROVIDERS 和 COMPILED SOLVE PLAN。pending entries标为 Legacy execution authority，legacy numerical input adapter明确标注。

## 10. Legacy boundaries

| 边界 | 状态 |
| --- | --- |
| density single-fluid / constant-density PISO/SIMPLE/PIMPLE | Implemented native authority chain；保留既有 numerical kernels。 |
| Eulerian phase transport/shared pressure | Legacy isolated：LegacySharedPressureTransformer、LegacyPlanFragment、specialized stepper 和原数值 execution；pending calls不独立执行。 |
| turbulence | common equation/contribution interface + pending HOW；Legacy numerical backend仍拥有 transport timing，未双重执行。 |
| IBM Ghost/ILW/forcing/projection/DLM/KKT | Legacy execution/backend 边界保留；未改数学、constraint ownership或 MPI。 |
| AssemblyPlan / TermKind / LegacyEquationRole | Legacy mathematical/backend compatibility；native single-fluid不创建平面副本。 |
| LegacyExecutionPolicy / LegacyPlanFragment | 原类型显式重命名并保留；真实 legacy用户仍需要，重命名不计入删除。 |
| variable-density pressure | Legacy::lowerPressure 和 legacyCouplingPlanFragment，明确旧 backend capability。 |
| Legacy FormulaCall/FormulaMode | legacy numerical proof/compatibility interface，不是 source HOW。 |
| arbitrary user equations/targets | contribution与provider接口可表达；通用 production kernel/任意 YAML parser未完整实现，不支持就 fail fast。 |
| nonorthogonal/mixed iterative temporal combinations | 现有 provider capability限制仍保留，未通过 fallback 掩盖。 |

## 11. Deleted authority

实际删除的 source authority：ExecutionProgram.temporal、TemporalMethodBinding、source Target 的 display/symbols/workspace/resources、EquationCall.requiresWorkspace、ExecutionScope.before/after、source consumedPolicies；pressure HOW 的 numerical-binding mutation；native数学 transformer的 numerical/runtime/resource declarations；generic compiler按方程/output/state flags选择 backend/defaults 的规则；temporal provider凭空创建 density physical Commit的位置决定。

没有宣称删除所有 legacy schedule、AssemblyPlan、TermKind、Eulerian backend、provider matching或旧 mathematical adapter。Legacy类型重命名与 implementations移动是职责隔离，不是功能删除。

## 12. Files changed

| 主要文件 | 最终职责 |
| --- | --- |
| core/system/SF_solveProgram.h | source Target/HOW 与独立 CompiledTarget/compiled plan。 |
| core/system/SF_numericalBinding.h | source NumericalBinding contract。 |
| solver/system/SF_numericalSelection.h | 完成 composition 的 source WHICH。 |
| solver/system/SF_pressureCoupling.* | pure pressure HOW 与独立 numerical contribution。 |
| solver/system/SF_transformation.* | native数学 transformation；明确 legacy shared边界。 |
| solver/system/SF_builtinProviders.cpp | 从 generic compiler 移出的 builtin numerical implementations、storage、capability与 lifecycle。 |
| solver/system/SF_builtinNumerics.cpp | numerical defaults/configuration composition。 |
| solver/system/SF_legacyNumerics.* | explicit legacy spatial selection / pressure numerical declarations。 |
| solver/system/SF_executionComposition.* | module/legacy composition 协调。 |
| solver/system/SF_methodObjects.* | deterministic binding / generic lowering；binding冲突和target invariants。 |
| solver/system/SF_numericalCompiler.* | neutral operator binding和 compiled result。 |
| solver/system/SF_systemBuilder.cpp / SF_singleFluidPreset.cpp / SF_presets.cpp | freeze boundary、source Commit、default occurrence。 |
| solver/system/SF_solvePlan.cpp / SF_providerResolver.cpp | compiled lowering和 builtin capability matching。 |
| solver/system/SF_formulaStorage.* / solver/run/SF_planExecutor.h | runtime compiled target interface。 |
| solver/system/SF_systemPrinter.cpp / SF_resolvedSimulationSystem.h | source WHICH与 compiled diagnostics分别展示。 |
| core/system/SF_planFragment.h / legacy contribution consumers | LegacyExecutionPolicy/LegacyPlanFragment命名适配。 |
| test_methodObjects / formulationArchitecture / pisoArchitecture / termRecipes | 新职责边界检查和旧断言迁移。 |
| tools/check_architecture.py / src/ARCHITECTURE.md / README 导航 | structural guards、ownership matrix和当前架构正文。 |

Ownership before/after：数学和physical state authority保持；sourceHOW不再拥有 compiled storage/provider/runtime decisions；sourceWHICH独立拥有数值选择；provider compilation拥有 storage与局部 numerical lifecycle。

Execution order before/after：同一 pressure操作序列从source hooks移入compiled provider fragment；Commit的source位置明确。公开API变为纯 HOW transformation、独立 numerical selection和 source/compiled target。依赖方向为 composition→compiler→compiled plan→executor→provider/kernel。

没有改写 pressure/RhieChow、WENO/TENO、RK coefficient、source数学或 MPI/halo/GlobalDof arithmetic；PressureOperators 的 numerical config 类型名和 Legacy interface 命名有同步适配。本轮属于“改怎么解 equation”的执行/绑定职责收口，不增加物理 term 或 equation。没有删除 A-type unsupported-combination checks，没有修改B-type numerical safety-net。除 numerical selection冲突规则和更早的具名capability error外，不有意改变算术或并行语义；numerical equivalence尚未验证。

## 13. Architecture checks

通过的实际命令：

- `cmake --build build -j 6`：现有 configured targets 全部编译/链接。
- `python3 tools/check_architecture.py --quiet`：dependency/authority checks。
- `git diff --check`：无 whitespace error。
- `./build/test_methodObjects`：generic IR、immutable sources、target/occurrence、precedence、fake pressure lifecycle。
- `./build/test_formulationArchitecture`：数学 composition / pressure provider / native registry。
- `./build/test_pisoArchitecture`：native与legacy capability、contribution、provider和plan contracts。
- `./build/test_termRecipes`：operator recipe binding、unused recipe、AST source terms。

Build log：`/private/tmp/sonic-strict-build.log`。Architecture checker 的 quiet 模式现在失败时也输出错误，避免掩盖实际失败。

这些检查覆盖：source冻结、同一方程不同target/occurrence、scope-local sorting、绑定优先级、typed target parser/storage、generic termination、PISO/SIMPLE/PIMPLE provider生命周期顺序、unsupported temporal capability和 source关系扩展。

CLI explain 的五个 case 均返回 runnable，并显示分开的 source/compiled sections：density `test/Sod/sodCase_weno7`、PISO `test/pressure/constantDensityPiso`、SIMPLE/PIMPLE 临时 copies `/private/tmp/sonic-how-simple` / `/private/tmp/sonic-how-pimple`、Eulerian `test/eulerianEulerianCase`。日志为 `/private/tmp/sonic-strict-explain-{density,piso,simple,pimple,eulerian}.log`。explain只解析/编译，不推进 CFD 时间。

Numerical regression intentionally not run in this architecture phase.

未运行 CFD timesteps、长算例或完整 numerical regression；不声称 numerical equivalence。结构测试有fake callbacks与小型model evaluator检查，它们不替代数值基线/MPI执行回归。

### Required self-review

| 问题 | 实际答案 |
| --- | --- |
| EquationRegistry可不依赖HOW构建？ | Yes，core EquationRegistry和数学 contribution独立。 |
| HOW可不知provider检查？ | Yes，只有semantic target与generic scope。 |
| WHICH可检查而不改math/order？ | Yes，NumericalSelection及独立pressureNumerics。 |
| compiler接收已composed三者？ | Yes，const math/provider接口、const HOW/WHICH；sorting用副本。 |
| compiler不选择PISO/SIMPLE/PIMPLE？ | Yes，module composition已提供topology。 |
| compiler不决定momentum predictor？ | Yes，source target.kind=Working。 |
| applyPressureExecution只修改HOW？ | Yes，signature/body静态检查。 |
| temporal在sourceHOW之外？ | Yes，NumericalSelection.recipes.time。 |
| sourceTarget不存compiled storage？ | Yes，仅symbol/kind。 |
| sourceHOW无raw runtime OpIds？ | Yes，无before/after；generic signal是控制语义。 |
| PressureConstraintTransformer不造runtime ops？ | Yes，native class仅math；LegacyShared仍为显式debt。 |
| native Builder不选择solver family？ | Yes，HOW来自module occurrences；StateRealization用于storage/report。 |
| 同一momentum可不同target重复调用？ | Yes，EquationRef独立于occurrence/target；实际执行需selected provider capability。 |
| provider选择与order独立？ | Yes，WHICH按稳定occurrence binding；不重排source。 |
| PlanExecutor只执行compiled decisions？ | Yes，generic tree、targets、signals和callbacks。 |
| legacy隔离明确？ | Yes，显式Legacy类型/adapter/pending解释；尚未完成全部Eulerian/turbulence/IBM runtime迁移。 |

Deferred：剩余模块runtime迁移、任意YAML数学解析、更多target/stage/provider capabilities，以及独立 numerical/MPI regression。
