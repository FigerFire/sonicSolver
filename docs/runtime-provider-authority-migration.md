# Runtime Provider Authority Closure + Active STATE Cleanup

状态：完成：两个原始 case 均完整复跑，前后 evidence exact equality。

本阶段只收口 native single-fluid provider ownership 和 active STATE；WHAT / STATE / HOW / WHICH 顶层定义不变。仓库开始时已有大量未提交修改，本报告的 phase.diff 以本阶段开始时源码 snapshot 为基准，不能用整个 Git diff 代表本阶段变更。

## 1. Baseline before refactor

先复制当前源码和原始输入，构建当前源码，保存 binary，再完整运行以下两个 case；baseline gate 通过后才开始 authority 源码修改。

| case | 原始目录 | 配置与终止时间 | 修改前运行 |
|---|---|---|---|
| Sod | test/Sod/sodCase_weno7_t0p2 | WENO7 / Rusanov / forwardEuler / CFL 0.5 / t=0.2，原始串行 | 454 steps，5 VTS，exit 0 |
| cylinderFlow | test/IBM/cylinderFlowGhost | WENO5 / characteristic LaxFriedrichs / CENTRAL2 / classicalRK4 / CFL 0.3 / t=5，原始串行 Ghost IBM | 11450 steps，51 VTS，exit 0 |

Sod 与已有 frozen baseline 的最终 VTS SHA 完全相同，积分通过已有 5e-12 容差；圆柱与已有 full Ghost baseline 的完整 time/dt 序列、全部 51 VTS SHA、stage closure diagnostics 完全相同。[基线门禁](evidence/runtime-provider-authority-20261002/baseline-gate.json)、[源码 manifest](evidence/runtime-provider-authority-20261002/source-manifest.json)、[binary identity](evidence/runtime-provider-authority-20261002/binary.json)。

before binary SHA256：`c3c3278ec84e125506f416853bfc336e7268d43234306fac4e19380490d718ae`。

after binary SHA256：`0b139c6c8f85387fdccda8bff5473194a3f96417829530438f754ce1fce1d653`。

原始 cylinderFlowGhost 实际没有压力循环、MPI 或 turbulence。已向用户说明这一配置与请求中希望覆盖的 Working(U)/Correction(p) 不吻合，未收到替代配置，因此保留原始 case，不修改算法来制造覆盖。压力 views 仅结构验证，不能以 Ghost CFD 结果证明其数值正确。

原始 full reference 位于 `/private/tmp/sonic-phase29-ghost-full`，日志为 `/private/tmp/sonic-phase29-baseline-20260930/phase29-ghost-full.log`。本阶段完整源码、日志、VTS、两个 binary 位于 `/private/tmp/sonic-provider-authority-20261002`。最初探索运行遇到旧 binary，停止后以 build 完成的 binary 重新正式运行；`.stale.log` 不进入本报告证据。

## 2. Before provider authority

```text
WHICH method selection
  -> compileExecutionProgram / numerical compile / SolvePlanner
  -> Builder::resolveOperationBindings
       -> StateRealization + capability/schedule predicates
       -> ProviderMatchContext
       -> ProviderCatalog matcher chooses backend
  -> runtime callback registration branches on capability flags
  -> PlanExecutor executes OpId
```

准确说，旧 `resolveOperationBindings` 也是 Builder 启动编译时执行，并不是 PlanExecutor 内逐步调用；问题是它在 WHICH method compilation 之后另行选择 backend，并且 runtime 注册阶段还有按 state capability 选实现的分支。

旧实际位置（冻结副本）：[Builder:429](evidence/runtime-provider-authority-20261002/before-source/SF_systemBuilder.cpp)、[resolver:226](evidence/runtime-provider-authority-20261002/before-source/SF_providerResolver.cpp)、[catalog:40](evidence/runtime-provider-authority-20261002/before-source/SF_providerCatalog.cpp)。resolver 的 conservative/constant pressure schedule 判断位于旧 41/61 行，catalog 的 family matchers 位于旧 76–105 行。

## 3. After provider authority

```text
WHAT + active STATE + HOW + WHICH
  -> selected IProvider declares runtimeProvider / operationProvider
  -> compiler validates equation, target, state/view/resource capabilities
  -> CompiledEquationCall.backendProvider
  -> owned numerical / temporal / lifecycle fragment leaves
  -> local SolvePlan lowering + compileOperationBindings validation
  -> returned CompiledSolvePlan: every executable production leaf owns provider
  -> runtime installs implementations for frozen IDs
  -> PlanExecutor validates owner / lookup / invoke
```

Numerical selection occurs once during compilation. Runtime provider lookup is not numerical selection. A compiled operation already owns its selected provider/backend. State realization may validate provider requirements, but it must not route solver families.

`ResolvedSimulationSystem::runtime.operationBindings` 是 compiled leaf ownership 的索引/诊断投影。只有显式 legacyAdapter leaf 可以在 compiler 内调用 legacy selector；返回 compiled system 后同样冻结。

## 4. Removed runtime routing

| 旧行为 | 本阶段结果 |
|---|---|
| ProviderMatchContext 的 conservativeState / phaseState / pressure schedule family match | 删除整个 Context 和 matcher API；catalog 只验证显式 selected ID + operation capability |
| Builder 在 numerical compile 后调用 resolveOperationBindings | 删除；统一进入 NumericalCompiler::compileSystem |
| Native pressure / conservative callback registration 根据 constantDensity/conservativeState 决定实现 | 改为读取冻结 flow.pressure-operators / flow.conservative ownership |
| Native constantPressureScheduleSupported 决定后端 | 改名 pressureOperatorRequirementsSatisfied，仅验证已选 pressure/rhie provider；失败不尝试另一个后端 |
| conservativePressureScheduleSupported | 移入 Legacy，未删除，保留旧 variable/shared pressure 数学路径 |
| phase/Eulerian、IBM compatibility family routing | 显式 Legacy::selectOperationProvider，仅 legacyAdapter leaf 可用，compile-time 冻结；未宣称全部迁移 |

State role summaries 和 signature flags 仍可用于 diagnostics、storage、selected-provider capability validation、composition applicability 和明确 legacy 兼容。它们不在 native PlanExecutor / method ownership selection 中决定 backend。

## 5. Provider selection

| WHICH method/fragment | Primary compiled backend | executable leaf owner |
|---|---|---|
| ConservativeResidual | flow.conservative | fused RHS stage、prepare/dt/begin、physical/time commit 均为 flow.conservative |
| PressureMomentum | flow.pressure-operators | momentum assembly / solve / predictor publication；pressure lifecycle 同 owner |
| PressureCorrection | flow.pressure-operators | pressure correction assembly / solve / correction publication |
| VelocityCorrection | flow.pressure-operators | velocity correction及publication |
| PressureUpdate | flow.pressure-operators | pressure update及publication |
| FluxCorrection | flow.rhie-chow | FluxCorrect 为 flow.rhie-chow；既有 state publication leaf 为 flow.pressure-operators |
| FixedTimeRelaxation | flow.pressure-operators | fixed-time iteration/relaxation/consistency/convergence callbacks 同 owner |
| explicit temporal recipe | 继承 residual body 的 owner | forwardEuler / SSPRK3 / classicalRK4 编译 fragment，recipe coefficient 和 stage 数不变 |

Fusion 要求 backendOperation 和 backendProvider 都一致；保留 continuity/momentum/energy 独立 call/provenance，然后编译为当前 packed ConservativeRHS kernel。已选 backend 的布局能力在编译时验证，不由 runtime 猜测。

## 6. STATE activation

`BuiltinStateCatalog` 的 private metadata vector 知道 14 个标准 symbol：rho、rhoU、rhoE、U、p、T、h、phi、alpha、k、omega、mu、nu、conductivity。构造 catalog 不向 active registry 插入 symbol，不创建 Field/StateBundle 或数值数组，也不启用模型。

`StateRegistry` 仅保存当前 composition 实际激活的 base symbols。`SystemCompositionBuilder::requireState(id)` / contribution.requiredStates 明确请求标准 metadata；`addState(customMetadata)` 保持自定义和显式覆盖 contract。`BuiltinStateCatalog::require(active,id)` 在 active 不存在时才补标准 metadata；未知 id 必须 addState，否则 fail fast。

`requireTargetStates` 激活 HOW 非 Workspace target 的 base symbol；不扫描任意 AST operand 名字来猜模型或添加 equation/HOW。closure/coefficient 依赖由 preset/module 显式请求。

Explain 实际 active：

- Sod：rho、rhoU、rhoE、U、p、T。
- cylinderFlowGhost：上述六项 + mu、conductivity。
- constantDensityPiso（只 explain）：rho 常量、U、p、nu、phi；rhoE/T/k/omega 不因 catalog 已知而激活。

结构测试验证单独 require h 只激活 metadata，不增加 equation/model；inactive k/omega/alpha/rhoE 不自动获得 physical storage。用户普通预制公式仍无需重复 U/p/rho 类型声明。

## 7. STATE views

| view | demand 来源 | storage authority / lifetime |
|---|---|---|
| Physical | HOW target + active STATE contract | 既有 physical storage；rho/rhoU/rhoE alias 同一 packed Q，U/p 为原有 physical authority |
| Working | HOW 的 U* 等 qualifier | compiler lazy manifest；既有 pressure predictor/provider workspace，provider publication写目标 |
| Correction | HOW 的 p' 等 qualifier | compiler lazy manifest；既有 correction workspace；不会改写 Physical(p) |
| OldTime | WHICH temporal method | 多 stage 时才请求原有 old-state workspace，solver execution lifetime |
| Stage | WHICH immutable recipe | stage-qualified manifest；当前 explicit backend stage alias 既有 physical packed authority，没有第二份 Q |
| Workspace | HOW/provider internal numerical target | numerical provider-owned temporary workspace；不预注册成 STATE |

STATE 仍只保存 base symbol，不保存 U*、p' 或 rho_stage_1。`realizeTarget` / `realizeTemporalViews` 的数学与 storage 实现未改。Physical(phi) 保持 canonical face physical contract；临时/修正 flux 属于 backend workspace。

## 8. SystemBuilder

Builder 收集 contribution 的 WHAT/STATE/HOW/WHICH，检查 composition applicability，应用既有显式 transformations，激活 requiredStates 和 HOW base targets，然后调用 compiler。移走 execution ordering、StateRealization compilation、operation provider resolution；Builder 不再根据 realization 给 native 选择 workflow/backend。

Builder 仍保留 EOS/transport/模型 applicability、执行 dependency requirement 汇总、runtime environment 能力检查和 validation/report 组织。这些是已有组合/环境职责，不能表述为 Builder 只剩四条容器 copy。

## 9. Compiler

真实函数链：

```text
SystemBuilder::build
  collect requiredStates -> BuiltinStateCatalog::require
  requireTargetStates(active STATE, HOW)
  NumericalCompiler::compileSystem
    compileStateRealization(state,constraints)        // storage/role summary only
    builtinTemporalMethods().at(selected time).compile
    compileExecutionProgram
      orderExecution(HOW)
      equation registry.at(EquationRef)
      state.at(non-Workspace base target)
      select explicit NumericalBinding (occurrence > equation > wildcard)
      ProviderRegistry::at(method)
      selectedProvider.compile
      backendProvider = selectedProvider.runtimeProvider
      stamp fragment via operationProvider
      realizeTarget(active STATE, compiled target)
      validate workspace/temporal capability
      compileMethodProgram / compatible temporal fusion
      temporalMethod.compileFragment
      realizeTemporalViews
    NumericalCompiler::compile                    // AST spatial recipe binding
    SolvePlanner::compile                         // lower selected fragments
    compileExecutionCapabilities                  // diagnostics/validation
    compileOperationBindings                      // validate frozen owners;
                                                   // legacy-only guarded selection
    reportOperationBindings
  return authoritative ResolvedSimulationSystem
```

HOW target parser 已产生 base symbol + TargetKind，compiler 将 U、U*、p' 分别 realize 为 Physical(U)、Working(U)、Correction(p)。当前真实实现交错 method/view/temporal/spatial compilation，不伪造与源码不符的严格 11 个独立 pass。

## 10. Runtime

`OpRegistry::Entry` 保存 OpId、provider owner、callback。生产 runtime 根据 frozen provider IDs 创建既有资源、注册同 owner callback。`PlanExecutor::validateBindings` 在任何 leaf 执行前检查每个 required OpId 和 owner；executeNode 调用 `operations.invoke(node.operation, execution, node.provider)`。

owner 不匹配立即失败，无 alternative provider search。结构负例验证错误 owner 产生异常且执行计数为零。生产 SystemValidator 要求 resolved leaf 的非空 provider 与 binding 索引一致。

Standalone IR/mock tests 仍允许 anonymous OpRegistry 和空 owner 的手工 plan；这不是可执行 production system 的通过路径。Direct/Linear 等尚无 production backend 的 generic method 不被补隐藏 fallback。Runtime 不读取 state shape 或 PISO/SIMPLE/PIMPLE 来重新定义 HOW/WHICH。

## 11. Remaining legacy

- variable/shared conservative pressure：保留 legacyCouplingPlanFragment/sharedPressurePlanFragment、transformer 和既有 ConservativeRHS 压力执行。
- Eulerian：保留 EulerianStepper、AssemblyPlan、TermKind 和 flow.eulerian-pressure；本阶段只冻结/核对 callback owner，不改 numerical bodies。
- turbulence：保留现有模型/transport provider；catalog 与 active k/omega 分开，不宣称完成 native four-module migration。
- IBM：Ghost/ILW boundary closure、forcing/constraint compatibility 均保留；ibm.constraint legacy ownership 显式标记，未重写 DFM/DLM/KKT。
- StateRealization 的既有 role summary flags 和 legacy policies 保留用于 storage/diagnostic/兼容，不是 native method selector。

Native constant-density PISO/SIMPLE/PIMPLE leaf 的 owner 由 selected method 提供，legacyAdapter=false，不调用 Legacy selector。Explain 在 legacy fragments/leaf 上标示 Legacy；具体数学能力缺失继续 Unsupported，禁止换 solver/scheme 来掩盖。

## 12. Numerical verification

修改前/后输入 SHA 相同；采用相同 build flags、原始 CFL/scheme/tolerance、完整原始终止时间，未新增第三个 CFD case。

| case | steps | endTime | VTS count | after 对 before |
|---|---:|---:|---:|---|
| sod | 454 | 0.2 | 5 | Passed identical（全量 evidence 相等） |
| cylinder | 11450 | 5.0 | 51 | Passed identical（全量 evidence 相等） |

比对内容：完整 printed time/dt 字符串序列、每一份 VTS SHA256、最终 field extrema、正 rho 点的 fluid extrema、质量/动量/总能量积分、rho/p/U/momentum/energy point L2 norms、全部 stage closure 行的 SHA256、全部输入文件 SHA256。相等判断为整个 evidence JSON dictionary equality，不只比较最后一帧。

| case | fluid rho min/max | fluid p min/max | mass | momentum xyz | energy |
|---|---|---|---:|---|---:|
| sod | [0.25924657327767214, 0.7622215092077242] | [13037.114603387012, 68970.12209014251] | 98.67807117165862 | [23619.052634273416, 1.9076916692407833e-08, 0.0] | 21770158.12505416 |
| cylinder | [1.2234377171537614, 1.2289806676854897] | [101290.53594640028, 101354.66376544903] | 236.5775588113999 | [561.0061515566838, -0.03355714201230428, 0.0] | 48915365.8856814 |

积分采用输出 uniform tensor grid 的 trapezoid weights（点值包含 IBM blanking 的零值）；fluid extrema 排除 rho=0 的 blank/solid 点。它们是固定可重复的 VTS 统计，不能误称为 fluid-cut-cell 几何守恒积分。L2 为点值 sqrt(sum(v²))，U 按全部分量统计，具体值见 evidence。

最终 Sod VTS SHA：`8bf55a1c4973f3a2322aeb8785305bf0589cc9eb4cf89505203309f0ccba0768`。

最终 cylinder before/after VTS SHA：`7ce841390cb215ea3b22fdddf910a0879e5eb5335e4f8a27a7948cb9dcec2c5f`。

[Before Sod](evidence/runtime-provider-authority-20261002/before/sod-evidence.json)、[After Sod](evidence/runtime-provider-authority-20261002/after/sod-evidence.json)、[Before cylinder](evidence/runtime-provider-authority-20261002/before/cylinder-evidence.json)。[After cylinder](evidence/runtime-provider-authority-20261002/after/cylinder-evidence.json)。

- Passed identical：已经完整 capture 并 equality 的 case，见上表。
- Failed with first difference：两个完整 case 均未发现 solver 首差异。
- Not tested：pressure Working/Correction CFD、MPI/turbulence CFD、DLM/KKT/forcing 数值；SIMPLE/PIMPLE 的数值收敛行为不在这两个原始配置中。
- Unsupported：缺 provider/capability、错误 packed target layout、unsupported pressure recipe/equation resource contract 均保留 fail fast；没有为了通过本阶段添加替代算法。

原始日志不打印每个 raw candidate/canonical flux 或 raw directional/global residual 数组，未另启内部 dump 造成额外 CFD run。全 VTS、dt 和 stage closure 一致，加上保护 kernel 源码未变，证明此次两条运行未观察到差异；不冒充逐面/raw residual checkpoint 测试。

## 13. Build/static checks

完整 build 通过，architecture checker 通过，git diff --check 通过。9/9 targeted tests：methodObjectCompilation、formulationArchitecture、pisoArchitecture、formulaExecutionContract、programDoesNotOwnRuntimeState、termRecipeAuthority、explicitStageMathematics、fixedTimePressureContract、providerCatalogResolution。

修改前结构测试 8/9；providerCatalogResolution 的一个源项 metadata 断言还期待旧 equationExtensions，当前源码已用 mathematicalExtensions.family=momentum。本阶段 API 迁移时修正断言，不改 source arithmetic。其余失败的探索/工具问题未计为 solver regression。

Architecture checker 对新检查先剥注释，再检查 method owner / compiled leaf / legacy guard / executor owner / catalog activation；保留原有 executor 不含 PISO/SIMPLE/PIMPLE routing 等约束，避免纯注释关键词误报。

| 必须汇报项 | 本阶段结果 |
|---|---|
| 影响层 | composition STATE activation、WHICH contract、compiler binding、runtime registry、explain |
| 数学作用分类 | 改“怎么解 equation”的 authority ownership；未增加数学 term/equation/constraint |
| ownership before/after | 后置 family matcher + runtime capability registration -> compiled leaf provider；catalog混入active -> 分离 |
| execution order before/after | native HOW/loops/stages/commit 顺序不变；provider freeze 移入 compiler |
| public API | requireState + catalog.require；IProvider runtimeProvider/operationProvider；selected-ID catalog resolve；owned OpRegistry bind |
| dependency | Builder 不依赖 resolver；compiler依赖 binding validation；executor只依赖已编译 plan/registry |
| MPI/parallel | COPY/SUM/owner/halo/kernel 源码不变；本次串行未数值验证 MPI |
| numerical semantics | coefficient、flux/residual、dt、pressure/IBM arithmetic 未改；两条 full CFD 对比按上表 |
| regression | 两条 full CFD + 9结构测试 + 静态检查；无新增 CFD case |
| deferred | native variable pressure、turbulence、Eulerian、IBM/user IO 后续迁移，未自动开始 |

自审清单：

- [x] WHAT / STATE / HOW / WHICH 顶层定义没有变化。
- [x] Native provider selection 只在 compile-time 按 selected method 一次完成；legacy-only compatibility 在 compiler 内冻结。
- [x] Runtime resolver 不按 solver-family/state booleans 选择 backend。
- [x] Executable production compiled leaf 包含最终 provider identity。
- [x] StateRealization 负责 role/storage/views，不选择 native algorithm。
- [x] BuiltinStateCatalog 与 active StateRegistry 分离。
- [x] Inactive builtin metadata 不自动分配 physical state 或启用模型。
- [x] U* / p' 不是 StateRegistry base symbols。
- [x] Working / Correction / Stage views 仍 lazy realize。
- [x] Native constant-density PISO/SIMPLE/PIMPLE 不使用 legacy provider routing。
- [x] Conservative fusion backend 由 method compile-time binding 决定。
- [x] PlanExecutor 只校验并执行 compiled decisions。
- [x] Sod/cylinderFlow 完整修改前后 frozen numerical behavior 保持一致。

证据目录：`docs/evidence/runtime-provider-authority-20261002`；该仓库既有规则忽略 docs，文件已保存到工作区但没有强行 git add/commit。Architecture 主文档位于 [src/ARCHITECTURE.md](../src/ARCHITECTURE.md)。
