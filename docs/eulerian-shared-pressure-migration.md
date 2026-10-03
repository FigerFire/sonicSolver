# Eulerian–Eulerian shared-pressure native execution — 2026-10-03

本阶段完成 Eulerian core 的执行权威迁移：逐相方程、共享压力关系、三层耦合循环及 Commit 全部来自原生 WHAT / STATE / HOW / WHICH。`EulerianStepper` 继续拥有工作区、绑定专用数值回调并执行冻结计划。现有数学、BC/halo、HYPRE 参数和相间闭式不变。

## 1. 基线来源与覆盖

使用现有 `test/eulerianEulerianCase`，没有新增 CFD 算例，也没有改动输入参数。原 binary、源码、输入、输出和轨迹在编辑前冻结，原入口直接运行成功，无诊断补丁或启动检查绕过。

该 case 为 water/air 双相、参考相 water、PIMPLE 1 outer / 1 pressure / 1 nonOrthogonal pass，Upwind、HYPRE，终止时间 0.0002。执行 phase continuity、momentum、共享压力、phase/face-flux correction、enthalpy，以及 SchillerNaumann drag、virtual mass、RanzMarshall heat transfer 和 MRF。初始两相温度相同，热交换驱动很小；不能把一次短时回归视为闭式物理精度验证。无 transported turbulence、Lee/RPI、wall heat、IBM；未做 Eulerian MPI/multi-patch 或三相 CFD。

前后均 6 steps / 7 VTS / 168 Ee operations。所有原始输入 hash、完整精度 plan clock/dt 和循环上下文、步诊断、所有输出字段/边界统计及 VTS SHA256 完全相同。`baseline.json` 与 `tools/check_eulerian_native.py` 保留为可重复回归。

## 2. Before 权威链

`SF_presets.cpp::addEulerianEulerianTemplate / addPhaseEquationPack` 注册 STATE、legacy EquationDescriptor/Definition 和 order 10/20/60 `legacyExecution`。`addSharedPressureConstraint` 请求 transformation。

`SF_pressureCoupling.cpp::contributePressureCoupling` 创建 `S_EE_PIMPLE` LegacyExecutionPolicy，另存 outer/pressure/nonOrthogonal counts；`SF_transformation.cpp::LegacySharedPressureTransformer` 生成 opaque E_SHARED_PRESSURE 和整个 Ee operation inventory。

`SF_executionComposition.cpp::composeContributions` 调用 `Legacy::sharedPressurePlanFragment`，片段决定整个 timestep topology。SolvePlanner lower 后，`SF_eulerianStepper.cpp::registerOperations` 绑定数值回调，`advance` 调用 PlanExecutor；`bindSolvePlan` 还依赖 legacy policyKind 判定是否允许运行。

## 3. 数学审计与 native WHAT

审计覆盖 `eulerian/equation/SF_phase.cpp`、`SF_pressure.cpp`、`SF_transport.cpp`、`SF_timeStep.cpp`、`SF_workspace.cpp`，PhaseSystem recovery/reference closure、phase sources/interphase、边界和 EulerianStepper callbacks。

| 数学对象 | 旧描述 / 问题 | 实际冻结内核与 native 表达 |
|---|---|---|
| 非参考相 continuity | ddt(phaseMass)+div(phaseMassFlux)=phaseMassSources | `massNew=massOld+dt*(-div(canonicalMassFlux)+massSource)`。连续更新暂时按质量比缩放 momentum/H，保持 U/h；随后仍使用旧时间 primary workspace 装配矩阵，不重复计入相变传递 |
| 参考相 continuity | 同样标记独立 ddt/div；**pre-existing authority mismatch** | 实际跳过独立 mass update，recovery 计算 `alphaRef=1-sum(alphaOther)`、`massRef=rhoRef*alphaRef`，重新发布 reference momentum/H。native E_CONTINUITY.reference 为此代数闭合关系，并保留 C_VOLUME_FRACTION provenance |
| momentum | `gradient(alphaPressure)` 且 diffusion 在左侧 | 实际为 `alpha*gradient(p)`，并非 `gradient(alpha*p)`；扩散为正贡献。native `ddt(M)+div(momentumFlux)+alpha*gradient(p)=diffusion(alpha*(mu+mu_t),U)+phaseMomentumSources` |
| momentum 数值实现 | diagonal/interphase 细节不明确 | 初始 diagonal=mass/dt；implicit Upwind/diffusion HYPRE 预测，momentumRelaxation，diagonal 除 relaxation。source ledger 中先减去待单独实现的相间 coupling contribution，随后原 semi-implicit all-phase correction 实现 drag/virtual-mass coupling。均保留原实现 |
| shared pressure | `constraint(sharedPressure)=0` 过于 opaque | native `(sumGas(alpha/p)/dt)*pPrime - diffusion(sumPhase(faceAverage(alpha²/A)),pPrime) = sum(massSource/rho)-div(sumVolumeFlux)-pressureHistory`。history=`sumGas(alpha/p)*(p-pOld)/dt`，参考行 pPrime=0；同一 matrix、geometry metrics、HYPRE solve |
| phase correction | 只有 generated operator 概述 | 各相 `M=mass*(predictedU-alpha/A*gradient(pPrime))`，shared `p=p+pressureRelaxation*pPrime`，接着原 primitive/EOS/reference recovery |
| phase face flux | 只有 OP_PHASE_FLUX_CORRECTION | 各相 paired relation 同时表达 volumeFlux 减 response*faceGradient(pPrime)、massFlux 减 faceDensity*同一 deltaFlux；canonical copy/halo 保留原运行语义 |
| enthalpy | diffusion 左侧；未单列 pressure work | native `ddt(H)+div(enthalpyFlux)=diffusion(alpha*(k/Cp+eddyDiffusivity),h)+energyLedger+alpha*((p-pOld)/dt+U·grad(p))`。HYPRE 解 h，发布 H=mass*h 与 T(h)，执行原 recovery |

`SF_eulerianRelations.*` 提供真实 AST；同一 factory 用于注册与 selected provider 的严格 capability matcher。没有引入更理想的新方程、改变通量、压力公式、系数、松弛、floors、相变或壁面公式。既有非法值 fail-fast 保留。

## 4. STATE 与实际 backing

用户选择 Eulerian PhaseSystem/phases 即显式贡献其 STATE，不要求把每相变量重复列进单流体 `state/state.yaml`。StateRegistry 仍是唯一 source authority。

| 对象 | role / storage | 唯一 backing |
|---|---|---|
| phaseMass.phase / momentum.phase / enthalpy.phase | model Transported / ProviderDistributed | PhaseState::primary 的现有 ScalarField/三分量数组 |
| shared p | Algebraic / ProviderDistributed，storage=pressure | PhaseSystem::sharedPressure |
| alpha/rho/U/h/T.phase | Derived / ProviderDistributed | PhaseState::primitive 原 recovery cache；不是 primary HOW target |
| Correction(p) | NumericalProvider view，storage=pressureCorrection | PhaseSolverWorkspace 原压力增量数组，不新增 pPrime base STATE |
| momentumDiagonal、previous primary/velocity/p、canonical phase face flux | numerical Workspace | EulerianStepper::PhaseSolverWorkspace；不注册为用户 base STATE |
| flux output view | NumericalProvider Workspace，3 components / EulerianFace / CanonicalFace | 原 phaseN.volumeFaceFlux；paired mass flux仍使用原数组 |
| time/dt/step | 原 StateBundle clock | 原 timestep authority，Ee callbacks 使用原 dt context，不另算时间循环 |

`registerState` 保留为原 backing 的非拥有注册，primary/derived storage keys 与 source metadata 一致。prepare 将 compiled physical views alias 到这些 DistributedFieldViews，并将 Correction(p)/face workspace 显式 bind 到已有数组。没有新增 phase/p/alpha/U 数组，没有架构态到旧态的 copy bridge。

单元测试使用真实 PhaseState 数组验证 Physical(phaseMass.air) 写入立即出现在原数组；Correction(p) 写入原 delta，物理 p 不改变。CFD 则经过实际 PhaseSystem 初始化和同一 binding。

## 5. After 权威链与 source HOW

`Eulerian PhaseSystem contribution -> native EquationRegistry / StateRegistry / default EquationCalls / NumericalBindings -> explicitly selected PIMPLE composition -> compiler bind/validate -> explicit N-member fragment fusion -> owned CompiledSolvePlan -> PlanExecutor -> frozen flow.eulerian-pressure -> existing assembler / PhaseSystem / HYPRE callbacks`。

HOW 的实际树如下；方程组保留独立成员，source 内没有 Ee OpIds：

```text
Sequence EE.step
  Loop EE.outer [outerCorrectors]
    Sequence phase.continuity
      E_CONTINUITY.water -> phaseMass.water
      E_CONTINUITY.air -> phaseMass.air
    Sequence phase.momentum
      momentum.water -> momentum.water
      momentum.air -> momentum.air
    Loop EE.pressure [pressureCorrectors]
      Loop EE.nonOrthogonal [nonOrthogonalCorrectors+1]
        E_SHARED_PRESSURE -> Correction(p)
        Sequence phase.correction
          correctPhase.water -> momentum.water
          correctPhase.air -> momentum.air
          correctSharedP -> p
        Sequence phase.fluxCorrection
          correctPhaseFlux.water -> Workspace(volumeFaceFlux.water)
          correctPhaseFlux.air -> Workspace(volumeFaceFlux.air)
    Sequence phase.enthalpy
      E_ENTHALPY.water -> enthalpy.water
      E_ENTHALPY.air -> enthalpy.air
  Commit EE.commit
```

三层次数只保存在 source HOW，不另存 execution policy。`SF_eulerianCoupling.cpp::applyEulerianExecution` 消费明确注册的 PIMPLE；当前其他 coupling 被明确拒绝，不从 Eulerian 或 STATE 自动选择 PIMPLE。WHICH capability matcher 检查该 backend 能实现的 authored topology，不更改顺序或循环次数。

## 6. WHICH、group 与 compiled micro-order

| 方法 | 独立 equation -> target 成员 | fusionKey / runtime owner | fragment |
|---|---|---|---|
| EulerianPhaseContinuity | 每相 E_CONTINUITY -> phaseMass | 同名 explicit group / flow.eulerian-pressure | interphase/source prerequisites、diagonal/flux、continuity、BC |
| EulerianPhaseMomentum | 每相 momentum -> momentum | 同名 explicit group / 同上 | all-phase predictors/interphase correction、BC/diagonal/flux publication |
| EulerianSharedPressureCorrection | E_SHARED_PRESSURE -> Correction(p) | 单 occurrence / 同上 | solve、publish、sync |
| EulerianPhaseCorrection | 每相 correctPhase -> momentum，加 correctSharedP -> p | 同名 explicit group / 同上 | correctAllPhases 一次 |
| EulerianPhaseFluxCorrection | 每相 correctPhaseFlux -> face workspace | 同名 explicit group / 同上 | paired flux correction、BC、canonical copy |
| EulerianPhaseEnthalpy | 每相 E_ENTHALPY -> enthalpy | 同名 explicit group / 同上 | all-phase energy、final BC、outer validation |

实际 compiled order 完全保留：

```text
ee.dt.compute -> ee.step.begin
Loop outer:
  ee.interphase.compute -> ee.sources.assemble
  [legacy turbulence.prepare if contributed]
  ee.sources.validate -> ee.momentum.diagonal -> ee.momentum.flux
  ee.faceFlux.canonical -> ee.continuity.assemble -> ee.boundary.prepare
  ee.momentum.solve -> ee.interphase.correct -> ee.boundary.afterMomentum
  ee.diagonal.sync -> ee.momentum.flux.after -> ee.faceFlux.canonical.after
  Loop pressure:
    Loop nonOrthogonal:
      ee.pressure.solve -> ee.pressure.publish -> ee.pressure.sync
      ee.phase.correct -> ee.faceFlux.correct -> ee.boundary.afterPressure
      ee.faceFlux.canonical.pressure
  ee.energy.solve -> [legacy turbulence.solve if implemented]
  ee.boundary.final -> ee.outer.validate
ee.step.commit -> ee.time.commit
```

包括 phase/flux correction 在内的所有七个压力操作，都处在原 nonOrthogonal pass 内。prepare/dt 来自 continuity method 的 compiled root lifecycle，state/time commit 来自 source Commit 的 method lowering。没有 fake source/interphase physics equation，也没有 source BC/MPI nodes。

## 7. Compiler、runtime 与旧组件去向

通用 fusion 由原 leaf 扩展到 N-member fragments：验证完整 explicit members、targets、相邻 Sequence siblings、无重复/跨 scope、同 method/backend、相同 executable fragment topology，包括 loop/termination/legacy 属性。lower 第一 member 的 fragment 并把全部 equation/target provenance 附在各 leaf；保留每一 source occurrence 和 STATE view。

`registerOperations` 的 callback 分类审计：InterphaseCompute/SourcesAssemble、MomentumDiagonal/MomentumFlux/MomentumFluxAfter、ContinuityAssemble/MomentumSolve/InterphaseCorrect、PressureSolve/PhaseCorrect/FaceFluxCorrect、EnergySolve 为数值内核或源项 prerequisite；全部 Boundary、DiagonalSync、PressureSync 和 FaceFluxCanonical variants 为 boundary/communication lifecycle；PressurePublish 为增量发布并累计线性诊断；SourcesValidate/OuterValidate 为检查/诊断；DtCompute/StepBegin/StepCommit/TimeCommit 为原 dt、准备、旧时间工作区和 StateBundle clock lifecycle。保留的 TurbulencePrepare/Solve 是未迁移兼容回调。它们没有选择 outer/pressure/nonOrthogonal 次数或重排方程；没有 legacy scheduler decision callback。

Provider 验证 inputs 与 PhaseSystem 声明 storage slot 顺序、每相 backing/layout、数学 AST、shared p 和 Correction(p)。排序/匹配不依赖 unordered map 或无关 STATE insertion。通用 compiler 没有 Eulerian/phase-name 路由分支。当前 recipe/host backend capability guard 保留；没有扩展 MPI、多 patch 或其他 phase convection。

| 旧组件 | 处理 |
|---|---|
| Legacy::sharedPressurePlanFragment | 删除声明、定义与生产调用；没有兼容 scheduler 副本 |
| kSharedPressureScheduleId / S_EE_PIMPLE policy | 删除；contributePressureCoupling 不再生成 schedule/count policy |
| phase legacyExecution 10/20/60 | 替换为 native EquationCalls 和 explicit numerical bindings，无重复 occurrence |
| LegacySharedPressureTransformer runtime inventory | transformer 改为数学关系与 constraint transformation；不再生成 Ee operation inventory |
| EulerianStepper policyKind/pimpleStageActive | 删除旧策略检查；bindSolvePlan 检查冻结计划 identity 与实际 callback capability |
| registerOperations / advance | 保留 callback binding/PlanExecutor；没有自行实现 outer/pressure/nonOrthogonal loop |
| Legacy::selectOperationProvider 的 Eulerian 分支 | 删除；core leaf 由 method 冻结 owner，legacyAdapter=false |
| AssemblyPlan / TermKind / legacy Definition | 保留为专用 assembler 的临时 term inventory adapter，无 order/count/provider authority |
| Eulerian turbulence | 保留 legacy entries及必要的 provider-local compatibility leaves；未迁移为 native turbulence HOW |
| IBM adapters | 本阶段未迁移 |

## 8. AssemblyPlan 边界与 Eulerian turbulence

最终 authored AST freeze 后，`projectEulerianBackendDefinitions` 从它重新投影 legacy term inventory，包括模型 source extensions。旧 Definition 不再是第二个可编辑数学输入。符号运算、乘积与正负号由 AST 保留，旧 assembler 仅消费原 term capability checks 和自身冻结实现；本阶段没有声称它已直接执行任意 AST。

参考相 AST 没有独立 ddt/div，旧 assembler 对参考相的两处 bind-time term 要求和 assembleContinuity 的 term check 改为只检查非参考相。实际 update、缩放、recovery、矩阵和源项没有修改。

Eulerian RAS 原有断点也如实保留：临时 variant 明确配置 phases 后，冻结旧入口报 `Cannot extend undeclared execution policy S_EE_PIMPLE`；旧 contribution 的 `k.phase/omega.phase` identity 又与 runtime 的 `E_TURB_*` 检查不一致。删除 obsolete policy extension 后，目前明确报告这些 legacy equations 没有 compiled HOW ownership；未静默关闭湍流，也未声称该路径 Runnable。旧真正 E_TURB_* backend adapter 所需 prepare/solve leaves仍可保留于 provider-local compatibility。实际 Eulerian turbulence CFD **Not tested / native execution deferred**。

下一阶段应是 **Eulerian Equation/Assembly Authority Migration**：让专用 assembler 从 authoritative AST/provider term contracts 装配，删除 AssemblyPlan/TermKind/legacy Definition adapter；本轮不自动启动。随后独立迁移 Eulerian turbulence，再处理 IBM 和一般自定义四模块 IO。

## 9. 数值结果、验证与限制

| 比较项 | before / after | 状态 |
|---|---|---|
| inputs / mesh / numerical parameters | 全部相同，无 tuning | Passed identical |
| steps/time | 6 / 0.0002 | Passed identical |
| plan context / clock / dt | 168 完整精度轨迹相同 | Passed identical |
| outer / pressure / nonOrth passes | 每步 1 / 1 / 1 | Passed identical |
| all phase fields / shared p / boundary stats | min/max/L2/sum 全部相同 | Passed identical |
| VTS | 7/7 SHA256 相同 | Passed identical |
| closure / total mass / H / pressure linear diagnostics | 6 条完整 step diagnostics 相同 | Passed identical |
| interphase/source/output diagnostics | VTS 中 source与mechanicalHeating字段相同；未打印的单独 ledger audit 不虚构 | Passed identical / unavailable diagnostics explicitly noted |
| Eulerian MPI/multi-patch、三相 CFD、多 pass CFD、phase change/RPI | 未运行 | Not tested |

最终 shared `p`（不是旧 Field placeholder `Pressure`）为 min 101324.99002629421、max 101325.01993900056、L2 1468338.6868276538。最终 total phase mass=235.8562、total phase enthalpy=295875700（原打印精度），maxAlphaSumError=0，pressure solve iterations=2，relative residual=3.413064e-16，HYPRE rebuilds/solves=1/6。

最后一份 VTS SHA256：`65743c2113d900e1c09b35b9601a1026a8069a5e3097d3ddc7f88faa00256519`。所有原输出数组及统计见 baseline/evidence；原 VTS 未直接输出 primary momentum/h，因此额外 momentum/h statistics 明确由 mass*U、H/m 重建，不冒充 raw primary dump。

完整 configured build 成功。完整 CTest 最终 **33/33 Passed**，包括旧 single-fluid pressure/Sod/turbulence、MPI pressure/halo 契约，以及新 nativeEulerianExecution / eulerianNativeNumericalRegression。曾有一次新增兼容测试错误要求 Eulerian RAS Runnable，现已用冻结旧入口证据纠正为明确 Unsupported；失败日志保留，不删除或伪称通过。

新增结构测试覆盖三相 separate WHAT/HOW/STATE，explicit group provenance，2 outer * 3 pressure * 2 passes 的 12 次 pressure/phase/flux operations，其他 all-phase kernels 每 outer 一次，step/time commit 一次；缺项、重复、wrong target、bad backing、scope splitting、math mutation、wrong Correction(p)、missing Commit、unsupported reordering和PISO拒绝。它们属于结构/执行次数验证，不代替三相或多 pass CFD。

`test_methodObjects.cpp` 还以无 Eulerian 名称的 synthetic provider 独立验证通用两成员、三成员 fragment fusion：逆序 registry insertion 不改变显式成员顺序；缺成员、重复、错误 target、重排、跨 scope 和不同 provider 都被拒绝。所有 source calls/targets 与每个 compiled leaf 的成员 provenance 保留。

architecture checker、consumer scope/12-field aggregate guard、git diff --check 均通过。170 个受保护数值源码文件与冻结源码逐字节相同；12 个关键函数中 11 个函数体完全相同，continuity 函数去掉参考相 term-capability 条件这一处差异后完全相同。protected source/function equality 与完整 before/after source/binary/inputs/output/logs/manifest 归档在 `docs/evidence/eulerian-native-20261003`。MPI owner/COPY/SUM/canonical face 规则无变化；本回归只证明串行 Eulerian 冻结结果。

## 10. AGENTS 非平凡修改报告

1. **层**：模型数学/STATE composition、native HOW、WHICH provider lowering、通用 fusion/compiled workspace layout、runtime binding、解释/静态守卫/回归。
2. **分类**：改怎么解已声明 equations；补充冻结算法实际使用的 correction relations。未增加新的物理模型或改变 PDE arithmetic。
3. **ownership**：旧 shared-pressure policy/fragment 的 schedule 转交 ExecutionProgram；PhaseSystem arrays、StateBundle clock、Stepper workspace 不变。
4. **order**：source authority 改变；compiled Ee 顺序与操作次数完全相同。
5. **public API**：去掉 shared-pressure fragment/schedule constant；增加 Eulerian relation/method/preset bindings；CompiledTarget 的 workspace shape/location/ownership 是 compiler output，source Target 不含 storage。
6. **dependency**：专用 provider 消费 core math/state/HOW IR；runtime 消费冻结 compile products，未引入对 aggregate 的依赖、reverse model->solver lifecycle 或第二个 SimulationPlan。
7. **parallel**：无 raw MPI、新平均、新 COPY/SUM 更改或 HYPRE partition 修改。
8. **numerical**：无公式/stencil/松弛/容差/CFL更改；参考相 term capability guard 纠正，recovery arithmetic不变。
9. **regression**：33/33、7 VTS/168operations identical，结构/alias检查和静态守卫通过。
10. **deferred**：Assembly AST直接装配、Eulerian turbulence启动/identity断点、Eulerian MPI/三相/多 pass CFD、IBM和一般用户四模块 IO。

## 11. 19 项自审

- [x] phase continuity/momentum/enthalpy occurrences 为 native HOW。
- [x] outer/pressure/nonOrthogonal loops 为 native HOW。
- [x] sharedPressurePlanFragment 不再是生产 authority。
- [x] shared-pressure policy 不再重复存 loop counts。
- [x] shared-pressure transformer 不声明完整 Ee lifecycle。
- [x] 每相物理变量一个 backing。
- [x] pressure correction 为 view，不是另一个 base p。
- [x] source HOW 无 Ee OpIds。
- [x] BC/halo/canonical copy 留在 provider/runtime lifecycle。
- [x] interphase/source prerequisites 没有伪装为物理 equations。
- [x] all-phase kernels 每个实际 group 只执行一次。
- [x] phase provenance 保留。
- [x] provider owners 在 compile 时冻结。
- [x] EulerianStepper 没有独立 scheduler。
- [x] 数值算术保持冻结。
- [x] operation order 与原 baseline 相同。
- [x] generic compiler 无 Eulerian-name routing。
- [x] Eulerian turbulence 明确 deferred。
- [x] AssemblyPlan/TermKind adapter 状态与下一阶段边界明确。
