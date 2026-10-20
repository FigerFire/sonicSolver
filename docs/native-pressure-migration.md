# Native Pressure Validation + Conservative Pressure Migration

2026-10-02。本阶段完成既有 native 常密度压力 CFD 验证，以及 single-fluid conservative / variable-density pressure 的 native execution authority 迁移。四模块定义保持原样，未开始 turbulence、Eulerian、IBM 或 user IO 后续迁移。

证据目录：[evidence](evidence/native-pressure-migration-20261002/)。完整源代码冻结、二进制、网格和各回归生成的 case 保留在 `/private/tmp/sonic-native-pressure-20261002`；长期证据目录保存 source/input/output hashes、全部 run/explain logs、字段统计、基线门禁、前后比较及保守压力的实际 YAML/VTS。现有工作树中的先前修改保持不动。本阶段修改清单：[phase-changed-files.json](evidence/native-pressure-migration-20261002/phase-changed-files.json)。

## 1. Native pressure validation

先完整运行当前实现，再开始 Phase B 源码修改。Phase A 门禁 [native-baseline-gate.json](evidence/native-pressure-migration-20261002/native-baseline-gate.json) 全部通过。使用冻结的既有 regression scripts，保留原输入、既有 test variants 与 tolerances；没有新增 CFD case。`SF_PLAN_TRACE=1` 只记录原操作执行，不修改数值配置。

| 既有 suite | 每个 before/after 的 solver run 数 | 实際覆盖 |
|---|---:|---|
| `check_constant_density_pressure.py` / constantDensityPiso | 2 | Working(U) predictor、Correction(p)、Physical(U/p/phi)、1/2 次 pressure corrector、gauge、checkerboard/continuity |
| `check_poiseuille_pressure.py` / constantDensityPoiseuille | 4 | 16²/32² 通道速度剖面、压降、连续性、边界/reference |
| `check_fixed_time_pressure.py` | 16 | PISO/PIMPLE outer=1 等价；SIMPLE/PIMPLE fixed-time base、outer min/max、early exit、relaxation、flux restore、iteration begin/end、correction/state/time commit、零 gauge |
| `check_distributed_pressure.py` | 32 | cavity/channel 的 211/121/221 分区；SIMPLE/PIMPLE 211；平衡态 early exit；pressure matrix/RHS、U/p halo、canonical phi/Rhie–Chow/owner COPY |

Working(U)、Correction(p)、Physical(phi)、pressure Loop、SIMPLE/PIMPLE outer iteration 和 pressure MPI 均实际进行了数值运行，不只是 explain。原 PressureOperators::bind 将原工作数组绑定为 views，assembleMomentum 通过 workingVelocityView_ 推进，correctVelocity/publishVelocity 发布 Physical(U)，压力求解/修正和面通量发布使用原有 provider backing。

既有 cases 的 nonOrthogonalCorrectors=0，因此数值验证的是 generic nonOrthogonal Loop 的 **一个 pass**。更多 non-orthogonal passes 没有 CFD coverage，仍由现有能力校验拒绝；没有更改输入来制造覆盖。

## 2. Before legacy pressure architecture

冻结源代码中的实际调用链：

```text
SystemBuilder::build
 -> Compose::addPressureConstraintFluid (SF_presets.cpp)
      legacy physical DSL; packed momentum 错用 U 名称
 -> PressureConstraintTransformer (SF_transformation.cpp)
      base pPrime + flat pSimple
 -> contributePressureCoupling
      S_PRESSURE 的第二份 outer/pressure/nonOrthogonal schedule
 -> composeContributions (SF_executionComposition.cpp)
      Legacy::lowerPressure + legacyCouplingPlanFragment
 -> SolvePlanner / compileOperationBindings
      legacyAdapter leaves -> Legacy::selectOperationProvider
      conservativePressureScheduleSupported 再检查 schedule
 -> SingleFluidStepper::bindSolvePlan / bindPressureProvider
 -> bindPressureOps -> OpRegistry -> PlanExecutor
```

`Legacy::lowerPressure` 声明操作/compiled equations，legacyCouplingPlanFragment 拥有旧顺序，Legacy selector 从状态与 schedule 决定 flow.conservative。`bindPressureOps` 中有效的 arithmetic/callback registration/state publication 必须保留。迁移前 explain：[before/conservative/explain.log](evidence/native-pressure-migration-20261002/before/conservative/explain.log)。

## 3. After migration architecture

```text
WHAT: EquationRegistry
STATE: active StateRegistry
HOW: native ExecutionProgram
WHICH: NumericalSelection / explicit ConservativePressure* bindings
  ↓ compileExecutionProgram / NumericalCompiler / SolvePlanner
  ↓ compileOperationBindings 验证 selected provider capability
native CompiledSolvePlan (所有 12 production leaves owner=flow.conservative)
  ↓ OpRegistry owner lookup/validation + generic PlanExecutor
原 ConservativeRHS + Time::Explicit::forwardEuler + PressureBased::Corrector
```

没有平行 SimulationPlan，没有新 solver family，没有 runtime state/schedule provider discovery。架构统一并不要求数值 kernel 统一：常密度仍使用 flow.pressure-operators / flow.rhie-chow，保守变量仍使用 flow.conservative。迁移后 explain：[after/conservative/explain.log](evidence/native-pressure-migration-20261002/after/conservative/explain.log)。

## 4. Equation mapping

物理方程复用 `Preset::installSingleFluid` 的原 identity：continuity（ddt(rho)）、momentum（ddt(rhoU)）、energy（ddt(rhoE)）。没有新增 predictorMomentum 等副本。旧的 packed U 实际是 rhoU，现使用正确数学名字，primitive U 为 derived view。

现有 pressure 数学以 `conservativePressureRelations()` 注册，并由 selected provider 对同一 factory 的 canonical AST 做 capability 验证：

| Equation | 原 numerical implementation 的数学意义 |
|---|---|
| pSimple | p′/(rho c² dt²) − div(faceAverage(1/rho) grad(p′)) = −div(U)/dt；原 gauge、BC 和矩阵结构不变 |
| correctU | rhoU ← rhoU − momentumRelaxation dt grad(p′) |
| correctP | prepared p ← physical p + pressureRelaxation p′ |
| correctFluxp | invalidate derived flux after conservative momentum correction |
| publishPressure | rhoE ← EOS.totalEnergyFromPressure(rho, rhoU, prepared p) |
| relaxIterate / restoreFlux / checkConvergence | generic fixed-time relation 描述；conservative provider 缺少该 execution capability，不能 Runnable |

修改 pressure relation AST 而没有实现对应 provider 会 fail fast。新增单元验证覆盖这一点；保守融合 flux operand 换成未知数学也不能静默调用原 flux kernel。continuity/energy 是 predictor 的显式 source-math fusion inputs，不会额外推进一遍。

## 5. STATE mapping

| Base symbol | 当前 storage / derivation |
|---|---|
| rho | packed conservative，offset 0，1 component |
| rhoU | packed conservative，offset 1，3 components |
| rhoE | packed conservative，offset 4，1 component |
| U | 从同一 Q/EOS 读取 velocity |
| p | 从同一 Q/EOS 读取 pressure |
| T | 从同一 Q/EOS 读取 temperature |

实际 compiled demands：Physical(rhoU)、Correction(p)、Working(p)、Workspace(fluxValidity)、Physical(rhoE)。rho 的 packed backing 同时被显式 fused math contract 读写，runtime base binding 使用原 physical authority。

Correction(p) 非 owning alias 到 Corrector::Workspace::correction；Working(p) alias 到 targetPressure；Workspace(fluxValidity) alias 到原 fluxCorrected flag。view backing 只在原 correction lifetime 内可读，越界/lifetime 错误 fail fast。StateRealizer 不分配第二份压力/Q。pPrime、U*、p′ 均不是 active base symbol；fluxValidity 是 numerical workspace，未注册成 physical STATE。没有激活 phi/alpha/k/omega 等未参与此 case 的 catalog symbol。

原预测器本来直接发布 Physical(Q)，所以 HOW 真实表达 Physical(rhoU) 及 rho/rhoU/rhoE fusion。没有为了理想图造出 Working(Q) 并改变 publication timing。该 conservative Corrector 没有独立 persistent phi，correctFlux 仅标记 derived flux invalidated；下一次 spatial assembly 从已修正 state 构造通量。

## 6. Actual HOW

pressureConstraintPiso 的源 HOW（与最终 explain 相同）：

```text
Sequence
  10 momentum@predictor -> Physical(rhoU)
  20 Loop pressure limit=1
    Sequence pressureCorrection
      10 Loop nonOrthogonal limit=1
        10 pSimple -> Correction(p)
      20 correctU -> Physical(rhoU)
      30 correctP -> Working(p)
      40 correctFluxp -> Workspace(fluxValidity)
      50 publishPressure -> Physical(rhoE)
  1000000 Commit physicalStep.commit
```

lower 后原 12 操作顺序完全保持：

```text
pressure.prepare (BC/halo + 原 CFL dt)
momentum.assemble (ConservativeRHS)
momentum.solve (forwardEuler + publish Q)
[pressure.boundary.prepare -> pressure.assemble -> pressure.solve
 -> velocity.correct -> pressure.update.prepare -> flux.correct
 -> pressure.correction.commit (EOS energy publication)]
pressure.step.commit (final BC/closure/diagnostics)
time.commit (一次 dt/time/step commit)
```

HOW Loop 是 pressure count 的唯一源 authority；改变 pressureCorrectors=2 的结构测试检查 native Loop/capability counts。Commit 不另起 lifecycle。WHICH 不增加源 calls/loops、不修改源 order；provider 的 local assemble/solve 和 scope lifecycle 由当前编译接口 lower。

SIMPLE generic HOW 为 outer Loop → predictor → 单 correction sequence → relaxIterate → restoreFlux → checkConvergence → Commit；PIMPLE 在 outer 内再嵌 pressure Loop。原常密度 implementation 支持并数值验证此流程。conservative 同一 HOW 可表达这些节点，但 FixedTimeIteration capability 未实现，所以明确 Unsupported。没有用 PISO 代替 SIMPLE/PIMPLE。

## 7. WHICH and frozen owners

| HOW occurrence | Explicit method | Compiled provider / leaf owner |
|---|---|---|
| momentum@predictor → Physical(rhoU) | ConservativePressureMomentum，fusion continuity+energy | flow.conservative |
| pSimple → Correction(p) | ConservativePressureCorrection | flow.conservative |
| correctU → Physical(rhoU) | ConservativeVelocityCorrection | flow.conservative |
| correctP → Working(p) | ConservativePressureUpdate | flow.conservative |
| correctFluxp → Workspace(fluxValidity) | ConservativeFluxCorrection | flow.conservative |
| publishPressure → Physical(rhoE) | ConservativePressurePublication | flow.conservative |

root/Commit leaves 也由该 selected predictor method freeze 为 flow.conservative。最终 explain 中 12 leaves 均显示 `legacyAdapter=false`，单元测试逐 leaf 验证。momentum occurrence 显式声明 writes=rho/rhoU/rhoE、三个 math references，explain 输出完整 writes，不把 side effects 隐藏为只有 rhoU。

conservative 的 fixed-time selected methods 为 ConservativeFixedTimeRelaxation、ConservativeFluxConsistency、ConservativeResidualConvergence；它们不会获得不存在的 backend capability。time recipe 仍为 immutable forwardEuler；unsupported stage treatment 不被自动降为一个 stage。

## 8. Legacy deletion / isolation

**删除：** Legacy::lowerPressure 声明/定义、legacyCouplingPlanFragment 声明/定义、conservativePressureScheduleSupported、couplingReportOperations（操作诊断改读最终 owned plan）。单流体 S_PRESSURE policy 和 pPrime base registration 不再生成；旧 conservative pressure compiledEquations 不再生成。

**不再可达：** migrated native leaves 不调用 Legacy::selectOperationProvider。Legacy selector 仍为未迁移 compatibility operations 存在，禁止用它匹配 pressure operation 并恢复 conservative schedule authority。

**仍需保留：** `Legacy::sharedPressurePlanFragment` 明确隔离在 Legacy namespace，只用于 C_SHARED_PRESSURE / Eulerian shared pressure policy。composeContributions 若缺少单流体 native predictor，明确报错，不提供 generic legacy pressure fallback。Eulerian 操作顺序和三层 Loop contract 的既有结构检查继续通过。没有把 relocation 称作删除。

## 9. bindPressureOps audit

函数保留，定位为 **callback registration only**。其完整函数体与冻结源完全相同。它保留 CFL dt、ConservativeRHS、Euler predictor、pressure matrix/solve、velocity correction、pressure preparation、derived-flux invalidation、EOS energy commit、BC/halo 与 clock callbacks。

它没有 timestep/pressure-corrector/outer loop，不拥有次数或全局顺序。CompiledSolvePlan 决定每个已注册 callback 在何时被调用。bindPressureProvider 只按冻结 operation binding 构造对应实现，不重新根据 state/formulation/coupling 选择 provider。其 `LegacyPressureInputs` 数据名仍作为原 numerical constructor inputs 保留，不代表 legacy scheduler；user-IO compatibility 命名不在本阶段迁移。

## 10. Conservative numerical before/after

先通过迁移基线门禁 [migration-baseline-gate.json](evidence/native-pressure-migration-20261002/migration-baseline-gate.json)，然后修改源码。完整 15 个 case input 文件的前后 SHA256 也完全相同，见 [conservative-all-inputs.json](evidence/native-pressure-migration-20261002/conservative-all-inputs.json)。原 scheme/CFL/time/linear/relaxation/corrector 设置保持不变，没有调参。

**结论：Passed identical。** [comparison.json](evidence/native-pressure-migration-20261002/comparison.json) 对输入、两个 VTS hashes、所有输出数组 min/max/L2、time/dt、math diagnostics 和 12 操作顺序做 exact comparison；全部 true。实际 VTS 也保存在 evidence 的 before/after/pressureConstraintPiso/result。

| 量 | before = after |
|---|---|
| Steps / pressure correctors | 1 / 1 |
| Final time / dt | 1.659315e-06 / 1.659315e-06 |
| rho extrema | 0.9999955959523205 … 1.0000037674640512 |
| U component extrema | −3.135893340993762e-05 … 0.050263138503177415 |
| p extrema | 101323.94804752528 … 101326.05196077126 |
| rhoE extrema | 253309.87011881324 … 253315.13112216647 |
| pressure linear solve | iterations=2，residual=4.3877e-11 |
| continuity maxDiv before/after | 2.01178 / 2.01052 |
| interface/jump | interfaceFaces=0，maxJump(target/correction)=0/0 |
| HYPRE rebuilds / solves | 1 / 1 |

字段 extrema/norm 是原 VTS PointData 的完整数组统计，不声称是 unique-owned-cell 积分。Euler final 和 flow correction 的 closure diagnostics、soundSpeed/temperature/rhoE 全数组输出一致；原 EOS.totalEnergyFromPressure commit 保持不变。

```text
initial VTS SHA256
4da441517514509c2f477076bef2685d1a562e5b666ee92d70ff6ce14b9c37d2
final VTS SHA256
4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511
```

既有 checker 的 frozen tolerance 5e-12 没有放宽，并且最终 VTS byte-identical。原算法的 maxDiv 只小幅下降是冻结行为，本阶段没有数值修复。

另对 12 项 frozen kernel 文件/函数体做 exact source comparison：ConservativeRHS、Time::Explicit、compressible fused backend、PressureOperators、Rhie–Chow 整个文件，以及 Corrector 的 assemble/solve/preparePressureUpdate/correctVelocity/correctFlux/commitPressureUpdate 与 bindPressureOps 函数体。全部 unchanged。

## 11. Constant-density regression

Phase A 和最终 Phase B 二进制都跑完同一四个 suites，共 54 solver runs/phase；另有 conservative 1 run/phase。所有既有 assertions/tolerances 通过。同配置 before/after 输入、VTS、字段统计、time/dt 和数值诊断全部 exact identical。

PISO 首两次 correction 的 maxDiv：(2.0, 3.39853e-08)，(0.64824, 4.63737e-08)；checkerboard/rms=4.594e-05，reference=101325。

Poiseuille 16²：L2=0.000746526192，Linf=0.000988599105，maxDiv=2.23908e-09；32²：L2=0.000354528614，Linf=0.00047749875，maxDiv=1.19914e-09；压降均为 0.1。

fixed-time cavity：clock 为 (0.000390625,0.000390625)、(0.00078125,0.000390625)；PIMPLE outer=1 与 PISO 一致；SIMPLE/PIMPLE max=4 实际每步于 2 次 outer 退出，relaxed flux continuity 与 zero gauge 检查通过。SIMPLE Poiseuille L2=0.00822508807、Linf=0.0112953895、maxDiv=1.37627e-09、压降=0.1，每步 outer=2，finalMaxDelta=9.26491e-06、finalConverged=false（保留原结果）。

MPI 同配置 before/after 的 VTS hashes 均 exact identical，但各 rank console trace 的交错顺序不同，不称为日志逐字节相同。serial vs MPI 遵循原 tolerance：matrix/RHS 1e-12、U 1e-8、pressure shape 1e-6，不能称 serial/MPI bitwise identical。原 socket sandbox 限制导致一次 baseline MPI 启动失败，随后以同一命令、相同 inputs 在允许本机通信的执行环境运行通过；失败尝试单独保留，没有改变 MPI kernel/ownership。

## 12. Remaining legacy / deferred capabilities

剩余 legacy 是 Eulerian/shared pressure、turbulence specialized execution、IBM forcing/DLM/KKT adapters、unmigrated Eulerian AssemblyPlan/TermKind，以及 explicit compatibility / user-defined IO 待迁移部分。explain 中 classification 的 legacy-configured 属于现有输入/诊断标签，不参与 provider selection；不是本路径 legacyAdapter。

conservative single-fluid pressure 的 production PISO 已 native；未实现的 conservative SIMPLE/PIMPLE fixed-time、额外非正交 pass、多 temporal stages、多 patch/MPI conservative pressure 与额外 IBM/phase/turbulence services 继续 fail fast，不能解读为新增支持。这次 pressure MPI 数值证据属于既有 constant-density distributed pressure，不扩大为 conservative pressure MPI。

后续仅报告顺序：turbulence → Eulerian/shared pressure → 剩余 IBM adapters → user-defined 四模块 IO。未自动启动。Sod/Ghost 数值证据来自上一阶段，这一阶段没有重跑或扩大其覆盖范围。

## 13. Build / static checks and change impact

完整 configured build 通过；architecture checker 通过（保留原两个 allowlisted methods→solver dependency debt，未新增）；git diff --check 通过。最终状态：[verification.json](evidence/native-pressure-migration-20261002/verification.json)。这些结论与上面的 CFD numerical validation 分开。

12/12 targeted tests：pisoArchitecture、formulaExecutionContract、rhieChowFaceCoupling、fixedTimePressureContract、explicitStageMathematics、methodObjectCompilation、providerCatalogResolution、termRecipeAuthority、programDoesNotOwnRuntimeState、formulationArchitecture、distributedPressurePrimitives、distributedPressureFailStop。新增/调整 contract 检查保护 native leaf ownership、correct conservative target/fusion、lazy pressure views、changed relation rejection、unknown fused operand rejection、HOW counts 和既有 Eulerian compatibility；没有删除数值断言来掩盖失败。

| AGENTS.md change dimension | 本阶段 before → after |
|---|---|
| Layer | composition + method/provider contract + state view binding + capability validation + explain/docs |
| Mathematical class | 改 execution authority / 怎么执行已有 equation；把已有 pressure relations 显式注册，未增加物理 equation 或改变数值 term |
| Ownership | legacy schedule/selector → native HOW + selected/frozen provider；Physical Q 与 Corrector workspace authority 保持 |
| Execution order | 原 12 OpIds → 相同 12 OpIds；time/BC/halo/publication 位置不变 |
| Public API | 删除 lowerPressure/legacyCouplingPlanFragment/couplingReportOperations；shared fragment 明确 Legacy namespace；pressure HOW/WHICH 接受 explicit primary target；Corrector 增加 bindStateViews；增加 selected ConservativePressure* method IDs |
| Dependency | numerical providers 使用现有 equation/coupling math factory；Corrector 通过现有 StateRealizer 接口绑定 non-owning views；generic compiler/runtime 不新增 domain/model dependency |
| Parallel | raw MPI、owner COPY、canonical face、GlobalDof SUM、HYPRE partition/backend 均未修改 |
| Numerical | frozen numerical files/function bodies unchanged；dt/stage/flux/residual/EOS/relaxation/tolerance 不变 |
| Regression | 55 runs/phase 的同输入压力回归 + 12 units + static checks |
| Deferred | 仅上述 capability gaps 与后续模块；没有 hidden fallback 或自动下一阶段 |

## Self-review — 19 items

- [x] Native pressure CFD 实际使用 Working(U)。
- [x] Native pressure CFD 实际使用 Correction(p)。
- [x] Physical(phi) publication 实际进行数值验证。
- [x] SIMPLE/PIMPLE fixed-time semantics 已测试。
- [x] 已运行既有 pressure MPI regression。
- [x] single-fluid conservative pressure 不再依赖 Legacy::lowerPressure。
- [x] single-fluid conservative pressure 不再依赖 legacyCouplingPlanFragment。
- [x] migrated production leaves 全部 legacyAdapter=false。
- [x] conservative pressure arithmetic 原样复用。
- [x] 未把 constant-density PressureOperators 强接入 conservative math。
- [x] HOW 是唯一 source scheduler。
- [x] STATE 不插入 calls 或排顺序。
- [x] WHICH 不修改 source HOW order/loops。
- [x] Runtime 不发现/选择 numerical provider。
- [x] provider owner 在编译期冻结。
- [x] state views 仍 lazy，仅绑定实际 demands。
- [x] 未引入 solver-family enum 或 hidden replacement。
- [x] Eulerian/shared pressure 只做明确 Legacy 隔离，未部分迁移数学。
- [x] ARCHITECTURE.md 已按最终源码与实际验证范围更新。
