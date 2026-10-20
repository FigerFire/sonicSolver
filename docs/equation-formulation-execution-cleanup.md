# Equation Formulation + Solver-Family Removal + P0–P2 Closure

日期：2026-09-24  
状态标记只使用：`Implemented` / `Interface-only` / `Legacy` / `Unsupported`。

---

## 1. Before architecture

```text
case input
  models/equations.yaml preset            (compressible / eulerianEulerian)
  models/multiPhase.yaml                  (physics template)
  solvers/algorithm.yaml: type            (densityBase / pressureBase)
  solvers/algorithm.yaml: algorithm       (SIMPLE / PISO / PIMPLE)
        ↓
CaseAdapter → CaseConfig{solver.numerics.solver, solver.pressure.workflow}
        ↓
System::build
  if (singleFluidPreset)        → conservative rho/rhoU/rhoE pack
  else if (composition.declared)→ rhoConst pack
  else if (EulerianEulerian)    → phase packs
  else if (numerics.solver == DensityBased)  ← solver-family dispatch
                                → conservative pack（多相/level-set 桥接）
  else                          → addPressureBasedFluid(...)  ← family dispatch
        ↓
每个 pack 自己硬编码 S_PRESSURE policy / transformation / schedule 映射
        ↓
ResolvedSimulationSystem{ formulation = SolverAlgorithm }
        ↓
runtime 重新读 raw config：numerics.solver / numerics.cfl / maxDeltaT
        ↓
CompressibleAlgorithm（单流体）/ EulerianEulerian::PressureStepper（Eulerian）
```

问题：

```text
1. densityBase/pressureBase 既选方程包，又当运行时身份（双重 authority）。
2. SIMPLE/PISO/PIMPLE 只被当成 runtime loop；注册在守恒系统上时被静默忽略。
3. pressure-coupling 的 formulation/policy/schedule 分散在 system builder 的三个
   pack 里，无法作为"注册模型"报告状态。
4. ResolvedSimulationSystem 暴露 solver-family 字段；state realization 由人读。
5. Eulerian stepper 从 raw config 读 dt/CFL 上限，并每步重建 OpRegistry。
6. solver/system 的 CMake 直接编译 models/**、run/** 的实现文件。
```

---

## 2. Final authority model

```text
Equation Composition (WHAT)
    builtin preset | declared equations | pressure-constraint request
    + model contributions (sources / level-set / turbulence / IBM)
    + user add/extend  (replace/disable = Interface-only, fail-fast)
        ↓
Registered models / coupling presets
    SIMPLE / PISO / PIMPLE  —— requirements 与 resolved equation/constraint 匹配
        ↓  active / inactive / invalid / unsupported
Formulation / Transformation
    pressure-constraint / shared-pressure / immersed-constraint transformers
        ↓
Executable Equation System (WHAT, frozen)
        ↓
State Realization (compiled)
    transported / derived / multiplier role groups + storage/index binding
        ↓
Numerical Compiler (HOW)
    term recipes + CompiledTimeIntegration.recipe + TimeStepPolicy{cfl,maxDeltaT}
        ↓
Solve Planner (ORDER)
    execution policies → structured control-flow IR
        ↓
Runtime Requirements
    providers / services / capabilities（capability 来自真实 binary）
        ↓
Runtime provider binding → SingleFluidStepper / EulerianStepper
        ↓
PlanExecutor（只解释 control flow）
```

`src/ARCHITECTURE.md`（Revision 2026-09-24 v2）已经描述该模型；本报告记录它
在本轮被**实现**到什么程度。

---

## 3. Removed solver-family dispatch

### 3.1 删除的类型 authority

```text
deleted   FDM::SolverAlgorithm { DensityBased, PressureBased }
deleted   FDM::toString(SolverAlgorithm)
deleted   NumericsConfig::solver
deleted   SolverPropertiesConfig::type
deleted   ResolvedSimulationSystem::formulation
```

替代：

```text
CaseConfig::compatFlowLabel          Legacy 兼容输入标签（"densityBase"/"pressureBase"）
CaseConfig::pressureCouplingDeclared 是否显式注册 coupling preset
BuildRequest::pressureConstraint     显式 equation-source 请求
BuildRequest::coupling               显式 coupling preset 注册
ResolvedSimulationSystem::realization  STATE REALIZATION（compiled）
ResolvedSimulationSystem::coupling     coupling preset 状态报告
```

`compatFlowLabel` 只在 composition root（`SF_inspection.cpp`）被翻译成
`pressureConstraint` / `singleFluidPreset` 请求；它不进入
`ResolvedSimulationSystem`，也不参与任何 runtime 分支。

### 3.2 删除的 dispatch 分支

| 位置 | Before | After |
| --- | --- | --- |
| `System::build` 方程包选择 | `else if (numerics.solver == DensityBased) … else addPressureBasedFluid()` | `singleFluidPreset` / `pressureConstraint` / declared composition / Eulerian template（equation source） |
| 多相/level-set 桥接 | 依赖 `numerics.solver` | 依赖 `templateOrigin != SingleFluid`（物理模板） |
| `strictCoreTermCoverage` | `formulation == DensityBased` | `realization.conservativeTransportedMass` |
| homogeneous 执行能力 | `numerics.solver != DensityBased` | `!realization.conservativeTransportedMass` |
| Eulerian 执行能力 | `numerics.solver == PressureBased` | Eulerian template + shared-pressure constraint + finite maxDeltaT |
| `PressureBased::Corrector` guard | `workflow.type != PressureBased` 抛错 | 只校验 typed coupling/linear 配置（调用点由 compiled plan 决定） |
| `EulerianStepper` guard | `workflow.type != PressureBased` 抛错 | `executable system` 必须含 `C_SHARED_PRESSURE` 约束 |
| `SF_report` | 打印 `toString(numerics.solver)` | 打印 equation/coupling 语义 |
| `SF_multiPatch` / `SF_singleFluid` 多相摘要 | 传 `numerics.solver` | 传 `compatFlowLabel`（纯展示） |

---

## 4. Runtime renames

```text
Implemented
  src/solver/algorithm/SF_compressible.h/.cpp
      -> src/solver/algorithm/SF_singleFluidStepper.h/.cpp
      class CompressibleAlgorithm -> class SingleFluidStepper

  src/solver/algorithm/SF_densityBasedRHS.h/.cpp
      -> src/solver/algorithm/SF_conservativeRHS.h/.cpp
      namespace SF::SolverAlgorithm::DensityBasedRHS
      -> namespace SF::SolverAlgorithm::ConservativeRHS

  src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.h/.cpp
      -> src/solver/algorithm/eulerian/SF_eulerianStepper.h/.cpp
      class EulerianEulerian::PressureStepper
      -> class EulerianEulerian::EulerianStepper
```

`CompressibleAlgorithm` / `DensityBasedRHS` / `PressureStepper` 现在出现在
`tools/check_architecture.py` 的禁止列表里：生产源码一旦重新引入这些名字即失败。

`SingleFluidStepper` 的职责被显式限定为：own OpRegistry/workspaces、bind
compiled operations、execute CompiledSolvePlan、commit state。

---

## 5. Pressure-coupling formulation changes

新增（`Implemented`）：

```text
src/solver/system/SF_couplingStatus.h     状态枚举、注册请求、报告、schedule id
src/solver/system/SF_pressureCoupling.h   匹配 + 贡献 API
src/solver/system/SF_pressureCoupling.cpp 实现
```

```text
matchPressureCoupling(request, raw)
   requirements:
     - momentum equation (E_MOMENTUM / E_MOMENTUM.<phase>)
     - incompressibility / pressure-multiplier constraint
       (C_INCOMPRESSIBILITY 或 C_SHARED_PRESSURE)
     - pressure unknown bound to that constraint
   结果：
     active        要求全部满足
     inactive      要求不满足（附原因，例如 "incompressibility /
                   pressure-multiplier constraint not present"）
     invalid       注册的 corrector 数非法等自身矛盾
     unsupported   数学成立但缺 provider/operation

contributePressureCoupling(...)
   active 时贡献：
     1. formulation/transformation 请求（pressureConstraint 或
        sharedPressureConstraint；已请求则不重复）
     2. execution policy（外层/压力/非正交重复次数 + derived operations）
   inactive/invalid/unsupported 一律不贡献，且状态进入 explain/validation。
```

schedule 映射的唯一 authority：

```text
correction.repeatCount        = request.outerCorrectors
correction.nestedRepeatCount  = request.pressureCorrectors
correction.innerRepeatCount   = request.nonOrthogonalCorrectors+1
```

原先分散在 `addPressureBasedFluid` / `addConstantDensityFluid` /
`addEulerianEulerianSolvePolicy` 中的三份 policy 构造代码已删除；system builder
不再出现 `S_PRESSURE` / `S_EE_PIMPLE` 字面量（guard 强制）。schedule id 由
`kPressureScheduleId` / `kSharedPressureScheduleId` 常量提供，planner、
turbulence contribution、printer 引用同一常量。

### 5.1 显式注册不再被静默忽略（Implemented）

`test/IBM/cylinderFlowGhost`（`preset: compressible` + `algorithm: PIMPLE`）：

```text
REGISTERED MODELS
  coupling PIMPLE : inactive
  reason          : incompressibility / pressure-multiplier constraint not present
  requirements    : momentum equation; incompressibility / pressure-multiplier
                    constraint; pressure multiplier unknown;
```

该 case 的数值结果与注册前逐位一致（§15），但"注册了却没生效"不再静默。

`test/pressureConstraintPiso`：`coupling PISO : active`，plan 由 coupling
贡献的 policy 展开，数值与冻结 baseline 字节一致。

---

## 6. State-realization changes

新增（`Implemented`）`src/solver/system/SF_stateRealization.{h,cpp}`：

```cpp
struct CompiledStateRealization {
    bool conservativeTransportedMass;   // rho transported
    bool conservativeMomentum;          // rhoU / U transported
    bool pressureMultiplier;            // p 由约束承担乘法子
    bool thermodynamicPressure;         // p 由 EOS closure 派生
    bool phaseTransportedState;         // 逐相 transported state
    bool multiplierConstraints;
    std::vector<std::string> constraints;
    std::vector<RealizedRoleGroup> transported / derived / multipliers;
};
CompiledStateRealization compileStateRealization(const ExecutableEquationSystem&);
```

它只读 executable equation system 的 `UnknownRole` / `StorageBinding` /
`ConstraintDescriptor`，不含任何 solver-family 开关，出现在
`explain` 的 `STATE REALIZATION` 段。

```text
conservative (sodCase / cylinderFlowGhost)
  conservative transported mass : yes
  transported : rho[1@0,conservative] rhoU[3@1,conservative] rhoE[1@4,conservative]

pressure-constraint (pressureConstraintPiso)
  pressure as multiplier : no（density formulation 的 p 由 EOS closure 派生）
  derived : p[1@0,thermodynamicClosure] pPrime[1@0,pressureCorrection]

constant-density composition（rho=rho0 + div(U)=0）
  pressure as multiplier : yes（p 是约束乘法子）

Eulerian（eulerianEulerianCase）
  phase transported state : yes
  transported : phaseMass.<phase> momentum.<phase> enthalpy.<phase>
```

`RealizedUnknown`（runtime 绑定，`SF_stateRealizer.cpp`）保持不变：编译期
realization 描述 role，运行时 realizer 把它绑定到 `StateBundle` 的 storage。

---

## 7. PIMPLE / rhoConst execution status

```text
pimplle + pressure-constraint（pressureConstraintPiso, PISO/PIMPLE）
    Implemented：coupling active，plan 展开 Predictor / PressureSolve /
    VelocityCorrect / FluxCorrect，numerical regression 字节一致

Eulerian shared-pressure PIMPLE（eulerianEulerianCase）
    Implemented：coupling active，`ee.*` plan 与之前一致

constant-density（rho=rho0 + div(U)=0 + p multiplier）
    mathematical structure: Implemented
    execution provider:       Unsupported
    reason（现在指名缺哪个 operation）:
      "Unsupported: the resolved system carries a pressure-multiplier
       constraint (rho=rho0, div(U)=0) but no active pressure-coupling
       schedule/provider is registered; derived operations [pressure.prepare,
       momentum.assemble, momentum.solve, ...] have no provider."

SIMPLE/PIMPLE fixed-point on single-fluid conservative state
    Unsupported：没有 dedicated predictor provider；不重复施加 full-dt
    momentum step 假装外层迭代（原有语义保持不变）
```

不再出现 "pressure-based solver not implemented" 这类 family 级拒绝。

---

## 8. Compiled HOW authority

```text
Already Implemented（上一轮）
  CompiledNumericalSystem::time.recipe / dt.{cfl,maxDeltaT}
  SingleFluidStepper 只读 numerics_.time.recipe / numerics_.dt.*
  装配一致性守卫（compiled vs raw config 不一致即 fail）

本轮 Implemented
  EulerianStepper 构造函数改为接收 CompiledNumericalSystem
  ee.dt.compute 使用 numerics_.dt.maxDeltaT / numerics_.dt.cfl
  stableTimeStep() 的 probe 上限使用 numerics_.dt.maxDeltaT
  新增同一一致性守卫（Eulerian）
```

physical model constants（gamma / gas constant / mu / Pr / ILW order /
IBM boundary scheme）**没有**被塞进 `CompiledNumericalSystem`；分类见 §16。

---

## 9. Model activation / incompatibility behaviour

```text
注册来源                     结果                        行为
无注册                       coupling inactive           正常执行（无 pressure schedule）
有约束 + 注册 PISO/PIMPLE     active + Unsupported        报告缺 provider/operation
有约束 + 注册且 provider 齐   active + Runnable           正常执行
无约束 + 注册 PIMPLE          inactive + reason           照常执行，不静默忽略
corrector 数非法              invalid                     配置校验期失败
显式注册且 invalid            —                           System::validate() 立即失败
monolithic IB k=KKT 覆盖      inactive + reason           显式说明被 KKT 约束取代
```

`test/test_formulationArchitecture.cpp` 覆盖上述 1–4。

---

## 10. SystemBuilder split

```text
before   solver/system/SF_systemBuilder.cpp   1017 lines（compose + build +
         requirement + explain formatter 混在一起）

after    solver/system/SF_systemBuilder.cpp    679 lines  build/orchestration
         solver/system/SF_presets.cpp          218 lines  built-in equation packs
         solver/system/SF_stateRealization.cpp  compile-time realization
         solver/system/SF_pressureCoupling.cpp  coupling preset matching/contribution
         solver/system/SF_systemDescription.cpp explain/IR enum 字符串
         solver/system/SF_numericalCompiler.cpp numerics compile（原有）
         solver/system/SF_solvePlan.cpp         solve planning（原有）
         solver/system/SF_systemValidator.cpp   invariants（原有）
         solver/system/SF_transformation.cpp    transformers（原有）
```

一个小的公开入口保持不变：`System::build(config, request)`。
没有引入 `BuildManager` / `CompilationContext`。

---

## 11. CMake / dependency cleanup

```text
before  solver/system/CMakeLists.txt 直接编译：
          ../../models/physics/SF_sourceContribution.cpp
          ../../models/ibm/SF_ibmSystemContribution.cpp
          ../../models/turbulence/SF_turbulenceSystemContribution.cpp
          ../../models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.cpp
          ../run/SF_planExecutor.cpp

after   contribution 实现回到 owning target：
          SF_physics     ← SF_sourceContribution.cpp + levelSet contribution
          SF_ibm         ← SF_ibmSystemContribution.cpp
          SF_turbulence  ← SF_turbulenceSystemContribution.cpp
        control-flow 解释器独立成 SF_run（src/solver/run/CMakeLists.txt）
        SF_solverSystem 只消费：link SF_physics / SF_ibm / SF_turbulence / SF_run
        PlanExecutor 不再需要 system 的 out-of-line symbol：
          toString(ExecutionPolicyKind/PlanNodeKind) 改为 SF_solveProgram.h 内联
```

guard：`tools/check_architecture.py` 现在拒绝 `solver/system` 下带 `..` 的
`.cpp` 条目。

Legacy：`SF_solverSystem ↔ SF_turbulence` 仍存在静态库互链（composition
contract 的 include 方向尚未下沉）。见 §19 #1。

---

## 12. Include cleanup

```text
renames 已消除两个重复 basename：
  SF_compressible.h   （algorithm 层）→ SF_singleFluidStepper.h
  SF_pressureStepper.h→ SF_eulerianStepper.h
本轮新增的 eulerian/SF_eulerian.h 目录转发头已删除（避免第三个 SF_eulerian.h）
```

`tools/check_architecture.py` 继续对 duplicate public basename 做 baseline
比较（新增/扩展即失败）——该检查就是 §13 要求的 "new ambiguous public
basename" guard。

---

## 13. P0 / P1 / P2 completion

| 项 | 状态 | 说明 |
| --- | --- | --- |
| P0-A solver-family removal | `Implemented` | enum、字段、dispatch 分支全部删除 |
| P0-B coupling model contribution | `Implemented` | `SF_pressureCoupling` + 状态报告 |
| P0-C constant-density path | `Implemented`(math) / `Unsupported`(provider) | 结构成立，provider 缺失具名报告 |
| P0-D runtime renames | `Implemented` | stepper/RHS 按 runtime domain 命名 |
| P0-E compiled HOW everywhere | `Implemented` | 两个 stepper 都读 compiled time/dt + 守卫 |
| P0-F planner/transformer boundary | `Implemented`（保持） | Planners 只排序；policy 由 coupling/IBM 贡献 |
| P0-G model contribution contract | `Implemented` | 复用 `SystemCompositionBuilder`；未复活 `ISystemContribution` |
| P1-A SystemBuilder contraction | `Implemented` | 见 §10 |
| P1-B CMake ownership | `Implemented` | 见 §11（一个 link 级 legacy 除外） |
| P1-C application execution | `Implemented` | composition root 只装配/绑定 |
| P1-D OpRegistry/workspace lifetime | `Implemented` | 两个 stepper 都用 solver-owned registry |
| P1-E include ambiguity | `Implemented` | renames + baseline guard |
| P2-A native decode split | `Legacy` | 仍单文件；见 §19 #2 |
| P2-B explicit override semantics | `Interface-only` | typed precedence 槽位 + fail-fast；见 §19 #3 |
| P2-C diagnostics extraction | `Legacy` | 见 §19 #4 |
| §17 explain sections | `Implemented` | REGISTERED MODELS / STATE REALIZATION / DERIVED OPERATIONS |
| §20 architecture guards | `Implemented` | 见 §11、§12 与 guard 列表 |

---

## 14. Remaining Unsupported providers

```text
Unsupported  constant-density predictor/corrector provider
             （rho=rho0 + div(U)=0 + p multiplier）
             报告：missing operations [pressure.prepare, momentum.assemble,
             momentum.solve, pressure.boundary.prepare, pressure.assemble,
             pressure.solve, pressure.update.prepare, velocity.correct,
             flux.correct, pressure.correction.commit, pressure.step.commit]
Unsupported  SIMPLE/PIMPLE fixed-point 单流体 predictor provider
Unsupported  distributed surface（Peskin/diffuse kernel）IBM MPI 覆盖缺口
```

---

## 15. Numerical invariance evidence

每个阶段的改动之后都在同一 binary 上复跑并逐位 diff（`Step time` 与
`State closure` 的 dt / min(rho) / min(p) / min(T) / min(c)）：

```text
阶段                                            case                    结果
renames + solver-family removal                 cylinderFlowGhost 20步  IDENTICAL
耦合模型迁移（policy 唯一 authority）            cylinderFlowGhost 20步  IDENTICAL
SystemBuilder split                             cylinderFlowGhost 20步  IDENTICAL
CMake ownership 迁移                            cylinderFlowGhost 20步  IDENTICAL
Eulerian compiled HOW + registry lifetime       cylinderFlowGhost 20步  IDENTICAL
同上一组                                        sodCase 4-rank 20步     IDENTICAL
```

额外的强等价证据：

```text
tools/check_piso_regression.py --case test/pressureConstraintPiso
  -> final VTS SHA-256 与 frozen legacy result 字节一致
```

没有修改任何 CFD 公式、CFL、tolerance、clamp 或 flux 选择。

---

## 16. Raw residual authority classification

| 值 | 当前位置 | 分类 | 结论 |
| --- | --- | --- | --- |
| `numerics.timeRecipe` / `cfl` / `maxDeltaT` | compiled `CompiledNumericalSystem` | numerical HOW | 已迁入 compiled authority |
| `numerics.termRecipes*` | compiled bound terms | numerical HOW | 已在上一轮迁入 |
| `idealGasGamma` / `idealGasConstant` | raw `NumericsConfig` | **EOS/model contract** | 保留在 model/EOS 侧；与 FluidStateModel 做一致性校验（不是 dt/recipe） |
| `dynamicViscosity` / `prandtl` | raw `NumericsConfig` | closure/transport contract | 保留；由 transport/thermo 模型消费 |
| `ibmBoundary` / `ilwOrder` | raw `NumericsConfig` + `BoundaryConfig` | boundary contract + numerical HOW | 仍双处可见；见 §19 #4 |
| `sources.enabled` | `SourceConfig` | model presence | 保留（模型存在性） |
| `pressure.workflow.*` | coupling preset 输入 | formulation 输入 | 只用于注册 coupling preset，不再当 solver identity |

---

## 17. Explain output

`sonicSolver explain <case>` 现在的段落：

```text
FLOW (WHAT)
  template origin / density behavior / thermo compressibility
REGISTERED MODELS
  coupling <PRESET> : active|inactive|invalid|unsupported
  reason / requirements / derived equations / derived operations
STATE REALIZATION
  conservative transported mass / momentum / pressure multiplier /
  thermodynamic pressure / phase transported state
  transported / derived / multipliers role groups（含 storage binding）
TIME RECIPE
COMPILED NUMERICAL SYSTEM
  required halo width / bound terms / workspace / provider requirements
CONTRIBUTIONS
RAW EQUATION SYSTEM
SYSTEM TRANSFORMATIONS
EXECUTABLE EQUATION SYSTEM（STATE VARIABLES / EQUATIONS / CONSTRAINTS）
SOLVE BLOCKS
COMPILED SOLVE PLAN
EXECUTION REQUIREMENTS / WORKSPACE REQUIREMENTS
NUMERICAL PROVIDERS
RUNTIME SERVICES
```

示例（守恒系统上注册 PIMPLE，显式 inactive，见 §5.1）。

---

## 18. Architecture tests

新增 `test/test_formulationArchitecture.cpp`（ctest：`formulationArchitecture`）：

```text
1. built-in pressure-constraint 与 manual Momentum+rhoConst 都产生
   C_INCOMPRESSIBILITY + E_MOMENTUM（同一约束/动量结构）
2. 注册 PIMPLE → coupling active 且 plan 含 outerCorrectors/pressureCorrectors 循环
3. 守恒 transported-rho 系统注册 SIMPLE → inactive（附原因），且不泄漏
   pressure.* 操作（不存在 density-solver 分支）
4. state realization 由 role 导出（conservative vs pressure multiplier）
5. 同一组方程换 template provenance 不改变 realization
6. compiled dt/time recipe 来自 CompiledNumericalSystem
7. Unsupported 报告必须指名 provider/operation + missing operations
8. preset 与手写 IR 的 equation definition 结构一致
9. 显式注册永不静默忽略（inactive 也必须带 reason）
```

既有测试同步更新为新 authority（`test_pisoArchitecture.cpp`、
`test_termRecipes.cpp`）：不再写 `numerics.solver` / `workflow.type`，改为
注册 typed coupling request。

---

## 19. Remaining architecture debt

| # | 项 | 状态 | 说明 / 建议 |
| --- | --- | --- | --- |
| 1 | `models ↔ solver/system` composition contract | `Legacy` | `SF_equationContribution.h`（builder）+ `RawEquationSystem`/`ExecutionPolicy` 类型仍在 solver/system，导致 `SF_solverSystem ↔ SF_turbulence` 静态库互链。建议把这些 WHAT/ORDER value types 下沉到 `core/system`，之后 models 只依赖 core |
| 2 | native decode split（P2-A） | `Legacy` | `SF_nativeDecode.cpp` 仍单文件（~3.3k 行，已按语义分节且不依赖 parser globals）。拆分前置条件：把 JSON 小工具抽到 `compatibility/private/` 头 |
| 3 | equation override lowering（P2-B） | `Interface-only` | 已建立 typed precedence 槽位（builtin → models → user add/extend → replace/disable → validation）；add/extend/replace/disable 的 lowering 尚未实现，当前一律 fail-fast，不静默 |
| 4 | Eulerian diagnostics extraction + `ibmBoundary`/`ilwOrder` 双处可见（P2-C） | `Legacy` | `ee.step.commit` 内仍有 cell traversal/phase 诊断；`ibmBoundary`/`ilwOrder` 同时存在于 numerics 与 boundaries |
| 5 | 高阶生产算例数值缺陷 | `Legacy` | 4-rank TENO5 Sod 在 t≈0.0746 negative pressure fail-fast（既有缺陷，见 docs/high-order-first-divergence.md），本轮未处理 |
| 6 | distributed surface IBM MPI | `Unsupported` | `cylinderFlowPeskinMPI4` 首步无法推进（既有覆盖缺口） |
| 7 | `SF_eulerianStepper` 的 `config_` 其余读取 | `Legacy` | phase source CFL / phase convection 仍在 raw workflow；它们属于 coupling preset 输入而非 time/dt authority |
| 8 | CMake 引用的空目录 | `Legacy` | `src/methods/numerics/time`、`src/solver/discretization/time` |

---

## 20. Acceptance criteria

```text
[x] `CompressibleAlgorithm` is gone from production code.
    guard 拒绝该 token；文件/类/命名空间已改名
[x] `DensityBasedRHS` is gone from production code.
    -> ConservativeRHS（RHS 描述的是守恒状态数学，不是密度求解器）
[x] single-fluid runtime stepper is named by runtime domain, not compressibility.
    -> SingleFluidStepper
[x] Eulerian runtime stepper is not named by pressure strategy.
    -> EulerianEulerian::EulerianStepper
[x] compressible/incompressible preset provenance does not select runtime paths.
    方程来源只声明 WHAT；guard 拒绝 preset provenance → stepper 分支
[x] DensityBased/PressureBased no longer select physical equation packs.
    enum/字段删除；compatFlowLabel 只在 composition root 翻译成 equation request
[x] state realization is compiled from executable equations/formulations.
    CompiledStateRealization（SF_stateRealization.cpp）
[x] SIMPLE/PISO/PIMPLE contribute formulation + solve-plan semantics.
    SF_pressureCoupling：transformation 请求 + execution policy
[x] no hidden SIMPLE/PISO/PIMPLE solver lifecycle exists inside a stepper.
    stepper 只 bind operations 并执行 CompiledSolvePlan
[x] rhoConst + incompressibility can be transformed into a truthful pressure system.
    C_INCOMPRESSIBILITY + p multiplier + pressureConstraint formulation
[x] missing pressure execution capability fails only as a named Unsupported provider.
    reason 列出 derived operations（见 §7 / §14）
[x] explicit incompatible coupling registration is reported, never silently ignored.
    explain: `coupling PIMPLE : inactive` + reason（cylinderFlowGhost / sodCase）
[x] all steppers consume compiled time/dt/numerical authority.
    SingleFluidStepper + EulerianStepper，均有 compiled-vs-raw 一致性守卫
[x] solver/system CMake no longer compiles foreign model/run `.cpp`.
    SF_physics / SF_ibm / SF_turbulence / SF_run 各自拥有实现
[x] PlanExecutor remains an ORDER interpreter only.
    只解释 Sequence/Loop/StageLoop 与 OpId；不派生方程
[x] no CFD formula was silently changed.
    逐位一致的回归证据（§15）
[x] Sod 20-step regression is unchanged.
    mpirun -np 4 (clean build) explain + --steps 20 -> IDENTICAL
[x] cylinderFlow 20-step regression is unchanged.
    explain + --steps 20 -> IDENTICAL
[x] clean build passes.
    cmake -S . -B /tmp/sonic-clean2 -G Ninja -DCMAKE_CXX_COMPILER=clang++
      -DBUILD_CLI=ON -DBUILD_GUI=OFF -DBUILD_TESTS=ON   -> HYPRE + MPI enabled
    cmake --build /tmp/sonic-clean2 -j 8                 -> 0 errors
[x] ctest passes.
    ctest --test-dir /tmp/sonic-clean2 -> 14/14 Passed
    （含新增 formulationArchitecture）
[x] final report is produced.
    docs/equation-formulation-execution-cleanup.md
```

---

## 21. Evidence index

---

## 22. Phase 25 — Executable Operation Authority

> 注：本轮收到的 spec 文本在 §4 开头（"Do not create a parallel …"）被截断，
> 以下实现对应 §0–§3 已明确的部分（唯一 operation authority + 责任边界），
> 以及 §4 已可推断的目标（formulation 声明 executable operations）。§4 之
> 后的 P0-A… 具体 API/命名/验收细则尚未收到。

### 22.1 收敛前的三份 authority

```text
1. formulation      SF_transformation.cpp  PressureConstraintTransformer
                    -> OP_PRESSURE_UPDATE / OP_VELOCITY_CORRECTION /
                       OP_FLUX_CORRECTION（symbolic）
2. coupling         SF_pressureCoupling.cpp couplingOperations()/
                    couplingProviders() 手工 operation/provider 名单
3. planner          SF_solvePlan.cpp compileGenericPiso()/
                    compileUnavailablePressureSchedule() 手工 operation 名与顺序
```

### 22.2 收敛后的唯一 authority

```text
Formulation / Transformation（唯一 authority）
    SF_transformation.cpp 声明 ExecutableOperation：
      {OpId, stage, provider, requiresLinearSolver, provenance}
    -> ExecutableEquationSystem::operations

Coupling preset（只引用 stage）
    SF_pressureCoupling.cpp couplingPlanFragment()
      顺序 + 重复次数 + OperationStage 引用
      formulation 未声明的 step 只能给"具名缺失标记"
    -> PlanFragment（SF_planFragment.h）

SolvePlanner（只 merge/order/解析）
    lowerFragmentNode(): stage -> OpId（通过 findExecutableOperation）
    删除 compileGenericPiso / compileUnavailablePressureSchedule
    planner 源码不再出现 pressure/momentum/velocity/flux operation 名字

Provider resolver
    couplingPlanResolved() 报告 formulation 未声明、需要 provider 的
    operation ids；builder 据此 require provider / runtime service，并把
    缺失清单写进 Unsupported 报告
```

关键文件：

```text
new  src/solver/system/SF_planFragment.h        片段 IR（Sequence/Loop/Leaf）
new  ExecutableOperation / OperationStage       SF_resolvedEquationSystem.h
new  addExecutableOperation / findExecutableOperation
del  couplingOperations() / couplingProviders()
del  compileGenericPiso() / compileUnavailablePressureSchedule()
```

### 22.3 行为等价性

```text
pressureConstraintPiso（PISO）
  coupling PISO : active
  plan node id / 顺序 / OpId 与收敛前完全一致
    PISO.prepare -> PISO.outerCorrectors -> PISO.pressureCorrectors ->
    PISO.nonOrthogonalCorrectors -> ... -> PISO.commit
  tools/check_piso_regression.py -> final VTS 与 frozen legacy 字节一致

固定点 schedule（SIMPLE/PIMPLE on pressure-constraint system）
  片段引用具名缺失标记 pressure.schedule.*（与收敛前同一组 id）
  runtime status = Unsupported，missingOperations 不变

monolithic KKT（FullyImplicitDLM / DFMAugmentedLagrangian / …）
  coupling PIMPLE : inactive（superseded by the monolithic KKT/DLM
  constraint contribution），plan 仍为 ibm.kkt.solve / constraint.project

守恒系统（sodCase / cylinderFlowGhost）
  coupling PIMPLE : inactive（incompressibility constraint not present）
```

### 22.4 本轮同时修正的问题

```text
strict core-term coverage 的 authority 推导过粗：
  phase 24 用 "realization.conservativeTransportedMass" 作为
  NumericalCompiler 的 strict 开关，导致 pressure-constraint 方程族
  （rho/rhoU/rhoE transported + C_INCOMPRESSIBILITY）被要求
  "每个选中的 recipe 都被 core equation 使用"，从而让
  cylinderFlowFullyImplicitDLM / cylinderFlowDFMAugmentedLagrangian
  在配置校验期失败（"diffusion TermRecipe was selected but the core
  Equation System contains no diffusion term"）。

now:
  strict = conservativeTransportedMass && !hasConstraint(
               executable, "C_INCOMPRESSIBILITY")
  即只有守恒 density formulation 才要求 recipe 使用；pressure-constraint
  方程族的 diffusion/pressure 项由 coupling provider 施加。这两个 case
  恢复为可构建并通过（plan 为 KKT/constraint-projection，runnable）。
```

### 22.5 新增 guard 与测试

```text
tools/check_architecture.py
  - coupling 层不得出现 couplingOperations/couplingProviders 或裸
    pressure/momentum operation id
  - SolvePlanner 不得出现 compileGenericPiso /
    compileUnavailablePressureSchedule 或裸 pressure operation id
  - formulation 必须声明 addExecutableOperation /
    OperationStage::PressureSolve / OperationStage::VelocityCorrect
  - fixed-point 缺失标记必须来自 coupling 片段

test/test_formulationArchitecture.cpp（第 9 组断言）
  - formulation 声明 PressureSolve -> "pressure.solve"（provider +
    requiresLinearSolver）
  - PISO 片段在 executable operations 上 100% 解析
  - plan 里每个 operation 要么由 formulation 声明，要么是
    pressure.schedule.* 具名缺失标记
  - fixed-point 片段必须 unresolved，且缺失清单含
    pressure.schedule.predictor.solve、不含 momentum.solve
```

### 22.6 验证（clean build）

```text
configure + build（/tmp/sonic-p25, Ninja, HYPRE+MPI）  -> 0 errors
ctest --test-dir /tmp/sonic-p25                         -> 14/14 Passed
tools/check_architecture.py                             -> 0 ERROR
cylinderFlowGhost --steps 20（clean binary）            -> 与 baseline 逐位一致
mpirun -np 4 sodCase --steps 20（clean binary）         -> 与 baseline 逐位一致
pressureConstraintPiso                                  -> final VTS 字节一致
cylinderFlow*（全部 serial）--steps 2                    -> exit 0
mpirun -np 4 cylinderFlowFictitiousDomainMPI4 --steps 2 -> exit 0（IBM 诊断正常）
eulerianEulerianCase --steps 2                          -> 与收敛前 summary 一致
```

### 22.7 尚未完成（等待完整 spec）

```text
Legacy   Eulerian `ee.*` plan 仍在 SolvePlanner（compileEulerianPimple）。
         它的 operation id 已是 formulation/provider 侧命名，但完整
         "shared-pressure formulation 声明 executable operations +
         fragment" 迁移尚未做。
Legacy   PlanExecutor 侧没有独立的 provider resolver 类型；当前解析结果
         由 builder 写入 RuntimeReport（missingOperations/reason）。
```

```text
clean build
  cmake -S . -B /tmp/sonic-clean2 -G Ninja -DCMAKE_CXX_COMPILER=clang++ \
        -DBUILD_CLI=ON -DBUILD_GUI=OFF -DBUILD_TESTS=ON
  cmake --build /tmp/sonic-clean2 -j 8

tests
  ctest --test-dir /tmp/sonic-clean2 --output-on-failure        -> 14/14
  python3 tools/check_architecture.py                           -> no ERROR

regression (clean-built binary)
  mpirun -np 4 /tmp/sonic-clean2/sonicSolver explain test/Sod/sodCase
  mpirun -np 4 /tmp/sonic-clean2/sonicSolver run --steps 20 test/Sod/sodCase
      -> dt / min(rho) / min(p) / min(T) / min(c) 与改动前逐位一致
  /tmp/sonic-clean2/sonicSolver explain test/IBM/cylinderFlowGhost
  /tmp/sonic-clean2/sonicSolver run --steps 20 test/IBM/cylinderFlowGhost
      -> 同上逐位一致（IBM + ILW + turbulence + PIMPLE/HYPRE + RK4/WENO5/LF）
  python3 tools/check_piso_regression.py --solver <build>/sonicSolver \
        --case test/pressureConstraintPiso
      -> final VTS 与 frozen legacy result 字节一致
  /tmp/sonic-clean2/sonicSolver run --steps 2 test/eulerianEulerianCase
      -> exit 0（shared-pressure coupling active）
```
