# SonicSolver Executable System Audit

日期：2026-09-16

本文记录本轮修改前的 runtime authority。审计范围包括 `solver/system`、
`solver/equation`、workflow/time/discretization/boundary/linear algebra、
`core/state`、models/IBM、application run composition 和 execution runtime。

## 1. Authority map

| 对象 | Resolved authority | Runtime authority | 当前差距 |
|---|---|---|---|
| density mass/momentum/energy | `equationDefinitions` + density `AssemblyPlan` | `Equation::Compressible::System` adapter + `DensityBasedRHS` | 已由 resolved term gate backend；具体 operator contract 仍在 adapter/RHS |
| Eulerian phase continuity/momentum/enthalpy | resolved definitions | `PhaseEquationAssembler` 固定方法与 `PressureStepper::stepImpl` | equation identity 已 resolved，term 与调用顺序仍硬编码 |
| turbulence transport | resolved unknown/equation/`S_TURBULENCE` | `Turbulence::EquationSystem` storage/kernel；PressureStepper 直接调用 | model-specific equation storage 与 solve 调用仍隐藏在 specialized objects |
| homogeneous phase mass | `E_PHASE_MASS` definition | `HomogeneousPhaseChangeCoupling` lifecycle hooks | resolved equation 尚无 AssemblyPlan/独立 executor |
| legacy alpha | `E_LEGACY_ALPHA` definition | `LegacyMultiphaseEquationCoupling` | equation execution 仍由 hook 隐式完成 |
| level set | `E_LEVEL_SET` definition | `InterfaceEquationCoupling`/level-set model | advection/reinitialization/commit 仍是 lifecycle hook |
| IBM Ghost/ILW | immersed descriptor + boundary solve block | boundary pipeline/IBM adapter | boundary requirement 仍未关联到具体 operator plan |
| IBM forcing/KKT | resolved unknown/constraint/solve block | immersed strategy、constraint port、linear algebra backend | row/unknown identity仍部分由 IBM private object构造；无 solve-block registry |
| source | resolved right-side terms + registration tables | density/Eulerian source registries | selection 已无 central switch；provider handle 尚未进入 AssemblyPlan |

## 2. 九个审计问题

### 2.1 已完全由 resolved definition 驱动的 equation

density `E_MASS`、`E_MOMENTUM`、`E_ENERGY` 已由 resolved definition 创建，
compressible adapter 只借用其 term 并建立 AssemblyPlan。Convection、diffusion、
source backend 是否执行由 plan 中的 `TermKind` 决定。具体 WENO/viscous/source
kernel 仍是 specialized backend，这是合法的“怎么解”。

### 2.2 仍隐藏在 specialized object 的 runtime equation identity

- `Turbulence::EquationSystem` 仍按 model object 返回 k/epsilon/omega transport
  equation objects；resolved system 虽已声明相同数学 equation，二者尚无 stable
  handle binding。
- `PhaseEquationAssembler` 的 continuity/momentum/enthalpy 方法隐含 phase equation
  分类。
- homogeneous phase mass、legacy alpha 和 level-set 的 execution identity 仍由
  `IEquationSystemCoupling` 实现类型决定。
- IBM constraint rows、J/JT/lambda storage identity 仍部分存在于 IBM private
  method state。

### 2.3 specialized assembler 硬编码的 term

Eulerian continuity 的 transient/face-flux/source，momentum predictor 的 transient、
pressure gradient、interphase/source，enthalpy transport/source，pressure correction、
turbulence solve 均由 `PhaseEquationAssembler`/`PressureStepper` 固定调用。Density
operator 的精确 boundary/halo/canonical barriers 仍由 `DensityBasedRHS` 固定实现。

### 2.4 physics-specific code 手工创建的 unknown/storage

- conservative storage 由 `Field` 与 EquationSet 初始化；
- per-phase primary/primitive storage 由 `PhaseSystem` 创建；
- turbulence scalar storage 由 turbulence manager/equation system 创建；
- level-set/legacy/homogeneous auxiliary storage由各 model 创建后手工注册；
- pressure correction、momentum diagonal、face flux、previous-time state 由
  `PhaseSolverWorkspace::setupLike` 创建；
- IBM multiplier/constraint storage 由 IBM method state 创建。

`UnknownDescriptor` 当前只表达 ID、名字、位置、分量和 ownership，不能声明
role、初始化、边界、restart/output 或 runtime-storage requirement。不存在统一
state realization/validation。

### 2.5 stepper 中硬编码的 solve order

- `CompressibleAlgorithm` 固定 dt -> transport correction -> RK/Euler -> optional
  flow correction -> commit/tail。
- `PressureStepper::stepImpl` 固定 outer PIMPLE、continuity、momentum、pressure
  correctors、energy、turbulence、commit 顺序。
- workflow stages 只在进入时间循环前验证；没有通用 solve-block executor registry。

### 2.6 隐式 workspace requirement

`PatchWorkspace`、`PhaseSolverWorkspace`、pressure matrix、turbulence matrix、IBM
constraint/KKT workspace 均由 specialized executor 自己分配。ExecutionRuntime
看到 canonical/global workspace access 后才在 timestep 中报
`workspace synchronization must be invoked with explicit borrowed storage`，说明
required borrowed view 没有在 composition 阶段收集和验证。

### 2.7 硬编码 output/diagnostic field

`SF_output.cpp` 明确枚举 Eulerian pressure、alpha/rho/T/U、phaseMass、enthalpy、
source ledger、wall boiling、turbulence；`SF_singleFluid.cpp` 明确枚举 momentum、
energy、thermodynamic cache、IBM mask/multiplier。新增 unknown 不能自动进入
restart/output/diagnostics。

### 2.8 仍依赖 solver path 的 boundary/numerics

Density 在 `CompressibleAlgorithm`/`DensityBasedRHS` 中选择 WENO、flux、viscous
与 boundary pipeline；Eulerian 在 PhaseEquationAssembler 中使用自己的 transport
与 pressure discretization。尚无 per-term resolved numerics 或 BoundaryRequirement。
高阶 characteristic path 仍明确只支持 five-variable PerfectGas，并正确 fail fast。

### 2.9 仍使用 lifecycle hook 的 module

`IEquationSystemCoupling::{beginStep,prepareRHS,assembleRHS,
preparePressureCorrection,commitStep}` 承载 homogeneous phase mass、legacy alpha、
level set 和 phase change。Turbulence transport 通过 transport/model correction 与
PressureStepper direct call 进入 lifecycle。它们尚未全部降为 term/equation/
constraint/solve-block contribution。

## 3. 当前 generic execution gap

```text
ResolvedSimulationSystem
  -> Workflow::Plan (same stage values)
  -> runFlow binds stages
  -> specialized INavierStokesStepper
       -> fixed internal lifecycle
```

缺少的启动阶段 contract：

```text
Unknown realization + stable state handles
AssemblyPlanRegistry + lowering status
WorkspaceRequirement collection/binding
SolveBlockExecutorRegistry
SystemExecutor
Output/Diagnostic providers
```

`SF_executionBuilder.cpp` 仍通过 `S_EE_PIMPLE` 选择两个 composition entry。该判断
基于 solve strategy 而非 template provenance，但仍是 central specialized dispatch。

## 4. 最小安全迁移顺序

1. 扩展 unknown contract，并建立只绑定/验证既有 storage 的 StateRealizer；不改变数组布局。
2. 建立 AssemblyPlanRegistry，显式标记 generic plan 或 specialized lowering；禁止 silently ignored equation。
3. 为 SolveBlock 建立 strategy/capability registry 和启动期完整性检查。
4. 让 Eulerian assembler 请求 resolved plan term，保留原公式与顺序。
5. 显式收集/bind workspace views，先消除 Eulerian runtime-late borrowed-storage failure。
6. 在上述 contract 有真实 consumer 后引入 SystemExecutor；不创建空壳 facade。
7. 最后迁移 lifecycle coupling、output/diagnostic 和 acceptance modules。

Phase D 以后若需要同时改变 phase/PIMPLE/IBM 数值公式才能完成，将停止并记录
blocker，不以 metadata 或 no-op executor 宣称完成。
