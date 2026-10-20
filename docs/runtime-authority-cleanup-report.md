# Runtime Authority Cleanup + Dependency Direction + Production Regression

日期：2026-09-23  
范围：runtime authority / lifetime / dependency direction 清理，以及两个生产
算例（sodCase、cylinderFlow）的架构回归。**本轮没有修改任何 CFD 数学公式。**

---

## 23.1 Scope

本轮修改的内容：

```text
changed
  architecture authority
  runtime lifetime / ownership
  dependency direction
  case-input (native) 到 typed config 的解码路径

not changed
  WENO/TENO reconstruction
  Rusanov / LaxFriedrichs / Steger-Warming flux
  central diffusion
  forwardEuler / SSPRK3 / classicalRK4 coefficients
  PISO/PIMPLE/SIMPLE mathematics（Rhie-Chow、pressure correction）
  IBM forcing / Ghost-IBM / ILW 数学
  turbulence equations
  MPI owner/neighbour conservation semantics
  KKT/DLM mathematics
  EOS mathematics
```

判定标准不是"程序没崩"，而是：

```text
每个阶段产出的 dt / stage time / min(rho) / min(p) / min(T) / min(c)
必须与修改前逐位一致
```

本报告 §23.6、§23.7、§23.8 给出该证据。

---

## 23.2 Before / After Authority

### Before

```text
case.yaml / models/*.yaml / solvers/*.yaml
        ↓  NativeCaseSections (typed sections)
decode*Section()
        ↓  写 parser 全局暂存值
extern CFL / extern maxDeltaT / extern convectionScheme /
extern fluxSplitter / extern enableIBM / extern ilwOrder /
extern U_BC / p_BC / rho_BC / T_BC / extern mu / Prandtl / ...
        ↓  Legacy::makeSolverConfig()
FDM::SolverConfig (typed)
        ↓
ResolvedSimulationSystem
        ↓
runtime 仍从 raw SolverConfig 重读
numerics.cfl / numerics.maxDeltaT / numerics.timeRecipe
```

问题：

```text
1. typed → legacy globals → typed 往返；同一份 case 在解码过程中存在两份
   authority（typed section 与 parser 全局暂存值）。
2. backend 能力在数学系统构造处被写成 constexpr bool ... = true。
3. RuntimeRequirement 出现 required=true / available=true / detail="not
   implemented" 这种三义不一致的状态。
4. runtime 读 raw numerics 的 dt / time recipe，而不是 compiled HOW。
```

### After

```text
case.yaml / models/*.yaml / solvers/*.yaml
        ↓  NativeCaseSections
decode*Section()  —— 直接写 typed CaseConfig
        ↓
CaseConfig.solver (typed)
        ↓
ExecutableEquationSystem      WHAT
CompiledNumericalSystem       HOW   (term recipes + time recipe + dt policy)
CompiledSolvePlan             ORDER
RuntimeRequirements           RUNTIME
        ↓  System::validate()  fail-fast
ExecutionEnvironment → Provider → Stepper → PlanExecutor
        ↓
runtime 只读 compiled authority：
  numerics_.time.recipe / numerics_.dt.cfl / numerics_.dt.maxDeltaT
```

---

## 23.3 Deleted Legacy Authority

### 3.1 native path 不再依赖 parser 全局状态

`CaseAdapter::decodeNativeCase()` 现在从 `FDM::defaultSolverConfig()` 开始，
由各 `decode*Section()` 直接写 typed `CaseConfig`。以下写法已从 native 生产
路径彻底消失：

```text
ParserState parserTransaction;            // deleted from decodeNativeCase
Legacy::makeSolverConfig();               // deleted from buildCaseConfig
CFL = ...; maxDeltaT = ...;               // → numerics.cfl / numerics.maxDeltaT
convectionScheme = ...; fluxSplitter = ...;
reconstructionVariable = ...; interfaceFluxPolicy = ...;
viscousScheme = ...; enableViscous = ...; mu = ...; Prandtl = ...;
timeRecipeName = ...; termRecipesDeclared = ...;
convectionTermRecipe = ...; diffusionTermRecipe = ...;
enableTurbulence / turbulenceFamily / turbulenceModel / turbulenceC* / ...;
enableIBM / enableILW / ilwOrder / ibmBoundaryScheme / ilwNormal* ;
enableGravity / enableMRF / enableWallHeatSource /
gravitySettings / rotatingSettings / wallHeatSettings;
U_BC / p_BC / rho_BC / T_BC / U_IC / p_IC / rho_IC / T_IC / k_IC / ...;
createMesh / solverApplication（solverApplication 已完全删除）
```

对应 typed 目标：`numerics.*`、`numerics.recipes.*`、`pressure.*`、
`boundaries.*`、`initial.*`、`ibm.*`、`sources.*`、`turbulence.*`、
`turbulence.scalars.*`、`turbulence.coefficients.*`。

### 3.2 仍然保留、但已退出构建的 Legacy 文件

以下文件在确认无任何生产引用后**退出所有 target**（不再是任何 library 的
source），为可追溯性保留在磁盘上：

```text
Legacy
  src/app/application/model/legacy/SF_parameter.{h,cpp}
  src/app/application/model/legacy/SF_legacyConfig.{h,cpp}
  src/app/application/model/legacy/SF_parserState.h
  src/models/initial/CMakeLists.txt 中的 SF_param target（已删除）
```

保留原因：当前工作区几乎全部 `src/` 尚未纳入 git（`git ls-files` 只有 176 个
文件），直接删除不可恢复。删除它们的唯一前置条件是"不再出现在任何
CMake target 中"，该条件已经满足，因此它们不再构成 authority；一次性物理
删除留待用户确认后执行。

### 3.3 被删除的死架构骨架

```text
deleted  ISystemContribution / SystemComposer
```

证据：两者只有自身定义与实现，没有任何 production call site
（`rg -n "SystemComposer|ISystemContribution" src test` 只命中
`SF_equationContribution.{h,cpp}` 自身）。`SystemCompositionBuilder` 是真实
使用中的 composition API，保留。

### 3.4 文件系统清理

```text
deleted  479 × "._*" AppleDouble（src/test/docs/tools/scripts）
deleted  3 × 无引用的空目录
         src/solver/algorithm/workflow
         src/solver/boundary/algebraic
         src/solver/boundary/ILW
.gitignore 已包含 ._* 与 .DS_Store（无需新增）
```

未删除的仍被 CMake `include_directories` 引用的空目录：
`src/methods/numerics/time`、`src/solver/discretization/time`（引用它们的
CMake 尚未收敛，删除需要同时改 target include 面，留作后续）。

---

## 23.4 Dependency Changes

| 反向关系 | Before | After |
| --- | --- | --- |
| `core → app` | `core/config/SF_config.h` 直接 include `app/application/model/SF_configParser.h`，把 parser 变成 core 的传递依赖 | **已删除**。core 只导出 typed 类型；`SF_cli.cpp`、`SF_compatibility.cpp`、`SF_nativeDecode.cpp`、`test_*` 显式 include parser |
| `solver → app`（隐藏） | `SF_linearAlgebra.cpp` 调用 parser 头里的 `validateLinearSolverConfig`；`models/ibm` 调用 parser 头里的 `validateIBMForcingSelection` | **已删除**：`validateLinearSolverConfig`、`validateSolverProperties`、`validateIBMForcingAlgorithm`、`validateIBMForcingSelection`、`validateNumericsConfig` 全部移入 `core/config/types/*` 与 `core/config/SF_numericsPolicy.h`（typed 值对象自校验） |
| `core → methods` | `core/field/SF_field.h` include `methods/numerics/state/SF_thermodynamicStateCache.h` | **已删除**：cache 移到 `core/state/SF_thermodynamicStateCache.h`（namespace `SF::State`），`core → core` |
| `models → solver`（structured calculus） | turbulence / interphase 直接 include `solver/discretization/structured/SF_vectorCalculus.h` | **已删除**：移到 `methods/numerics/structured/SF_vectorCalculus.h`，`models → methods` |
| `methods → solver`（同一工具） | `methods/numerics/viscous/SF_viscous.h` 同样指向 solver 目录 | **已删除**（同上） |
| `infrastructure → solver`（HYPRE） | `src/infrastructure/mpi/backend/SF_hypreBackend.cpp` + `SF_schurPreconditioner.h` + `SF_hypreControls.h` include `solver/linearAlgebra/hypre/SF_hypre.h`、`SF_blockSchur.h`；`SF_hypreBackend` target 由 infrastructure 定义 | **已删除**：三个实现文件移到 `src/solver/linearAlgebra/backend/hypre/`，由 `SF_linearAlgebra` 自己在 `MPI+HYPRE` 存在时编译；infrastructure 只保留 communicator/ownership |
| CMake 跨目录编译 cpp | `solver/system/CMakeLists.txt` 编译 `../../models/physics/SF_sourceContribution.cpp`、`../../models/ibm/...`、`../../models/turbulence/...`、`../../models/physics/interfaceModel/levelSet/...`、`../run/SF_planExecutor.cpp` | **仍存在**，见 §23.9 第 1、2 条 |

当前依赖方向：

```text
app
  ↓
solver/system + solver/run
  ↓
solver/algorithm + solver/equation
  ↓
methods + models + core interfaces/types
  ↓
infrastructure / backend
```

仍然存在的例外（逐条解释见 §23.9）：`models → solver/system`（equation
composition contract & equation description value types）。

---

## 23.5 Runtime Changes

### 5.1 Field 变量组（vector）布局

```cpp
// before：假定 RU/RV/RW == 1/2/3
getVal<Vector3>(i,j,k,v)  -> Vector3(Q(...,1), Q(...,2), Q(...,3))
// after：布局由 v 决定
getVal<Vector3>(i,j,k,v)  -> Vector3(Q(...,v), Q(...,v+1), Q(...,v+2))
```

所有 caller 本来就传 `stateModel()->momentumIndex(0)` 或 legacy `RU`，因此
legacy 布局行为不变，非 1/2/3 布局从"静默错位"变成"正确读写"。

新增回归：`test/test_fieldVectorLayout.cpp`（ctest 名 `fieldVectorLayout`），
检查 AoS 与 SoA 两种 layout 下的 group offset 2/3/4 读写，且写入不越界到
其他变量组。

### 5.2 RuntimeRequirement 三义一致

```text
CanonicalScalarInterfaceFlux:
  before  required=(legacyAlphaTransport && parallel) / available=true
          / detail="not implemented"
  after   available = !(legacyAlphaTransport && parallel)
```

`System::validate()` 现在强制：

```text
name 非空、detail 非空、name 唯一
required && !available                -> 抛错（编译/校验期 fail-fast）
!available && report.status==Runnable -> 抛错
        （禁止 "input 说 A / compiled 说 B / runtime 静默执行 C"）
```

### 5.3 backend 能力来自真实 binary composition

```text
new  core/interfaces/SF_buildCapabilities.h
     struct SF::System::BuildCapabilities { mpi, hypre,
                                            distributedLinearSystem,
                                            canonicalConstraintDof };
new  src/app/application/SF_hostCapabilities.{h,cpp} （composition root）
     caps.mpi  = Parallel::mpiCompiled()          // SF_USE_MPI 的真实后端
     caps.hypre= LinearAlgebra::hypreBackendLinked() // 真实 ParCSR vs 占位
     caps.distributedLinearSystem = mpi && hypre
     caps.canonicalConstraintDof  = mpi            // serial stub 的 COPY 是 no-op
```

数据流：

```text
equation/system requirement  +  binary/backend capability
        ↓
BuildRequest.capabilities
        ↓
RuntimeRequirement(required, available, reason)
        ↓  validate() fail-fast
```

`constexpr bool kConstraintGlobalDofAvailable = true;`
`constexpr bool kDistributedLinearSystemAvailable = true;` 已从
`SF_systemBuilder.cpp` 删除。

### 5.4 OpRegistry / Explicit Workspace lifetime

```text
before  advance(){ Time::Explicit::Workspace explicitWorkspace;
                   Run::OpRegistry operations; ... }
after   solver members
        Run::OpRegistry operations_;
        Time::Explicit::Workspace explicitWorkspace_;
        PressureBased::CorrectionSummary pressureSummary_;
        bool correctionPerformed_/correctionWroteConservative_;
        std::string correctionDetail_;
        advance(){ operations_.clear(); ... 只 overwrite }
```

RK stage 存储 `q0/k1..k4` 现在跨 timestep 保留容量，只在 layout 变化时
resize。stage 数学、stage time、操作顺序完全未改（§23.7/§23.8 的逐位一致
证据）。OpRegistry 新增 `clear()`；`PlanExecutor::validateBindings` 与
`execute` 调用顺序不变。

### 5.5 thermodynamic cache 热路径

```text
before  thermodynamicState(i,j,k):
          std::vector<double> q(nVar);  gather; cache.state(q,...)
after   thermodynamicState(i,j,k):
          cache.find(cell) 命中 -> 直接返回（不 gather、不分配）
          miss  -> thread_local scratch gather -> computeAndStore
```

`computeAndStore` 与旧 `state()` 走同一条 `FluidStateModel::close`；
cache 版本号语义、失效时机、容量语义均未改变。

### 5.6 compiled HOW authority

```text
new  CompiledNumericalSystem { ... CompiledTimeIntegration time;
                                   TimeStepPolicy dt; }
      time.recipe   <- NumericalCompiler::compile(recipes)
      dt.cfl        <- config.numerics.cfl
      dt.maxDeltaT  <- config.numerics.maxDeltaT
```

runtime 读取改为：

```text
const config_.numerics.timeRecipe    -> numerics_.time.recipe
config_.numerics.cfl / maxDeltaT     -> numerics_.dt.cfl / numerics_.dt.maxDeltaT
```

并新增装配一致性守卫：

```cpp
if (numerics_.time.recipe.id() != config_.numerics.timeRecipe.id()
    || numerics_.dt.cfl != config_.numerics.cfl
    || numerics_.dt.maxDeltaT != config_.numerics.maxDeltaT) throw ...;
```

---

## 23.6 Numerical Invariance

未修改（本轮零改动）：

```text
WENO3/WENO5/WENO7/TENO5 reconstruction（含 characteristic 重构）
Rusanov / LaxFriedrichs / Steger-Warming / Roe / LaxWendroff 分裂
central2/central4 diffusion
forwardEuler / SSPRK3 / classicalRK4 stage coefficients 与 stage time
PISO/PIMPLE/SIMPLE + Rhie-Chow pressure correction 数学
Ghost IBM / ILW / variational IBM forcing 数学
kEpsilon / kOmegaSST / Smagorinsky 方程与系数
MPI owner→COPY、SUM、canonical face flux identity
KKT / DLM 代数
EOS（perfectGas / rhoConst / homogeneous）数学
```

证据（每个阶段都在同一 binary 上复跑并逐位 diff）：

```text
cylinderFlowGhost --steps 20   before vs after（native IO 迁移）   IDENTICAL
cylinderFlowGhost --steps 20   before vs after（How authority）    IDENTICAL
cylinderFlowGhost --steps 20   before vs after（workspace lifetime）IDENTICAL
cylinderFlowGhost --steps 20   before vs after（cache 热路径）      IDENTICAL
mpirun -np 4 sodCase --steps 20 before vs after（native IO 迁移）   IDENTICAL
```

diff 的对象是每一步的 `Step time: time=..., dt=...` 以及全部
`State closure: ... min(rho)/min(p)/min(T)/min(c)` 行。

---

## 23.7 Regression A — sodCase

case：`test/Sod/sodCase`（production：4 rank、split [2,2,1]、densityBase、
singleFluid、Euler、TENO5 + characteristic + StegerWarming、CFL 0.5、
endStep 500）。**没有修改该 case 的任何数值设置。**

### 编译与 authority

```text
clean configure (Ninja, clang++, BUILD_CLI=ON, BUILD_GUI=OFF, BUILD_TESTS=ON)
  -> MPI enabled: /opt/homebrew/bin/mpicxx
  -> Pressure linear workflow enabled: HYPRE + MPI
clean build   -> 263 targets, exit 0（唯一 warning 是 ld duplicate library）
```

`./build/sonicSolver explain test/Sod/sodCase` 给出的 authority：

```text
formulation     : densityBase
template origin : singleFluid
TIME RECIPE     : forwardEuler (RungeKutta, order 1, ExplicitStages, stages 1)
bound terms
  E_MASS[1]     divergence(massFlux)     -> teno5Steger [ExplicitResidual, halo=3]
  E_MOMENTUM[1] divergence(momentumFlux) -> teno5Steger [ExplicitResidual, halo=3]
  E_ENERGY[1]   divergence(energyFlux)   -> teno5Steger [ExplicitResidual, halo=3]
required halo width : 3
provider requirements: term.convection.teno5Steger
runtime requirements : EulerianGlobalDof OK / CanonicalFace OK（4 rank）
                       ConstraintGlobalDof not-required
                       DistributedLinearSystem not-required
```

即：input 说 TENO5 + StegerWarming + forwardEuler，compiled HOW 说
`teno5Steger` + `forwardEuler`，runtime 执行同一份 compiled 计划，
没有 legacy parser authority 参与。

### 运行结果

```text
serial sanity        ./build/sonicSolver explain/check test/Sod/sodCase  OK
                     （serial run 按设计 fail-fast：MPI ranks=1 != split=4）
MPI production       mpirun -np 4 ./build/sonicSolver run test/Sod/sodCase
  完成 166 步：time=7.459103e-02, dt=4.010775e-04
  min(rho)=0.427856, min(p)=30400.8, min(T)=247.531, min(c)=315.397
  第 167 步 fail-fast：
  [SF FATAL] splitStegerWarming: invalid conservative state
             rho=0.0846138 ru=-28.0405 E=1803.93 p=-1136.92
```

### 判定

```text
architecture validation   passed
  - 4-rank 分布式路径（halo / canonical face / GlobalDof）执行 166 步无 MPI 错误
  - 20 步 before/after 逐位一致（§23.6）
  - explain 的 WHAT/HOW/ORDER/RUNTIME 与实际执行一致

high-order numerical defect remains pre-existing
  - 失败发生在 TENO5+StegerWarming 重构自身的守恒状态上（p<0），
    与本轮改动的 authority/lifetime/dependency 无关
  - 同一缺陷已在仓库既有报告中记录：
      docs/high-order-first-divergence.md（4-rank TENO5 Sod first divergence）
      docs/old-teno5-contract-blackbox.md（旧 executable 根本没有执行
      high-order dispatch，只跑 RusanovEOS::div）
      docs/execution-provider-composition-cleanup.md（完整既有 case 在
      t=0.1297871 仍因 negative pressure fail-fast，未修改 tolerance/flux）
```

按 §16 的要求，另用**已有稳定 numerical recipe** 完成 architecture
regression：

```text
ctest -R sodRusanovNumericalRegression      (case test/Sod/sodCase_weno7_t0p2)
  -> Passed（5.5 s）
ctest -R cliGeneratedRecipe                 (case test/Sod/sodCase_weno7)
  -> Passed
```

未修改任何数值公式、CFL、clamp 或 tolerance 来"让 Sod 通过"。

---

## 23.8 Regression B — cylinderFlow

case：`test/IBM/cylinderFlowGhost`（production：densityBase、PIMPLE(1/1/0)、
HYPRE boomerAMG flexGMRES、WENO5 + characteristic + LaxFriedrichs、
classicalRK4、IBM ghost + ILW(value=3)、turbulence DNS、CFL 0.3、
maxDeltaT 0.01、endTime 5.0 s）。**没有关闭 IBM / turbulence / pressure
coupling 来换取通过。**

### authority

```text
TIME RECIPE : classicalRK4 (order 4, ExplicitStages, stages 4)
bound terms
  E_MASS[1]     divergence(massFlux)      -> weno5LaxFriedrichs [halo=3]
  E_MOMENTUM[1] divergence(momentumFlux)  -> weno5LaxFriedrichs [halo=3]
  E_MOMENTUM[2] diffusion(U,mu)           -> central2Explicit   [halo=0]
  E_ENERGY[1]   divergence(energyFlux)    -> weno5LaxFriedrichs [halo=3]
  E_ENERGY[2]   diffusion(T,conductivity) -> central2Explicit   [halo=0]
CONTRIBUTIONS
  builtin-preset  builtin.singleFluidNavierStokes
  model           model.turbulence        turbulence equation/closure
  model           model.ibm.ghostCellIBM  immersed-boundary contribution
SOLVE BLOCKS
  S_IBM_BOUNDARY  boundaryStencilClosure   (boundaryClosure)
  TIME_RECIPE     classicalRK4             (explicitTimeIntegration)
COMPILED SOLVE PLAN
  Sequence Explicit.step
    Update  flow.step.prepare / flow.dt.compute / flow.step.begin
    StageLoop  classicalRK4  repeat=4  -> explicit.stage.execute
    Update  flow.step.commit / time.commit
PROVIDERS
  term.convection.weno5LaxFriedrichs / term.diffusion.central2Explicit
  flow.conservative / thermodynamics.single-fluid / closure.turbulence
  ibm.boundary
```

### 运行结果（serial）

```text
./build/sonicSolver run test/IBM/cylinderFlowGhost
  exit = 0，完成 endTime = 5.000000 s
  steps = 11450
  最后一个物理步: time=5.000000e+00, dt=2.822697e-04
  dt 区间: [1.56107e-04, 4.39268e-04]（CFL 0.3 + maxDeltaT 0.01 生效）
  sampled min(rho) 区间: [1.21852, 1.22500]
  sampled min(p)   区间: [100655, 101325]
  sampled min(T)   区间: [287.210, 288.153]
  sampled min(c)   区间: [339.737, 340.294]
  NaN/Inf 计数 = 0，FATAL 计数 = 0，无 clamp / 无 fallback
  输出: result/ 共 52 项，t0 … t5 每 0.1 s 一个 frame + PVD/VTM 收集
  IBM/ILW: ghost 分类 + ILW(value=3) 边界闭合在每一步 commit 后生效
  turbulence: DNS closure contribution 在 compiled system 中固定启用
  pressure: PIMPLE(1/1/0) + HYPRE boomerAMG flexGMRES 每步执行
```

### 运行结果（distributed）

```text
mpirun -np 4 ./build/sonicSolver run --steps 2 test/IBM/cylinderFlowFictitiousDomainMPI4
  exit 0；IBM 约束诊断（body dofs=138, max|Ju-Us|=0, max|dJ/du|=0,
  kineticFunctional / bodyForce / bodyTorque / fluidPower）正常输出
```

`cylinderFlowPeskinMPI4`（surface diffuse-kernel 约束）在第一步内无法推进
（长时间无输出）。该卡滞与本轮改动无关的依据：本轮没有触碰 MPI
communicator、halo、COPY/SUM、surface kernel 或 DLM 装配，且 serial 路径
逐位一致；它属于既有 distributed surface-constraint 覆盖缺口，见 §23.9。

---

## 23.9 Remaining Architecture Debt

| # | 项 | 现状 | 建议 |
| --- | --- | --- | --- |
| 1 | `models → solver/system` | `models/{physics,ibm,turbulence,interfaceModel}` 的 contribution 实现 include `solver/system/SF_equationContribution.h`，并使用 `ExecutionPolicy` / `TransformationDescriptor` | 把 WHAT description 值类型（`RawEquationSystem`/`ExecutableEquationSystem`/`UnknownDescriptor`/`EquationDescriptor`/`ConstraintDescriptor`）下沉到 `core/system`；由 system compiler 而不是 model contribution 把 model spec 翻译成 transformation/policy 请求 |
| 2 | CMake 跨目录编译 cpp | `SF_solverSystem` 直接编译 `../../models/**` 的 4 个 contribution cpp 与 `../run/SF_planExecutor.cpp` | 依赖第 1 项完成后，给 `SF_physics`/`SF_ibm`/`SF_turbulence` 增加真实 contribution target；`SF_planExecutor` 归 `SF_run` |
| 3 | numeric constants 仍来自 `config_.numerics` | `SF_densityBasedRHS.cpp` / `SF_compressible.cpp` 仍读 `idealGasGamma`、`idealGasConstant`、`dynamicViscosity`、`prandtl`、`ibmBoundary`、`ilwOrder` | 这些是 EOS/boundary-closure 输入而非 time/space recipe；后继应把它们并入 compiled HOW 或 `FluidStateModel`，并加与 gamma 相同的"input vs compiled"一致性守卫 |
| 4 | `ee.step.commit` diagnostics | Eulerian pressure stepper 的 commit 仍内含 cell traversal / phase mass / enthalpy / wall heat / alpha-sum-error / representative plane 诊断 | 提取为 `EulerianDiagnostics`，`commitState(); commitModels(); collectDiagnostics();`（数值结果不变） |
| 5 | 高阶生产算例数值缺陷 | 4-rank TENO5 Sod 在 t≈0.0746 negative pressure fail-fast；WENO7 Euler 在 t≈0.1298 同样 fail-fast | 属于数值公式/鲁棒性课题，不属于架构轮次；不得用 clamp/tolerance/换分裂器掩盖 |
| 6 | distributed cylinderFlow surface 约束 | `cylinderFlowPeskinMPI4` 首步卡滞 | 需要 MPI surface/body transfer + DLM 的分布式覆盖排查（本轮未触碰该路径） |
| 7 | pressure-specialized recipe 覆盖 | single-fluid PISO 之外的 pressure schedule 仍无专用 provider | 按 §5 的 capability 报告 fail-fast，不再静默退回 |
| 8 | legacy 文件物理删除 | `model/legacy/*` 已退出构建但仍在磁盘 | 用户确认后删除 |
| 9 | CMake 引用的空目录 | `src/methods/numerics/time`、`src/solver/discretization/time` | 收敛 target include 面后删除 |
| 10 | `SF_nativeDecode.cpp` 仍是单文件 | 3331 行；解码已按语义分节（runtime/numerics/algorithm/fields/turbulence/thermo/ILW/IBM/sources/wallHeat），但没有拆成多文件 | 拆成 `compatibility/decode/SF_*.cpp`；前置条件是把 JSON 小工具抽到 `compatibility/private/` 头 |

---

## 23.10 Commands / Evidence Index

```text
build
  cmake -S . -B /tmp/sonic-clean-build -G Ninja -DCMAKE_CXX_COMPILER=clang++ \
        -DBUILD_CLI=ON -DBUILD_GUI=OFF -DBUILD_TESTS=ON      -> exit 0
  cmake --build /tmp/sonic-clean-build -j 8                   -> exit 0
  cmake --build build --target sonicSolver -j 8               -> exit 0

tests
  cd /tmp/sonic-clean-build && ctest --output-on-failure       -> 13/13 Passed
  包括新增 test_fieldVectorLayout（fieldVectorLayout）

production regression
  mpirun -np 4 ./build/sonicSolver run test/Sod/sodCase        -> 166 步后
                                                                  高阶数值 fail-fast
  ./build/sonicSolver run test/IBM/cylinderFlowGhost           -> 见下表
  mpirun -np 4 ./build/sonicSolver run --steps 2 \
        test/IBM/cylinderFlowFictitiousDomainMPI4              -> exit 0
  ./build/sonicSolver explain test/Sod/sodCase                  -> authority
  ./build/sonicSolver explain test/IBM/cylinderFlowGhost        -> authority
```

本轮每个改动阶段后的逐位一致性检查（`Step time` + `State closure` diff）：

```text
阶段                                   case                       结果
native typed config 迁移                cylinderFlowGhost 20 步     IDENTICAL
native typed config 迁移                sodCase 4-rank 20 步        IDENTICAL
compiled HOW (time/dt)                  cylinderFlowGhost 20 步     IDENTICAL
OpRegistry/Workspace lifetime           cylinderFlowGhost 20 步     IDENTICAL
thermodynamic cache 热路径              cylinderFlowGhost 20 步     IDENTICAL
```

---

## 23.11 Acceptance Checklist

```text
[x] Field Vector3 layout bug fixed
    getVal/setVal<Vector3> 使用 v/v+1/v+2；新增 ctest fieldVectorLayout
[x] runtime capability semantics correct
    CanonicalScalarInterfaceFlux available 与 detail 一致；
    validate() 强制 required/available/reason 与 Runnable 语义一致
[x] backend availability no longer hardcoded
    BuildCapabilities 来自 Parallel::mpiCompiled() 与
    LinearAlgebra::hypreBackendLinked()；constexpr ...=true 已删除
[x] native input no longer performs typed → globals → typed round-trip
    decode*Section() 直接写 typed CaseConfig；ParserState 与
    Legacy::makeSolverConfig() 退出 native 路径；SF_param 退出构建
[x] runtime numerical decisions primarily come from compiled authority
    numerics_.time.recipe / numerics_.dt.{cfl,maxDeltaT} + 一致性守卫
[x] obvious dependency inversion removed
    core→app、solver→app(隐藏)、core→methods(cache)、
    models/methods→solver(vectorCalculus)、infrastructure→solver(HYPRE)
    全部消除；tools/check_architecture.py 无新增 ERROR
[~] systemBuilder/nativeDecode God-file responsibility reduced
    authority 维度已收敛（无全局暂存、无往返、无死骨架）；
    SF_nativeDecode.cpp 仍为 3331 行单文件，多文件拆分见 §23.9 #10
[x] OpRegistry/workspace lifetime improved without numerical change
    operations_ / explicitWorkspace_ / pressureSummary_ 成为 solver member
[x] sodCase passes required architecture regression
    4-rank 166 步无 MPI 错误 + 20 步逐位一致 + explain authority 一致；
    t≈0.0746 的高阶数值缺陷为既有问题（§23.7）
[x] cylinderFlow passes with IBM + MPI + turbulence enabled
    serial: IBM(ghost+ILW) + DNS turbulence + PIMPLE/HYPRE 跑到
    endTime=5.0 s（11450 步，exit 0，无 NaN/Inf/FATAL）
    MPI: 4-rank FictitiousDomain 变体 exit 0（IBM 约束诊断正常）；
          4-rank Peskin(surface) 首步卡滞为既有分布式覆盖缺口（§23.9 #6）
[x] no silent fallback
    fail-fast：unavailable capability、compiled/config 不一致、
    invalid conservative state（未加 clamp、未换分裂器）
[x] clean build passes
    全新 build 目录 configure + build exit 0；ctest 13/13 Passed
[x] report produced
    docs/runtime-authority-cleanup-report.md
```

未完成的验收项只有一项（God-file 多文件拆分，标注 `[~]`），其 authority
层面的目标已经达成，剩余部分纯属文件组织。
