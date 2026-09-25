# Executable Equation Authority 重构报告

日期：2026-09-16

## 1. 目标与范围

本轮把现有 equation-centric description 推进为 runtime authority，迁移链为：

```text
CaseConfig
  -> ResolvedSimulationSystem::equationDefinitions
  -> Workflow::Plan::stages
  -> AssemblyPlan
  -> specialized equation/solve-block executor
```

本轮属于“改变怎么组织和执行 equation”的架构迁移。没有修改物理项、离散公式、时间系数、PIMPLE 校正公式、IBM/KKT 数学或并行归约规则。既有 specialized executor 继续承担算法 lowering；resolved system 负责声明“解什么”。

修改前的 authority 图见 `docs/executable-equation-authority-audit.md`。

## 2. Equation ownership：before / after

### Before

| 对象 | description authority | runtime authority |
|---|---|---|
| density mass/momentum/energy | `ResolvedSimulationSystem::equations` | `Equation::Compressible::System` 内部重新创建 definition |
| Eulerian continuity/momentum/enthalpy | resolved descriptor | `PhaseEquationAssembler` 固定成员函数 |
| turbulence transport | 无 resolved equation | `Turbulence::EquationSystem` 根据 model 隐式创建 |
| source term | 配置枚举 | density/Eulerian 各自 central `switch(SourceKind)` |

### After

`ResolvedSimulationSystem::equationDefinitions` 是 executable equation 的唯一 definition authority。每个 `EquationDescriptor` 必须存在同 ID 的 `Equation::Definition`，validator 在 composition 阶段 fail fast。

`Equation::Compressible::System` 保留为 transitional assembly adapter，但只读取 resolved definitions 并创建 non-owning `AssemblyPlan`；它不再创建 mass/momentum/energy definition。`sonicSolver explain` 也直接读取同一组 definition，因此 explain 和 density runtime 不再维护两套 equation 文本。

Eulerian specialized assembler 仍保留现有数值实现，但 continuity、momentum、enthalpy、shared-pressure 以及 turbulence transport equation 的 identity 和 term definition 已进入 resolved authority。完整 term-by-term lowering 仍是后续工作。

## 3. Unified equation contribution contract

新增窄接口：

```text
EquationSystemBuilder
  addUnknown
  addEquation
  addTerm
  addConstraint
  addSolveBlock
  addClosure
  require

IEquationContribution::contribute(EquationSystemBuilder&)
```

该接口只表达数学贡献，不提供 timestep hook，不拥有 state、workspace、MPI barrier 或时间循环。它直接演化现有 `ResolvedSimulationSystem`，没有创建第二套 Equation IR、Context、Manager 或 ServiceLocator。

## 4. Source contributions

density `SourceTerm::Sp` 和 Eulerian source composition 的 central `switch(SourceKind)` 已删除，改为显式 registration table。Gravity、MRF、WallHeat 仍按 `config.sources.enabled` 原顺序执行，原 source kernel 和 application order 不变。

resolved equation definition 同时记录：

- gravity / MRF -> momentum equation；
- wallHeat -> energy 或各相 enthalpy equation。

不支持的 `SourceKind` 会 fail fast，不会静默忽略或替换。

## 5. Turbulence proves “module adds equations”

transported RAS model 现在注册真实 unknown、equation、closure 和 solve block：

| Model | Unknowns | Equations | Closure | Solve block |
|---|---|---|---|---|
| k-epsilon | `k`, `epsilon`（Eulerian 路径按启用 phase 后缀） | `E_TURB_K`, `E_TURB_EPSILON` | `mu_t` | `S_TURBULENCE` |
| k-omega SST | `k`, `omega`（Eulerian 路径按启用 phase 后缀） | `E_TURB_K`, `E_TURB_OMEGA` | `mu_t` | `S_TURBULENCE` |
| DNS / Smagorinsky | 不新增 transport equation | 无 | model closure | 无 |

Eulerian runtime 校验 `S_TURBULENCE` 与实际 transport-equation capability 一致，并由绑定后的 workflow stage 控制原有 `solveTurbulence(dt)`。RPI case 的 explain 显示 liquid phase 的 k/epsilon equations，不再错误地为未启用 turbulence 的 vapor phase 注册 equation。

## 6. AssemblyPlan

新增 `Equation::AssemblyPlan`。它只保存 authoritative `Equation::Definition` 及其 `Term` 的 non-owning reference，不复制 physics metadata。

density mass、momentum、energy adapter 通过 AssemblyPlan 查询 convection、diffusion 和 source term，再调用原有 WENO/TENO、flux splitter、CENTRAL 和 source backend。没有修改任何 kernel、索引、符号或 clear timing。

## 7. Workflow::Plan runtime execution

`Workflow::Stage` 与 runtime 使用同一个 `FDM::SolveStage` 值对象。`runFlow` 在创建 `Time::Driver` 前执行一次：

```text
stepper.bindSolveStages(workflow.stages)
```

两个 stepper 都逐项验证 stage 的顺序、ID、equations 和 constraints 与 resolved solve blocks 完全一致。绑定后不能在 timestep 中替换。

- density/pressure adapter 根据已绑定的 `S_FLUID` 或 `S_PREDICTOR + S_PRESSURE` 选择现有执行路径；
- Eulerian stepper 要求 `S_EE_PIMPLE`，并由已绑定的 `S_TURBULENCE` 决定是否执行 turbulence transport。

`S_TURBULENCE` 仍按原数值 contract 嵌在每个 PIMPLE outer corrector 中。没有把平面 stage 列表机械执行成 “PIMPLE 完成后再解 turbulence”，因为那会改变现有 coupling mathematics。这里的 specialized stepper 是 solve-block lowering executor，而不是第二份 workflow authority。

## 8. Public API：before / after

| API | Before | After |
|---|---|---|
| `CompressibleAlgorithm` construction | config + independently-owned runtime equation definitions | config + resolved system；adapter 引用 authoritative definitions |
| `PressureStepper` construction | config/model 驱动隐式 equation selection | config + resolved system，并验证 turbulence solve block |
| `INavierStokesStepper` lifecycle | `bindServices -> advance` | `bindServices + bindSolveStages -> advance` |
| `runFlow` | 只驱动 time，未消费 plan | 稳定绑定 plan 后驱动 time |
| explain equation source | descriptors / 独立打印 representation | authoritative executable definitions |

## 9. 量化指标

| 指标 | Before | After |
|---|---:|---:|
| density mass/momentum/energy definition authorities | 2 | 1 |
| resolved descriptors 缺少同 ID executable definition | 可存在 | 0（validator 强制） |
| central runtime `SourceKind` switches | 2 | 0 |
| transported turbulence equations 注册到 resolved system | 0 | 每个启用 phase 2 个 |
| density equation AssemblyPlan | 0 | 3 |
| `Workflow::Plan` runtime bind points | 0 | 1（统一 `runFlow`） |
| specialized stepper stage validators/consumers | 0 | 2 |
| timestep 内 workflow/service rebind | 0 | 0 |
| 新 Equation IR | 0 | 0 |
| 新 physical-state authority | 0 | 0 |

## 10. Ownership 与 execution order

### Ownership

- physical state 仍由 `StateBundle` 持有；
- flux/residual/RK workspace 仍由 solver execution lifetime 持有；
- equation definition 由 `ResolvedSimulationSystem` 持有；
- AssemblyPlan 和 compressible adapter 只借用 definition；
- runtime source/turbulence objects 仍持有自己的 numerical workspace，不再拥有 equation identity authority；
- workflow stage 由 `Workflow::Plan` 持有，stepper 仅稳定绑定并消费。

### Execution order

以下顺序均未改变：

- density physical BC -> halo -> Ghost/ILW -> RHS；
- RK/Euler stage time、stage coefficient 和 publication；
- candidate flux -> canonical F* -> COPY -> `+F* / -F*` residual；
- local GlobalDof contribution -> SUM -> owner；
- PIMPLE outer/correction/non-orthogonal loops；
- Eulerian phase energy 后、outer corrector 尾部的 turbulence solve；
- state commit -> derived refresh -> diagnostics/output。

## 11. Dependency 与并行语义

架构检查结果：

```text
source files                    602
allowlisted dependency edges    10
new dependency violations        0
duplicate header basenames      17（既有基线）
```

未引入 core -> solver、methods -> solver 或新的 infrastructure -> solver edge。Canonical face 仍是唯一 F* 后 COPY；GlobalDof 仍为 SUM；state/geometry/labels 仍为 owner COPY。没有新增 raw MPI 调用。

## 12. Build 与 regression

| 验证 | 结果 |
|---|---|
| `cmake --build build --parallel 4` | Passed，306/306 targets；Ninja 报既有 metadata recovery warning 后完整成功 |
| `python3 tools/check_architecture.py` | Passed，allowlist=10 |
| serial WENO7 Sod，1 step | Passed；`dt=6.681531e-04`，`rho_min=0.125`，`p_min=10000` |
| serial Ghost IBM，WENO5/LF/RK4，1 step | Passed；`dt=4.392680e-04`，`rho_min=1.22383`，`p_min=101208` |
| 4-rank TENO5 Sod，1 step | Passed；`dt=6.681531e-04`，canonical multi-patch path 正常 |
| RPI wall-boiling `check/explain` | Passed；`S_EE_PIMPLE -> S_TURBULENCE`，liquid k/epsilon equations 可见 |
| host CTest 可用 targets | 5/5 passed，包括 3 个 MPI/HYPRE tests |

完整 `build-phase1a-tests` 无法重新生成，因为当前工作树缺少 CMake 已引用的五个 test source：`test_registryIO.cpp`、`test_distributedConstraintLayout.cpp`、`test_sparseCanonicalCopy.cpp`、`test_distributedKKT.cpp`、`test_distributedSurfaceSchur.cpp`。因此新增的 `numericalFluxContract` executable 也无法由该 build tree 生成；这是测试源树完整性问题，不是本轮 compile/assertion failure。现存 5 个可执行测试在 host 环境全部通过。RPI 一步执行仍在既有 `workspace synchronization must be invoked with explicit borrowed storage` fail-fast 处停止，本轮没有绕过或修改该基线问题。

## 13. Numerical semantics self-review

| 检查项 | 结论 |
|---|---|
| WENO/TENO/SW/Rusanov/LF formula changed | No |
| RK/Euler coefficient or stage time changed | No |
| pressure/PIMPLE mathematics changed | No |
| boundary/halo order changed | No |
| source application order changed | No |
| residual sign/clear timing changed | No |
| canonical COPY semantics changed | No |
| GlobalDof SUM semantics changed | No |
| IBM/KKT mathematics changed | No |
| hidden fallback introduced | No |

## 14. Deferred issues

- Eulerian `PhaseEquationAssembler` 的 term-level AssemblyPlan lowering；
- homogeneous phase mass、legacy alpha、level-set 与 phase-change contribution 的完整 executable binding；
- IBM boundary/forcing/KKT solve block 的统一 lowering；
- pressure/density specialized executor 对更多 solve-block 组合的 capability；
- high-order WENO/TENO + characteristic flux 与 generic thermodynamics interface consolidation；
- pressure-based workspace borrowed-storage synchronization 基线问题；
- 当前缺失的 CTest source tree 恢复；
- 10 条 allowlisted dependency debt 与 17 个 duplicate-header baseline。

这些项目均未在本轮顺手重构。
