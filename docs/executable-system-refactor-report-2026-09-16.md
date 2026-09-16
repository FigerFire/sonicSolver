# Executable System 重构阶段报告（2026-09-16）

## 1. 范围与结论

本轮严格保持数值公式、RK 系数与 stage time、boundary/halo 顺序、IBM
校正位置、canonical face COPY 及 GlobalDof SUM 语义。完成了 executable
architecture 的安全前置部分：resolved unknown storage contract、完整 term
vocabulary、全 equation AssemblyPlan 建表、Eulerian equation/plan 绑定，以及
Eulerian solver workspace 的启动期验证和显式 borrowed synchronization。

本轮没有宣称完成最终 `SystemExecutor` 架构。审计发现当前
`S_TURBULENCE` 在 resolved workflow 中是顶层 block，但真实 Eulerian 数值生命
周期把 implicit turbulence solve 放在每个 PIMPLE outer corrector 的 energy solve
之后；density auxiliary equation 又嵌在每个 explicit RK stage 中。当前
`SolveBlock` 不表达这种 nested placement。直接让一个 generic executor 顺序遍历
现有 blocks 会改变时间离散与耦合顺序，因此按任务的“若需同时改变数值行为则停
止”规则，本轮停在该边界，没有用 forwarding `SystemExecutor` 假装完成 Phase G–X。

## 2. Audit

详细审计见 `docs/executable-system-audit-2026-09-16.md`。主要 hidden authority：

- density mass/momentum/energy 已从 `ResolvedSimulationSystem::equationDefinitions`
  建立 plan；homogeneous、legacy alpha、level-set 仍通过 lifecycle coupling 执行；
- Eulerian continuity/momentum/enthalpy 的方程身份已 resolved，但原 assembler 仍
  隐含 operator 集合；
- turbulence equation 已声明，执行位置仍由 `PressureStepper` 内部生命周期决定；
- IBM descriptor 声明 unknown/constraint/solve block，但具体 KKT/forcing backend
  仍由 specialized IBM strategy 持有；
- output/diagnostics 仍在 application runner 中按物理类型枚举。

## 3. Phase A — Runtime State Contract

`UnknownDescriptor` 新增：value shape、primary/transported/algebraic/multiplier/
derived role、storage binding、storage key/component offset、初始化/边界/restart/
output eligibility、runtime-storage requirement 与 namespace。

新增 `StateRealization`。`runFlow` 在进入 `Time::Driver` 前调用 executor
`prepare()`；density 和 Eulerian executor 在该阶段把 resolved unknowns 绑定到
`StateBundle::distributed` 的稳定 view，并校验 patch 数量与 component layout。
数值数组仍由 Field、phase state、turbulence state 或 IBM specialized storage
持有；没有复制 Field，也没有创建第二份 physical state。

IBM multiplier 与尚未迁移完的 homogeneous/legacy auxiliary storage 被明确标为
`SpecializedExecutor`，而不是伪装成已统一的 named storage。这是保留的 transitional
binding，报告中不把它计作最终 state realization 完成。

## 4. Phase B/C — Term 与 AssemblyPlan Authority

`Equation::Term` 现在可表达 gradient、algebraic relation、explicit/implicit/
algebraic evaluation mode、time level、closure provider、boundary requirement 和
linearization requirement。它仍是有限 operator vocabulary，没有引入 symbolic CAS。

新增启动期 `AssemblyPlanRegistry`：每个 `EquationDescriptor` 必须从同一个
`equationDefinitions` 建立稳定 plan；validator 同时要求每条 equation 至少属于一
个 solve block。`E_PHASE_MASS`、`E_LEGACY_ALPHA`、`E_LEVEL_SET` 已加入现有 primary
fluid block，消除了 resolved equation 无 block owner 的静默状态。

## 5. Phase D — Eulerian Lowering

`PressureStepper` 在构造时为全部 resolved equations 建立 registry，并把同一 registry
绑定给 `PhaseEquationAssembler`。Assembler 根据 phase equation ID 取得 continuity、
momentum、enthalpy 与 shared-pressure plan，并在执行对应 kernel 前 fail fast 验证所需
term 存在。

Eulerian authoritative definitions 已补全真实执行中的 pressure gradient、momentum
diffusion 与 enthalpy diffusion。因此 `explain` 展示的 operator 与 runtime 请求一致，
但具体 sparse assembly kernel 与原公式完全未改。

## 6. Phase L — Workspace Requirement 与 Canonical Face

`ResolvedSimulationSystem` 新增显式 `WorkspaceRequirement`。Eulerian composition 声明：

- 每相 `momentumDiagonal`；
- 每相 `volumeFaceFlux` 与 `massFaceFlux`；
- `pressureCorrection`。

`StateRealizer` 在 timestep 前确认这些 solver-owned workspace 已注册并且 component
layout 一致。旧路径把 canonical-face write 交给不带 storage 的 `finalize()`，运行时
必然抛出 `workspace synchronization must be invoked with explicit borrowed storage`。
新路径借用已注册的 phase face-workspace view 并调用
`ExecutionRuntime::synchronizeTransient()`。它继续执行 canonical owner → COPY；没有
第二个 face flux、没有 average，也没有 full-field copy。

## 7. Explain 与 Architecture Guard

`sonicSolver explain` 现在显示 unknown role、storage binding/key/offset，以及 solver
workspace requirements。`tools/check_architecture.py` 新增检查：

- `runFlow` 必须在时间循环前调用 state preparation；
- executor interface 必须提供 prepare hook；
- validator 必须为全部 equations 建立 `AssemblyPlanRegistry`；
- state realizer 必须验证 workspace requirements；
- Eulerian face workspace 不得退回隐式 `writeCanonicalFace` finalize。

依赖 allowlist 仍为 10，没有新增 dependency edge。

## 8. Ownership / Execution Order / API

| 项目 | Before | After |
|---|---|---|
| resolved unknown storage | description only；runtime module 各自注册 | description 带 binding contract；run 前统一 realize/validate |
| equation lowering | density 三方程局部 vector；Eulerian hidden | 全 equation registry；Eulerian kernel 绑定 authoritative plans |
| auxiliary equation block ownership | 3 条 equation 可无 block owner | validator 要求至少一个 owner |
| Eulerian face workspace | solver owned，但 sync 调用没有 view | solver owned，runtime 显式借用已注册 view |
| public transient API | bind services/stages → advance | bind services/stages → prepare state → advance |

数值执行顺序保持：

```text
Eulerian:
boundary/halo
→ interphase/source/turbulence preparation
→ momentum diagonal + face flux
→ canonical COPY
→ continuity
→ momentum predictor/interphase
→ canonical COPY
→ pressure/non-orthogonal correction
→ energy
→ turbulence
→ boundary/halo
→ commit
```

## 9. Verification

| 检查 | 结果 |
|---|---|
| `cmake --build build --parallel 4` | Passed，307 targets；Ninja 每次报告 metadata `premature end of file; recovering` 后自动 regenerate，未出现 compile/link failure |
| `python3 tools/check_architecture.py` | Passed；615 source files，allowlist=10 |
| single WENO7 + Steger-Warming，一步 | Passed；`dt=6.681531e-04`，`rho_min=0.125`，`p_min=10000`，与既有基线一致 |
| serial Ghost IBM WENO5/LF/RK4，一步 | Passed；`dt=4.392680e-04`，`rho_min=1.22383`，`p_min=101208` |
| 4-rank TENO5/Steger-Warming Sod，一步 | Passed；`dt=6.681531e-04`；canonical multi-patch path 正常 |
| Eulerian representative 完整 case | Passed 至 `2e-4`；末步 pressure residual `3.413064e-16`，phase mass/enthalpy 保持既有量级 |
| RPI wall boiling 完整 case | Passed 至 `2e-4`；不再出现 borrowed-workspace failure |
| RPI 一步基线 | `dt=1e-5`，`pIter=5`，`pResidual=2.660242e-10`，`phaseMass=3.213304e+03`，`phaseEnthalpy=1.062648e+09` |
| variational DFM IBM，一步 | Passed；`max|Ju-Us|=0`，`dt=4.392680e-04` |
| Architecture CTest binaries | 5/5 existing binaries passed，包括 MPI/HYPRE tests |

完整 CTest 不能重建。当前工作树缺少 CMake 引用的
`test_registryIO.cpp`、`test_distributedConstraintLayout.cpp`、
`test_sparseCanonicalCopy.cpp`、`test_distributedKKT.cpp` 和
`test_distributedSurfaceSchur.cpp`。旧 build tree 的五个现存 executable 全部通过；
`numericalFluxContract` executable 不存在，因此 CTest 将它标为 **Not Run**。没有凭空
重写这些测试。

完整 serial WENO7 在 `t=0.1297871` 复现已记录的真实 high-order
Steger-Warming invalid-state fail-fast；一步基线通过。本轮没有 clamp、降阶、自动换
flux 或修改 CFL。

`mixtureCase` 继续按既有 capability validation 拒绝 generic six-variable EOS 与
characteristic WENO3/Steger-Warming 组合；没有退回 Rusanov。

## 10. 尚未完成与停止边界

以下工作未在本轮伪装为完成：

- Phase E/F：把 homogeneous/level-set/phase-change/IBM specialized lifecycle 全部
  lowering 为 contribution/executor；
- Phase G–K：`SystemExecutor`、`SolveBlockExecutorRegistry`、PIMPLE nested block 与
  turbulence generic executor；
- Phase M/N：boundary freshness 与 per-operator numerics handle；
- Phase O：generic characteristic thermodynamics capability；
- Phase P：output/diagnostic registry；
- Phase Q/R/W：删除 physics-specific runner、transitional adapters 与
  `INavierStokesStepper` 命名；
- Phase S/T/U：PassiveScalar、constraint、source 三个 end-to-end acceptance；
- Phase X：最终 central-code physics-name-free 验收。

继续前必须先把现有 nested numerical placement 编码为 resolved execution metadata，
并用当前 PIMPLE/RK checkpoint 固定顺序；否则 generic stage iteration 会改变数值生命
周期。该问题应作为下一安全子阶段处理，不能用只转发旧 stepper 的 Manager/Adapter
绕过。
