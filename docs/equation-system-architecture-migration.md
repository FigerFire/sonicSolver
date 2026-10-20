# Equation System Architecture Migration Report

日期：2026-09-16

## 1. Previous architecture

迁移前 `ResolvedSimulationSystem` 同时保存 unknown、physical equation、algorithmic pressure
equation、constraint、flat solve block、closure、workspace 和 capability。`SF_systemBuilder.cpp`
在一个阶段内完成 composition、pressure/IBM transformation-like expansion 和 solve-order
construction。`Workflow::makePlan` 再把 solve block 一对一转成 flat stage，最终由 density 或
Eulerian specialized stepper 执行自己的固定 lifecycle。

主要问题是：raw mathematical system、transformed executable system 与 execution plan 没有
类型边界；pressure constraint 甚至由 explain 根据 `E_PRESSURE` 名称合成，而不是 raw system
中的正式 constraint。

详细证据见 `docs/system-composition-transformation-audit.md`。

## 2. New architecture

启动期现在形成唯一链：

```text
SystemCompositionBuilder
        ↓
RawEquationSystem
        ↓
TransformerRegistry / TransformationPipeline
        ↓
ExecutableEquationSystem
        ↓
SolvePlanner
        ↓
CompiledSolvePlan
        ↓
transitional legacy execution adapter
```

`ResolvedSimulationSystem` 只聚合这些冻结阶段对象以及 runtime capability/workspace requirement。
Runtime assembly、state realization、validator 和 specialized steppers 已改为读取
`executableSystem`；workflow flat view 来自 `solvePlan.legacyBlocks`。

## 3. Composition layer

新增/演化的 contract：

- `SystemCompositionBuilder`：统一注册 unknown、equation、term extension、constraint、closure、
  dependency、transformation request 和 execution policy；
- `ISystemContribution`：builtin/model/user 共用的 contributor 接口；
- `SystemComposer`：按显式顺序应用 contributor；
- `ContributionRecord` / `Provenance`：记录 `BuiltinDefault`、`BuiltinPreset`、`Model`、`User`、
  `Generated`；
- `SystemModification`：区分 `Add`、`Extend`、`Replace`、`Disable`。

当前 production composition 已统一通过 `SystemCompositionBuilder`。Source term 使用
`extendEquation`，因此 explain 能显示它是显式 extend modification，而不是同名覆盖。

## 4. Built-in presets

Single-fluid NS 现在以 `builtin.navierStokes` contribution 注册 rho/rhoU(or U)/rhoE 与
Mass/Momentum/Energy。Eulerian-Eulerian 以 `preset.eulerianEulerian` 注册每相 unknown、
continuity/momentum/enthalpy 和 raw shared-pressure constraints。

Pressure coupling 会额外记录 `preset.coupling.SIMPLE/PISO/PIMPLE`。Preset 的结果是普通 raw
objects、transformation request 和 execution policy；runtime 不读取 provenance 选择路径。

## 5. Model contributions

Level-set、homogeneous multiphase、legacy mixture、turbulence、source 与 IBM 通过同一个
composition builder 写入 raw system。Model 不创建 timestep loop。

Ghost IBM 只贡献 boundary-closure execution policy；variational IBM 贡献 multiplier、constraint、
constraint equation、transformer request 和 solve policy。两类 IBM 的数学语义没有被合并。

## 6. User overrides

接口层已建立显式 modification vocabulary：

```text
add
extend
replace
disable
```

Typed add 与 equation-term extend 已实现。Replace/disable 尚无完整配置 payload 和 lowering，
当前明确抛出异常；不会使用“同名后写覆盖”或静默忽略。用户 equation DSL/config reader 尚未
接入该接口，列入 implementation pending。

## 7. RawEquationSystem

`RawEquationSystem` 现在正式保存：

- unknowns；
- physical/constraint equations 与 `Equation::Definition`；
- constraints；
- closures/dependencies；
- contributions/modifications/provenance。

它不保存 solve block、RK/PIMPLE loop、MPI schedule 或 runner identity。Validator 明确拒绝 raw
system 中的 `AlgorithmicDerivedEquation`。

## 8. System Transformation layer

新增：

- `IEquationSystemTransformer`；
- `ExecutableEquationSystemBuilder`；
- `TransformerRegistry`；
- `TransformationPipeline`；
- `TransformationDescriptor/Match/Record`；
- `NotRegistered / RegisteredAndActive / RegisteredButNotApplicable /
  RegisteredButInvalid / Applied` 状态。

Pipeline 按 priority 稳定排序。未注册 transformer、显式请求但不适用、或 invalid contract
都会 fail-fast。Transformer 不 include/call MPI。

## 9. SIMPLE/PISO/PIMPLE repositioning

Single-fluid pressure composition 现在先在 raw system 注册：

```text
Momentum physical equation
pPrime unknown
C_INCOMPRESSIBILITY constraint
```

`PressureConstraintTransformer` 只检查该数学 contract，不检查 density/pressure solver identity。
匹配后生成：

```text
E_PRESSURE
OP_PRESSURE_UPDATE
OP_VELOCITY_CORRECTION
OP_FLUX_CORRECTION
```

SIMPLE/PISO/PIMPLE 名称和迭代语义由 execution policy/plan 保留。当前具体 pressure numerical
lifecycle 仍由 `LegacyPressureExecutionAdapter` 执行。

## 10. IBM transformation integration

`ImmersedConstraintTransformer` 注册在同一个 registry，并按 `C_IBM*` mathematical constraint
匹配。Fractional/monolithic IBM 的现有 constraint equations 继续由 IBM descriptor 提供；
transformer 建立统一 architecture identity，numerical projection/KKT lowering 仍由明确报告的
legacy IBM adapter 执行。

Ghost/ILW 没有注册 constraint transformer，仍是 boundary/stencil closure。

## 11. ExecutableEquationSystem

`ExecutableEquationSystem` 是 runtime equation authority。它从 raw system 继承 physical objects，
并接收 transformer 生成的 algorithmic equations 与 correction operators。

迁移后：

```text
CompressibleAlgorithm equations_
PressureStepper AssemblyPlanRegistry
StateRealizer unknown binding
SystemValidator
SystemPrinter
```

都读取 `executableSystem`，不再读取 `ResolvedSimulationSystem` 上的第二套 flat equation fields。

## 12. SolvePlanner

`SolvePlanner` 输入 executable system、execution policies 和 time-integrator identity。它不新增或
替换物理方程。Pressure predictor/correction policy 在 structured view 中合成为一个 loop；
legacy blocks 仍保持原顺序，供现有 stepper compatibility binding 使用。

## 13. CompiledSolvePlan

`CompiledSolvePlan` 支持以下 node kind：

```text
Sequence, Loop, StageLoop, Subcycle, BlockSolve,
Assemble, Solve, Correct, Update, Synchronize,
Reduction, Commit, ConvergenceCheck
```

当前已实际生成：

- density explicit `StageLoop`；
- pressure/Eulerian `Loop`；
- fractional IBM `Sequence(Solve, Correct)`；
- monolithic IBM `BlockSolve`；
- Ghost IBM boundary `Update`，位于 density `StageLoop` 内、`Assemble` 之前。

Subcycle、Reduction、Commit 等 node 已有类型，但 generic runtime lowering 尚未实现。

## 14. Runtime authority

Production 数值执行仍使用现有 algorithms。`CompiledSolvePlan` 明确记录：

- `LegacyDensityExecutionAdapter`；
- `LegacyPressureExecutionAdapter`；
- `LegacyEulerianExecutionAdapter`；
- `LegacyImmersedExecutionAdapter`。

这些 adapter 消费 executable equations 与 compiled legacy blocks，不创建 raw/executable/plan
authority。`PhysicsTemplateKind` 没有进入 application runner dispatch。

## 15. Legacy adapters

保留 adapter 是为了保护已验证的 RK、Eulerian PIMPLE、IBM、WENO/TENO、boundary/halo 和 MPI
lifecycle。它们是 transitional implementation boundary，并已在 explain 中逐项列出。

尚未删除的 specialized class 不再被文档描述为最终 solver family。

## 16. Empty/stub implementations

### Interfaces complete but implementation pending

| 接口/能力 | 当前状态 |
|---|---|
| User field/equation config lowering | Contribution/modification API 已完成；case DSL reader 未实现。 |
| User replace/disable | 语义类型与 fail-fast 已完成；typed payload/lowering 未实现。 |
| Custom transformer config registration | C++ registry/interface 可用；用户配置 loader 未实现。 |
| SIMPLE/PISO/PIMPLE generic execution | Transformation 与 structured policy 已完成；数值执行仍用 legacy pressure adapter。 |
| IBM generic transformation lowering | Registry/match/provenance 已完成；projection/KKT kernel 仍由 legacy IBM adapter 执行。 |
| Generic structured-plan executor | Plan IR 已完成；runtime 尚未逐 node 执行。 |
| Dependency-ordered transformer pipeline | 固定 priority 已实现；general dependency ordering 未实现。 |
| Complete `CompiledEquation` binding | 现有 `AssemblyPlan` 继续使用；component slice、全部 storage/sync binding 仍需完善。 |

所有 pending 路径均没有 silent fallback。现有 production configuration 继续进入已验证 adapter；
新的不完整 user modification/transformer request 会 fail-fast。

## 17. Numerical behavior preserved

本轮没有修改 numerical formula、RK coefficient、stage time、boundary/halo order、convection
dispatch、residual sign、source clear timing、IBM correction order 或 storage layout。

验证结果：

- 完整 Debug build：通过，303 targets；无 compile/link failure。Ninja 每次先报告
  `premature end of file; recovering`，CMake/Ninja 自动恢复后完整成功；
- density thermal production case：完整运行到 `t=0.02`，通过；
- Eulerian-Eulerian PIMPLE/HYPRE production case：完整运行到 `t=2e-4`，通过；
- Ghost IBM：完成 setup 和多步 RK4 smoke；每个 stage 的 state closure 通过；因原 case
  `endTime=5` 很长，smoke 后人工终止，没有 assertion/fail-fast；
- 4-rank TENO5/Steger-Warming Sod：复现既有记录中的同一已知 failure，写出 step 100 后在
  `splitStegerWarming` 对 `rho=0.104189, p=-872.514` fail-fast。该数值问题已记录于
  `docs/high-order-first-divergence.md`，不是本轮首次出现，也没有被 fallback 掩盖。

CTest：host MPI 下 `registryIO`、`distributedConstraintLayout`、`sparseCanonicalCopy`、
`distributedMinimalKKT`、`distributedSurfaceSchur` 五个现有 executable 全部通过。
`numericalFluxContract` 无法运行，因为工作区缺少 `test/test_numericalFlux.cpp` 和对应 executable；
测试 build 重新生成还同时报告四个其他缺失 test source。这是当前 repository export/test-tree
完整性问题，不是测试断言失败。

## 18. MPI invariants preserved

本轮没有修改 MPI/infrastructure 文件和 numerical synchronization code：

```text
state/geometry owner COPY     unchanged
canonical face F* COPY        unchanged
residual/source/load SUM      unchanged
constraint canonical owner    unchanged
```

Host MPI 的 sparse canonical copy、minimal KKT 和 surface Schur tests 全部通过。Transformer
source 中没有 MPI include/call，plan 只声明 Synchronize/Reduction node semantics。

## 19. Remaining implementation work

1. 将 user fields/equations/overrides 从 case model lowering 到 `SystemCompositionBuilder`。
2. 为 replace/disable 增加 typed payload，并继续拒绝名称覆盖。
3. 让 pressure transformer 真正生成完整 predictor/correction compiled bindings。
4. 让 IBM transformer 生成 KKT/fractional executable block representation，而不依赖 descriptor
   已经展开的 equations。
5. 建立 generic plan-node executor，并逐个替换 legacy adapters。
6. 完成 `CompiledEquation` 的 storage、component slice、read/write 与 synchronization binding。
7. 修复/恢复缺失测试源和 `test_numericalFlux` executable。
8. 单独处理已经记录的 high-order TENO5/Steger-Warming numerical divergence。

## 20. Recommended next step

下一步应选择一条可独立验证的 transformation 做完整 vertical lowering。推荐 pressure-constraint
PISO 小 case：

```text
Raw constraint
  -> generated pressure equation/operators
  -> compiled bindings
  -> structured Loop execution
  -> remove LegacyPressureExecutionAdapter for that capability
```

该阶段必须建立专门的 small regression，比较 predictor、pressure residual、velocity/flux correction
和 commit state。不要同时迁移 IBM KKT 或 high-order convection。
