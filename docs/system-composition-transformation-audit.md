# Equation System Composition / Transformation 架构审计

日期：2026-09-16

## 1. 审计范围与结论

本审计覆盖 `CaseConfig` 到 runtime stepper 的完整启动链。当前代码已经建立了
`Equation::Definition`、`ResolvedSimulationSystem`、类型化 `SolveStrategyKind` 与
`AssemblyPlan`，但 `ResolvedSimulationSystem` 仍同时承担数学组合结果、算法派生方程、
执行顺序和 backend capability 四类职责。它因而不是严格意义上的 raw system，也不是
严格意义上的 executable system。

当前 production 数值链可以保留，但启动阶段必须明确拆成：

```text
SystemContribution
  -> RawEquationSystem
  -> TransformationPipeline
  -> ExecutableEquationSystem
  -> SolvePlanner
  -> CompiledSolvePlan
  -> transitional legacy executor
```

## 2. 当前 equation authority

| Authority | 文件 | 当前关系 | 问题 |
|---|---|---|---|
| 符号方程 | `src/solver/equation/SF_expression.h` | `Equation::System` 保存 `Equation::Definition` | 可作为 composition/transformation 共用的方程定义容器。 |
| 方程描述 | `src/solver/system/SF_resolvedSimulationSystem.h` | `ResolvedSimulationSystem::equations` 与 `equationDefinitions` 必须按 id 保持一致 | descriptor 与 definition 是同一 authority 的两个互补视图，但当前没有阶段边界。 |
| 执行 lowering | `src/solver/equation/SF_assemblyPlan.*` | 由 resolved definitions 编译 term/operator binding | 当前只编译 term，尚未编译 structured control flow。 |
| specialized execution | `src/solver/algorithm/SF_compressible.cpp`、`src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp` | 复制 equation system，按 solve block 校验并运行固定 lifecycle | 仍是 production executor；必须明确标记为 transitional adapter，不能成为新的 composition authority。 |

## 3. 当前 model authority

`CaseConfig` 和 `BuildRequest` 决定启用哪些物理模块；
`src/solver/system/SF_systemBuilder.cpp` 中的自由函数直接向最终 resolved object 写入未知量、
方程、约束、closure、workspace 与 solve block。source、turbulence、level-set、multiphase
和 IBM 因而没有统一的 contribution provenance，也没有统一的 add/extend/replace/disable
修改语义。

IBM 的 `ImmersedAlgorithmDescriptor` 已经接近 contributor input，但 builder 同时把其内容
解释成方程、约束和 solve block。Ghost/ILW 是 boundary closure；variational constraint
才应进入 transformation infrastructure。这两个数学角色必须继续区分。

## 4. 当前 transformation-like logic

| 位置 | 当前逻辑 | 实际架构含义 |
|---|---|---|
| `addPressureBasedFluid` | composition 时直接生成 `E_PRESSURE`、`S_PREDICTOR`、`S_PRESSURE` | 把 pressure-constraint transformation 和 solve planning 混入 builtin NS composition。 |
| `addSharedPressureConstraint` / `addEulerianEulerianSolveBlock` | 直接生成 `E_SHARED_PRESSURE` 与 PIMPLE block | shared-pressure constraint、derived equation 和 PIMPLE schedule 混在一起。 |
| `addImmersed` | descriptor 直接生成 unknown/equation/constraint/block；monolithic 时修改 block 并删除 `S_PRESSURE` | IBM physical constraint、KKT transformation 与 execution policy 混在一起。 |
| system printer | 根据 `E_PRESSURE`/`E_SHARED_PRESSURE` 合成 `C_INCOMPRESSIBILITY` 文本 | constraint 只存在于 explain 推断中，不是 raw system authority。 |

## 5. 当前 solve-order authority

`ResolvedSimulationSystem::solveBlocks` 是 flat order。`Workflow::makePlan` 将其一对一映射成
`SolveStage`，随后 specialized stepper 再按自身固定 lifecycle 执行。它不能表达 loop、
stage loop、subcycle、block solve、synchronization、reduction、commit 或 convergence check
之间的层级关系。

因此当前存在两个互相校验但语义层次不同的执行 authority：flat solve block 和 specialized
stepper control flow。迁移后 `CompiledSolvePlan` 是 explain/validation 的计划 authority；在
generic executor 完成前，它显式记录 legacy adapter，production lifecycle 仍由该 adapter
执行。

## 6. 当前 runner dispatch

应用层已主要使用 resolved equation/strategy 查询，而非 physics template id 选择 density
生命周期；但仍保留两类 specialized executor：density explicit algorithm 与 Eulerian
pressure stepper。它们是实现能力边界，不应反向定义 mathematical system。

禁止在新层引入：

- 根据 `PhysicsTemplateKind` 选择 timestep lifecycle；
- 根据 `SolverAlgorithm::DensityBased/PressureBased` 判断 SIMPLE/PISO/PIMPLE applicability；
- 未匹配 transformation 时静默使用另一种算法。

## 7. 当前 preset behavior

Single-fluid NS、Eulerian-Eulerian、pressure workflow 与 IBM descriptor 当前由 builder 分支
直接展开到最终对象。它们行为上是 preset，但没有普通 contribution、transformation 和
execution-policy 表示，用户无法用同一种对象显式构造等价系统。

迁移目标是让 built-in preset 只注册普通 architecture objects；provenance 可用于 explain、
override 和 debug，runtime 不按 provenance 分派。

## 8. 重复或隐式 authority

1. `config.numerics.solver`、`ResolvedSimulationSystem::formulation` 与 specialized stepper
   都能暗示 execution family。
2. `equations`/`equationDefinitions` 同时包含 physical 与 algorithmic-derived equation，
   没有 raw/executable 边界。
3. pressure constraint 在 printer 中被合成，但不属于 raw system。
4. `solveBlocks` 描述顺序，stepper 实现真实 loop/stage order。
5. IBM descriptor 与 system builder 都解释 IBM solve semantics。

## 9. 本轮迁移边界

本轮建立并接入：

- unified contribution + explicit modification pipeline；
- `RawEquationSystem` 与 provenance；
- transformer registry、match/status/validation/application；
- `ExecutableEquationSystem`；
- structured `CompiledSolvePlan`；
- layered validator 与完整 explain；
- production legacy adapter 的显式标记。

本轮不修改现有 RK 系数、stage time、boundary/halo/IBM 顺序、WENO/TENO/Steger-Warming、
canonical face COPY、GlobalDof SUM、residual 符号或 storage layout。generic execution lowering、
完整用户 DSL、通用 pressure/KKT executor 可以保持 unsupported，但必须在 capability validation
或 explain 中显式可见。

