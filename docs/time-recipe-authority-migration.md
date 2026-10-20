# Time Recipe Authority Migration Report

日期：2026-09-20

## 1. Scope

本轮只迁移时间离散 authority，并将 `Equation::Term` 收敛为纯数学描述。修改属于“改变怎么解 equation”的架构迁移；没有增加物理项或 governing equation，也没有修改时间推进公式、空间离散、边界、IBM、PISO 或 MPI 数值语义。

状态：**Implemented**

## 2. Previous time authority

迁移前，时间方法同时存在于输入字符串、`TimeScheme`、`explicitStageCount()`、physics composition 生成的 `ExplicitStages` policy、SolvePlanner 的字符串解析、`Time::Explicit` 以及旧 `ddtDispatch`/Euler/RK4 文件链。stage 数和 lifecycle 因而有多个可写或可推导来源。

迁移后，输入字符串只在配置边界解析一次；`ResolvedSimulationSystem::timeRecipe` 是 resolved authority。Planner 只消费 recipe 的 topology 和 stage count，provider 只消费 recipe ID 执行局部公式。

状态：**Implemented**

## 3. Equation::Term before/after

| 项目 | Before | After |
|---|---|---|
| 数学类别 | transient/divergence/gradient/diffusion/source/constraint | 保持数学类别 |
| source 类别 | `ExplicitSource`、`ImplicitSource` | `Source` |
| execution metadata | `EvaluationMode`、`TimeLevel` | 删除 |
| lowering hints | closure/boundary/linearization requirement | 删除 |
| 职责 | 数学描述和执行 hint 混合 | 只描述 equation 中的数学项 |

状态：**Implemented**

## 4. Removed EvaluationMode / TimeLevel status

`EvaluationMode`、Equation AST 内的 `TimeLevel` 及其字段已删除。仓库中与 multiphase physical old-time storage 有关的同名方法不属于 Equation AST，本轮未改变其物理状态语义。

状态：**Implemented**

## 5. ExplicitSource/ImplicitSource migration

`TermKind::ExplicitSource` 与 `TermKind::ImplicitSource` 已合并为 `TermKind::Source`。`implicitSource()` 已删除，callers 和 printer 均迁移到 `source()`/`Source`。Equation definition 不再选择 source 的显式或隐式 treatment。

状态：**Implemented**

## 6. Physics/time decoupling

Density、pressure、constant-density 和 Eulerian-Eulerian 的 physics composition 只创建 unknown、equation、term、constraint 与 closure。它们不再接收 time name，不再根据 Euler/SSPRK3/RK4 创建 stage policy。相同 physics 分别使用三个 built-in recipe 时，RawEquationSystem signature 完全相同。

状态：**Implemented**

## 7. SystemBuilder before/after

| 项目 | Before | After |
|---|---|---|
| physics helper 参数 | 携带 `timeIntegrator`/time name | 无 time 参数 |
| explicit schedule | physics helper 创建 | SolvePlanner 从 typed recipe lowering |
| resolved system | 保存 time string | 保存 `TimeRecipe` |
| mathematical system | 可能随 time policy 改变 | 与 time recipe 无关 |

状态：**Implemented**

## 8. TimeRecipe design

`TimeRecipe` 是内建、不可由用户拆分配置的值对象，包含：

- typed `TimeRecipeId`
- `TimeFamily`
- `TimeTopology`
- order
- stage count

用户只选择完整 recipe 名称，不能单独设置 `implicit=true`、stage count 或 Butcher tableau。无效名称直接报错，不会降级到另一个算法。

状态：**Implemented**

## 9. TimeRecipeId list

| 输入名 | ID | Family | Topology | Order | Stages | Status |
|---|---|---|---|---:|---:|---|
| `forwardEuler` | `ForwardEuler` | RungeKutta | ExplicitStages | 1 | 1 | Implemented |
| `SSPRK3` | `SSPRK3` | RungeKutta | ExplicitStages | 3 | 3 | Implemented |
| `classicalRK4` | `ClassicalRK4` | RungeKutta | ExplicitStages | 4 | 4 | Implemented |

## 10. Euler -> ForwardEuler migration

输入、内部 ID、日志和测试均使用明确的 `forwardEuler`/`ForwardEuler`。旧 `Euler` 输入 alias 未保留。Forward Euler 更新数学仍调用原有 divergence update，并在相同位置失效 thermodynamic cache、publish 和 validate。

状态：**Implemented**

## 11. RK4 -> ClassicalRK4 migration

输入、内部 ID、日志和测试均使用 `classicalRK4`/`ClassicalRK4`，避免把任意四阶 RK 与经典四级 RK 混为一谈。旧 `RK4` 输入 alias 未保留。

状态：**Implemented**

## 12. stageCount authority

stage count 只存于 resolved `TimeRecipe`。`explicitStageCount()` 已删除。SolvePlanner 用 `recipe.stageCount` 编译 `StageLoop` repetition；`Time::Explicit` 用同一 recipe 校验 stage index 和选择局部 stage 公式。不存在第二份 stage-count switch。

状态：**Implemented**

## 13. SolvePlanner before/after

迁移前，Planner 接收 time string，重新 parse scheme，再通过独立 helper 推导 stage count。迁移后，Planner 接收 `const TimeRecipe&`，只检查 topology 并读取 stage count；它不比较 ForwardEuler、SSPRK3 或 ClassicalRK4，也不解析这些名字。

显式 plan 顺序保持：

```text
flow.step.prepare
flow.dt.compute
flow.step.begin
StageLoop x TimeRecipe.stageCount
    explicit.stage.execute
optional IBM operation
flow.step.commit
time.commit
```

状态：**Implemented**

## 14. Time::Explicit remaining responsibility

`Time::Explicit` 保留 workspace、begin 和单次 `executeStage()`，并在 provider 内依据 `TimeRecipeId` 执行 Forward Euler、SSPRK3 或 Classical RK4 的系数与组合数学。它不拥有全局 stage loop；stage index 由 `PlanExecutor` 的 `StageLoop` context 提供。

状态：**Implemented**

## 15. ddtDispatch audit

全局 caller 审计确认旧 `ddtDispatch`、solver/discretization time bridge 和 methods/numerics/time Euler/RK4 lifecycle 不再有 production 或 meaningful test caller。显式时间 production path 为：

```text
CompiledSolvePlan
  -> StageLoop
  -> explicit.stage.execute
  -> Time::Explicit::executeStage
```

状态：**Implemented**

## 16. Deleted time dispatch files/functions

已删除：

- `src/solver/discretization/time/SF_time.h`
- `src/methods/numerics/time/SF_time.h`
- `src/methods/numerics/time/SF_Euler.h/.cpp`
- `src/methods/numerics/time/SF_RK4.h/.cpp`
- 对应 CMake targets、link entries 和 umbrella includes
- `ddtDispatch()` 与旧 stage-count helper

状态：**Implemented**

## 17. Dead interfaces/helpers deleted

已删除 `TimeScheme`、`parseTimeScheme()`、`explicitStageCount()`、Equation execution metadata、旧 source helper，以及 physics contributor 的 time 参数。未增加 alias、deprecated forwarding wrapper、`ExplicitSolver`、`ImplicitSolver` 或空 NumericalCompiler facade。

状态：**Implemented**

## 18. Numerical invariants

以下内容保持不变：Forward Euler、SSPRK3 与 Classical RK4 系数和 stage time；RHS assembly；WENO/TENO；flux splitting；diffusion；CFL；boundary/halo 顺序；canonical face COPY；GlobalDof SUM；IBM projection/KKT insertion；PISO matrix/correction；MPI state ownership。

本轮只改变 authority、命名、typed resolution、dead dispatch 和 Equation metadata。

状态：**Implemented**

## 19. Regression results

| 验证 | 结果 | Status |
|---|---|---|
| normal build | 全部 targets 编译和链接完成 | Implemented |
| test build | 编译和链接完成 | Implemented |
| architecture checker | allowlisted edges = 10；无新增依赖 | Implemented |
| Physics independent of TimeRecipe | 三种 recipe 的 RawEquationSystem signature 相同 | Implemented |
| explicit stage mathematics | stage count、constant RHS、publish/validate 检查完成 | Implemented |
| Forward Euler WENO7 one-step | `1ab9ccf4cfb3083459de82fa2312e9eaafc0dbe4be896cdf354d9aafd9d55d07` | Implemented |
| Classical RK4 Ghost IBM | `b9d70a4855121102dfb351841332986ae12080862f0facc2aa696d194f329b0e` | Implemented |
| PISO regression | `11cba96d719f1b00ced3a01237ba489c636ffec13de3c0c516144af47ca88040` | Implemented |
| Velocity Forcing IBM | serial one-step 完成 | Implemented |
| fractional DFM IBM | serial one-step 完成 | Implemented |
| surface KKT IBM | 72 iterations，relative residual `8.04775e-10` | Implemented |
| MPI-2 forcing IBM | closure `max|Ju-Us| = 1.9613e-12` | Implemented |
| `numericalFluxContract` | 迁移前已存在的 reconstructed Rusanov constant-stencil mismatch；本轮未改公式、容差或期望 | Legacy |

指定的 `build-phase1a-tests` build tree 当前不存在；最终 CTest 在已重新配置为 `BUILD_TESTS=ON` 的 `build-tests` 上执行。四项中三个本轮相关测试完成，唯一未完成项为上表既有数值通量问题。

## 20. Architecture guards

`tools/check_architecture.py` 现禁止：

- Equation AST 出现 execution/time/Newton/HYPRE authority
- physics composition 接收或解释 time recipe/scheme/string
- Planner parse 时间字符串、比较具体 RK 名称/ID或恢复独立 stage-count helper
- 恢复 `ddtDispatch`、旧 time dispatcher 文件或 global solver-family lifecycle
- `ExplicitSolver`、`ImplicitSolver`、`RK4Solver`、`SDIRKSolver` 类型

状态：**Implemented**

## 21. Updated runtime graph

```text
Input
  ├── Physics -> Equation System
  └── Numerical selection -> resolved TimeRecipe
                  │
ExecutableEquationSystem + TimeRecipe
                  ↓
             SolvePlanner
                  ↓
          CompiledSolvePlan
                  ↓
       ProviderRequirements
                  ↓
          Concrete Providers
                  ↓
             OpRegistry
                  ↓
            PlanExecutor
```

状态：**Implemented**

## 22. Remaining Legacy

- OpenFOAM-like compatibility reader 仍识别其历史配置 key，并在 IO 边界立即解析为 typed recipe；该字符串不进入 Planner 或 runtime dispatch。
- Eulerian pressure coupling 继续使用已迁移的 structured pressure plan 和现有 numerical providers；本轮没有改变其 PISO/PIMPLE 数值算法。
- `numericalFluxContract` 的 constant-stencil mismatch 属于迁移前数值问题。

状态：**Legacy**

## 23. Remaining Unsupported

`backwardEuler`、`SDIRK2`、`ARK2`、runtime-defined tableau、runtime-defined Term/Time recipe 和隐式 stage provider均未注册。选择未知 recipe 会 fail fast；不存在 Forward Euler fallback，也没有空 placeholder implementation。

状态：**Unsupported**

## 24. Global symbol search

生产源码与测试的全局搜索确认以下旧 symbol 无 caller：

```text
TimeScheme
timeIntegrator
EvaluationMode
ExplicitSource
ImplicitSource
implicitSource
explicitStageCount
ddtDispatch
SF_Euler
SF_RK4
ExecutionPolicyKind::ExplicitStages
```

剩余命中只位于 architecture obsolete-symbol guard，或 level-set compatibility deprecation 文本中的历史配置 key，不构成 execution authority。

状态：**Implemented**

## 25. Final authority table

| 问题 | Authority / Answer |
|---|---|
| 谁决定 governing equations？ | Equation System |
| Equation::Term 是否知道 explicit / implicit？ | 否 |
| 谁决定当前使用 Classical RK4？ | Resolved TimeRecipe |
| 谁决定 RK4 有几个 stage？ | TimeRecipe |
| 谁决定 stage coefficient mathematics？ | `Time::Explicit` numerical provider |
| 谁决定 StageLoop 执行几次？ | CompiledSolvePlan，次数来自 TimeRecipe |
| 谁解析用户字符串？ | IO / recipe resolution，只解析一次 |
| SolvePlanner 是否重新 parse `RK4`？ | 否 |
| Physics builder 是否知道 time recipe？ | 否 |
| ddtDispatch 是否仍存在第二套 time authority？ | 否 |
| 用户是否能够自己定义新 time recipe？ | 否 |
| 用户是否单独设置 `implicit=true`？ | 否 |

状态：**Implemented**

## 26. Recommended next task

下一阶段应处理 **Term Recipe Authority**：让 WENO/TENO、flux splitter、central diffusion 和 source evaluation 逐步成为 built-in TermRecipe contract，并由数学 Term 与 recipe 共同生成 numerical operation requirements。继续复用现有 kernel，先完成 term binding 与 capability validation，再考虑真正承担 lowering 的 NumericalCompiler。第一个 implicit TimeRecipe 应在该 authority 完成后实现。

状态：**Interface-only**
