# Explicit Time Plan Authority 迁移报告

日期：2026-09-19

## 1. Scope

本轮只迁移 density formulation 显式时间推进的 stage control-flow authority。数值公式仍由 `Time::Explicit` 与 `DensityBasedRHS` 实现；方程、空间离散、IBM 数值实现、MPI backend、application/environment composition 均未重构。

本轮属于“改变怎么解 equation”，不增加物理 term、equation 或 constraint。

## 2. Before lifecycle

迁移前的 production 链为：

```text
CompressibleAlgorithm::stepDensity
  -> Time::Explicit::advance
       -> 内部选择 Euler / SSPRK3 / RK4
       -> 内部循环全部 stages
       -> RHS / publish / validate
  -> PlanExecutor::executeOperation(IBM)
  -> finishDensityStep
  -> clock commit
```

`Time::Explicit::advance` 同时拥有 stage 数学和全局 stage 顺序，IBM 通过主 Plan 之外的单操作入口插入。

## 3. After lifecycle

production 链现为：

```text
CompiledSolvePlan
  -> flow.step.prepare
  -> flow.dt.compute
  -> flow.step.begin
  -> StageLoop x N
       -> explicit.stage.execute
  -> optional ibm.constraint.project / ibm.kkt.solve
  -> flow.step.commit
  -> time.commit
```

`CompressibleAlgorithm::stepDensity()` 只注册现有 numerical callbacks 并执行 `resolved_.solvePlan`。它不再包含 scheme switch 或 stage loop。

## 4. StageLoop semantics

`Run::PlanExecutor` 现在解释 `PlanNodeKind::StageLoop`：按 `repetitions` 顺序执行 child，并在进入前保存、结束后恢复外层 `ExecutionContext`。空 StageLoop 和非正 repetition 均 fail fast。`Subcycle` 仍显式 unsupported。

`SF_PLAN_TRACE=1` 会输出 `stage=1/N` 至 `stage=N/N`，默认仍关闭。

## 5. ExecutionContext

新增的最小 context 只有：

```text
stageIndex   0-based 当前 stage
stageCount   本 StageLoop 的总 stage 数
```

`OpRegistry` 同时接受无参数 callback 和 context-aware callback；旧 pressure/IBM callback 自动包装，只有 `explicit.stage.execute` 读取 stage context。Executor 不检查 Euler、SSPRK3、RK4、density 或 IBM identity。

## 6. Explicit Workspace

`Time::Explicit::Workspace` 持有单个 physical step 跨 stage 的 `Q_n`、`k1..k4` 和 registered-variable snapshots。`begin()` 只初始化 workspace；不调用 RHS、publish、validate、IBM 或 commit。Euler 不创建不需要的 RK snapshot。

stage count 的唯一映射位于公共 `FDM::explicitStageCount(TimeScheme)`：Euler=1、SSPRK3=3、RK4=4，未知枚举 fail fast。Planner 和 `Time::Explicit` 共同读取该映射，没有两份可能漂移的 switch。

## 7. Euler lowering

Plan 生成一个 stage。stage 0 的顺序保持为 RHS(`t_n`) → conservative/registered-variable update → cache invalidation → publish → validate。生产 WENO7 Sod 回归的最终 VTS 与冻结 baseline 字节一致。

## 8. SSPRK3 lowering

三次 `executeStage()` 保持既有 Shu-Osher 系数与 stage time：

| Stage | time fraction | base weight | Euler weight |
| --- | ---: | ---: | ---: |
| 0 | 0 | 0 | 1 |
| 1 | 1 | 3/4 | 1/4 |
| 2 | 1/2 | 1/3 | 2/3 |

stage-level test 从 `time=2.0, dt=0.1` 观察到 `2.0, 2.1, 2.05`，常量 RHS 的最终值为 `9.8`，publish/validate 各三次。

## 9. RK4 lowering

四次 Plan stage 保留旧实现的精确位置：k1 后 half update、k2 后 half update、k3 后 full predictor、k4 后 final `1/6(k1+2k2+2k3+k4)` combination。stage time 为 `t, t+dt/2, t+dt/2, t+dt`。旧实现中 stage-1 RHS 后的 validation 以及每次 publish/validation 的位置均保留，没有新增 `explicit.finalize` 假操作。

## 10. Fused RHS decision

`DensityBasedRHS::assembleAllPatches` 继续一次完成 validated boundary → halo → fused spatial RHS → canonical face handling → source/reduction。Planner 将所有兼容的 `ExplicitStages` policies 聚合成一个 StageLoop。`BoundaryClosure` policy 被该 fused stage 消费，不生成第二次 boundary operation。

## 11. IBM insertion

Ghost/ILW 仍在 fused RHS 的 boundary closure 位置。forcing/constraint IBM leaf 直接属于主 `CompiledSolvePlan`，位于 StageLoop 之后、state commit 之前。production callback 仍调用现有 `flowAlgorithm_->correct()` numerical adapter，未改变 IBM 数学。

实际 Velocity Forcing FTS 与 fractional DFM trace 均为：四次 `explicit.stage.execute` → `ibm.constraint.project` → `flow.step.commit` → `time.commit`。

## 12. Removed legacy functions

已删除：

- `Time::Explicit::advance`
- `PlanExecutor::executeOperation`
- `CompressibleAlgorithm::executeFlowCorrection`
- `CompressibleAlgorithm::finishDensityStep`
- 依赖单操作执行入口的测试

全局搜索没有剩余 production caller 或 forwarding wrapper。

## 13. Numerical invariants

以下语义未改变：Euler/SSPRK3/RK4 系数、stage time、stage state combination、registered scalar update、RHS 公式、boundary/halo 顺序、canonical face COPY、residual/GlobalDof SUM、IBM projection/DFM/KKT、StateBundle ownership、MPI ownership。

ownership 变化仅为 control flow：全局 stage order 从 `Time::Explicit::advance` 移到 `CompiledSolvePlan::StageLoop`；workspace 仍由 solver execution lifetime 持有。

## 14. Architecture tests

- StageLoop test 收集 `0/4, 1/4, 2/4, 3/4` context。
- Planner test 验证 Euler/SSPRK3/RK4 repetitions 为 `1/3/4` 且只生成一个 explicit StageLoop。
- numerical stage test 验证三种 scheme 的 stage time、常量 RHS 最终状态及 publish/validate 次数。
- `tools/check_architecture.py` 增加 guard：StageLoop 必须实现；Executor 禁止 scheme/solver identity；`Explicit::advance`、`executeOperation`、显式 helper 内全局 stage loop、`stepDensity` scheme switch/stage loop 均禁止重新出现。
- architecture checker 通过，allowlisted dependency edges 仍为 10。

## 15. Euler regression

single Sod WENO7，1 step：

```text
dt = 6.6815310478106089e-04
final VTS SHA-256 = 1ab9ccf4cfb3083459de82fa2312e9eaafc0dbe4be896cdf354d9aafd9d55d07
```

与迁移前冻结值完全一致。

## 16. SSPRK3 stage regression

`explicitStageMathematics` 通过。SSPRK3 的 RHS 次数、stage time、Shu-Osher combination 和 publish/validate 位置均由 unit test 冻结。

## 17. RK4 regression

single Ghost IBM，1 step：

```text
dt = 4.3926797327883271e-04
final VTS SHA-256 = b9d70a4855121102dfb351841332986ae12080862f0facc2aa696d194f329b0e
```

与迁移前冻结值完全一致。Pressure PISO control case 的最终 VTS SHA-256 仍为 `11cba96d719f1b00ced3a01237ba489c636ffec13de3c0c516144af47ca88040`。

## 18. IBM regression

- Ghost IBM：serial RK4 通过，冻结 hash 一致。
- Velocity Forcing FTS：serial 1-step 通过，IBM leaf 位于四个 stages 之后。
- fractional DFM self-propelled：serial 1-step 通过，state-closure diagnostics 正常，IBM leaf 顺序正确。
- Velocity Forcing MPI-2：host 环境 2-rank 1-step 通过；distributed mesh/halo/canonical face 初始化和 constraint projection 完成，`max|Ju-Us|=1.9613e-12`。

## 19. Remaining Legacy

`IFlowAlgorithm`、`makeFlowAlgorithm` 和 `flowAlgorithm_` 仍作为现有 IBM numerical adapter。它们不再拥有 stage/commit lifecycle。`Subcycle` 仍为明确 Unsupported。

CTest 共四项：`pisoArchitecture`、`explicitStageMathematics`、`pisoNumericalRegression` 通过；`numericalFluxContract` 仍因迁移前已知的 reconstructed Rusanov constant-stencil mismatch 失败，本轮未修改 tolerance 或 flux 数学。

## 20. Next task

下一任务为 **Flow Algorithm Adapter Removal**：单独审计并决定是否删除 `IFlowAlgorithm`、`makeFlowAlgorithm`、`flowAlgorithm_` 以及 `CompressibleAlgorithm` 中剩余 formulation/lifecycle branch。本轮未实施。

## Final authority table

| Concern | Before | After | Status |
| --- | --- | --- | --- |
| explicit global stage order | `Time::Explicit::advance` | `CompiledSolvePlan::StageLoop` | Implemented |
| stage index | internal loop | `ExecutionContext` | Implemented |
| Euler math | `Time::Explicit` | `Time::Explicit` | Implemented |
| SSPRK3 math | `Time::Explicit` | `Time::Explicit` | Implemented |
| RK4 math | `Time::Explicit` | `Time::Explicit` | Implemented |
| RHS | `DensityBasedRHS` | `DensityBasedRHS` | Implemented |
| IBM position | legacy tail + `executeOperation` | main `CompiledSolvePlan` | Implemented |
| step commit | `CompressibleAlgorithm` hidden lifecycle | Plan Op `flow.step.commit` | Implemented |
| clock commit | `CompressibleAlgorithm` hidden lifecycle | Plan Op `time.commit` | Implemented |
| `IFlowAlgorithm` | legacy adapter | legacy adapter | Legacy |
