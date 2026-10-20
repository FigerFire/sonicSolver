# PISO Pressure-Constraint Vertical Lowering Report

日期：2026-09-17。

## 1. Previous effective authority

迁移前，`RawEquationSystem` 已能描述 pressure constraint，`PressureConstraintTransformer` 和
`CompiledSolvePlan` 也能打印 pressure correction，但真实顺序仍由
`CompressibleAlgorithm::stepPressure()`、`IFlowAlgorithm` 和
`PressureBased::Corrector::correct()` 共同拥有。`CompiledSolvePlan` 最终只选择
`LegacyPressureExecutionAdapter`，没有驱动 predictor、pressure solve 和 correction。

本轮属于“改变怎么解 equation”。物理状态继续由 `StateBundle` 中唯一 patch `Field` 持有，
solver workspace 继续由 algorithm/corrector 持有，MPI backend ownership 未改变。

## 2. Raw equation cleanup

pressure composition 不再创建 raw `pPrime` 或 raw `E_PRESSURE`。Raw system 只保留 production
predictor 已经实际推进的 `E_MASS`、`E_MOMENTUM`、`E_ENERGY`，authoritative pressure `p`，以及
`C_INCOMPRESSIBILITY`。

这里必须准确说明当前数值实现：现有 single-fluid pressure production path 是 PerfectGas
conservative predictor 后接 pressure correction，并不是纯 `U/p` 的 constant-density、isothermal
离散系统。因此本轮没有为了匹配理想化示意图而删除 production 实际推进的 mass/energy，或新写
另一套 PISO 公式。严格的 constant-density `U/p` composition 和 operator 尚未实现，状态列为
**unsupported**。

Validator 新增 fail-fast：raw unknown 中出现 `pPrime`，或 raw equation 中出现任何
`AlgorithmicDerivedEquation`，启动即失败。

## 3. Transformation lowering

`PressureConstraintTransformer` 只匹配：

- `C_INCOMPRESSIBILITY`；
- compatible `E_MOMENTUM`；
- pressure multiplier capability `p`。

它不检查 density/pressure solver family，也不调用 MPI。匹配后记录
`origin=Generated, source=pressureConstraint`，生成 pressure-correction workspace、derived
equations、correction operators 和 compiled bindings。

## 4. Generated algorithmic equations

Transformer 生成：

- `pPrime`：`SpecializedExecutor:pressureCorrection`，不初始化、不 restart、不 output，也不要求
  permanent runtime storage；
- `E_MOMENTUM_PREDICTOR`；
- `E_PRESSURE`；
- `OP_PRESSURE_UPDATE`；
- `OP_VELOCITY_CORRECTION`；
- `OP_FLUX_CORRECTION`。

这些对象只存在于 `ExecutableEquationSystem`。`sonicSolver explain` 会同时显示 raw、generated
provenance 和 executable system。

## 5. CompiledEquation bindings

新增的类型化 binding 如下：

| Equation | Operator | Storage/read-write contract | Matrix | RHS |
|---|---|---|---:|---:|
| `E_MOMENTUM_PREDICTOR` | `ExplicitMomentumPredictor` | conservative `Q[0:5]`、momentum residual slice `[1:3]`、boundary freshness、ReadHalo/WriteOwned | No | Yes |
| `E_PRESSURE` | `PressureCorrection` | conservative `Q`、transient `pPrime`、`pressureMatrix`、`pressureRhs` | Yes | Yes |

Validator 检查 compiled equation 必须对应 executable descriptor、必须具有 generated provenance、
必须有合法 storage/component/access binding。Executor 不根据 equation ID 重新解释数学意义。

## 6. Structured PISO plan

满足当前最小 capability 时，SolvePlanner 生成并冻结：

```text
Sequence
  PreparePressureStep
  AssembleMomentumPredictor
  SolveMomentumPredictor
  Loop PressureCorrectors(1)
    PreparePressureCorrection
    AssemblePressureCorrection
    SolvePressureCorrection
    CorrectVelocity
    PreparePressureUpdate
    CorrectFlux
    CommitPressureCorrection
  CommitPressureStep
```

`CorrectVelocity -> PreparePressureUpdate` 保持现有 production 的真实数学顺序：旧实现先用
`pPrime` 修正 momentum，再从修正后的 momentum 与旧 energy 恢复当前 pressure，最后形成目标
pressure 并重建 total energy。交换两步会改变 kinetic-energy contribution，因而本轮没有按概念
示意图重排。

## 7. Generic plan executor

`GenericPlanExecutor` 仅实现本 vertical slice 需要的 `Sequence`、`Loop`、`Assemble`、`Solve`、
`Correct`、`Update` 和 `Commit`。其他 node kind 明确 fail-fast。

Executor 递归遍历 `SolvePlanNode::children` 并调用 `PlanOperationRegistry` 的类型化 callback；源文件
不包含 PISO operation 顺序、不检查 equation ID，也不查询 turbulence、phase 或 IBM flag。

## 8. Existing operator reuse

没有重写 pressure mathematics。现有实现被拆成同一 `PressureBased::Corrector` 的阶段接口：

- `assemble()`：原 row map、fixed-pressure boundary、compressibility、Laplacian coefficient、jump、
  matrix 与 RHS；
- `solve()`：原 HYPRE solve 和 transient `pPrime` synchronization；
- `correctVelocity()`：原 pressure-gradient momentum correction；
- `preparePressureUpdate()`：原 pressure relaxation；
- `correctFlux()`：声明当前 FDM path 没有 persistent pressure face-flux storage；下一次 spatial
  assembly 从 corrected cell state 重新构造 face flux；
- `commitPressureUpdate()`：原 EOS/PerfectGas energy reconstruction 与 thermo-cache invalidation。

旧 `correct()` 仍用于 legacy capabilities，但只按同一阶段 API 调用，因此不再拥有另一份公式。

## 9. State/workspace ownership

| Object | Before | After |
|---|---|---|
| `U/rho/rhoE` physical state | `Field` | unchanged |
| authoritative clock | `StateBundle` | unchanged |
| `pPrime` | local temporary inside monolithic corrector | `Corrector::Workspace`, generated compiled binding |
| pressure matrix/RHS | local temporary inside monolithic corrector | `Corrector::Workspace` |
| residual/face-flux workspace | `PatchWorkspace` | unchanged |

没有新增 physical-state copy、第二份 `Field` authority 或 persistent face-flux copy。

## 10. Boundary behavior

执行顺序保持：step-start boundary/closure freshness，predictor assembly，predictor publish，target-time
boundary refresh，pressure assembly/solve/correction，thermodynamic state commit，step-end boundary
refresh，clock commit。pressure fixed-value boundary row handling和 interface jump relation未改变。

当前 migrated capability 是 serial/single-patch。Plan/compiled binding 只声明 ReadHalo/WriteOwned；
MPI implementation 仍只位于 `ExecutionRuntime`/backend。Canonical COPY 与 residual SUM 语义没有修改。

## 11. Numerical A/B regression

在修改 staged corrector 与 runtime routing 前，先用同一 case、mesh、initial state、CFL、Euler、
WENO5/Steger-Warming 和 HYPRE 配置冻结 legacy 输出；修改后运行 generic plan。

| Quantity | Legacy | Generic plan |
|---|---:|---:|
| `dt` / final time | `1.659315e-06` | `1.659315e-06` |
| pressure iterations | 2 | 2 |
| pressure residual | `4.38769e-11` | `4.38769e-11` |
| max divergence before | `2.01179` | `2.01179` |
| max divergence after | `2.01052` | `2.01052` |
| HYPRE rebuilds/solves | 1 / 1 | 1 / 1 |
| final VTS SHA-256 | `11cba96...a88040` | `11cba96...a88040` |

最终 VTS 字节完全一致。其 `U`、`p`、`rho`、`TotalEnergyDensity` 等全部输出数组因此逐值一致；
baseline 文件还保存这些数组的 min/max/L2。matrix/RHS 公式没有替换，而是从原函数原样移动到
`assemble()`；pressure residual、correction 后 continuity 和最终 state 提供了运行期闭环。

`tools/check_piso_regression.py` 可重复运行该比较。它同时验证 backend 为
`GenericPlanExecutor`、adapter section 不含 legacy pressure adapter、step diagnostics 和最终数组。

## 12. Architecture tests

`test_pisoArchitecture` 直接检查对象结构，不比较 explain 字符串：

- raw 中没有 `pPrime`/`E_PRESSURE`；
- transformer 生成 predictor、pressure equation 和 operators；
- matrix/RHS/storage/component bindings 完整；
- plan root/loop/operation 顺序正确；
- minimal capability backend 是 `GenericPisoPlan`，legacy adapter 列表为空；
- 两次 pressure corrector 明确路由到 legacy backend。

验证结果：

| Check | Result |
|---|---|
| full application build | Passed，305 targets；最终 `sonicSolver` rebuild 262 targets |
| `python3 tools/check_architecture.py --quiet` | Passed；633 source files，allowlist=10 |
| `test_pisoArchitecture` | Passed |
| `check_piso_regression.py` | Passed；legacy/generic final VTS byte-identical |
| existing `numericalFluxContract` | Failed：constant-stencil reconstructed Rusanov mismatch；与本 PISO lowering 无关，未改 tolerance/WENO |

Build tree 仍偶发报告 `ninja: premature end of file; recovering`；Ninja 自动恢复 metadata 后完成
全部编译和链接，没有 source compile/link failure。

## 13. Capability routing

Generic path 当前要求 executable system/policy 同时满足：Euler、PISO、一个 pressure corrector、零
non-orthogonal corrector、单 incompressibility constraint，以及精确的 production equation set。
Runtime 再验证 exactly one patch，且无 transported registry、transport model、auxiliary equation
system 或 IBM service。条件不满足时不会改变数学算法：plan 会显式标记对应 legacy backend，或
generic runtime validation fail-fast。

`explain test/pressureConstraintPiso` 显示 `GenericPlanExecutor` 和 `(none)` legacy adapters；
`levelSetLaplaceCase` 显示 `LegacyPressureExecutionAdapter`，证明 routing 可见。

## 14. Legacy adapter removal for migrated capability

对上述最小 production capability，`CompressibleAlgorithm` 不再构造或调用 pressure
`IFlowAlgorithm`/`LegacyPressureExecutionAdapter`。它创建 staged `Corrector` operator provider，
并把 callback 注册到 plan executor。Plan 是唯一顺序 authority。

Public API before/after：`Corrector::correct()` 保留给未迁移 capability；新增 staged operation API
供 typed plan binding。未新增 solver family、ServiceLocator 或第二套 timestep lifecycle。

## 15. Unsupported pressure capabilities

状态分类：

| Capability | Status |
|---|---|
| 当前 single-patch single-fluid Euler/PISO(1/0) production pressure path | **implemented** |
| Compiled storage/operator/sync descriptors | **implemented** for this slice |
| `Sequence/Loop/Assemble/Solve/Correct/Update/Commit` executor | **implemented** |
| other structured node kinds | **interface-only**, executor fail-fast |
| multiple pressure correctors、non-orthogonal、SIMPLE/PIMPLE、level set、turbulence、IBM、multiphase | **legacy** |
| exact constant-density/isothermal `U/p` raw system and discretization | **unsupported**；当前 repo 没有可提取的 verified implementation |
| distributed generic PISO | **unsupported** |

## 16. Remaining work

下一步应独立建立并验证真正 constant-density/isothermal Momentum + `div(U)=0` production
operator，再让同一 transformer/compiled-plan infrastructure消费它。之后可逐项迁移 non-orthogonal
loop、multiple correctors、SIMPLE/PIMPLE 与 distributed pressure solve。它们都不得通过 hidden
fallback、solver-family判断或修改本轮冻结的 density high-order/IBM/KKT 数值路径来完成。

