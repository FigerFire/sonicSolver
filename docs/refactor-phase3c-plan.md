# Phase 3C 计划——证明 Final Tail Contract 并合并 Density-Based Algorithm

## 范围与预检

本阶段只合并 density-based 的 outer timestep lifecycle。保留
`CompressibleAlgorithm` 作为最终类名，删除 `MultiPatchAlgorithm`；不触及
pressure-based timestep、RK algebra、DensityBasedRHS 的空间步骤、WENO/TENO、
Steger–Warming、EOS、IBM 数学、Field storage 或 BoundaryPipeline API。

源代码修改前的验证已完成：

| 检查 | 结果 |
| --- | --- |
| `cmake --build build --parallel 4` | 成功，297/297 target 完成 |
| `ctest --test-dir build-phase1a-tests --output-on-failure --parallel 1` | 6/6 通过 |
| `python3 tools/check_architecture.py` | 10 条 allowlisted dependency edge |
| `mpirun -np 4 ./build/sonicSolver run test/sodCase` | 完成 500 步，尾部为 `time=2.435715e-01`、`dt=5.352234e-04` |

## 现有生命周期

### CompressibleAlgorithm（single density）

```text
prepareBoundaryState(t)
→ CFL / beginStep / transport correction
→ DensityBasedTime::advance({field})
→ flowAlgorithm.correct(target=t+dt)
→ publish corrected conservative owner state (wroteConservativeState)
→ correction validation / observer message
→ equationSystem.commitStep(field)
→ final prepareBoundaryState(t+dt)，但 forcing correction 时跳过
→ bundle time/step commit
→ observer time-step message
```

### MultiPatchAlgorithm（multi density）

```text
prepareBoundaryState(t)
→ all-patch CFL / beginStep
→ DensityBasedTime::advance(fields)
→ flowAlgorithm.correct(target=t+dt)
→ publish corrected conservative owner state (wroteConservativeState)
→ correction validation / observer message
→ final prepareBoundaryState(t+dt)
→ equationSystem.commitStep(fields)
→ bundle time/step commit
→ observer time-step message
```

两条路径在 `DensityBasedTime` 之后的差异只有 final conservative closure 与
`commitStep` 的相对顺序，以及 single forcing correction 的 final-closure skip。

## Final Tail Read/Write 表

| Operation | Reads | Writes | 需要 fresh replica | 产生 stale replica | 需要 derived closure |
| --- | --- | --- | --- | --- | --- |
| `flowAlgorithm.correct` | predicted conservative owner state；forcing 时 constraint/system；通过 callback 准备其所需 boundary state | forcing/constraint 的 conservative owner state；可能 IBM internal state | correction 自己通过 `prepareBoundaryState` callback 取得 | `wroteConservativeState=true` 时 conservative replica 失效 | pressure correction 可需要；density correction 本身不读取 thermo cache |
| `publishCorrectedConservativeState` | `wroteConservativeState` | Runtime `WriteOwned("conservative")` publication | 否 | 否；仅声明 owner state 已更新 | 否 |
| `equationSystem.commitStep`（homogeneous） | model diagnostic state | diagnostic output | 否 | 否 | 否 |
| `equationSystem.commitStep`（legacy multiphase） | owned field、registered auxiliary state | previous time level、auxiliary state；随后自身 `WriteOwned`/`ReadHalo` auxiliary field，并重建 auxiliary boundary/coupling | auxiliary halo 由 commit 内部准备 | 否 | 不读取 conservative ghost/thermo cache |
| `equationSystem.commitStep`（single interface） | `phi`、interface model state | reinitialization/complete-time-step state、curvature；内部 prepare `phi` halo 与 interface state | 其 callback 自行 apply boundary/prepare phi halo | conservative state 否；auxiliary writes 由自身发布 | 不读取 conservative ghost/thermo cache |
| `equationSystem.commitStep`（multi interface） | `phi`、interface model state | complete-time-step state | 自行 apply `phi` boundary/prepare phi halo | conservative state 否 | 不读取 conservative ghost/thermo cache |
| final `prepareBoundaryState` | corrected conservative owner state、target `t+dt`、IBM boundary port | physical conservative boundary/ghost values；IBM ghost 时 `WriteOwned("conservative")` 后 `ReadHalo` | pipeline/manual implementation owns required conservative halo | 否，完成后 replica 可读 | 使下一 space read 与 output-side derived reads ready |
| clock commit | committed state、`dt` | `bundle.time += dt`、`bundle.step++` | 否 | 否 | 否 |
| observer | committed bundle clock、final read-ready field state | observer sink only | 对 field-consuming observer/output：是 | 否 | 是，来自 final boundary closure |

### commitStep 的结论

所有当前 concrete `commitStep` 实现都不读取 final conservative ghost/boundary
value 或 conservative thermodynamic cache。interface 与 legacy multiphase commit
需要的 auxiliary boundary/halo 由它们自身完成。因此 conservative final closure
不必位于 `commitStep` 之前。

## 选择的 Final Tail 顺序

```text
DensityBasedTime complete: Q^(n+1) predicted
↓
flowAlgorithm.correct(target=t+dt)
↓
publish corrected conservative owner state iff wroteConservativeState
↓
correction validation / diagnostic message
↓
equationSystem.commitStep(all patches)
↓
prepare final conservative boundary/halo/IBM ghost state at t+dt
↓
bundle time/step commit
↓
observer / Time::Driver output
```

理由：`commitStep` 只提交 model/history 或自行准备其 auxiliary data，因而不需要
final conservative closure。将 final closure 放在 commit 后使 observer 和紧随
`advance` 的 output 都看到 read-ready `Q^(n+1)`、边界/halo 与 derived state。该顺序
与现有 non-forcing single 路径一致，消除 multi 的 pre-commit closure 与 forcing
special skip。forcing 路径的额外 final preparation 使用和下一 timestep 开始相同的
`t+dt` target time；它不改变 owned conservative update、RK、canonical COPY 或
residual SUM。

`DensityBasedRHS::publishCorrectedConservativeState` 继续只由
`FlowAlgorithmResult::wroteConservativeState` 驱动；不以 IBM method 分类。

## 保持既有 convection thermodynamic contract

审查发现现有 single composition 使用 `EquationSetRusanov`，multi composition
使用 `HighOrderPerfectGas`。把 multi 改为 Rusanov 明确违反 Phase 3A 决策；把 single
改成 high-order 会改变其既有数值基线。因此 final `CompressibleAlgorithm` 将在构造时
接收现有 enum `ConvectionThermodynamicContract`，由 application composition 显式
选择，且不依据 patch count 分支：

| 既有 composition | 保留 contract |
| --- | --- |
| single runner | `EquationSetRusanov` |
| multi runner | `HighOrderPerfectGas` |

这不是新的 service、state 或 lifecycle；它保留现有空间离散能力。density mainline
对 one patch 和 N patches 均只执行同一 lifecycle。generic EOS/high-order flux
统一继续 deferred。

## levelSetInjectionCase 处置

`OneFluidInterface` 是 density-based compressible workflow：它注册 `phi` 作为
transported state，同时仍求解 conservative density system。single/multi runner 目前
只在 `SingleFluid` workflow 绑定 PerfectGas EquationSet，导致
`StateBundle::equations` 为空。若 field 已由既有 initialization 完成，本阶段只在
application composition 中创建并附着同一 existing PerfectGas binding 到所有 patch，
不重新初始化 conservative storage、不加入 gamma fallback，也不改变 level-set
equation ownership。

## 计划修改

1. 把 `CompressibleAlgorithm` 的 density mainline 改为接收
   `StateBundle::patches`，保留 pressure route 的 exactly-one-patch guard 与既有
   `ddtDispatch` 路径。
2. 在普通内部 density tail helper 中实现上文证明的 correction、commit、final
   readiness、clock、observer 顺序；不放入 CFL/RHS/RK。
3. 把 multi application caller 迁移到 `CompressibleAlgorithm`，在构造时传入其既有
   high-order contract；single caller 继续使用其既有 Rusanov contract。
4. 为 OneFluidInterface 的 single/multi application composition 附着已有
   PerfectGas EquationSet。
5. 删除 `SF_multiPatch.{h,cpp}`，移除其 CMake 与 umbrella-header references；不留
   alias、wrapper 或继承壳。
6. 用现有 test facilities 增加最小 final-tail regression，至少验证
   `wroteConservativeState=false/true` 时 publication、commit、final readiness、
   clock 与 observer 的顺序。

## 不变量

- Euler / SSPRK3 / RK4 coefficient 与 stage target time 不变。
- DensityBasedRHS 的 boundary/halo、canonical face COPY、GlobalDof residual SUM
  及 IBM stage order 不变。
- `StateBundle::{time,dt,step}` 只在 successful tail 完成后更新一次。
- `ExecutionRuntime` 继续拥有 reduction、halo、canonical COPY 与 SUM；algorithm
  不读取 MPI rank 或 patch-count mode。
- pressure solver 及 `methods/numerics/time` pressure route 不变。
