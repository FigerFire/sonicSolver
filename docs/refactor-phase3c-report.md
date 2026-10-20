# Phase 3C 报告——证明 Final Tail Contract 并合并 Density-Based Algorithm

## 1. Final Tail 读写分析

| 操作 | 读取 | 写入 | 需要新鲜 replica | 产生 stale replica | 需要 derived closure |
| --- | --- | --- | --- | --- | --- |
| `flowAlgorithm.correct` | predicted conservative owner state、constraint/system；通过 callback 准备需要的 boundary state | forcing/constraint 的 conservative owner state；可能的 IBM internal state | correction callback 自行准备 | `wroteConservativeState=true` 时 owner state 需要 publication | density correction 不读取 thermo cache |
| `publishCorrectedConservativeState` | `wroteConservativeState` | `WriteOwned("conservative")` | 否 | 否 | 否 |
| homogeneous `commitStep` | model diagnostics | diagnostics | 否 | 否 | 否 |
| legacy multiphase `commitStep` | owned field、auxiliary state | previous level 与 auxiliary state；自身完成 auxiliary `WriteOwned`/`ReadHalo` 和 boundary/coupling refresh | 其内部准备 auxiliary halo | conservative 否 | 不读取 conservative ghost/thermo cache |
| interface `commitStep` | `phi` 与 interface model state | reinitialization / complete-time-step / curvature state；自身准备 `phi` halo | 其内部 apply boundary 并准备 `phi` halo | conservative 否 | 不读取 conservative ghost/thermo cache |
| final `prepareBoundaryState` | corrected conservative owner state、`t+dt`、IBM boundary port | physical conservative boundary/ghost；ghost IBM 时 publication 与 halo refresh | BoundaryPipeline/manual path 自行准备 conservative halo | 否，结束后可读 | 为下一空间读取和 field-consuming output 提供 read-ready state |
| clock commit | 已完成的 model commit 与 final readiness | `bundle.time`、`bundle.step` | 否 | 否 | 否 |
| observer | committed clock；其消息不直接读取 field | observer sink | output 若读取 field，则 final readiness 已完成 | 否 | 需要的 field state 已准备 |

## 2. commitStep Contract

实际 concrete `commitStep` 实现表明，commit 负责 model/history：homogeneous
只报告 diagnostics；legacy multiphase 写 previous time level 后自行同步 auxiliary
field；single/multi interface 都自行闭合 `phi` 并准备其 halo。没有实现读取 final
conservative ghost/boundary value 或 conservative thermodynamic cache。

因此，conservative final boundary closure 不必位于 `commitStep` 前。这个结论来自
`SF_equationCoupling.cpp` 和 `SF_interfaceCoupling.cpp` 的实现，不是按惯例推断。

## 3. Final Boundary / Closure Contract

`DensityBasedRHS::prepareBoundaryState` 负责 existing BoundaryPipeline/manual
physical boundary、conservative halo 和 ghost IBM preparation。它在 target
`t+dt` 后运行，使 timestep 返回时的 conservative read path 已就绪。interface
与 legacy multiphase 的 auxiliary closure 仍由各自 `commitStep` 负责，未移动到
boundary pipeline。

## 4. 选择的 Final Ordering 及理由

统一后的顺序为：

```text
DensityBasedTime complete: predicted Q^(n+1)
↓
flowAlgorithm.correct(target=t+dt)
↓
publish corrected conservative owner state iff wroteConservativeState
↓
correction validation / message
↓
equationSystem.commitStep(all patches)
↓
prepare final conservative boundary / halo / IBM ghost state at t+dt
↓
bundle time/step commit
↓
observer 与 Time::Driver output
```

`commitStep` 位于 final conservative closure 前，因为它不读取 conservative
ghost/thermo，而自身处理所需 auxiliary state。final closure 位于 commit 后，使
observer 及紧随 `advance` 的 output 都看到 read-ready `Q^(n+1)`。clock 最后提交，
避免 tail 抛异常时出现 clock 已前进而 model commit 未完成的状态。

原先 forcing correction 会跳过 single-path final closure；该 special tail 已删除。
forcing 的 final preparation 使用下一 timestep 开始也会使用的 `t+dt` target time，
不改变 owned conservative update、RK algebra、canonical COPY 或 residual SUM。

## 5. Observer Contract

time-step observer 只消费成功提交的 `StateBundle::{time,dt,step}` 并发送结构化
message，不直接读取 field。Time::Driver 在 `advance` 返回后才调用 output callback；
因此 output 也发生在 final boundary/halo readiness 与 clock commit 之后。

## 6. Conservative Correction Path

`FlowAlgorithmResult::wroteConservativeState` 仍是唯一 publication 判据：

```text
wroteConservativeState = false → 不增加 conservative WriteOwned publication
wroteConservativeState = true  → WriteOwned("conservative")
```

algorithm 不检查 Ghost、direct forcing、DLM 或其他 IBM 类型。RK4
velocity-forcing smoke 实际覆盖了 `true` path；ghost IBM smoke 覆盖了无
conservative correction 的 path。

## 7. Lifecycle Before

### CompressibleAlgorithm

```text
single field → density RHS/time → correction → commit →
final boundary（forcing 时跳过）→ clock → observer
```

### MultiPatchAlgorithm

```text
field vector → density RHS/time → correction → final boundary →
commit → clock → observer
```

两条 lifecycle 的 final closure/commit order 不同，且前者存在 forcing special
branch。

## 8. Lifecycle After

```text
StateBundle::patches (one or many)
→ stable CFL reduction / beginStep / optional bound transport port
→ DensityBasedTime
→ flow correction
→ shared finishDensityStep
→ committed observer/output
```

patch count 不参与 density lifecycle 选择。`DensityBasedRHS` 仍负责 spatial stage
orchestration；`DensityBasedTime` 仍负责 explicit tableau；`ExecutionRuntime` 仍
负责 MPI reduction、halo、canonical COPY 和 residual SUM。

## 9. Algorithm Classes：2 → 1

| 指标 | 改造前 | 改造后 |
| --- | ---: | ---: |
| Density-based Algorithm class | 2 | 1 (`CompressibleAlgorithm`) |
| Density public timestep lifecycle | 2 | 1 |
| Density RHS orchestration | 1 | 1 |
| Density explicit time integration | 1 | 1 |
| Final tail implementation | 2 | 1 (`finishDensityStep`) |
| single/multi density lifecycle branch | 1 | 0 |
| `MultiPatchAlgorithm` source references | N | 0 |

`CompressibleAlgorithm` 仍保留 pressure-based single-field route，且该 route 未被
本次 density merge 修改。

## 10. 删除的文件

- `src/solver/algorithm/SF_multiPatch.h`
- `src/solver/algorithm/SF_multiPatch.cpp`

同时从 algorithm CMake target 与 umbrella header 删除它们。应用层的
`run/SF_multiPatch.cpp` 是 mesh/application runner，不是旧 Algorithm，实现仍保留。

## 11. 删除的 Compatibility Path

- 删除 `MultiPatchAlgorithm` class、constructor、`advance`、private helper 与所有
  application caller。
- 未保留 inheritance shell、type alias、forwarding wrapper 或 `advanceMulti` API。
- 删除 single forcing correction 的 final-boundary skip branch。

## 12. Public API 改造前/后

application-level density API 现在只有 `CompressibleAlgorithm::advance`。为保持
Phase 3A 已冻结的空间离散，构造函数接收已有的
`ConvectionThermodynamicContract`：single runner 保留
`EquationSetRusanov`，multi runner 保留 `HighOrderPerfectGas`。这是 immutable
composition policy，不是第二个 state/service/timestep entry，也不按 patch count
选择 lifecycle。

| API | 改造前 | 改造后 |
| --- | --- | --- |
| density algorithm | `CompressibleAlgorithm` + `MultiPatchAlgorithm` | `CompressibleAlgorithm` |
| density `advance` | 两个 class 各一个 | 一个 class 一个 |
| legacy multi compatibility wrapper | 无 | 无 |

## 13. 数值回归

| Case | 结果 |
| --- | --- |
| 4-rank Sod / Euler | 通过；尾部仍为 `time=2.435715e-01`、`dt=5.352234e-04`，与 Phase 3A 基线一致。 |
| single WENO7 SSPRK3 short | 通过；至 `0.002`，stage 1/2/final 全部完成，末步 `dt=1.376924e-04`。 |
| serial ghost IBM RK4 short | 通过；至 `0.001`，四 stage 与 final closure 完成。 |
| serial velocity-forcing BP RK4 short | 通过；至 `0.001`，每步有 flow correction，证明 conservative correction path 可完成 shared tail。 |
| RPI wall boiling smoke | 有界通过；到达 `1e-5`、`2e-5` 的既有 pressure Eulerian summary，未见新失败。 |

未改变 RK coefficient、stage target time、WENO/TENO、Steger–Warming、EOS formula、
viscous/phase-change formula 或 IBM mathematical method。single/multi 既有的
convection thermodynamic contract 也未互相替换，避免把 multi high-order path 改为
Rusanov 或反向改变 single baseline。

## 14. Full Build / CTest

| 检查 | 结果 |
| --- | --- |
| 修改前完整 build | 通过，297/297 target |
| 修改后完整 build | 通过，296/296 target；少一个已删除 algorithm source target |
| 修改后 CTest | 6/6 通过 |

## 15. MPI Ownership 确认

没有修改 canonical face owner、COPY、GlobalDof residual SUM、MPI halo contract 或
rank/partition ownership。最终 algorithm 不检查 MPI/rank；`ExecutionRuntime` 和
`DensityBasedRHS` 继续拥有这些执行语义。4-rank Sod 与 pre-change baseline 一致。

## 16. Architecture Allowlist

`python3 tools/check_architecture.py` 仍报告 10 条 allowlisted dependency edge。
本阶段未新增 reverse dependency，也没有处理 scalar transport、viscous、models 或
HYPRE 的既有 debt。

## 17. levelSetInjectionCase 处置

该 case 是 one-fluid-interface density-based compressible workflow。single/multi
application composition 现在都会为它附着 existing PerfectGas EquationSet，而不
重新初始化 conservative storage、引入 null EquationSet 或 gamma fallback。

运行已越过原先 `StateBundle::equations missing` fail-fast，随后因 case 输入的
phase thermal conductivity 为 0 触发既有 EOS viscous heat-flux validation：
`EquationSet viscous heat flux requires finite positive phase thermal conductivity`。
这不是 lifecycle merge 错误；本阶段未添加物理 fallback。该 case 需要独立修正其
物性配置或在后续 thermodynamics/interface ownership 阶段处理。

## 18. 延后问题

- pressure-based time integration 与 methods/numerics/time ownership；
- scalar transport、multiphase temperature advance、model derivative 与
  `SF_viscous` ownership；
- BoundaryPipeline redesign；
- level-set HJ-WENO 与 thermodynamics/flux generic EOS interface；
- Field responsibility split；
- IBM architecture、distributed monolithic KKT、HYPRE ownership；
- Equation DSL / AssemblyPlan。

Self review：仍有两个 density Algorithm class **否**；仍有两条 density outer
lifecycle **否**；single/multi special lifecycle branch **否**；
`MultiPatchAlgorithm` compatibility wrapper **否**；RK coefficient、stage time、
DensityBasedRHS barrier、EOS capability、canonical COPY、residual SUM、pressure-based
solver **均未改变**；新增 Manager/Adapter/Context **否**。
