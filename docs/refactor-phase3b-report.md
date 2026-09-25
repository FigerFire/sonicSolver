# Phase 3B 报告——统一 Density-Based 显式时间积分

## 1. 时间积分改造前

改造前，density-based timestep 有两套独立维护的实现：
`CompressibleAlgorithm` 通过 `ddtDispatch` 调用 `methods/numerics/time`，而
`MultiPatchAlgorithm` 调用 `SF_multiPatchTime.cpp` 中的显式 Euler / SSP-RK3 /
RK4。两条路径的数学公式相同，但 stage storage、更新循环和 callback 调用位置
相互独立。

改造后，两个 algorithm 都调用内部实现
`solver/algorithm/SF_densityBasedTime.cpp`。二者唯一的输入差别是参与计算的
field vector 长度。

| 指标 | 改造前 | 改造后 |
| --- | ---: | ---: |
| Density Euler 实现数 | 2 | 1 |
| Density SSP-RK3 实现数 | 2 | 1 |
| Density RK4 实现数 | 2 | 1 |
| single/multi density stage orchestration | 2 | 1 |
| 应用层公开时间 API | 0 | 0 |
| 生效中的 density 时间积分源码 | 约 927 行 | 461 行 |

历史 `methods/numerics/time` 源码仍会为冻结不动的 pressure-based 路径编译；
它已不再被 density-based timestep 选用。

## 2. Euler 公式改造前/后

旧路径与统一路径均保留已有 residual 符号约定：

```text
Q1 = Q0 - dt R(Q0, t)
q1 = q0 + dt rhs(q0, t)
```

`DensityBasedTime::stepEuler` 在 committed `StateBundle::time` 对所有 patch
计算一次 RHS，更新每个参与的 field 及其已有的 registered variable，之后调用
既有的 publication 和 closure check。系数与正负号均未改变。

## 3. SSPRK3 公式改造前/后

统一实现保留已有 Shu–Osher 形式：

```text
Q1 = Q0 - dt R(Q0)
Q2 = 3/4 Q0 + 1/4 (Q1 - dt R(Q1))
Q3 = 1/3 Q0 + 2/3 (Q2 - dt R(Q2))
```

registered transported variable 使用对应的 `q + dt rhs` 形式，并保持相同的
`3/4, 1/4` 与 `1/3, 2/3` 权重。既有的非自治 RHS time node 原样保留；本阶段没有
把它改写成另一种 SSPRK 表述。

## 4. RK4 公式改造前/后

统一代码保留 classical four-stage state construction：

```text
k1 = R(Q0, t)
Q2 = Q0 - dt/2 k1
k2 = R(Q2, t + dt/2)
Q3 = Q0 - dt/2 k2
k3 = R(Q3, t + dt/2)
Q4 = Q0 - dt k3
k4 = R(Q4, t + dt)
Qn+1 = Q0 - dt/6 (k1 + 2k2 + 2k3 + k4)
```

transported variable 保留既有的正 RHS 形式，仍使用四份 RHS snapshot 与最终的
`1, 2, 2, 1` 组合。storage layout、residual layout 和 update policy 均未改变。

## 5. Stage Time 改造前/后

`StateBundle::time` 仍是 committed physical clock。integrator 只构造临时 RHS
target time，不会在 stage 中写入 physical clock。

| Scheme | 改造前 RHS target-time sequence | 改造后 RHS target-time sequence |
| --- | --- | --- |
| Euler | `t` | `t` |
| SSPRK3 | `t`, `t + dt`, `t + dt/2` | `t`, `t + dt`, `t + dt/2` |
| RK4 | `t`, `t + dt/2`, `t + dt/2`, `t + dt` | `t`, `t + dt/2`, `t + dt/2`, `t + dt` |

因此 time-dependent boundary condition 和 prescribed IBM motion 仍会看到与
改造前完全相同的 target time。

## 6. Conservative State 的 Stage 生命周期

`StateBundle::patches` 持有 committed `Qn`。integrator 创建临时 `q0` 和 residual
snapshot，将 temporary `Qstage` 写入已有 field storage，通过 Phase 3A contract
发布 stage，最后将同一 field storage 写为 `Qn+1`。RK snapshot 是 workspace，
不是第二份 authoritative physical state；physical clock 不属于该 workspace。

## 7. Transported Variable 的 Stage 生命周期

single patch 使用 `StateBundle::transported`。multi-patch equation-system
coupling 则查询该 patch 既有的 registry。实现不会复制或合并 registry。update
policy 为 Explicit、BoundedExplicit 或 HamiltonJacobi 的变量与 conservative state
使用相同 Euler / SSPRK3 / RK4 stage 和 `dt`；其他 policy 保持既有行为。

这保持当前 `phaseMass` 及 custom scalar 的 policy 语义。

## 8. Shared Time Integrator 设计

新增 [SF_densityBasedTime.cpp](/Volumes/SSD_LPF/sonicSolver/sonicSolver/src/solver/algorithm/SF_densityBasedTime.cpp)
及其小型内部头文件。它接收 participating field vector、既有 bundle、选择的
scheme、可选的 equation-system registry provider，以及三个显式 callback：
all-patch RHS、publication 和 closure validation。

stage 结构仍然是 stage outermost：每个 RHS callback 均在同一 target time 下先
看到全部 patches，之后才发布任何 stage update。integrator 不知道 convection、
WENO/TENO、EOS、boundary、halo、IBM 或 MPI 细节。没有新增 Manager、Adapter、
Context、Strategy 或应用层 public API。

## 9. DensityBasedRHS 集成

`CompressibleAlgorithm` 向 shared integrator 提供
`DensityBasedRHS::assembleAllPatches`、`publishIntegratedState` 和
`validateStateClosure`。`MultiPatchAlgorithm` 提供其已有的 all-patch RHS、
publication 和 validation operation。因此 `DensityBasedRHS` 仍是 density spatial
RHS order、canonical-face COPY、GlobalDof residual SUM、boundary/halo preparation
及 IBM stage behavior 的唯一 owner。

## 10. 已删除的时间实现

`src/solver/algorithm/SF_multiPatchTime.cpp` 中重复的 density orchestration 已在
迁移后删除，没有 forwarding wrapper。`SF_multiPatch.h` 也不再暴露三个 private
multi-patch Euler / SSPRK3 / RK4 helper。

## 11. SF_multiPatchTime.cpp 状态

已删除。其 all-patch storage 和 tableau 行为迁入中立的
`SF_densityBasedTime.cpp`；文件名不再暗示该实现只能由 multi-patch caller 使用。

## 12. methods/numerics/time 状态

这些 methods 被刻意保留，因为 pressure-based branch 仍调用 `ddtDispatch`。
density-based single patch 已不再调用它们。移动或重构 pressure-based time path
超出 Phase 3B 范围，因此没有引入 compatibility forwarding layer，也没有进行
大规模目录迁移。

## 13. Public API 改造前/后

应用仍配置既有 time scheme 并调用 `advance`，没有改变 algorithm public API。
`DensityBasedTime::advance` 仅为两个 algorithm implementation 使用的 solver
内部组合。

| API 类别 | 改造前 | 改造后 |
| --- | ---: | ---: |
| Public density timestep entry | 2 个既有 algorithm entry | 2 个，未改变 |
| Public density time-integrator entry | 0 | 0 |
| Private multi-patch scheme helper | 3 | 0 |

## 14. 数值回归

改造前 architecture checker 报告所需的 10 条 allowlisted edge。改造后仍报告
`Known dependency debt: 10 allowlisted entries`。

| 验证项 | 结果 | 证据 / 限制 |
| --- | --- | --- |
| 修改的 algorithm translation unit | 通过 | `SF_densityBasedTime.cpp`、`SF_compressible.cpp`、`SF_multiPatch.cpp` 均通过 Ninja 编译。 |
| CMake link command | 通过 | 更新 static library 并执行 CMake 生成的 `sonicSolver` link command 成功。 |
| CTest | 通过 | 6/6，零失败。 |
| 4-rank Sod | 通过 | 结束于 `time=2.435715e-01`、`dt=5.352234e-04`，与 Phase 3A frozen tail 一致。 |
| Single-patch WENO7 SSPRK3 short run | 通过 | 运行至 `0.002`；每一步都有三个 closure label。 |
| Single-patch WENO7 RK4 short run | 通过 | 运行至 `0.002`；每一步都有四个 stage label 和 final closure。 |
| Serial IBM RK4 short run | 通过 | 运行至 `0.001`；含四个 stage 与 final closure，prescribed IBM 仍在 stage RHS path。 |
| RPI wall-boiling smoke | 有界执行 | 到达 `1e-5` 和 `2e-5` 的 pressure Eulerian step，未见新失败；长运行受交互命令时限停止。 |

完整 `cmake --build build --parallel 4` 在当前交互执行时限内没有结束。因此改为直接
编译修改过的 object，并执行最终 linker command；未遗留 compile 或 link diagnostic。
release 前仍应在无时限 CI/terminal session 中执行完整 serial WENO7 及完整 RPI
long regression。

已有的 `levelSetInjectionCase` SSPRK3 composition 仍会在时间积分开始前失败，原因是
Phase 3A 必需的 `StateBundle::equations` binding 缺失。这是既有 fail-fast
configuration 问题，不是本次 tableau 改动；本阶段没有加入 EOS fallback。

## 15. Stage-Time Unit Test

Phase 3B plan 中的 source-level truth table 已与统一代码逐项核对。运行时 stage
regression 覆盖了精确的 target-time sequence：short SSPRK3 density run 观察到
stage 1 / 2 / final；short RK4 density 和 IBM run 观察到 stage 1 / 2 / 3 / 4 /
final。integrator 从 committed bundle clock 计算这些 node，且不赋值
`StateBundle::time`。

没有新增独立 ODE CTest target：构造 synthetic `Field` 会引入另一条 test-only
lifecycle，反而无法覆盖 all-patch RHS 和 IBM target-time contract。两组 short
density case 是本次 consolidation 实际加入的回归覆盖。若测试基础设施今后提供最小
field factory，仍可补充 pure algebra test。

## 16. Single/Multi Integrator 一致性

现在 one-field `CompressibleAlgorithm` 与 complete-vector `MultiPatchAlgorithm`
调用完全相同的 `DensityBasedTime::advance`。4-rank Sod 验证 multi-patch 路径；
short WENO7 SSPRK3/RK4 case 验证 single-patch 路径。patch 数量只影响 shared stage
callback 所见 field 数量。

## 17. 剩余 Algorithm 差异

两个 algorithm class 有意保留各自的 pre/post-time composition 与 final tail。
特别是 final closure 和 `commitStep` order 未改变，application composition detail、
pressure-based stepping 及 flow/IBM correction placement 也未改变。`beginStep`
仍在 RK stage 之外，correction 仍在所有 stage 之后。

## 18. Architecture Allowlist

allowlist 保持为 10。Phase 3B 未增加 reverse dependency，也未涉及 scalar transport、
viscous ownership、models 或 linear backend。新实现在 `solver/algorithm`，通过
callback orchestration 既有 equation/discretization 工作。

## 19. 延后的 Final Closure / Commit Contract

以下问题明确留给 Phase 3C 或后续阶段：

- single/multi final closure 与 `commitStep` ordering；
- BoundaryPipeline redesign；
- pressure-based time integration ownership；
- scalar-transport equation ownership 与 multiphase temperature advance；
- model derivative 与 viscous ownership；
- Field responsibility split；
- IBM architecture、distributed monolithic KKT、HYPRE ownership 与 Equation DSL / AssemblyPlan。

Self review：RK coefficient **未改变**；stage time **未改变**；DensityBasedRHS、
boundary/halo order、canonical barrier、residual SUM、flow correction、final
boundary/commit order **均未改变**；删除 algorithm class **否**；新增
adapter/manager/context/strategy **否**。
