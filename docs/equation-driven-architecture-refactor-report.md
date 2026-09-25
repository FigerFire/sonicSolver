# Equation-driven 架构重构报告

日期：2026-09-15

## 1. 目标与边界

本次重构把启动和执行路径从物理类型驱动的 runner 分派，收敛为：

```text
CaseConfig
  -> physics/equation template composition
  -> ResolvedSimulationSystem
  -> formulation resolution
  -> solve-block Workflow::Plan
  -> execution builder
  -> generic runFlow
  -> ExecutionRuntime
```

本次没有修改 WENO/TENO、数值通量、RK/Euler、pressure correction、
Eulerian PIMPLE、IBM、EOS、phase change 或 MPI 数学语义。

## 2. Phase 0 基线

生产代码修改前完整 build 成功。记录的代表性一步结果为：

| Case | 数值路径 | 基线结果 |
|---|---|---|
| single WENO7 Sod | WENO7 + Steger-Warming + Euler | `dt=6.681531e-04`，`rho_min=0.125`，`p_min=10000` |
| Ghost IBM | WENO5 + Lax-Friedrichs + RK4 | `dt=4.392680e-04`，`rho_min=1.22383`，`p_min=101208` |
| variational IBM | WENO3 + Lax-Friedrichs + Euler | `dt=8.564235e-05`，3882 body dofs，`max|Ju-Us|=0` |
| 4-rank Sod | TENO5 + Steger-Warming + Euler | `dt=6.681531e-04` |

当前源码基线还包含三项既有 fail-fast 行为：homogeneous 六变量 EquationSet
不满足 PerfectGas 高阶通量 contract；level-set injection case 的热导率为零；
Eulerian 一步执行缺少显式 borrowed workspace。重构没有用 fallback、改输入或
放宽检查来隐藏这些问题。

## 3. 系统解析与 ownership

### Before

```text
CaseConfig
  -> Workflow::Request / Workflow::Kind
  -> physical workflow kind
  -> System::build
  -> physics-specific runner dispatch
```

Workflow 同时参与物理模板选择和执行计划选择，`SF_run.cpp` 根据
SingleFluid、Homogeneous、OneFluidInterface、EulerianEulerian 选择 runner。

### After

`SF_inspection.cpp` 先从 case 解析 `System::BuildRequest`，构造并验证
`ResolvedSimulationSystem`，随后才由 resolved system 生成 Workflow 和
ModuleGraph。`ResolvedSimulationSystem` 成为 unknown、equation、constraint、
solve block、formulation 和 execution requirement 的权威来源。

字段语义调整如下：

| Before | After | 语义 |
|---|---|---|
| `PhysicsStateKind` | `PhysicsTemplateKind` | 只记录模板来源，不参与运行时分派 |
| `system.flow` | `system.formulation` | primary-state formulation |
| `system.physics` | `system.templateOrigin` | composition provenance only |

新增的只读查询包括 `hasUnknown`、`hasEquation`、`hasEquationPrefix`、
`hasConstraint`、`hasSolveBlock`、`hasRequirement` 和 `requiresCapability`。
这些查询直接读取 resolved representation，没有复制第二套 Equation IR。

## 4. Physics、formulation 与 coupling

System builder 分别处理模板 composition、formulation、coupling、model
contribution 和 execution requirement。

Eulerian-Eulerian 由 phase 列表循环组装 N 份 phase equation pack，每相声明
phase mass、momentum、enthalpy 及其方程，并显式声明：

- shared pressure relationship；
- `sum(alpha.phase) = 1` volume-fraction closure；
- `S_EE_PIMPLE` solve block 中的方程与约束执行关系。

模板 composition 不再用“Eulerian 必然等于 pressure”来定义数学系统。
尚未实现的 formulation/equation 组合由 execution capability validation 拒绝。

## 5. Workflow 与 ModuleGraph

`Workflow::Kind` 和 `Workflow::Request` 已删除。`Workflow::Plan` 只保存有序的
solve-block 引用与 time-integrator 信息，由 `ResolvedSimulationSystem` 构造。

ModuleGraph 现在读取 equations、constraints、formulation、IBM enforcement 和
requirements。它不再读取 physical workflow kind。Eulerian 的普通代数约束
不会被误判成 IBM forcing；IBM module 由 resolved immersed enforcement 判断。

## 6. Generic FlowRunner 与 execution builder

应用层 `runFlow` 是唯一创建 `Time::Driver` 的位置。它统一处理：

- `INavierStokesStepper::advance`；
- accepted-step 检查；
- step/time 输出 hooks；
- after-step diagnostics；
- MPI finished agreement。

原 `runSingleOrHomogeneous` 与 `runEulerianEulerian` 已收敛为内部 equation
execution assembler。公开入口 `runSingleField` 依据 resolved solve block 选择
已经准备好的 executor，不依据 `templateOrigin` 或 physics booleans。
`SF_run.cpp` 不再分派物理类型 runner。

MultiPatch 仍保留独立的 state/domain 组装入口，但时间推进复用同一 `runFlow`。
它表示 patches、halo、canonical face 与 GlobalDof 的 execution domain，不表示
另一种 physics workflow。

## 7. Public API 与依赖变化

| 项目 | Before | After |
|---|---:|---:|
| `Workflow::Kind` | 1 enum | 0 |
| `workflowPlan.kind` consumer | 多处 | 0 |
| `Workflow::Request` | 1 struct | 0 |
| physics-specific public runner | 2 | 0 |
| application `Time::Driver` construction sites | 多条 lifecycle | 1 |
| resolved mathematical authority | 1 + workflow shadow decisions | 1 |
| runtime physics-template dispatch | 存在 | 0 |

最终搜索以下符号均为零：

```text
Workflow::Kind
workflowPlan.kind
runEulerianEulerian
runSingleOrHomogeneous
homogeneousEquationSet
PhysicsStateKind
```

架构守卫通过：580 个源文件，allowlisted dependency edges 保持 10，没有新增
依赖债务。重复 header basename 仍为既有 17 组，本次没有通过搬文件扩大范围。

## 8. Explain / check 验收

`sonicSolver check test/eulerianEulerianCase` 返回 `CONFIGURATION VALID`。
`sonicSolver explain` 从运行时使用的同一份 resolved system 输出：

- formulation：`pressureBase`；
- template origin：`eulerianEulerian`；
- water/air 各自 phaseMass、momentum、enthalpy unknown；
- 两相 continuity、momentum、enthalpy equations；
- incompressibility、shared pressure、volume-fraction constraints；
- `S_EE_PIMPLE` 的 9 行 algebraic system；
- `S_EE_PIMPLE` solve stage 和 execution requirements。

## 9. 数值与并行验收

最终完整 build 成功，共完成 304 个 Ninja 步骤，CLI `sonicSolver` 与 GUI
`sonicGui` 均成功链接。Ninja 输出过 `premature end of file; recovering`，但其
恢复后完整重建成功，未出现 compile 或 link failure。

| 验证 | 结果 |
|---|---|
| architecture checker | Passed，allowlist 10 |
| host CTest | 6/6 Passed |
| single WENO7 Sod | 与基线相同 |
| serial Ghost IBM | 与基线相同 |
| serial variational IBM | 与基线相同 |
| 4-rank Sod | Passed，`dt=6.681531e-04` |

受限 sandbox 内首次 CTest 的三个 MPI 用例因 PRTE 无权绑定 socket 而未启动；
host 环境复跑后 `sparseCanonicalCopy`、`distributedMinimalKKT`、
`distributedSurfaceSchur` 全部通过。这不是 assertion failure。

并行语义保持不变：physical state 使用 owner-to-replica COPY；shared face 使用
单一 canonical flux 后 COPY 并以正负号装配；residual/source 使用 SUM；
constraint row 继续使用唯一 GlobalDof ownership。

## 10. Self review

| 检查 | 结果 |
|---|---|
| 新增第二套 Equation/System authority | No |
| physics type 决定 runner | No |
| density/pressure 决定 physics template | No |
| MultiPatch 作为 physics workflow | No |
| 修改 stage timing / boundary order / halo order | No |
| 修改 flux / residual / pressure correction | No |
| 修改 canonical COPY / GlobalDof SUM | No |
| 使用 hidden fallback | No |

## 11. Deferred issues

- generic EOS 与高阶 characteristic flux 的 thermodynamic interface；
- homogeneous multiphase execution capability；
- level-set thermal-property 输入有效性；
- Eulerian borrowed workspace 的显式 execution contract；
- scalar transport 与 viscous numerics 的最终模块归属；
- Field 剩余职责、IBM architecture、distributed KKT 与 HYPRE ownership；
- Equation DSL / AssemblyPlan 的长期 executable authority 收敛。
