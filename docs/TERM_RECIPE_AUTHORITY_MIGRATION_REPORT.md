# Term Recipe Authority Migration Report

日期：2026-09-21

## 1. Scope

本轮完成 density-formulation production path 的 term recipe authority 与
physical term-presence 迁移：

```text
Equation System
    + NumericalRecipeSet
    -> NumericalCompiler
    -> CompiledNumericalSystem
    -> density RHS / boundary / mesh runtime
    -> existing numerical kernels
```

已迁移 convection、core-fluid diffusion、gravity/MRF/wallHeat source 的
binding、halo、workspace 和 provider requirements。未修改 WENO/TENO、
Rusanov/Steger-Warming、中心扩散、RK、PISO、IBM/KKT 或 MPI 数学。

Pressure/Eulerian specialized term execution 尚未迁入严格 recipe coverage；
它们继续使用既有 numerical provider，以保持本轮边界和数值行为。

整体状态：**architecture implementation complete；full Sod numerical
validation blocked by a pre-existing reconstructed-Rusanov defect**。按任务规则，
没有修改该数值公式，也没有降低 CFL 或增加 clamp。

## 2. Before architecture

旧 density execution 同时存在多套 authority：

```text
SolverConfig
  -> ConvectionScheme + FluxSplitter
  -> AssemblyContext carries SolverConfig
  -> divDispatch reinterprets selection

scheme string -> requiredGhostLayersForConvection() -> halo

Equation always contains diffusion
  + runtime viscousEnabled / transport pointer -> execute or skip

Equation contains source
  + SourceTerm::Sp loops config.enabled -> select again
```

因此 equation description、halo selection 和 runtime execution 可以漂移。

## 3. After architecture

```text
Input / legacy compatibility input
    -> one-time built-in recipe resolution

ExecutableEquationSystem + NumericalRecipeSet
    -> NumericalCompiler
    -> CompiledNumericalSystem
         BoundTerm
         requiredHaloWidth
         workspaceRequirements
         providerRequirements
    -> DensityBasedRHS
    -> existing WENO/TENO + flux kernel
```

`AssemblyContext` 不再携带完整 `SolverConfig` 来选择 convection、flux 或
diffusion presence。它只接收已经编译的 convection/diffusion recipe、source
bindings，以及 kernel 所需物理参数。

`sonicSolver explain` 新增 `COMPILED NUMERICAL SYSTEM`，逐项显示 equation、
term ordinal、recipe、temporal role、halo、workspace 与 provider requirement。

## 4. Physical term presence

### Before

`addDensityBasedFluid()` 总是声明 momentum/energy diffusion；runtime 再用
`viscousEnabled || transport != nullptr` 静默决定是否跳过。

### After

composition boundary 把旧 `viscousEnabled` 和已有 transport/model composition
转换一次为 `BuildRequest::fluidDiffusion`：

- Euler/inviscid raw equations 不包含 diffusion；
- viscous/model-required raw equations明确包含 diffusion；
- NumericalCompiler 对存在的 diffusion 要求 recipe；
- density system 上选择未使用的 diffusion recipe 会 fail fast；
- `laplacianDispatch` 不再接收 enabled flag。

这是 ownership/API 迁移，没有改变数组布局或扩散公式。

## 5. TermRecipe catalog

当前 catalog 中每一项均映射现有 production kernel，没有 placeholder：

- Convection：WENO3、WENO5、TENO5、WENO7 与 Steger-Warming、Rusanov、
  Lax-Friedrichs、Roe、Lax-Wendroff 的既有组合；
- Diffusion：`central2Explicit`、`central4Explicit`；
- Source：`gravityExplicit`、`mrfExplicit`、`wallHeatExplicit`。

本轮 primary closure 实测 `weno7Rusanov` 与 `weno7Steger`。Recipe 是不可逐项
修改的 value/spec；native input 不能重新组合 scheme/flux flags。旧输入只在
IO boundary 解析一次并映射到 catalog。

## 6. NumericalCompiler

NumericalCompiler 实际遍历 executable core-fluid definitions：

- `div(...)` -> selected convection recipe；
- `diffusion(...)` -> selected diffusion recipe；
- `source(gravity/MRF/wallHeat)` -> 对应 source recipe；
- `source(zero)` 保持为方程零右端语义，不生成 provider。

编译结果汇总：

- WENO3 halo 2，WENO5/TENO5 halo 3，WENO7 halo 4；
- face flux、residual、characteristic reconstruction workspace；
- `term.convection.*`、`term.diffusion.*`、`term.source.*` provider requirements；
- explicit TimeRecipe topology compatibility。

缺少 convection/diffusion recipe、density equation 不含 diffusion 但选择
diffusion recipe、未知 recipe 名都会 fail fast。

## 7. Removed authority

已删除或停止使用：

- `requiredGhostLayersForConvection(string)`；
- density runtime 对 `numerics.convection/flux/reconstruction` 的再次选择；
- equation assembly 的 `viscousEnabled` skip；
- diffusion dispatcher 的 enabled 参数；
- `SourceTerm::Sp` 对 `config.enabled` 的第二次 physics selection；
- `AssemblyContext` 中仅用于数值选择的完整 `SolverConfig`；
- native input 中 `terms` 与 legacy convection/diffusion 同时出现的歧义。

Native template 现在生成 `forwardEuler`，不恢复 `Euler` alias。

## 8. Remaining Legacy

- Pressure/Eulerian numerical term execution 尚未使用 strict compiled-term
  coverage；PISO/PIMPLE/IBM provider 保留现有实现。
- Level-set、turbulence transport、phase-change 与 Eulerian phase term 的
  provider composition 仍有 EquationId/专用 executor mapping。
- `divDispatch` 内保留一个对已经 resolved recipe 的局部 kernel switch；它不再
  读取用户输入或 raw SolverConfig。
- Compatibility input 仍接受旧 convection/flux/diffusion 字段，但只在 IO
  boundary 映射一次。
- 2026-09-20 17:44 的并发 config-type ownership 拆分使当前 architecture
  checker 报告 7 条新 dependency edges；详见 Regression。它不属于本轮
  numerical authority 修改，未通过扩大 allowlist 隐藏。

## 9. Unsupported

以下仍显式 Unsupported，未添加 fallback：

- implicit TermRecipe；
- backwardEuler、SDIRK2、ARK、IMEX；
- runtime-defined recipe；
- runtime-defined mathematical term；
- generic EOS high-order characteristic flux；
- pressure/Eulerian 的 strict generic term lowering。

## 10. Regression

| Check | Result |
|---|---|
| Normal build `cmake --build build --parallel 4` | PASS |
| Test build `cmake --build build-tests --parallel 4` | PASS |
| `pisoArchitecture` | PASS |
| `explicitStageMathematics`（含 RK4 stage mathematics） | PASS |
| `termRecipeAuthority` | PASS |
| `cliGeneratedRecipe` | PASS；CLI template -> production reader |
| `pisoNumericalRegression` | PASS |
| `numericalFluxContract` | FAIL；修改前已存在的 reconstructed-Rusanov constant-stencil mismatch |
| Architecture checker（本轮 guard 加入后、并发拆分前） | PASS，allowlist=10 |
| Architecture checker（最终共享工作区） | FAIL，7 条并发 config ownership dependency edges；allowlist 未增加 |
| Rusanov one-step serial | PASS；迁移前后 VTS SHA-256 完全相同 |
| Steger one-step serial | PASS；迁移前后 VTS SHA-256 完全相同 |
| Rusanov vs Steger dispatch | PASS；两个 one-step hash 不同 |
| Full configured Rusanov Sod | FAIL at `t=0.06478187`；旧输入复跑在同一步同样失败 |
| 4-rank existing Sod one-step | PASS；4 partitions，`dt=6.681531e-04` |
| IBM config validation | PASS，10/10 built-in method cases |
| IBM production one-step | PASS，10/10 built-in method cases |

IBM one-step 覆盖 Ghost/ILW、Peskin、FTS、BP、prescribed/self-propelled
fractional DLM、fully implicit prescribed/self-propelled KKT 与 augmented
Lagrangian。代表性诊断：FTS `max|Ju-Us|=1.9613e-12`、CG 9 iterations；
fully implicit prescribed KKT 2 iterations、residual `1.07973e-09`。

MPI 回归在 host socket 环境执行。4-rank case 仍使用 4 个 partition，未修改
canonical face COPY、owner/neighbour residual assembly 或 GlobalDof SUM。

## 11. Sod result

### One-step exact regression

| Recipe | Time | dt | Before SHA-256 | After SHA-256 | Status |
|---|---:|---:|---|---|---|
| `weno7Rusanov` | 6.681531e-04 | 6.681531e-04 | `844036d46ff1119aa7cf3340328621f008b42b958cb1f2eb457148930ece0dad` | same | bit-identical |
| `weno7Steger` | 6.681531e-04 | 6.681531e-04 | `1ab9ccf4cfb3083459de82fa2312e9eaafc0dbe4be896cdf354d9aafd9d55d07` | same | bit-identical |

### Full configured run attempt

- Recipe: `weno7Rusanov`
- Time method: `forwardEuler`
- Original end time: `5.0`
- Completed steps before fail-fast: 149
- Last committed time: `0.06478187`
- Last dt: `3.889387e-04`
- Failure: `PerfectGasEOS produced an invalid state`

迁移后 native recipe 与迁移前 legacy `WENO7 + Rusanov` 运行具有相同 dt/state
序列、相同失败点，且 `t=0.05` artifact hash 同为：

```text
016a301439972238965d5dd09169f9c21a87e9ef323dcf31ea50c9573ec82c91
```

`t=0.05` 有限场统计：

| Field | Min | Max |
|---|---:|---:|
| rho | 0.1249965240 | 1.000314316 |
| p | 9999.538674 | 100044.1838 |
| |U| | 6.39394e-10 | 361.3572902 |
| rhoE | 24998.84669 | 250110.4666 |

`numericalFluxContract` 同时报告 constant stencil 上 reconstructed Rusanov 与
physical flux 不一致（例如 mass flux 36.0551 vs 36）。这是 full-run stability
未闭环的真实 numerical blocker；本轮禁止修改该公式，因此没有更新 tolerance
或 baseline。

## 12. Authority table

| Question | Authority after migration |
|---|---|
| equation 是否有 convection？ | Equation System |
| equation 是否有 diffusion？ | Equation System / physics composition |
| convection 是否用 WENO7/Rusanov？ | TermRecipe |
| halo width？ | TermRecipe requirements -> CompiledNumericalSystem |
| Forward Euler？ | TimeRecipe |
| StageLoop 次数？ | TimeRecipe -> CompiledSolvePlan |
| 实际 flux kernel？ | compiled recipe binding 消费的既有 provider/kernel |
| runtime 是否重读 scheme string？ | 否 |
| runtime 是否按 viscousEnabled 跳过 diffusion？ | 否 |
| Source provider 是否重新决定 MRF presence？ | 否 |

## Completion status

A-I、K 已满足；L 的现有 PISO、IBM 和 MPI smoke/invariants 未回归。
J（完整 Sod 到原 end time）未满足，原因是已通过 before/after 同步复跑证明的
既有 reconstructed-Rusanov numerical defect。因此本轮不能把整体状态标为
fully validated，也不能在架构迁移中偷偷修复或掩盖该缺陷。
