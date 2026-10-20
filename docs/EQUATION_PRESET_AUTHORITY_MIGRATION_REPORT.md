# Equation Preset Authority Migration Report

日期：2026-09-21
阶段状态：**Complete（以标准短时 Sod `t = 0.2` 作为冻结回归基线）**

## 1. 范围与结论

本阶段依次完成了两个 Gate：

1. 修复 characteristic reconstruction 不能保持常量守恒状态的问题，并建立可重复的 WENO7 + Rusanov + Forward Euler Sod 数值基线；
2. 将 density single-fluid 的核心 Navier–Stokes 方程定义从巨型 `SF_systemBuilder.cpp` 移到唯一的 built-in equation preset。

本阶段只改变：

- Gate A：修正经 contract test 证明错误的 PerfectGas 特征左特征向量系数；
- Gate B：改变方程 composition 的 ownership 和调用路径。

未改变 RK 系数、stage time、CFL、Rusanov 定义、WENO/TENO 标量重建、边界/halo 顺序、canonical face COPY、GlobalDof SUM、IBM/KKT 数学或压力修正算法。

## 2. Gate A：Reconstructed-Rusanov 数值闭环

### 2.1 根因

错误位于：

```text
src/methods/numerics/convection/SF_riemann.cpp
Riemann::buildCharacteristicMatrix()
```

原 acoustic left eigenvector 的密度分量使用：

```text
b2 * |u|^2 ± un / (2c)
```

其中：

```text
b2 = (gamma - 1) / (2 c^2)
```

正确 kinetic coefficient 是：

```text
(gamma - 1) |u|^2 / (4 c^2)
= 0.5 * b2 * |u|^2
```

错误系数使 `L * R != I`，因此即使输入 stencil 完全为常量，`R * W` 也不再等于原守恒状态。Rusanov 本身对相同左右状态是一致的，偏差产生在 characteristic projection/inverse projection，而不是 Rusanov flux。

### 2.2 修改前诊断

测试状态：

```text
Q0 = [1.2, 36, -2.4, 1.2, 253855.50000000006]
```

修改前 characteristic round trip 得到约：

```text
R * (L * Q0)
= [1.201837369, 36.055121074, -2.403674738,
   1.201837369, 254399.3314]
```

因此常量面的质量通量从正确值 `36` 变为约 `36.0551`。

修改后：

```text
L * Q0
= [0.3428571428571429, 0, 0,
   0.4285714285714285, 0.4285714285714285]

R * (L * Q0)
= [1.2, 35.99999999999999, -2.4,
   1.2, 253855.50000000003]

max |L*R - I| = 1.76e-16
```

对应物理通量和 equal-state Rusanov 通量均为：

```text
[36, 102405, -72, 36, 10655415.000000002]
```

### 2.3 Contract tests

`test/test_numericalFlux.cpp` 现在覆盖：

- WENO3、WENO5、TENO5、WENO7 的常量 conservative stencil preservation；
- Rusanov、Steger-Warming、Roe、Lax-Friedrichs 的 equal-state consistency；
- characteristic WENO reconstruction 后的 Rusanov、Steger-Warming、Lax-Friedrichs 常量面通量一致性；
- 配置中的 `Rusanov` 枚举解析与输出。

所有比较继续使用 `1e-12` 相对尺度容差；未放宽容差。

### 2.4 冻结 Sod baseline

保留原 `test/Sod/sodCase_weno7` 的 `endTime = 5.0`，并新增标准短时回归：

```text
test/Sod/sodCase_weno7_t0p2
WENO7 + Rusanov + Forward Euler
endTime = 0.2
```

冻结结果：

| 指标 | 结果 |
|---|---:|
| step count | 454 |
| final time | 0.2 |
| rho min/max | 0.25924657327767214 / 0.7622215092077242 |
| p min/max | 13037.114603387012 / 68970.12209014251 |
| \|U\| min/max | 95.968361480711 / 391.43134094748274 |
| mass integral | 98.67807117165903 |
| momentum integral | [23619.052634273423, 1.907691669240853e-08, 0] |
| energy integral | 21770158.12505413 |
| final VTS SHA-256 | `8bf55a1c4973f3a2322aeb8785305bf0589cc9eb4cf89505203309f0ccba0768` |

`baseline.json` 冻结全部 454 个 `dt`，回归脚本还验证 step count、final time、有限性、`rho > 0`、`p > 0`、范围、积分和输出 SHA。

原 `endTime = 5.0` case 未删除。修复后它从原先约 `t = 0.06478187` 的失败推进到约 `t = 0.3758227`，随后仍因 PerfectGas 非法状态 fail fast。这是 WENO7 + Forward Euler 在长时间运行中缺少 positivity-preserving capability 的独立限制。本阶段没有通过降 CFL、clamp、density/pressure floor 或一阶 fallback 隐藏它；稳定数值回归的 required end time 明确采用标准 Sod 时间 `t = 0.2`。

## 3. Gate B 修改前的数据流

修改前 density single-fluid core equations 的真实 authority 链为：

```text
case input
  -> CaseConfig
  -> Application::inspectCase()
  -> BuildRequest::templateOrigin + BuildRequest::fluidDiffusion
  -> System::build()
  -> addDensityBasedFluid()
  -> rho/rhoU/rhoE + Mass/Momentum/Energy DSL
  -> NumericalCompiler
  -> SolvePlanner
  -> runtime
```

`SF_systemBuilder.cpp` 同时决定 physics equations、transformations、solve policy、numerical compile、runtime requirements 和 capability reporting。`fluidDiffusion` 还是 application 写入、builder 再读取的第二个 physical term presence flag。

## 4. 修改后的数据流

现在 density single-fluid 使用：

```text
typed CaseConfig
  -> resolve SingleFluidPresetSpec
  -> Preset::installSingleFluid()
  -> RawEquationSystem
  -> ExecutableEquationSystem
  -> NumericalCompiler + NumericalRecipeSet
  -> CompiledNumericalSystem
  -> CompiledSolvePlan
  -> runtime
```

`Preset::installSingleFluid()` 是以下数学内容的唯一注册点：

- unknowns：`rho`、`rhoU`、`rhoE`；
- equations：Mass、Momentum、Energy；
- physical term presence：viscous spec 下 Momentum/Energy 含 diffusion，inviscid spec 下不含。

Preset 只决定 WHAT。它不读取或选择 WENO、TENO、Rusanov、Steger-Warming、diffusion stencil、Euler、RK 或 execution order。

Preset 放在现有 `solver/system` composition 层，而没有放入 `solver/equation`：它需要写入 `SystemCompositionBuilder`。若放入 equation 层，会形成 `equation -> system` 的反向依赖。没有为此新增 Manager、Adapter、Context、Facade 或 ServiceLocator。

## 5. 删除的 authority

实际删除或停止 production 使用的 authority：

- `addDensityBasedFluid()`；
- `BuildRequest::fluidDiffusion`；
- `SF_systemBuilder.cpp` 内 density core NS 的 `Equation::ddt/div/diffusion` DSL；
- density single-fluid 根据 `templateOrigin` 隐式生成 core equations 的路径；
- numerical diffusion recipe 根据 `BuildRequest` flag 决定 term presence 的路径。

NumericalCompiler 现在从 executable equations 的实际 `TermKind::Diffusion` presence 得出是否需要 diffusion recipe。Source term 仍由 Equation System 注册，再由 NumericalCompiler lower；没有重新引入 runtime physics selection。

## 6. Authority tests 与架构守卫

新增或加强的验证包括：

1. contribution provenance 必须包含 `builtin.singleFluidNavierStokes`；
2. inviscid preset 不含 diffusion；
3. viscous preset 的 Momentum/Energy 含 diffusion；
4. 同一 raw Equation System 可分别绑定 WENO7-Rusanov 和 WENO7-Steger，而不重建数学方程；
5. `SF_systemBuilder.cpp` 不得再次出现 density core `rhoU` DSL；
6. 禁止恢复 `addDensityBasedFluid()`；
7. 禁止恢复 `BuildRequest::fluidDiffusion`。

## 7. BuildRequest 与文件收缩

| 指标 | Before | After |
|---|---:|---:|
| `SF_systemBuilder.cpp` LOC | 1223 | 1213 |
| 独立 preset implementation LOC | 0 | 62 |
| 独立 preset header LOC | 0 | 18 |
| BuildRequest data members | 17 | 17 |
| BuildRequest boolean feature flags | 11 | 10 |
| density core equation construction functions in giant builder | 1 | 0 |
| density core NS DSL occurrences in giant builder | 1 family | 0 |

BuildRequest 总字段数保持 17，是因为删除 `fluidDiffusion` 后加入了 typed、optional 的 `SingleFluidPresetSpec`。它不再用一个布尔 feature flag让 giant builder 重建方程；preset spec 是对选定 mathematical preset 的显式输入。

新文件只有 62 行实现，包含单流体 preset 的全部 unknown/equation 注册，没有把旧 giant builder 整体搬到另一个大文件。

## 8. Dependency contraction

阶段开始时 architecture checker 报告 7 条新增 violation：

| Source layer | 旧 target ownership |  offending type/include |
|---|---|---|
| discretization | algorithm | time recipe |
| core/config | models/IBM | IBM config value types |
| core/config | models/turbulence | turbulence config value types |
| core/config | solver/pressure | pressure policy value types |
| core/config | solver/time | time recipe value types |
| core/config | solver/discretization | term recipe value types |
| core/config | solver/linearAlgebra | linear solver config value types |

这些类型都是只保存配置值的 DTO，现移到最低公共配置层：

```text
src/core/config/types/
  SF_ibmConfigTypes.h
  SF_turbulenceConfigTypes.h
  SF_pressureConfigTypes.h
  SF_timeRecipe.h
  SF_termRecipe.h
  SF_linearSolverConfigTypes.h
```

结果：

| 指标 | Before | After |
|---|---:|---:|
| 新 architecture violations | 7 | 0 |
| allowlisted dependency edges | 10 | 10 |

没有扩大 allowlist、修改 baseline、屏蔽目录或删除 checker 规则。

## 9. Native IO 边界

新 preset 只消费 typed `SingleFluidPresetSpec`，不读取磁盘、不解析 dictionary，也不访问 runtime loop。`Model::Description` 中的 equation/algorithm/thermodynamic composition 已在 application compatibility binding 中直接写入 typed `EquationCompositionConfig`，不会通过 dictionary text 再解析这些 semantic objects。

完整 native case 仍有剩余 compatibility debt：`ModelLoader::build(Model::Description)` 对 control、solver、numerics、field 与若干 built-in model 参数仍调用 `CaseAdapter::build()`，其中会构造内存 dictionary text 并复用 legacy `readCase()`。因此本阶段只建立了 preset 的 typed consumption boundary，没有宣称完成全量 native IO 迁移。该 bridge 只能产生 `CaseConfig`，不控制 runtime 或重新选择 preset。

## 10. Remaining Legacy

以下内容保持现状，仍由 giant builder 或 compatibility path 参与：

- pressure single-fluid composition 与 pressure execution；
- Eulerian multiphase template；
- level-set / one-fluid interface composition；
- turbulence equations/contributions；
- phase change；
- IBM contribution、constraint 与 KKT；
- homogeneous/legacy mixture bridge；
- native input 中 control、solver、numerics、field 和部分 model 的 dictionary compatibility translation；
- `BuildRequest::templateOrigin`，仍服务上述 legacy paths 和 reporting。

本阶段没有新增 `EulerianEulerianPreset`。未来 Eulerian 应由独立 EquationSystem instances 先各自构造，再 resolve coupling。

## 11. 回归结果

### Build 与静态检查

| Check | Result |
|---|---|
| `cmake --build build --parallel 4` | PASS |
| tests build | PASS |
| `python3 tools/check_architecture.py` | PASS；allowlist = 10 |
| scoped `git diff --check` | PASS |

### CTest

`ctest --test-dir build-tests --output-on-failure --parallel 1`：**7/7 PASS**。

```text
numericalFluxContract
pisoArchitecture
explicitStageMathematics
termRecipeAuthority
cliGeneratedRecipe
sodRusanovNumericalRegression
pisoNumericalRegression
```

### 数值/并行 smoke

| Case | Result |
|---|---|
| WENO7 + Rusanov + Forward Euler Sod, t=0.2 | PASS，454 steps，冻结全 dt 序列和 VTS SHA |
| WENO7 + Steger-Warming one step | PASS，dt = 6.681531e-4，rho/p 保持正值 |
| 4-rank TENO5 + Steger-Warming Sod one step | PASS，4 partitions，canonical runtime initialized |
| RK4 stage mathematics | PASS |
| PISO architecture + corrected numerical regression | PASS |
| 10 built-in IBM serial one-step smoke | PASS |

IBM smoke 覆盖：Ghost Cell、Peskin Original、DFM explicit、DFM fractional prescribed/self-propelled、DFM implicit prescribed/self-propelled、DFM augmented Lagrangian、Velocity Forcing FTS、Velocity Forcing BP。此次没有修改 IBM，因此按 IBM validation policy 使用一阶段 smoke 检查跨路径污染，没有把 Ghost test 当作 forcing/KKT correctness 的证明。

## 12. Numerical、execution 与 MPI invariants

| Invariant | Changed? |
|---|---|
| TermRecipe / NumericalCompiler authority | No |
| TimeRecipe / SolvePlan authority | No |
| RK coefficients / stage times | No |
| Boundary / halo order | No |
| Canonical face owner -> COPY | No |
| Residual/source/load -> SUM | No |
| GlobalDof ownership/layout | No |
| IBM/KKT numerical implementation | No |

Gate A 唯一数值公式修改是把错误的 acoustic left eigenvector kinetic coefficient 改为数学上一致的 `0.5 * b2 * q2`。Gate B 的 preset 迁移对冻结的 density baseline 为 bit-identical。

## 13. Self Review

```text
numericalFluxContract passed?                       Yes
stable WENO7/Rusanov Sod t=0.2 passed?              Yes
stable Sod baseline frozen?                         Yes
density single-fluid core equations have one owner? Yes
giant builder still contains density NS DSL?        No
BuildRequest::fluidDiffusion remains?                No
TermRecipe authority changed?                        No
TimeRecipe / SolvePlan authority changed?            No
pressure/Eulerian/IBM regression introduced?         No
architecture allowlist increased?                    No
positivity clamp or numerical fallback added?        No
```

下一步应继续按 equation family 逐个收缩 Remaining Legacy。原 `endTime = 5.0` WENO7 + Forward Euler case 的长期正性能力需要独立数值阶段处理，不能在 architecture refactor 中加入隐式 safety net。
