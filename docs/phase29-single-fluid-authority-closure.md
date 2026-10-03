# Phase 29 — Single-Fluid WHAT / HOW / HOW NUMERICALLY Authority Closure

> 历史报告：后续 Equation Registry / Ordered Execution 迁移已替换本报告中的 FormulaGroup、ProgramStep 和 pressure schedule authority。当前状态见 [迁移报告](equation-execution-migration.md)。

```text
WHAT
  authored Formula / FormulaGroup
HOW
  ExecutionProgram.root：Formula -> explicit Output + structured controls
HOW NUMERICALLY
  subject-bound TemporalMethod + EquationMethod + Formula occurrence bindings
                         ↓
  CompiledExecutionProgram（method references / source inputs / workspace contract）
                         ↓
  Formula AST spatial compilation + compiled numerical providers
                         ↓
  CompiledSolvePlan
                         ↓
  requiredOperations -> resolved providers -> PlanExecutor / OpRegistry
                         ↓
  existing numerical kernels -> committed StateBundle
```

日期：2026-09-30。本轮修改基于已有 dirty 工作树继续演化，未 reset 前轮改动，未上传 GitHub。根目录没有 `ARCHITECTURE.md`，实际架构文件是 `src/ARCHITECTURE.md`，已更新。报告保留在本地忽略的 `docs/`。

## 1. 已删除的生产 authority

| Before | After | 主要文件 |
|---|---|---|
| 常密度 pressure PlanFragment 直接持有宏观 schedule；布尔参数选择旧/新引用 | `pressureProgram()` 生成唯一 structured HOW；删除 `bindConstantDensityFormulas` 分支 | `solver/system/SF_pressureCoupling.*`、`SF_systemBuilder.cpp` |
| EquationMethod 生成 FormulaCall.mode；Plan 携带执行模式 | 生产引用为 `CompiledMathRef {formula,target}`；方法编译自己的 numerical fragment | `core/system/SF_solveProgram.h`、`SF_methodObjects.cpp` |
| BoundTerm / Equation::Term / EquationRole 决定 single-fluid spatial binding | `CompiledSpatialBinding` 由 compiled method 的 residual references/source inputs、Output 和 Formula occurrence 生成 | `SF_numericalCompiler.*`、`SF_numericalSystem.h` |
| 源项 provider 根据 EquationRole、旧 Term 匹配 | provider 根据 AST operator、明确 Output、selected recipe 匹配；捕获的参数仍由 model-owned kernel 持有 | `core/system/SF_sourceProvider.h`、`models/physics/{gravity,mrf,heat}/` |
| density Formula 来自 legacy definition；两份数学依靠 equality guard 同步 | builtin 直接注册 authored AST；源项扩展修改 AST；删除 equality guard | `SF_singleFluidPreset.cpp`、`SF_presets.cpp`、`SF_equationContribution.cpp`、`SF_singleFluidStepper.cpp` |
| fused backend 再用 AssemblyPlan.contains() 判定 convection/diffusion/source 是否运行 | 删除 backend 的 Definition 指针、AssemblyPlan、contains 和旧构造参数；只消费编译后的 recipe/source kernels | `solver/equation/compressible/SF_compressible.*` |
| numerical compiler 再解析 TemporalMethod | composition 解析并编译一次；同一冻结结果供 method fragment、numerics、Plan 和 stage kernel 使用 | `SF_systemBuilder.cpp`、`SF_methodObjects.*`、`SF_numericalCompiler.cpp` |
| compiler 自动补 Output workspace provider | EquationMethod 必须声明 workspace provides；requires 按编译顺序验证 | `SF_methodObjects.cpp` |

这些旧入口从迁移路径中删除，不是新增 forwarding wrapper。`Equation::Definition` 的 legacy metadata mirror 和未迁移系统仍存在，但不再控制迁移后的 single-fluid spatial mathematics。旧通用 FormulaCompiler 的 FormulaMode API 仅保留在明确的能力证明/旧接口范围。

## 2. 实际调用链与职责

`SystemBuilder` 先组合数学系统，生成 structured HOW 和 EquationMethod bindings，解析 TemporalMethod，再编译 EquationMethod。空间编译器只读 compiled method 声明的 residual references 和 source-math inputs，遍历相应 Formula AST；它不重新按方法名推断 target 或流动物理，也不读取旧 Definition 的项。

常密度压力的五个 step：

| HOW | EquationMethod | numerical fragment / workspace |
|---|---|---|
| `E_MOMENTUM_PREDICTOR -> U` | PressureMomentum | MomentumAssemble、MomentumSolve；显式输入 E_MOMENTUM Formula；提供 HbyA/rAU/faceResponse/correctedFlux |
| `E_PRESSURE -> pPrime` | PressureCorrection | boundary prepare、assemble、solve；产出 pressureCorrection |
| `F_VELOCITY_CORRECTION -> U` | VelocityCorrection | 保留原 velocity correction 与 boundary callback |
| `F_PRESSURE_UPDATE -> p` | PressureUpdate | 保留 pressure update/gauge/boundary callback |
| `F_FLUX_CORRECTION -> phi` | FluxCorrection | 保留原 Rhie–Chow/canonical flux 运算；产出 pressureFaceFlux |

`pressureProgram()` 声明 outer/pressure/non-orthogonal loops、fixed-time begin/end、relaxation、flux restore、convergence 和 commit。EquationMethod 仅产生局部操作片段；PressureOperators 仅执行原 numerical helpers。没有新的 SimpleSolver/PisoSolver/PimpleSolver。

PISO 顺序仍是 prepare → predictor → pressure loop（non-orthogonal pressure passes → velocity → pressure update → flux → correction commit）→ step/time commit。SIMPLE/PIMPLE 仍先冻结 physical time base，进入 outer iteration，再执行 predictor/corrections、relaxation、flux restore、convergence 和 iteration end；PIMPLE 保留 inner pressure loop。loop count、early-exit signal、minimum iterations 与 commit 位置保持原语义。

Density 的 `Flow -> Q` 由 ConservativeResidual 产生明确公式/Output references，TemporalMethod 生成已有 stage topology。single/multi patch 使用同一 Plan/RHS 链；没有恢复两套 density lifecycle。

## 3. Spatial binding 与 fail-visible 能力边界

`TermMatchContext` 只提供 Formula、AST expression、occurrence、Output、EquationMethod 和 selected recipe。`CompiledSpatialBinding` 保存 Formula/执行 subject、Output、方法、occurrence、operator/operands、side、recipe/provider/owner 和 compiled kernels；不携带旧 TermKind/EquationRole。

优先级为：global operator default < formula default < specific occurrence override。生产 TermProviderCatalog 与通用 FormulaOperatorCatalog 共用 `selectFormulaBinding()`，选择结果不受注册顺序影响。最高优先级重复、选定 provider 缺失或不兼容、地址没有被选中的 Formula occurrence 消费，都明确失败；不会自动换一个 provider。

本轮没有实现任意公式的 fused lowering。已有 transport backend 支持 additive、unit-coefficient terms 及已支持的 operands/side；未知 spatial operator、未支持的正负号/系数算术、未绑定 residual data、错误/重复 ddt Output 或 term side 明确 fail-fast。fused provider 同时验证实际支持的 mathematical operands：massFlux/momentumFlux/energyFlux 与 mu/U、conductivity/T；不能把未实现的 userFlux 当作 NS flux。编译器不理解微积分，也不静默忽略 AST 数学。

声明了方程但没有对应 compiled HOW ownership 时，runtime report 明确标记 Unsupported，并指出缺失 equation step。这样压力湍流案例仍可 explain 缺少的 `closure.turbulence` / `equation.turbulence-transport`，但不会被伪装成只计算 stress 的 Runnable case。

## 4. Temporal / Workspace / State ownership

TemporalMethod 绑定到 transient subject，例如 Flow。compileExecutionProgram 接收已经冻结的 CompiledTimeRecipe 和同一个方法对象；NumericalCompiler 不再解析第二个时间方法。配置的时间 ID 保留为选择及一致性检查数据，不成为另一组 runtime stage 系数。

StateBundle/Field 继续拥有原 physical storage 和 clock；PatchWorkspace、pressure workspace、RK storage 仍由原 solver execution 持有。本轮未改变 conservative array、flux/residual stride、halo buffer、GlobalDof、lambda layout，也未增加 Field/Q 或 numerical workspace 的完整副本。

workspace validation 位于 method compilation。PressureMethods 声明具体 numerical workspace；DirectEvaluation/LinearEquation 明确声明自己产出的 workspace Output。编译器不得自动把任意 Output label 视为已提供资源。

## 5. Public API 与依赖

- 常密度 pressure composition 改为 `pressureProgram(request)`；旧路径明确命名 `legacyCouplingPlanFragment(request)`，不再有布尔切换。
- NumericalCompiler 输入改为 CompiledExecutionProgram 和已编译时间方法，输出 operators，删除 BoundTerm/terms API。
- Compressible fused equation backend 删除 Definition 构造参数、definition() view、AssemblyPlan 成员和 contains()。
- PressureOperators、ConservativeRHS、ProviderResolver 消费 Formula-derived binding。
- 新旧状态/API 未加 forwarding wrapper。FormulaOperatorBinding 属于 solver numerical catalog；Formula WHAT 不保存 numerical provider 选择。
- Application → compiled system → Plan/algorithm → existing numerical kernels 的方向保持；core 未新增 solver implementation 依赖。architecture checker 保留 2 条原有 methods → solver allowlist。

## 6. 数值和 MPI 不变量

未修改 WENO/TENO、Steger–Warming/Rusanov/LF、特征结构、RK 系数和 stage times、CFL、fixed-time predictor 公式、pressure matrix、gauge、Rhie–Chow 或 HYPRE 参数。

boundary → halo → Ghost/ILW → spatial reads 的顺序不变。各 patch candidate flux 完成后仍只有一次 canonical finalize；canonical F* 仍 COPY，随后 +/- residual assembly。GlobalDof residual/source/load 仍 SUM。residual sign、workspace clear、source clear 和 diagnostics-after-commit 时机未改。没有新的 fallback/clamp、算法降阶或 tolerance 调整。

## 7. 验证

最终 Debug binary 与 build/CTest/checker/A-B 日志另保存在 `/private/tmp/sonic-phase29-baseline-20260930`，solver SHA-256 为 `136558b2b392f02bbd990f43820a7418940f0fde488e4bcc6f8bca0b0fe12a2d`。这用于后续回归，不是新的运行时 state。

迁移前 Phase 28D.3 的独立记录：完整构建、host CTest 27/27、2 条 architecture allowlist。最终验证使用 build（Debug）和 build-tests（CLI + 27 CTest）；host 执行真实 MPI/PRTE/HYPRE 测试，未把 sandbox launch failure 记作数值通过。

| 项目 | 结果/证据 |
|---|---|
| 完整 `cmake --build build --clean-first --parallel 4` | 298 个 Ninja build actions 成功，无 compile/link failure |
| `cmake --build build-tests --parallel 4` | 成功 |
| architecture checker / diff whitespace | 通过；2 条既有 allowlist，无新增依赖 debt |
| 完整 host CTest | 27/27 通过；日志 `/private/tmp/phase29-ctest-final.log` |
| WENO7/Rusanov/FE Sod t=0.2 | 454 步；全部 frozen dt、time、积分、extrema 及 VTS SHA 校验通过，未调整阈值 |
| legacy PISO、constant-density PISO | frozen 数值及 continuity/gauge/boundary 验证通过 |
| PISO Poiseuille | 完整既有 profile/pressure-drop regression 通过 |
| SIMPLE/PIMPLE fixed-time | 同一 clock、PIMPLE outer=1 等价 PISO、relaxed flux continuity、early exit、zero gauge、SIMPLE Poiseuille 验证通过 |
| pressure MPI | cavity/channel 的 2×1×1、1×2×1、2×2×1，以及 SIMPLE/PIMPLE 分布式配置通过；matrix/RHS 和 U/p 使用原容差 |
| source/legacy A/B | conservative gravity/MRF/wallHeat、PISO gravity/MRF、SIMPLE gravity、PIMPLE MRF、Eulerian 8 组与 Phase 28C 独立冻结结果的 VTS/time-dt 完全一致 |
| Ghost Debug 两步 | dt=4.392680e-4、4.392260e-4；最终 min(rho)=1.22285、min(p)=101109，与 Phase 28D 记录一致 |
| Ghost 完整原配置 5 秒 | Passed：11450 步到 t=5；51 个 VTS（含初始和最终状态）与仓库已有 full Ghost baseline 逐字节一致；PVD 覆盖 0→5，进程 exit=0。Release 构建，不替代上面的 Debug 同配置短基线。长程 binary 使用已删除旧 AssemblyPlan gate 的 kernels；随后补充的 capability guards 不改变有效案例的编译结果，最终 Debug 与长程 Release 的 explain 完全一致 |

Sod frozen VTS SHA-256：`8bf55a1c4973f3a2322aeb8785305bf0589cc9eb4cf89505203309f0ccba0768`。最终 time=0.2；mass=98.67807117165903，momentum=[23619.052634273423, 1.907691669240853e-08, 0.0]，energy=21770158.12505413；rho/p 范围为 [0.25924657327767214, 0.7622215092077242] / [13037.114603387012, 68970.12209014251]。以上全部由回归脚本对最终实际输出验证，dt 序列逐项比较而非只比较末值。

首轮完整 CTest 有两处架构测试失败：DirectEvaluation 未声明 output workspace；pressure turbulence 在 Unsupported explain 前触发 solve-block ownership 异常。已分别移入方法提供声明、显式无 HOW capability reporting；随后完整 CTest 通过。没有为恢复绿色重建旧 term 或 pressure schedule authority。

Ghost full 最终 VTS SHA-256 为 `7ce841390cb215ea3b22fdddf910a0879e5eb5335e4f8a27a7948cb9dcec2c5f`；rho/p 范围为 [0.0, 1.2289806676854897] / [0.0, 101354.66376544903]，last dt=0.0002822697。以上 raw VTK 范围包括 138 个被 blank 的 IBM 非流体点，它们的 rho/p 导出值为零，是原有输出行为；IBMCellType=0 流体点的 rho/p 范围为 [1.2234377171537614, 1.2289806676854897] / [101290.53594640028, 101354.66376544903]。完整日志和 JSON 证据已存入最终 baseline 目录，输出保留于 `/private/tmp/sonic-phase29-ghost-full/result`。caffeinate 仅覆盖长程进程，结束后随之退出。

### 原配置 4-rank TENO5/Steger–Warming Sod 的已验证失败

新版本实际运行完整原输入，在 166 个已提交物理步之后、t=0.07459103 的下一次装配触发 `splitStegerWarming` negative-pressure fail-fast：rho=0.0846138、ru=-28.0405、E=1803.93、p=-1136.92。

用保存的独立 Phase 28C binary `/private/tmp/sonic-phase28c-baseline-20260928/sonicSolver` 运行完全相同输入：失败位置和两个 rank 的全部报错状态相同；166 条 time/dt 序列相同；initial/50/100/150 步共 16 个 VTS 逐字节一致。日志为 `/private/tmp/phase29-sod-mpi4{,-reference}.log`，输出保存在对应 `sonic-phase29-sod-mpi4{,-reference}/result`。

因此该完整 TENO/Steger–Warming case **Failed — independently reproduced pre-existing numerical instability**；不是 Unsupported，也不是 MPI launch failure，更不能写作 500 步通过。本轮没有观察到架构迁移首差异；保持 fail-fast，不修改格式、CFL 或阈值。稳定的 frozen WENO7/Rusanov 完整数值回归另列为 Passed。

## 8. 实现边界

| 分类 | 范围 |
|---|---|
| Implemented | authored single-fluid Formulas、structured constant-density HOW、production mode-free refs、Formula spatial bindings、shared precedence、stable temporal binding、method workspace validation、truthful runtime reporting、architecture guards/explain |
| Implemented backend | PressureOperators 的既有 arithmetic、ConservativeRHS 的 fused numerical orchestration、Time::Explicit 单 stage kernels、现有 Rhie–Chow/HYPRE/MPI helpers |
| Interface-only / proof | DirectEvaluation/LinearEquation 在通用单 patch Cartesian tests 中可运行，但没有宣称 native 任意公式和任意 pressure Formula 已 generic production lowering |
| Legacy | Eulerian Equation::Definition/专用 transport、turbulence specialized transport、部分 IBM equation/plan、无 HOW 系统的 role/target selection；后者隔离在 named legacy helper，迁移 Flow 禁止走该 fallback |
| Unsupported | 任意新数学结构的 fused lowering、generic pressure matrix 取代专用 callbacks、native 任意 formula/replace/disable、通用隐式 high-order flux、RK2/implicit stage 和未实现的 stage/constraint 组合 |

Ghost 验证只证明边界闭合及共享流动链；本轮没有改 J/JT、lambda、DLM/KKT 或 rigid-body algebra，不能把 Ghost 通过写成所有 IBM 方法已验证。分布式 pressure CTest 也不等于 distributed monolithic KKT 完成。

## 9. Self review 与后续范围

| 检查 | 结论 |
|---|---|
| migrated WHAT 是否直接 authored Formula | Yes |
| constant-density pressure macro HOW 是否唯一来自 root | Yes |
| method references 是否无 FormulaMode | Yes |
| single-fluid spatial binding 是否不读旧 Definition/Term | Yes |
| fused backend 是否删除旧 AssemblyPlan 数学开关 | Yes |
| migrated target 是否来自 ProgramStep Output | Yes |
| numerical compiler 是否自行理解微积分 | No |
| 第二 physical state/clock/registry 或全 workspace copy | No |
| numerical formula/stage/boundary/halo/IBM/COPY/SUM 是否改动 | No |
| hidden fallback 或 tolerance 调整 | No |

后续独立任务才处理 native custom Formula IO、通用存储/数值 lowering、legacy Eulerian/turbulence/IBM 剩余数学入口，以及 high-order 数值稳定性调查。本轮不重新设计 IBM/KKT，不新增 solver family，也不自动进入下一阶段。
