# SonicSolver：WHAT / HOW / STATE / WHICH / Compiler

SonicSolver 组合可解释、可执行的数值系统。主未知量、数学方程、模块贡献、耦合、时间方法、空间离散和运行时并行语义保持独立；case 不选择一个预打包的 solver family。

> WHAT defines equations.
>
> STATE defines base variables and storage/closure contracts.
>
> HOW invokes equations.
>
> WHICH selects numerical implementations.
>
> Numerical selection occurs once during compilation.
>
> Runtime provider lookup is not numerical selection.
>
> A compiled operation already owns its selected provider/backend.
>
> State realization may validate provider requirements, but it must not route solver families.
>
> Compiler binds and lowers; it does not make solver decisions.
>
> The compiler combines the four authorities without making their choices.
>
> The same equation may appear multiple times in HOW without duplicating its mathematical definition.
>
> Solver validates whether a configuration is executable, not whether the user's numerical choice is wise.

```text
MODULES
   ├── WHAT  ── EquationRegistry ─────┐
   ├── STATE ── StateRegistry ───────┤
   ├── HOW   ── ExecutionProgram ────┼── Compiler
   └── WHICH ── NumericalSelection ──┘       │
                                        Storage Binding
                                             │
                                      CompiledSolvePlan
                                             │
                                         PlanExecutor
                                             │
                                           Kernel
```

## Authority ownership

| Concept | WHAT | STATE | HOW | WHICH | Compiler |
| --- | --- | --- | --- | --- | --- |
| equation exists / mathematical AST | owns | | | | reads |
| base symbol / shape / physical storage / closure | | owns | | | resolves |
| semantic target / occurrence | | | owns | | resolves |
| order / generic Loop / physical Commit position | | | owns | | lowers |
| Euler/RK, WENO/flux, Central, HYPRE selection | | | | owns | validates |
| PISO/SIMPLE/PIMPLE execution topology | | | owns | | reads |
| Working/Correction view demand | | | requests | | binds |
| OldTime/Stage and backend workspace demand | | | | requests | binds |
| runtime operation provider/backend identity | | | | declares final implementation | freezes |
| provider capability validation | | | | declares capability | checks |
| pressure arithmetic implementation | | | | numerical provider | lowers |

禁止依赖：HOW transformer 写 NumericalBinding；Equation 保存 target/order；NumericalCompiler 选择 PISO；StateRealization 选择 topology；Builder 按 state flags 选择 solver family；source ExecutionProgram 持有 compiled workspace；WHAT transformer 创建 runtime OpIds。

## WHAT — Equation

WHAT defines mathematical equations. `SF::System::Equation` 保存稳定的 `id`、左右两侧 `FormulaExpr` AST 和 provenance。`EquationRegistry` 是数学定义 authority；它检查重复注册、未注册 lookup、缺失 replacement 和非法 AST/operator occurrence。Equation 不保存 execution target、order、循环归属或 numerical method。

实际 single-fluid density registry 是 `continuity`、`momentum`、`energy`。常密度 NS registry 是 `momentum`、`continuity`；后者是 `div(U)=0`，由压力一致性方程满足，不额外伪造一次独立 continuity solve。压力耦合注册 `pSimple`、`correctU`、`correctP`、`correctFluxp`；固定时间迭代另有 `relaxIterate`、`restoreFlux`、`checkConvergence` 数学关系。

Builtin 与 C++ user-defined equations 经 `SystemCompositionBuilder::addEquation(Equation)` 进入同一 registry。编译器/provider 根据 AST、target 和已实现的 capability 绑定，不能按 builtin/user 来源切换执行路径。native YAML 的任意数学表达式解析仍未实现，不能把 C++ 注册接口误说成任意 YAML 方程已可执行。

Formula/Expression AST 是 equation mathematical implementation detail，不是顶层 execution authority。`SF_formula.h/.cpp` 实现 AST、Equation 与 registry；保留此文件名不意味着恢复 Formula solver。FormulaGroup/FormulaGroupRegistry 已移除，守恒向量 fusion 属于 numerical provider 优化。

未迁移模块仍有 `SF::Equation::Definition` in `legacyDefinitions`、`EquationDescriptor` in `legacyEquations`、`LegacyEquationRole` 和 `TermKind` 平面 DSL。它们用于 legacy numerical backend 及诊断，不能重新控制已经迁移的 single-fluid 方程或执行顺序。density、常密度 NS 和 native Eulerian core 不再创建平面数学副本。AssemblyPlan 已没有数值消费者，相关类及实现已删除。

## STATE — Base variables and realized views

`core/system/SF_stateRegistry.h` 的 `StateRegistry` 与 `StateSymbol` 是变量 namespace authority。符号保存 shape/components、location、parallel ownership、primary/derived role、physical storage contract、closure derivation、initialization/BC/restart 属性。没有 equation id、target binding、order、Loop、PISO、RK recipe。

原生 case 必须分别选择 WHAT、solution STATE 和 HOW。`EquationCompositionConfig::stateDeclared/solutionVariables` 是输入；冻结后 `StateRegistry::solutionVariables()` 是成员身份 authority。`StateRole` 仍描述数学/存储职责，不决定是否属于 solution；不增加 `Conserved` 互斥角色。物理守恒密度及其 balance 属于 WHAT，一个主未知量可以同时表示守恒密度。

```yaml
# state/state.yaml
SonicFile:
  object: state
  type: registry
use: [U, p]                 # 或 [rho, rhoU, rhoE]
```

`BuiltinStateCatalog::solution(id)` 提供已选择变量的 writable metadata：U 为 packed velocity、p 为 named algebraic pressure、rho/rhoU/rhoE 为 packed conservative components。`require(id)` 只激活 dependency；保守 STATE 下的 U/p/T 保留 EOS/closure view。p 是否参与 multiplier constraint 由已声明数学关系决定，不能永久在 catalog 中设为 Multiplier。

EOS 只贡献 closure。rhoConst 在未选择 writable rho 时激活常数 rho；既不添加 solution U/p，也不选择 SIMPLE/PISO/PIMPLE。`validateSolutionProviderContract` 在不支持的数学实现进入 lowering 前报告缺失 provider：variable-density primitive evolution、rhoConst 下 writable rho evolution、primitive energy executor 等；不改选 STATE/HOW/WHICH。

`installSingleFluid` 不注册 solution pack，WHAT 物理 identity 保持 continuity/momentum/energy，默认 HOW target 消费显式 STATE。`addConstantDensityFluid` 只在显式选择 U 后，用 rho=rho0 的数学条件归约 continuity 与 momentum；不会把 variable-density ddt(rhoU) 直接换成 ddt(U)。U/rhoU 同时被选为默认 Momentum target 时 fail fast。

旧 densityBase/pressureBase 仅由 `CaseAdapter::build` 输出标明 compatibility provenance 的 WHAT/STATE/HOW tuple。`inspectCase` 和 builder 不读取旧 family 标签。原生缺 STATE 或 HOW 即报错；完整原生输入的 HOW 优先于旧控制文件里的选择标签，线性/迭代容差保持各自嵌套作用域。

Raw/ExecutableEquationSystem 是 composition snapshot，包含独立的 `registry` (WHAT) 和 `state` (STATE) channel。这不是两份运行时 physical state。`compileExecutionProgram` 主接口显式接收 WHAT、STATE、HOW、WHICH；保留 snapshot convenience overload。`Program.state` 引用同一 registry，runtime `realizeState` 只接收 STATE 与编译需求，不读取方程排序。

Solver 的 `SF_builtinState.h/.cpp` 中 `BuiltinStateCatalog` 只保存 known metadata：rho、rhoU、rhoE、U、p、T、h、phi、alpha、k、omega、mu、nu、conductivity。它不是 active STATE，没有 StateBundle/Field、物理数组或模型生命周期。`StateRegistry` 只包含当前 composed system 请求的 base symbols。

模块用 `SystemContribution::requireState(id)` 或 `SystemCompositionBuilder::requireState(id)` 请求 builtin metadata；HOW 的非 Workspace target 通过 `requireTargetStates` 请求 base symbol。Composition 在 freeze 前按需激活；显式 `addState(customState)` 的 storage/closure contract 优先。未知 symbol 必须给出自定义 metadata。激活变量不会增添 Equation、HOW、WHICH 或模型。

例如本阶段 Sod 的 active STATE 为 rho/rhoU/rhoE/U/p/T；圆柱 Ghost 另请求 mu/conductivity。常密度压力 case 为 U/p/phi，以及实际使用的 rho 常数 closure、可选 nu coefficient。闲置 k/omega/alpha/rhoE 不因目录可用而进入这些压力 case；k/omega 由已启用 turbulence module 的明确声明激活。用户仍仅为新自定义变量注册 metadata，为预制模型提供初始化、边界和参数。

`CompiledExecutionProgram.stateViews` 是 compiler 输出的 lazy manifest：HOW qualifiers 要求 Working/Correction；WHICH temporal lowering 要求 OldTime/Stage。未引用的版本没有数组。相同 base/kind 只有一份 authority；不同 provider 要求不同 backing 时明确冲突失败。

`StateRealization` 绑定 StateBundle 的 physical authority。Generic DirectEvaluation/LinearEquation 的临时数组归 solver execution lifetime；pressure backend 绑定已有 predictor、pressure correction 和 canonical face arrays。RK OldTime 引用 q0，Stage 是当前 published physical state 的时间限定 view，仅在对应 stage 内可读，不创建第二套 Q。Thermodynamic U/p/T/h/transport values 使用 readonly closure views，写入 derived physical target 明确拒绝。

`StateBundle.stateModel` 是 thermodynamic/layout closure 引用，不是 equation-to-target binding；`validateThermodynamicBinding` 验证它与各 patch 一致。clock 仍只存在于 StateBundle。

详细实现与边界见 [four-module-state-migration.md](../docs/four-module-state-migration.md)。

## HOW — Execution

公共控制词汇只有 `EquationCall`、`Loop`、`StageLoop`、`Commit` 和 `Sequence` 容器，定义在 `core/system/SF_solveProgram.h`。

```cpp
EquationCall { EquationRef equation; Target target; occurrence; }
ExecutionScope { kind; children; Order order; repetitions;
                 minimumIterations; terminationSignal; origin; }
```

HOW 回答求哪个方程、结果写到哪里、先后顺序和迭代结构。`ExecutionProgram.root` 是唯一 source execution authority；没有第二份 flat steps。一个 equation 可以重复引用：

```text
EquationRegistry: momentum -> one mathematical definition

HOW:
10 momentum -> U*
50 momentum -> U
```

这是 one mathematical equation / two execution occurrences。数值 provider 若尚不支持特定 target/stage 组合，明确拒绝该 capability，不能复制 momentum definition 来掩盖缺口。

### Scope-local order

`order is local to an execution scope`。`orderExecution` 对每个 scope 的 children 使用 `std::stable_sort`，再递归处理子 scope。相同 order 按稳定的模块注册/插入顺序打破平局。没有全局 flatten、依赖 DAG 或猜测 CFD 方程先后。排序后的 scope path 是默认 occurrence address；作者也可提供显式 occurrence ID。

### Typed targets

Source `Target {symbol, kind}` 仅表达 `Physical / Working / Correction / Workspace`。它不持有 workspace 名字、resource、offset、pointer 或 MPI binding。`targetFromSyntax` 在 composition 时把 `U*` / `p'` 转为 symbol 与 kind。

`CompiledTarget` 由 provider 提出语义 target 与可选 backend storage demand，再由 compiler 的 STATE binding 填入 storage、components、offset、access 与 synchronization bindings。所有 Physical/Working/Correction target 都查询基础 StateRegistry；provider 不能因 target 是 working 而绕过基础符号验证。Plan leaf / ExecutionContext 只消费冻结后的 CompiledTarget。

| 表达 | Storage semantics |
|---|---|
| `U` | StateBundle 的 physical/current target |
| `U*` | execution-owned predictor/working storage |
| `p'` | Correction(p)；native STATE 不注册 pPrime，数学 AST 可以使用修正量记号 |
| `phi` | Physical(phi)，canonical face authority，由 pressure backend 提供 storage |

`* does not prescribe a particular numerical formula`，不规定只能读 old-time。lagging、relaxation 和 stage semantics 由数值 provider 决定。runtime 不解析 `*` 或 `'` 字符串；Plan leaf 携带 typed Target，generic ExecutionContext 将其传给 numerical callback。

常密度 `PressureOperators` 的 `workingVelocity_`/view 保存 predictor；`correctU -> U` 把边界和 halo 已就绪的结果 COPY 到 physical velocity。它不注册第二个 U field、clock 或 StateBundle。`-> U` 的写入不等于 physical-time commit。

### Generic loops and commit

Loop 是 nonlinear/correction iteration；StageLoop 是 temporal topology。通用 executor 只读取 repetitions、minimumIterations、terminationSignal 和 children，不知道 SIMPLE/PISO/PIMPLE、pressure、alpha 或 IBM。convergence signal 在每次完整 loop body 后检查。

不可变 time recipe 的 coefficients/stage times 仍由编译后的 time recipe 提供。HOW 不创建 RK4Loop/EulerLoop，也不单独定义 stage 数。

Source scope 没有 before/after，也没有 raw OpIds。pressure preparation 和 fixed-time snapshots 由 PressureMomentum provider 的 root lifecycle lowering 生成；outer iteration begin/end 仅附着于包含对应 predictor 的 termination Loop；correction publication 位于 FluxCorrection provider 的 compiled fragment；Commit 的 state/clock operations 由 provider 实现。density 与常密度 pressure 的 source HOW 都显式包含最终 Commit。

这些均为 compiled lifecycle，不重建 coupling topology，不移动 sibling calls。数学 relaxation/flux consistency/convergence 仍是注册关系和 EquationCall。

## WHICH — Numerics

`NumericalSelection`（`solver/system/SF_numericalSelection.h`）是独立 source WHICH，保存 numerical bindings、recipes.time、空间 recipes、线性配置、provider numerical inputs 和 operator overrides。`selectNumerics` 在 composition 阶段完成 builtin 默认值，compiler 不补猜 diffusion 或 pressure numerical method。

WHICH 是 temporal method、spatial method、flux method、linear solver 与 equation provider 的选择。它不增加 equation、不改变 order、不创建 correction Loop、不选择 solver family。

`NumericalBinding` 绑定 provider 到 equation default 或具体 occurrence；优先级为 global `*` default < equation default < execution occurrence override。最高优先级重复、悬空 binding、缺失 provider 均 fail fast；低优先级 defaults 不影响更高优先级的选择，注册顺序不是 fallback。operator binding 继续使用 global operator default < equation default < mathematical operator occurrence；WENO/TENO、Rusanov/Steger-Warming、Central 与 source kernel 的数学保持原样。

`IProvider::runtimeProvider/operationProvider` 是 selected method 的 implementation ownership contract。ConservativeResidual 固定绑定 flow.conservative；常密度压力 methods 固定绑定 flow.pressure-operators，ConservativePressure* methods 固定绑定 flow.conservative；FluxCorrection 的 FluxCorrect leaf 绑定 flow.rhie-chow，而 correction publication 仍绑定 flow.pressure-operators。`ProviderCatalog::resolve(operation, selectedProvider)` 只验证这个 ID 的实现能力，不搜索候选、不读取 solver-family flags。

`IProvider::compile` 只生成某个 EquationCall 的数值 realization，例如 assemble/solve 或 evaluate。`ProviderRegistry` 中的 PressureMomentum/PressureCorrection/VelocityCorrection 等是 numerical implementations，不是 HOW node、另一套数学 Equation 或 solver lifecycle。

当前守恒 provider 将 ordered `continuity -> rho`、`momentum -> rhoU`、`energy -> rhoE` calls 编译融合为一个 packed conservative stage kernel。WHAT 仍有三个方程，HOW 仍有三个 calls 与显式 Commit。fusion 不得吞掉其他 calls 或 loops；没有实现的混合/迭代 stage topology 明确报 capability missing。

现有时间 recipe 为 forwardEuler、SSPRK3、classicalRK4。常密度压力 backend 当前只覆盖单 stage、现有 spatial/constraint capability。PISO + RK4 等组合若不可执行，原因是缺少相应 stage/predictor implementation，不是用户选择“不聪明”。不修改 CFL、RK coefficients、linear tolerances 或 flux 数学来绕过 capability 缺口。

## Pressure coupling

Pressure module 有三个独立通道：PressureConstraintTransformer 注册数学 relations（WHAT）；`applyPressureExecution(program, request)` 只变换 HOW；`pressureNumerics(request)` 只贡献 WHICH。HOW 函数不接收 NumericalBinding，不知道 PressureMomentum/PressureCorrection provider IDs。

Base equations：`momentum`、`continuity`。SIMPLE/PISO/PIMPLE module 请求所需 pressure relation registration，并用 `applyPressureExecution` 消费 base momentum occurrence。它们不提供第二个 momentum definition，也不重建其他模块的默认 HOW。

实际 PISO：

```text
10 momentum -> U*
20 Loop pressure [pressureCorrectors]
   Sequence pressureCorrection
      10 Loop nonOrthogonal [nonOrthogonalCorrectors + 1]
         10 pSimple -> p'
      20 correctU -> U
      30 correctP -> p
      40 correctFluxp -> phi
1000000 Commit
```

实际 SIMPLE：

```text
20 Loop outer [outerCorrectors, min=1, convergence]
   10 momentum -> U*
   20 Sequence pressureCorrection [same relations and nonOrthogonal scope]
   30 relaxIterate -> iterate*
   40 restoreFlux -> phi
   50 checkConvergence -> converged
1000000 Commit
```

实际 PIMPLE：

```text
20 Loop outer [outerCorrectors, min=1, convergence]
   10 momentum -> U*
   20 Loop pressure [pressureCorrectors]
      Sequence pressureCorrection [same relations and nonOrthogonal scope]
   30 relaxIterate -> iterate*
   40 restoreFlux -> phi
   50 checkConvergence -> converged
1000000 Commit
```

Predictor 保留原 NS AST，由 PressureMomentum provider 调用已有 PressureOperators assemble/solve。pSimple 保留原 pressure matrix / continuity-defect 数学。velocity correction、pressure update、canonical face flux correction、relaxation、flux restore 与 convergence arithmetic 不重写。

transformation provenance 显示 `default momentum -> U/rhoU @ 10` 被哪个 preset 变为 working predictor；其他模块 entries 保留其 origin/order。所有 single-fluid pressure 不再注册独立 PressureSchedule LegacyExecutionPolicy；loop/count authority 只在 HOW。Eulerian shared pressure 也由 native HOW 保存三层循环次数；旧 sharedPressurePlanFragment/policy 已删除。

## Single-fluid conservative / variable-density pressure

常密度与保守变量压力共享 WHAT/STATE/HOW/WHICH → Compiler → owned CompiledSolvePlan → Runtime。架构统一不要求数值 kernel 统一：前者使用 flow.pressure-operators / flow.rhie-chow，后者使用 flow.conservative 与原 ConservativeRHS、Time::Explicit::forwardEuler、PressureBased::Corrector。

原生保守压力示例显式选择 WHAT `[Continuity, Momentum, Energy, PressureConstraint]`、STATE `[rho, rhoU, rhoE]` 和 HOW `PISO`。PressureConstraint 是已有数学约束的库入口，不是 solver family。WHAT 复用同一 continuity/momentum/energy identity（ddt(rho/rhoU/rhoE)），新注册现有校正器的 pressure compliance、inverse-density diffusion、动量修正、prepared pressure、derived-flux invalidation 和 EOS energy publication 关系。conservativePressureRelations() 同时服务注册、explain 与 selected method 的 AST capability 校验；修改关系而没有相应 kernel 就 fail fast。

STATE 为 packed rho/rhoU/rhoE 与 EOS-derived U/p/T。原预测器已直接发布 Physical(Q)，因此不伪造独立 Working(Q)。Correction(p) 和 Working(p) 分别非 owning 地绑定原 Corrector::correction/targetPressure；Workspace(fluxValidity) 绑定已有 fluxCorrected flag。此算法不存独立 canonical phi；下一次 spatial assembly 从已修正状态重建 derived flux。常密度 Physical(phi) 的 canonical face COPY 语义保持原样。

实际 runnable 保守压力 HOW：

```text
10 momentum@predictor -> Physical(rhoU)
   WHICH: ConservativePressureMomentum
   explicit fusion inputs: continuity -> rho, momentum -> rhoU, energy -> rhoE
20 Loop pressure [pressureCorrectors]
   Sequence pressureCorrection
      10 Loop nonOrthogonal [nonOrthogonalCorrectors + 1]
         10 pSimple -> Correction(p)
      20 correctU -> Physical(rhoU)
      30 correctP -> Working(p)
      40 correctFluxp -> Workspace(fluxValidity)
      50 publishPressure -> Physical(rhoE)
1000000 Commit
```

PressureConstraintTransformer 读取 composition 解析出的默认 momentum occurrence target（U 或 rhoU），不再用 rho.constantValue 选择关系。closure 只参与当前实现能力的校验。

HOW 是唯一 loop/count/order authority。WHICH 的局部 fragment 提供原 assemble/solve/correction kernels；root/Commit 的 pressure.prepare、pressure.step.commit、time.commit 按 selected method lifecycle lower。所有生产 leaves `legacyAdapter=false` 且 owner=flow.conservative。bindPressureOps 保留为 callback registration；无主时间或 corrector loop。

现有 conservative provider 仅支持单 patch、单 stage、PISO、无额外非正交 pass 的 laminar single-fluid contract。conservative SIMPLE/PIMPLE fixed-time、RK stage treatment、多 patch/MPI conservative pressure、IBM/turbulence/phase services 尚未实现；编译或启动 capability validation 明确拒绝。保留 generic outer/pressure/nonOrthogonal HOW 不代表这些 kernel 已支持。

## Module composition

`SystemContribution` 统一贡献 equations/AST extensions、unknowns/closures、execution entries 与 numerical defaults。source contribution 不拥有全局 timestep 或 MPI barrier。

| Module | Equation contribution | Execution contribution | Numerical contribution |
|---|---|---|---|
| Navier–Stokes | momentum, continuity；density 另含 energy | ordered default calls | conservative / pressure momentum providers |
| MRF、gravity、wall heat | extends momentum/energy AST | none | existing compiled source kernels |
| Single-fluid SST / k-epsilon | separate native k, omega/epsilon AST；provider-owned mu_t closure | order 1/2 generic EquationCalls before flow calls | TurbulenceTransport；explicit pair fusion -> flow.turbulence |
| Eulerian turbulence | native phase-mass weighted k/epsilon/omega AST、implicit sink 与 mu_t algebraic closure；原数组 ProviderDistributed aliases | 显式 outer closure group → flow groups → transport group，完整 phase/model fusion | EulerianTurbulenceClosure / EulerianTurbulenceTransport → flow.eulerian-turbulence，原矩阵内核 |
| VOF / level set | alpha/phi transport and closures, current module definitions | contribution interface available | Legacy numerical backend / migration pending |
| PISO | pSimple, correctU, correctP, correctFluxp；conservative 另有 publishPressure | consumes momentum call; inserts pressure Loop | constant / conservative relation providers |
| SIMPLE | pressure and fixed-time relations | outer Loop with single correction | constant density: original PressureOperators；conservative fixed-time Unsupported |
| PIMPLE | pressure and fixed-time relations | nested outer/pressure Loop | constant density: original PressureOperators；conservative fixed-time Unsupported |
| Eulerian–Eulerian | namespaced native Formula equations | native PIMPLE EquationCalls / Loops / Commit | AST-validated compiled contracts / specialized kernels |
| Ghost IBM | native boundary read/write/stage/order contract；不增加虚假 equation | every spatial evaluation 的原 BC/halo/ILW 链 | explicit ibm.boundary port；原 kernel |
| forcing / DLM / KKT | native impulse/work、lagged force、projection 或 coupled constraint AST；original STATE aliases | explicit complete post-predictor group / block | immutable native Immersed method contract；原 kernel；KKT production capability 仍有缺口 |

`legacyEntries` 明确记录尚未独立执行的模块 calls，explain 显示 migration pending。不能同时把这些 entries 和已存在的 legacy local backend 执行一次，造成湍流/相方程重复推进。Eulerian core 的 shared-pressure fragment 和 AssemblyPlan 已删除；Eulerian turbulence 已使用原生 closure/transport occurrences，不再进入 legacyEntries。

## Compiler and execution

```text
Collect WHAT/STATE/HOW/WHICH       SystemCompositionBuilder / module contributions
Select Solution + Dependencies    BuiltinStateCatalog::solution + require / requireTargetStates
Freeze Source Composition         EquationRegistry / StateRegistry / ExecutionProgram / NumericalSelection
Order Each HOW Scope              compileExecutionProgram -> orderExecution (independent copy)
Resolve WHAT + STATE Targets       EquationRegistry::at / StateRegistry::at
Lookup Selected WHICH Method      ProviderRegistry::at
Bind Method Implementation        IProvider::runtimeProvider / operationProvider
Realize Required Views            realizeTarget / realizeTemporalViews
Bind Operator Numerics            NumericalCompiler::compile
Lower Owned Executable Operations SolvePlanner::compile / compileMethodProgram
Validate + Freeze Final Owners    compileOperationBindings -> selected-ID ProviderCatalog validation
Execute Frozen Decisions          PlanExecutor / owned OpRegistry / existing kernels
```

`SystemBuilder::build` 收集 modules 的独立 contributions。`composeContributions` 做 module applicability matching 和显式 legacy adapters；native pressure 匹配已贡献的 mathematical/default occurrence，不读取 StateRealization flags。composition 完成后保存 source NumericalSelection，再调用 `NumericalCompiler::compileSystem`。

`compileExecutionProgram` 后校验每个显式 solution variable 是否由 method 的 `writes` 或 writable resource 覆盖；保守 predictor 的 fusion 同时写 rho/rhoU/rhoE，不能只看外层 rhoU target。未使用的 solution T 等明确失败。

Compiler 通过 const source HOW/WHICH 和 const equation/provider interfaces 接收输入。scope sorting 发生在 compiler 的独立 HOW 副本，source 不被修改；provider 不得改变 equation/occurrence/target semantics。编译器不增添 equation、不选择 preset、不替换 numerical selection。compiled operation/resource metadata 是输出，不是另一份数学 authority。

`SF_methodObjects.cpp` 只做 deterministic occurrence binding、capability validation、generic lowering；domain numerical implementations 位于 `SF_builtinProviders.cpp`。`SF_numericalCompiler.cpp` 按 provider 声明的 spatial terms 和 kernel/workspace requirements 绑定，不按 momentum/U/constantDensity 识别 solver。旧 flat DSL 选择由 `Legacy::spatialInputs` 在 composition 显式提供，generic compiler 不导入 legacy numerical adapter。不新建平行 SimulationPlan、manager 或 forwarding planner。

`SolvePlanNode.provider` 是每个 executable leaf 的最终 implementation identity。`RuntimeRequirements.operationBindings` 只是该 frozen leaf ownership 的校验/index 投影；`SystemValidator` 拒绝两者不一致。`OpRegistry::bind(id, provider, callback)` 注册实现的 owner，`PlanExecutor` 在任何 operation 执行之前检查所有 ID 与 owner，然后直接 invoke；它不读取 state/formulation 或 coupling preset。SingleFluidStepper 按 frozen provider ID 构造/绑定资源，不再由 constantDensity flag 选择 backend。

只有显式 `legacyAdapter` plan leaf 允许在编译时调用 `Legacy::selectOperationProvider`。IBM constraint 已由 native Immersed method 冻结 owner/operation，不进入兼容 selector；native Eulerian shared pressure 同样不进入该 selector。Eulerian turbulence 使用独立原生 numerical provider，不进入 compatibility selector。所有 native single-fluid pressure 不走该 adapter。不存在 compile 后的 provider matching。

`CompiledSolvePlan.root` 是冻结后的 executable operations，不再是独立 authored HOW；source program、compiled calls、numerical bindings 与 explain 都来自同一编译链。compiler 不猜 equation target/order，也不选择 coupling preset。

`sonicSolver explain` 首先分别显示 WHAT / EQUATIONS、STATE / BASE VARIABLES、STATE VIEWS、HOW / EXECUTION、WHICH / NUMERICS 的 source choices；之后显示 COMPILED OCCURRENCES / STORAGE / PROVIDERS 和 COMPILED SOLVE PLAN，明确区分 source 和 lowered outputs。Legacy execution authority 与 pending entries 明确展示。unsupported 的真实边界必须保留，不能用解释文字冒充运行实现。

## Native single-fluid turbulence and mixed temporal lowering

显式 primary flow STATE selection 不要求用户把每个模型辅助变量加入 `state/state.yaml`。用户显式选择 kOmegaSST 后，模型在同一个 StateRegistry 贡献 k/omega；选择 kEpsilon 则贡献 k/epsilon。用户 primary flow selection 仍是 rho/rhoU/rhoE，U/p/T 是 flow closure views，mu_t 是模型 derived closure output。没有第二个 STATE registry。

RAS 的 k/omega/epsilon/mu_t backing 仍由 Turbulence::Manager 内的 ScalarFields 持有；`ProviderDistributed` 描述这一事实。Manager::registerDistributed 发布原数组的非 owning views，StateRealization 和 Physical(target) 直接引用它们，compiled view owner 为 NumericalProvider。没有为架构另分配数组；ScalarFields 原先分配的四个槽位仍保留。当前没有 turbulence restart reader，因此这些 STATE 不宣称 restart capability。

single-fluid native WHAT 如实描述现有 production RAS 内核：ddt(k/omega/epsilon)，显式 production/dissipation、除以 rho 的固定局部扩散，以及 SST cross diffusion。旧注册形式包含实际未执行的 div；本阶段删除这个数学描述偏差，没有给内核增加对流或隐式汇。kEpsilon 扩散仅 mu_t/sigma；SST 为 laminarMu + sigma*mu_t。既有 positivity floors、production limits、wall-distance 逻辑和 in-place sweep 完整保留，不代表数值缺陷已解决。Eulerian 的质量加权输运和隐式汇是另一个旧内核，未与此实现合并。

source HOW 仍是通用 EquationCalls：order 1 k->k、order 2 omega/epsilon->自身，随后 flow 的 continuity(10)、momentum(20)、energy(60)，最后 terminal Commit。WHICH 为两个 RAS occurrences 显式选择 TurbulenceTransport，并声明完整 pair 输入。编译出的 fusionKey/fusionMembers 合约要求同一 Sequence 中两个相邻 sibling、正确目标、正确 model/state/math、无缺失或重复。两份 WHAT 和两份 source HOW 保留，只有数值 backend operation 合成一次 turbulence.advance，owner=flow.turbulence；共享 provider ID 本身不会触发这种 fusion。

通用 temporal compiler 支持最小 mixed pattern：flat independent once-step prefix -> one contiguous fused temporal group -> terminal Commit。ITemporalMethod 接收编译后的 physicalStepPrefix；explicit recipe 在原 prepare/dt 后、begin/snapshot 前安放它。因此实际顺序为 flow.step.prepare -> flow.dt.compute -> turbulence.advance -> flow.step.begin -> StageLoop(flow RK) -> flow.step.commit -> time.commit。temporal -> independent -> temporal、嵌套 mixed scopes、缺失 terminal Commit，以及尚未指定 temporal placement 的额外 root lifecycle decorations 均明确 Unsupported，不猜测、重排或丢弃。此 prefix 是降阶输出，没有新增 source scheduler 或第五模块。

TurbulenceTransport 保留当前一次/物理步的显式原位更新，不采用 flow RK4 stages。它消费 StateBundle::dt；不计算第二个 dt 或推进 clock。FlowStepBegin 已移除 correctTransportModel 调度；CompiledSolvePlan 中的 turbulence leaf 绑定旧 correction callback。其局部 BC -> WriteOwned -> ReadHalo -> correct -> BC -> WriteOwned -> ReadHalo 顺序保持不变，MPI 实现仍由 execution runtime/backend 承担。application 只构造 Manager 和参数，不再通过 transportedTurbulence predicate 选择执行拓扑。

DNS 不贡献输运方程；Smagorinsky 仅贡献 mu_t 代数闭合及 TurbulenceClosure refresh，保持原缓存刷新位置，不产生 k/omega/epsilon 输运。全局 term recipes 只由声明 usesSpatialRecipes 的 provider 消费；RAS 的固定局部 Central2 不会误触发 flow diffusion recipe 默认选择。

RAS 当前执行能力仅 single-fluid、serial、single patch、density explicit flow；pressure、MPI/multi-patch、IBM、Eulerian 和 implicit RAS 未放行。旧 single-fluid pressure service guard 保留；Eulerian 的 EeTurbulenceSolve 已是独立原生 provider callback；Eulerian core scheduler 与 assembly authority 均已原生化，AssemblyPlan 已删除。详细审计、诊断对照的基线来源及 18 项自审见 [single-fluid turbulence migration report](../docs/single-fluid-turbulence-migration.md)。

## State, workspace and parallel invariants

Physical state/clock/field composition 的唯一 authority 仍是 StateBundle。预测量、p'、flux、RK k 值和 linear/constraint buffers 属于 execution/provider lifetime，不塞回 Field，也不创建新 physical registry。

Runtime interfaces 隔离 MPI；raw MPI 仍仅在 infrastructure/backend。state、geometry、metrics、labels：owner -> COPY -> replicas。每个 physical face 只有 canonical F*，COPY 后 owner +F* / neighbour -F*。residual/source/load/constraint contributions：local -> SUM -> canonical owner。不得平均 state/geometry/face identity，或把 SUM 改为 COPY。

本轮 architecture migration 不改 WENO/TENO、flux split、pressure/Rhie–Chow arithmetic、RK coefficients、CFL、HYPRE tolerance、IBM mathematics 或 distributed ownership。编译和 explain/static checks 不代表数值基线验证。

本阶段先完整运行原始 Sod 与 cylinderFlowGhost，再在 provider authority/active STATE 修改后运行同一份输入。详细 frozen evidence、数值验证状态和未覆盖项见 [runtime provider authority migration report](../docs/runtime-provider-authority-migration.md)。原始 cylinderFlowGhost 是 density/RK4/Ghost，不能作为 pressure Working/Correction 或 MPI/turbulence CFD coverage 的证据。

## Native Eulerian shared-pressure execution

Eulerian PhaseSystem 明确贡献每相 phaseMass/momentum/enthalpy、shared p 和派生 alpha/rho/U/h/T，用户不必重复列为单流体 primary STATE。ProviderDistributed views alias 原 PhaseSystem/PhaseState 数组；Correction(p) alias 原 pressureCorrection，canonical face arrays、diagonal、previous state 保持 numerical workspace。

`SF_eulerianRelations` 描述冻结数学：参考相 mass 由 volume closure 重建；动量压力项为 alpha*grad(p)，扩散为 RHS；shared pressure 包含 sumGas(alpha/p)/dt reaction、sum faceAverage(alpha²/A) diffusion、volume source/divergence/history；energy 包含 material pressure work。旧 phantom/不准确描述已纠正，原算术未升级。

明确注册的 PIMPLE 贡献 generic Sequence/EquationCall/Loop/Commit：outer 内 continuity group -> momentum group -> pressure/nonOrthogonal 下 shared-pressure correction、phase correction、paired phase flux correction -> enthalpy group，最后 Commit。三层次数只在 source HOW；当前其他 Eulerian coupling 明确 Unsupported，不自动替换。source 无 raw Ee operations。

六个 Eulerian methods 显式选择 flow.eulerian-pressure，完整 phase inputs/group targets 按声明的 PhaseSystem storage order 验证。通用 compiler 支持 N-member fragment fusion，保留每相数学 provenance，core leaves legacyAdapter=false。Provider-local prepare/source/BC/halo/canonical-copy lowering 保留原 Ee 微操作顺序。`sharedPressurePlanFragment`、S_EE_PIMPLE policy 与 transformer Ee inventory 已删除；EulerianStepper 不从 policyKind 选择算法，只执行 frozen plan。

native Eulerian 数学只在 Formula EquationRegistry 中表示一次。所选 provider 直接校验准确 AST、源项扩展和 STATE/backing，产生不可变的 CompiledEulerianAssemblyContract：方程/输出、实现关系、phase slot/name、method/provider 与支持的扩展。generic compiler 只携带不透明 provider 编译产物，不解释 Eulerian 类型。

```text
WHAT Formula AST -> selected WHICH provider -> strict AST capability matcher
    -> immutable compiled provider contract -> HOW lowered CompiledSolvePlan
    -> PlanExecutor -> flow.eulerian-pressure -> existing numerical kernels
```

EulerianStepper 启动时只验证冻结相槽与实际 PhaseSystem；不再解析方程名或遍历 AST。PhaseEquationAssembler 是专用数值实现，不是数学 authority 或通用 AST interpreter；不再绑定 Equation、Definition、TermKind 或 AssemblyPlan。旧投影 helpers 和没有真实数值消费者的 AssemblyPlan infrastructure 均已删除。未知数学结构明确 Unsupported。

架构法则：native 数学方程由 Formula-based EquationRegistry 唯一表示。数值 provider 可以从它编译 typed implementation contracts，禁止再创建另一份可编辑数学方程。

Eulerian turbulence 的 identity/HOW/provider 断点已完成原生迁移；其真实数学和执行契约见下一节。core assembly 三阶段门槛、源码保护和数值证据见 [Eulerian assembly authority migration](../docs/eulerian-assembly-authority-migration.md)；调度迁移历史见 [Eulerian shared-pressure migration](../docs/eulerian-shared-pressure-migration.md)。

## Native Eulerian turbulence

模型贡献 `k.phase` 与 `epsilon.phase`/`omega.phase` 方程，以及真实代数闭合 `mu_t.phase`。WHAT 对应原矩阵：ddt(phaseMass*phi) + div(canonicalMassFlux*phi) + implicitSink*phi = diffusion(D,phi) + explicitSource。kEpsilon 与 SST 的相质量、alpha 加权、生产限幅、耗散汇和 SST cross diffusion 保留；closure coefficients 在预测器前计算、输运时滞后使用，求解后刷新。没有把 single-fluid 的显式原位 RAS 内核复用于此矩阵路径。

STATE 继续由 EquationSystem::PhaseEquationState 的原数组持有；model contribution 用完整 PhaseSystem 顺序解析真实 slot，selected phases 可以是 subset 或反序。registerState 与 synchronization 使用 phaseIndex，不用 selected-array ordinal。没有第二份 k/omega/epsilon/mu_t、Field 或 clock；暂不宣称 restart reader 已实现。

source HOW outer 明确包含 closure group → continuity → momentum → pressure/nonOrthogonal corrections → enthalpy → transport group。EulerianTurbulenceClosure/Transport 选择独立 owner flow.eulerian-turbulence、显式 N-phase/model group 与不可变 providerContract；每组融合一次原 prepare/solve。flow provider 不再通过字符串或 legacy identity 插入湍流。closure provider 的局部 prerequisites 实现 interphase → sources → prepare，continuity 不再重复前两项；outer final BC/validate 放在最后一个 authored group 后。source sibling 顺序不重排。

所选 provider 的支持签名只用于严格 AST/backing/phase/model 能力校验，不注册第二份 WHAT；system 不 include concrete model。runtime 只校验冻结实现绑定与实际模型数组，不重新解析方程名。generic compiler 不识别 turbulence，explain 展示每个 occurrence 的模型、slot、group、owner 及数值方法。

保留完整生命周期：initial prepare/commit；dt probe 中 prepare/sourceTimeStep；每 outer 在预测器前 prepare/同步、energy 后 solve、solve 内 BC/validate/reprepare；outer final BC/同步；成功物理 step commit(m*phi)，随后 diagnostics、clock commit。这些 prepare 不可去重。

当前能力为 serial single block、既有 Upwind/diffusion/HYPRE backend；MPI/multi-patch Eulerian turbulence 在 capability validation 明确 Unsupported。实测 2-rank host 输入因缺少分布式 phase transport provider fail-fast，是 implementation capability gap，不是 sandbox launch failure。已有 MPI pressure regressions 保留。

数值基线来源、独立矩阵 harness、输入/输出 hashes、phase subset/reorder、多层迭代与旧内核对照见 [Eulerian turbulence native migration](../docs/eulerian-turbulence-native-migration.md)。原 production turbulence 入口无法运行，报告不把恢复接线后的诊断执行冒称历史 production baseline。全局 compatibility DSL 仍为 level-set/mixture 等 consumer 保留；IBM 的 native contribution 接线见下一节。

## Current implementation and verification status

实现状态与数值验证范围分别记录；已有 executable path 或结构测试不能代替该路径的 CFD baseline。

| Path / contract | 当前源码状态 | 已有数值验证 |
|---|---|---|
| Single-fluid conservative explicit | WHAT/STATE/HOW/WHICH 编译为 owned fused stage；runtime 执行冻结 provider | Sod 454 steps / 5 VTS，Ghost cylinder 11450 steps / 51 VTS；前后 time/dt、全部输出 hash、积分、极值、范数、closure diagnostics 完全一致 |
| Constant-density PISO/SIMPLE/PIMPLE | native HOW + method-owned pressure/Rhie–Chow operations；无 legacy provider routing | 本阶段 PISO、16²/32² Poiseuille、fixed-time SIMPLE/PIMPLE、既有 MPI 回归前后全部通过；同配置前后 VTS hash 一致 |
| Active STATE / lazy views | catalog metadata 与 active registry 分离；HOW/WHICH 请求 view，runtime 绑定 storage | inactive symbol、qualified target、storage/owner contract 检查通过；本阶段 CFD 实际覆盖 Working(U)、Correction(p)、Physical(U/p/phi) |
| Single-fluid conservative/variable-density pressure | native six occurrences + explicit three-equation predictor fusion；bindPressureOps 仅注册已有 kernels | pressureConstraintPiso 1 step，time/dt=1.659315e-06；输入、诊断、操作顺序、字段统计与两个 VTS 均前后完全一致；fixed-time conservative 等 capability 仍 Unsupported |
| Single-fluid transported RAS | native WHAT/STATE/HOW/WHICH；完整相邻 pair 显式 fusion；CompiledSolvePlan 调用 flow.turbulence，一次/物理步 | 新增 SST 串行覆盖：6 steps / 4 VTS 与旧内核诊断对照逐字节一致；kEpsilon 为结构、存储及更新单元覆盖，未做 CFD baseline |
| Eulerian turbulence | 原生 mass-weighted WHAT、单 registry STATE aliases、显式 closure/transport HOW 和独立 WHICH；无 E_TURB/legacy routing | SST/kEpsilon：每项 6 steps / 7 VTS / 180 operations；subset/reversed phases 及多 PIMPLE passes 严格一致；多 pass 588 operations；Smagorinsky closure-only 174 operations；全部与隔离旧内核接线对照比较，非历史 production baseline |
| Eulerian shared pressure | native phase/pressure WHAT/STATE/HOW/WHICH；explicit N-member fragment fusion；EulerianStepper 仅绑定/执行；AST -> immutable compiled assembly contract，零 native flat DSL 依赖 | 现有双相/MRF/PIMPLE：6 steps、7 VTS、168 Ee operations 和诊断前后逐字节一致；三阶段 CFD 均完全一致；最终完整回归见本次迁移报告 |
| IBM / distributed execution | native boundary/impulse/projection/block contracts；原 numerical kernels 和 COPY/SUM semantics | 各方法分别验证；不能用 Ghost 代表 forcing/KKT；本次完整时域结果见 IBM migration report |

## Next migration sequence

上一阶段 native pressure CFD 验证与 single-fluid conservative pressure authority 迁移已完成。完整 13 项报告、19 项自审、数值与结构证据见 [native pressure migration report](../docs/native-pressure-migration.md)。MPI 同配置前后输出完全一致；串行与 MPI 的比较遵循既有 tolerance，不称为 bitwise identical。既有案例的 nonOrthogonalCorrectors=0，覆盖单 pass；更多 pass 仍 Unsupported。

本轮显式 solution STATE 阶段已完成：原生 Sod、常密度 PISO、fixed-time SIMPLE/PIMPLE 与保守 PISO 每侧 20 次 run 调用（含 mesh 生成），32 个 VTS 比较项全部 hash 一致；字段统计、clock/dt、diagnostics、operation/stage context 序列一致。22 个契约/输入/并行测试、11 个输入组合、完整 build 和架构检查通过。详细 13 项报告及 12 项自审见 [explicit solution STATE migration](../docs/explicit-solution-state-migration.md)。本轮没有新增 ResolvedSimulationSystem 聚合字段；结构计数守卫同步到冻结源码中已存在的 12 个字段，NumericalSelection 是独立 source WHICH，未放宽 CFD 容差。

后续顺序仅作为路线，不在本阶段执行：

Eulerian turbulence 原生迁移已完成，历史完整 build、host 38/38 CTest 见 [中文迁移报告](../docs/eulerian-turbulence-native-migration.md)。IBM native contribution 的 authority 迁移见下节。

1. 补齐 IBM 组合能力：RAS stage/壁距/边界、pressure predictor/block、Eulerian phase ports 与 distributed KKT；不以开放 enabled 伪装 numerical capability。
2. level-set、mixture 等剩余 compatibility consumers；之后仅删除无生产消费者的全局 DSL。
3. 用户自定义四模块 IO：复用同一 compiler、validator、explain。

## Native IBM contributions

IBM 使用现有四契约与同一 compiler。Ghost/ILW 贡献 boundary closure 的 read/write/stage/order contract，physical BC → halo COPY → reconstruction → publication → halo COPY → operator，不伪造 EquationCall。

BP 是 post-predictor momentum penalty 与 midpoint work，没有乘子；Peskin / explicit DFM 消费旧 lagged force 后计算下一层乘子，不宣称当前步精确 no-slip；FTS/fractional DLM 贡献原 mass-norm projection。显式/fractional自推进使用 virtual-fluid generalized mass projection，并区分 active translation/rotation；未改为实体刚体质量时间增量。

原 multiplier/laggedForce/solid velocity 通过 SpecializedExecutor STATE aliases 和 original-storage port 绑定；surface lambda 仍是原数值内核的 transient workspace。native complete mathematical group 经 adjacent fusion lower 一次原 project/KKT callback，所有数学 members 保留 provenance/write-set。compiled immutable contract 校验 AST/backing/target/algorithm/rigid DOFs；native correction 当前只支持完整 implemented predictor 之后、commit 之前，不支持任意 stage treatment。

IBM 的 legacy Definition/HOW/policy/provider routing、空 ImmersedConstraintTransformer、C_IBM 前缀 matching 与独立 descriptor inventory 已删除；descriptor 只保留实际 method options/KKT capability 与 display summary。generic temporal compiler 不按 IBM 字符串选择 execution。

原 ILW/J/Jᵀ、body/surface projection、rigid solve、KKT/HYPRE numerical helpers 保留。三种 KKT 有 structural native contracts，当前 canonical production inputs 新旧均因 recipe/lowering 缺口 fail-fast；不能把它们标为已验证 Runnable。Eulerian＋IBM、RAS/LES＋IBM、pressure/implicit＋IBM 的缺失 mathematical/stage/storage contracts 在编译/capability 阶段明确拒绝。Peskin lagged surface MPI 允许有限非负的 local partial diagonal（包括空支撑的零），Runtime SUM 后严格要求正性；单独验证合法空 rank、canonical lambda COPY 后移除该编译限制。force/mask 是 Eulerian owner 已完成的诊断量，输出 replica 使用 Runtime COPY，不是 SUM 或 average。物理 ILW 的法向和模板遵循 source mesh 的 communication mask；ILW 在重建前声明 conservative halo read，之后保持 physical BC → halo → immersed closure 的相对顺序。IBM 的 mass/spreading/load 使用统一的物理 nodal volume，内部 partition endpoint 不引入半体积，EMPTY 端点共同表示真实厚度。47/47 host CTest 与首差异定位已经完成，7 个串行冻结案例逐字节一致，8 个 MPI 案例全部运行至 endTime=5 且逐帧 output replica 冲突为零；完整时域差异、剩余 kernel 一阶矩/力矩性质与 RAS boundary 缺口见 [IBM parallel consistency validation](../docs/ibm-parallel-consistency-validation.md)。历史迁移 preservation 结果见 [IBM native contributions migration](../docs/ibm-native-contributions-migration.md)。
