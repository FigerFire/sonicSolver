# 显式 solution STATE 收口报告 — 2026-10-02

本轮实现原生 WHAT + STATE + HOW + WHICH 的显式组合。主未知量成员身份与数学角色分开；物理守恒密度及 balance 留在 WHAT。没有增加第五个模块、第二个物理状态 authority 或 solver family。

## 1. 修改前的实际架构

`SF_systemBuilder.cpp::validateSemanticComposition` 曾把 rhoConst/PerfectGas 的实现限制写成方程族规则；`build` 根据 equation-source/EOS 分支安装状态包。`SF_singleFluidPreset.cpp::installSingleFluid` 无条件注册 rho/rhoU/rhoE。`SF_presets.cpp::addConstantDensityFluid` 无条件注册 U/p 和常数 rho。`PressureConstraintTransformer` 用 rho.constantValue 是否存在选择 primitive/conservative 校正关系。

调用关系因此包含 `EOS/preset → state pack → physical equations/default HOW`。旧标签还在 `Application::inspectCase` 选择 pressure-constraint 或 conservative preset。

## 2. 新 STATE 输入

实际文件为 `state/state.yaml`，通过现有 CaseIO/IORegistry/FactoryRegistry 读取：

```yaml
SonicFile:
  object: state
  type: registry
use: [U, p]
```

保守表示使用 `use: [rho, rhoU, rhoE]`。相邻输入分别是 `equations/equations.yaml` 的 `use` 和 `algorithms/algorithms.yaml` 的 Explicit/SIMPLE/PISO/PIMPLE 对象。闭合关系仍在 models/thermoDynamics.yaml。

输入落到 `EquationCompositionConfig::stateDeclared/solutionVariables/stateSelectionOrigin`；冻结后 `StateRegistry::solutionVariables()` 保存成员身份。原生输入缺 WHAT、STATE 或 HOW 会报错；不从 EOS、fields 文件或旧 type 标签补 solution。

旧控制文件中的 linearSolvers、relaxation、reference 等参数继续由已有 decoder 读取。原生 HOW 对象决定算法；旧文件中的算法名仅帮助读取既有参数块。两者都提供参数时原生参数覆盖，线性容差与 outer convergence 容差保持独立嵌套作用域。C++ BuildRequest 的显式算法与 coupling contribution 不一致也会失败。

## 3. Builtin catalog

`BuiltinStateCatalog` 保持 solver-owned。`solution(id)` 与 `require(active,id)` 分别处理选中的可写变量和依赖：

| 符号 | solution metadata | dependency metadata |
|---|---|---|
| U | Primary，packed velocity，3 components，无 derivation | Derived velocity view |
| p | Algebraic，named pressure，无 derivation | Derived EOS pressure view |
| rho/rhoU/rhoE | Primary，conservative packed offsets 0/1/4 | 按方程/closure 请求的依赖 |
| T/h 等 | 有对应 builtin metadata；执行还需实际 provider | 按需 EOS/closure view |

用户选择名称并提供现有初始化/边界/模型参数，不填写 shape、storage key、component offset 或 ownership。C++ `SystemContribution::addState` / `SystemCompositionBuilder::addState` 继续支持真正自定义变量；本轮没有实现任意自定义 STATE YAML。

p 不在 catalog 中永久等同于 multiplier。`compileStateRealization` 从实际 constraint relationship 判断 algebraic p 是否承担 multiplier 职责。solution 成员身份不从 Primary/Algebraic/Multiplier/Derived 推断，也没有新增互斥的 Conserved role。

## 4. STATE 激活

`System::build` 首先用显式选择激活 writable builtin metadata，并记录 selection provenance。model/user 的 requireState 仅激活 dependency；非 Workspace HOW target 的 requireTargetStates 仍按需请求 base metadata。

rhoConst 在未选择 writable rho 时激活 constant/derived rho：rho0 直接冻结为 constantValue，没有可写数组，也不作为 restart authority。选中 writable rho 时不会覆盖成常数或换成 U/p，而是由缺失的 closure-consistent evolution provider 拒绝。

Raw/ExecutableEquationSystem 的 STATE 是编译快照；运行时物理数组和 clock 仍由 StateBundle/已有 Field authority 持有。Working、Correction、OldTime、Stage 和 Workspace 继续由 HOW/WHICH 请求，不能成为另一个 solution registry。

## 5. WHAT 的数学表示

`installSingleFluid` 不再注册 solution state pack，消费显式 STATE 和 equation selection。canonical identities 仍是 continuity/momentum/energy。

保守表示仍是 ddt(rho)、ddt(rhoU)、ddt(rhoE)。选中 U 后，当前 `addConstantDensityFluid` 在 rho=rho0 的真实数学条件下构造 ddt(U) 与 div(U)=0；这包含对 momentum 的常密度归一化和 continuity 的归约，不是字符串替换 ddt(rhoU) → ddt(U)。variable-density primitive 表示尚无对应 momentum/continuity/pressure-response provider，进入不支持的 lowering 前就报 capability 错误。

原生保守压力示例 WHAT 显式列出 `[Continuity, Momentum, Energy, PressureConstraint]`。PressureConstraint 是已有 pressureVelocityConsistency 数学约束的库入口，不是 solver family。HOW 不静默增加或重定义物理约束。

默认 Momentum target 从显式 solution U/rhoU 解析。双选在当前 preset 下属于歧义，明确失败；原生 YAML 当前没有实现任意 target override。缺少合法 target 同样失败，不能依赖 catalog 中恰好存在同名 dependency。

## 6. HOW 独立性

Explicit/SIMPLE/PISO/PIMPLE 必须显式选择。`composeContributions` 继续从默认 momentum occurrence 的 target 取得 U 或 rhoU，再调用 `applyPressureExecution(program,coupling,target)` 和 `pressureNumerics(coupling,target)`。

`makePressureConstraintTransformer(momentumTarget)` 接收同一个已解析 target，选择相应 pressure relation set；不再检查 rho.constantValue 来决定表示。STATE 本身没有 equation destination、order、Loop、corrector count 或 time recipe。

压力算法的 loop、correction、relaxation、flux restore、convergence 与 Commit 顺序维持旧实现。保守 predictor 同时写 rho/rhoU/rhoE，compiled ownership 校验读取 method 声明的完整 writes/resource contract，不能只看外层 rhoU target。

## 7. Closure 独立性

删除原 native validator 中 rhoConst→无 Energy/压力算法、PerfectGas→Energy/仅 Explicit 的族式判断。保留已声明 EOS/caloric/transport 标识的存在性校验、Energy 所需 caloric closure 和实际 provider 能力要求。

`validateSolutionProviderContract` 只拒绝当前缺失的实现，不选择任何 replacement。原生目前执行已证明的 PerfectGas/hConst 保守表示和 rhoConst/const primitive 表示。注册表能识别 JANAF/Sutherland 名称，但当前单流体 provider 未实现它们，因此明确拒绝，不执行默认模型冒充所选模型。

本轮没有扩展通用 EOS/thermodynamic property binding。既有参数和初始化通道保持，不能把名称注册或框架表示误认为任意 thermophysical provider 已可执行。

## 8. WHICH 独立性

NumericalSelection 与 CompiledNumericalSystem 继续分别保存 source selection 和 lowered numerical bindings。ConservativeResidual、PressureMomentum、pressure/Rhie–Chow 与 conservative pressure methods 的数学/kernel 未改。

WHICH 不改 source HOW、STATE 或 EOS。时间仍选择不可变 built-in recipe；没有恢复 TimeScheme、explicitStageCount 或 ddtDispatch。原有 WENO/primitive reconstruction、stage treatment、nonOrthogonal、single-patch 等能力校验继续生效，不能隐式换方法。

## 9. 实际编译调用链

```text
ModelLoader::read
  → CaseIO::read / IORegistry / SonicFile serialization
  → CaseAdapter::build
      equationRegistry + stateRegistry + algorithmRegistry + thermoDynamics
      → decodeNativeCase（复用数值控制、初始化和 BC decoder）
  → Application::inspectCase → BuildRequest
  → System::build
      validateSemanticComposition / 显式 HOW consistency
      BuiltinStateCatalog::solution → provider representation capability check
      installSingleFluid / addConstantDensityFluid / explicit constraint WHAT
      model/user contributions → requireState → selectSolution
      TransformationPipeline::apply（明确 momentum target）
      composeContributions → requireTargetStates → selectNumerics
  → NumericalCompiler::compileSystem
      compileStateRealization
      temporal method compile
      compileExecutionProgram
      校验 selected solution 的完整 write owner
      numerical operator/provider compile
      SolvePlanner::compile
      compileExecutionCapabilities → compileOperationBindings → runtime report
  → System::validate
  → owned OpRegistry / PlanExecutor → existing kernels
```

表示本身尚无数学 provider 的组合在 `validateSolutionProviderContract` 提前失败；不会伪造完整 executable plan。支持表示的 method/operator/stage 能力继续在后续 compiler/provider validation 处理。

## 10. 已移除推断与保留的兼容边界

| 原逻辑 | 本轮处理 |
|---|---|
| installSingleFluid 无条件添加 Q | 删除；STATE 独立激活 |
| addConstantDensityFluid 无条件添加 U/p | 删除；仅构造已选择表示的数学 |
| EOS 分支选 state/equation pack | native builder 改读显式 STATE；EOS 仅 closure/capability |
| constant density 选择 pressure relations | transformer 改为显式默认 momentum target |
| inspectCase 读取 densityBase/pressureBase | 删除；消费 adapter 已翻译的 typed requests |
| native HOW 缺失时靠旧标签/默认算法 | 明确失败 |
| solution 成员身份等同 StateRole | 增加独立 membership；role 保留数学职责 |

旧 case 无原生语义声明时，CaseAdapter 一次性翻译其显式 legacy 标签，生成带 compatibility provenance 的 WHAT/STATE/HOW tuple。冻结前的旧 Sod 和旧 pressureConstraintPiso 用最终 binary explain 均通过，native leaves 不走 legacyAdapter。原生部分声明缺 STATE/HOW 时不会进入兼容补齐。

原生 tuple、model prerequisites、NumericalSelection 和 immutable time recipe 没有新 Manager/Context/Facade；没有新增 ResolvedSimulationSystem 聚合字段。旧结构计数守卫为 11，冻结源码已经为 12（已有 source WHICH NumericalSelection）；已按这份实际输入同步 baseline 并记录明确架构决定。聚合头文件逐字节未变，CFD 容差与 frozen numerical baseline 未改。

## 11. Regression 与验证

数值基线在修改前冻结，比较使用原有 checker/tolerance。fixed-time/distributed 测试 helper 只补充原生 HOW 选择；数值参数、断言和 CFD 容差不放宽。

| 回归 | 每侧 run 调用数（含 mesh） | VTS hash 比较项 | 结果 |
|---|---:|---:|---|
| Sod WENO7/Rusanov，t=0.2 | 1 | 5 | 全部一致；原 checker 的 mass/momentum/energy 积分及 rho/p/speed ranges 通过 |
| constantDensityPiso | 2 | 3 | 全部一致；checkerboard/gauge/continuity 通过 |
| fixed-time SIMPLE/PIMPLE cavity + SIMPLE Poiseuille | 16 | 22 | 全部一致；相同时钟、outer/convergence/relaxation/zero gauge 通过 |
| pressureConstraintPiso | 1 | 2 | 全部一致；保留原 12-operation sequence |

Sod 最终积分为 mass=98.67807117165903、momentum=[23619.052634273423, 1.907691669240853e-08, 0]、energy=21770158.12505413；前后完全相同。

全部 32 个 VTS hash 比较项、字段 min/max/L2、clock/dt、closure/continuity diagnostics、operation sequence 与逐行 SF_PLAN stage context trace 相同。SIMPLE Poiseuille L2=0.00822508807、Linf=0.0112953895、maxDiv=1.37627e-09、Δp=0.1。fixed-time cavity 时钟仍为 (0.000390625,0.000390625)、(0.00078125,0.000390625)。保守 PISO 最终 VTS SHA-256 为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`。

完整 build、22 个 unit/input/parallel contracts、architecture dependency/scope/fieldcount 检查、git diff --check 通过；原生输入检查覆盖 11 个组合。113 个 numerical/boundary/time/MPI 相关源码文件（含原 equation-method kernels 和 pressure relation formulas）与本轮冻结源码逐字节一致。

本轮四组 CFD 为串行基线；并行验证覆盖 MPI primitives、fail-stop 和 halo 契约，未新跑完整 MPI CFD、Ghost/forcing/DLM/KKT 或 turbulence/Eulerian CFD。没有新增每个 candidate/canonical flux 或内部 stage Q 的 dump；不把已有测试覆盖范围扩称为所有数学组合已验证。

最终 binary hashes、逐次输入/输出、stage trace、比较结果及 build/unit/static logs 见 [evidence manifest](evidence/explicit-state-20261002/manifest.json)。第一次验证中的参数作用域、fusion write-owner 和空 mapping serialization 问题均已修复；失败尝试保留在临时冻结目录，不覆盖最终结果。

## 12. 已验证的 Unsupported/Invalid 边界

`check_solution_state_selection.py` 从实际原生 case 经过 production reader/compiler 验证：

| 组合/缺项 | 实际错误依据 |
|---|---|
| rhoConst + [U,p] + SIMPLE | Runnable；无 type family label 也可解释 |
| rhoConst + [rho,rhoU,rhoE] + Explicit | 缺 closure-consistent density evolution / conservative energy provider |
| PerfectGas + [U,p] + PISO | 缺 variable-density primitive momentum/continuity/pressure-response provider |
| [U,p,rhoU] | Ambiguous Momentum target；不选其中一个 |
| 缺 STATE / 缺 HOW | 必须显式声明；不从 EOS 或旧标签补齐 |
| 重复 STATE | Duplicate solution STATE |
| [U,p,T] 无 T equation owner | No compiled HOW output owns selected solution STATE 'T' |
| 未声明 custom metadata | Solution STATE has no metadata: custom |
| JANAF / Sutherland | 当前单流体 caloric/transport provider 缺失；不默认执行 hConst/const |

现有 RK-pressure stage treatment、WENO 对 primitive contract、conservative SIMPLE/PIMPLE、额外 nonOrthogonal、MPI conservative pressure、auxiliary physical modules 等限制仍按实际 method/provider 能力报告。没有扩大这些执行能力。

## 13. 剩余迁移与本轮影响核对

后续仍为 turbulence specialized execution、Eulerian/shared-pressure fragment、IBM boundary/source/constraint adapters，以及完整自定义 WHAT/STATE/HOW/WHICH IO 和 thermophysical provider binding。下一阶段未开始。

AGENTS 要求的影响核对：

1. 层：input/composition、STATE metadata、compiler/provider validation、explain/serialization。
2. 分类：改变方程表示的配置/编译方式；显式公开已有 constraint 库入口，没有新物理项或新数值算法。
3. ownership：原 EOS/preset 决定 state pack → 用户 STATE 决定 membership；物理数组/clock 仍 StateBundle，workspace 仍 execution/provider。
4. execution order：旧 kernel 顺序 → 同一冻结 HOW/recipe/Plan 顺序；全部 stage trace 验证相同。
5. public API：BuildRequest 必须携带显式 STATE；installSingleFluid/pressure preset 消费 composition；pressure transformer 接收 target；StateRegistry 增加 membership 查询。
6. dependency：复用既有 core config/state 与 solver compiler/provider；新增 stateRegistry 使用既有 IO/factory，无新反向依赖。
7. parallel：COPY/SUM/canonical face identity 未改；并行契约测试通过。
8. numerics：kernel、RK、CFL、pressure/Rhie–Chow、EOS arithmetic 与 CFD 容差未改。
9. regression：上表四组 frozen CFD、22 contracts、11 native combinations 与 build/static 检查。
10. deferred：上述未迁移模块、任意 target/representation、通用 EOS/property binding 和对应 CFD 验证。

### 12 项自审

1. rhoConst 是否仍自动选 U/p？**原生否**；只在已选 U 后做常密度数学归约，缺 p 失败。
2. PerfectGas 是否仍自动选 Q？**原生否**；Q 来自 stateRegistry；旧标签 tuple 只在 compatibility adapter。
3. SIMPLE/PISO/PIMPLE 是否选 STATE？**否**；它们变换 HOW 并验证 target/constraint。
4. STATE 是否产生 Loop/order？**否**；membership/catalog/state views 不调度。
5. WHICH 是否改 source HOW？**否**；lowering 编译独立副本，source choices 保留。
6. compiler 是否发明 state representation？**否**；未实现表示提前 capability fail，支持表示绑定显式 target。
7. explain 能否区别 solution 与 dependencies？**能**；在 STATE 首段分别打印 use、require 和 origin，再列 compiled views。
8. U/p/rho/rhoU/rhoE metadata 是否 solver-owned？**是**；BuiltinStateCatalog 提供 writable/dependency contracts。
9. 用户是否只需选择 builtin names？**是**；另提供现有 initialization/BC/model parameters，无 storage/shape/ownership 注册。
10. rhoConst + U/p + SIMPLE 是否复现旧结果？**是**；fixed-time cavity/Poiseuille 全部 VTS、时钟和 diagnostics 相同。
11. Q + PISO 是否复现旧结果？**是**；两个 VTS 和原 12-operation sequence 相同，最终 SHA 相同。
12. 跨组合是否因 capability 失败而非 solver family？**是**；实际 negative tests 的错误对应缺失数学/provider 或明确输入歧义，无 STATE/HOW/WHICH 替换。
