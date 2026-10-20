# Equation Registry + Ordered Execution migration report

> 2026-09-30 阶段历史记录。后续严格职责收口已修改 source Target、temporal、lifecycle 与 legacy API，当前状态见 [严格解耦报告](strict-what-how-which-migration.md) 和 [架构正文](../src/ARCHITECTURE.md)。

本报告描述本轮 WHAT / HOW / WHICH 收口后的实际实现。工作区开始时已有大量未提交修改；这里列出本轮涉及的职责，不将工作区所有新增、删除文件归因于本轮。

## 1. Before

WHAT 同时存在 FormulaRegistry、FormulaGroup、平面 Equation::Definition 和 EquationDescriptor/EquationRole。HOW 同时存在 MethodProgram steps、pressureProgram、压力专用 schedule 和 SolvePlan。WHICH 由 method objects、数值 recipe 和 provider 分别绑定。

重复 authority 包括：momentum 与 momentum predictor 的数学定义、FormulaGroup 的向量调用身份、flat steps 与结构化执行树、native 方程与平面 DSL 副本。Builder 同时组装贡献和执行 numerical compilation，压力 schedule 与实际执行结构分别保存循环次数。

## 2. New architecture

```text
MODULES
   │
   ├───────────────────┐
   ▼                   ▼
EQUATIONS           EXECUTION
 WHAT                  HOW
   │                   │
   └────────┬──────────┘
            ▼
         COMPILER ◄──── NUMERICS (WHICH)
            │
            ▼
    executable operations
            │
            ▼
          kernels
```

EquationRegistry 保存数学定义；ExecutionProgram.root 保存调用、target 和局部顺序；NumericalBinding/provider/recipe 实现数值方法。compiler 绑定并 lowering。完整设计正文见 [src/ARCHITECTURE.md](../src/ARCHITECTURE.md)。

## 3. Equation Registry

实际 native density 方程为 `continuity`、`momentum`、`energy`。常密度 NS 为 `momentum`、`continuity`；continuity 的 `div(U)=0` 由压力一致性路径满足。

压力耦合注册 `pSimple`、`correctU`、`correctP`、`correctFluxp`；固定时间关系为 `relaxIterate`、`restoreFlux`、`checkConvergence`。momentum 始终只有一个数学 definition，predictor 引用该 definition。

`Equation {id,lhs,rhs,origin,authored}` 不含 target/order。`EquationRegistry` 提供 add/replace/disable/at/contains/entries，builtin 与 C++ 自定义注册使用同一入口。任意 YAML 数学表达式解析尚未实现。

## 4. Execution model

公共结构位于 `core/system/SF_solveProgram.h`：

- `EquationCall`：EquationRef、typed Target、workspace requirement 和 occurrence。
- `Target`：id、symbols、workspace、kind、编译后 resources；kind 为 Physical/Working/Correction/Workspace。
- `Order`：scope 内整数优先级。
- `ExecutionScope`：Sequence/EquationCall/Loop/StageLoop/Commit，children、order、循环次数、最小迭代次数与 termination signal。
- `ExecutionProgram.root`：唯一 source execution tree；legacyEntries 明确记录尚未迁移贡献。

`orderExecution` 在每个 scope 内递归 stable_sort。相同 order 保留贡献插入顺序，不建立全局 DAG，也不 flatten 后猜测物理顺序。默认 occurrence 使用 scope path；显式 occurrence 可指定 numerical override。数值绑定按全局默认、方程默认、occurrence override 选择；缺失或不支持组合 fail fast。

Native density 默认调用 continuity→rho @10、momentum→rhoU @20、energy→rhoE @60。守恒向量 fusion 由 provider/time lowering 完成，每个 equation call 仍保留来源。

## 5. Pressure coupling

`applyPressureExecution` 消费已有默认 momentum occurrence，将其 target 改为 Working U*，保留 EquationRef 和来源；不会复制 momentum AST。其他模块调用保留。

实际公共 correction body 为：

```text
Sequence correction
  Loop nonOrthogonal (configured count + 1)
    pSimple -> p' (Correction)
  correctU     -> U   (Physical)
  correctP     -> p   (Physical)
  correctFluxp -> phi (Workspace)
  after: correctionCommit
```

实际 PISO：

```text
momentum -> U*
Loop pressure (configured correctors)
  correction body
Commit: pressure stateCommit -> timeCommit
```

实际 SIMPLE：

```text
Loop outer (maximum / minimum iterations, convergence signal)
  momentum -> U*
  correction body
  relaxIterate     -> iterate*
  restoreFlux      -> phi
  checkConvergence -> converged
Commit: pressure stateCommit -> timeCommit
```

实际 PIMPLE：

```text
Loop outer (maximum / minimum iterations, convergence signal)
  momentum -> U*
  Loop pressure (configured correctors)
    correction body
  relaxIterate -> iterate*
  restoreFlux -> phi
  checkConvergence -> converged
Commit: pressure stateCommit -> timeCommit
```

准备、iterationBegin/iterationEnd、correctionCommit 使用 scope before/after lifecycle hooks，降低到原有 operation；没有新增算法专用 HOW node。循环数量和执行顺序来自这棵树。Native 常密度 capability 从 compiled loops 读取次数。原有 relaxation、flux restore、收敛和 commit 的 arithmetic 保留；数值等价性本轮未验证。

`PressureOperators` 的 Working predictor 使用 solver-owned velocity workspace，correctVelocity 后发布到 physical state。Physical target 直接发布。已有 boundary/halo-ready values 通过本地 COPY 发布，没有增加 MPI 通信、Field authority 或 clock。

## 6. Compiler

真实调用链：

```text
SystemCompositionBuilder::applyContribution / addEquation
  collect registry + execution + numerical bindings + pending legacy entries
SF_systemBuilder.cpp
  compose ExecutionProgram
  applyPressureExecution                 [SF_pressureCoupling.cpp]
  compileSystem                          [SF_numericalCompiler.cpp]
    orderExecution                       [SF_solveProgram.h]
    compileExecutionProgram              [SF_methodObjects.cpp]
      resolve typed target resources
      select numerical binding / provider
      compile each EquationCall
      compile immutable time recipe / temporal fusion
    NumericalCompiler::compile
      bind spatial recipes and term providers
    SolvePlanner::compile                [SF_solvePlan.cpp]
      compileMethodProgram               [SF_methodObjects.cpp]
      lower generic scopes to operation plan
    freeze sourceProgram + compiledProgram
PlanExecutor                             [SF_planExecutor.cpp]
  pass each typed Target to operation callbacks
```

Builder 负责 composition，`compileSystem` 负责 numerical compilation 和 lowering。旧模块 fragments 明确作为 legacyFragments 传入，不替代 native execution authority。

## 7. Deleted authority

本轮实际删除/替换：

- FormulaGroup / FormulaGroupRegistry 及公开 group forwarding API。
- native single-fluid density、constant-density NS 的平面数学副本。
- 独立 momentum predictor equation definition。
- pressureProgram 函数，改为变换既有 execution 的 applyPressureExecution。
- source ExecutionProgram 的 flat steps authority，改用 root。
- native 常密度 pressure 的独立 PressureSchedule policy authority；循环 capability 从 compiled HOW 读取。
- 公开 HOW 的 FormulaCall/FormulaMode；它们仅在 formula compiler 的 Legacy numerical compatibility 接口中保留。

没有删除所有旧 pressure fragments、所有平面 DSL 或所有 AssemblyPlan。以下边界仍明确存在。

## 8. Remaining legacy

| 边界 | 实际状态与原因 |
| --- | --- |
| Eulerian | phase legacy equations 和 stepper 保留；phase continuity/momentum/enthalpy 的 ordered entries 明确标记 pending，避免重写其数值 lifecycle。 |
| Turbulence | k/omega/epsilon equation 注册及 ordered pending entries 接入 common contribution；旧 backend 仍执行输运，不重复执行 pending entries。 |
| IBM | Ghost/ILW closure、forcing/projection/KKT 旧执行机制保留；constraint/workspace 的 pending entries 和旧 policy 标明未迁移。没有修改 IBM 算法。 |
| AssemblyPlan | 未迁移 numerical backend 的内部 lowering/兼容入口保留。 |
| EquationRole | 更名 LegacyEquationRole，用于 legacy metadata/capability，不控制 native pressure HOW。 |
| TermKind / flat DSL | 在 legacyDefinitions、legacyEquations 中保留；native equations 已不创建这些数学副本。 |
| legacy provider | FormulaCall/FormulaMode 在 Legacy namespace 作为数值 proof/compatibility 输入保留；不是 public execution API。 |
| variable/shared pressure | 原有 legacyCouplingPlanFragment/sharedPressurePlanFragment 和能力限制保留，未冒充已完成通用迁移。 |

当前 fused conservative backend 明确拒绝尚未实现的混合独立 provider body / iterative temporal 组合。StageLoop 可 lowering，但不是承诺任意 equation、target、stage 组合均已具备 kernel capability。

## 9. File changes

| 主要文件 | 本轮职责 |
| --- | --- |
| core/system/SF_formula.h/.cpp | Equation/EquationRegistry 与数学 AST，删除 group authority。 |
| core/system/SF_solveProgram.h | typed target、occurrence、generic scopes、局部排序、compiled calls。 |
| core/system/SF_equationIR.h、SF_systemContribution.h | native registry、显式 legacy metadata、common contributions 与 mathematical extensions。 |
| solver/system/SF_equationContribution.* | 收集 equations、HOW、numerics、legacy entries。 |
| solver/system/SF_singleFluidPreset.cpp、SF_presets.cpp | native 方程和默认调用。 |
| solver/system/SF_transformation.*、SF_pressureCoupling.* | pressure relations、变换已有 momentum occurrence、通用循环树。 |
| solver/system/SF_methodObjects.*、SF_numericalCompiler.* | provider binding、target resolution、编译 pipeline 与 lowering。 |
| solver/system/SF_systemBuilder.cpp、SF_solvePlan.cpp、SF_providerResolver.* | composition 简化、plan lowering、读取 compiled capability。 |
| solver/run/SF_planExecutor.*、solver/algorithm/SF_singleFluidStepper.cpp | runtime 传递并消费 typed target。 |
| solver/algorithm/pressure/SF_pressureOperators.* | Working predictor workspace 与 physical publication。 |
| models/physics/SF_sourceContribution.cpp | AST source extension，不拥有 HOW。 |
| models/turbulence/SF_turbulenceSystemContribution.cpp、models/ibm/SF_ibmSystemContribution.cpp | common contribution 与明确 pending legacy 边界。 |
| solver/system/SF_systemPrinter.cpp 等 | explain 展示 registry、execution 和 legacy 边界。 |
| test/test_methodObjects.cpp、相关架构测试、tools/check_architecture.py | API 适配、结构断言和静态 authority 检查。 |
| src/ARCHITECTURE.md、src/README.md、README.md、相关历史 audit | 主文档更新与历史状态标注。 |

本轮删除的是上述 authority/type/function，不额外删除生产 numerical kernel 文件。工作区已有的大量文件删除和移动不在这里重复认领。

Ownership：EquationRegistry 拥有 native math，ExecutionProgram.root 拥有 source order/target，compiled plan 拥有执行绑定，StateBundle 仍拥有唯一 physical state/clock，predictor workspace 属于 solver execution lifetime。

Execution order：压力专用 schedule 改由 generic scope 显式表达；kernel 内原有 arithmetic 顺序保留。Public API：Formula/group/program step 改为 Equation、EquationCall、Target、ExecutionScope、NumericalBinding。Dependency：composition→compiler→plan→runtime callbacks→原 kernel；无新增 raw MPI 依赖。

并行语义：owner/COPY/SUM、canonical face 和 halo kernel 未在本轮改变。Numerical semantics：没有有意修改数学公式、RK coefficients、CFL、tolerance；Working storage/publication 改动尚未以 numerical baseline 证明等价。

Deferred：迁移剩余模块 runtime、任意 YAML equation parsing、更多 target/stage/provider capability，以及后续独立 numerical regression。

## 10. Build status

- 现有 build 配置下全部 configured targets 编译/链接通过，包括 sonicSolver 和测试 binaries。
- `test_methodObjects` 结构单元检查通过，包含 repeated equation occurrences、scope-local order、numerical override、typed target 和 fake operation callback；没有 CFD 时间推进。
- `python3 tools/check_architecture.py --quiet` 通过。
- `git diff --check` 通过。
- native density、constant-density PISO、SIMPLE、PIMPLE 和 Eulerian 的 CLI explain 解析为 runnable。SIMPLE/PIMPLE 使用临时复制 case；没有改动原 case 来执行计算。

Build log：`/private/tmp/sonic-what-how-build-all.log`；explain logs：`/private/tmp/sonic-explain-{density,piso,simple,pimple,eulerian}.log`。

Numerical regression intentionally not run in this architecture phase.

本轮未运行 CFD timesteps、长算例或完整 numerical regression，也不声称 numerical equivalence。其他测试 targets 的编译通过不代表其所有 runtime assertions 已执行。
