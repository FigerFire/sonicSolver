# WHAT / HOW / STATE / WHICH 框架收口

## 两份建议的取舍

两份建议都正确地把 HOW 定为唯一 scheduler。采用第一份的轻量输入：STATE 只注册 base symbols；采用第二份的运行时分层：StateRealization 绑定实际存储、旧值、工作值和修正量。输入注册表不预声明所有版本，编译输出才列出实际使用的 views。

STATE 可以保存变量的物理 storage/closure contract，但不能保存 equation → target binding。`momentum → U*` 完整地属于 HOW。WHICH 声明 provider 与额外 numerical storage demand，compiler 统一完成状态 binding。Views 不改变 HOW 顺序。

## 唯一 authority

| 模块 | 类型/入口 | 决定什么 |
|---|---|---|
| WHAT | EquationRegistry / Equation | 数学左右侧与方程 identity |
| STATE | StateRegistry / StateSymbol | base variable、类型、physical storage、closure |
| HOW | ExecutionProgram / EquationCall | EquationRef、target qualifier、order、Loop、Commit |
| WHICH | NumericalSelection / NumericalBinding | equation occurrence 的 method、operator/time recipe |
| Compiler | compileExecutionProgram / realizeTarget | lookup、验证、target binding、lazy view manifest、lowering |
| Runtime STATE | StateRealization | 绑定 physical authority 与 execution-owned storage |

Raw/ExecutableEquationSystem 保留独立 WHAT/STATE channel 的 composition snapshot；不是第二份 runtime state。Program 增加 STATE 的只读引用，仍不拥有物理数组或 clock。

## 基础变量目录

Solver 内置 typed `BuiltinState` 目录，注册 rho、rhoU、rhoE、U、p、T、h、phi、alpha、k、omega、mu、nu、conductivity。Preset/model/user 的明确声明优先，目录只补缺。变量存在不自动增加方程、激活模型、分配数组，也不推断 formulation。

例如 rhoConst case 注册表知道 rhoE，但 rhoE 是未激活的可用符号；该 case 不增加 energy equation，不申请 rhoE 的 physical state。k/omega 的数值存储仍由启用的湍流模块贡献。用户不需要重复声明这些常用符号，仍需提供适用的初始化、边界和模型参数。

STATE 的 `addState(StateSymbol)` 是自定义变量接口，替代原 `addUnknown`。不存在 `if (name == "phi") autoRegister...` 的 cell-loop 或用户字段解析分支。

## 编译链

```text
preset/model/user contributions
    ├── WHAT: equations / mathematical transforms
    ├── STATE: base symbols / physical contracts
    ├── HOW: ordered equation occurrences / scopes
    └── WHICH: numerical selection
             ↓
WHAT lookup + STATE lookup + WHICH lookup for every HOW call
             ↓
CompiledTarget + CompiledStateView manifest + compiled numerical fragment
             ↓
Plan execution + runtime state/view binding
```

编译器先验证 WHAT 与基础 STATE symbol，再调用 provider，最后统一绑定 storage。Physical、Working、Correction 都不能绕过 STATE lookup。Workspace 是数值 bookkeeping 的显式内部 target，不作为 physical base symbol；SIMPLE 的 iterate/converged 属于此类。

同一 base/kind 只有一个 backing；冲突 provider demand 会 fail fast。Source HOW 不保存 storage key、component offset、pointer 或 MPI layout。

## Views 与实际运行时

| View | 来源 | storage authority |
|---|---|---|
| Physical(U/p/rho/...) | 变量契约 + HOW | 原 StateBundle/Field/ScalarField |
| Working(base) | HOW `*` | generic compiler workspace 或 provider 原有 predictor array |
| Correction(base) | HOW `'` | generic compiler workspace 或 provider 原有 correction array |
| OldTime(base) | 多 stage temporal recipe | explicit backend 原 q0 |
| Stage(base,s) | temporal recipe | backend published stage state 的限定别名 |
| Physical(phi) | pressure HOW | backend 原 canonical face flux |

Generic Working 初始化为可用 physical view 的快照；Correction 初始化为零增量。Direct/Linear cell binding 根据 typed target 读写对应 view；写入临时 view 不发布到 physical state。发布/commit 必须由 HOW 中的数学关系和 provider 执行。

PressurePrepare 将原 workingVelocity、correction、correctedFlux arrays 暴露为 Working(U)、Correction(p)、Physical(phi)，不复制 authority。Native transformation 不注册 pPrime 基础变量；AST 中 pPrime 是 backend correction 数学记号。Legacy pressure adapter 仍保留自己的数学 auxiliary descriptor。

RK stage 数仍由 immutable built-in recipe 唯一决定。编译器按实际 transported calls 列出 stage views；runtime 绑定已有 backend arrays，不另外分配一套 RK lifecycle。Stage view 只能在对应 active stage 内读取；它不承诺保存所有历史 stage 快照。ForwardEuler 无多 stage q0 需求。

U/p/T/h 等 derived values 通过 Field 的真实 thermodynamic closure 读取。没有合法 closure 时读取失败；不会自动选择 PerfectGas，也不允许把 derived physical target 当成可独立写入的变量。

## 普通符号无需用户注册的 C++ 接口

```cpp
ExecutableEquationSystem what;
StateRegistry state;
installBuiltinState(state); // composition normally does this
what.registry.add({"pressureDelta", FormulaExpr::symbol("p"),
    FormulaExpr::constantValue(7.0), {}});
ExecutionProgram how{EquationCall{"pressureDelta", targetFromSyntax("p'")}};
std::vector<NumericalBinding> which{{"pressureDelta", "DirectEvaluation"}};
auto program = compileExecutionProgram(what, state, how, which, builtinProviders());
// program.stateViews contains only Correction(p); there is no addState(p') call.
```

Production SystemBuilder 自动安装目录，无需用户调用 installBuiltinState。单独的 compiler API 保留显式 STATE 输入，便于用户自定义数学系统和扩展 providers。

## 影响报告

1. **层**：core system specification、solver composition/compiler、runtime state realization、explain。
2. **类型**：改变 equation 的组织/求解实现方式；不增加物理 term 或物理 equation。
3. **Ownership**：旧 unknown vector → StateRegistry；provider 自行绑定 target → compiler STATE binding；temporary arrays 归 execution，physical authority 仍为 StateBundle。
4. **Execution order**：仍由同一 HOW scope tree 排序；STATE 不插入 equation call。Native phi target 从 Workspace 明确为 Physical；数值操作顺序不变。
5. **Public API**：StateSymbol/StateRole/addState/StateRegistry；显式 compiler STATE 输入；runtime realizeState 接收 registry 与 view manifest；Program 增加只读 STATE 引用。
6. **Dependency**：WHAT EquationRegistry 不依赖 scheduler/storage；HOW 引用 base symbol；WHICH 提需求；compiler 负责 STATE binding；runtime 消费冻结结果。
7. **MPI**：不新增 raw MPI，不改变 COPY/SUM、canonical face 或 GlobalDof 规则；view binding 不做同步。
8. **Numerical semantics**：不改 RHS/flux/pressure algebra、RK coefficients、dt 或 stage times。Working/Correction 严格隔离，stage lifetime 明确验证。
9. **Verification**：见下方本轮实际结果；结构/unit 与 explain 不是 frozen CFD 数值基线证明。
10. **边界**：任意 YAML equation/HOW parser、通用 multi-patch cell GlobalDof provider、未实现的任意方法组合仍需相应 backend capability。Eulerian/turbulence/IBM 的 explicit legacy adapters 保留；未宣称其所有局部算法已重写为 generic cell FormulaCompiler。

没有新增 fallback、clamp、降阶或 solver family。仓库原有未提交改动保留；本轮前快照位于 `/private/tmp/sonic-four-module-before`。

## 本轮实际验证结果

- `cmake --build build -j 6`：全部配置 targets 构建成功；日志 `/private/tmp/sonic-four-state-build.log`。
- `test_methodObjects`、`test_formulationArchitecture`、`test_pisoArchitecture`、`test_formulaExecution`、`test_programDecomposition`、`test_termRecipes`、`test_explicitStages`、`test_fixedTimePressure`：全部通过。
- 新增实际 Field/storage 检查：同一 equation 写入 Working 与 Correction，不改 physical；缺失 base symbol 被拒绝；STATE 拒绝注册带 qualifier 的 U*；builtin p 无用户 addState 可编译并写入 Correction(p)；显式独立 WHAT/STATE compiler 输入；RK4 自动要求四个 Stage view。
- `python3 tools/check_architecture.py --quiet` 与 `git diff --check`：通过。
- Density、PISO、SIMPLE、PIMPLE、Eulerian 的 `explain`：均退出成功并报告 runnable，输出包含 STATE 基础符号与实际 view manifest。日志 `/private/tmp/sonic-four-{density,piso,simple,pimple,eulerian}.log`；SIMPLE/PIMPLE 使用此前 `/private/tmp/sonic-how-{simple,pimple}` 的 PISO case 配置变体。

上述验证覆盖编译、存储隔离、接口、recipe/plan 结构与可执行 capability；没有运行完整 frozen CFD 数值基线对比，也没有据此宣称 serial/MPI 数值等价已经重新验收。
