# SonicSolver：双数学入口、统一方程驱动的 CFD 求解器架构

## Math / Typed WHAT 的职责更新（2026-10-21）

本节更新核心设计与实施顺序，优先于较早的“独立数学分析层”和先推进内置控制的规划。2026-10-21 的 Phase A1 已在生产构造、注册和选定数值绑定路径接入统一 Math；不据此声明通用离散编译器或所有数学表达式已有执行能力。实现、逐例设置与验证边界见 [A1 报告](../docs/math-what-numerics-a1.md)。

**Math 定义数学，WHAT 组织方程，Numerics 实现数值算法，WHICH 选择实现；STATE 提供符号与状态契约，HOW 定义求解过程。** Math 是公共基础库，不是第五个 authority。

| 部分 | 职责 | 不承担 |
|---|---|---|
| Math | 数学算子身份、签名传播、量纲、张量规则、纯向量/张量代数 | 网格、离散格式选择、Field lifecycle |
| WHAT | 引用 Math 构造 Typed AST、Equation、Constraint；保留身份、引用、occurrence、来源 | 第二套类型规则、数值求值 |
| STATE | 注册符号签名、位置、存储/视图与版本契约 | 算子类型规则、方程排序 |
| Numerics | 离散算子、融合算法及计算 kernel；声明支持的签名和资源 | 重定义数学规则、替用户选择方法 |
| WHICH | 显式选择 Numerics 实现及绑定 | 重写离散算法、静默替换选择 |
| Compiler | 解析符号、调用 Math 校验入口、规范化和能力匹配、lowering | 独立 MathAnalyzer、推断用户未选的算法 |
| HOW / Runtime | HOW 的调用/耦合/时间声明经编译后由 runtime 执行 | runtime 再决定数学或数值方法 |

### 数学签名与注册

签名核心为 `Signature = (rank, Dimensions)`；七维量纲顺序固定为 **质量、长度、时间、温度、物质的量、电流、发光强度**，元素是指数。另保留必要的指标尺寸、分量布局、对称性信息；rank 相同不是充分的类型兼容条件。STATE 的空间位置、ownership 与版本独立于量纲，仍参与编译/binding 校验。

例如 T 为 rank 0、[0,0,0,1,0,0,0]；U 为 rank 1、[0,1,-1,0,0,0,0]；rho 为 rank 0、[1,-3,0,0,0,0,0]；p 为 rank 0、[1,-1,-2,0,0,0,0]。

注册/解析边界必须取得已知签名；未知量纲与明确无量纲严格区分。常数也有量纲，不把带单位系数或零值默认为可匹配任意量纲。未知签名不能进入已验证 Typed AST；前置引用先解析再构造，不在每个 timestep 重新推导。签名只存于符号/表达式元数据，不为每个 cell 增加元数据。

### Math 唯一定义传播规则

以下规则描述合法输入的输出，仍须检查指标尺寸与结构：

| 运算 | rank | 量纲 |
|---|---|---|
| A ± B | 输入兼容，输出保持 | 输入必须相同 |
| scalar × A / A ÷ scalar | 保持 A | 指数相加 / 相减 |
| outer(A,B) | r+s | 指数相加 |
| 单指标收缩 dot(A,B) | r+s−2 | 指数相加 |
| cross(a,b) | 支持的三维向量输入输出 rank 1 | 指数相加 |
| grad(A) | r+1 | 除以长度 |
| div(A) | r−1，输入 r≥1 且收缩合法 | 除以长度 |
| laplacian(A) | r | 除以长度平方 |
| ddt(A) | r | 除以时间 |
| 当前固定空间方向偏导 Partial | r | 除以长度 |

普通乘法、点积、外积、双点积/多指标收缩必须有明确不同的数学语义，不用一个 Multiply 猜测。`diffusion(D,A)` 应根据其实际定义传播，例如标量 D 的 `div(D grad(A))` 输出量纲为 [D][A]/L²；不能当作无系数 Laplacian。

方程两侧与同一加减表达式的各项最终签名必须兼容；乘积中的每个操作数不要求同阶。同样，连续向量通量的 div 与已投影/积分的离散 face flux balance 不强套一个签名，离散投影、测度和位置变化要明确表达。

内置算子枚举不是封闭数学全集。扩展算子通过明确身份、签名规则及所需约束接入 Math；相应 Numerics 可以尚未实现，此时数学成立但不可执行。不要用匿名不透明回调绕过签名检查。

### 构造期校验与数值绑定

`Expr` 运算构造调用 Math 的签名规则；Equation 构造调用同一兼容检查。WHAT 不另写 MathAnalyzer，provider 不重新实现通用量纲/张量检查。Compiler 可在反序列化、数学改写和入口防御校验时复用这些规则，不维护另一套扫描算法或把未校验 AST 标成 typed。

例：`ddt(rho*U) + div(rho*outer(U,U)) = -grad(p)` 各项最终均为 rank 1、[1,-2,-2,0,0,0,0]。这只证明数学类型兼容，不证明存在可运行离散。

WHICH 按算子 occurrence 与用户选择绑定实现：例如 momentum.pressureGradient → Central2。**不是全局要求 rank=1**；rank=1 仅在所选实现接受向量输入/输出的相应位置成立。完整匹配还包括算子、输入/输出结构、位置、layout、时间、边界与并行能力。允许融合 provider 覆盖多个方程/算子，不强制逐算子拆 kernel。

数学错误为 Invalid；数学合法但无兼容实现/资源为 Unsupported；数学、数值、STATE、HOW 和资源均满足后才为 Executable。Runtime 继续处理除零、非有限值、迭代失败等动态问题，不能由静态签名替代。

`dot(Vector3,Vector3)` 等纯代数可在 Math 直接求值；`grad(Field,Mesh)` 属于 Numerics。表达式 dot 构造 AST，数值 dot 计算分量，二者共享语义但不强迫数值热路径经过 AST。

### A1 的实际落点与边界

- `methods/math/SF_signature.h`：`Dimensions`、`Signature`、`Math::Op`、`resultSignature()`；只依赖标准库。量纲指数统一按绝对容差 1e-12 比较。
- 现有 `core/system/SF_formula` 是唯一 AST；私有结构和 const getters 禁止直接更改操作数。`withArguments()` 经原工厂重新推导类型；Equation 构造及注册/替换校验两侧。canonical key 包含签名，不包含 provenance。
- `FormulaSymbols` 和 STATE 的 `expression()` 提供已声明类型。普通非零数字常量明确无量纲；带单位常数必须声明。只有语法上的 literal zero 在加减/方程两侧获得对方的具体类型，随后固定为 typed zero；显式声明了量纲的零不具有此特例。
- 物理符号和扩展关系归属实际模块：IBM 的 `SF_immersedMathematics.h`、湍流的 `SF_transportMathematics.h`、level-set 的 `SF_levelSetSystemContribution.h`。具体 numerical provider 可读取这些数学声明并独立核对支持范围；通用 compiler 不创建模型，不按模型名调度。公共 core 不再维护 domain 清单；纯 Math 不认识 rho、k、IBM 等物理名字。
- 连续 `div` 要求 rank≥1；积分面率使用明确的 `faceVolumeDivergence`、`faceMassDivergence` 等既有离散关系。异质 iterate、几何测度采用类型 product，不能冒充无量纲 scalar；当前 product 最多四个非嵌套成员，函数签名最多三个输入。这是现有数学契约的受限元数据，尚非通用类型语言。
- `NumericalCompiler` / `TermProviderCatalog` 和 `FormulaOperatorCatalog` 复用 Math identity + occurrence + 原 provider matcher。AST pattern binding 由所选 provider 的 lowering 承担，公共 compiler 只组织调用与绑定校验；`solver/discretization` 的 Central2/feedback kernels 不 include WHAT。原融合 flow、pressure、turbulence、IBM kernels 保留。
- native 显式数学字段须声明 `dimensions`；省略不再解释为无量纲。保留现有有限语法与单目录注册身份，不自动升级历史输入。已声明数学合法不保证存在 provider。
- STATE 的 Q、clock、storage、view/halo 生命周期及 HOW topology 没有迁移。本轮没有新增通用 implicit/Jacobian、非线性源项 kernel、任意多 provider stage 或模块组合能力。

### A2：结构提取、能力匹配与显式非线性 scalar（2026-10-21）

本轮实际实现见 [A2 报告](../docs/math-what-a2.md)。Typed WHAT → `inspectFormula()` 的不可编辑 AST 索引/依赖/target-relative degree → 所选 provider 的能力检查 → `compileFormulaValue()` → 原 StageGroup/PlanExecutor。Math 继续是唯一签名规则来源。

`explicitFormulaRhs()` 仅提取一个加性、正号、单位系数的 `ddt(target)`；保留其他加减/负号树的求值次序，重建仍用 Math 工厂。它不会把代数关系、约束、`ddt(rho*c)` 或质量矩阵自动转换为显式 ODE，不做线性化或算子恒等变换。现有 conservative/shared/SST provider 保留真实融合形式。

ScalarDiffusionCentral2 对现有串行单 patch scalar 支持原 Central2 加局部 Add/Subtract/Negate/Multiply/Divide/Pow；反应 `-k*c*c` 与 `-k*pow(c,2)` 已经实际执行和独立验证。新增合法局部组合无需增加 ScalarMethod 的 AST 模板分支。A2 初版要求一个 target diffusion；A6 已将普通命名场路径扩展为零个或多个常系数 Central2 occurrence、交叉 Stage 输入和三分量向量。空间需求只施加到实际 stencil 输入；transport/shared/SST 的原融合数学与限定范围保留。其他算子、非线性 implicit 或未实现布局继续准确拒绝，不换算法。

扩展数学身份在每个 EquationRegistry/OperatorCatalog 的编译范围内检查 name/hash/rule 一致性，无进程全局可变 registry。STATE 同名引用在公共编译入口核验签名。现有 core domain 声明未扩充；模块扩展可声明规则和注册现有 catalog 的 numerical implementation，无需修改 Math 枚举、中央物理清单或 PlanExecutor。

Explain 的 normalization、每节点实现和 target-relative 非线性说明来自实际 compiled scalar contract；原数学、occurrence、来源、STATE Stage 版本及融合成员覆盖继续保留。Legacy FormulaCall/FormulaMode/CompiledFormulaCall/FormulaCompiler::compile 已删除；有效矩阵/值核验证迁到现役 provider/GlobalDofSystem，`compileFormulaValue` 与数值内核保留。

这不等于通用非线性隐式求解、通用资源装配、Subcycle 或任意模块组合已经实现。后续 STATE/HOW/资源任务保持独立。

### A3：可扩展派生 STATE 与资源视图（2026-10-21）

`StateSymbol::producer` 引用不可变 `StateProducer`：模块声明数学关系/签名、依赖、支持的 view kinds 与现有计算的 binder。它不是第二个 STATE registry。`compileDerivedViews()` 校验声明/依赖环/签名/版本/输出 slice 并冻结 manifest；只绑定同一请求版本的依赖，不推导新的 HOW 顺序。缺实现或版本、布局不支持时明确拒绝。

Velocity 在 A3 已退出原中央物理 switch；2026-10-22 的核心单路径收口进一步删除整个 `StateDerivation` 机制。builtin contribution 绑定原 `Field::thermodynamicState().velocity`；其原 `(cell,version)` EOS memoization、writer invalidation 与 cache miss 时机保留。普通 view 查询不增加刷新或同步；本轮没有把全部 EOS lazy closure 改成 eager Update。U 为 primary 时仍直接引用原 velocity storage，不从 conserved state 强行派生。homogeneous 的 U 声明在 composition 引用实际 partial-density sum，执行仍调用原 closure。

物化生产者经 `StateRealization::bindDerivedViews()` 绑定原输入/输出数组，经显式 HOW Update callback 调用 `updateView()`。读前核验已有 distributed generation、依赖可读性、输出是否被覆盖、原 Stage epoch 和要求的 halo publication；读取不更新、不通信。revision 句柄在启动阶段冻结；版本快照是有效性 metadata，既有数组、registry 与 owner 不变，没有第二个数值 cache 或时钟。所有 distributed 注册必须先完成，绑定后保持原有句柄稳定的生命周期约束。

测试专用 `q=2a+b` 经同一 contribution/typed registry/compiler/provider/PlanExecutor 运行 Current update，并用原 StageLoop/Stage reader 验证 Stage 和 Current 分离；没有新增 q 分支或 enum。该测试不代表任意物化闭合已经具有 RK stage treatment。Working(U) 继续是原 compiler-owned snapshot，不把 snapshot 伪装成重新派生另一个时间层。

非 Field 的 `ResourceViewRequirement/ResourceView` 是小型 typed borrowed handle，检查 id/type/owner/access、真实 owner lifetime 和 binder 提供的 availability/layout 条件。具体方法可通过原 providerContract 携带需求；Residual 测试直接引用原 PatchWorkspace 对象，不构造假 StateSymbol、不复制数组。不宣称 application 所有资源已自动装配，也不提供动态插件 ABI。

A3 的原始范围和证据见 [A3 报告](../docs/derived-state-resource-views-a3.md)。当前 p/T/h/μ/ν/导热率和 immersed wallDistance 已迁到同一 producer 入口，见下方核心单路径收口；原 EOS/cache 和具体资源 owner 保留。该迁移不等于通用 MPI view 生命周期已经证明。A4 的限定耦合贡献与 HOW 生命周期见下节；Subcycle、通用 Newton/monolithic 与整体资源装配仍后置。

### 核心单路径收口（2026-10-22）

Equation 的 provenance 仅描述来源，不改变数学。`Equation::authored` 和裸 RHS symbol 自动包装 Source 的解释已删除。融合残差 provider 只接受明确的 `Source(...)`；原普通 scalar/赋值中的代数 RHS 保留自己的 AST 语义，不一律改成 source。真实的 signed/coefficient、EOS、时间、布局和边界能力校验继续存在。

所有已登记的派生 STATE 经 `StateSymbol::producer` → `compileDerivedViews` → `StateRealization` 绑定。p/T/h/μ/ν/导热率的 producer 只投影原 `Field::thermodynamicState()`，wallDistance 只读取原 Field geometry。无额外数组、物理 clock、EOS binding 或 cache；绑定不触发 closure。p/T/h 声明 density/momentum/energy 依赖；transport 投影显式命名 Field-owned thermodynamic cache，不为 primitive-pressure case 注册不存在的 rhoU/rhoE。cache 的全 Q 失效仍由原 Field writer/version 负责，不能把这一 owner-bound closure 宣称为任意可替换 STATE 输入的通用热物性实现。

Physical packed 输入需要 componentOffset；数值 provider 已切片的 Stage 输入不再二次偏移。Stage 投影只在原 stage epoch 可读。Working/Correction/Workspace 不自动继承 Physical producer：compiler snapshot、PressureOperators/Corrector 的修正数组继续由原 owner 持有。U/p/T 作为独立主变量时移除 producer，绑定原 writable storage。读取不插入 Update、halo 或 MPI。

架构守卫禁止恢复旧来源解释和派生枚举；solver/system 仍禁止 include concrete model implementation，仅允许上述三个模块数学声明头，并检查其不引入 Field/MPI/solver/backend 实现头。provider 的独立 AST/capability 校验不因数学声明归属调整而删除。

核心单路径阶段未实施 A5；后续首个 A5 闭环已删除 scalar 专用执行分流，flow/Eulerian 的其他具体装配仍保留，见下方 A5 范围。既有 serial/MPI pressure shape 差 `1.0851072147488594e-6` 仍超过 `1e-6` 原阈值，不改容差、不因同配置前后相等而称跨配置一致性通过。实施与本轮验证见 [核心单路径报告](../docs/core-single-path-closure.md)。

### A4：耦合贡献与 HOW 生命周期（2026-10-22）

`PressureCouplingInputs` 是 composition 的方法输入，复用 EquationCall、Target 与既有 constraint 引用；不持有状态、矩阵或 clock。默认方案声明这些引用，显式组装可覆盖它们。Transformation、HOW 与 WHICH 消费同一个已选输入。control 层只选择 typed preset/counts；方法展开保留原 Sequence/Loop/EquationCall/Commit，不创建全局 Coupling driver。

Pressure methods 的 `viewReads` 显式声明 Correction(p) 和适用的 Working(p)，compiler 冻结到同一 NumericalProvider-owned CompiledStateView。读取不创建新数组、不执行刷新。rAU/faceResponse 在原 predictor 产生，correction sequence 读冻结资源；下一 predictor 或 correctedFlux 更新通过 `workspaceInvalidates` 明确结束上一资源版本，不能把旧 correction 带到下一次消费。`commitRequires` 要求显式 correction/publication/flux-consistency 已完成。普通 loop-local workspace 不自动逃逸；已失效资源不能由 scope 恢复或 stepWorkspaces 名单重新复活。

原 PressureOperators 的 baseVelocity_ 仍是 fixed-time U_n snapshot，iterationVelocity_ 是 outer 起点快照；没有新增 clock 或第二份物理 Q。VelocityCorrect 在原位置发布 U，PressureUpdate 更新 p，FluxCorrect/原 correction tail 发布通量或 EOS energy；lazy EOS cache 继续由原 writer 失效。上述校验覆盖声明过的 source HOW、版本与资源，不宣称所有 kernel 读写都有通用运行时证明。OldTime snapshot 的通用借用接口和整体 binder 装配保留到后续。

同配置短时串行/两进程压力结果、首个矩阵/RHS与 trace 保持。串行/MPI 压力形状差仍有迁移前即存在的约 1.085e-6（原阈值 1e-6），没有扩大容差。改名、错误版本、缺 publication、lazy invalidation、scope 负例及独立双方程 Loop 见 [A4 报告](../docs/coupling-contributions-how-lifecycle-a4.md)。没有实现 Subcycle、通用 Newton/BlockSolve、新 SIMPLE/PISO/PIMPLE 数学或 A5。


### A5 首个闭环：所选 lowering 与声明驱动的标量资源（2026-10-22）

公共 `NumericalCompiler` 通过 `IProvider::lowerSpatial` 调用冻结的 WHICH 实现；不再解释融合后端的裸符号、符号/系数或项左右侧限制。`SF_fusedSpatial` 供现有 ConservativeResidual、primitive/conservative pressure 实现共享，保留原 matcher 与 kernel 限制。普通 scalar 继续使用其局部表达式 lowering；没有放宽扩散、时间、layout、边界或 EOS 能力。

`CompiledEquationCall::physicalResources` 复用已有 typed `ResourceViewRequirement`。所选 ordinary scalar 实现声明 immutable `ScalarEquationData` 端口，application 的 `bindScalarResources` 核验 owner/type/occurrence、字段初值/边界和唯一存储，完成全部校验后才向 StateBundle 注册原数组 views。数组与强端口由局部 `ScalarResources` 保活，SolverServices 只借用端口；不分配第二份 Q/clock 或 stage workspace。Explain 区分编译需求与实际绑定日志。

独立多标量和 flow+多标量消费同一 binder；删除 `executeScalarEquations` 及中央 scalar-provider 装配分支。已有流体主变量由原 BuiltinStateCatalog 声明，native added equations 不再把它们重新注册为 NamedDistributed 标量。控制层仅解析已有 metadata，没有调用 kernel。

验收组合是既有显式 flow + 两条交叉反应/扩散标量，共用 RK4；这不是一般资源自动装配、任意方程支持或新模块数值耦合。Eulerian/flow 的其他资源、模型自有数组、MPI/phase/face ports 仍由现有 owner/binder 负责。原 stage、BC/halo、canonical COPY 和 residual SUM 不变。

压力原 serial/MPI shape 阈值失败保留；独立收紧 absolute linear tolerance 的诊断支持其来源为线性求解精度，未更改正式输入或原阈值。旧 level-set isolated-kernel 参考仍未完成首差异归因；当前同配置回归及选定 D1 serial/MPI publication 对照单独验收。完整范围和逐例证据见 [A5 报告](../docs/resource-binding-a5-report.md)。

### 更新后的顺序与文件约束

A1：Math 签名 + Typed WHAT + 至少一个现役 Numerics/WHICH 绑定闭环。
A2：数学规范化、结构/依赖提取与 Numerical capability matching；这不等于重新建立类型分析器。
随后：派生 STATE 扩展 → HOW 高层求解语义 → 通用资源装配；再完善结构化 Explain、内置控制、外部公式及远期变分入口。

现有 Explain 和支持路径持续保护，不因顺序调整移除。Math/WHAT 的概念边界不要求按概念各建文件；复用当前 `methods/math`、`core/system` 及数值模块，优先短名 `SF_signature`、`SF_formula`、`SF_contribution`，检查同名冲突后再命名。不增加 `SF_analysis`/MathAnalyzer 平行类型体系，不为图示创建新顶层目录。本节为设计更新，不授权自动重命名现有文件。

## 历史实施优先级（2026-10-20）

本节保留 2026-10-20 的规划；最新顺序以首页 2026-10-21 Math / Typed WHAT 章节为准。不是新增能力完成声明。**四轴耦合契约 → 完整执行解释 → 内置模型/预制方案控制 → 外部自定义方程**。变分入口仍为远期规划。旧章节、旧路线中的物理模块扩展顺序不覆盖本节；数值支持范围仍以源码与逐例验收为准。

当前旧 Equation DSL、影子方程容器、旧 execution policy/fragment 已退役；705 个 C/C++ 文件、96,886 物理行的统计见 [退役报告](../docs/legacy-code-retirement-report.md)。已有有限 native scalar 输入继续维护，但近期不扩大外部公式语言。详细计划与冻结门槛见 [四轴优先工作计划](../docs/four-axis-first-work-plan.md)；可执行任务说明见 [长任务提示词](../docs/four-axis-first-long-task-prompt.md)。

### 分层与控制方向

```text
GUI / CLI
    ↓
IO → typed user settings
    ↓
内置选择与组装 ← 模型、算法、预制方案库
    ↓
WHAT + STATE + HOW + WHICH
    ↓
Compiler / capability validation
    ↓
Compiled Program ───→ Explain
    ↓
Runtime resource binding → PlanExecutor → kernels / backend

后续用户公式/变分前端 → 同一贡献入口 → 同一四轴与执行链
```

GUI 只编辑设置、展示解释和结果；IO 只负责配置身份、schema 和 typed decoding。内置选择与组装将用户选择展开为已有 SystemContribution、执行作用域及数值绑定，不拥有 timestep、迭代、MPI 或第二份 state。库内的模型、算法、预制方案可以由不同代码实现，不要求新增一个统一“大控制器”类。

运行时只绑定和执行冻结结果，不根据模型名字重新选择数学或算法。Explain 关联同一 authored/resolved/compiled 对象及 provenance，不另维护展示用方程库。配置选择、展开来源、实际执行方法和未参与执行的设置必须区分。

### 四轴核心契约的冻结目标

| 轴 | 唯一职责 | 需要明确的边界 |
|---|---|---|
| WHAT | 方程、项、代数关系、约束及其数学引用 | 不持有 storage、全局循环或隐式模型调度 |
| STATE | 符号身份、主变量/派生量角色、形状、位置、存储/闭合契约 | 不决定方程顺序；保守形式与存储变量分开；workspace 不成为第二份 physical authority |
| HOW | equation occurrence、target、顺序、迭代、stage placement、读版本与提交关系 | 不重定义数学；耦合与时间推进独立；不静默重排用户流程 |
| WHICH | 方程/算子的具体实现、recipe、能力和实现需求 | 不从 formulation 或模型名称推断另一轴选择；不支持就明确拒绝 |

Compiler 验证引用、读写版本、系数生产/消费、资源 ownership、作用域和 provider 能力，并降低为现有 compiled program。耦合未知量的依赖环允许由显式耦合方法处理；派生闭合环、无定义读取、重复 authority 仍应拒绝。一般非线性/隐式系统需要实际装配与求解实现，不因存在 Loop/BlockSolve 类型就宣称支持。

冻结的是经过代表性实现验证的语义、ownership、最小接口和错误行为，不冻结全部 provider 清单、内存 ABI 或所有数学能力。后续新增模型原则上只贡献数学、资源需求及实现/注册；若确需新增核心语义，应记录缺口、兼容影响和新验收后显式修订，不以万能 payload 绕过契约。

### 内置模块与预制方案

模型选择与完整数值方案分开：k–ε 增加 k/ε 方程、状态需求和湍流闭合，不能因此隐式选择 PIMPLE 或 RK4；laminar 不增加湍流输运方程。MRF 贡献项；Ghost/ILW 贡献边界闭合；DLM/KKT 贡献约束及乘子；SIMPLE/PISO/PIMPLE 提供明确的耦合流程及所需派生关系。它们共享接入机制，但数学作用不同。

完整预制方案可以显式组合四轴选择。默认值必须携带来源，显式用户选择优先；两个同级显式贡献冲突应报错，不按加载顺序覆盖。模块关闭后不应留下该模块专属的活动方程、回调或资源需求；共享资源由仍在使用它的贡献决定，不能误删。

### Explain 的完整性要求

每个活动方程和约束均可追溯：选择来源 → WHAT → STATE → HOW occurrence → WHICH → compiled operation → resource/parallel requirements。融合 kernel 可以对应多个数学成员，不强制拆融合；没有对应的数学成员或操作需说明其边界、同步、诊断等职责。

时间解释来自实际 compiled temporal treatment。2026-10-20 已删除只用于摘要的 SolveBlock/TIME_RECIPE：configured recipe/tableau 与实际 temporal root、provider 固定时间处理分别显示。Eulerian 的现有 segregated update 只接受 forwardEuler 意图，不执行通用 RK stages；显式 RK4 输入现在 fail-fast。固定时间 corrector/outer loop 不提交第二份物理时间，没有为修解释更改数值方法。

Explain 完整覆盖当前活动系统；Unsupported 输入应报告阻断链，不制造可执行计划。未来尚未实现的方程不属于“全部已解释/支持”。

### 内部契约 v1 与解释覆盖（2026-10-20，选定能力已验收）

A0/A1/B 已完成限定范围；C1 已在 laminar、单流体 SST、gravity 和两种完整方案的限定范围完成后端闭环；C2 GUI 与 D 扩展仍为计划。详细审计、逐例输入、冻结对照及未支持范围见 [内部四轴契约与 Explain 收口报告](../docs/four-axis-contract-closure-report.md) 和 [C1 报告](../docs/builtin-composition-closure-report.md)。

### C1 内置选择与有来源的默认

`CaseIO → ModelLoader::build → CaseAdapter → inspectCase → SystemBuilder → 原 concrete binder → PlanExecutor`
保持唯一生产链。注册的 `compositionPreset` 只在 ModelLoader 中展开普通 Description 值；
`singleFluidLaminar` 和 `singleFluidSST` 提供明确的四轴默认，不拥有另一条执行路径。
单独选择 turbulence 不强绑时间、EOS 或 coupling；完整方案的默认与模型本身的要求分开。
EOS、物性、网格、初值和边界仍需实际输入，缺失时按既有能力契约拒绝。

显式用户值覆盖适用默认；冲突默认只有在其目标被明确覆盖时才可消解。同级重复显式角色和
未覆盖的歧义默认报来源及目标，不按文件顺序选赢家。原始输入 presence 决定显式性，不能从
NumericsConfig 的 typed 缺省值反推。Description::selectionOrigins 只记来源，不保存第二份数学
或计划，也不参与 runtime dispatch。authored HOW 仍须显式非空 WHICH，方案默认不能补救。

optional turbulence/source block 缺失为 disabled，存在默认 enabled，`enabled:false` 保留输入但
退出贡献。laminar 无额外湍流 PDE。gravity/MRF 的实际活动性由注册模型决定，不再需要
numerics.sources 第二次启用；因此旧输入中“注册 gravity 但没写数值源开关”应明确写
`enabled:false` 才保持旧 inactive 数学。关闭专属消费者不删除用户字段文件，共享 mu 继续由
活动扩散消费者绑定。隐式默认输出只绑定实际 STATE；显式 output 请求没有 provider 时仍报错。

编译后的 STATE required initialization/boundary 约束用于检查直接命名的 provider 输运数组；
重建相状态及专用端口继续由原 binder 校验，不新增中央字段名分支。CLI explain 和启动解释共用
Application::describeCase，展示来源/覆盖、实际活动性及原四轴数学覆盖链。C2 未来复用这些接口，
本轮未实现或验收 GUI；有限 native scalar 和现有 Unsupported 能力范围不变。

- 数学仍只在 `EquationRegistry/FormulaExpr`。`SystemModification` 记录 RHS extension 的 canonical identity 和贡献来源；重复 contribution、同一方程的重复 RHS contribution、缺失 target 明确拒绝。
- `ConstraintDescriptor::equations` 与 `ClosureDescriptor::{equations,states,origin}` 只引用现有数学/STATE。闭合说明不是另一份可执行公式；常量、lazy view 不因此获得虚假 PDE 或刷新操作。PISO 不再注册未参与求解的 fixed-point 关系。
- `equationConsumers/equationOperations` 只读查询原 compiled calls/Plan leaves，供编译器、validator、Explain 和测试共享。可运行系统的每个活动数学成员必须有实际 consumer 和操作；先使用直接引用，只有无直接 occurrence 的约束关系才沿显式 constraint realization 引用追溯。融合成员可以共享 leaves。
- 数组、cache、epoch、OldTime/Stage/Frozen/Lagged、系数 lifetime 和 runtime binding 沿用既有实现。没有新增 STATE、workspace registry、clock 或时间循环；声明检查不冒充对所有 provider 内部读写的运行时证明。
- 代表性验证包含双向标量解析/时间收敛、改名与第三实例、源项、stage/once-step SST、Eulerian RAS、PISO/固定时间约束与 no-slip Ghost/SST 短时保护。不能推导任意 PDE、任意模块组合、production KKT 或通用隐式/Jacobian 能力。

后续修改核心语义须说明消费者、ownership/版本/时序影响，新增有独立参考的正负验收并保护现有冻结集；不建立插件 ABI 或另一个版本分派系统。外部输入、GUI 或新的 preset 只能生成同一贡献与四轴，不拥有另一执行链。


> 2026-10-14 实现更新：生产 IO 仅 SonicFile native registry；CLI/GUI 共享创建、
> reader/writer 和 typed inspection。Level-set RHS/pseudo workspace 与 thermal
> numerical execution 已迁入 solver；mixture flow/thermo、advanced KKT 等仍明确 Unsupported。
> 历史章节中的验收范围不自动升级为本轮证据。详见
> [全模块收口报告](../docs/native-architecture-closure.md) 与
> [Native IO/GUI 契约](../docs/native-io-gui-contract.md)。

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

## 数学入口与统一执行主干（目标设计，2026-10-08）

正式目标是 **双数学入口、统一方程驱动**：用户可直接贡献方程，也可贡献泛函、耗散关系和约束，由变分前端生成方程。两种输入汇入同一 WHAT，并组合独立的 STATE / HOW / WHICH。这里的“变分驱动”指数学生成，不引入第二个全局求解器、时间循环或 physical-state authority。

```text
Mathematical Input / Model Contributions
    ├─ Equation Input: PDE / source / constraint ─────────────┐
    └─ Variational Input: Action / FreeEnergy /              │
                          Dissipation / Constraint          │
                  ↓                                         │
          Variational Compiler                              │
          declared rules + assumptions + boundary terms     │
                  ↓                                         │
          equations / constraints / closure relations ───────┤
                                                            ↓
                                          WHAT / EquationRegistry
                                                            +
                                               STATE / HOW / WHICH
                                                            ↓
                                          Existing System Compiler
                                                            ↓
                                  Compiled Program (equations, numerics,
                                                    solve, runtime)
                                                            ↓
                                  Runtime STATE realization + binding
                                                            ↓
                                              OpRegistry / PlanExecutor
                                                            ↓
                                            numerical kernels / backend
```

图中的 Variational Input / Compiler 是 **Planned**，当前没有承诺通用泛函解析或符号变分能力。Equation Input 已有 C++ contribution；任意 YAML PDE 输入仍受实现限制。runtime binding 主干已完成职责拆分，见本文末尾及 [2026-10-08 报告](../docs/runtime-binding-migration-report.md)。后续阶段、门槛和验收统一维护在 [双数学入口实施路线](../docs/dual-mathematical-input-roadmap.md)。

六个职责层保持独立：Composition → Semantic IR（WHAT/STATE/HOW/WHICH）→ Compiler → Compiled IR → Runtime Binding → Execution。变分前端位于 Composition 到统一方程语义的入口；现有 compiler 负责能力检查、绑定与 lowering。Provider 是 WHICH 选择到具体实现的桥梁，不是第五个数学 authority；具体 binder 理解自己的数学和 workspace，通用 PlanExecutor 只执行已编译控制流。

### 变分数学契约（规划）

| 对象 | 数学职责 | 必须声明的条件 |
|---|---|---|
| Action | 由驻定作用量生成动力学关系 | 变量、允许变分、时空域、端点与边界条件 |
| FreeEnergy | 生成化学势、界面力或其他能量导数 | 积分测度、梯度依赖、边界项、变量约束 |
| Dissipation | 通过选定耗散原理生成耗散关系 | 速率/通量变量、配对、符号和原理；不能一律套用 δS=0 |
| Constraint | 生成约束及乘子耦合关系 | primal/multiplier、作用域、配对、约束方程和可解性要求 |

Constraint 复用已有 WHAT/STATE 语义；泛函引用统一符号身份，不创建独立变量注册表。前端生成的数学通过现有 SystemContribution / EquationRegistry 接入；需要新增符号时贡献 STATE metadata，不分配数组、不替用户选择主未知量、时间 recipe 或 coupling。模型可以同时提供直接方程、泛函、约束和可选数值实现；WHICH 仍显式选择实现。直接与生成贡献必须检查重复 equation/term/constraint，禁止双重施加同一物理作用。

首版仅实现明确列举的解析变分规则，不要求通用符号微分。每条生成方程保留 functional/rule/assumption/boundary provenance；缺少推导规则、边界条件或下游 numerical provider 时明确失败，不以占位公式或隐式替代宣称支持。explain 应能从来源泛函追溯到生成 WHAT，再沿已有六级解释到数值与并行执行。

连续泛函推导后离散，与先构造离散泛函再求导，是两种显式路径；不假定二者交换或具有相同守恒性质。连续弱式必须保留分部积分边界项；当前 FDM provider 不支持的弱式不能直接当作可执行强式。离散规则必须记录 measure、mass/weight、adjoint 和边界处理。

### KKT 的归属与首个离散示例（规划）

对有限维速度 u，选择对称正定质量/度量矩阵 M、冻结 predictor u*、线性约束 D u=0 与 J u=Ub，定义：

```text
L(u,p,lambda) = 1/2 (u-u*)^T M (u-u*) + p^T D u + lambda^T (J u-Ub)

stationarity:
[ M  D^T  J^T ] [ u      ]   [ M u* ]
[ D   0    0  ] [ p      ] = [ 0    ]
[ J   0    0  ] [ lambda ]   [ Ub   ]
```

这是指定离散坐标配对下的约束投影，不是 Hamilton 时间作用量的直接表达，也不是所有 CFD pressure/IBM 的通用矩阵。加权场空间使用相应 adjoint 或已吸收权重的代数矩阵；不能把现有 IBM spreading 无条件替换为裸转置。物理压力与此乘子的符号、单位和 dt 缩放必须由具体推导给出。压力 gauge、约束兼容性/冗余、边界项和 nullspace 必须显式处理；生产隐式动量块也不必等于此示例的 M。

WHAT 持有 momentum/continuity/no-slip 数学，STATE 持有所选 U/p/lambda，HOW 描述 BlockSolve/迭代，WHICH 选择 KKT block 实现、Schur 或增广拉格朗日策略及线性后端。变分前端不强制整体 KKT；SIMPLE/PISO/PIMPLE、显式 forcing、fractional projection 和 Ghost 保留各自适用的数学与执行方法。当前 production KKT 缺口不会因为新增前端而自动消失。

### 模块扩展与代码边界（规划）

Euler 流体作用量、黏性耗散、固体弹性能、界面自由能和相场可各自提供推导规则；热传输需明确热力学变量及能量/熵平衡。VOF 守恒输运、level-set 重初始化和经验湍流闭合仍可直接贡献方程与数值方法，不要求全部变分化；压力也并非在所有 formulation 中都是不可压缩乘子。

建议未来在 `src/solver/variational/` 建立独立数学前端，按实际需要实现 `SF_functional`、`SF_derivative`、`SF_compiler` 和规则注册。模型特定泛函由 models 贡献；纯数学运算可复用 methods。Constraint 使用现有数学类型，不另建不兼容系统。Registry 仅在需要多条可扩展规则时承载规则查找，不能成为第二份 STATE/WHAT authority；不预建空目录、通用 Manager 或额外大体系。

前端依赖 core 数学 IR 和无 Field 的数学规则，不能依赖 runtime binding、MPI/HYPRE、mesh 实例或具体求解器生命周期。数值 workspace 仍属于运行时实现，physical state、clock、phase registry 与 canonical COPY / contribution SUM 语义保持唯一。

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
| coefficient versions / publisher capabilities | physical dependencies only | physical symbols only | scope / occurrences | declares allowed treatments and provided set | freezes unique producer/resource |
| pressure arithmetic implementation | | | | numerical provider | lowers |

禁止依赖：HOW transformer 写 NumericalBinding；Equation 保存 target/order；NumericalCompiler 选择 PISO；StateRealization 选择 topology；Builder 按 state flags 选择 solver family；source ExecutionProgram 持有 compiled workspace；WHAT transformer 创建 runtime OpIds。

## WHAT — Equation

WHAT defines mathematical equations. `SF::System::Equation` 保存稳定的 `id`、左右两侧 `FormulaExpr` AST 和 provenance。`EquationRegistry` 是数学定义 authority；它检查重复注册、未注册 lookup、缺失 replacement 和非法 AST/operator occurrence。Equation 不保存 execution target、order、循环归属或 numerical method。

实际 single-fluid density registry 是 `continuity`、`momentum`、`energy`。常密度 NS registry 是 `momentum`、`continuity`；后者是 `div(U)=0`，由压力一致性方程满足，不额外伪造一次独立 continuity solve。压力耦合注册 `pSimple`、`correctU`、`correctP`、`correctFluxp`；固定时间迭代另有 `relaxIterate`、`restoreFlux`、`checkConvergence` 数学关系。

Builtin、C++ user contribution 与有限 native typed `equationRegistry.add` 经 `SystemCompositionBuilder::addEquation(Equation)` 进入同一 registry。编译器/provider 根据 AST、target 和已实现的 capability 绑定，不能按 builtin/user 来源切换执行路径。native 首版开放具名标量的 ddt/diffusion 和显式算术 AST；结构可解析与 provider 可执行分别校验，未知符号/operator、非法参数、不支持的数学明确拒绝。来源文件、方程与表达式路径进入原 provenance/operator occurrence，不建立第二份数学 IR。任意 PDE、extend/replace/disable、通用符号微分仍未实现。

Formula/Expression AST 是 equation mathematical implementation detail，不是顶层 execution authority。`SF_formula.h/.cpp` 实现 AST、Equation 与 registry；保留此文件名不意味着恢复 Formula solver。FormulaGroup/FormulaGroupRegistry 已移除，守恒向量 fusion 属于 numerical provider 优化。

2026-10-20 退役：Raw/Executable system 不再保存 legacyEquations / legacyDefinitions，
旧 descriptor/Definition 注册重载、role extension、执行 policy、flat→AST conversion 与
legacy plan/provider fallback 已删除。唯一数学 authority 是 EquationRegistry/FormulaExpr；
HOW/WHICH 编译产生 method/temporal 控制流，SolvePlanner 直接冻结它，不推断旧 target。
兼容链退出后，Definition/Expression/TermKind 已无生产消费者，其定义与实现也已退役；
Equation::Compressible::System 与 mixture 的旧 RHS/workspace/commit 是独立现役实现，保持不变。
AssemblyPlan 已没有数值消费者，相关类及实现在此前阶段删除。

## STATE — Base variables and realized views

`core/system/SF_stateRegistry.h` 的 `StateRegistry` 与 `StateSymbol` 是变量 namespace authority。符号保存 shape/components、location、parallel ownership、primary/derived role、physical storage contract、模块 producer 或剩余 legacy derivation、initialization/BC/restart 属性。没有 equation id、target binding、order、Loop、PISO、RK recipe，也不保存 coefficient publication resource 或具体 Stage/Frozen 选择。

原生 case 必须分别选择 WHAT、solution STATE 和 HOW。`EquationCompositionConfig::stateDeclared/solutionVariables` 是输入；冻结后 `StateRegistry::solutionVariables()` 是成员身份 authority。`StateRole` 仍描述数学/存储职责，不决定是否属于 solution；不增加 `Conserved` 互斥角色。物理守恒密度及其 balance 属于 WHAT，一个主未知量可以同时表示守恒密度。

```yaml
# solvers/state.yaml
SonicFile:
  object: solver
  type: state
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

`StateRealization` 绑定 StateBundle 的 physical authority。Generic DirectEvaluation/LinearEquation 的临时数组归 solver execution lifetime；pressure backend 绑定已有 predictor、pressure correction 和 canonical face arrays。RK OldTime 引用 q0，Stage 是当前 published physical state 的时间限定 view，仅在对应 stage 内可读，不创建第二套 Q。U 的模块 producer 与剩余 thermodynamic p/T/h/transport legacy derivations 使用 readonly closure views。普通赋值写入 derived physical target 明确拒绝；物化结果只由冻结的同一 producer Update 写回原输出。

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

Pressure module 有三个独立通道：PressureConstraintTransformer 注册数学 relations（WHAT）；`applyPressureExecution(program, request, inputs)` 只变换 HOW；`pressureNumerics(request, inputs)` 只贡献 WHICH。HOW 函数不接收 NumericalBinding，不知道 PressureMomentum/PressureCorrection provider IDs。

内置方案声明 base equations 及 `PressureCouplingInputs`；SIMPLE/PISO/PIMPLE 方法请求 pressure relations，并用显式 EquationCall/Target 引用消费 predictor、融合 balances 与替换的 Commit。方法不查找 `momentum.default` 或默认 equation ID；自定义 C++ composition 可传入不同 equation/occurrence/relation identity。当前 provider 的 U/rhoU/p layout 与数学能力限制独立校验，不等于支持任意未知量或任意压力数学。

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

transformation provenance 显示实际 predictor equation/输入 occurrence、target 与原 contribution 来源及哪个方法执行变换；其他模块 entries 保留其 origin/order。所有 single-fluid pressure 不再注册独立 PressureSchedule LegacyExecutionPolicy；loop/count authority 只在 HOW。Eulerian shared pressure 也由 native HOW 保存三层循环次数；旧 sharedPressurePlanFragment/policy 已删除。

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

HOW 是唯一 loop/count/order authority。WHICH 的局部 fragment 提供原 assemble/solve/correction kernels；root/Commit 的 pressure.prepare、pressure.step.commit、time.commit 按 selected method lifecycle lower。该流程生产 leaves 的 frozen owner=flow.conservative；无兼容selector标记。bindPressureOps 保留为 callback registration；无主时间或 corrector loop。

现有 conservative provider 仅支持单 patch、单 stage、PISO、无额外非正交 pass 的 laminar single-fluid contract。conservative SIMPLE/PIMPLE fixed-time、RK stage treatment、多 patch/MPI conservative pressure、IBM/turbulence/phase services 尚未实现；编译或启动 capability validation 明确拒绝。保留 generic outer/pressure/nonOrthogonal HOW 不代表这些 kernel 已支持。

## Module composition

`SystemContribution` 统一贡献 equations/AST extensions、unknowns/closures、execution entries 与 numerical defaults。source contribution 不拥有全局 timestep 或 MPI barrier。

| Module | Equation contribution | Execution contribution | Numerical contribution |
|---|---|---|---|
| Navier–Stokes | momentum, continuity；density 另含 energy | ordered default calls | conservative / pressure momentum providers |
| MRF、gravity、wall heat | extends momentum/energy AST | none | existing compiled source kernels |
| Single-fluid SST / k-epsilon | separate native k, omega/epsilon AST；provider-owned mu_t closure | order 1/2 generic EquationCalls before flow calls | TurbulenceTransport；explicit pair fusion -> flow.turbulence |
| Eulerian turbulence | native phase-mass weighted k/epsilon/omega AST、implicit sink 与 mu_t algebraic closure；原数组 ProviderDistributed aliases | 显式 outer closure group → flow groups → transport group，完整 phase/model fusion | EulerianTurbulenceClosure / EulerianTurbulenceTransport → flow.eulerian-turbulence，原矩阵内核 |
| Level set / homogeneous / mixture | level-set 原生 advection、pseudo-time、geometry；homogeneous/mixture 原生声明 | level-set 原生 physical/pseudo-time/Commit；其余按实际能力 | level-set 原有 kernel；mixture flow/thermo 仍 Unsupported，不能将声明当执行 |
| PISO | pSimple, correctU, correctP, correctFluxp；conservative 另有 publishPressure | consumes momentum call; inserts pressure Loop | constant / conservative relation providers |
| SIMPLE | pressure and fixed-time relations | outer Loop with single correction | constant density: original PressureOperators；conservative fixed-time Unsupported |
| PIMPLE | pressure and fixed-time relations | nested outer/pressure Loop | constant density: original PressureOperators；conservative fixed-time Unsupported |
| Eulerian–Eulerian | namespaced native Formula equations | native PIMPLE EquationCalls / Loops / Commit | AST-validated compiled contracts / specialized kernels |
| Ghost IBM | native boundary read/write/stage/order contract；不增加虚假 equation | every spatial evaluation 的原 BC/halo/ILW 链 | explicit ibm.boundary port；原 kernel |
| forcing / DLM / KKT | native impulse/work、lagged force、projection 或 coupled constraint AST；original STATE aliases | explicit complete post-predictor group / block | immutable native Immersed method contract；原 kernel；KKT production capability 仍有缺口 |

源 HOW 不再包含 legacyEntries；模块只贡献正式 occurrences。没有生产者的
LegacyExecutionPolicy/LegacyPlanFragment、其 lowering 及 provider 推断已退役。
未知/缺失的编译调用、provider 或能力仍明确失败，不能回到旧执行器。
SolveBlock/TIME_RECIPE 及其重复校验、打印已删除。时间意图、实际 temporal treatment 和 Plan root 分别来自原有冻结对象，不保存第二份摘要计划。

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

`SystemBuilder::build` 收集 modules 的独立 contributions。`composeContributions` 做 module applicability matching；native pressure 匹配已贡献的 mathematical/default occurrence，不读取 StateRealization flags。composition 完成后保存 source NumericalSelection，再调用 `NumericalCompiler::compileSystem`。

`compileExecutionProgram` 后校验每个显式 solution variable 是否由 method 的 `writes` 或 writable resource 覆盖；保守 predictor 的 fusion 同时写 rho/rhoU/rhoE，不能只看外层 rhoU target。未使用的 solution T 等明确失败。

Compiler 通过 const source HOW/WHICH 和 const equation/provider interfaces 接收输入。scope sorting 发生在 compiler 的独立 HOW 副本，source 不被修改；provider 不得改变 equation/occurrence/target semantics。编译器不增添 equation、不选择 preset、不替换 numerical selection。compiled operation/resource metadata 是输出，不是另一份数学 authority。

`SF_methodObjects.cpp` 只做 deterministic occurrence binding、capability validation、generic lowering；domain numerical implementations 位于 `SF_builtinProviders.cpp`。`SF_numericalCompiler.cpp` 按 provider 声明的 spatial terms 和 kernel/workspace requirements 绑定，不按 momentum/U/constantDensity 识别 solver。空间输入由所选 provider 的 spatialTerms/sourceMathInputs 显式声明，generic compiler 不导入旧 numerical adapter。不新建平行 SimulationPlan、manager 或 forwarding planner。

`SolvePlanNode.provider` 是每个 executable leaf 的最终 implementation identity。`RuntimeRequirements.operationBindings` 只是该 frozen leaf ownership 的校验/index 投影；`SystemValidator` 拒绝两者不一致。`OpRegistry::bind(id, provider, callback)` 注册实现的 owner，`PlanExecutor` 在任何 operation 执行之前检查所有 ID 与 owner，然后直接 invoke；它不读取 state/formulation 或 coupling preset。SingleFluidStepper 按 frozen provider ID 构造/绑定资源，不再由 constantDensity flag 选择 backend。

生产 leaf 必须在编译时拥有明确 provider/operation。IBM、Eulerian pressure/turbulence、single-fluid pressure 均由所选方法冻结绑定，OpRegistry/PlanExecutor 只消费冻结结果。旧兼容 selector 与其空标记已删除；不存在 compile 后的 provider matching。

`CompiledSolvePlan.root` 是冻结后的 executable operations，不再是独立 authored HOW；source program、compiled calls、numerical bindings 与 explain 都来自同一编译链。compiler 不猜 equation target/order，也不选择 coupling preset。

`sonicSolver explain` 首先分别显示 WHAT / EQUATIONS、STATE / BASE VARIABLES、STATE VIEWS、HOW / EXECUTION、WHICH / NUMERICS 的 source choices；之后显示 COMPILED OCCURRENCES / STORAGE / PROVIDERS 和 COMPILED SOLVE PLAN，明确区分 source 和 lowered outputs。Legacy execution authority 与 pending entries 明确展示。unsupported 的真实边界必须保留，不能用解释文字冒充运行实现。

## Native single-fluid turbulence and mixed temporal lowering

显式 primary flow STATE selection 不要求用户把每个模型辅助变量加入 `solvers/state.yaml`。用户显式选择 kOmegaSST 后，模型在同一个 StateRegistry 贡献 k/omega；选择 kEpsilon 则贡献 k/epsilon。用户 primary flow selection 仍是 rho/rhoU/rhoE，U/p/T 是 flow closure views，mu_t 是模型 derived closure output。没有第二个 STATE registry。

RAS 的 k/omega/epsilon/mu_t backing 仍由 Turbulence::Manager 内的 ScalarFields 持有；`ProviderDistributed` 描述这一事实。Manager::registerDistributed 发布原数组的非 owning views，StateRealization 和 Physical(target) 直接引用它们，compiled view owner 为 NumericalProvider。没有为架构另分配数组；ScalarFields 原先分配的四个槽位仍保留。当前没有 turbulence restart reader，因此这些 STATE 不宣称 restart capability。

single-fluid native WHAT 如实描述现有 production RAS 内核：ddt(k/omega/epsilon)，显式 production/dissipation、除以 rho 的固定局部扩散，以及 SST cross diffusion。旧注册形式包含实际未执行的 div；本阶段删除这个数学描述偏差，没有给内核增加对流或隐式汇。kEpsilon 扩散仅 mu_t/sigma；SST 为 laminarMu + sigma*mu_t。既有 positivity floors、production limits、wall-distance 逻辑和 in-place sweep 完整保留，不代表数值缺陷已解决。Eulerian 的质量加权输运和隐式汇是另一个旧内核，未与此实现合并。

source HOW 仍是通用 EquationCalls：order 1 k->k、order 2 omega/epsilon->自身，随后 flow 的 continuity(10)、momentum(20)、energy(60)，最后 terminal Commit。WHICH 为两个 RAS occurrences 显式选择 TurbulenceTransport，并声明完整 pair 输入。编译出的 fusionKey/fusionMembers 合约要求同一 Sequence 中两个相邻 sibling、正确目标、正确 model/state/math、无缺失或重复。两份 WHAT 和两份 source HOW 保留，只有数值 backend operation 合成一次 turbulence.advance，owner=flow.turbulence；共享 provider ID 本身不会触发这种 fusion。

通用 temporal compiler 支持最小 mixed pattern：flat independent once-step prefix -> one contiguous fused temporal group -> terminal Commit。ITemporalMethod 接收编译后的 physicalStepPrefix；explicit recipe 在原 prepare/dt 后、begin/snapshot 前安放它。因此实际顺序为 flow.step.prepare -> flow.dt.compute -> turbulence.advance -> flow.step.begin -> StageLoop(flow RK) -> flow.step.commit -> time.commit。temporal -> independent -> temporal、嵌套 mixed scopes、缺失 terminal Commit，以及尚未指定 temporal placement 的额外 root lifecycle decorations 均明确 Unsupported，不猜测、重排或丢弃。此 prefix 是降阶输出，没有新增 source scheduler 或第五模块。

TurbulenceTransport 保留当前一次/物理步的显式原位更新，不采用 flow RK4 stages。它消费 StateBundle::dt；不计算第二个 dt 或推进 clock。FlowStepBegin 已移除 correctTransportModel 调度；CompiledSolvePlan 中的 turbulence leaf 绑定旧 correction callback。其局部 BC -> WriteOwned -> ReadHalo -> correct -> BC -> WriteOwned -> ReadHalo 顺序保持不变，MPI 实现仍由 execution runtime/backend 承担。application 只构造 Manager 和参数，不再通过 transportedTurbulence predicate 选择执行拓扑。

DNS 不贡献输运方程；Smagorinsky 仅贡献 mu_t 代数闭合及 TurbulenceClosure refresh，保持原缓存刷新位置，不产生 k/omega/epsilon 输运。全局 term recipes 只由声明 usesSpatialRecipes 的 provider 消费；RAS 的固定局部 Central2 不会误触发 flow diffusion recipe 默认选择。

RAS 当前执行能力仅 single-fluid、serial、single patch、density explicit flow；pressure、MPI/multi-patch、Eulerian 和 implicit RAS 未放行；IBM 仅开放后述显式静止黏性 Ghost/SST 组合。旧 single-fluid pressure service guard 保留；Eulerian 的 EeTurbulenceSolve 已是独立原生 provider callback；Eulerian core scheduler 与 assembly authority 均已原生化，AssemblyPlan 已删除。详细审计、诊断对照的基线来源及 18 项自审见 [single-fluid turbulence migration report](../docs/single-fluid-turbulence-migration.md)。

## State, workspace and parallel invariants

Physical state/clock/field composition 的唯一 authority 仍是 StateBundle。预测量、p'、flux、RK k 值和 linear/constraint buffers 属于 execution/provider lifetime，不塞回 Field，也不创建新 physical registry。

Runtime interfaces 隔离 MPI；raw MPI 仍仅在 infrastructure/backend。state、geometry、metrics、labels：owner -> COPY -> replicas。每个 physical face 只有 canonical F*，COPY 后 owner +F* / neighbour -F*。residual/source/load/constraint contributions：local -> SUM -> canonical owner。不得平均 state/geometry/face identity，或把 SUM 改为 COPY。

本轮 architecture migration 不改 WENO/TENO、flux split、pressure/Rhie–Chow arithmetic、RK coefficients、CFL、HYPRE tolerance、IBM mathematics 或 distributed ownership。编译和 explain/static checks 不代表数值基线验证。

本阶段先完整运行原始 Sod 与 cylinderFlowGhost，再在 provider authority/active STATE 修改后运行同一份输入。详细 frozen evidence、数值验证状态和未覆盖项见 [runtime provider authority migration report](../docs/runtime-provider-authority-migration.md)。原始 cylinderFlowGhost 是 density/RK4/Ghost，不能作为 pressure Working/Correction 或 MPI/turbulence CFD coverage 的证据。

## Native Eulerian shared-pressure execution

Eulerian PhaseSystem 明确贡献每相 phaseMass/momentum/enthalpy、shared p 和派生 alpha/rho/U/h/T，用户不必重复列为单流体 primary STATE。ProviderDistributed views alias 原 PhaseSystem/PhaseState 数组；Correction(p) alias 原 pressureCorrection，canonical face arrays、diagonal、previous state 保持 numerical workspace。

`SF_eulerianRelations` 描述冻结数学：参考相 mass 由 volume closure 重建；动量压力项为 alpha*grad(p)，扩散为 RHS；shared pressure 包含 sumGas(alpha/p)/dt reaction、sum faceAverage(alpha²/A) diffusion、volume source/divergence/history；energy 包含 material pressure work。旧 phantom/不准确描述已纠正，原算术未升级。

明确注册的 PIMPLE 贡献 generic Sequence/EquationCall/Loop/Commit：outer 内 continuity group -> momentum group -> pressure/nonOrthogonal 下 shared-pressure correction、phase correction、paired phase flux correction -> enthalpy group，最后 Commit。三层次数只在 source HOW；当前其他 Eulerian coupling 明确 Unsupported，不自动替换。source 无 raw Ee operations。

六个 Eulerian methods 显式选择 flow.eulerian-pressure，完整 phase inputs/group targets 按声明的 PhaseSystem storage order 验证。通用 compiler 支持 N-member fragment fusion，保留每相数学 provenance，core leaves 使用冻结的 provider/operation。Provider-local prepare/source/BC/halo/canonical-copy lowering 保留原 Ee 微操作顺序。`sharedPressurePlanFragment`、S_EE_PIMPLE policy 与 transformer Ee inventory 已删除；EulerianStepper 不从 policyKind 选择算法，只执行 frozen plan。

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

所选 provider 的支持签名只用于严格 AST/backing/phase/model 能力校验，不注册第二份 WHAT；system 不 include concrete model implementation；具体 provider 只消费模块的数学声明。runtime 只校验冻结实现绑定与实际模型数组，不重新解析方程名。generic compiler 不识别 turbulence，explain 展示每个 occurrence 的模型、slot、group、owner 及数值方法。

保留完整生命周期：initial prepare/commit；dt probe 中 prepare/sourceTimeStep；每 outer 在预测器前 prepare/同步、energy 后 solve、solve 内 BC/validate/reprepare；outer final BC/同步；成功物理 step commit(m*phi)，随后 diagnostics、clock commit。这些 prepare 不可去重。

当前能力为 serial single block、既有 Upwind/diffusion/HYPRE backend；MPI/multi-patch Eulerian turbulence 在 capability validation 明确 Unsupported。实测 2-rank host 输入因缺少分布式 phase transport provider fail-fast，是 implementation capability gap，不是 sandbox launch failure。已有 MPI pressure regressions 保留。

数值基线来源、独立矩阵 harness、输入/输出 hashes、phase subset/reorder、多层迭代与旧内核对照见 [Eulerian turbulence native migration](../docs/eulerian-turbulence-native-migration.md)。原 production turbulence 入口无法运行，报告不把恢复接线后的诊断执行冒称历史 production baseline。全局 compatibility DSL 仍为 level-set/mixture 等 consumer 保留；IBM 的 native contribution 接线见下一节。

## 同步多实例与多 provider 时间组（2026-10-09）

独立方程可以共同参与一个 compiled temporal group。Equation occurrence、numerical provider type、融合 numerical owner、temporal participant 是不同身份：两个 scalar 共用 provider 类型仍有独立操作和 workspace；NS 的 continuity/momentum/energy 仍合成一个 `flow.conservative` participant。`fusionKey` 只描述原数值融合，不承担跨 provider 的时间组合。

`CompiledExecutionProgram::temporalParticipants` 冻结 identity、calls、targets、provider 和阶段操作；scalar 的 OpId 由操作类型与 occurrence 长度/内容组成，OpRegistry 的重复绑定检查不取消。WHAT/STATE/HOW/WHICH 仍是输入 authority：HOW 决定成员顺序，WHICH 必须提供 split-phase temporal capability。参与者集合不自动推断物理算法、不重排 RHS，也不把 Current/OldTime 静默替换为 Stage。

```text
Application::Execution::executeEquations
  -> existing Time::Driver / SingleFluidStepper
  -> realizeState
  -> bindSingleFluidOperations
       FlowOperations::temporalParticipant (only when required)
       scalarParticipants (one numerical owner per scalar occurrence)
       FlowOperations::bindStageSources (after all Stage views exist)
       Time::bindStageGroup (bind callbacks in compiled participant order)
  -> unchanged Run::PlanExecutor / OpRegistry
```

共同计划只有一个 StageLoop 和一次 clock commit：

```text
participant preparation
one dt = min(driver limit, config maxDeltaT, participant proposals)
all old-time snapshots
StageLoop(recipe.stageCount):
  Open(epoch/index/time)
  all prepareStage (including each provider's boundary readiness)
  Ready: Stage reads available
  all evaluateRHS (no Stage writes)
  RhsReady: Stage reads closed, advance available
  all advanceStage (existing recipe mathematics)
  Close
validate all final publications
publish all
time.commit
```

`validateTemporalPlan()` 在 capability report 和 runtime binding 校验冻结 leaf/provider 顺序、stageCount、全部成员阶段及唯一 Commit；缺失、错位、重复 leaf 不是 Runnable。runtime epoch 检查 snapshot、index、recipe、RHS readiness、advance 和 publication。StageReader 在绑定时解析 storage/component view，不在 cell loop 查找字符串；只在共同 Reading phase 可读。复用 buffer 每次读表示当前 epoch，不是持久旧 stage snapshot。

scalar physical arrays 由 application/STATE 持有，ScalarOperations 只持有 old/stage/RHS/increments。新非拥有式 scalarInstances 表按 occurrence 绑定 immutable BC/source data，并独立检查 target；weak_ptr 用于检测数据已过期。旧单 scalar raw pointer 仅兼容一个 occurrence，仍要求调用者保证生命周期。binding 表和 reader 不复制别的 provider 的 physical arrays 或 workspace。

普通 native case 现可显式注册 STATE、有限 WHAT、HOW 与 WHICH，运行多个独立常系数扩散标量。`CaseAdapter` 将 typed add lower 到原 SystemContribution/FormulaExpr；fields registry 的声明提供自定义符号的 metadata，已注册 internalField/boundaries 提供初值和边界。每个方程有独立 equation identity、每次调用有独立 occurrence 和 target；compiler 决定 OldTime/Stage/RHS demands，不由用户配置数组。

没有流体 preset 请求时，通用 composition 不要求 Momentum、Continuity 或 EOS；这些检查仍由真正的流体 preset 保留。inspection 从冻结的 flow provider requirements 得到 mesh 是否需要 conservative storage。geometry-only loader 不分配 Q、不执行 fluid initialization；声明的 scalar default 覆盖整个几何，不套用 fluid 的 set-only IC 覆盖推断。`adapters/SF_scalarResources` 只为编译端口需求创建唯一 NamedDistributed ScalarField 数组及 immutable boundary ports；geometry-only 和 flow 入口共同消费它，随后由已有 `executeEquations`/统一 Plan 执行。它不拥有新的时间循环，不依据字段名或数量选择算法。flow/Eulerian 的既有资源装配保留，并未宣称所有 provider binder 已统一。

独立 scalar 的正式输出沿用 registeredStateFields/ResultWriter/fields registry；没有伪 rho/rhoU/rhoE 或 FluidStateModel。当前 native 闭环限制为 serial 单 patch、均匀轴对齐网格、完整面 constant fixedValue/zeroGradient、ScalarDiffusionCentral2、显式 forwardEuler/classicalRK4。native source callback、通用初始化函数、restart 文件持久化、scalar MPI、新模块 stage 组合不因此开放。`check/explain` 校验声明/编译能力，完整 BC 面覆盖和实际 layout 在 runtime binding 校验；尚未加载 mesh 时不宣称已验证全部 runtime geometry。

flow 保留原 `Time::Explicit` workspace、stage 系数、RHS 和 Q publication。`evaluateRHS()` 与 `executeStage()` 拆出可选预装配屏障，但旧 standalone fused-stage 路径仍有效。混合组中 Q 是原实现发布的 time-qualified Stage storage；所有 RHS 完成后才允许 advance。未新增永久 Q，未承诺任意运行异常的全状态事务回滚。scalar 最终 publication 前先校验所有成员，已能确定的 binding/plan/数据生命周期错误不提交半份 scalar 或 clock。

现已实现的 WHICH：

- `ScalarDiffusionCentral2`：多个均匀正交网格 scalar；一个常系数 diffusion，加 known source/有限常量/有限常系数乘共同 Stage scalar。真实双向耦合允许代数读环，因为它们都读取冻结的同一 Stage，不能被当作顺序更新。
- `ScalarTransportUpwindCentral2`：`ddt(rhoC)+div(rhoU,rhoC/rho)=diffusion(rho*D,rhoC/rho)+source(Sc)`；读取 Stage rho/rhoU/rhoC，真实一阶 upwind convection、面平均 rho*D diffusion。此选择使用中心面 mass flux，未借用高阶 NS canonical mass flux；不能宣称任意可变密度下浓度恒定性与高阶 flow 严格兼容。
- 显式 `linearScalarForce` / `linearScalarWork` source：compiled Stage source kernel 经稳定 reader 进入原 ConservativeRHS；动量力密度为 kappa*rhoC，能量功使用同一 Stage 速度。不是 application 的 step-end 修正。

共同组当前支持 Euler/classicalRK4、serial single patch、无额外 auxiliary stage treatment 的原 conservative NS。纯 scalar 不安装 NS、Q 或 EOS。MPI/multi-patch、scalar implicit/SSPRK3、任意 nonlinear coupling、超出有限 typed native 范围的 AST/YAML、turbulence/IBM/interface 的新共同 stage treatment 均 Unsupported。已有其他模型自身原执行路径不改为此组，也不失去原有能力。新 provider 需提供真实阶段实现及 binder，不以空 callback 宣称支持。

通用 PlanExecutor/Time::Driver 不认识 scalar、feedback 或 flow family；compiler 与 execution 不新增模型专属 runner、第二套 tableau/clock/physical registry。kernel 不调用 MPI，原 canonical COPY、GlobalDof SUM 及 FluxField/Residual solver ownership 不变。原多实例数学与范围见 [共同 stage 迁移报告](../docs/multi-equation-stage-migration-report.md)；正式 native 入口、独立 C++ 参考、错误配置与生产保护见 [native 方程执行闭环报告](../docs/native-equation-execution-closure.md)。

## Current implementation and verification status

实现状态与数值验证范围分别记录；已有 executable path 或结构测试不能代替该路径的 CFD baseline。

| Path / contract | 当前源码状态 | 已有数值验证 |
|---|---|---|
| Single-fluid conservative explicit | WHAT/STATE/HOW/WHICH 编译为 owned fused stage；runtime 执行冻结 provider | Sod 454 steps / 5 VTS，Ghost cylinder 11450 steps / 51 VTS；前后 time/dt、全部输出 hash、积分、极值、范数、closure diagnostics 完全一致 |
| Constant-density PISO/SIMPLE/PIMPLE | native HOW + method-owned pressure/Rhie–Chow operations；无 legacy provider routing | 本阶段 PISO、16²/32² Poiseuille、fixed-time SIMPLE/PIMPLE、既有 MPI 回归前后全部通过；同配置前后 VTS hash 一致 |
| Active STATE / lazy views | catalog metadata 与 active registry 分离；HOW/WHICH 请求 view，runtime 绑定 storage | inactive symbol、qualified target、storage/owner contract 检查通过；本阶段 CFD 实际覆盖 Working(U)、Correction(p)、Physical(U/p/phi) |
| Single-fluid conservative/variable-density pressure | native six occurrences + explicit three-equation predictor fusion；bindPressureOps 仅注册已有 kernels | pressureConstraintPiso 1 step，time/dt=1.659315e-06；输入、诊断、操作顺序、字段统计与两个 VTS 均前后完全一致；fixed-time conservative 等 capability 仍 Unsupported |
| Single-fluid transported RAS | native WHAT/STATE/HOW/WHICH；完整 pair fusion；Plan 每物理步调度 flow.turbulence，solver 执行原位 sweep/BC/diffusion，model 保留物理闭合与原数组 | 2026-10-16：非均匀 kEpsilon/SST、no-slip Ghost/SST、LES、stage SST、Eulerian RAS 和 Sod 共 8 对短时 CFD 与迁移前一致；不代表长期湍流物理验收，见 once-step ownership 报告 |
| Eulerian turbulence | 原生 mass-weighted WHAT、单 registry STATE aliases、显式 closure/transport HOW 和独立 WHICH；物理数组/闭合属于 model，old-time、节点 BC、装配/线性求解属于 solver；无 model commit | Phase C：10 个短时案例各 6 steps / 7 VTS；矩阵、解、物理输出与迁移前当前工作区逐字节一致；单 outer RAS 180→198 leaves，多 outer 870→924（只展开原后处理）；serial single-block，MPI/multipatch Unsupported。详见 execution ownership 报告 |
| Eulerian shared pressure | native phase/pressure WHAT/STATE/HOW/WHICH；explicit N-member fragment fusion；EulerianStepper 仅绑定/执行；AST -> immutable compiled assembly contract，零 native flat DSL 依赖 | 现有双相/MRF/PIMPLE：6 steps、7 VTS、168 Ee operations 和诊断前后逐字节一致；三阶段 CFD 均完全一致；最终完整回归见本次迁移报告 |
| IBM / distributed execution | native boundary/impulse/projection/block contracts；原 numerical kernels 和 COPY/SUM semantics | 各方法分别验证；不能用 Ghost 代表 forcing/KKT；本次完整时域结果见 IBM migration report |

## Next migration sequence

截至 2026-10-08，四契约、跨模块能力/状态版本检查、静止黏性 Ghost/SST 受限组合和 runtime binding 职责迁移已有实现；这些不等于完整 SST、所有 IBM、多相或分布式组合已通过物理验收。历史数值证据继续由各迁移报告保存，不改写为新能力证明。

当前顺序为：先完成 turbulence、IBM、level-set、Eulerian–Eulerian 的选定能力与数值验证，再启动变分前端。模块完成以明确 capability matrix 与冻结基线为门槛，不以所有可能组合均支持为前提。随后依次推进有限规则前端、pressure/IBM 离散约束应用和其他变分物理模块。详细任务、依赖和验收见 [双数学入口实施路线](../docs/dual-mathematical-input-roadmap.md)，本节不维护第二份独立待办顺序。

## Native IBM contributions

IBM 使用现有四契约与同一 compiler。Ghost/ILW 贡献 boundary closure 的 read/write/stage/order contract，physical BC → halo COPY → reconstruction → publication → halo COPY → operator，不伪造 EquationCall。

BP 是 post-predictor momentum penalty 与 midpoint work，没有乘子；Peskin / explicit DFM 消费旧 lagged force 后计算下一层乘子，不宣称当前步精确 no-slip；FTS/fractional DLM 贡献原 mass-norm projection。显式/fractional自推进使用 virtual-fluid generalized mass projection，并区分 active translation/rotation；未改为实体刚体质量时间增量。

原 multiplier/laggedForce/solid velocity 通过 SpecializedExecutor STATE aliases 和 original-storage port 绑定；surface lambda 仍是原数值内核的 transient workspace。native complete mathematical group 经 adjacent fusion lower 一次原 project/KKT callback，所有数学 members 保留 provenance/write-set。compiled immutable contract 校验 AST/backing/target/algorithm/rigid DOFs；native correction 当前只支持完整 implemented predictor 之后、commit 之前，不支持任意 stage treatment。

IBM 的 legacy Definition/HOW/policy/provider routing、空 ImmersedConstraintTransformer、C_IBM 前缀 matching 与独立 descriptor inventory 已删除；descriptor 只保留实际 method options/KKT capability 与 display summary。generic temporal compiler 不按 IBM 字符串选择 execution。

原 ILW/J/Jᵀ、body/surface projection、rigid solve、KKT/HYPRE numerical helpers 保留。三种 KKT 有 structural native contracts，当前 canonical production inputs 新旧均因 recipe/lowering 缺口 fail-fast；不能把它们标为已验证 Runnable。Eulerian＋IBM、除后述静止黏性 Ghost/SST 外的 RAS/LES＋IBM、pressure/implicit＋IBM 的缺失 mathematical/stage/storage contracts 在编译/capability 阶段明确拒绝。Peskin lagged surface MPI 允许有限非负的 local partial diagonal（包括空支撑的零），Runtime SUM 后严格要求正性；单独验证合法空 rank、canonical lambda COPY 后移除该编译限制。force/mask 是 Eulerian owner 已完成的诊断量，输出 replica 使用 Runtime COPY，不是 SUM 或 average。物理 ILW 的法向和模板遵循 source mesh 的 communication mask；ILW 在重建前声明 conservative halo read，之后保持 physical BC → halo → immersed closure 的相对顺序。IBM 的 mass/spreading/load 使用统一的物理 nodal volume，内部 partition endpoint 不引入半体积，EMPTY 端点共同表示真实厚度。上一冻结版本 47/47 host CTest 与首差异定位已经完成，7 个串行冻结案例逐字节一致，8 个 MPI 案例全部运行至 endTime=5 且逐帧 output replica 冲突为零；完整时域差异、剩余 kernel 一阶矩/力矩性质与 RAS boundary 缺口见 [IBM parallel consistency validation](../docs/ibm-parallel-consistency-validation.md)。历史迁移 preservation 结果见 [IBM native contributions migration](../docs/ibm-native-contributions-migration.md)。

## Surface transfer 的一阶矩契约

`partitionOfUnity` 与 `linearReproducing` 是 WHICH 的两个显式 numerical selections；WHAT、STATE、HOW、原 J/Jᵀ/Schur/projection lifecycle 共用。新选择在原 geometry support 上由 owner raw moments 经 Runtime SUM 形成一个完整 affine Gram matrix，解出同一组 J/spreading 权重。typed selection 被冻结在 CompiledImmersedContract，runtime 不得另选 kernel；未知策略或不满秩支撑直接失败，不能 clip、正则化、扩大支撑或降阶 fallback。旧选择不增加 moment communication。

物理坐标 x/X 与同一个力矩中心、marker area m、Eulerian physical nodal volume V 定义载荷与功率。常量重现和 `Σ w(x-X)=0` 一起保证 marker/Eulerian force、torque 一致；带权 adjoint 保证 power 一致。Peskin 的 consumed lagged force 与 next force 必须分别使用对应时间层。新权重可有符号，当前仅支持完整 3D affine basis（EMPTY 双平面可满足；coplanar 单层支撑拒绝）；不由力矩性质推导单调性或通用 KKT 可运行性。[专项报告](../docs/ibm-first-moment-report.md)分别记录数值实现、算子指标与完整 CFD 验证状态。

## Native interface 与剩余多相能力

Level-set WHAT 分离 `ddt(phi)+U·grad(phi)=0`、frozen-reference pseudo-time 与 normal/curvature algebraic relations；CSF 的 force/energy 数学成员绑定原融合 RHS 写集合，Ghost-fluid 仅声明 jump query。STATE alias 原模型数组。HOW 显式冻结 phi0、伪时间 Loop、geometry publication、physical commit。WHICH 在 reference preparation 固定伪时间 recipe，stage 消费该 workspace，不拥有第二份 dt/time；界面 port 不拥有循环。

D1 的实际发布链为：物理 BC → phi owner 写失效声明 → HJ 所需 phi halo → 第二次物理 BC/必要时重新发布 → normal 计算 → normal owner COPY/一层 halo → curvature、mask、物性与 band 更新 → curvature owner COPY/一层 halo → owner-only 物理体积积分 SUM → 消费/输出。BC 的 zeroGradient donor 使用真实 set 的法向轴，不能从局部 patch endpoint 推断物理壁面。共享 physical replica 可参与同一局部几何计算，但 authoritative normal/curvature 始终由 owner 发布；几何计算不得在 normal COPY 后再清零 normal。material cache 与 normal/curvature publication 分开标记有效，失效状态不得被 CSF、jump 或输出消费。

初始化输出在 StateBundle attach 与上述发布之后；physical commit 使用同一准备链，writer 不执行隐式通信或副本平均。每 rank 一个 patch 的 MPI 也必须归并 diagnostics，不能将 single-field storage layout 当作串行。原 HJ-WENO、RK、CSF、Young–Laplace、canonical face COPY 与 GlobalDof SUM 数学保持。倾斜平面、圆界面、CSF 与不同 x/y/2×2 分区及实际多 source-zone MPI 的短时证据见 [D1 修复报告](../docs/level-set-mpi-publication-repair.md)。MPI 重初始化和 massCorrection 仍 Unsupported，串行多 source-zone solve 当前也 Unsupported；sharp indicator 对 phi=0 的舍入敏感性仍保留，不由 D1 扩展数学能力。D2 的伪时间 workspace/执行与 D3 的独立 thermal 执行 ownership 已在 2026-10-14 收口（见下文）；mixture flow/thermo binding 仍 Unsupported。

Homogeneous 按 partial-density component layout 声明原生 balances；mixture 推进 phaseMass，alpha 保持派生。它们没有内置 flat Definition/HOW producer，但无 generic-EOS flux / native mixture closure binding 时不能报告 Runnable。共享 flat DSL 保留于真实 compatibility consumers；死 prefix/forwarding API 已删除。

`executionProgram` 输入完整 typed source HOW，`providerBindings` 提供明确 WHICH；整份 authored HOW 不继承默认 binding，unknown/conflict/unused binding 与 unsupported topology 必须明确失败。`targetKind: Workspace` 与 working `*`、correction `'` 不混淆。native preset 与 authored 已知契约经同一 compiler/runtime，回归比较 operation sequence 和完整字段。任意 mathematical YAML/custom storage capability 不由该入口自动推导。详见[中文迁移报告](../docs/native-remaining-migration-report.md)。

## Scope、版本与跨模块能力验证

HOW 是唯一调度 authority。显式 ExecutionProgram 按作者的 sibling 数组顺序编译；默认 recipe 主干保持既定顺序，模块使用 PlacementRequirement 在其指定 scope 内声明 before/after，composer 仅插入这些已声明节点。多个数值含义不同的可行顺序没有明确关系时要求显式 HOW。依赖不能自行选择 coupling，也不能展平 Loop/StageLoop 后跨 snapshot/publication/Commit/synchronization 排序。RAS/LES、IBM、level-set contribution 不再拥有全局绝对 order。完整 ordered fusion group 仍是部分现有 provider 的实现能力。

STATE 描述真实 storageBinding/storageKey/location/component slice 与派生依赖；StateUse 声明 Current/Stage/OldTime/Frozen/Lagged，StateEffect 声明写入、cache invalidation、materialized publication 与 halo availability。owner current 不代表 halo current。U/p/T/h 的 Lazy view 依赖原版本失效机制，不引入 refresh operator；normal/curvature 等物化数组需要 publication。Eulerian coefficients 每 outer 冻结，Peskin force 合法滞后，level-set phi0 在整个 pseudo-time scope 冻结。普通 loop workspace 不泄漏，真实 timestep 数组由 provider 声明 stepWorkspaces。派生缓存环与 coupled unknown block 分开检查。

这些是 provider 明确声明的数据流编译验证；未声明的 specialized accesses 不自动获得完整运行时 freshness 证明。实际 Field/model/storage/cache ownership 保持，验证不复制数值数组，不执行刷新或通信。

WHICH provider 声明 capabilities 与 conditional requirements；通用 resolver 报告真实缺失 contract，不写模块配对黑名单。IBM descriptor 从 composition 冻结到 executable/compiled contract，provider 支持表及 AST 校验保持独立；runtime 只校验冻结 fluid-port 与实际 storage，不再根据 AST 猜 variant。BP/Brinkman 保持 penalty/work，不伪装 no-slip constraint。除后述已绑定静止黏性 Ghost/SST 外的 RAS＋IBM、phase IBM ports、production KKT 等缺口继续准确 Unsupported。数值内核与 canonical COPY/GlobalDof SUM 不变；实现与验收见[跨模块组合报告](../docs/generic-cross-module-composition.md)。


## 静止黏性 immersed wall 的数值契约（2026-10-06）

`ImmersedWallClosure` 是一次解析后冻结的 numerical selection：默认 EulerSlip 数学不变；StationaryNoSlipAdiabatic 是独立的约束二次 Ghost 重建，不能作为 Euler ILW 的自动降阶路径。它要求静止单流体、serial single patch、PerfectGas 五分量 conservative state、Central2 viscous gradients，Euler ILW 显式关闭。U 在局部壁面切平面为零，rho/p 零法向导数给出绝热 T；full-rank donor support 不成立或重建非物理时 fail-fast。曲壁是局部平面近似，不宣称 generic viscous ILW。

IBM geometry lifetime 保存 donor indices/weights；physical Q 和 clocks 不复制。壁面上的数值零距离点属于边界，不作为 SST fluid-cell 计算点。`immersed.wallDistance` 是 Field 原数组的只读 derived STATE view，不触发刷新、不复制几何。已有 turbulence arrays 仍由 Manager 持有，`IImmersedTurbulenceBoundary` 只写原 k/omega/mu_t ghost arrays，没有 solve/schedule authority。

SST port 的 k=0、omega=60 nu/(beta1 d1²)、mu_t=0 以真实第一流体支撑距离为依据。signed polynomial ghost 不是独立 physical unknown；原 positivity floors 只在 fluid 点执行，新 Ghost 值不被 floors 改写。既有 source/diffusion/in-place transport 数学、once-step 更新和 flow RK topology 保留。新组合读取 prepare 后的 physical velocity，发布 step-frozen mu_t，flow stage 只读取缓存；原未绑定壁面时的 dynamic-viscosity 查询方式保留。该原内核没有 convection，不能由此次 boundary capability 推断完整 SST RANS 已实现。

WHICH 的 bound capability 使用 `(name, equation, fluidPort, boundary)` 精确匹配。带对象要求不能由无绑定字符串或另一个 wall/port/equation 满足。当前 wall=`immersed.ghost` 对应现有唯一 IB 实例，不是多壁面 registry，也不新增调度层。provider 仍独立验证 AST、原数组及方法支持范围。

Plan 顺序继续为 prepare(boundary/halo/Ghost) -> dt -> once-step RAS -> flow snapshot -> RK stage boundary/operator -> commit/time commit；边界闭合不变成 WHAT 的伪输运方程。默认和 authored HOW 使用相同 lowering/provider。负例覆盖 old slip、混用 Euler ILW、MPI、顺序错误，以及 capability 存在但对象错误。局部解析 profile、实际黏性算子网格收敛与短时 CFD 是本轮证据；长期充分发展湍流、移动壁、MPI 新壁面和三维 CFD 验收没有被这些测试替代。详见 [中文报告](../docs/viscous-immersed-wall-sst-native-composition.md)。


## Runtime binding 的稳定职责边界（2026-10-08）

编译后由同一个 `RuntimeRequirements` 和 frozen `leaf.provider` 驱动 implementation binding。`SingleFluidStepper` 仅保持 bundle/compiled products、generic STATE realization、稳定 `OpRegistry`、`PlanExecutor` 与唯一 physical clock commit；不解析 equationMethod/providerContract，不持有模块数值 workspace，不在每 advance 重绑 callbacks。

`bindSingleFluidOperations` 是一次性 composition entry，不是 scheduler。原 flow/RK/pressure workspace 由 `FlowOperations` 持有，transport、level-set、immersed 使用窄 bind functions 和原 ports。具体 binder 可以理解自己的数学实现、参数与 backing；通用 host/executor 不按模型名选择执行。原 boundary/publication/finalize/geometry/clock 的操作位置保持独立。

callbacks 借用稳定 bundle/services/compiled objects，强引用实际 numerical owner；动态 dt limit 从 host 的稳定输入读取。phi0/temporal view 仍 alias 原 workspace，physical state 不复制。host 不可 copy/move，patch identity 在绑定之后不可替换。compiler 与 binder 共享纯 compiled-stage 校验；geometry/storage/service validation 不伪装为静态证明。

设计与生命周期表见 [runtime binding architecture](../docs/runtime-binding-architecture.md)，具体验证范围与剩余能力见 [migration report](../docs/runtime-binding-migration-report.md)。这次收口不实现 SST convection、level-set MPI 修复、production KKT 或新组合能力。

## 普通物理壁面条件与重构的独立选择（2026-10-10）

速度 `noSlip` 声明静止壁面 U=0；速度 `slip` 声明静止平面壁面的 U_n=0、dU_t/dn=0。这两个名字是边界定律的便捷输入，分别展开成已有 FixedValue(0) / Symmetry 中立约束，由同一 Linear 或显式 ILW3/5/7/9 重构执行。它们不创建 solver、时间方法、额外状态或特殊 runtime path；相同展开输入必须产生相同数值执行与结果。

rho、p、T 仍独立声明各自的边界条件。速度壁面不自动选择绝热/固定温度/热流，也不修改压力条件、EOS、黏性、coupling、时间或 IBM 算法。slip/noSlip 用于标量/温度或 initial condition，以及给静止壁面附非零速度，必须显式拒绝；移动壁面应使用已实现且明确的其他数学契约。通用数学表达与某个 numerical provider 的实现能力不同：例如当前 constant-density pressure provider 支持展开的 noSlip，尚未实现 slip 的法向/切向装配，继续准确报 Unsupported。

多 patch 的全局 thermal set 可以有零本地成员，该 patch 不进行法向推断或写入；不存在的 set 身份仍报错。pressure→thermal 和有交点的 ordered thermal patches 之间必须失效 EOS 缓存，让消费者读取已发布的原守恒数组，不创建第二份温度 authority。physical boundary→halo→IBM→halo→operator 的调度顺序与 canonical COPY / residual SUM 不变。

这条契约服务普通结构化物理边界，不将 IBM 曲面 Euler ILW 声称为通用黏性壁面 ILW。现有曲面壁面选择及支持范围继续由其独立数学契约校验。实际测试设置与证据见 [slip/noSlip 解耦报告](../docs/slip-noslip-boundary-composition.md)。

## 热物性选择、参数与执行绑定（2026-10-10）

EOS、caloric closure 和 transport 属于 `models/`，与 WHAT 中的 Energy equation、STATE 的能量变量、HOW 的时间推进、WHICH 的空间离散分别选择。`hConst` 不增加 Energy equation，也不推进 T；`const` 不选择 CENTRAL2。用户仍须显式提供方程与离散。rhoConst 的 U/p provider 不需要 PerfectGas 能量模型；geometry-only scalar 不需要 EOS。

```text
models/thermoDynamics.yaml
  selection + required parameters + provenance
        ↓ strict native validation / once-only historical normalization
shared_ptr<const ThermophysicalContract>
        ├─ MeshRuntimeConfig / Init (before conservative energy initialization)
        ├─ FluidStateModel factory (single and participating patches)
        ├─ boundary / ILW / IBM initialization
        ├─ molecular transport / RAS coefficients
        └─ CompiledNumericalSystem → explain / numerical consumers / output
```

这份冻结值是 case 热物性的唯一可编辑输入所产生的 authority。SolverConfig、BoundaryConfig 和 CompiledNumericalSystem 共享同一个 const 对象；NumericsConfig、mesh、IBM 和 turbulence 中的标量是 kernel 所需的只读投影。进入 compiler/runtime 时校验投影不能漂移，native boundary/output 在缺少所需 FluidStateModel 时 fail-fast。Field/StateBundle 的 physical state、clock 与 PatchWorkspace ownership 不变。

当前 native 单流体支持 perfectGas + hConst + const，gamma/R 显式必填；选择 const 时 mu/Pr 显式必填。cp/cv 由 gamma/R 推导，可选输入仅用于一致性校验。热导率为 `mu * cp / Pr`。rhoConst 显式提供 rho；未选择 caloric 时不要求 cp/cv/gamma/R。JANAF、Sutherland、未知模型和不支持的 EOS/STATE/energy 组合在 kernel 前拒绝，不能改用 PerfectGas。当前黏性热通量核要求正 conductivity；mu=0 可用于无 diffusion 的系统，不能伪装成该核已支持零导热。

旧 `numerics.transport.enable=false` 的显式停用语义仍保留并显示来源；声明了备用 diffusion scheme 不等于强行启用该项。不含旧开关的 native 输入由实际 diffusion 选择派生 enable 投影。

旧配置仍可声明 `numerics.transport.mu/Pr`；解析一次后冻结，gamma=1.4、R=287.05、mu=0、Pr=0.72 的历史默认来源会显示在 explain 中。native 与 legacy 同时声明同一参数时必须相等，并保留两个来源；冲突报错。历史短参数 diffusion overload 仅用于显式低层兼容调用，native production 使用完整参数或已绑定模型。

标准 native 文件示例：

```yaml
SonicFile:
  object: models
  type: thermoDynamics
thermoDynamics:
  equationOfState: perfectGas
  thermo: hConst
  transport: const
properties:
  equationOfState:
    gamma: 1.32
    R: 310
  transport:
    mu: 2.1e-5
    Pr: 0.81
```

`properties: {}` 不是有效的 native 热物性声明；需要显式补齐参数。CLI init 生成完整的默认声明，而不是让 reader 猜测参数。Thermophysical / thermophysical 的同义 selection 入口使用相同校验；历史 temperature-only contribution 仍独立存在，不等同于 EOS selection。

### Canonical case IO

|语义|目录|SonicFile.object|
|---|---|---|
|STATE / WHAT / algorithm registry|solvers/|solver|
|HOW/WHICH/runtime、旧 flow equationSystem、普通边界 ILW|solvers/|solver|
|EOS/caloric/transport、IBM、turbulence、物理 source|models/|models|
|field declaration、初值、边界|fields/|fields|
|mesh declaration|mesh/|mesh|

温度等 derived view 可以只声明边界而不声明独立初值；configuration decoding 不得因缺少 T 初值而丢弃 thermal boundary。p 初始化和 T 初始化应共享同一套显式边界。

## Phase B：共享守恒通量与 single-fluid stage SST

共享质量通量是原 NS `PatchWorkspace::convectiveFlux` 中密度分量的只读借用，不是另一个 flux provider 的计算结果。`StageFluxPublication` 只保存 patch、FluidStateModel、component、layout、numerical recipe、physical step/stage/time 和 generation 元数据；`StageMassFluxView` 不持有数值数组。RHS 清理、resize、move 和 stage advance 使旧 view 失效。候选面通量完成现有 canonical finalize/COPY 后才能发布，同 epoch 重复发布拒绝。面按正计算方向编号，读取 cofactor 加权 `F_hat`；散度是 `Jac*(F_right-F_left)`，不重复乘面积或除 dx。canonical COPY 与 residual SUM 的实现不变。

`ScalarTransportSharedFluxCentral2` 通过 `equation.scalar-shared-flux` 消费同 stage canonicalMassFlux；浓度采用 Upwind1，扩散采用面平均 `rho*D` 的 Central2。WHAT 为 `ddt(rhoC)+div(rhoU,C)=diffusion(rho*D,C)+source`，STATE 的原 rhoC 是唯一 physical authority。D 可以是常数或明确的 Stage STATE。独立几何标量和旧中心 rhoU 构造仍是不同的显式数值选择，不作为新 shared method 的 fallback。第一版共享标量限串行单 patch、均匀正交网格；MPI/new multipatch 支持不能由 NS 的并行能力推断。

模型显式选择 `transportMethod: stageConservative` 时，single-fluid kOmegaSST 贡献 `ddt(rho*k)` / `ddt(rho*omega)`、真实 convection、rho 加权 diffusion、production、sink 和 cross diffusion，WHICH 为 `SSTStageTransport`。primitive physical k/omega/mu_t 继续 alias 原 Manager 数组；solver 标量 participant 持有 rho*k/rho*omega 的 OldTime/Stage/RHS/RK workspace。模型只计算系数，不通过 `correct()` 推进该路径；solver 持有七个 SST 系数 workspace。mu_t 发布到原数组，并被真实 NS 黏性算子读取。不存在第二份 physical k/omega、第二个 clock 或 SST runner。

具体 coefficient 数学符号仍在 STATE/WHAT：mu_t 和 `SST.nu.k/omega`、`SST.production.k/omega`、`SST.sink.k/omega`、`SST.cross`。
其 physical alias、抽象 TransientWorkspace 属性和闭合依赖不决定更新方式。SSTStageTransport 显式声明完整八成员 Stage 发布集合、prepare/publish 操作及其数值读取需求；once-step RAS/LES 显式声明既有 step-Frozen mu_t 发布。
ConservativeResidual 声明可消费的已实现版本，compiler 只匹配本次 WHICH 已选择的唯一 publisher，不根据基础 STATE 属性选择 Stage/Frozen。

`CompiledEquationCall::coefficientUses` 保存 provider 接受版本与 compiler 冻结的 version、producer identity 和 resource；同一绑定派生 stateUses 与 temporal participant。
这些是既有编译 IR 的数值输入/输出契约，不是第二份 STATE 或资源注册表。缺失、冲突、不兼容版本/halo、缺少完整输出集合和错误 participant 镜像明确拒绝。
StageTransportClosure 在绑定阶段独立核对实际八成员输出，分配原七个 workspace，并将 mu_t Stage view 绑定原物理数组；不按模型名或扫描第一个同名资源重选 publisher。

公共 temporal 编译契约具有 optional coefficient resources：全部 participant boundary/prepare -> Stage ready -> 明确 coefficient publication -> 全部 ordered RHS -> RHS ready -> 全部 advance。RHS 不写其他 participant 的 Stage 状态。NS canonicalMassFlux producer 必须位于 shared scalar consumer 之前，coefficients 必须在 flow/SST RHS 之前；compiler 校验 authored 顺序，不自行重排。全部 physical publish 后刷新最终系数，再执行唯一 `time.commit`。recipe coefficients、stage abscissae、CFL、NS flux、residual 和 MPI 语义不改变。

StageGroup 对冻结绑定和 callbacks 做一致性检查；每 stage 清除 coefficient-ready，RHS 前要求所有已声明发布完成。
下一 stage / 下一 physical step 不能沿用上一 epoch；全部 physical publication 后的最终 mu_t 刷新位置不变。初始化仍用原物理闭合。
`explain` 可显示 producer、resource、version、八成员集合与 flow/k/omega consumer；原 HOW occurrence/provider/provenance 保持可追溯。
契约迁移、同 STATE 的版本能力测试、真实 epoch 负例与 Euler/RK4 精确数值证据见[专项报告](../docs/stage-coefficient-contract-ownership-migration.md)。

stage SST 第一版仅支持 single-fluid kOmegaSST、serial/single patch、forwardEuler/classicalRK4、明确完整面 fixedValue/zeroGradient 边界；继承 mesh halo 标签不等于允许内部/局部点集冒充完整面。新 stage SST＋IBM、pressure、MPI/multipatch、Eulerian、kEpsilon 或未实现的 stage recipe 明确 Unsupported。已有 `legacyOnceStep`（省略选择时保持当前兼容默认）有真实 Ghost/SST、kEpsilon 和旧 single-fluid 消费者；它仍使用前述 once-step 数学，不得称为共享守恒 SST。Eulerian 保留自己的相质量加权矩阵内核。两种 SST 数学不会静默互换。

## Phase C：Eulerian turbulence execution ownership（2026-10-12）

模型 `Turbulence::EquationSystem` 保留原 k/epsilon/omega、mu_t、wallDistance 和 diffusivity/explicitSource/implicitSink 缓存；模型贡献是 WHAT authority，provider 中的 AST signature 只验证已实现数学。物理数组仍由唯一 STATE alias 暴露，不复制为 solver physical state。模型不应用通用节点 BC，不拥有线性装配/求解或 numerical time commit。

`PhaseSolverWorkspace::turbulenceOldTime` 每个被输运 unknown 分配一份 `phaseMass_n * phi_n`；条目同时绑定实际 phaseIndex 和原 physical array。选中相的顺序不能被当作 phase storage slot。初次 attach 在原 initialization BC/prepare 后初始化 snapshot；所有 outer 都读取同一 snapshot，唯一 compiled `ee.step.commit` 后更新下一物理时间层。`StateBundle.time/step` 仍是物理 clock，workspace 的 step/outer epoch 只是消费有效性标记。失败的装配/求解不更新 old-time；不承诺物理数组的整步回滚。

节点 scalar BC 移入 `solver/boundary/SF_scalarBoundary`：保留固定值的镜像 ghost、zeroGradient/symmetry 的边界相邻点复制和全部 ghost 常值延拓、EMPTY 的原 stencil。它不同于普通 scalar provider 的 Neumann 模板，不能因 API 名称相近替换。`PhaseEquationAssembler` 在已有 transport 实现中组装和求解相质量加权矩阵，拥有每 unknown 的可复用 SolverSession，workspace 复用一个 diagonal；不创建第二个 runner 或 outer loop。

编译 fragment 显式显示原尾部：`ee.turbulence.solve`（求解并写回原数组）→ `ee.turbulence.boundary` → `ee.turbulence.validate` → `ee.turbulence.refresh` → 原 final boundary/outer validation → physical commit。predictor 前 prepare/publication、transport 消费 Frozen 系数、solve 后 prepare、dt probe 的 prepare/sourceTimeStep 均保留原频率。系数 epoch 必须匹配当前 physical step/outer，不能以读取“最新值”为由改成 fully coupled turbulence。

本 provider 明确限定 closure gradient=Central2、diffusion=Central2、canonical mass-flux convection=Upwind1；显式选择其他 order 或未知参数 fail-fast。原 CENTRAL2 闭合 kernel 暂留 model 内，公开为 fixed implementation capability，不声称已支持任意梯度。Eulerian RAS/LES 仍仅串行单块；MPI 在 capability 阶段拒绝，多 patch 在单 patch state binding 拒绝。COPY/SUM、canonical face orientation、隐式 sink、线性容差和 PIMPLE 三层循环不变。数值证据和逐例设置见 [Phase C 中文报告](../docs/eulerian-turbulence-execution-ownership-migration.md)。

验收与范围见 [共享守恒输运迁移报告](../docs/shared-conservative-flux-transport-migration.md)。离散守恒预算使用原散度控制体权重，输出统计使用物理梯形积分，两者须分别解释。新增 SST 对流及 common-stage 推进是数学变化，不能以旧 once-step SST 的冻结结果作为新模型的正确性证明。

CaseIO writer 在写入任何文件前计算全部目标路径，拒绝 safeFileName 或大小写碰撞、保留 registry 名、重复 category/name、重复 registry、非法父目录和目标目录冲突。solver/model 的同名对象位于不同目录，不覆盖彼此。Reader 只按当前 object registry 的 reference 解析，不搜索猜测文件；旧分散布局须显式转换，生产 reader 不再读取历史固定路径。重复 referenced file（含规范化路径或 symlink 指向同一文件）和 object/type 不匹配会拒绝。短 registry selection 写出后保持，浮点参数按 max_digits10 写出，避免读回改变系数。


### Native IO 注册身份与单目录约束（2026-10-11）

语义 WHAT / STATE / HOW / WHICH 独立不要求创建独立 case 目录。所有 solver 配置使用 `SonicFile.object: solver`，集中在 `solvers/`，由唯一 `solvers/solvers.yaml`（type: registry）注册。models/fields/mesh 同样保持自己的 object。

双重核验区分 **文件身份** 与 **内部角色**：

|Registry 键 / 子文件 type|注册项内部 type|子文件|职责|
|---|---|---|---|
|equations|equationSystem|solvers/equations.yaml|方程系统 preset/实例|
|state|stateRegistry|solvers/state.yaml|solution STATE 选择|
|equationRegistry|equationRegistry|solvers/equationRegistry.yaml|WHAT 方程选择|
|execution|algorithmRegistry|solvers/execution.yaml|HOW 算法选择|
|how|executionProgram|solvers/how.yaml|显式 HOW（按需注册）|
|which|providerBindings|solvers/which.yaml|显式 WHICH（按需注册）|
|internalField|internalField|fields/internalField.yaml|初值（fields registry 注册）|
|boundaries|boundary|fields/boundaries.yaml|边界（fields registry 注册）|

注册键 `state` 对应子文件 `type: state`；`stateRegistry` 仅作为母文件条目内部角色。母子 object、键与子 type 逐字核验；保留键必须匹配指定角色，单例角色不能换名创建第二 authority。自定义实例（例如 models 的 surface: {type: geometry}）保持 surface 文件身份与 geometry 角色独立，不用几何角色覆盖实例身份。

默认 `flow` 键已改为 `equations`。原有 equationRegistry 选择表使用同名独立键，不与方程系统共用键。HOW/WHICH、runtime、numerics、algorithm 和 ILW 都在 solvers/；热物性、IBM、湍流、物理源属于 models/。

`case.yaml` 是 case 身份/版本元数据，固定的四个 object registry 是引导入口。其子文件只能由显式注册读取。已登记的内置 compact selection 可以解析到登记过的默认文件位置；这不授权扫描目录。fields 的内联字段声明仍是 STATE IO metadata，不是额外的 document registry；初值/边界 child reference 显式注册，并在字段声明加载后分派。

加载链：registry -> parseEntry（保留键/角色校验）-> resolveFile -> object 目录/重复引用校验 -> 子文件 object/type 校验 -> 内部角色 handler -> typed Description。数值 compiler 仍消费原内部角色及 payload；文件身份不成为算法 dispatch authority。

拒绝注册的缺失文件、跨 object/嵌套目录引用、symlink 越界、重复 canonical 文件和第二份 registry。未注册文件即使格式错误也不读取。删除 state、equationRegistry 或 execution 注册后，不能通过同名文件存在恢复；后续能力校验照常报告必需配置缺失。旧式头部/分散布局须显式转换，生产 reader 不保留自动发现 fallback。

CLI init、CaseIO writer 与正式算例遵守同一布局。Writer 先检查整个输出集和注册契约，再写任何文件。物理 state、clock、numerical workspace、provider、stage/boundary/halo order 与 MPI COPY/SUM 均不因 IO 迁移改变。验证证据见 [注册身份 IO 收口报告](../docs/registry-identity-io-closure.md)。

## Native IO / D2 / Thermal ownership closure（2026-10-14）

2026-10-15 source启用语义补充：phaseChange block存在时默认启用，显式`enabled:false`
不得因模型名称再次启用；inactive RPI不要求wall-boiling服务，也不产生沸腾质量源。独立启用的
普通wallHeat仍提供热贡献。dt probe和各outer允许重算源率，但必须先清空重装配；physical state
只由编译数值操作推进，物理clock每步只commit一次。选定Lee/RPI的source-only预算和短时CFD
证据见[执行与预算报告](../docs/source-execution-budget-closure.md)，不代表任意热力学组合已验证。

### 生产输入与 GUI

`Application::createNativeCase → CaseIO → ModelLoader/inspectCase` 是 CLI 与 GUI 唯一创建/校验链。
模型与 solver 配置按 object registry 显式注册；GUI 原生文本编辑先 staging validation，失败
rollback，optional 模块默认 inactive。不存在 OpenFOAM dictionary reader/editor/backend 或
vtkOpenFOAMReader。application 的 native 数值转换 helper 不是旧格式 filesystem parser。
历史输入只能显式转换，不能因文件存在自动加载。

mesh registry 可选注册 `generation: {type: structuredMesh, file: mesh/generation.yaml}`，
child header 为 `object: mesh, type: generation`。CaseIO 验证typed参数，CaseAdapter 冻结
MeshParameters；MeshGen 保留原trilinear/grading/TFI算法。face classification 分别使用每个显式
表面region，不能把相对面的point-set union当成一个表面。生成SFM须可被正式loader重新读取。
GUI四角挤出草图消费同一参数；更复杂编辑UI明确Unsupported，保留真实几何kernel。

GUI用typed native parallel配置启动QProcess；single-step preview 是 test/ 下独立副本的
`sonicSolver steps 1`，正常结束后读取结果。按step文件编号与PVD physical time分开：
PVD标签来自committed StateBundle.time，不拥有时间状态。

### Level-set

phi、normal、curvature仍由原LevelSetModel持有。`solver/algorithm/interface/SF_levelSetOperations`
持有RHS、Reinit::Workspace与frozen phi0，`SF_advection`和`SF_reinitialization`负责数值执行。
原模型不再有primaryScalarRHS/advance/reinitialize转发入口。geometry、CSF、jump数学留在models。

HOW/Plan决定physical RK、pseudo Loop、geometry publication与唯一physical Commit。
D1的BC/invalidation/halo/donor/normal COPY/curvature COPY/diagnostics SUM链不变。
37个单source配置和6个MPI多source配置before/after result/checkpoint exact。
MPI/multipatch reinit仍因缺真实scheduler而Unsupported；不隐式降级。

### 独立温度与 mixture

thermal贡献普通WHAT/STATE/HOW/WHICH；`thermal.T` alias原MultiPhaseModel温度数组，
不能与EOS T/rhoE混成两个writable authority。存在EOS T时选择独立fieldName，mean-flow T
输入/BC与该target分开。`ThermalEulerTransport`仅声明已实现的Forward Euler、serial单patch。
convection要求真实五变量PerfectGas velocity closure；unrelated U alias、primitive pressure port、
未知time treatment、MPI均拒绝。pure diffusion不要求convection port。old/stage/RHS属于solver，
通过现有StateRealization绑定。模型不再有advanceTemperature或beginStep完成Euler的隐藏路径。

mixture仍由phaseMass主量推进，alpha/T派生；真实native flow/thermo及publication binding未实现，
保持Unsupported，不能根据RPI/相物性自动选pressure/density family。Auto/family配置与
print-only conservativeLayout镜像已删除。独立温度不开放为mixture的第二份energy积分。

### 清理与能力边界

native 模块不新增 legacy Definition/role/flat plan；AssemblyPlan 已删除，仍有真实消费者的全局兼容 DSL 按模块保留。
single-fluid legacyOnceStep RAS 的原位 sweep、dt 写回、标量扩散和通用 BC 已迁入 solver；
模型仅提供局部只读物理闭合、原物理数组、floors/壁面关系和代数 mu_t，见下一节。

LES无配置rho clamp已改为finite-positive fail-fast。Scalar Project仍是显式bounding选择；
本轮不删除合法limiter/regularization。TransferLedger是源项预算authority；homogeneous total-energy
mapping只计external energy一次，不重复加phase wall allocations；不支持的phase mechanical conversion
在写source之前拒绝。generic-EOS homogeneous flux、mixture binding、JANAF/Sutherland、primitive energy、
Subcycle、density KKT/advanced rigid augmentation/restart仍按真实实现能力报告。

本轮最终验证为45个serial C++ tests、11个host MPI/HYPRE tests、72组短时保护对照与GUI实际运行，
不替代历史完整CFD或长期物理验收。逐例配置、终点、hash与限制见收口报告和test/证据。


## Single-fluid once-step RAS execution ownership（2026-10-16）

`TurbulenceTransport → turbulence.advance → bindTransportOperations` 仍是原 HOW/WHICH 入口。
solver/algorithm/transport 的 SF_onceStep 执行 k/j/i 原位 sweep、两个原标量数组的扩散与 dt 写回；
每个 cell 调用模型只读闭合计算 production/dissipation、SST blending/cross diffusion。
`OnceStepTransportData` 只借用原 ScalarFields 的数组、BC 配置和物理系数，不分配或复制 physical state。
通用标量 fixedValue/zeroGradient/symmetry/empty 实现在 solver/boundary；保留旧单流体模板的
角点方向与 ghost 写入顺序，不换成语义不同的 Eulerian 节点 BC。

实际叶操作内部顺序保持：BC → Runtime finalize → Runtime prepare → 原位更新 → BC →
mu_t 代数刷新 → BC → finalize → prepare。新 BC 表示通用模板后再执行原 floors/immersed wall 数学。
初始化由 application 在原数组创建/IC 后立即请求 solver BC/closure；没有新增 timestep 或 Commit。
模型不再暴露 correct(field,dt) 或数值推进转发；LES/DNS 用不接受 dt 的 refreshClosure。
绑定数值 advance 时缺少原两标量端口明确失败，不能以 closure-only 刷新冒充 transport。

本实现仍是原 legacyOnceStep source/diffusion 数学，未增加 convection，未转成 rho*k/rho*omega
stageConservative。mu_t 在 sweep 内冻结；无壁面时原 dynamicViscosity 查询计算、Ghost/SST 的
step-frozen 缓存、stage SST 的 publishStageViscosity 以及 Eulerian 矩阵路径分别保留。
SST F1/F2/cross/strain 等物理闭合所用既有 Central2 kernel 仍属 model 的真实消费者。
不开放 pressure/MPI RAS、移动壁、其他 IBM 组合或新时间策略。
数值证据、API 删除与完整逐例设置见[迁移报告](../docs/single-fluid-ras-execution-ownership-migration.md)。


## A5 后续：统一贡献及模型资源契约（2026-10-23）

延用SystemContribution，不新增模块控制器或另一Registry。ContributionRecord.module表示实现类别，id表示实例；默认单实例id=module，显式实例为module@instance。这些元数据仅用于冲突与Explain，runtime不据此调度。

所有方程/STATE注册完成后、transformation之前，builder.finalizeMathematics按原贡献顺序合并显式加性RHS。不得按相同AST去重，亦不自动合并系数。families是作者显式提供的目标引用/namespace选择，不用于猜数学。closure.outputs显式声明闭合输出并拒绝冲突；states仍是依赖引用。缺输出的历史描述不能冒充一个已声明owner。

numerics默认来自module/preset；用户显式WHICH优先，精确重复的显式key拒绝，未覆盖的同级default冲突拒绝。key为equation+occurrence，equation级显式覆盖其默认项，细粒度选择沿原occurrence优先级。authored HOW仍要求显式非空WHICH，不偷偷继承默认。NumericalBinding.origin保留选择/覆盖来源，不参与数学或provider匹配。

CompiledEquationCall.physicalResources继续复用ResourceViewRequirement。普通scalar与single-fluid模型输运分别请求ScalarEquationData/ITransportModel端口；application具体binder按需求保留原数组owner，validatePhysicalResources逐项核验真实view的身份/类型/权限/lifetime，拒绝遗漏、重复及未消费端口。Manager不因configured但inactive而分配数组；已有startup初始化、BC、once-step source/diffusion、stage SST和publication顺序不变。资源校验不插入刷新、同步或时间操作。

本轮没有重建k–ε：其现有single-fluid once-step仅实现冻结的source/diffusion、无对流；支持范围和现有floors保留。native k–ε+gravity、禁用/覆盖、flow+双标量及SST真实运行完成。新增内部双标量模块经同一贡献/绑定/Plan执行，并以独立解析解验证；不等于任意native新公式或任意资源已实现。

无生产写入的SolvePlanNode旧兼容标记及其fingerprint/打印已删除。LegacyPressureInputs、DistributedRowMap、native decode仍有真实消费者，保留到相应模块替换。详见 [本轮报告](../docs/contribution-a5-report.md)；原压力阈值失败及历史Level-set参考失配不因本轮通过而消失。


## A6：通用显式局部方程与多 occurrence Central2（2026-10-23）

普通 `ScalarDiffusionCentral2` 是现役 WHICH 方法标识，现在支持 NamedDistributed 的标量及三分量向量；名称不表示向量被转换成速度 U。没有新增 solver family、provider-family dispatch、时间循环或 Plan 节点。纯反应、多个 diffusion、交叉 diffusion 和向量 diffusion 均使用原 `FormulaExpr → 所选方法 → FormulaOperatorProvider → 资源 binder → 同步 temporal group → PlanExecutor`。

编译的 `diffusions` 每项保存原 occurrence（无显式 occurrence 时使用 AST path）、输入、常系数和分量数；实际 kernel 各自捕获自己的系数与 Stage slot。相同表达式不去重，交叉输入不替换成 target。不再有一个全局 diffusivity 代表所有普通扩散。signed 常系数按作者数学计算，不自动判断扩散矩阵稳定性或修改系数；已有 conservative/shared/SST 的 rho-weighted 数学与系数限制不变。

局部式没有 stencil 时不要求 BC、halo、正交 spacing 或 active axis。被其他 diffusion 消费的局部方程仍必须提供其自身真实 Stage 边界。当前 native SFM reader 保留 NG>0 及非退化几何要求，这是网格入口范围，不是反应算子的 halo 需求；内部零 ghost layer 运行另行验收。

向量以三个原 ScalarField 分量组成唯一多分量 physical view；旧 scalar 和 vector 使用同一 binder。OldTime/Stage/RHS/RK 数组采用 component-major，仅属 solver workspace。Math 声明的标量系数在逐分量值求解中广播；标量 STATE 读取固定分量 0，向量 STATE/邻域读取实际分量。fixedValue 必须给全三个分量，zeroGradient 各分量独立。registeredStateFields/原 writer 输出 NumberOfComponents=3，不造第二份输出数组。绑定失败在任何 distributed 注册之前退出。

每 stage 保持全部 participant prepare/boundary → Stage ready → ordered RHS → RHS ready → 全部 advance；physical publication 后唯一 clock Commit。cross diffusion 读取 partner 已准备的原 Stage workspace，不复制 partner physical state。公共 NumericalCompiler、PlanExecutor 与中央模块装配选择未因这些数学能力增加分支。

本轮数值能力限 serial single patch、单位质量 additive ddt、forwardEuler/classicalRK4、Cartesian constant-coefficient Central2；不包括通用变系数/各向异性、dot/cross/tensor execution、向量 source callback、隐式/Jacobian/DAE/Subcycle、MPI 新资源或任意 PDE。正交 stencil、完整边界和实际 storage 缺失仍 fail-fast。独立解析/半离散参考、native 六例、现役前后对照与已知历史失败见 [A6 报告](../docs/explicit-a6.md)。
