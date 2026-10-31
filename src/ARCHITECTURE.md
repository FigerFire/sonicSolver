# SonicSolver 架构：四轴声明、统一编译、能力驱动执行

更新：2026-10-31。本文定义当前架构边界与下一阶段目标；实现状态另行标明。历史迁移报告保留证据，不作为相互叠加的新规范。stage SST 的局部闭合组已接通普通方程调用及明确的空间系数发布；闭合独立化的构建与未测试边界见 [prepare.md](../docs/prepare.md)，此前系数发布见 [closure.md](../docs/closure.md)，首条关系的历史收口见 [ordinary.md](../docs/ordinary.md)。

## 1. 目标与主干

SonicSolver 组合一个可解释、可执行的数值系统。用户明确指定数学关系、执行过程、状态版本和数值实现；编译器验证并翻译这些选择，不自动寻找求解方法，也不纠正用户的耦合策略。

```text
GUI / CLI
    ↓ 用户设置
IO：读取注册文件与 typed configuration
    ↓
控制层：内置方程、模型、算法声明的选择与组合
    ↓
WHAT + HOW + STATE + WHICH
    ↓
Lowering：数学检查 → 过程/数据检查 → 已选实现的能力检查
    ↓
CompiledSolvePlan + 资源需求
    ↓
Runtime Binding → PlanExecutor → Numerics / Boundary / Linear Algebra
    ↓
Runtime interfaces → Infrastructure backend
```

Math 提供纯数学符号、签名及代数运算；它不拥有场存储、IO、MPI 或物理模型。Numerics 使用 Math 实现离散与求解，不拥有全局时间循环。

主未知量形式、物理系统、耦合算法、时间推进、空间离散各自独立。不能从 density/pressure 推断 compressible/incompressible、EOS 或默认算法；不以 solver family 重新捆绑这些维度。

## 2. 四轴唯一职责

| 轴 | 唯一权威 | 不应承担 |
|---|---|---|
| WHAT | 方程、约束、显式声明的数学关系 | 选择执行顺序、储存第二份执行公式 |
| HOW | 方程调用、循环、作用域、提交顺序 | 猜测物理模型、自动重排耦合策略 |
| STATE | 符号签名、调用使用的状态版本与存储映射 | 重定义方程、另选求解目标 |
| WHICH | 方程/算子/方法及唯一时间 recipe 的选择 | 更改数学关系、补入未声明算法 |

Equation ID、HOW step ID、注册场符号是不同身份。HOW 的第一个引用必须解析为方程，目标必须解析为可写的注册量；两者不要求同名。STATE 和 WHICH 引用 typed ID，不凭字符串位置猜身份。

下例是设计说明，**不是可直接复制运行的配置语法**。Ca 的量纲应与 D/b 一致：

```text
WHAT:
    E_D: D = Ca*b
    E_a: ddt(a) = diffusion(D,a)
HOW:
    Stage:
        S1: E_D -> D
        S2: E_a -> a
    Commit(a)
STATE:
    S1: b(Stage)
    S2: D(Stage), a(Stage)
WHICH:
    E_D: LocalAffine
    E_a.diffusion: Central2Arithmetic
    Time: classicalRK4
```

示例中的 b(Stage) 必须由真实的阶段输入或发布者提供；仅声明一个名字不产生数据。

### WHAT 与 Math

方程是一般等式，`D=f(b)` 与 `D-f(b)=0` 可以表达同一关系。由 HOW 指定求解哪个目标，再由选定方法判断能否求解；不新增 Define 数学权威或强制区分赋值式方程类别。

Math 检查量纲和形状。符号具有七维量纲向量与张量签名；加减要求兼容，算子按自己的签名传播形状和量纲。合法签名不等于存在离散、导数、线性化或求解实现。

普通显式方法接受其明确支持的 `ddt(target)=RHS` 形式。用户负责规范化 `C*ddt(a)=R`，编译器不自动除以 C、不自动推导一般质量矩阵。若写为 `ddt(a)=R/C`，对应输入版本必须有效，除法必须有定义。C 是另一方程的未知量，不妨碍其已发布版本在当前调用中成为已知输入。

守恒映射如 `q=rho*a` 属于 WHAT 数学关系。STATE 可引用该关系提供派生视图，但不得独立再保存另一份 q/rho 公式。

### HOW

保持精简：单方程目标调用 `E -> a`、整体块调用 `[E1,E2] => [a,b]`、Sequence/Loop/Stage 作用域及 Commit。本文不要求新增公开 Update/Solve/Evaluate 指令。一次调用是代数求值、RHS 评价还是隐式求解，由方程、目标和 WHICH 的契约确定。

循环条件引用统一残差能力和显式参数；绝对/相对残差及归一化由残差实现解释。嵌套保持树结构，编译器不得通过扁平化改变作用域、停止条件和状态生命周期。

严格保留用户顺序。先读后写不必然非法：读取显式 OldTime/预测版本可以合法；要求本轮尚未产生的新版本必须拒绝。标准 RK 的阶段快照、全部 RHS 就绪后推进等规则仍需检查。HOW 不能通过覆盖阶段输入暗中改变所选时间方法。

### STATE

对外简洁选择版本，对内保留真实生命周期与存储保障：

| 表达 | 含义 |
|---|---|
| a | 当前物理版本，按调用契约解析 |
| a* | 工作版本 |
| a' | 上一时间层/本步入口快照 OldTime |
| a(Stage) | 当前 RK 阶段版本，必须有真实绑定 |
| deltaP | 独立修正量符号，可借用内部 correction workspace |

`a'` 不表示修正量。HOW 已决定目标身份，STATE 只决定输入版本及目标允许的视图；不重复选择另一个输出变量。字段、版本、内存不是三套 physical authority。

纯 Lazy 关系可以按请求版本读取，不能暗中更新其他状态、推进时间或触发隐藏模型执行。Materialized 数学关系通过普通 HOW 方程调用执行，并发布对应版本。更新后依赖失效、阶段结束和物理提交都必须有明确语义。

### WHICH

默认按方程与算子 occurrence 绑定，必要时以 step 作显式覆盖；不能靠同名算子猜选中哪一项。时间方法只保留一份完整 recipe authority，不另设 stage 数或互相冲突的 explicit/implicit 开关。

选择某方法仅表示用户意图。方法必须实际支持所需算子、场形状、输入版本、边界、网格及并行资源。禁止为了接受组合而静默换方法。

## 3. 控制层和内置模型

内置 C++ typed 声明是主要可读维护入口。方程和模型提供基础数学，SIMPLE/PISO/PIMPLE 等提供 HOW 过程及所需关系；套餐组合这些声明，不能成为第二套运行器。

模型可以增加 term、equation、constraint 或 closure；算法改变求解过程。`SystemContribution` 是现有声明汇集容器，不意味着模型必须另建专用 provider。

用户选择 kEpsilon，应由模型提供未知量、方程和闭合依赖，算法声明将调用插入明确位置。不能只按一个数值优先级解决所有相互依赖；冲突、循环作用域和所读版本必须明确。用户自定义 HOW 与套餐 HOW 的替换也应显式声明。

内置声明恢复已完成，见 [builtin.md](../docs/builtin.md)。外部模板库及展开语言已退役；普通 native 四轴输入保留。不要重复建设模板引擎，也不要为了隐藏数学把全部声明塞进巨型控制类。

## 4. 编译与 provider 边界

三类检查分别给出诊断：

1. 数学合法性：引用、量纲、形状、算子签名。
2. 过程合法性：顺序、作用域、版本可用性、目标写入、stage/commit 契约。
3. 实现可用性：所选数值实现、离散操作、边界、halo、布局和并行能力。

统一 Lowering 协调 HOW 树、WHAT 方程、STATE 使用和 WHICH 选择，生成现有 CompiledSolvePlan；不新建平行计划或第二套执行器。Runtime Binding 绑定已声明资源，PlanExecutor 执行已编译操作，不重新选模型、方程或数值方法。

HOW 生命周期由过程编译统一解释。`IProvider::lowerLifecycle` 及全部生产实现已退役；局部方法不再扫描根 HOW、装饰其他调用或追加物理时间提交。`CompiledEquationCall` 只声明数值输入准备、迭代工作量、局部组前置动作和物理结果发布的操作身份，不携带新控制树。编译器在已有 Step/Loop/Sequence 作用域和显式 Commit 处实现这些要求，按操作身份去重并拒绝 owner 冲突；不创建数学循环或移动方程调用。固定时间 Working predictor 的快照属于其显式带信号祖先循环，不能落到压力内修正循环。

三种发布分别处理：普通代数/块的局部输出仍在调用内发布；Working/Correction COPY 与压力修正保持原位置；物理 Commit 发布完整数值结果后，才执行唯一 `time.commit`。显式 recipe 仍生成快照、Ready、全部 RHS 就绪屏障、advance 和 Close。全部物理量发布后独立闭合 Current 刷新及最终普通关系重求，均先于 clock。非积分的 Stage 关系不获得 RK 增量或独立 clock。

时钟操作由编译器声明 `PhysicalTimeCommit` 能力，宿主只绑定原 StateBundle 的 dt/time/step 更新。Eulerian 的旧 `ee.time.commit` 身份已合并为 `time.commit`，其实际 callback 算术及 owner 不变。重复/Stage 作用域内推进与重复 clock 被拒绝；普通静态关系或显式变量 COPY 不因此推进时间。同步显式组保留其终端 Commit/相邻融合组限制；Eulerian 固定拓扑与完整湍流组校验为只读实现能力检查，不能将其放宽成通用求解支持。Explain 从实际计划区分结果发布与物理 clock；本轮未执行 Explain。迁移审阅与生产构建范围见 [how.md](../docs/how.md)。

Provider 是待退役的迁移期结构，不是最终架构。最终由四轴编译解释生成操作、资源绑定提供输入、所选数值方法执行数学，PlanExecutor 运行冻结计划；不再保留横跨准入、数学、选择、资源与生命周期的大接口。迁移期间普通方法必须消费完整 AST，不凭模型/字段身份准入；未替代的 fixed kernel 仍保留完整能力匹配。

**整式匹配不能直接删掉。** 当前固定 kernel 的匹配保护了“未实现项不能被悄悄忽略”。只有新路径实际消费了全部方程项、系数、符号和版本，才能退役对应匹配。让 validator 和模型共用公式生成器，本身不证明 runtime 能消费新增项。

恒密度 `PressureMomentum` 已按项编译实际 WHAT：一个正的 `ddt(U)`、一个 `div(momentumFlux)`、至多一个 `diffusion(nu,U)` 与右侧显式 `source(expr)`。遍历完整 AST，逐 occurrence 绑定原 PrimitiveUpwind1/Central2 与 primitive acceleration 回调；重复空间项、缩放/符号变换、未知算子或无实际 kernel 的选择明确拒绝。`CompiledPressureMomentum` 只引用原 CompiledSpatialBinding 索引和时间 occurrence；PressureOperators 消费这些冻结绑定，不再重做 AST/recipe/固定操作数匹配。这里 nu 仍是该后端的原固定 transport 投影，执行权重为已解析 mu/rhoConst，不支持任意系数字段或 producer。

压力过程只读验证消费已编译调用及真实 HOW 作用域，已删除 resolver 中 pressureOperatorRequirementsSatisfied 的第二份动量模板和 additionalContributions 整体拒绝。现有压力计划可以显式加入 Current 动态输入或常量输入、独立 TransientWorkspace 输出的 LocalAffine；其本地求值/发布由原 binder 实现，无额外 clock。普通非积分方法仅核验自己的操作，不能因外围压力操作或时间提交而整体拒绝。写 U/p/phi/压力工作量、Stage/OldTime 附加读取、其他时间组、块求解及未实现边界端口仍具体拒绝；非正交限制按实际未校正的压力求解循环检查，插入无关代数不改变该限制。固定 U/p/rho/phi 布局、压力修正关系和其他后端能力保护继续保留。动量逐项编译见 [momentum.md](../docs/momentum.md)，动态反馈与方法退役增量见 [loop.md](../docs/loop.md)；未运行测试、未数值验收。

资源布局限制、块结构限制、真实支持的边界限制仍可保留。诊断必须说明缺少的具体能力，而不是仅报告“模型不匹配”。

## 5. 普通算子、融合实现和数学更新

普通路径的目标是 AST → 已选算子实现 → 数值贡献组合。融合路径消费同一数学声明，在自身确实支持的结构内提供优化；两者不应各保存一份物理方程。

融合 kernel 的完整成员、相邻性与单次 RHS 条件只约束选中该实现的调用。普通独立方法不继承这些限制。禁止暗中重排、重复执行融合 RHS，或在不匹配时自动换成另一实现。

Central2/Central4 的可选性需要真实算子、边界闭合、halo 与面通量实现，不能仅添加枚举。暂不为抽象完整性重写全部 CFD kernel。

数学上的系数更新（如 D、mu_t）逐步进入 WHAT/HOW；分配、绑定、halo 准备仍是执行基础设施。壁距、梯度、滤波和限幅各有真实算法，不能冒充纯局部代数。退役 coefficientPrepare 的某部分之前，必须列清其数学输出及所有消费者。

stageConservative SST 的四条局部关系（sink.k、nu.k、nu.omega、sink.omega）由模型唯一 WHAT 声明，通过普通 LocalAffine 执行。HOW 明确 Stage 调用，位置在原剩余 prepare 之后、完整 flow RHS 之前。两个输运 RHS 仍消费原 rho*nu 扩散及 production-sink 数学，不拆共享通量或改质量加权。步末先发布全部物理未知量、刷新剩余闭合，再重求四条关系，最后唯一 TimeCommit；不复制最后一个 Stage 充当最终关系。

F1 来自原 stageCoefficients 非局部算法，与 mu_t/production/cross 一起由剩余 prepare 发布。F1 使用 solver workspace，mu_t 仍借用原模型数组；实际 Stage/Current 视图有发布 epoch 检查。LocalAffine 和闭合 prepare 使用同一配置的层流动力黏度，实际端口不一致时拒绝，不从 EOS 或默认值替换。

普通 LocalAffine 可借用原标量 STATE 或物化 TransientWorkspace。workspace 默认只发布内部节点；只有声明 WorkspaceBoundary::NearestInterior 的输出，才在内部最终值计算后执行原最近内部节点延拓。当前仅两条迁出的扩散系数请求此准备，局部 sink 不增加边界处理。该准备是所选 workspace 的串行读取区域契约，不是用户新增 HOW 命令、物理 BC 或 MPI halo，也不是在 ghost 输入上重算数学。空间消费者核对实际 publisher 和准备契约；未准备区域或未发布版本不得读取。

普通 Stage 输入的合法性由实际 publisher、版本及 HOW 顺序确定，不以是否为演化未知量判断。已声明的闭合输出读取绑定到实际 prepare；迁出的三项不再被模型准备列表按名称认领。Explain 消费同一编译计划、系数输入与空间准备契约，不维护另一份数学。

SST 的 mu_t、F1/F2、production 和 cross diffusion 仍由原梯度/壁距/限制器算法准备，完整 k/omega 端口、serial 单 patch 与未迁移 fixed-kernel 保护继续存在。剩余闭合已由独立 SSTNonlocalClosure 的五输出数学组与 Stage HOW 调用执行，独立拥有最终 Current 刷新；它不是被积分的 temporal unknown。CapacityTransportSharedFluxCentral2 消费 rho 加权输运 AST、共享通量和独立 ScalarEquationData 资源，不认识 SST 或 k/omega 配对。模型的数组/边界由应用资源端口借出，普通输运不持有模型。闭合内部的原梯度/壁距与配对限制仍是所选非局部实现的实际边界；没有实现任意容量场或任意闭合。once-step RAS、supportedRasMathematics、StageGroup 和通量内核保持。Lowering/STATE 对独立的非积分 Stage 输出作最小适配，保留原屏障与显式顺序。以上是源码接通和构建结论，未运行数值验收。


## 6. 存储、并行与边界

- physical state、clock、equation binding 各有唯一 authority；残差、通量、RK 增量和线性系统属于 solver workspace。
- state/geometry/labels：owner → COPY → replicas。
- 一个共享面只有一个 canonical flux；COPY 后两侧按正负装配，不平均独立 flux。
- residual/source/load/constraint contributions：SUM 到 canonical owner。
- Raw MPI 只存在 infrastructure/backend；上层使用并行语义接口。
- Ghost/ILW 是边界闭合；forcing 按实际方法是源项或预测后修正；DLM/KKT 增加未知量/约束和块代数。它们均不拥有第二套 CFD 时间循环。
- 不改变既有边界公式、近壁阶数或壁面物理来掩盖架构问题。IBM 壁面验收遵守现有 no-slip 要求。

## 7. IO、GUI 与 Explain

GUI 表达用户选择；IO 读取已注册设置，不据文件存在自动启用模块。solver 文件直接位于 solvers/，models、fields、mesh 使用自己的 object 目录；每个目录一份 registry，母键与子文件身份严格一致，禁止路径逃逸和第二份 authority。

EOS、thermo、transport 属于 models；能量输运方程和热物性关系不同。不能从 formulation 或套餐偷偷选择 PerfectGas。可选模型 block absent 为禁用，present 默认启用，enabled:false 不进入执行路径。

Explain 从同一份 authoritative description 和实际 compiled plan 输出：全局流、活动方程、耦合/约束、时间方法、每项离散、边界/并行语义；同时显示调用引用、读取版本、默认/覆盖来源、真实实现限制。禁止一份展示 DSL 配另一份独立执行数学。

## 8. 当前实现与剩余距离

截至本次源码/报告审阅：

| 内容 | 状态与边界 |
|---|---|
| Math/WHAT AST、签名、方程目标引用 | 已有基础；不等于任意 PDE 可解 |
| HOW 树、统一 lowering、typed STATE/WHICH 引用 | 已有实现；方法仍有真实能力限制 |
| 普通代数调用与部分 Stage 系数链 | 已有基础；不能泛称所有 materialized producer 已支持 |
| Step/Stage 作用域及融合 owner 单次 RHS | 已实现迁移；最新修改未数值验收 |
| 内置 typed 声明 | 已恢复；声明来源变化不增加求解能力 |
| 普通 scalar/vector 表达式路径 | 已有；版本、源回调、布局、边界和并行范围有限 |
| 通用隐式/非线性/多步/任意块 | 未完成；已有线性块能力不代表全部支持 |
| 模型无关 provider | 普通 LocalAffine→物化系数→输运已接入一个真实 SST 切片；其余 RAS 整式、SST 非局部系数和边界端口仍专用 |

关键证据：SF_builtinProviders.cpp 的 TurbulenceTransport 与 supportedRasMathematics 保护未迁移 once-step 内核；SF_scalarMethod.cpp 的容量加权输运已不含 SST helper/配对，非局部配对与数学保护属于独立 SSTNonlocalClosure。四条局部系数已退出 prepare 数学与专属输出集合，完整 AST 由 LocalAffine 执行；F1 等剩余输出保留原 prepare，空间准备见 closure.md，不能宣称整套湍流已通用化。

近期顺序：三条主要应用装配、HOW 生命周期与恒密度压力动量逐项编译已在限定范围接通；继续处理其余具体方法的固定数学准入与全部 AST 消费，迁出剩余 provider 编译/绑定职责，再删除 registry/resolver/catalog 和旧接口。数值验收由用户指定算例；mixture、IBM production KKT 和一般隐式分别后置。

## 9. 后续双数学入口

未来 Equation Input 与 Variational Input 汇入同一 WHAT。独立变分前端处理 Action、FreeEnergy、Dissipation、Constraint 及各自规则，生成方程/关系，不重新建立求解器。约束语义复用已有 WHAT/STATE。

KKT 是一类约束数值求解方法；离散约束投影不等同 Hamilton 时间作用量。变分前端不自动生成所有 VOF/Level-set/湍流算法。该方向待现有模块接入与四轴执行能力收口后推进，目前不是已实现能力。

## 10. 验收与历史

默认不运行测试、回归、诊断或架构检查脚本，由用户指定算例。构建与数值验收分开；构建成功不证明方程组合正确。授权的实验全部放 test/，逐例报告实际数学、配置、时间、边界、并行和输出。不得覆盖历史基线。

下一步工作见 [problem20261020.md](../docs/problem20261020.md) 和 [provider-plan.md](../docs/provider-plan.md)。历史实现与验证范围见 [lowering.md](../docs/lowering.md)、[scopes.md](../docs/scopes.md)、[stages.md](../docs/stages.md)、[builtin.md](../docs/builtin.md)；[templates.md](../docs/templates.md) 是已退役外部模板路线的历史记录。

独立闭合调用占据原 Ready 后、普通局部关系前的位置；最终刷新占据全部 physical publication 后的位置。它通过既有 stageLocal 计划语义执行，不寄生于任何输运参与者，不增加 snapshot/RK advance/TimeCommit。旧 SSTStageTransport/provider 与输运自动附加闭合已退役；详见 prepare.md。

## Provider 退役第一阶段：conservative 输入

FlowOperations/ConservativeRHS 不再保存或读取 SolverConfig。已有 CompiledNumericalSystem.conservative 冻结 IBM stencil policy、ILW order、界面通量策略、分子兼容参数及来源，thermophysical/numerics 合法投影只在组合入口解析一次。对流 gamma 来自绑定 FluidStateModel 的能力；T、mu、导热率继续来自模型/缓存与已选 transport 资源，不冻结动态场值。单/多 patch 的模型绑定、halo 与界面策略均消费该编译结果。time/dt 继续使用已有 compiled recipe/policy，不另加配置权威。

守恒边界定义由原 Applicator/Pipeline 持有；恒密度压力边界由其编译 PressureNumericalConfig/PressureOperators 持有。FlowOperations 只额外保留未迁移 conservative corrector 的参数与 gamma；SingleFluidStepper 仅在非恒密度压力路径提取这些 legacy 输入，不存整份配置。对应 conservative pressure/KKT binder 仍使用窄参数，没有修改其数学。输入消费矩阵、构建与未验收边界见 [authority.md](../docs/authority.md)。

## Provider 退役第二阶段：资源装配与操作绑定

conservative 的应用入口不再查询 flow.conservative：由冻结 compiledProgram.stateViews 及 producer 依赖查询所需 conservative backing/分量，借用 Field 原 packed Q。withSingleFieldResources 构造或借用 EOS、模型、边界和字段端口，并在 owner 存活期交给公共 executeEquations；它不绑定操作或拥有时间驱动。独立普通字段也交给同一公共入口。单/多 patch 的重复 runConservative 入口已删除，拓扑特有资源准备保留。

bindSingleFluidOperations 仅在计划请求 FlowStepPrepare/PressurePrepare 时创建对应操作 workspace。FlowOperations 保留实际请求的非时间组操作，snapshot/prepare/RHS/advance/publication 由同一 StageGroup 绑定；不再按 provider 家族集合取得操作。数值实现身份仍随冻结 binding 校验，PlanExecutor 不重新选择算法。新增已有资源能力的普通方程只贡献数学与四轴选择，不修改应用家族判断。

冲突按实际操作与 scope 检查：同一 StageLoop 的不同 temporal owner 不可重叠写同一 storage slice，不可重复 RHS/advance 或发布同一 stage 资源；顺序数学更新仍由原别名/版本/失效契约判断，不一律禁止多写者。packed 路径的物理 clock commit 不可重复或位于重复/stage scope。完整操作覆盖、OpRegistry 重复绑定拒绝及计划 owner 核验保留。

恒密度压力的 U/p、rhoConst、预测/修正及持久面通量已由实际方法/资源需求接入同一应用装配边界，详见下节。Eulerian 的 phase/old-time/shared-pressure 资源也已按冻结契约装配，数值 workspace 仍为独立实现，见本节末尾。压力资源与 packed/phase 的布局、发布、通量和 clock 组合能力仍未实现时明确拒绝；不对 conservative＋普通方程/独立闭合施加家族互斥。multi-patch 普通字段/闭合资源端口仍未实现。

已有独立 SST 闭合、四条局部关系、NearestInterior 延拓和 Stage/Current 最终刷新顺序不变；本轮只完成源码装配并构建，不运行测试或数值验收。实际资源矩阵、剩余分派与增量见 [assembly.md](../docs/assembly.md)。后续继续迁移其他装配边界，再逐方法迁出 provider 编译/绑定，最终删除 IProvider/registry/resolver/catalog；不维持长期新旧双轨。

## Provider 退役第三阶段：恒密度单流体压力资源

requiredPressureInputs 查询实际编译 PressureMomentum 调用，返回已有冻结 PressureNumericalConfig 的引用；核验 U[3]/velocity、独立 p[1]/pressure、正 rhoConst、deltaP Correction 与 canonical pressureFaceFlux views，以及完整预测/修正/发布/clock 操作组。不因 U/p 同名字段、PressurePrepare 单个 OpId 或未消费的数值默认值而启用算法。该查询放在已有 built-in numerical selection 实现中，不新增 IR、selector 或资源中心。

executeSingle、withSingleFieldResources、Environment 和 inspection 不再按 flow.pressure-operators 选择压力存储或布局。网格先提供几何，单场装配一次建立三分量速度 backing 与独立 p；原先随后会被 resize 清除的五变量/PerfectGas 初始化不再执行。Field 的 resize 原本即清零全部分量；原速度初值、p 默认填充值和逐 set 初值覆盖仍在原位置进行。默认 p 填充值沿用已解析参考值，运行中的参考 cell/value 是另一职责，不合并为初值或数学约束的第二份来源。

原 PressureOperators 是真实 Working/Correction、rAU/HbyA、持久 corrected face flux、矩阵/线性求解和压力边界 owner；恒密度路径不再创建无消费者的 conservative Applicator/Pipeline，也不读取 conservative ILW 数值输入。Environment 的单 source zone/单 patch-per-rank 限制和 U/p empty-dimension 规则消费同一资源契约与 compiled 边界；密度/热边界的几何集合输入仍属原应用职责。未新增多 patch/MPI 能力。

压力 binder 由完整 Program 的方法/资源契约构造 PressureOperators；conservative pressure corrector 保留独立 legacy 窄输入。OpId/provider 标签仍验证冻结实现身份，不再决定应用 family。压力与 conservative/phase 组合的拒绝检查实际布局/发布需求；Rhie-Chow 必须具有真实 canonical face view 和对应操作 owner，不能以 provider 已注册作为充分条件。唯一 TimeCommit 检查也覆盖本切片，原 alias-aware/Working/Correction/发布检查保留。

固定符号 U/p/rho/phi/deltaP、pressureOperatorRequirementsSatisfied 的完整数学保护、原边界、松弛、非正交能力与既有时序保持。未实现任意变量重命名、任意动量方程、任意普通字段耦合或 KKT/phase 组合。Explain 只读取实际契约，绑定成功仍需 runtime 核验；本轮未执行 Explain。源码/构建及增量见 [pressure.md](../docs/pressure.md)，未运行测试、未数值验收。

后续进入逐方法迁出 provider 编译/绑定 → 删除 registry/resolver/catalog 和旧接口；phase 资源与输入的后续收口见下节。Provider 尚未全部退役。

## Phase 资源与执行输入收口

requiredPhaseInputs 从实际 CompiledEulerianAssemblyContract 的 phase slot/name、reference closure 和关系组确定相布局，核验 STATE 主状态/primitive alias、共享压力、pressureCorrection、canonical volume/mass flux workspace 与操作绑定。配置存在或 provider 标签不启用资源。execution、inspection、multiPatch 和 validator 已删除 flow.eulerian-pressure 身份分派；与 conservative/恒密度布局的未实现组合继续按实际资源/发布/clock 契约拒绝。

withPhaseResources 持有 PhaseSystem、Turbulence::EquationSystem、PhaseSourceRegistry、PhaseBoundaryApplicator、runtime 和输出 callbacks。源对象只在装配阶段从 typed 参数创建，并与每相实际数学 source extension 对齐，不在数值 cell loop 再判断 enable。输出前核验真实 phase 顺序/reference 与模型分组，保留 initialPressure/referencePressure 的原相等检查及阈值：参考压力修正行设 deltaP=0，随后在当前共享压力上累加修正，不能在本轮擅自解除这项初始压力锚定。

独立 runPhasePlan 借用上述资源，创建数值 owner EulerianStepper、注册原数组/workspace view、绑定服务，再调用原公共 flowLoop；不要求 phase 使用 SingleFluidStepper。Stepper 持有 old-time/outer/flux/矩阵 workspace 与操作 registry，只借用模型和边界，不保存 SolverConfig。输出模型数组由应用保证存活，summary 由数值 owner 保证存活；原归约顺序保留。

新增的窄 PhaseNumericalConfig 只有实际消费的 momentum/pressure relaxation、reference、四类 linear 配置和 shared-pressure BC，由已有 selectNumerics 解析并冻结。方程装配只读该值；dt/CFL/phase transport 已有 compiled 值继续唯一，raw/compiled dt 相等断言删除。模型参数归实际模型对象，动态相状态不冻结成常数。没有将整个 SolverConfig/PressureCorrectionConfig 搬到新对象。

N-1 alpha、每相 mass/momentum/enthalpy 主状态、primitive caches、old-time、interphase/source/湍流准备频率、压力修正发布、paired canonical flux 与最终 commit 均沿用原实现。Workspace 目标可以不注册物理 STATE；packed backing 查询已按这项既有契约区分 numerical workspace，未知物理符号仍拒绝。完整数学匹配和 scope/alias 检查保留；phase 物理发布与统一 TimeCommit 分开，callback 仍更新同一个 StateBundle。没有开放任意 phase/普通字段/IBM 或 multi-local-patch 组合。

Explain 显示同一资源需求、窄输入与仍有限的实现能力，不宣称实际绑定已完成。本轮生产构建成功，未运行测试、算例、诊断或数值验收，见 [phase.md](../docs/phase.md)。三条主要应用身份分派已在限定范围收口，provider 方法编译、完整数学匹配、冻结操作标签与部分专用 binder 仍有真实消费者，尚未全部退役。


## 四轴动态反馈切片（2026-11-01）

本轮结构完成，数值闭环待验收。实际链为 native typed WHAT/HOW/STATE/WHICH → Lowering → 已选数值函数 → 编译资源 → 原 binder/PlanExecutor → 原 RHS。`LocalAffine`、`ConservativeResidual` 及恒密度压力方法已退出 IProvider/ProviderRegistry；按 WHICH 的精确方法身份直接编译，不从模型名称选择实现。已迁移操作按方法/过程实际声明的 operation/owner 冻结，不依赖旧 ProviderCatalog 的模型准入。未迁移消费者才走旧 registry/catalog 验证；不是编译失败后回退。

HOW 的 recipe scope 依据所选实现返回的 temporalResidual 分类，不能仅因 WHAT 有 ddt 就把固定时间压力 predictor 包进 RK。共同阶段仍由原 recipe、一次融合 RHS、Ready/RhsReady/advance/Close 屏障生成，普通关系按明确 Stage/Step 位置求值；相邻融合成员限制不变。统一时间计划在编译后与实际 StageGroup 绑定时均校验。

`component(vector,index)` 与三标量 `vector(x,y,z)` 使用同一 Math 签名和 AST 求值。LocalAffine 可借用真实三分量输入，以独立标量 workspace 为输出。Current reader 借用原物理数组与 generation；materialized workspace 在本步显式发布后可读，Current 依赖写入即失效，不能隐式刷新。Stage/OldTime 继续读取真实阶段/q0，并受原 epoch 约束；没有用 Current 数组代替它们。halo、命名 snapshot、perIteration 与未经实现的 Frozen/Lagged 端口仍明确拒绝。

`source.field` 编译 RHS 的唯一 source(expr)，冻结输入版本和纯局部计算，并在启动时绑定 reader。密度动量使用力密度、能量使用功率密度，直接累加原 Residual.source；恒密度动量使用加速度，保留原 dt 权重，不补除 rho。HOW 负责显式计算其输入 workspace；源读取不触发模型 correct/prepare。新增同数值形式的关系改名不要求新增 provider 或 application 分支。局部源内部算子当前不可独立替换，显式选择未实现的内部实现会被拒绝。

密度路径仍固定五变量 rho/rhoU/rhoE 与原可用 EOS/边界/空间方法；恒密度压力仍固定 U/p/rho/phi/deltaP、PrimitiveUpwind1/Central2、rhoConst 和持久 Rhie-Chow 通量。压力动量按项消费；压力修正/速度/压力/通量/松弛等剩余固定 kernel 以结构及扩展规则检查保护，不将新公式只打印后忽略。源项可重复显式相加；重复对流/扩散在融合方法尚无累加实现时拒绝。

两条真实数据流、准备输入和未验证范围见 [loop.md](../docs/loop.md)。普通关系/反馈绑定当前限串行单 patch；密度 field source 需要实际 synchronous stage group，旧 unsplit fused executor 未提供其 reader 端口时拒绝；不宣称压力可耦合任意时间输运。once-step RAS、Eulerian、IBM、Level-set、thermal/mixture、普通显式输运/线性块等仍有旧 provider 消费者，后续按实际数值方法迁移。旧四条 SST 关系、独立非局部闭合、NearestInterior 延拓和最终重求顺序没有重做。生产构建不代替数值验收，历史压力 MPI 阈值与 Level-set 参考问题保留。
