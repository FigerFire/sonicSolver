# SonicSolver 架构：控制层展开声明，核心执行通用数学与数值操作

## 从哪里读源码

以 `predictor: E_U -> U` 为例。`E_U` 保存动量方程，`viscU` 标识其中一项。启用 SST 动量反馈时，控制组合按该引用替换应力项，保留其他项及其权重。下面的入口也适用于默认预制生成的普通调用。

| 想做什么 | 从哪里读 | 用 `E_U` 追踪 |
| --- | --- | --- |
| 添加由既有算子组成的方程 | [SF_nativeBinding.cpp](app/application/model/SF_nativeBinding.cpp) 读取 native 声明。[SF_model.cpp](app/application/model/SF_model.cpp) 展开模型，[SF_equationContribution.cpp](solver/system/SF_equationContribution.cpp) 组合 WHAT | `momentumStress` 引用 `E_U/viscU`。`finalizeMathematics` 替换唯一的黏性项，保留 `appliedForce` 等其他项 |
| 改求解顺序 | [SF_declarations.cpp](solver/system/SF_declarations.cpp) 解码用户 HOW，[SF_pressureCoupling.cpp](solver/system/SF_pressureCoupling.cpp) 生成压力预制。[SF_lowering.cpp](solver/system/SF_lowering.cpp) 解析调用并生成操作树 | HOW 决定先发布黏度与 `SST.muEff`，再运行 `predictor`。过程契约检查依赖，[PlanExecutor](solver/run/SF_planExecutor.cpp) 执行冻结的顺序 |
| 选输入版本 | [SF_stateRealization.cpp](solver/system/SF_stateRealization.cpp) 解析调用的读写需求，[SF_stateRealizer.cpp](solver/system/SF_stateRealizer.cpp) 绑定实际视图 | 预测调用读取已发布的 `SST.muEff: Current`。PISO 的预测目标是 Physical U，SIMPLE/PIMPLE 可通过 HOW 指定 Working U |
| 更换离散方法 | [SF_builtinNumerics.cpp](solver/system/SF_builtinNumerics.cpp) 解析 WHICH，[SF_builtinProviders.cpp](solver/system/SF_builtinProviders.cpp) 将该次调用交给选中的方法 | [SF_cell.cpp](solver/system/SF_cell.cpp) 在同一入口核对 Cell 方程、逐项权重和所选 kernel。`viscU` 的 `stress.central2.arithmetic-faces` 绑定实际算术面牵引 |
| 追踪计算结果 | [SF_singleFluidBinding.cpp](solver/algorithm/SF_singleFluidBinding.cpp) 先挂输出、再接读取、最后装操作。[SF_cellBinding.cpp](solver/algorithm/scalar/SF_cellBinding.cpp) 绑定该次预测工作区 | [SF_cellEquation.h](solver/discretization/SF_cellEquation.h) 装配面牵引。[SF_faceBinding.cpp](solver/algorithm/SF_faceBinding.cpp) 从 `predictor` 借用预测速度和 `dt/rho` 响应，压力校正后发布实际 `phi` 数组 |

局部仿射、非局部 SST 闭合、普通标量输运和线性块的编译入口分别是 [SF_affine.cpp](solver/system/SF_affine.cpp)、[SF_sstMethod.cpp](solver/system/SF_sstMethod.cpp)、[SF_scalarMethod.cpp](solver/system/SF_scalarMethod.cpp) 和 [SF_block.cpp](solver/system/SF_block.cpp)。方法编译核对本次数学，冻结后的跨调用兼容性由 [SF_executionContract.cpp](solver/system/SF_executionContract.cpp) 核对。迁移细节与验证范围见 [method.md](../docs/method.md)。

## 公式、物理模型、流程算法与数值方法（2026-11-17）

本节记录用户最新设计，优先用于后续重写。各模块消费自己需要的输入，提供明确结果。用户负责选定耦合策略、调用顺序、迭代次数及停止准则。开发优先补齐执行操作并删除固定身份准入；高级组合诊断后置。

| 对象 | 内置内容 | 进入四轴的职责 |
| --- | --- | --- |
| 公式 | NS、mu_t、MRF、输运与约束等数学关系，以及可绑定符号和追加项位置 | 提供 WHAT 表达式，声明实际依赖 |
| 物理模型 | 选择和组合 k/omega、mu_t 等可复用公式，提供常数、壁面关系及应力/热输运贡献 | 主要生成 WHAT，声明所需 STATE 与可调度的调用片段 |
| 流程算法 | SIMPLE/PISO/PIMPLE、显式推进等预制过程，引用基础公式与活动模型的调用 | 主要生成 HOW 和 STATE，提供可覆盖的 WHICH 默认选择 |
| 数值方法 | WENO3/WENO5、中心差分、通量构造、矩阵装配和线性求解等具体实现 | 接收 WHICH 绑定的操作与数据，返回数值贡献或结果 |
| 控制层 | 用户选择、内置声明、显式覆盖和模型接入位置 | 展开并组合一份 WHAT/HOW/STATE/WHICH，交给统一编译入口 |

mu_t 闭合独立保存为可复用的内置公式及其实际计算实现。SST、k-epsilon 等模型选择各自的闭合关系，绑定输入、参数和输出符号，通过控制层加入 WHAT。公式可以继续放在物理模型目录下的独立函数或短文件中，复用现有公式机制。

HOW 通过 `E_mu -> mu_t` 明确安排计算位置，STATE 指定各输入版本并绑定共享输出，WHICH 选择表达式求值及其中空间算子的实现。梯度、壁距等前置量由声明的关系或实际数值输入提供。模型提供默认调用片段与接入位置，控制层完成组合；用户显式 HOW 决定最终顺序。

动量、能量、其他输运或闭合方程可读取同一次发布的 mu_t。热扩散等消费者通过自己的关系将 mu_t 转换成所需系数，例如结合 cp、Pr_t 形成湍流导热系数。更新时机由 HOW 决定，消费者按 STATE 读取选定版本。每个独立物理系统/相可绑定自己的输出符号，多次模型引用同一公式时由控制层按实际实例复用。

不同模型的湍流黏度公式具有独立数学身份，共享输出量的物理含义。数学关系与其执行实现保持一致，内部临时值和跨调用结果按真实生命周期分配存储。现役 SST 用独立的 `SSTCentral2` 调用发布 mu_t。其余非局部关系与普通系数关系读取同一次发布。公式参数可绑定实际符号，原物理模型的数组与边界仍由原 owner 提供。

算法预制可调用 NS 公式，安排源项与模型调用的执行位置。MRF 数学由 MRF 模块提供，控制层将其追加到所选方程；预制使用声明的接入位置及依赖完成组合。用户显式 HOW 按原顺序编译。时间 recipe 的阶段系数仍由唯一时间方法持有。

数值方法按算子语义共用输入输出接口。WENO3/WENO5 可接受同类重构输入，各自声明实际 stencil/halo 宽度。字段名、模型名和预制名仅在声明选择时使用，普通数值操作消费绑定后的数据。需要的几何、边界、矩阵与并行资源由对应实现提供。

### 检查归属与删除方向

核心保持解析、符号解析、数学签名和实际数据绑定所需的局部契约。取消根据模型身份、惯用组合或收敛偏好拒绝用户过程的规则。数值方法只报告本次操作需要的输入和实际实现，避免跨层遍历整个模型组合。反复出现的同一检查应收敛到拥有该契约的一处。

固定点更新可以由一次调用、固定次数 Loop 或用户选定的 While 控制；松弛因子由用户选择。残差定义可独立注册，HOW 按需要引用。声明中存在其他模型常用的符号，不构成当前方程的冲突。

现役部分 provider 仍调用写死数学的整组内核。退役任务同时迁出其数学、调用顺序和数据读取，接通实际消费 WHAT/WHICH 的操作后删除原入口、整式匹配与专用流程。旧类和接口可以直接重写。工作项以实际缺少的计算描述，例如变系数向量应力、质量通量输运、阶段压力校正，避免为每个模型组合新增适配宿主。

### 现役绑定边界

单流体启动时为每次调用绑定工作区，验证并封存操作表。单实例标量端口按实际借用核对。固定区间按明确的准备、步长提议、快照和接受动作接线。仍重绑的 Eulerian 继续逐次验证。历史改造见 [kernel.md](../docs/kernel.md)。

已删除 once-step RAS 对无关 omega/epsilon 注册的排斥，固定点方法强制离散残差 AND 循环的限制，固定点松弛因子的 (0,1] 限制，以及所有残差必须被 Loop 引用的要求。固定点算术、数据快照、发布和时间提交使用原实现。

恒密度 fixedVolume 由控制层将具体方程与黏性 occurrence 传给模型，模型提供 muEff 和显式偏应力数学。普通 CellExplicit 按 WHICH 装配算术系数面牵引，压力按实际 predictor occurrence 绑定响应。模型数组所有者按冻结调用的 STATE 发布效果和物理存储切片提供视图，逻辑输出名可变。fixedVolume Step 可以把现有黏度调用绑定到另一条已选关系。原保守 Stage 发布、Current 刷新和旧 once-step 扫描仍有不同的执行契约。代码路径和边界见 [bind.md](../docs/bind.md)，原应力算法见 [stress.md](../docs/stress.md)。能量反馈、变密度压力输运、压力 RK、MPI 应力面交换及其他未迁移 provider 仍待实现。

## 时间作用域、残差与共享结果（2026-11-10）

本节为最新设计契约。Push/While、独立残差及常密度冻结系数输运已接通生产源码；实际范围见 [push.md](../docs/push.md)。生产构建与数值验收分开，未运行测试。

HOW 使用 `Push(dt)` 表示一个物理时间步作用域。入口确定入口时间、dt、目标时间及必需快照；内部按已选 recipe 执行；成功结束时接受结果并通过唯一内部提交更新 clock。`Loop(n)` 重复 n 次，`While(res(p))` 通过独立残差模块决定继续。循环本身不推进物理时间。Stage 由 recipe 确定状态及求值时刻；隐式/固定时间迭代消费同一 dt 与入口快照。

`res(p)` 返回继续迭代布尔值，未收敛为 true，已收敛为 false。残差定义、容差、范数、归一化、绝对/相对策略和数值记录由残差模块管理，HOW 仅引用结果。每个残差绑定明确方程调用及数据版本；同一字段存在多个残差定义时通过配置明确选择。首次求值、检查时机、最大迭代及异常行为须具有明确配置和实现契约，缺失或非有限残差不能视为收敛。

跨调用消费的扩散系数、生产/耗散项和通量由 STATE 注册并绑定真实存储，生产调用写入并发布，消费调用读取同一版本。绑定阶段分配可复用数组；局部临时值按实际生命周期安排。相同数学含义的通量允许保持同一符号并由 WHICH 替换构造方法，体积/质量/动量通量保留各自量纲、形状和位置。预测/校正版本可在旧值无消费者时复用存储。

模型贡献源项、系数、通量/应力以及附加方程。通用方程的源接口追加实际表达式；未知量依赖及隐式装配由所选方法处理。规范化形式由用户或模型声明提供，普通显式方法继续消费受支持的单位时间质量形式。

编译器检查数学、过程和实现能力并冻结契约；绑定器连接输入输出存储；数值核消费数据及必要几何资源。普通 scalar/transport 方法以算子、数据形状、版本和资源为契约，模型身份与预制流程由控制层处理。

## 字段布局与模型输入契约（2026-11-09）

本节记录最新设计决定。以下为目标契约；SST 当前接通范围及未验收状态见 `docs/sst.md`。

模型生成方程、闭合、源项表达式及依赖。k、omega、epsilon 分别具有独立方程身份和调用目标，模型组合选择所需方程及闭合。源项通过实际表达式追加到目标方程；仅在需要跨调用发布时物化源场。gravity/MRF 的贡献形式按具体目标数学确定，速度形式的 g 与守恒动量形式的 rho*g 分别进入 WHAT。

字段需求由 WHAT 的数学依赖、HOW 的求解目标、STATE 的注册及版本、WHICH 的资源要求共同确定。编译与绑定阶段按需求建立存储。常量、分量视图、Lazy 关系和局部临时量可以无需独立数组。第一阶段在 case 编译和绑定时确定布局；运行中动态增删字段属于后续能力。

STATE 保持数学符号、量纲、张量形状、位置和逻辑版本权威。Field/ScalarField/工作区承载实际数据，packed AoS/SoA、独立标量/向量和面存储均可作为实际布局。复用已有 StateRegistry、StateRealization、DistributedFieldView 与资源绑定，统一物理字段 owner。现有 Field 的几何与热力学缓存依赖需要逐消费者迁移。

普通数值核通过编译后的输入输出读取器和明确几何资源访问数据。速度三分量属于数学形状；RHO/RHOU/RHOE 偏移和五变量布局属于具体存储后端。选用特殊融合内核时允许声明实际布局能力。更换守恒量还需对应方程、闭合、通量、边界、初始化及持久化能力。

U、U*、U' 表示逻辑值及生命周期，U' 继续表示 OldTime。旧内容没有后续读取者时允许原位覆盖，同一数组可依次承载入口、预测和接受结果。空间 stencil、RK 入口快照或固定时间迭代仍需旧值时必须保存。RHS、阶段缓存和线性代数工作区归执行实例。别名复用必须同时维护版本、边界准备和发布代次。

时间区间保存入口时刻、dt、入口快照和目标时刻。Stage 按已选 recipe 给出阶段状态及求值时刻；Loop 在所属时间区间内重复；唯一内部 Commit 接受结果并提交 clock。dt 可由用户固定或已实现的稳定性策略计算。隐式时间项同样消费 dt 和入口快照。公开 Push 作用域遵循2026-11-10契约；子步分裂和多速率需另行明确语义。

本轮已在限定范围接通声明驱动存储与真实读取端口：SST使用显式rhoU/rho关系选择原packed/cache梯度，
或使用注册named速度/密度及其真实边界reader；纯Step闭合可读Current/Working，不积分、不继承Stage一次调用限制。
字段需求复用原FieldDataRequirement与binder，模型公式引用的已注册输入也进入native STATE；不新建字段registry。
源项按真实目标等式及显式Target/Form生成，删除全局U名称单位判断，精确方程引用不重复扩展family。
原模型数组owner和复杂Field/thermo布局仍保留；新增显式fixedVolume客户端消费冻结系数与校正phi；恒密度固定时间动量反馈已接线，专用湍流墙面、变密度压力容量与 MPI 应力仍有限。
详见[fields.md](../docs/fields.md)：生产构建完成，未运行测试、未数值验收。

更新：2026-11-07。本文定义当前架构边界与下一阶段目标；实现状态另行标明。历史迁移报告保留证据，不作为相互叠加的新规范。stage SST 的局部闭合组已接通普通方程调用及明确的空间系数发布；闭合独立化的构建与未测试边界见 [prepare.md](../docs/prepare.md)，此前系数发布见 [closure.md](../docs/closure.md)，首条关系的历史收口见 [ordinary.md](../docs/ordinary.md)。

## 设计契约与审查基准（2026-11-05）

**目标是控制层把预制方程、模型和算法展开成四轴，核心按实际数学与已选数值操作执行。不是让一个压力或其他物理专用执行器逐步接纳更多附加方程。** 本节是后续功能设计和任务规划的基准；本文的迁移记录说明实现状态，不将过渡限制升级为架构规则。

### 用户入口与默认展开

用户可选择 NS、PISO、SIMPLE、PIMPLE、湍流、源项和 IBM 等预制，也可明确提供四轴。控制层收集所选模块的方程、字段需求、闭合、算法过程及数值选择，生成唯一 typed WHAT/STATE/HOW/WHICH。预制是声明生成函数，不持有 CFD 时间循环或另一个执行器。

默认 HOW 只在声明组合阶段根据作用域、数据版本和前后约束构造；完整 authored HOW 不自动插入或移动数学操作。无法确定位置、存在依赖环或未实现版本时准确诊断，不私自滞后、加迭代或换方法。普通数值方法不认识 PISO 外循环、SST 等套餐身份。

模型可能同时提供方程、闭合和数值 kernel，但两种职责必须分开：控制层贡献数学；Numerics 通过明确输入输出实现对应数学。不允许只声明一套用于 Explain 的公式，运行时再依据模型名调用另一套固定公式。

### 核心执行单位与边界

编译器按 HOW occurrence 解析 WHAT、STATE、WHICH。每次调用具有自己的已编译贡献、输入视图、输出、数值资源需求和发布契约；同一 equation ID 被多次引用时不能共享可变调用状态。不得用全程序唯一动量/cell 方程槽限制受支持的多个方程。

数学项及其系数、符号、重复贡献必须被真实消费。普通方程不能因为不叫 momentum/k/omega 或不属于某个 provider 而拒绝。选择某种方法仍需检查真实算子、边界、布局、版本和实现能力。具体方法名是选择真实实现的合法键，物理模型名不是通用执行准入条件。

Runtime Binding 只将编译资源连到实际存储和数值实现。PlanExecutor 执行冻结过程，不再解释模型名来决定方程或循环。可以保留真实的线性求解器、约束块、通量构造或优化融合内核；它们不得拥有隐含的完整 CFD 求解过程。

### 状态、时间与通量

U/U* 是逻辑视图，不等于强制分配两份物理数组。WHAT 仍使用数学符号 U，STATE 指定读取哪份值；HOW 指定目标。旧值不再被读取时允许原位复用；空间 RHS 尚需旧邻点、RK 或固定时间关系需要快照时保留真实缓冲。先满足数学读取语义，再简化内存。

调用结束必须完成自己声明的输出发布；边界准备作用于实际选中的版本。编译器可生成已选数值方法需要的技术准备，不能隐式重求闭合或改变用户耦合顺序。

RHS/残差或矩阵评价、状态推进/固定时间迭代更新、结果接受、物理 clock 提交分别有明确职责。循环不自动代表时间推进；固定时间迭代不能反复累加完整 dt。唯一时间 recipe 和唯一 clock 保持权威。

物理通量属于 WHAT 数学；SW/Roe 等通量构造及 WENO/TENO 等重构由 WHICH 选择，实际组合须兼容。共享 phi 的构造/修正由 HOW 调用、STATE 发布，消费者读取同一份结果。Rhie–Chow 不专属于 PISO，也不是任意 div 项的通用替代。仅在内部使用的临时通量无需强制公开注册。

### 当前源码与目标的差距

以下以现役源码调用链为据；生产构建与未运行数值验证的边界见本轮报告。数值支持范围不能由架构推导。

| 边界 | 当前状态 | 是否已达长期目标 |
| --- | --- | --- |
| 普通 Cell 方程实例 | SF_cell 编译与空间 lowering；CompiledEquationCall 内持有 CompiledCellEquation；SF_cellBinding 每 occurrence 绑定实例 | 独立实例主链已接通；仅覆盖所选方法的有限数学 |
| 普通向量与压力 | CellEquationExecution 不持有 PressureOperators；压力显式消费所引用调用的响应 | 宿主反向依赖已解除，不能再把旧压力宿主依赖列为现役事实 |
| scalar/vector | 共用字段端口唯一性、Current/Working shape读取与物理区间戳；cell 与 nodal/transport 保留不同边界和数值路径 | 应继续消除重复机制，不能为合并代码更改离散 |
| 时间接口 | 普通cell Stage接入ForwardEuler/SSPRK3/RK4，与普通scalar共享StageGroup；Step/压力响应保持Euler；普通scalar新增显式partitionedFixedEuler终点迭代及缺陷判据 | 核心时间接入仍有限，不等于所有算法可组合 |
| 面资源与压力 | FaceFluxExecution从实际调用绑定预测面资源，压力借同一数组；独立面不创建矩阵/压力宿主；仍有限 U/p/phi 布局 | 单个面布局接通；多系统、任意面布局尚未实现，未数值验收 |
| 模型 | Stage SST 已退出整体模型闭合/专用provider编译，使用独立黏度调用与可共享的其他非局部关系；once-step RAS、Eulerian、IBM等仍有固定内核及旧provider消费者 | 必须逐真实数学迁移；能力矩阵见 equation.md，cell轮没有新增这些组合能力 |
| 数值证据 | cell/equation/working/face报告以构建与源码证据为主 | 未运行不能宣称默认/显式、多实例和数值保持已验收 |

关键源位置：`solver/algorithm/scalar/SF_cellBinding.cpp`、`solver/system/SF_scalarMethod.cpp`、`solver/system/SF_numericalSystem.h`；实例轮历史见 [cell.md](../docs/cell.md)，面资源及剩余职责收口见 [cell-closure.md](../docs/cell-closure.md)。具体检查及数值范围以现役代码和实际运行证据为准。

### 重写、精简与验收原则

旧代码没有保留义务。保留数学算法、数值语义、并行不变量和历史证据；允许替换类、接口、目录及调用链。优先删除重复 authority、重复数学、过渡转发和固定身份分派，而不是在旧框架外新增 Adapter。新路径接通后删除被替代生产路径，不长期维护双体系。

不按总行数承诺砍掉算法。四轴、编译和绑定应精简；网格、边界、并行、IO及实际数值 kernel 有独立价值。每次重构报告删除了什么职责、仍有什么专用执行骨架，不能只统计改名或文件拆分。

功能验收基准：受支持的方程改名、增添已支持源项、实例化多个普通方程、选用合法输入视图，不需要增加中央模型分支；预制与等价显式声明进入同一链；真实不支持的算法准确说明缺失能力。结构接通、构建通过和数值验证分别报告。

### 后续任务提示词规则

每次生成或更新开发任务提示词前，必须重新阅读本架构文档，并核对最新相关迁移报告和源码；不能仅沿用聊天中的旧结论。提示词明确哪些属于设计契约、哪些是现役能力、哪些是本轮待修改内容。若用户改变设计，先同步本文，再据此规划任务。历史提示词/报告不是额外长期规范。

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

保持精简：单方程目标调用 `E -> a`、整体块调用 `[E1,E2] => [a,b]`、Sequence/Loop/While、根Push与Stage频率。变量Commit保留，物理Commit由Push内部生成。本文不要求新增公开 Update/Solve/Evaluate 指令。一次调用是代数求值、RHS 评价还是隐式求解，由方程、目标和 WHICH 的契约确定。

Loop只重复次数，While通过res(id)引用残差模块；绝对/相对残差及归一化由残差实现解释。嵌套保持树结构，编译器不得通过扁平化改变作用域、停止条件和状态生命周期。

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

HOW 生命周期由过程编译统一解释。`IProvider::lowerLifecycle` 及全部生产实现已退役；局部方法不再扫描根 HOW、装饰其他调用或追加物理时间提交。`CompiledEquationCall` 只声明数值输入准备、迭代工作量、局部组前置动作和物理结果发布的操作身份，不携带新控制树。编译器在已有 Step/Loop/Sequence 作用域和显式 Commit 处实现这些要求，按操作身份去重并拒绝 owner 冲突；不创建数学循环或移动方程调用。固定时间 Working predictor 的入口旧值由普通方程区间准备保存一次；迭代入口值由引用它的压力外循环请求，不能在压力内修正循环覆盖旧时间。

三种发布分别处理：普通代数/块的局部输出仍在调用内发布；Working/Correction COPY 与压力修正保持原位置；物理 Commit 发布完整数值结果后，才执行唯一 `time.commit`。显式 recipe 仍生成快照、Ready、全部 RHS 就绪屏障、advance 和 Close。全部物理量发布后独立闭合 Current 刷新及最终普通关系重求，均先于 clock。非积分的 Stage 关系不获得 RK 增量或独立 clock。

时钟操作由编译器声明 `PhysicalTimeCommit` 能力，宿主只绑定原 StateBundle 的 dt/time/step 更新。Eulerian 的旧 `ee.time.commit` 身份已合并为 `time.commit`，其实际 callback 算术及 owner 不变。重复/Stage 作用域内推进与重复 clock 被拒绝；普通静态关系或显式变量 COPY 不因此推进时间。同步显式组保留其终端 Commit/相邻融合组限制；Eulerian 固定拓扑与完整湍流组校验为只读实现能力检查，不能将其放宽成通用求解支持。Explain 从实际计划区分结果发布与物理 clock；本轮未执行 Explain。迁移审阅与生产构建范围见 [how.md](../docs/how.md)。

Provider 是待退役的迁移期结构，不是最终架构。最终由四轴编译解释生成操作、资源绑定提供输入、所选数值方法执行数学，PlanExecutor 运行冻结计划；不再保留横跨准入、数学、选择、资源与生命周期的大接口。迁移期间普通方法必须消费完整 AST，不凭模型/字段身份准入；未替代的 fixed kernel 仍保留完整能力匹配。

**整式匹配不能直接删掉。** 当前固定 kernel 的匹配保护了“未实现项不能被悄悄忽略”。只有新路径实际消费了全部方程项、系数、符号和版本，才能退役对应匹配。让 validator 和模型共用公式生成器，本身不证明 runtime 能消费新增项。

`CellExplicit` 使用普通向量 target、严格 `ddt(target)` 左侧和完整加性 RHS。`CompiledCellEquation` 冻结实际 term 索引/权重、梯度 kernel、显式 capacity 常量和数值响应；重复空间/源项按 AST 顺序累加。PrimitiveUpwind1 自输运、Central2 常系数扩散、串行正交网格的算术面变系数应力、实际源 reader 和受支持的 `-grad(scalar)/positiveCapacity` 均进入同一 RHS。预制把常 transport 的 mu/rhoConst 声明为 nu 常量。资源与区间由每次普通 Cell 调用绑定，支持无压力向量、双实例及与普通标量共用 recipe 的 Stage。其他空间格式、一般质量矩阵和多个不同 recipe 仍超出实现范围。

压力过程只读验证消费已编译调用及真实 HOW 作用域，已删除 resolver 中 pressureOperatorRequirementsSatisfied 的第二份动量模板和 additionalContributions 整体拒绝。现有压力计划可以显式加入 Current 动态输入或常量输入、独立 TransientWorkspace 输出的 LocalAffine；其本地求值/发布由原 binder 实现，无额外 clock。普通非积分方法仅核验自己的操作，不能因外围压力操作或时间提交而整体拒绝。写 U/p/phi/压力工作量、Stage/OldTime 附加读取、其他时间组、块求解及未实现边界端口仍具体拒绝；非正交限制按实际未校正的压力求解循环检查，插入无关代数不改变该限制。固定 U/p/rho/phi 布局、压力修正关系和其他后端能力保护继续保留。动量逐项编译见 [momentum.md](../docs/momentum.md)，动态反馈与方法退役增量见 [loop.md](../docs/loop.md)；未运行测试、未数值验收。

资源布局限制、块结构限制、真实支持的边界限制仍可保留。诊断必须说明缺少的具体能力，而不是仅报告“模型不匹配”。

## 5. 普通算子、融合实现和数学更新

普通路径的目标是 AST → 已选算子实现 → 数值贡献组合。融合路径消费同一数学声明，在自身确实支持的结构内提供优化；两者不应各保存一份物理方程。

融合 kernel 的完整成员、相邻性与单次 RHS 条件只约束选中该实现的调用。普通独立方法不继承这些限制。禁止暗中重排、重复执行融合 RHS，或在不匹配时自动换成另一实现。

Central2/Central4 的可选性需要真实算子、边界闭合、halo 与面通量实现，不能仅添加枚举。暂不为抽象完整性重写全部 CFD kernel。

数学上的系数更新（如 D、mu_t）逐步进入 WHAT/HOW；分配、绑定、halo 准备仍是执行基础设施。壁距、梯度、滤波和限幅各有真实算法，不能冒充纯局部代数。退役 coefficientPrepare 的某部分之前，必须列清其数学输出及所有消费者。

stageConservative SST 的四条局部关系（sink.k、nu.k、nu.omega、sink.omega）由模型唯一 WHAT 声明，通过普通 LocalAffine 执行。HOW 明确 Stage 调用，位置在原剩余 prepare 之后、完整 flow RHS 之前。两个输运 RHS 仍消费原 rho*nu 扩散及 production-sink 数学，不拆共享通量或改质量加权。步末先发布全部物理未知量、刷新剩余闭合，再重求四条关系，最后唯一 TimeCommit；不复制最后一个 Stage 充当最终关系。

Stage SST 在同一 HOW 中先发布独立 mu_t，再用其实际 Stage 视图计算 F1、生产与交叉扩散的已选输出。
`CompiledSstClosure` 按 occurrence 保存已编译数学、参数与输入输出端口。模型的动黏度数组由原 owner 提供，其他工作量由调用持有。
Packed/cache 与独立速度 reader 继续选择各自原 Central2 实现。末阶段物理接受之后，黏度、其余闭合与 LocalAffine 依序以 Current 重求，最后由唯一 TimeCommit 提交。

普通 LocalAffine 可借用原标量 STATE 或物化 TransientWorkspace。workspace 默认只发布内部节点；只有声明 WorkspaceBoundary::NearestInterior 的输出，才在内部最终值计算后执行原最近内部节点延拓。当前仅两条迁出的扩散系数请求此准备，局部 sink 不增加边界处理。该准备是所选 workspace 的串行读取区域契约，不是用户新增 HOW 命令、物理 BC 或 MPI halo，也不是在 ghost 输入上重算数学。空间消费者核对实际 publisher 和准备契约；未准备区域或未发布版本不得读取。

普通 Stage 输入的合法性由实际 publisher、版本及 HOW 顺序确定，不以是否为演化未知量判断。已声明的闭合输出读取绑定到实际 prepare；迁出的三项不再被模型准备列表按名称认领。Explain 消费同一编译计划、系数输入与空间准备契约，不维护另一份数学。

普通输运 `CapacityTransportSharedFluxCentral2` 消费 rho 加权 AST、共享通量和独立边界端口。SST 独立闭合只接入串行单 patch，壁距、限制器和空间梯度按其显式实现获取。Once-step RAS 的整式匹配和原位扫描仍服务其真实消费者。新动量路径按面产生一次定向牵引，普通装配按相反符号使用。MPI 的 canonical face COPY 与残差 SUM 原协议继续适用，新应力尚缺分布式牵引实现。


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

关键证据：SF_builtinProviders.cpp 的 TurbulenceTransport 与 supportedRasMathematics 保护未迁移 once-step 内核；SF_scalarMethod.cpp 的容量加权输运已不含 SST helper/配对，非局部数学和签名校验属于独立 SSTCentral2。四条局部系数已退出 prepare 数学与专属输出集合，完整 AST 由 LocalAffine 执行；F1 等剩余输出保留原 prepare，空间准备见 closure.md，不能宣称整套湍流已通用化。

近期顺序：三条主要应用装配、HOW 生命周期与恒密度压力动量逐项编译已在限定范围接通；继续处理其余具体方法的固定数学准入与全部 AST 消费，迁出剩余 provider 编译/绑定职责，再删除 registry/resolver/catalog 和旧接口。数值验收由用户指定算例；mixture、IBM production KKT 和一般隐式分别后置。

## 普通方程残差与时间recipe（2026-11-06）

CellExplicit是普通cell空间方法，不再绑定Euler时间身份。显式frequency:Stage声明同一WHAT进入现有
ForwardEuler/SSPRK3/classicalRK4 participant；Step及原压力预测保持Euler关系；新增partitionedFixedEuler时压力仍按t_old采样，普通scalar按t_new迭代。CellEquation::evaluate
产生完整不含dt的RHS，梯度、源项及几何归一化均在其中；Step Euler和StageGroup使用同一空间入口。
scalar/cell共用小型recipe组合函数；StageGroup保持全snapshot、全prepare、全RHS屏障、全advance、最终发布及唯一clock。
LocalAffine Stage关系按HOW显式执行，不获得积分权。每occurrence拥有真实OldTime/Stage/RHS，Physical不被Stage覆盖。

Stage空间读取使用实际版本的prepareSpatial与stageTime；原cell/nodal边界差异保留。补齐cell的boundaries
资源声明，修复被公共需求整理清掉boundaryRequired的接线；该修复和RHS合并的浮点结合变化需数值验收。
RhieChow不得借RK participant的阶段系数伪造Euler响应，压力固定时间工作量、面资源及COPY/SUM仍独立。
旧方法名CellExplicitEuler明确报退役错误，输入需显式迁移，无fallback。

普通Step与Stage混用、普通vector固定终点迭代、多版本任意空间输入、分布式Stage与任意packed切片仍有限制。
这不是压力/IBM/湍流RK能力扩展。生产构建成功，未运行测试/数值验收；完整映射、浮点变化和候选见
[time.md](../docs/time.md)及test/time/；普通scalar固定时间接线见[iterate.md](../docs/iterate.md)，两轮均待用户指定数值验收。

## 普通固定时间迭代（2026-11-06）

用户明确选择完整recipe `partitionedFixedEuler`，不从Loop推断时间数学。普通标量的
`FixedPointCentral2`/`FixedPointVolumeFluxUpwind1`复用原ScalarOperations与RHS：
`candidate=a_old+dt*R(a_iter,t_new)`，`a_next=a_iter+omega*(candidate-a_iter)`，omega 为用户指定的有限松弛系数。
HOW显式Step调用与Loop决定顺序；每个物理区间共用一次dt准备、入口OldTime及最终接受，唯一TimeCommit。
普通ForwardEuler Step的重复推进禁令保持；原Stage FE/SSPRK3/RK4不改成固定点方法。

压力是该区间/快照机制的客户端：CellEquationExecution持有U_old和迭代U采样，PressureOperators
借用而不再注入base指针；保留p/phi采样、矩阵/reference、修正、松弛与压力收敛。
原压力在t_old准备/评价，普通标量在t_new评价并借该轮发布的校正phi；这是用户批准的分区关系，
不是完整耦合backwardEuler。RHS、cell/nodal边界、dt/rho响应、canonical COPY/SUM没有新增算法。

固定标量按 HOW 调用次数更新，支持固定次数或用户选择的残差控制。
DiscreteEquation 缺陷为 `a-a_old-dt*R(a,t_new)`，StateChange 采样迭代增量；HOW 选择具体判据及组合。
残差通过原OpRegistry采样，由独立ResidualEvaluation评价，支持既有绝对/相对/归一化设置；max退出默认失败，
显式failOnMaximum:false才接受未收敛退出。恒密度压力默认双判据与普通res共用残差模块，可用While及平面组合判据；保守压力/Eulerian旧信号仍保留。
默认HOW contribution可携带残差和明确scope placement；只在默认声明组合阶段安放，authored HOW不重排。

固定标量限单patch串行、独立named或明确借用的标量owner及既有Central2/Upwind1能力；显式算术面occurrence可读已发布变系数。
Current读取实际视图；OldTime限真实入口快照的局部读取，未实现终点OldTime stencil明确拒绝。
Working保持原已发布Current端口限制，不能把t_old空间准备冒充终点能力；未发布phi/过期区间运行时拒绝。
生产链、native示例、各owner、删除及未验证范围见[iterate.md](../docs/iterate.md)，候选test/iterate/均未运行。
本轮未扩展普通vector fixed-point、压力RK、RAS/IBM/Eulerian、一般质量矩阵或通用隐式求解。

## 首个模型完整普通执行迁移：独立热输运（2026-11-07）

模型contributeThermal只贡献原ddt(T)+div(U,T)=diffusion(D,T)、字段初边值与默认四轴；
普通ScalarMetricUpwindCentral2按实际AST/operator occurrence编译，复用ScalarOperations与共同FE StageGroup。
它保留原metric上风/扩散、集合BC/ILW、显式bounds及前后顺序，未把几何离散换成Cartesian。
热物理数组由普通字段binder唯一持有；不再创建MultiPhaseModel热宿主、ThermalOperations或专用服务端口。
旧ThermalEulerTransport明确退役，无兼容别名；D与对流仅来自WHAT，数值参数不存重复数学。

速度由STATE明确Current/所支持视图读取：原flow借真实thermodynamic U producer，named vector借自己的边界端口，
不要求固定五变量速度身份。没有对流的热扩散无需EOS、压力或phi；独立热不反馈EOS T/rhoE。
模型初边值ScalarFieldInput、编译集合边界选项及无分配ScalarArrayView围绕真实binder/kernel复用，
不新增registry、宿主、clock或Lowering/PlanExecutor模型分支。完整authored HOW不自动插入热调用。
原方法目前仅FE、serial单patch；不扩热RK/固定点/MPI/任意metric系数或cell边界能力。

全部保留功能的源码去向、native/typed界限、删除与构建记录见[model.md](../docs/model.md)，
全模型生产接入清单及待迁移顺序见[models.md](../docs/models.md)。模型注册与可求解分开判断，
JANAF/Sutherland选择身份、局部EOS类、KKT底层或Lagrangian接口不等于完整生产算法。
本轮源码与生产构建完成；test/model/只是待用户指定的输入，未运行测试、未数值验收。

### 公共快照与模型公式入口（2026-11-07）

StageGroup 继续在所有 snapshot 完成后发布 OldTime。公开旧值在本区间不得由 prepare/边界/bounds 改写。
metric FE 保留原准备后更新基值，但使用原 ScalarOperations 实例的独立 updateBase_；old_ 只在区间入口保存。
普通/capacity/cell/fixed-time 的原快照与全 RHS 屏障保持。这个修复可能改变原先误读准备后 OldTime 的消费者，
热自身更新仍用原 BC/ILW/bounds 后基值；未数值验收，不能声称逐位保持。

ScalarMethod 只保存有限 ScalarImplementation，替代 form 和连续 bool 配置；原十个用户方法名、实际算子、
版本/边界与时间限制保留。metric/固定点选项分别解析为现有 typed 值；完整数学仍只经一个普通编译入口。

控制层 ModelLoader::modelDeclarations 提供实际活动模型的 typed contributions，native 元数据与 inspection
组装共用此函数。native 先 decode 启用状态/参数，再收集 contribution.states 的真实签名并解析用户公式。
不按模型字段字符串补注册，不创建物理数组。原 extend 可在模型方程 RHS 添加局部源，最终 AST 仍由原编译链执行；
重复 equation/STATE、冲突签名和缺资源继续拒绝。完整 authored HOW 不自动加模型调用，replace/disable 修改仍未实现。
模型贡献尚未声明的 phase/插件符号、任意版本/MPI 不因签名入口收口而自动支持。
详见 [snapshot.md](../docs/snapshot.md)：生产构建完成，未运行配置、测试或 Explain；下一轮只评估而未实施完整 RAS 迁移。

## 9. 后续双数学入口

未来 Equation Input 与 Variational Input 汇入同一 WHAT。独立变分前端处理 Action、FreeEnergy、Dissipation、Constraint 及各自规则，生成方程/关系，不重新建立求解器。约束语义复用已有 WHAT/STATE。

KKT 是一类约束数值求解方法；离散约束投影不等同 Hamilton 时间作用量。变分前端不自动生成所有 VOF/Level-set/湍流算法。该方向待现有模块接入与四轴执行能力收口后推进，目前不是已实现能力。

## 10. 验收与历史

默认不运行测试、回归、诊断或架构检查脚本，由用户指定算例。构建与数值验收分开；构建成功不证明方程组合正确。授权的实验全部放 test/，逐例报告实际数学、配置、时间、边界、并行和输出。不得覆盖历史基线。

下一步工作见 [problem20261020.md](../docs/problem20261020.md) 和 [provider-plan.md](../docs/provider-plan.md)。历史实现与验证范围见 [lowering.md](../docs/lowering.md)、[scopes.md](../docs/scopes.md)、[stages.md](../docs/stages.md)、[builtin.md](../docs/builtin.md)；[templates.md](../docs/templates.md) 是已退役外部模板路线的历史记录。

独立闭合调用占据原 Ready 后、普通局部关系前的位置；最终刷新占据全部 physical publication 后的位置。它通过既有 stageLocal 计划语义执行，不寄生于任何输运参与者，不增加 snapshot/RK advance/TimeCommit。旧 SSTStageTransport/provider 与输运自动附加闭合已退役；详见 prepare.md。

## Provider 退役第一阶段：conservative 输入

FlowOperations/ConservativeRHS 不再保存或读取 SolverConfig。已有 CompiledNumericalSystem.conservative 冻结 IBM stencil policy、ILW order、界面通量策略、分子兼容参数及来源，thermophysical/numerics 合法投影只在组合入口解析一次。对流 gamma 来自绑定 FluidStateModel 的能力；T、mu、导热率继续来自模型/缓存与已选 transport 资源，不冻结动态场值。单/多 patch 的模型绑定、halo 与界面策略均消费该编译结果。time/dt 继续使用已有 compiled recipe/policy，不另加配置权威。

守恒边界定义由原 Applicator/Pipeline 持有；恒密度 U/p 使用已编译字段需求绑定的 cell 边界端口，PressureNumericalConfig 不再复制边界值，压力校正消费同一端口。FlowOperations 只额外保留未迁移 conservative corrector 的参数与 gamma；SingleFluidStepper 仅在非恒密度压力路径提取这些 legacy 输入，不存整份配置。对应 conservative pressure/KKT binder 仍使用窄参数，没有修改其数学。输入消费矩阵、构建与未验收边界见 [authority.md](../docs/authority.md)。

## Provider 退役第二阶段：资源装配与操作绑定

conservative 的应用入口不再查询 flow.conservative：由冻结 compiledProgram.stateViews 及 producer 依赖查询所需 conservative backing/分量，借用 Field 原 packed Q。withSingleFieldResources 构造或借用 EOS、模型、边界和字段端口，并在 owner 存活期交给公共 executeEquations；它不绑定操作或拥有时间驱动。独立普通字段也交给同一公共入口。单/多 patch 的重复 runConservative 入口已删除，拓扑特有资源准备保留。

bindSingleFluidOperations 仅在计划请求 FlowStepPrepare/PressurePrepare 时创建对应操作 workspace。FlowOperations 保留实际请求的非时间组操作，snapshot/prepare/RHS/advance/publication 由同一 StageGroup 绑定；不再按 provider 家族集合取得操作。数值实现身份仍随冻结 binding 校验，PlanExecutor 不重新选择算法。新增已有资源能力的普通方程只贡献数学与四轴选择，不修改应用家族判断。

冲突按实际操作与 scope 检查：同一 StageLoop 的不同 temporal owner 不可重叠写同一 storage slice，不可重复 RHS/advance 或发布同一 stage 资源；顺序数学更新仍由原别名/版本/失效契约判断，不一律禁止多写者。packed 路径的物理 clock commit 不可重复或位于重复/stage scope。完整操作覆盖、OpRegistry 重复绑定拒绝及计划 owner 核验保留。

恒密度压力的 U/p、rhoConst、预测/修正及持久面通量已由实际方法/资源需求接入同一应用装配边界，详见下节。Eulerian 的 phase/old-time/shared-pressure 资源也已按冻结契约装配，数值 workspace 仍为独立实现，见本节末尾。压力资源与 packed/phase 的布局、发布、通量和 clock 组合能力仍未实现时明确拒绝；不对 conservative＋普通方程/独立闭合施加家族互斥。multi-patch 普通字段/闭合资源端口仍未实现。

已有独立 SST 闭合、四条局部关系、NearestInterior 延拓和 Stage/Current 最终刷新顺序不变；本轮只完成源码装配并构建，不运行测试或数值验收。实际资源矩阵、剩余分派与增量见 [assembly.md](../docs/assembly.md)。后续继续迁移其他装配边界，再逐方法迁出 provider 编译/绑定，最终删除 IProvider/registry/resolver/catalog；不维持长期新旧双轨。

## Provider 退役第三阶段：恒密度单流体压力资源

requiredPressureInputs 查询实际 PressureCorrection、已编译 cell response 和 face prediction，返回已有冻结 PressureNumericalConfig 的引用；核验 U[3]/velocity、独立 p[1]/pressure、正 rhoConst、deltaP Correction 与 canonical pressureFaceFlux views，以及完整预测/修正/发布/clock 操作组。不因 U/p 同名字段、PressurePrepare 单个 OpId 或未消费的数值默认值而启用算法。该查询放在已有 built-in numerical selection 实现中，不新增 IR、selector 或资源中心。

executeSingle、withSingleFieldResources、Environment 和 inspection 不再按 flow.pressure-operators 选择压力存储或布局。网格先提供几何，单场装配一次建立三分量速度 backing 与独立 p；原先随后会被 resize 清除的五变量/PerfectGas 初始化不再执行。Field 的 resize 原本即清零全部分量；原速度初值、p 默认填充值和逐 set 初值覆盖仍在原位置进行。默认 p 填充值沿用已解析参考值，运行中的参考 cell/value 是另一职责，不合并为初值或数学约束的第二份来源。

预测 Working、candidate 与响应由显式引用的 CellEquationExecution 实例持有，PressureOperators 仅借用；Correction、矩阵/线性求解、p/phi迭代采样及压力判据仍归压力数值对象；U旧时间与迭代快照归普通方程实例；持久 corrected face flux 已由独立 FaceFluxExecution 持有并借给压力客户端；恒密度路径不再创建无消费者的 conservative Applicator/Pipeline，也不读取 conservative ILW 数值输入。Environment 的单 source zone/单 patch-per-rank 限制仍消费实际资源契约，empty-dimension 几何标签来自应用边界输入；实际数值边界只从已绑定字段端口消费。未新增多 patch/MPI 能力。

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

结构完成，随后两组明确获授权的串行闭环及限定指标/trace 对照已完成，见 loop.md；不代表任意组合或 MPI 验收。实际链为 native typed WHAT/HOW/STATE/WHICH → Lowering → 已选数值函数 → 编译资源 → 原 binder/PlanExecutor → 原 RHS。`LocalAffine`、`ConservativeResidual` 及恒密度压力方法已退出 IProvider/ProviderRegistry；按 WHICH 的精确方法身份直接编译，不从模型名称选择实现。已迁移操作按方法/过程实际声明的 operation/owner 冻结，不依赖旧 ProviderCatalog 的模型准入。未迁移消费者才走旧 registry/catalog 验证；不是编译失败后回退。

HOW 的 recipe scope 依据所选实现返回的 temporalResidual 分类，不能仅因 WHAT 有 ddt 就把固定时间压力 predictor 包进 RK。共同阶段仍由原 recipe、一次融合 RHS、Ready/RhsReady/advance/Close 屏障生成，普通关系按明确 Stage/Step 位置求值；相邻融合成员限制不变。统一时间计划在编译后与实际 StageGroup 绑定时均校验。

`component(vector,index)` 与三标量 `vector(x,y,z)` 使用同一 Math 签名和 AST 求值。LocalAffine 可借用真实三分量输入，以独立标量 workspace 为输出。Current reader 借用原物理数组与 generation；materialized workspace 在本步显式发布后可读，Current 依赖写入即失效，不能隐式刷新。Stage/OldTime 继续读取真实阶段/q0，并受原 epoch 约束；没有用 Current 数组代替它们。halo、命名 snapshot、perIteration 与未经实现的 Frozen/Lagged 端口仍明确拒绝。

`source.field` 编译 RHS 的唯一 source(expr)，冻结输入版本和纯局部计算，并在启动时绑定 reader。密度动量使用力密度、能量使用功率密度，直接累加原 Residual.source；恒密度动量使用加速度，保留原 dt 权重，不补除 rho。HOW 负责显式计算其输入 workspace；源读取不触发模型 correct/prepare。新增同数值形式的关系改名不要求新增 provider 或 application 分支。局部源内部算子当前不可独立替换，显式选择未实现的内部实现会被拒绝。

密度路径仍固定五变量 rho/rhoU/rhoE 与原可用 EOS/边界/空间方法；恒密度压力仍固定 U/p/rho/phi/deltaP、PrimitiveUpwind1/Central2、rhoConst 和持久 Rhie-Chow 通量。压力动量按项消费；压力修正/速度/压力/通量/松弛等剩余固定 kernel 以结构及扩展规则检查保护，不将新公式只打印后忽略。源项可重复显式相加；重复对流/扩散在融合方法尚无累加实现时拒绝。

两条真实数据流、准备输入和未验证范围见 [loop.md](../docs/loop.md)。普通关系/反馈绑定当前限串行单 patch；密度 field source 需要实际 synchronous stage group，旧 unsplit fused executor 未提供其 reader 端口时拒绝；不宣称压力可耦合任意时间输运。once-step RAS、Eulerian、IBM、Level-set、mixture、普通显式输运/线性块等仍有旧 provider 消费者，后续按实际数值方法迁移。旧四条 SST 关系、独立非局部闭合、NearestInterior 延拓和最终重求顺序没有重做。生产构建与数值结果分别记录。获授权的 uniform RK4/source 与 PISO/damping 四配置已到1 s/10 s，解析/预算/过程及有限通量一致性指标通过；历史压力 MPI 阈值与 Level-set 参考问题仍保留。


授权验收补齐了 native reader 的一个真实缺口：STATE 已显式选择或读取的 builtin unknown，其签名现在可进入同一公式符号表，不要求为 rhoU 等 packed 分量伪造 fields 数组。只借用 BuiltinStateCatalog 的元数据，不创建存储或默认 EOS，不新增名字分支。压力 traceState 为环境开关控制的只读快照，宿主传入原 step/time；不更新数学、clock 或同步。数值证据与范围见 loop.md；四配置之外未运行回归。

## PressureMomentum 声明化与向量数值执行（2026-11-04）

`PressureMomentum` 仅作为内置声明来源；原专用 lowerPressureMomentum、CompiledPressureMomentum 和 assembleMomentum/solveMomentum 已退役。控制层内置压力模板选择 CellExplicit（旧轮名称CellExplicitEuler已退役），并保留显式目标/occurrence 引用；原 ConservativePressureMomentum 的未迁移固定守恒 kernel 继续保护。SIMPLE/PIMPLE 的固定基值/松弛/收敛仍为实际压力方法能力，不能用循环内多次完整 Euler 推进替代。共同 cell 区间准备合并实际 dt 需求；压力 correction/FixedTimeRelaxation 只声明专门准备与固定时间工作量。结果接受由各真实调用需求产生，过程编译器只生成一个 TimeCommit，没有局部主循环。

RhieChow::construct 接受明确速度、压力、梯度、响应和面几何，独立完成原 predict 公式的面遍历；cell-closure 轮把几何绑定、数组 lifetime、canonical COPY 和发布接到独立 FaceFluxExecution，压力只借同一面数组；predict/correct 算法未替换。完整变更、分层组合矩阵、残留核心绑定职责、构建与未数值验收说明见 [equation.md](../docs/equation.md)。本轮未运行测试、Explain 或架构脚本。

## 方程实例与压力客户端（2026-11-05）

普通 Cell 向量经编译字段端口和 `bindSingleFluidOperations → bindCell` 进入每次 occurrence 的 `CellEquationExecution`。`CompiledCellEquation` 存于对应的 `CompiledEquationCall`，空间 lowering 完成后冻结。候选值、Working 和响应工作区依 occurrence 建立，压力 WHICH 用 `response` 引用具体预测调用。独立向量只申请自身需要的资源。

Cell与原scalar共享普通字段资源、读版本与occurrence lookup、Step dt需求、Working变量接受及唯一clock；保留cell几何加权边界和nodal reflection的真实数值差异。梯度按实际项消费并在写目标前完整采样，不再强制作为压力私有potential更新。Standalone RhieChow面资源已在后续cell-closure轮接通；cell Step与recipe Stage任意混用仍拒绝；普通cell/scalar共同Stage及有限scalar固定时间迭代分别接通，后者不进入Stage。SST/IBM/Eulerian等其余固定内核能力不扩大。源码/构建与未数值验收边界、逐项owner和API变化见 [cell.md](../docs/cell.md)，候选在test/cell/。

## 面资源与普通机制收口（2026-11-05）

FaceFluxExecution消费实际CompiledEquationCall.faceFlux及WHICH显式引用的cell响应。独立预测与压力校正
共享同一面资源owner；预测只发布Working phi*，校正过程才可发布Physical phi。PressureOperators不再拥有
面数组/发布代次、面预测叶或重复重构，只保留压力矩阵/reference、校正、p/phi迭代采样及压力专用数学。原边界、响应、
RhieChow公式、canonical COPY及唯一clock的设计不变；mask的metadata COPY移至面资源绑定阶段。

普通scalar/cell/面调用共用EquationInterval，普通Current/Working布局检查通过equationField，Block也共用
occurrence端口唯一性/弱生命周期。空间贡献标记callOccurrence，避免同方程不同调用的选择信息混用；
不重写nodal/cell边界或RK/capacity更新内核。Native scalar/face derived元数据可声明实际数值面输出，
没有初值/边界、不选择算法；geometry-only装配不要求它预先拥有物理数组。压力默认和authored流程
仍进入同一编译/绑定链，完整任意native packed压力资源映射尚未实现。范围、增量、生产构建与未测试说明
见[cell-closure.md](../docs/cell-closure.md)。不代表provider全部退役或新的湍流/IBM/多相组合可运行。
