# A5 后续收口：统一贡献与实际资源绑定

日期：2026-10-23。完成选定模块的基础契约与真实执行闭环；不代表全部物理算法或资源已经通用化。未修改正式算例、未提交/推送GitHub。

## 1. 起点和本轮增量

起点HEAD `761712ef1d4f3f140178c08c24954702ed15bb39`，包含上一轮未提交A5工作。`test/contributionA5/evidence/start-status.txt/start.patch/start-sources.json`保存实际起点。原先6份tracked历史docs改动及全部A5修改保留，不算成本轮增量。

- before：`test/contributionA5/baseline/sonicSolver`，SHA256 `77d83bb057449006367e404ef3ef0d8329504a6b54e23d3172d2a05eb0a8707e`。
- final：`build/sonicSolver`及不可变副本`test/contributionA5/candidate/sonicSolver`，SHA256 `26da6e5bde9ddd94aa03d0148cfb1d163fb76690dacfa8a2f0f6fab3e2fb9ae4`。
- 相对起点快照：**20份C++源码/头变化，+193/-67，净+126行；无新增生产源码文件**。CMake、测试、文档不计入该数字，逐文件见`evidence/metrics.json`。

本轮是composition、编译元数据与具体资源绑定的职责收口，属于“改怎么组合/绑定已有方程”；没有新增物理数学、kernel、时间循环或Registry。测试中的内部模块使用现有算子和scalar实现。

## 2. 完成项与删除项

1. ContributionRecord区分module类别与实例id：默认单实例id=module，显式实例id=module@instance。现役模型不因增加实例元数据获得不存在的数值多实例能力。
2. 显式加性RHS等所有声明注册后、transformation前合并，保留来源/原顺序；**删除按相同AST拒绝重复贡献的逻辑**。不推断系数相加，不从依赖图猜算法。
3. closure新增显式outputs，与原states依赖引用分开；重复闭合输出owner拒绝。注册完成后重新核验方程引用签名，防止后注册STATE绕过已有表达式检查。
4. 默认WHICH和显式选择沿原builder收口：显式优先、同级重复显式拒绝、未被覆盖的默认冲突拒绝，NumericalBinding保留实际来源和覆盖来源。authored HOW仍要求显式非空WHICH。
5. single-fluid模型输运与named scalar均声明typed物理端口。复用ResourceViewRequirement/ResourceView，统一验证需求与实际绑定，不创建ResourceManager。原scalar binder只消费自己负责的需求，模型端口由现有flow装配处绑定；ScalarEquationData为只读配置，ITransportModel为明确可写的更新/发布端口。
6. 不活动模型不构造Manager、不初始化模型数组；活动Manager仍是原数组owner，稳定借用给原CompositeTransportProvider。每个compiled occurrence绑定同一个实际输运实例，禁止重复/缺失/未消费资源。
7. **删除 SolvePlanNode::legacyAdapter、fingerprint字段、打印和旧断言**；保留真正的provider/operation、BlockSolve和执行顺序检查。架构文档中不存在的旧selector说明同步删除。

`LegacyPressureInputs`仍被conservative pressure binder/corrector消费，`DistributedRowMap`仍被Eulerian装配消费，native decode仍是生产IO；未因名字删除它们。合法输入拼写`legacyOnceStep`仍在使用，未假称已经退役。

## 3. 贡献规则和边界

| 情况 | 本轮规则 |
|---|---|
| 方程预制/模型/用户 | 同一SystemContribution、EquationRegistry、STATE和编译链；来源仅解释 |
| 模块/实例 | module是实现类别，id是实例唯一键；重复id拒绝，独立标量实例合法 |
| 显式加项 | 向作者明确提供的equation/namespace目标增加RHS，保留每一次声明；无目标拒绝 |
| 相同表达式 | 不自动去重；新增测试将两个相同`.2*b`实际计算为`.4*b` |
| coefficient/closure | 不默认相加；组合必须在WHAT中明确写出。closure.outputs冲突拒绝，states不冒充输出 |
| STATE/producer | 重复输出声明/符号拒绝；既有producer依赖、版本、scope和只读检查保留 |
| default/explicit | equation+occurrence是选择键；显式覆盖默认，重复显式和未消歧默认报目标/来源 |
| authored HOW | 使用自身的显式非空WHICH，不能暗中继承preset bindings |
| disabled | inspection不产生专属贡献；无执行operation、输运端口请求或模型数组 |
| 缺能力 | 保留具体Unsupported/缺STATE/缺资源等原因，不改模型、recipe或后端补救 |

显式equation级WHICH覆盖该方程默认项；occurrence精确选择保留原匹配优先级。HOW顺序始终来自声明/所选方法展开，没有为了让模块运行而自动拓扑排序。

`finalizeMathematics()`是原builder的组合完成步骤，不是第二份运行时数学IR。直接使用builder追加项的调用者需要在冻结前显式完成此步骤；不保留旧即时扩展转发入口。

closure.outputs只对明确声明的输出提供owner检查。旧的说明性closure若未声明outputs，不能据此宣称其完整读写已被证明。物理数学仍以真实AST/所选方法为准；本轮不实现任意闭合、系数或方程合成。

## 4. ownership、顺序、API和依赖

| 对象 | before → after |
|---|---|
| physical Q、clock、EOS/cache | 原Field/StateBundle → 原owner；无layout或cache复制 |
| named scalar | 原ScalarResources唯一数组 → 同一owner，增加实际typed view验证 |
| turbulence数组 | flow里无条件构造原Manager/按config初始化 → 按compiled端口需求构造原Manager，共享强owner保活；仍只有原数组 |
| 服务 | 原CompositeTransportProvider借用Manager → 同一借用，类型/身份/lifetime需求单独验证 |
| solver workspace | 原RK/Residual/FluxField/pressure workspace → 原owner，不进入Field |
| stage/BC/halo | 原startup初始化、once-step或stage流程 → 同一顺序；validator不插刷新/通信 |
| 数学合并 | 遇到贡献立即扩展 → 声明齐备后按同一贡献顺序扩展，旧输入加法顺序不变 |
| API | contribution实例/closure输出/WHICH来源为现有类型增量；新增builder finalize及资源验证函数，删除死flag，无forwarding overload |
| dependency | core中立类型 → solver组合/方法 → application具体binder；无IO反向依赖solver，无新增反向边 |

公共NumericalCompiler的selected lowering、既有Math/AST和PlanExecutor本轮未增加任何新模块名称分支。控制层仍只选择、传递typed参数和收集贡献，没有运行模型时间循环。

Manager的初始化及原`initializeTransport`仍在原startup位置显式执行；它们不是timestep stage。模型once-step的in-place source/diffusion、边界与发布仍由已有Plan callback控制。没有将once-step k–ε/SST改成RK stages，也没有删去合法Frozen/Lagged读取。

## 5. 真实资源盘点

| 资源 | 数学/编译需求 | 实际owner与binder | 消费/生命周期/并行语义及范围 |
|---|---|---|---|
| packed flow state | builtin WHAT/STATE、CompiledTarget/STATE views | Field、StateBundle；Environment和原single/multi初始化 | flow/pressure kernels；case lifetime；owner COPY，原支持范围 |
| named scalar | scalar method physicalResources: ScalarEquationData | ScalarResources数组/强端口，bindScalarResources | scalar Stage/BC；运行期借用，solver持workspace；当前serial单patch |
| single-fluid turbulence | RAS/LES贡献STATE/WHAT，所选方法请求ITransportModel | 原Manager自有k/epsilon/omega/mu_t，flow具体装配按需求绑定 | once-step或stage/closure；运行期；原WriteOwned/readHalo/publication，MPI能力未扩大 |
| phase state | Eulerian预制/逐相method slices | 原PhaseSystem、Eulerian application binder | phase equations/transport；case lifetime，原canonical state COPY、源/残差SUM；未通用化 |
| face flux | numericalWorkspaces/pressure与共享通量需求 | PatchWorkspace/原pressure face workspace | 阶段内候选→canonical F*→COPY→正负装配；未建空通用binder |
| geometry/IBM | boundary closure/producer/constraint声明 | Field geometry、原IBM对象/adapter | ghost读取或transfer；几何/标签COPY，贡献SUM；仅相应已实现方法 |
| multiplier | constraint WHAT/STATE/block需求 | 原IBM约束数组与现有binder | projection/KKT kernels，唯一owner/COPY；本轮不扩production KKT |
| thermo/closure | STATE producer/所选EOS能力 | 原Field stateModel和cache | Q writer失效、原lazy读取；不另算一份热物性；JANAF等缺口保留 |
| solver workspace | numericalWorkspaces、STATE workspace views | 原PatchWorkspace、ScalarOperations、PressureOperators等 | solver执行生命周期；residual/source SUM；无第二份物理状态 |

实际`validatePhysicalResources`核验端口身份/type/access、活着的owner、遗漏/重复/未消费项；不选择binder、不分配数组、不执行同步。Scalar binder接收其他具体binder已绑定的views，在发布任何自身注册前验证全局需求；其他端口缺失/权限不符时同样回滚局部未发布数组，失败不留下悬空注册。运行时实际STATE别名、边界和所选方法仍有独立检查，不能只相信一份自报descriptor。

本轮只把选定scalar及single-fluid模型端口闭合，不声称phase/face/IBM/thermo的全部资源均已转成通用资源对象。直接C++手工装配仍须提供真实services，原数值binder的独立能力校验继续保留。

## 6. 构建与集中测试

迁移前完整build为ninja no work to do。最终完整build成功，最后生产变更后的23个actions compile/link失败均0；测试增强后另完成incremental build。既有duplicate-library warning保留。架构检查通过：709源文件、2条原有依赖债务、11组同名头；新增守卫防止死flag恢复或模型/scalar端口绕过需求核验。

定向CTest最终**44/44通过**，不是全部项目CTest。新增`moduleContributionClosure`测试随后增强canonical四轴/Plan对照并独立重跑通过；另43项涵盖A1/A2/A3/A4、scalar耦合、物理/时间意图、source、pressure、Eulerian、RAS、SST、nativeIBM等。真实命令在`test/contributionA5/ctest/CTestTestfile.cmake`，stdout在`evidence/ctest-final.log/ctest-contract-final.log`。

中间nativeImmersed测试失败的原因是退役flag时误把结构节点当作operation leaf；已修复测试为仅对真实operation验证frozen provider，并保留BlockSolve、数量与顺序检查。不是IBM数值失败，也未放宽PlanExecutor。

### 新内部模块的实际数值验收

`test/contributionA5/test_contribution.cpp`：两个`internal.diffusion@a/@b`实例声明原scalar STATE和扩散方程；三个`internal.reaction@first/@second/@reverse`只提供显式RHS。模型源项先注册、目标方程后注册，组合完成后可解析真实目标。

```
ddt(a) = diffusion(.03,a) + .2*b + .2*b
ddt(b) = diffusion(.08,b) + .3*a
```

无量纲a/b，D量纲L²/T、reaction系数T⁻¹。9×1×1节点、NG1，x∈[0,1]，y/z不活动；初值a=1/b=2，x两端zeroGradient；无flow Q、EOS、IBM、MPI。显式HOW固定a→b→Commit；RK4/maxDt.001/start0/end.01，实际10步。通过typed contribution→原compiler→scalar资源binder→原Plan/数值kernel执行，不新建provider或CFD solver。

均匀解扩散严格为零但D真实存在。独立解析解使用w=sqrt(.12)：a=cosh(wt)+.8/w*sinh(wt)，b=2cosh(wt)+.3/w*sinh(wt)，误差<1e-13。两个相同`.2*b`都参与计算，没有去重。

只改变两个独立方程模块的注册顺序，不改HOW；canonical WHAT、按id比较的STATE、WHICH和Plan均一致，最后逐元素相等。Provenance/注册显示顺序不要求逐字相同。本测试不叠加冲突湍流闭合。

集中负例验证重复实例、symbol/producer输出、缺目标、默认歧义、显式冲突/覆盖、未知方法、闭合owner冲突、错误STATE形状及资源缺失/重复/owner/lifetime；A3/A4的producer依赖、失效、错误版本、scope和未publication负例本轮重跑通过。没有把这些断言通过算作CFD时域验收。

## 7. Native代表性控制链：8组实际前后运行

路径`test/contributionA5/{before,after}/<name>`。来源为上一轮现役stageSST的独立副本，保留native单目录/注册身份；改为以下明确的均匀初值、原once-step RAS和短时控制。不是修改正式case。`cases.py`保存生成规则，`evidence/cases-before.json/cases-after.json`保存binary SHA和命令，逐例有explain.log/run.log/result。

共同WHAT/STATE：NS的rho/rhoU/rhoE、原EOS closure；启用k–ε时新增k/epsilon及mu_t闭合，启用SST时为k/omega。PG gamma1.4/R287.05、hConst、const mu1e-5 Pa·s/Pr.72。域[0,100]×[0,20]×[0,.1]，100×20×1cells/101×21×2nodes、NG4、scale1/SI，serial单patch，无物体/IBM/ILW（wallClosure及热壁面不适用）。

IC rho1.2、U(.3,0,0)、p100000，k=.01、epsilon=.001、omega=1（未活动变量的输入仍保留，但不注册模型数组/执行PDE）。flow nativeRight/nativeSides zeroGradient、frontBack empty，nativeLeft未显式声明；k/omega在nativeLeft也有zeroGradient，epsilon只声明Right/Sides zeroGradient及z empty，Left未显式声明。这里保留实际支持的原BC行为，不把此短时控制回归说成标准近壁k–ε物理验证。

WHICH：WENO5/Rusanov、characteristic conservative shared flux、Central2黏性；没有IBM近壁重构。RAS输入明确为`transportMethod: legacyOnceStep`：沿用原source/diffusion in-place sweep与closure/floor数学，**无k–ε对流、无新linearization，不能冒称完整标准k–ε输运**。SSTOnce同样不是stageConservative。

共同start0/end.001、无endStep、CFL.5/maxDt.0001、初始及每步输出；实际每例10步/.001、11帧。所有比较使用同一输入的before/final二进制，不把不同物理选择比成相同。

| case | 控制选择与实际贡献 | 时间/验证 |
|---|---|---|
| absent | models registry不注册turbulence；gravity显式false，无RAS | Euler；验证未加载模型不产生操作/端口 |
| laminar | 显式laminar；无RAS、无gravity | Euler；与absent输出一致 |
| gravity | laminar；g=(.01,0,0)开启，原momentum source及融合保守kernel的mechanical-work行为 | Euler；中心Ux=.30000999999999983，相对off增量1e-5符合g*end |
| kepsilon | block存在未写enabled，默认启用kEpsilon及原k/epsilon/mu_t | Euler；实际模型推进与输出，不是注册名字 |
| kepsilonGravity | 同时启用上述kEpsilon和gravity | Euler；现有相容范围内组合执行 |
| disabled | 保留kEpsilon配置但enabled:false；gravity false | Euler；没有turbulence操作、端口或数组，结果与laminar一致 |
| override | kEpsilon；显式WHICH分别选择TurbulenceTransport、inputs[k,epsilon]覆盖默认，显式time为RK4 | RK4只推进flow，RAS仍once-step；Explain保留覆盖来源，未偷偷换成stage SST |
| sstOnce | kOmegaSST，原once-step，非k–ε替代验收 | RK4 flow；仅补充模型端口保护 |

8组同配置VTS/VTM、dt/Plan trace和全部统计一致；absent/laminar/disabled也彼此一致。source开关改变物理输入，仅检验预算/贡献，不要求gravity等于off。中心速度预算通过不等于整个开边界域是封闭守恒系统。

`verification.json`逐例保存物理体积加权质量、动量、能量、rho/p及湍流场min/max、L2；用实际节点坐标的梯形权重，不使用裸nodal sum。例：laminar质量240、x动量72；gravity末态质量239.99999999783532、x动量72.00239982034522。开放/未显式指定边界意味着不能把全域动量变化直接等同无通量预算。所有这些量在同配置迁移前后相等。

## 8. 原执行链保护：5组实际前后运行

路径`test/contributionA5/protection/{before,after}/<name>`，真实命令/binary见`evidence/protection-before.json/protection-after.json`。每个case都有run.log/result，最终均与before字段及clock/trace一致。

### 8.1 flowPair / scalarPair

来源上一轮`test/resourceA5/composition/<name>`，输入不改。

- 域[0,1]×[0,.25]²，33×5×5nodes/32×4×4cells，serial。flowPair NG3、scalarPair NG1；后者geometry-only，Q=0/EOS=none。
- 双标量数学：ddt(dye)=diffusion(.03,dye)+.7*salt；ddt(salt)=diffusion(.08,salt)+.7*dye。scalar无量纲，D为L²/T，.7为T⁻¹。
- dye默认0、x≤.5初值1；salt默认.5、x>.5初值2。dye Left固定1/Right固定0，其余四面zeroGradient；salt六面zeroGradient。
- flowPair另有NS Q，PG1.4/R287.05、hConst/const mu0/Pr.72；rho1、U(.1,0,0)、p1，六面rho/U/p zeroGradient；WENO5/StegerWarming/characteristic/shared flux，无viscous operator、无模型/IBM/ILW。两标量用Central2；HOW为原flow三条调用后双标量再统一Commit。
- RK4，start0/end.005/maxDt.0005/CFL.5，无endStep；初始/每步输出，实际10步/.005、11帧。保护真实交叉Stage读取和同一binder，不扩散/反应数学。

### 8.2 stageSST

来源`test/resourceA5/after/stageSST`；域/网格、PG和均匀IC同第7节，但WHICH为RK4/WENO7/Rusanov/Central2、stageConservative k/omega、共享质量通量Upwind1。flow BC同第7节；k/omega Left/Right/Sides zeroGradient、z empty，无IBM。start0/end.001/endStep1/CFL.5，无显式maxDt；初始/每步输出，实际1步7.316037672350e-4、2帧。没有把这次首步保护称作长期SST验证。

### 8.3 eulerianSST

来源`test/resourceA5/after/eulerianSST`。域[0,1]×[0,.2]×[0,.01]，20×4×1cells/21×5×2nodes，serial；原逐相mass/momentum/enthalpy、shared p、相质量加权SST。alpha .9/.1、U0/T300/p101325/k.01/omega1。水perfectLiquid rho998.2/mu.001002/cp=cv4182/kappa.598；空气PG输入rho1.1766/R287.05/mu1.846e-5/cp1007/cv718/kappa.0262。

Left固定上述相量，Right zeroGradient，top/bottom U0其余zeroGradient，z empty；k/omega Left固定.01/1其余zeroGradient，walls bottom/top。MRF中心(.5,.1,0)、z轴omega1；SchillerNaumann d.001、virtualMass.5、lift0、RanzMarshall，无Lee/RPI/IBM。

Euler intent/PIMPLE3/2/nonOrth1即2passes、U/p relax .7/.3，原phaseUpwind与turbCentral2/Upwind1；HYPRE FlexGMRES/BoomerAMG max100/rel1e-8/abs1e-12/Krylov30/refresh10。start0/end.0002/endStep1/CFL.5/maxDt.0001，实际1步3.624299498372e-5、2帧。它保护既有phase owner，未宣称将phase资源通用化。

### 8.4 no-slip SST Ghost

来源正式`test/IBM/cylinderFlowSSTGhost`的独立副本，网格、壁面、物性和终止条件均不改。实际域[0,100]×[0,20]×[0,.1]，100×20×1cells/101×21×2nodes、NG3，serial单patch、SI scale1。STL截面中心约(20,10)、半径约5，穿过薄z域；实际bounds保存在`evidence/ghost-geometry.json`。

NS Q+原once-step SST，PG gamma1.4/R287.05、hConst/const mu1.7894e-5 Pa·s/Pr.72。IC rho1.225、U(2.36643,0,0)、p101325、k1e-6、omega1；epsilon .0005输入保留但无epsilon PDE。Left固定这些量，Right/sideWalls zeroGradient、frontBack empty。

**Ghost stationaryNoSlipAdiabatic，静止无滑移、绝热，Euler ILW关闭（0）**；原quadratic ghost/viscous及SST壁距/壁面关系不变。RK4/WENO5/LaxFriedrichs/characteristic/shared flux+Central2；mu_t每步冻结，读取原step-entry boundary-ready速度，不改成stage SST。

start0/end.02/CFL.3/maxDt1e-5，无endStep，runTime每.01输出；before和final实际均2000步/.02、3帧。全部payload、clock/trace与原结果一致。因本轮触及输运端口而运行对应SST Ghost，**没有改IBM算法**；它不替代forcing/DLM/KKT验收，也不是长期绕流/湍流精度证明。总量统计包含VTS中的原nodal mask/solid状态，不声称cut-cell精确守恒积分。

## 9. 生产输入负例与错误分类

`test/contributionA5/negative/<name>`从第7节override独立复制，除下述单一因素外，WHAT/STATE、网格、IC/BC、recipe和终止设置不变。均实际调用candidate的`explain`，exit1、无物理步或VTS输出；命令与完整原因在`evidence/negative.json`。

| case | 改动 / 实际阻断链 |
|---|---|
| stageKepsilon | transportMethod改为stageConservative；所选实现只支持single-fluid SST，准确Unsupported，不换模型 |
| initial | 删除epsilon初值；required native STATE initialization缺失 |
| whichConflict | 增加同key显式k绑定；Conflicting explicit WHICH，指出目标和文件 |
| unknownMethod | 显式选择MissingNumericalImplementation；没有registered implementation，无fallback |
| modelConflict | 第二个turbulence对象、独立文件；typed explicit target Turbulence冲突，列出两个来源 |

数学签名错误在Math/STATE入口拒绝；组合/默认/owner冲突在composition拒绝；所选方法能力与后端不足由原provider/host检查；缺实际资源由binder验证。CLI外层仍可能显示CONFIGURATION INVALID，但内部保留具体原因，没有统一改写成“公式不合法”。

没有新建展示用方程表。Explain使用原贡献、修改、STATE及NumericalSelection，显示模块/实例、加项来源、闭合输出、显式覆盖、compiled端口和实现限制；运行成功后的Bound physical resource另记日志。不能仅凭Explain宣称端口就绪。

输入准备阶段曾有未活动omega请求输出和重复which注册导致拒绝，保留原日志；修正的是实验副本的注册/输出控制，没有借此修改正式输入或数值算法。最终8组模型及5组保护使用冻结完整配置。

## 10. host MPI与原失败保留

### 10.1 原PISO serial/MPI阈值

路径`test/contributionA5/pressure/{before,after}/{original-serial,original-MPI2}`，来自现役PISO输入，不改变容差。

U momentum+div(U)=0/p multiplier，rhoConst1、mu.01/Pr.72。域[-.5,.5]²×[0,.01]，80×80×1cells/81×81×2nodes、NG3。IC U0、p101300、输入rho1；top U(1,0,0)，其余三边U0，四边p/rho zeroGradient，z empty。Euler intent/PISO outer1 pressure2 nonOrth0、relax1；primitiveUpwind1/Central2/RhieChow，HYPRE FlexGMRES/BoomerAMG max1000/Krylov50、rel1e-10/abs1e-12、refCell0/p101325。

serial单patch；MPI2split[2,1,1]。start0/end1/endStep2/CFL.1/maxDt.001、初始/每2步输出，实际2步/.00078125。before/final四次均exit0，首矩阵/RHS按坐标相同，同配置payload保持。serial/MPI maxU差5.548720974878074e-10，压力shape差**1.0851072147488594e-6 > 原1e-6阈值，仍Fail**。本轮未调容差、未替换输入，之前tight诊断的通过不拿来代替这项失败。

### 10.2 当前D1 level-set publication

路径`test/contributionA5/levelSetMpi/{before,after}/{serial,mpi2x}`，复制D1的planeAdvection，输入不改。NS Q+phi advection/normal/curvature/material；域[0,100]×[0,20]×[0,.1]、100×20×1cells/101×21×2nodes、NG4；PG1.4/R287.05/hConst、mu.001/Pr.72；material water/air rho1/.125、mu.001/1.48e-5。IC rho1、p1e5、U(.2,.1,0)，phi为过(50,10,0)、normal(1,.37,0)的距离；flow/phi Left/Right/sideWalls zeroGradient、frontBack empty。k/epsilon配置保留但不活动，无IBM/RAS/表面张力/质量修正。

SSPRK3/WENO5/StegerWarming/Central2，HJ advection5/epsilon1e-6/power2，interfaceThickness1.5/narrowBand.02、reinitSteps0（pseudoDt.01未执行）；serial1patch/MPI2split[2,1,1]。start0/end.001/endStep3/CFL.3/maxDt1/3000，初始/每步输出；实际3步/.001、4帧。四次exit0、同配置全部payload保持；本轮没有重新定义Level-set算法或重初始化。

旧`test/interface/nativeLevelSet/frozen-kernel.json`未改、其历史失配仍保留。该参考自身标明不是可运行production，缺原始场；本轮没有重跑旧isolated哈希工具或重新冻结后宣称解决。当前D1保护与历史参考首差异是两件事。

### 10.3 并行基础保护

实际host执行`mpirun -np 4 build/test_distributedConstraintLayout`、`-np 4 build/test_sparseCanonicalCopy`、`-np 2 build/test_distributedKKT`，均exit0。它们保护唯一ownership、COPY和最小代数/HYPRE路径，不是production IBM KKT CFD。日志和真实命令/binary见`evidence/mpi-before.json/mpi-after.json/pressure-before.json/pressure-after.json`。

无sandbox launch失败；运行时AppleDouble退出清理warning若存在保留在原日志，不能当作数值assert失败。未改变raw MPI、communicator、GlobalDof SUM、canonical单F* COPY及正负残差装配。

## 11. 证据与复现

所有测试、脚本、日志和checkpoint在`test/contributionA5/`，数值输出各case/result，未使用系统临时测试目录。构建/测试工具TMPDIR指向本目录evidence。正式case、旧参考和原始长时输出未覆盖。

主要命令（从仓库根目录；before读取已冻结binary，after读取不可变candidate）：

```
DEVELOPER_DIR=/Library/Developer/CommandLineTools TMPDIR="$PWD/test/contributionA5/evidence" cmake --build build --parallel 4
python3 -B tools/check_architecture.py
python3 -B test/contributionA5/cases.py before
python3 -B test/contributionA5/cases.py after
python3 -B test/contributionA5/protect.py before
python3 -B test/contributionA5/protect.py after
python3 -B test/contributionA5/run_tests.py
python3 -B test/contributionA5/negative.py
python3 -B test/contributionA5/pressure.py before
python3 -B test/contributionA5/pressure.py after
python3 -B test/contributionA5/mpi.py before
python3 -B test/contributionA5/mpi.py after
python3 -B test/contributionA5/verify.py
```

MPI需host允许PRTE socket。脚本依赖本地保留的冻结源case和binary，不能假称只克隆Git即可重现全部历史证据；新增C++集中测试不依赖这些CFD输出目录。

最终索引：`build-atomic.log/build-final.log`、`architecture-final.log`、`ctest-final.log/ctest-contract-final.log`、`cases-before-validated.log/cases-after-validated.log`、`protect-before.log/protect-final.log`、`verification.json`、`negative.json`及MPI manifests。字段对照、weighted统计及每例实际终点均在verification.json。

## 12. 未完成项与后续顺序

**Implemented**：统一typed贡献的实例/来源、显式加项、闭合输出owner、默认/显式WHICH冲突规则；selected scalar/模型端口实际绑定与lifetime核验；native k–ε及gravity启停/组合/覆盖；内部多实例双标量解析验收；死flag退役；真实前后CFD保护。

**保留的真实实现**：专用融合flow/源项、原k–ε/SST内核与既有floors、Manager自有数组、Eulerian phase binder、pressure snapshots/HYPRE、IBM几何和transfer、现役native decode。原gravity保守source kernel包含mechanical-work行为，本轮未拆融合源项；这不是任意source完整read/write-set数学证明。

**未完成**：全部phase/face/MPI/constraint资源的通用需求驱动装配；legacy pressure输入和row-map消费者迁移；任意closure、producer/layout/MPI runtime证明；完整k–ε对流/壁面/线性化；不支持的RAS/pressure/IBM组合；production KKT；mixture全面迁移；通用Newton/DAE/Subcycle/外部公式/变分；原压力阈值和历史Level-set首差异。

建议按真实资源消费者逐类替换剩余binder并同轮删除旧入口；再按明确数学/参考扩展数值能力。压力精度策略和历史参考恢复独立处理。GUI继续消费同一生产IO/inspection，未实施C2功能。

**基础契约的证据**是新增内部模块仅提交自身typed贡献与已有资源需求即可通过同一控制/编译/binder/Plan执行，公共NumericalCompiler和PlanExecutor无需增加该模块分支；**不是所有物理算法、资源和任意组合已经完成**。
