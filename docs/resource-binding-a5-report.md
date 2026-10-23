# A5：所选数值实现边界与首个声明驱动资源组合

日期：2026-10-22。本轮完成选定资源类型的闭环，不宣称所有资源已自动装配，不实现新 PDE、隐式算法或模块物理耦合。不修改用户正式算例、不发布 GitHub。

## 1. 起点与实际增量

HEAD/分支：`761712ef1d4f3f140178c08c24954702ed15bb39` / `codex/core-single-path-closure-20261022`。起点已有6份 tracked docs 修改全部保留，未算成本轮修改。`test/resourceA5/evidence/start-status.txt/start.patch/start-sources.json` 记录起点。

- before binary：`test/resourceA5/baseline/sonicSolver`，SHA256 `9ca256fc0ef6df62b006c4aed92b99e50b6b02b3f0b5031a6e5c48e7e3b772af`。
- final binary：`build/sonicSolver`，SHA256 `77d83bb057449006367e404ef3ef0d8329504a6b54e23d3172d2a05eb0a8707e`；同一工作树未提交状态，不能只用 commit 标识。
- 与起点源码快照相比：20份C++源码/头发生变化，其中新增3份；**+448/-322，净+126行**。不含CMake、脚本、文档。详见 `evidence/src-metrics.json`。

## 2. 层次、数学类别与所有权

本轮属于“改已有 equation 的编译/绑定方式”，不增加物理 term、方程或数学能力。

| 项目 | before → after |
|---|---|
| 数学/AST/签名 | 原 Math/Typed WHAT → 不变；无新 MathAnalyzer |
| 融合空间 lowering | 公共 NumericalCompiler 遍历/限制AST → 所选 IProvider::lowerSpatial；共享 FusedSpatial 实现 |
| term 匹配 | Catalog 同时含融合操作数知识 → Catalog 负责选择，具体融合 matcher/注册归 FusedSpatial |
| 标量物理数组 | scalar-only专用执行装配，flow入口无命名数组 → 同一 bindScalarResources 消费编译端口需求 |
| 端口生命周期 | 专用入口局部数组/强端口 → ScalarResources 明确拥有数组/强端口，services借用 |
| 物理 Q/EOS/clock | 原 Field/StateBundle → 同一owner，未复制或改layout |
| RK、Residual、FluxField、压力workspace | 原 solver → 同一owner，不回塞Field |
| 执行顺序 | 原 CompiledSolvePlan → 不变；binder不运行stage/刷新/同步/求解 |
| 依赖 | application解析追加方程时重声明所有主变量 → 复用solver BuiltinStateCatalog；CaseAdapter显式链接现有solverSystem |
| MPI | 原 owner COPY / contributor SUM → 不变；不新增raw MPI |

调用链：native typed选择 → 原contributions/EquationRegistry/STATE → 所选数值provider编译 → CompiledEquationCall.physicalResources → 具体资源binder → 原distributed views/稳定services → 原Plan/OpRegistry/kernels。

没有第二份STATE/clock、全局资源控制器、新Registry、solver family或fallback。3个新生产文件是共用融合lowering的`.h/.cpp`及既有scalar资源实现的公开头，均有真实消费者。

## 3. 数值编译器边界收口

新增 `IProvider::lowerSpatial`，缺实现时明确拒绝；公共compile显式接收选定ProviderRegistry，所有调用者直接迁移，没有旧重载转发。

`SF_fusedSpatial` 被 ConservativeResidual、primitive pressure、conservative pressure 三类具体实现复用。原符号/系数、项左右侧、裸operator操作数、ddt(Output)、source kernel/EOS/flux结构要求原样保留。TermProviderCatalog不再自己拥有这些融合操作数规则。

公共编译仍负责统一recipe消费、显式binding是否被使用、halo/workspace汇总和计划组织。普通scalar继续使用已有FormulaCompiler的局部表达式lowering，没有强行通过融合结构检查。

测试专用provider收到带Negate的原AST，证明公共compiler调用所选lowering；它不是新的生产signed-flow数值实现。原融合signed/coefficient负例、未知flux、重复ddt和无用recipe负例继续通过。

## 4. 资源绑定闭环及删除项

普通scalar method在编译结果中通过现有 `ResourceViewRequirement` 声明只读 `ScalarEquationData`：owner、type、occurrence；target和storage仍引用原CompiledTarget/STATE，不复制配置authority。

`bindScalarResources`：

1. 核验具体owner/type/access/occurrence及对应CompiledScalarEquation。
2. 仅分配实际consumer要求的NamedDistributed标量数组；引用已有field声明/初值/集合，不按字段名特判。
3. 绑定已有fixedValue/zeroGradient边界；缺失、重叠或非完整轴对齐face明确拒绝。known-source callback仍须真实实现，native binder不造空回调。
4. 同一storage只有一个初始化owner；与已有owner冲突拒绝。同一symbol多次消费共享数组。
5. 全部数组/端口核验成功后才注册distributed views。失败不留下悬空STATE注册。
6. 局部ScalarResources持有数组和immutable强端口，覆盖整个运行；services只存weak端口。stage快照/RHS仍由原solver创建。

独立多标量及flow+多标量复用此函数。**删除 `executeScalarEquations` 的声明/实现与中央 `ScalarOps::Provider` 装配分支**；geometry-only尾部按已绑定资源运行唯一executeEquations。flow保留原资源/时间执行链，仅接入声明数组/端口和原STATE输出视图；single/multi helper调用直接迁移，无兼容overload。

native added-equation旧逻辑把rho/rhoU/rhoE也注册为NamedDistributed，导致组合输入即使合法仍无法装配。现在已有builtin由现有catalog负责，新增named字段保留其原注册检查；引用“数学符号存在”仍不意味着资源已满足，后续编译/binder继续核验。

Explain显示的是**编译资源需求**，实际binder成功后输出 `Bound physical resource`。不能仅凭Explain的模型名称声称数组已绑定。

本轮没有统一Eulerian/flow的全部物性、phase/face/IBM等资源；这些旧的真实binder继续有消费者。完整A5仍需逐类型迁移。

## 5. 构建、架构及集中测试

迁移前完整build为 `ninja: no work to do`。最终完整build成功，48个actions，compile/link failures均0；既有duplicate-library warning保留。开发期间的缺include/链接依赖已修正，旧失败日志保留，没有改kernel绕过。

架构检查：**709份源文件、2条原有反向依赖、11组同名头**，通过。新增守卫要求公共compiler调用selected lowering、禁止恢复融合AST限制和scalar家族装配分流。

最终定向CTest **43/43通过**，52.35秒；并非整个项目全部CTest。包含原42项Math/WHAT、A2/A3/A4、真实耦合scalar、fluid/scalar、shared flux、stage/once-step SST、Eulerian、pressure/source、nativeIBM/viscous wall等保护，及新增compiledResourceComposition。选择和真实命令在 `evidence/test-selection.json` / `ctest/CTestTestfile.cmake`，最终日志 `evidence/ctest-final.log`。

新增集中测试设置（独立于CFD）：

- `test/resourceA5/test_resources.cpp`。两个无量纲标量a/b，`ddt(a)=diffusion(.03,a)`、`ddt(b)=diffusion(.08,b)`；9×1×1节点、NG1、x∈[0,1]，没有流体Q/EOS/IBM/MPI。初值1/.5，x两端fixedValue0；其他方向不活动。RK4、maxDt.0005、start0/end.005，实际10步，未输出CFD云图。
- application binder与独立手工ScalarField/端口装配使用同一冻结系统执行，逐元素、dt/time/step完全一致。对照方没有调用binder。
- 验证错误owner/type/access/occurrence、重复需求/已有owner、缺field/初值/边界、失败不发布注册、强owner释放后weak失效。
- 测试贡献ID曾重复，初次新增单测失败；修正测试自身后重跑通过，随后最终整组43项全过。

## 6. 13组现役serial前后对照：共同规则

路径 `test/resourceA5/{before,after}/<name>`；由 `regression.py` 从当前可执行输入复制，不覆盖正式case。命令为对应binary的 `explain <case>` 和 `run <case>`，cwd为case；配置、explain.log、run.log、result完整保留。来源和binary逐例见 `evidence/cases-before.json/cases-after.json`。

SI量纲/scale1；本节所有case无活动IBM/ILW、无物体，IBM wallClosure/热边界不适用。输出在各自result。比较全部VTS/VTM、dt/Plan trace；`verification.json`逐例保存rho/p extrema、可用质量/动量/能量、各场L2，采用实际坐标梯形体积权重，不用nodal sum冒充积分。不输出的能量不伪造预算。

### 6.1 PISO / SIMPLE / PIMPLE

来自 `test/coreSinglePath/before/<name>`，不改变输入。U momentum+div(U)=0、p multiplier；rhoConst1、mu.01 Pa·s/Pr.72，无密度/能量/RAS输运。域[-.5,.5]²×[0,.01]，实际80×80×1cells/81×81×2nodes、NG3；未注册的16² blockMeshDict不控制网格。

IC U0、p101300、输入rho1。top U(1,0,0)，left/right/bottom U0；四周p/rho zeroGradient，z empty。ForwardEuler intent、primitiveUpwind1、Central2、Rhie–Chow；HYPRE FlexGMRES/BoomerAMG max1000/Krylov50、relative1e-10/absolute1e-12、referenceCell0/p101325。

共同start0/end1/endStep2/CFL.1/maxDt.001，初始及每2步输出。实际均2步，dt两次.000390625、终点.00078125，2帧。

| case | 方法配置 |
|---|---|
| PISO | outer1/pressure2/nonOrth0，U/p relax1 |
| SIMPLE | outer最多3/pressure1/nonOrth0，relax .8/.8，迭代relative.001/absolute1e-8 |
| PIMPLE | outer最多3/pressure2/nonOrth0，relax及迭代容差同SIMPLE |

### 6.2 conservativePressure

来源 `test/coreSinglePath/before/conservativePressure`，同上域/网格/NG3。NS Q+原压力修正；PG gamma1.4/R287.05、mu0/Pr.72。IC rho1、U(.05,0,0)、p101325；四边U0、rho/p zeroGradient，z empty。Euler/PISO1/1/0、WENO5/StegerWarming、无活动黏性。start0/end1e-5/endStep1/CFL.1/maxDt1e-5，实际1步t=dt=1.659315106641e-6，2帧。

### 6.3 sod / stageSST / eulerianSST

均来自 `test/coreSinglePath/before/<name>`，配置不变；实际只跑配置中的1步，不称完整物理验收。

| case | WHAT/STATE、网格、物性与IC/BC | HOW/WHICH及实际终止 |
|---|---|---|
| sod | NS Q，PG1.4/R287.05、mu0/Pr.72；[0,100]×[0,20]×[0,.1]，100×20×1cells/101×21×2nodes、NG4。左右rho1/.125、p1e5/1e4、U0；nativeRight/nativeSides zeroGradient，frontBack empty，nativeLeft未显式声明，无RAS | SSPRK3/TENO5/StegerWarming；start0/end.001/endStep1/CFL.5，无显式maxDt；每步输出，实际1步6.681531047811e-4，2帧 |
| stageSST | 同域/网格；NS Q+stageConservative k/omega/mu_t；PG同上、mu1e-5/Pr.72。rho1.2、U(.3,0,0)、p1e5、k.01、omega1。flow BC同sod；k/omega还有nativeLeft zeroGradient，right/sides同样zeroGradient、z empty。epsilon/phi配置保留不代表活动PDE | RK4/WENO7/Rusanov/Central2、共享质量通量Upwind1；start0/end.001/endStep1/CFL.5、无显式maxDt；实际1步7.316037672350e-4，2帧 |
| eulerianSST | [0,1]×[0,.2]×[0,.01]，20×4×1cells/21×5×2nodes。逐相mass/momentum/enthalpy、shared p、相质量加权SST。alpha .9/.1、U0/T300/p101325/k.01/omega1。水perfectLiquid rho998.2/mu.001002/cp=cv4182/kappa.598；空气PG输入rho1.1766/R287.05/mu1.846e-5/cp1007/cv718/kappa.0262。Left固定这些值，Right zeroGradient，top/bottom U0其余zeroGradient，z empty；k/omega Left固定.01/1其余zeroGradient，walls bottom/top。MRF中心(.5,.1,0)、z轴omega1；SchillerNaumann d.001、virtualMass.5、lift0、RanzMarshall，无Lee/RPI | Euler intent/PIMPLE3/2/nonOrth1即2passes、U/p relax .7/.3；phaseUpwind、turbCentral2/Upwind1。HYPRE FlexGMRES/BoomerAMG max100/rel1e-8/abs1e-12/Krylov30/refresh10；start0/end.0002/endStep1/CFL.5/maxDt.0001，实际1步3.624299498372e-5，2帧 |

### 6.4 五组独立named scalar

来源 `test/nativeEquationClosure/cases/<name>` 的独立副本。旧fixture未声明量纲，当前reader拒绝；只在本轮副本补齐field无量纲、D的L²/T及reaction系数T⁻¹，先用before binary实际冻结，未修改正式fixture或数学。

geometry-only，Q=0/EOS=none/无物理模型。域[0,1]×[0,.25]²，33×5×5nodes/32×4×4cells、NG1，serial1patch。首标量默认0、x≤.5初始化1；第二默认.5、x>.5初始化2；第三默认.75、x≤.5初始化.2。首标量Left固定1/Right固定0，其他四面zeroGradient；其他标量六面zeroGradient。

WHAT：`ddt(c)=diffusion(D,c)`，D分别.03/.08/.11（有第三个时）。start0/end.005/endStep未配置/maxDt.0005/CFL.5，初始及每步输出；实际全部10步/.005、11帧。没有速度对流、热边界、pressure coupling或隐式数学。

| case | 变量/方程 | 时间及额外数学 |
|---|---|---|
| pairEuler | dye/salt，spreadDye/spreadSalt | ForwardEuler |
| pairRK4 | 同上 | classicalRK4 |
| renamedRK4 | colourant/dopant，diffuseColour/diffuseDopant | RK4，仅改名 |
| tripleRK4 | colourant/dopant/signal，diffuseColour/diffuseDopant/diffuseSignal | RK4，增加第三独立扩散方程 |
| linearTermsRK4 | dye/salt | RK4；dye RHS额外+.2*dye，salt额外+.15 |

### 6.5 levelSet（现役Sod+伪时间）

来源 `test/coreSinglePath/interface/after/final/preset`，不改变输入。NS Q+physical phi advection、冻结phi0、3次pseudo reinitialization、normal/curvature；PG1.4/R287.05、mu.001/Pr.72，water/air属性rho1/.125、mu.001/1.48e-5，不把它们当作新EOS。无RAS/IBM/表面张力/质量修正。

域[0,100]×[0,20]×[0,.1]，101×21×2nodes/100×20×1cells，serial单块。U0、左右rho1/.125、p1e5/1e4；phi为x50法向(-1,0,0)距离。flow Right/sideWalls zeroGradient、z empty、Left未显式声明；phi Left/Right/sideWalls zeroGradient、z empty。

SSPRK3/WENO5/StegerWarming/Central2；phi advection/reinit5、epsilon1e-6/power2、pseudoDt.01/signFactor1、narrowBand.02。start0/end.003/endStep3/CFL.3/maxDt.001，runTime每.001输出；实际3步/.001、2帧，dt为.0004008919/.0003624616/.0002366466。伪时间未改变物理clock。

**13组同配置前后所有比较payload与clock/trace完全一致**；这不意味着历史参考全部匹配或所有跨配置阈值通过。

## 7. 首个实际flow+双标量组合：逐例设置及对照

路径 `test/resourceA5/composition/<name>`；由 `composition.py/references.py` 生成。共有设置：域[0,1]×[0,.25]²、33×5×5nodes/32×4×4cells，serial单patch、scale1，未加IBM/RAS。标量的初值/边界同6.4节两标量；其真实数学为：

```
ddt(dye)  = diffusion(.03,dye)  + .7*salt
ddt(salt) = diffusion(.08,salt) + .7*dye
```

scalar无量纲，D量纲L²/T、.7量纲T⁻¹。不是通过虚假diffusion开放纯反应能力；扩散真实参与空间计算。RK4的两个RHS读取同一Stage，两者在全部RHS完成后推进/发布。

| case / 来源与修改 | 实际配置与目的 |
|---|---|
| flowPair | pairRK4的输入结构，增加NS WHAT/原rho/rhoU/rhoE STATE及3个ConservativeResidual调用，增加上述交叉反应。NG3；PG gamma1.4/R287.05、hConst、const mu0/Pr.72；rho1、U(.1,0,0)、p1，全部六面rho/U/p zeroGradient。WENO5/StegerWarming、characteristic/shared flux、无viscous operator、ILW0；scalar Central2。HOW先3个flow calls，再dye/salt，最后同一Commit |
| reversedPair | 与flowPair完全同数学/物性/网格/IC/BC，只把两个scalar calls放在3个flow calls之前；检验共同Stage而非顺序偷写 |
| scalarPair | 同两条scalar数学、IC/BC/域，移除flow声明/物性/Q；NG1，geometry-only，用相同binder验证无流体配置 |
| referenceScalar | scalarPair的独立输入副本，**before binary** 的旧scalar资源装配运行，作为已有数值能力冻结对照；不调用新binder |
| referenceFlow | flowPair的独立副本，去掉两个scalar的字段、公式、HOW/WHICH，仅保留原NS；NG3，由**before binary**运行，作为流场保持对照 |

五例共同start0/end.005、无endStep、CFL.5/maxDt.0005、初始和每步输出，实际均10步/.005、11帧。flow热边界不另加壁面，全部zeroGradient；没有物体、no-slip/slip不适用。结果保留各case的explain.log/run.log/result，binary/命令见脚本和证据。

比较结果：flowPair与reversedPair全部输出逐元素一致；两标量均与scalarPair/referenceScalar逐元素一致；flowPair的所有纯流动输出与referenceFlow逐元素一致。没有改变flow来迁就scalar。内部集中测试另以手工数组/端口对照binder，覆盖独立装配一致性，不能把以上数值对照冒称“通用资源已自动生成”。

## 8. 输入/资源负例

`test/resourceA5/negative/<name>` 从flowPair独立复制，只改表中因素。域/网格、数学、初始/BC、HOW/WHICH及终止设置除表中改动均同第7节。所有负例exit1，**没有VTS输出/实际物理步**；它们不算CFD时域样例。

| case | 唯一改变 / 实际阻断 |
|---|---|
| noWhich | bindings置空；explicit HOW要求显式非空WHICH，不从preset继承 |
| missingField | 删除dye field declaration；registered internalField的unknown field dye首先拒绝，不伪造更晚的资源错误 |
| missingInitial | 删除dye初值；STATE initialization backing缺失 |
| missingBoundary | 删除dye的Left；实际资源绑定报axis0/side0缺边界，在初始输出之前拒绝 |
| unsupportedImplicit | 把salt method改为LinearEquation；实际首先被Stage STATE use不在声明StageLoop的契约拒绝。没有进入Jacobian/implicit solve，不能据此声称验证了其数值算法 |
| unknownMethod | 改为UnimplementedA5Method；没有registered implementation，拒绝且无fallback |

typed owner/type/access/重复/lifetime负例在集中测试中，不要求通过native格式暴露C++类型。实际错误、命令和binary见 `evidence/negative.json`。

## 9. 已知压力异常：保留原失败，区分线性精度

独立路径 `test/resourceA5/pressure/{before,after}/{original-serial,original-MPI2,tight-serial,tight-MPI2}`，均来自6.1的PISO。每组保持80²网格、U/p/rho、BC、PISO1/2/0、Euler/Upwind1/Central2/RhieChow、start0/end1/endStep2等设置。MPI2唯一拓扑变化为2ranks/split[2,1,1]；tight唯一数值控制变化为algorithm.yaml的absolute linear tolerance **1e-12→1e-16**，relative1e-10不变。tight是明确不同输入的诊断，未修改正式case或原比较阈值。

before和after八次实际host运行均exit0，实际2步/.00078125。每个case保留run.log/result/assembly；首个矩阵按物理坐标重建。两种容差下首矩阵/RHS均逐项相同，首RHS L2约6.25e-6；同配置前后payload完全一致。

| 输入 | max U差 serial/MPI | 压力去常量后shape差 | 原1e-6压力阈值 |
|---|---:|---:|---|
| original | 5.548720974878074e-10 | 1.0851072147488594e-6 | **Fail，保留** |
| tight诊断 | 3.674784001400644e-14 | 1.3096723705530167e-10 | Pass，仅针对诊断输入 |

此证据支持既有差异由线性求解绝对收敛精度引起，不能将original阈值失败写成已修复。HYPRE backend允许绝对或相对残差停止；tiny RHS的第二次correction可因absolute tolerance直接0次迭代。该停止条件和压力算法本轮未改。没有以tight结果替换original断言，也不推广为所有压力case的误差归因。

## 10. Level-set MPI与历史参考：分开验收

### 10.1 当前可复现D1样例

路径 `test/resourceA5/levelSetMpi/{before,after}/{serial,mpi2x}`，来源 `test/levelSetMpiD1/closeout/planeAdvection/<name>`，配置未改。NS Q+level-set physical advection及normal/curvature/material；无reinitialization、表面张力、质量修正、IBM或RAS。

域[0,100]×[0,20]×[0,.1]、100×20×1cells/101×21×2nodes，capacity NG4；PG1.4/R287.05、hConst/const mu.001 Pa·s/Pr.72。water/air的material rho1/.125、mu.001/1.48e-5。IC flow rho1、p1e5、U(.2,.1,0)；phi为过(50,10,0)、normal(1,.37,0)、scale1平面距离。rho/U/p/phi的Left/Right/sideWalls zeroGradient，frontBack empty；保留k/epsilon0及边界但无RAS。

SSPRK3/WENO5/StegerWarming/Central2；HJ advection5、WENO epsilon1e-6/power2、interfaceThickness1.5/narrowBand.02、reinitSteps0，pseudoDt.01等未被执行。serial1patch；MPI2split[2,1,1]。start0/end.001/endStep3/CFL.3/maxDt1/3000、初始及每步输出；实际3步/.001、4帧。

四次运行均exit0。同配置前后serial的4个payload、MPI的12个payload一致。MPI合并日志存在跨rank交错，比较trace多重集/clock而非伪称可从未标rank日志证明每rank完整顺序；serial trace完全相同，冻结Plan未变。

最终serial/MPI物理坐标对应：phi、curvature、U、p、rho、momentum、energy等场最大差均0；shared replica最大差0；PhaseIndicatorIntegral/PhaseMassIntegral差5.684341886080802e-14，相对诊断差2.842176702864442e-16。该结论仅覆盖倾斜平面均匀平流的3步D1 publication样例，不代表任意界面/MPI长期精度。

### 10.2 历史isolated-kernel参考仍不匹配

`test/interface/nativeLevelSet/frozen-kernel.json`明确记录其来源是临时绕过旧HOW guard的isolated kernel，**不是可运行的旧production baseline**，binary SHA为6dcdac…e61d。当前6.5样例在before/after的t0 hash和第二步dt都与其不同，说明不能简单说成末态输出路径变化；仍保留historicalMatch=false。

D1之后已发生真实BC donor、normal/curvature publication和诊断SUM修复，历史参考也没有封存完整原场可逐cell重建，因此仅靠哈希无法把历史差异归因到某一具体修改。结论为：**参考可复现性不足、首差异归因未完成**，不是本轮新增回归。未替换其哈希、容差或原失败测试。本轮用可执行起点binary及当前D1输入单独证明迁移保持。

## 11. Host MPI、数值不变量与证据

实际host命令除上面的MPI pressure/level-set外，还运行：

- `mpirun -np 4 build/test_distributedConstraintLayout`：exit0。
- `mpirun -np 4 build/test_sparseCanonicalCopy`：exit0。
- `mpirun -np 2 build/test_distributedKKT`：exit0；这是最小代数/HYPRE保护，不是production KKT CFD。

命令在脚本中使用绝对binary/case路径，cwd全在test/resourceA5；实际退出码/binary SHA存 `evidence/mpi-before.json/mpi-after.json`、`pressure-before.json/pressure-after.json`。无sandbox launch失败。外置盘AppleDouble清理warning可能出现，exit仍0；不当成数值断言失败，也不掩盖原始日志。

Stage coefficient/time、累加顺序、workspace clear、source clear、residual sign、BC→halo→operator、canonical单F* COPY、GlobalDof SUM的数值实现均未改。未新增通信或重排。13组同配置字段/clock/trace、首压力矩阵和MPI保护支持所选范围保持；没有对每个kernel增加全checkpoint instrumentation，不能声称所有模块组合都已证明。

## 12. 可复现命令与未覆盖范围

从仓库根目录：

```
DEVELOPER_DIR=/Library/Developer/CommandLineTools TMPDIR="$PWD/test/resourceA5/evidence" cmake --build build --parallel 4
python3 -B tools/check_architecture.py
python3 -B test/resourceA5/run_tests.py
python3 -B test/resourceA5/regression.py before
python3 -B test/resourceA5/regression.py after
python3 -B test/resourceA5/composition.py
python3 -B test/resourceA5/references.py
python3 -B test/resourceA5/negative.py
python3 -B test/resourceA5/mpi.py before
python3 -B test/resourceA5/mpi.py after
python3 -B test/resourceA5/pressure.py before
python3 -B test/resourceA5/pressure.py after
python3 -B test/resourceA5/verify.py
```

MPI命令需要允许PRTE socket的host。首次regression before先生成/校验typed副本；脚本不会覆盖正式case。所有数值输出各自result，无系统临时目录测试。最终关键日志：build-final.log、architecture-final.log、ctest-final.log、regression-validated.log、composition-validated.log、negative-validated.log、mpi-final.log、pressure-final.log、verification-final.log；结构化比较与加权统计verification.json。

本轮没有修改IBM算法，不重跑其全部CFD，nativeIBM/viscous wall/surfaceSchur集中保护通过不能代替forcing/DLM/KKT完整验收。未宣称任意EOS、用户PDE或资源类型支持。

后置工作：Eulerian/flow其余资源装配；一般phase/face/MPI资源与restart；纯反应、多个/变系数/交叉扩散、质量矩阵、非线性隐式/DAE；Subcycle/通用Newton；任意模块组合；历史level-set完整首差异恢复；压力original配置的精度策略。它们分别需要真实实现和独立验收，不重建WHAT。

**完成范围：公共/具体lowering边界已收口；flow+双标量与独立多标量通过同一声明驱动binder实际运行；被替代scalar执行分支已删除。完整通用资源装配仍未完成。**
