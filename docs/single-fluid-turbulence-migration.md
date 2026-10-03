# Single-fluid turbulence native execution — 2026-10-03

冻结阶段开始于 2026-10-02，完成于 2026-10-03（Asia/Shanghai）。

本阶段完成 single-fluid kEpsilon / kOmegaSST 的 native WHAT、STATE、HOW、WHICH 贡献，以及通用 mixed temporal lowering。现有 RAS 数值内核不变；两个方程通过显式 provider fusion 合成一次 physical-step 更新。pressure/MPI/multi-patch/IBM/Eulerian turbulence 未放行。

## 1. 原调用链与真实启动状态

原 composition：配置 -> Turbulence::contribute -> NamedDistributed k/omega 或 k/epsilon -> legacy EquationDescriptor/Definition -> order 50/51 addLegacyExecution。flow 的三个 native conservative occurrences 编译为 RK StageLoop，但两个 RAS legacy entries 没有 compiled HOW ownership。

**原可执行文件实际在启动检查处拒绝 RAS**：`Equation 'k' has no compiled HOW execution step. Equation 'omega' has no compiled HOW execution step.` 该失败保存在 evidence/original-entry-failure.log。不能称为已有可运行的历史湍流基线。

被启动检查挡住的旧 production 内核调用链是：application 创建 Manager -> CompositeTransportProvider -> FlowStepBegin -> correctTransportModel -> CompositeTransportProvider::correct -> Manager::correct -> KEpsilonModel / KOmegaSSTModel::correct -> Time::Explicit::begin。每个 FlowStepBegin 只调用一次 correction；RK stages 不再次推进 RAS。

审计读取了 turbulenceSystemContribution、RAS 两个 correct 内核、kEpsilonEquation/kOmegaSSTEquation、equationSystem、equationModelOps、turbulence Manager、CompositeTransportProvider、SingleFluidStepper 和 compiler。后面几个 EquationSystem/prepare 文件属于 Eulerian kernel，不能拿它们的对流/隐式输运能力来描述 single-fluid RAS。

## 2. 数学审计：原 WHAT 与实际内核

| 项 | 原注册描述 | single-fluid 实际算术 / 本轮 native WHAT |
|---|---|---|
| transient | ddt(k/second) | primitive scalar 的一次 dt 更新；没有 rho*k conserved unknown |
| convection | div(flux.k/second) | **未执行**；native WHAT 删除 phantom div，属于 pre-existing mathematical-authority mismatch |
| diffusion | 通用左侧 diffusion | kEpsilon：右侧 div((mu_t/sigma) grad(phi))/rho；SST：右侧 div((max(laminarMu,0)+sigma*max(cached mu_t,0)) grad(phi))/rho |
| production | opaque source | kEpsilon min(mu_t*S²,50*rho*epsilon)；SST min(mu_t*S²,20*0.09*rho*k*omega)，保留原限制 |
| dissipation | opaque source | k：-epsilon 或 -0.09*k*omega；epsilon：C1*epsilon*P/max(rho*k,kFloor) - C2*epsilon²/max(k,kFloor)；omega：gamma*P/max(mu_t,1e-30) - beta*omega²；全部显式 |
| cross diffusion | source 未区分 | SST omega：2*(1-F1)*rho*sigmaW2*grad(k)·grad(omega)/max(omega,omegaFloor)，与 diffusion 一起除以 rho |
| wall / blending | 未具体声明 | 既有 wallDistance、F1/F2、vorticity limiter；invalid/small distance -> 1e20 的原行为保留 |
| floors / BC | 未具体声明 | 原 positivity floors、rho/denominator guards、标量 BC、interior floor、post-update viscosity refresh 不变 |
| update realization | 未具体声明 | k 与 second 在原 i/j/k sweep 中原位更新；不改成 simultaneous/Jacobi 或 flow RK stages |

native AST 在模型贡献中声明数学；TurbulenceTransport 的 capability matcher 要求 AST 与其冻结内核完全一致，修改/添加未实现项会明确拒绝。matcher 不创建第二个 WHAT registry。AST 中 limitedProduction、gamma/beta/sigma、maxRhoK/maxK/maxMuT、crossDiffusion 等符号表示原闭合计算及 guards，不宣称新数值项。

SST 内核内部使用原固定 SSTConstants；部分同名配置 coefficient 在 single-fluid 路径原本就未被使用。kEpsilon diffusion 原本不含 laminarMu；两个内核都缺少物理对流、使用原位邻点更新，并保留 floors。均作为后续数值研究问题记录，本轮不补项、不改系数/迭代/壁面公式、不引入 HYPRE RAS。

Eulerian prepare 另有 alpha/rho/mass 加权系数、包含 laminar viscosity、显式生产和隐式对角汇；其 AssemblyPlan/shared-pressure 执行未迁移。不能据此把 single-fluid 的缺项解释为已实现。

## 3. 修改后的权威调用链

模型贡献独立的 WHAT/STATE/HOW/WHICH -> compiler bind / validate -> provider-declared adjacent pair fusion -> generic mixed lowering -> StateRealization aliases -> owned CompiledSolvePlan -> PlanExecutor -> frozen flow.turbulence leaf -> 原 correctTransportModel / CompositeTransportProvider / Manager / RAS 算术。

ResolvedSimulationSystem 仍是冻结的 12 字段结构；未增加第五模块、另一份 SimulationPlan、全局湍流 solver family 或模型时间循环。

## 4. STATE 和 backing

| STATE | 来源 / role | backing owner / realization |
|---|---|---|
| rho,rhoU,rhoE | 用户 primary flow selection | 原 Field packed Q；StateBundle 非 owning 引用 |
| U,p,T 等 | 实际 flow dependencies / derived | 原 EOS-backed flow views |
| k + omega/epsilon | 显式模型贡献 / Transported | Manager::ScalarFields；ProviderDistributed；Physical(target) 的 owner=NumericalProvider |
| mu_t | 模型 Derived closure output | 原 EddyMu slot；不加入 primary flow selection |
| time,dt,step | 原 StateBundle clock | 不新增 clock；RAS 使用已计算的 dt |

Manager::registerDistributed 接管原 application 的手工 name->slot 注册，仍用原 workspaceView/COPY contracts，直接引用同一 vector。StateRegistry、StateRealization、DistributedFieldRegistry 都指向原数组；没有架构副本。ScalarFields 原先分配四个槽位，包括未激活模型的 unused slot，这个历史分配没有改动。

模型数组初始化/BC/refresh lifecycle 不变；在初始化注册后绑定 native Physical views，物理步内 RAS 更新后保持可读，flow stages 继续读取现有 viscosity closure。当前没有 turbulence restart reader，metadata 不再虚报 restartEligible。单元测试通过 STATE view 写入后直接读到 Manager 数组同一值；实际 kEpsilon correction 更新两个 transported arrays 和 mu_t，并证明全部 flow Q components 未被直接写入。

真实新算例只在 state/state.yaml 选择 `[rho,rhoU,rhoE]`；k/omega 来自模型贡献。其 input/output 描述沿用现有 fields 输入格式，本轮未扩展任意自定义四模块 IO。

## 5. 精确 HOW 顺序

source HOW（SST；kEpsilon 将 omega 换为 epsilon）：

```text
Sequence
    1         k          -> Physical(k)
    2         omega      -> Physical(omega)
    10        continuity -> Physical(rho)
    20        momentum   -> Physical(rhoU)
    60        energy     -> Physical(rhoE)
    1000000   Commit
```

compiled / 实际执行：

```text
flow.step.prepare
flow.dt.compute
[ k->k, omega->omega ] -> turbulence.advance   owner=flow.turbulence
flow.step.begin                              snapshot
StageLoop classicalRK4
    explicit.stage.execute × 4               owner=flow.conservative
flow.step.commit
time.commit
```

RAS 一次/物理步显式原位更新使用自身冻结 method contract，flow RK4 不应用到 k/omega。dt 的计算、CFL、flow old-time snapshot、stage times、物理提交和 clock 提交保持不变。

## 6. WHICH、fusion 和 compiler

两个 RAS occurrences 各选择 `TurbulenceTransport`，并显式声明 pair inputs `[k,omega]` 或 `[k,epsilon]`。编译产品的 fusionKey 为 TurbulenceTransport/kOmegaSST 或 TurbulenceTransport/kEpsilon，fusionMembers 保留两个 equation->target references；operation=turbulence.advance，provider=flow.turbulence。其能力只有 single-fluid physical-step explicit RAS，不宣称 MPI、压力、Eulerian 或隐式输运。

Fusion 验证完整 equation set、正确目标/模型 STATE/数学、相邻 sibling、同一 Sequence、相同方法/backend/合约、无重复及无跨 scope。源 WHAT 和 HOW 保留两份；compiled leaf 保留两份 mathematical provenance；不存在 RASBlockEquation，也不因共享 provider ID 自动 fusion。

原 compiler 在 transient.size()!=steps.size() 时拒绝 mixed fused-stage/independent providers。现在通用编译器支持最小 topology：独立 prefix -> 一个连续 fused temporal group -> terminal Commit。ITemporalMethod 接收编译后的 physicalStepPrefix；explicit recipe 将它放在自己的 prepare/dt 后、begin 前。generic compiler 不读取 k/omega/epsilon、模型名或 turbulence provider ID。

temporal -> independent -> temporal、independent tail、nested mixed scope、缺少 terminal Commit，以及尚未明确指定 temporal placement 的额外 root lifecycle decorations 均明确 Unsupported，不静默重排或丢弃。synthetic Independent/A/B compiler test 完全独立于湍流名，验证一次 prefix、四个 flow stages、原 commits、交错布局失败和 root decoration 拒绝。

provider 的 usesSpatialRecipes 合约让全局 diffusion defaults 只匹配实际 spatial consumers；RAS 内核固定的局部 Central2 不误触发 main-flow recipe。Conservative provider/callback 的 tuple 检查现在只检查其 temporalResidual occurrences，仍严格要求 rho/rhoU/rhoE layout。

## 7. 旧路径去向和边界/halo

- **删除 single-fluid authority**：RAS legacy EquationDescriptor/Definition、order 50/51 addLegacyExecution；改为 native AST、order 1/2 EquationCalls 和 explicit numerical bindings。
- **删除 scheduling call**：FlowStepBegin 中的 correctTransportModel；begin 保留现有 flow begin/snapshot、workspace、temporal views 职责。
- **删除 application execution predicate**：transportedTurbulence 配置推断及旧 provider 比较；compiled requirement 决定需绑定的更新，配置仅创建模型/参数和输出。
- **移动 registration**：application 手工标量注册变为 Manager::registerDistributed；无新的数组或数据 COPY。
- **保留 provider arithmetic**：correctTransportModel、CompositeTransportProvider::correct、Manager::correct 和两个 RAS correct 内核。它们并未被“删除”；生产 RAS 只能通过 compiled turbulence leaf 调用。
- **保留 legacy Eulerian island**：其 legacy equations/execution/policy extensions、EulerianStepper、AssemblyPlan/TermKind、sharedPressurePlanFragment、EeTurbulenceSolve。

原 correction 内部精确顺序保留：applyBoundary -> execution finalize(WriteOwned) -> prepare(ReadHalo) -> transportModel->correct(dt) -> applyBoundary -> finalize(WriteOwned) -> prepare(ReadHalo)。RAS correct 内部原有 post-update BC 和 EddyMu refresh 也保留。它们分别属于 boundary closure、publication/halo contract 和局部数值更新；没有把 raw MPI 或新的同步 node 放入 source HOW。

Smagorinsky 的旧缓存刷新随 scheduling call 一起迁出：mu_t 代数 closure + TurbulenceClosure -> turbulence.closure.refresh，位置仍在 dt 与 snapshot 之间；不新增 transported state/equation。DNS 没有 transported equation update。保留 ITransportModel 的 viscosity/boundary/data/correct 接口，没有 manager/coordinator/scheduler 包装。

## 8. 冻结基线的来源和数值结果

没有现成 enabled、可运行的单流体 RAS case。**仅新增一个持久算例**：test/turbulence/singleFluidSST，从已有 WENO7/Rusanov Sod 模板派生，serial，RK4，Central2，endTime=0.003，k=0.01、omega=1。input hashes 包括 mesh、配置、初始化和边界。

原入口启动失败。为比较被挡住的现有算术，先冻结 source/binary，并构建 reference-only diagnostic binary：只给 `k`/`omega` 的启动 ownership 检查增加豁免；没有修改 HOW、FlowStepBegin 调度、任何算术、CFL、BC 或时间算法。reference-only.patch、构建日志、原二进制、诊断二进制、完整原源码和失败日志均保存。此豁免已撤销，最终生产源码仍严格检查 ownership。初始 legacy-source 输入对照另保存在 legacy-reference；最终 native 显式输入用同一个旧内核诊断二进制冻结，两者输出相同。

| 对比项 | 旧内核 reference / native after | 分类 |
|---|---|---|
| steps / final time | 6 / 0.003 | Passed identical |
| time/dt、flow plan stage contexts | 全序列相同；after 只增加一次/步 turbulence.advance | Passed identical |
| k、omega、mu_t、rho/U/p/T/energy | 四个输出每个数组的 min/max/L2/sum 相同 | Passed identical |
| active XI/ETA 边界统计 | 各输出边缘数组的 min/max/L2/sum 相同 | Passed identical |
| VTS | 4/4 SHA256 逐字节相同 | Passed identical |
| closure / high-order checkpoints | 完全相同；仅规范化进程地址文本 | Passed identical |
| execution count | 6 次 pair advance、24 次 flow stage | Passed |
| protected numerical sources | 76 个实际源文件、9 个关键 function bodies 完全相同；hash 清单另含 66 个 macOS 附属文件 | Passed identical |
| kEpsilon CFD | 没有另造第二个 CFD case；native structure/storage/correction unit 验证通过 | Not tested (CFD) |
| RAS pressure / MPI / multi-patch / IBM | explicit capability rejection | Unsupported |
| Eulerian / 近壁 SST CFD | 本轮未运行 | Not tested / deferred |

SST 最终物理输出统计（reference 和 after 相同）：

| 数量 | min | max | L2 |
|---|---:|---:|---:|
| k | 0.009997300577664637 | 0.010039358010026641 | 0.6513132720960351 |
| omega | 0.9997516509311345 | 1.0016004117970956 | 65.12248685924914 |
| mu_t | 0.0012499324932319267 | 0.010000329319918413 | 0.4598380715148193 |
| rho | 0.12499047708619548 | 1.00000600071477 | 45.95238832725648 |
| p | 9999.062742917338 | 100000.95781952726 | 4574669.868980194 |
| T | 278.6922017113756 | 348.3726101581729 | 20559.947407368614 |

最终 VTS SHA256：`40d59779c9645dfc3b5314c3284d9f4ec0736728ca256bf5b61ed12c85d1d0c5`。完整 hashes、U/vector statistics、dt、边界值和 diagnostics 在 baseline.json 与 evidence before/after/evidence.json。

生产内核不打印 turbulence floor-hit 计数；报告不虚构此项。原 guards/floors 源码逐字相同，输出正值；这不证明其物理合理性。小型 Sod 变体不覆盖近壁 SST。没有调整任何 frozen 参数来取得相等。

同一输入模板的 temporary Smagorinsky 选择变体另外验证 closure refresh：旧入口与新入口均成功，4 份 VTS hash 相同；不作为额外持久 RAS regression case。初次变体因仍要求输出不存在的 omega 被两侧同样拒绝，随后只关闭该输出 request，输入在两侧一致；没有修改数值参数。

## 9. 构建和 regression tests

完整 `cmake --build build -j 8` 成功。新增 nativeTurbulenceExecution、mixedTemporalComposition 和 singleFluidTurbulenceNumericalRegression；architecture checker 根据去注释 function/lambda bodies 检查 native 四模块、FlowStepBegin 去调度、frozen leaf owner 和 generic temporal compiler。

完整 CTest 首轮 **30/31** 通过；唯一失败 pressureTermProviders 仍要求旧 equation.turbulence-transport 名称，实际 runtime 已报告 Unsupported。将其 assertion 更新为 flow.turbulence/turbulence.advance 后该项单独复验通过。最终三项新测试及 methodObjectCompilation 再次通过；故 **31 个测试各自的最终结果全部 Passed**，不将首轮有失败的日志伪称一次全绿。

既有 Sod、constant-density PISO/Poiseuille、fixed-time SIMPLE/PIMPLE、conservative PISO、distributed pressure regression 和 MPI/halo contracts 均通过。这些 MPI 测试保护既有 pressure/runtime 行为，**不证明 RAS MPI 能力**。

check_architecture.py、check_program_decomposition.py（scope/fieldcount）、git diff --check 均通过。新测试、算例、runner 和报告增加精确 gitignore 例外；result、诊断二进制和 evidence 保持忽略，不改变其他未提交工作。真实 pressure+RAS explain 显示 Unsupported，run 在既有测试中被拒绝；parallel+RAS 与有效 Ghost IBM geometry+RAS 在 SingleFluidTurbulenceScope 处拒绝。pressure service guard 未删除。

复验入口：

```sh
cmake --build build -j 8
ctest --test-dir build -R 'nativeTurbulenceExecution|mixedTemporalComposition|singleFluidTurbulenceNumericalRegression' --output-on-failure
python3 tools/check_single_fluid_turbulence.py --solver build/sonicSolver --case test/turbulence/singleFluidSST
python3 tools/check_architecture.py --quiet
```

## 10. AGENTS 影响报告

1. **层**：模型数学贡献、STATE backing contract、numerical method/provider、通用 temporal lowering、runtime binding、application setup、explain/文档/测试。
2. **类型**：已有模型 equations/closure 迁入 native contribution；改调度/数值实现选择；没有给物理内核增加 term。
3. **ownership**：旧 legacy scheduling / FlowStepBegin -> HOW + explicit fusion + owned compiled leaf；数组依然属于 Manager，clock 依然属于 StateBundle。
4. **execution order**：RAS 的数值位置不变；旧隐式调用变为 dt 后、snapshot 前的显式 compiled operation。
5. **public API**：新增 ProviderDistributed、CompiledEquationCall fusionKey/fusionMembers、ITemporalMethod physicalStepPrefix、IProvider usesSpatialRecipes、Manager::registerDistributed；不新增 source HOW node 或 STATE authority。
6. **dependency**：模型只提供 neutral AST/state/scope/binding 值；generic compiler 不 include concrete turbulence model。provider 的数学 capability matcher 与真实 kernel 保持一致；未增加 reverse dependency。
7. **parallel**：COPY/ReadHalo/WriteOwned 和 SUM/canonical face invariants 未改；RAS 能力限制明确保留。
8. **numerical semantics**：RAS、flow、RK、flux、diffusion、BC、guards、viscosity 算术保持；WHAT 删除未执行的 phantom convection 是描述修正。
9. **tests**：31 个最终结果 Passed，SST 4/4 VTS 和 diagnostics identical；Smagorinsky closure 4/4 VTS identical；完整 build 和静态守卫通过。
10. **deferred**：RAS numerical defects/accuracy/近壁、kEpsilon CFD、pressure/MPI/IBM RAS、Eulerian/shared pressure、AssemblyPlan/TermKind、用户四模块 IO、thermophysical generalization。

## 11. 18 项自审

- [x] single-fluid transported RAS 不再使用 addLegacyExecution；Eulerian 明确保留。
- [x] k/omega 或 k/epsilon 仍是独立 WHAT equations。
- [x] 模型 STATE 在同一 StateRegistry 内显式贡献。
- [x] 用户 primary flow selection 保持 rho/rhoU/rhoE。
- [x] HOW 只使用通用 EquationCall/Sequence/Commit。
- [x] TurbulenceTransport WHICH binding 显式。
- [x] provider owner 在编译时冻结。
- [x] pair kernel 一次/物理步；SST 6 次、flow RK 24 stages。
- [x] FlowStepBegin 无 transport correction scheduling。
- [x] mixed temporal lowering 通用；synthetic test 无湍流 routing。
- [x] 未将 flow RK recipe 应用于 RAS。
- [x] 现有 numerical arithmetic 保持；protected sources/bodies 相同。
- [x] WHAT/kernel mismatch 已审计、修正描述并报告。
- [x] runtime arrays 无副本，alias/publication 实测。
- [x] pressure+RAS 未虚报支持。
- [x] MPI/multi-patch RAS 未虚报支持。
- [x] Eulerian turbulence 显式 deferred。
- [x] 真实 SST CFD 与冻结旧内核诊断对照一致，基线来源限制明确。

证据索引：[evidence manifest](evidence/turbulence-native-20261002/manifest.json)。
