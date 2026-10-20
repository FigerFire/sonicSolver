# Phase 4A：Field 与 Solver Data Responsibility 审计

审计日期：2026-09-11。本文依据当前源码重新检查，而不是复述 Phase 0 结论。未修改
`src/`、CMake 或 public API。`python3 tools/check_architecture.py` 仍报告 10 条
allowlisted dependency edge。

## 1. Current Field Responsibility Map

当前 `Field` 的最小、不可替代核心是**一个结构 patch 上的守恒存储及其索引布局**：
`data`、含 halo 的尺寸、`NVar` 和 storage layout。`StateBundle` 不是该数组的第二份
拷贝；它以非拥有的 `patches` 指针组成参与求解的 patch 集合，并持有 EOS binding、
registry、时钟和 distributed view 注册表（`src/core/state/SF_stateBundle.h:18-30`）。
因此目前的关系已经是：

```text
StateBundle：求解状态组合、时钟、EOS/registry 契约（不拥有大数组）
    └─ patches[] ──borrow──> Field：patch-local Q 与历史上聚合的其它数据
                                  ├─ geometry / boundary / distributed identity
                                  ├─ thermodynamic cache
                                  ├─ residual / face-flux stage workspace
                                  └─ IBM geometry / classification
```

这符合 Phase 2 后的“一个 physical state contract”，却仍不符合
`ARCHITECTURE.md:570-610` 对 patch state 的职责分离目标：后者要求区分 structured
patch、flow state、solver workspace，并要求 IBM-specific state 优先由 IBM subsystem
持有。问题不是 `Field` 与 `StateBundle` 同时拥有 Q；问题是单个 patch 对象同时承担
至少 11 类、生命周期不同的数据。

量化口径：把 `mx/my/mz/ng` 计为四个直接成员，重载分别计为一个 public callable
declaration。`Field` 有 **16 个直接语义数据成员**（15 个声明行）、约 **118 个 public
callable declarations**、**28 个返回可写指针/引用的直接入口**；再计入 setup、clear、
set、mark 和 stage 写入操作，共 **51 个可写 API family**。现有 14 个类别中，A/B/C/D/F/H/I/J/K/L/N 存在；
E（辅助物理数组）、G（独立 reconstruction 临时数组）和 M（诊断/输出数组）不在
`Field` 内。

`GeometryStorage`、`BoundaryMetadata` 与 `IBMGeometryStorage` 已是有用的内部分类，
但它们仍被 `Field` 唯一拥有，外部通过 `Field` 的 mutable accessor 写入；它们尚未
表达 subsystem owner，也没有访问权限边界。

## 2. Field Member Ownership Table

| Member | Category | Created by | Written by | Read by | Lifetime | Authoritative? |
| --- | --- | --- | --- | --- | --- | --- |
| `data` | C physical conservative state | `Field::setup` / `resizeConservedVariables` | initial condition、boundary/halo、DensityBasedTime、forcing/correction | equation、CFL、boundary、output | setup 后至 patch 销毁；每 RK stage 改写 | 是，patch-local `Q` |
| `mx,my,mz,ng,total_size` | A mesh topology | `Field::setup` | mesh setup | 所有索引、halo、workspace allocation | patch mesh lifetime | 是，local structured extent；不是全局 connectivity |
| `nVar_`、`storageLayout_` | N layout/configuration | constructor/setup/EOS factory | construction 或显式 resize 前 | Q indexing、flux/residual、distributed view | allocation lifetime | 是，Q layout contract |
| `equationSet_` | D thermodynamic binding | application/EOS factory | composition/setup | `thermodynamicState`、StateBundle consistency validation | workflow binding lifetime | Field 的 local binding；density timestep 的逻辑权威仍是 `StateBundle::equations` |
| `thermodynamicCache_`、`stateVersion_` | D thermodynamic cache | `setEquationSet` / setup | lazy `thermodynamicState` fill；各 Q writer 后 `invalidateThermodynamicCache` | boundary closure、viscous、EOS utilities | Q version lifetime | 否，完全由 `data + equationSet_` 派生 |
| `geometry_.x/y/z` | B geometry | mesh loader/generator | mesh generation、decomposition、geometry halo fill | all geometry-aware discretization、IBM、output | mesh lifetime；halo copy 时更新 replica | 是，local coordinate storage |
| `geometry_.jac, xi*, eta*, zeta*` | B metrics | mesh metric build / SFM load | mesh construction、metric halo fill | CFL、div/grad/laplacian、ILW、pressure path | mesh lifetime | 是，patch metric storage |
| `geometry_.canonicalFaceMetrics`、`canonicalFaceMetricMask`、`canonicalFaceOwnerMask` | J stable distributed ownership metadata | communication-plan construction | mesh/communication plan | canonical-face flux assembly | partition/topology lifetime | 是，canonical face geometry/owner metadata；不是 runtime MPI request |
| `boundary_.sets` | F boundary metadata | mesh import/normalization | mesh setup；`getSetWritable` 也允许任意 caller | physical BC、output、model boundary queries | mesh/config lifetime | 是，patch-local set-to-index map |
| `boundary_.faceMarkers` | F boundary metadata | setup allocation only | 当前 `src/` 未找到写/读调用点 | 无当前 reader | 未证实 | 否；疑似 dormant storage，P2 |
| `boundary_.cellType`、`solverMask` | K IBM classification / F boundary metadata | `Field::setup` | mesh import、CompositeIBM、ghost WeightBuilder；mask 再由 BC setup 标记 | traversal、boundary、IBM、multiphase/turbulence | mesh + IBM classification lifetime | `cellType` 是 current patch classification；`solverMask` 是 BC-derived cache |
| `boundary_.communicationHaloMask` | J distributed metadata | `Field::setup` | mesh/CompositeIBM partition mapping | MPI/IBM availability filtering | partition mapping lifetime | 是，mapped halo identifier；Field 不执行通信 |
| `boundary_.globalDofId`、`globalDofOwnerRank`、`globalDofOwnerMask` | J MPI/GlobalDof identity | communication plan | `buildHaloExchangePlan` | GlobalDof residual/IBM topology | partition topology lifetime | 是，stable identity and canonical owner metadata |
| `ibm_.ghostLayer` | L IBM runtime/ghost data | `Field::setup` | mesh import、CompositeIBM、ghost builder | ghost closure and cell traversal | IBM setup/reclassification lifetime | 是，当前 ghost-layer classification |
| `ibm_.fluidMask`、`signedDistance` | K IBM geometry/classification | `Field::setup` | mesh import、CompositeIBM、ghost builder | IBM, turbulence wall distance, output | IBM geometry lifetime | 是，current IBM classification/geometry cache |
| `ibm_.wall*`、`image*`、`normal*`、`wallVelocity*` | K IBM geometry | `Field::setup` then static/preclassified IBM setup | CompositeIBM / WeightBuilder / mesh preclassified data | ghost ILW and ghost closure | IBM geometry lifetime | 是，对 ghost closure 的 patch-local geometry cache；不是 rigid-body equation state |
| `convectiveFlux_` | I face flux | `Field::setup` / resize | convection face assembly；runtime canonical COPY overwrites replicas | canonical-face finalization、residual assembly | one RHS stage | 否，stage workspace |
| `residual_` | H residual | `Field::setup` / resize | convection/diffusion/source; local residual staging; runtime GlobalDof SUM | explicit update | one RHS stage，clear 后失效 | 否，stage workspace |

`ScalarField`、`VariableRegistry`、phase/multiphase/level-set state 并不保存在
`Field`；它们以独立存储和 `StateBundle::transported/distributed` view 注册
（`src/core/state/SF_distributedField.h:27-43,93-154`）。这说明 E 类没有被悄悄塞回
`Field`，也说明未来拆分不得复制 auxiliary array。

## 3. Physical State and Thermodynamic Cache

守恒量的 authoritative source 是 `Field::data`。`StateBundle::patches` 只借用它，
`registerConservativeState()` 也只构建可读写 view（`SF_stateBundle.h:81-88`）；
Runtime 不保留第二个 Q 数组。Density-based RK 将 `Qstage` 写回同一 Field，RK 的
`q0/k1...` 是 `SF_densityBasedTime.cpp` 的函数内 workspace，而非第二 physical
state；Phase 3B 的该项结论仍与当前实现一致。

primitive/thermodynamic state没有独立可写数组。`thermodynamicState()` 从当前 Q 调用
`EquationSet::close`，以 `stateVersion_` 命中小型 direct-mapped cache
（`src/core/field/SF_field.cpp:136-172`）。这属于 D 类 derived cache，不是
authoritative primitive state。当前风险在于 Q 有很多直接可写入口：普通
`operator()`、`conservativeData()`、`componentData()` 以及 halo/IBM/forcing 路径。
cache 失效由调用者协议负责；现有 RHS、time integrator、boundary pipeline、MPI 和
Peskin path 均有显式 invalidation 调用，但 raw pointer 写入不可能由类型系统自动发现。

`StateBundle::validateDensityEquationBinding()` 强制每 participating patch 的
`Field::equationSet()` 与 `StateBundle::equations` 指针相同
（`SF_stateBundle.h:50-66`）。因此 Phase 2 的多 EOS authoritative path 已收敛为一个
logical binding；保留 Field 的 local pointer 是为 `thermodynamicState` 提供不携带
bundle 的 patch-local closure。将它移动或删除会触及边界、viscous 和 cache 时序，属于
数值敏感改动，不是本阶段动作。

## 4. Mesh / Geometry Ownership

`Field` 当前确实**拥有** local coordinate and metric arrays，而不只是转发方便访问；
`SF_mesh.cpp:193-205,336-345` 创建 Q grid、coordinates 和 metrics，
`SF_MultiBlockMesh.cpp:171-180` 可从 SFM 写 metric。网格 halo 通过
`SF_mesh.cpp:486-495` 复制 metrics。Mesh topology/connectivity、source-block
descriptors、communication mappings 则在 `MultiBlockMesh` / communication plan，
不在 Field；Field 仅保存 structured extent、boundary set indices 和 local geometry。

因此 geometry 不能被误当作 physical solver state，也不能靠复制 `Field` 获得第二份
geometry。它的合理长期归属是 patch-local structured geometry；现有
`GeometryStorage` 足以作为未来提升的起点，不能再造 `GeometryStorage2`。但 geometry
mutable refs (`X/Y/Z/Jac/Xi*/Et*/Ze*`) 对所有 include 者开放，使 mesh-only writer
无法在接口上与 solver reader 区分。

## 5. Solver Workspace

### residual

`Residual` 已是单独的 core container，却被 `Field` 作为成员持有。它含方向 residual、
source、local staged residual、GlobalDof result 和 mask（`src/core/residual/SF_residual.h:15-82`）。
`DensityBasedRHS::assembleAllPatches` 的顺序是：每 patch 生成 face candidates → 一次
canonical-face barrier → diffusion/source/coupling → local strong residual staging → 一次
GlobalDof SUM barrier（`src/solver/algorithm/SF_densityBasedRHS.cpp:99-145`）。
`clearResidual()` 同时清残差、source、local/global staging 并清 `convectiveFlux_`
（`SF_field.cpp:216-226`）。这些值不跨 timestep 传递物理历史，属于 H 类
`SolverWorkspace candidate`。

### face flux

`FluxField` 也是独立 core container。convection 写 `ConvectiveFlux`，MPI canonical
owner COPY 后将同一 `F*` 回写每个 participant，再以 `+F*/-F*` 装配结构 residual
（`src/infrastructure/mpi/SF_haloExchange.cpp:1371-1502`；
`src/methods/numerics/structured/SF_faceAssembly.h:58-82`）。它是 I 类 one-stage
workspace，不是 committed physical state。移动它时必须保持“所有 patch candidate
完成后才进入一次 barrier”；不能把 workspace 拆分误写成每 patch 单独 finalize。

### temporary storage

Field 内没有 RK snapshot、WENO reconstruction buffer、ILW plan 或 matrix buffer。
RK `q0/k*` 在 `SF_densityBasedTime.cpp` 内局部拥有；ghost ILW plan 在
`WeightBuilder::ilwStorage_`；KKT surface/multiplier 在 IBM method state。它们是 G 类
workspace，但不应因 “Field cleanup” 被搬入 Field 或 StateBundle。

## 6. MPI / GlobalDof Metadata

`communicationHaloMask`、GlobalDof id/rank/owner mask 和 canonical-face masks 的写入
发生在 mesh/communication-plan setup；例如 GlobalDof owner 的稳定选择和写入在
`SF_communicationPlan.cpp:1685-1731`。它们描述的是 mesh-distribution identity：
谁是 owner、哪个 local cell 是 replica、哪个 face 可以计算唯一 `F*`。

MPI execution state不在 Field。communicator、neighbor payload、collective scheduling 和
actual COPY/SUM 都在 `Execution::Runtime` / `HaloExchange`；Runtime 只借用
`StateBundle*`（`SF_executionRuntime.cpp:13-16`），并且对 WriteOwned、canonical face
和 GlobalDof accumulation 执行不同 contract（`SF_executionRuntime.cpp:54-85`）。
这符合“patch 可以保存稳定身份，不能拥有 MPI request/algorithm”的原则。未来把
metadata 从 Field 移走时必须只移动 stable arrays，不能把 COPY/SUM 或 canonical
selection规则复制进新对象。

## 7. IBM Data Inside Field

`IBMGeometryStorage` 包含 15 个按 cell 分配的数组；无论 IBM 是否启用，
`Field::setup` 都分配 ghost layer、fluid mask、SDF、wall/image/normal/wall velocity
（`SF_field.cpp:65-79`）。这是“所有 Field 付费、只由 IBM feature 使用”的最清楚证据。

| Field data | Writer / update point | Consumers | Required by | Moving-body observation |
| --- | --- | --- | --- | --- |
| `cellType` | mesh preclassified import、`CompositeIB`、`WeightBuilder` | traversal、all IBM-aware models | generic flow traversal also uses it | classification may be rebuilt by IBM setup; not a physical Q variable |
| `ghostLayer` | same; CompositeIB grows layers | ghost/ILW traversal | Ghost/ILW; not Peskin/DFM/DLM by itself | no timestep writer found |
| `fluidMask`、`signedDistance` | CompositeIB/WeightBuilder SDF classification | ghost setup, k-omega wall-distance, output | Ghost plus some wall-aware closures | no per-stage field refresh found |
| wall point/image point/normal | CompositeIB, WeightBuilder, preclassified mesh | ghost WeightBuilder and ILW closure | Ghost/ILW | no per-stage field refresh found |
| `wallVelocity*` | same setup paths | ILW plan/ghost closure | moving Ghost/ILW only | Field holds setup-time cache; rigid-body `targetTime` motion and forcing surface state live in IBM subsystem |

证据表明 Ghost/ILW 需要这些 cache；Peskin、velocity forcing 和 DFM/DLM 的 target-time
surface/constraint state主要来自 `BodyModel` 和 `surfaceSystem`，并非依赖 Field 的
wall/image arrays（例如 `src/models/ibm/method/SF_methodState.cpp:18-26`）。因此 IBM
storage 不能默认视作 core flow state，也不能因为名字都叫 IBM 而把 ghost geometry、
forcing state、constraint multiplier 合成一个服务。当前也没有找到把 moving body 的
stage targetTime 回写 `Field::wallVelocity*` 的路径；这不是本报告认定的 bug，而是
Phase 4B 若触及该 cache 时必须先验证的语义边界。

## 8. Writable API Surface

| Writable access API | Data | Call sites / writer class | Legitimate writer? | Risk |
| --- | --- | --- | --- | --- |
| `conservativeData`, `componentData`, `operator()(…,v)`, `setVal` | Q | initialization、halo, RK, BC, IBM forcing | 是，但 raw pointer 不能自动失效 thermo cache | P1 cache freshness |
| `X/Y/Z/Jac/Xi*/Et*/Ze*` mutable refs | coordinates/metrics | mesh load/generation/decompose/metric halo; boundary metric fill | 仅 mesh/partition setup 应写 | P1，solver/model也能写 |
| `setCanonicalFaceMetrics/Owner` | face ownership metadata | communication-plan creation | 仅 plan construction | P1，canonical identity 若漂移会改变 interface flux |
| `getSetWritable`, mutable `getAllSets`, `setBoundarySets`, mask markers | sets and boundary mask | mesh setup; mutable getters currently have no external writer match except mesh | setup-time writer 合法 | P2，任何 consumer 可改 topology metadata |
| `CellFlag`, `IBMGhostLayer`, `setIBMFluidMask`, `IBMSignedDistance`, `setIBMGeometry` | IBM classification/geometry | CompositeIB, WeightBuilder, mesh import | IBM/mesh preprocessing | P1，Feature-specific data暴露给所有 callers |
| `setCommunicationHalo`, `setGlobalDofOwnership` | distributed identity | mesh/communication plan | 仅 plan setup | P1，影响 COPY/SUM owner |
| `ConvectiveFlux`, `ResX/Y/Z`, `Source`, stage/global residual API, clear APIs | stage workspace | face/viscous/source assembly、Runtime | RHS/Runtime only | P1，workspace 与 physical Q 混在 public Field surface |

建议不是为每个标量增加 setter。科学计算 kernel 需要高效 direct access；下一阶段应收缩
**owner namespace and aggregate**，不是把 28 个引用入口机械包装成 100 个 setter。

## 9. Dependency Hotspots

直接 include `core/field/SF_field.h` 或当前同目录 `SF_field.h` 的 source/header 文件为
**85** 个（85 个 include directives；同名 IO header 仍使裸 include 有歧义）。按目录统计：

| Subsystem | Direct includer files | Main use |
| --- | ---: | --- |
| `solver/algorithm` | 2 | lifecycle state access |
| `solver/equation` | 1 | compressible equation implementation |
| `solver/discretization` | 7 | CFL, WENO/div, diffusion/source, vector calculus |
| `solver/boundary` | 12 | BC and ILW reconstruction |
| `models` total | 27 | 其中 IBM 11；physics/turbulence/level-set 16 |
| `methods` | 16 | legacy Field-aware flux/time/viscous helpers |
| `infrastructure` | 8 | mesh, MPI, VTK, communication plan |
| `app` | 3 | composition and output |
| `core` | 8 | state views and interfaces |

高热区不是单纯 `solver/algorithm`，而是 boundary/discretization、IBM/model 和 mesh/MPI
共同经由一个 header 接触不同类别的数据。最小低风险切点因而是 residual/face-flux
workspace：它已有 `Residual`、`FluxField` 两个独立 storage type，直接 reader/writer
集中于 RHS、structured assembly、runtime，远少于 geometry 或 Q。不要把此结论用于
顺手处理 methods/viscous/scalar 或现有 10 条 allowlist debt。

## 10. Large Implementation Files

阈值按当前文本行数 `>=800`。长度不是数值错误；表中的 implementation-unit 建议均保持
现有 public abstraction，不增加 Manager/Adapter/Interface。

| File | Lines | Current internal stages | Risk |
| --- | ---: | --- | --- |
| `src/infrastructure/mesh/communication/SF_communicationPlan.cpp` | 2241 | topology discovery、donor spatial search、trilinear inversion、Lagrange/tensor weights、conformal map、dual volume、GlobalPoint/GlobalDof owner、canonical-face plan、validation | P1；owner/index/order 改动可改变 MPI/numerics |
| `src/models/ibm/method/ghost/SF_ilwClosure.cpp` | 2222 | geometry/frame、samples、2D/3D WENO fit/smoothness、precomputed plan、characteristic/wall first derivative、higher derivative、Taylor ghost、diagnostics | P1；ILW order/weights/failure behavior 敏感 |
| `src/infrastructure/mpi/SF_haloExchange.cpp` | 1510 | neighbor payload, Q/scalar pack/unpack、direct/tensor halo, GlobalDof COPY、canonical face flux COPY and residual assembly | P1；collective order、COPY/SUM 不可变 |
| `src/infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.cpp` | 1354 | SFM parsing、set normalization、topology validation、metrics/source canonical data、assembly/load | P2 |
| `src/infrastructure/mesh/createMesh/SF_meshGen.cpp` | 954 | structured generation、global-point dedup、boundary/interface report、curved edge/TFI、output | P2 |
| `src/infrastructure/mesh/decompose/SF_meshDecompose.cpp` | 832 | extrusion decision、logical cuts、source slicing、coverage validation、composite partition、coordinate COPY | P2 |
| `src/models/physics/multiphase/SF_multiphase.cpp` | 863 | config/init、thermal Euler、mixture properties、mass/alpha/T closure、conserved init、RHS/commit/BC/source/diagnostics | P1；time split/derived-state ordering 敏感 |
| `src/models/ibm/method/ghost/SF_weightBuilder.h` | 930 | SDF/classification、ghost layers、sample cache、ILW plan and donor weights | P2；header implementation exposes IBM complexity |
| `src/methods/numerics/scalar/SF_scalarTransport.h` | 802 | scalar BC/ILW、bounds、convection/diffusion RHS、Euler update | P1；deferred ownership issue |
| `src/core/config/SF_config.h` | 987 | cross-module configuration, parse/validation/defaults | P2 |
| GUI: `SF_widget.cpp`, `SF_controller.cpp`, `SF_project.cpp`, `SF_mainWindowTree.cpp`, `SF_editor.cpp`, `SF_configDocument.cpp` | 2053/1529/1300/1221/933/910 | view/controller/project/text parsing stages | P2；不属于本 solver data phase |

`SF_communicationPlan.cpp` 仍同时承担用户要求列出的全部 communication-plan 阶段；如未来
仅拆实现，应按“donor/interpolation”、“stable interface ownership”、“validation/reporting”
分 translation unit，保留一个 `CommunicationPlan` public abstraction。`SF_ilwClosure.cpp`
可按“precompute plan”与“runtime derivative/Taylor closure”划分实现单元；数学 kernel 和
IBM orchestration 的移动必须是后续独立决策。

## 11. Candidate Minimal Splits Option A/B/C

### Option A — Extract SolverWorkspace

将现有 `Residual` 与 `FluxField` 从 `Field` 成员移为 solver-owned patch workspace，
保留 Field 的 Q、geometry、stable metadata 和现有 numerical order。复用已有两个 storage
type；只需一个普通 data aggregate，不创建 Manager/Adapter/Context。预计触及 `Field`、
`DensityBasedRHS`、structured assembly、MPI canonical flux/GlobalDof assembly、部分
viscous/source readers，约 15–25 个直接实现/头调用点。

### Option B — Separate Structured Geometry from Flow State

在 Option A 之外，将 `GeometryStorage`、boundary set/mask、communication/GlobalDof
metadata 提升为 patch geometry/metadata owner，Flow Field 只留 Q/EOS/cache。会影响约
85 个直接 includer 中的大多数，特别是 mesh、MPI、boundary、models 和 output；需要
迁移 28 个 direct geometry references and view contracts。它可显著降低长期耦合，但这不是
小步重构。

### Option C — IBM-first extraction

将 `IBMGeometryStorage`、`cellType`、ghost layer 的 ownership 提升到 IBM patch state，
并让 ghost/forcing/constraint 通过明确 patch association 读取。可去除无 IBM case 的
15-array allocation，却触及 11 个 IBM header includer、mesh preclassification、turbulence
wall distance、output、halo availability 和 moving-body语义。需要先区分 static
preclassification、Ghost/ILW cache 和 target-time body state；不能合并这些语义。

## 12. Risk Comparison

| Option | Files / caller migration | Public API impact | Numerical risk | MPI risk | IBM risk | Reduced Field responsibilities |
| --- | --- | --- | --- | --- | --- | --- |
| A | low–medium；workspace readers concentrated | medium，Field workspace access迁移 | P1；必须保留 face barrier、residual sign、clear timing | P1；canonical COPY/GlobalDof SUM | P2，ghost仍使用同一 RHS order | H/I，及 Field 的工作区混合职责 |
| B | high；多数 Field includer | high | P1；metric/halo freshness/CFL sensitive | P1 | P1，IBM大量读 geometry | A/B/F/J 与 C/D 分离 |
| C | medium–high；IBM、mesh、models | high | P1；Ghost/ILW geometry and wall state | P1；partition halo filtering | P0/P1，静态/运动/forcing/constraint 分别验证 | K/L，消除 feature-specific core allocation |

Option A 风险最低，不因它是最小就可随意移动：它仍必须以 Phase 3A 的 canonical-face
barrier 和 GlobalDof SUM contract 为不可变前置条件。Option C 的内存收益明显，但目前
moving ghost cache 的更新时间语义尚未证明，不应先做。

## 13. Recommended Phase 4B

唯一建议：**Phase 4B — Extract SolverWorkspace（先抽取现有 `Residual` 与 `FluxField`）**。

理由是这两项已经有独立 storage implementation、数学寿命明确只到 RHS stage、读写者
集中，并且不改变 `StateBundle`/Field 的 authoritative Q、EOS binding、mesh geometry
或 IBM patch state。实施前必须把每个 API 的 owner 定为 RHS/discretization/runtime，保持：

```text
candidate convective flux for all patches
→ one canonical-face COPY/finalize
→ residual assembly
→ local residual stage
→ one GlobalDof SUM
→ explicit update
```

不得借此统一/改写 WENO、Steger-Warming、RK、boundary/halo timing、MPI collective、
IBM correction 或 storage layout。也不建议新增 “WorkspaceManager” 或适配器层。

## 14. Explicitly Deferred Work

- 完整重写或一次性拆分 Field；
- mesh geometry / physical state split（Option B）；
- IBM geometry/runtime extraction、Ghost/ILW、CompositeIBM、DFM/DLM/KKT（Option C）；
- single/multi density lifecycle（已在 Phase 3C 收敛）之外的 pressure-based lifecycle；
- methods/numerics time ownership、scalar transport、`SF_viscous`、model derivative ownership；
- phase boundary、level-set HJ-WENO、multiphase derived-state ownership；
- distributed monolithic KKT、HYPRE backend、Equation DSL/AssemblyPlan；
- architecture allowlist reduction。

自检：

| Question | Result |
| --- | --- |
| 是否修改源码？ | No |
| 是否提出整体重写 Field？ | No |
| 是否新增 Manager/Adapter/Context？ | No |
| 是否把 algorithm workspace 误认为 authoritative state？ | No |
| 是否混淆 MPI canonical ownership 与 Field physical-state ownership？ | No |
| 是否把 IBM-specific state 默认当作 core Field state？ | No |
| 是否基于旧代码而不是当前源码下结论？ | No |
