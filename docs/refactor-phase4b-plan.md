# Phase 4B 计划：从 Field 分离 Residual / Face-Flux Workspace

日期：2026-09-14。本计划在任何源码修改前依据当前源码编写。范围只包含
`Field::convectiveFlux_`、`Field::residual_` 和直接读写它们的 numerical data path。

## 1. 修改层、调用链与边界

本修改位于 `solver/algorithm → solver/equation → solver/discretization` 的数值数据
通路。现状是 `Field` 既拥有 Q/geometry/stable metadata，又拥有一个 RHS stage 的
`FluxField` 与 `Residual`。`StateBundle` 只组合 patch 指针、EOS、时钟、registry 和
non-owning distributed views；`ExecutionRuntime` 只调度 Runtime contract。

拟修改的直接调用链为：

```text
CompressibleAlgorithm::bindState
  → one vector<PatchWorkspace>, index-aligned with StateBundle::patches
  → DensityBasedTime / legacy ddt callback
  → DensityBasedRHS::assembleAllPatches(fields, workspaces)
  → Equation::Compressible::System
  → convection / diffusion / source assembly
  → FluxField + Residual
  → one canonical-face COPY barrier / one GlobalDof SUM barrier
  → explicit time update reads Residual
```

不新增 Manager、Context、Adapter 或 service。若需新增类型，只新增一个
solver/algorithm 侧的数据聚合 `PatchWorkspace { FluxField convectiveFlux; Residual
residual; }`；它不含 timestep、MPI、boundary、IBM 或 dispatch 行为。

## 2. 修改前基线

- `python3 tools/check_architecture.py`：574 source files、10 条 allowlisted edge。
- `cmake --build build --parallel 4` 已启动完整基线构建；本机 build tree 曾报告
  `ninja: premature end of file; recovering`，后续应在干净的增量构建中记录最终结果。
- `ctest --test-dir build-phase1a-tests --output-on-failure --parallel 1`：
  registryIO、surfaceSchurProjection、distributedConstraintLayout 通过；三个 MPI
  test 在 sandbox 中因 PRTE 无权绑定 socket 失败，不是测试断言失败。修改后将以同一
 方式复测，并在允许的执行环境运行 MPI regression。
- Phase 3C 数值基线：4-rank Sod 末尾 `time=2.435715e-01`、
  `dt=5.352234e-04`；single WENO7 SSPRK3、serial ghost IBM、velocity-forcing IBM
  与 RPI smoke 均有既有通过记录。

## 3. Complete Residual / Flux Call-Site Matrix（修改前）

| Caller | Operation | Flux / Residual | Density / Pressure / Common | Lifetime requirement |
| --- | --- | --- | --- | --- |
| `Field::setup/resizeConservedVariables` | allocation | `Residual::setup`、`FluxField::setup` | common | layout construction / explicit NVar resize |
| `Field::clearResidual` | clear residual, source, local/global, flux | both | common | every spatial assembly stage entry |
| `Equation::Compressible::System::begin` | calls clear | both | density and legacy pressure route | before convection |
| `RusanovEOS::div` | candidate flux + residual assembly | flux then directional residual | EquationSet Rusanov | one RHS stage |
| `divDispatch` / WENO/TENO | candidate flux + directional residual | both | high-order density | one RHS stage, before canonical barrier |
| `HaloExchange::assembleCanonicalInterfaceFluxes` | canonical owner flux COPY, then residual assembly | flux read/write then residual write | multi-patch / MPI | exactly once after all candidates |
| `Viscous::computeCentralViscousRHS` / `laplacianDispatch` | subtract diffusive face contribution | directional residual | common compressible | after canonical convection |
| `SourceTerm::Sp`, gravity, MRF, wall heat | source accumulation / source clear | residual source | common | after diffusion, before local residual stage |
| equation/interface coupling, multiphase, level-set, homogeneous phase change | model source accumulation | residual source | density coupling | same source accumulator; no separate IBM/model residual |
| `stageLocalSpatialResidual` | local strong residual snapshot | residual local/global mask reset | density | after all flux/source writes |
| `GlobalDofResidualAssembler` | SUM then write owner complete value | local/global residual | multi-patch / MPI | exactly once after all patch local stages |
| `spatialResidual`, `applyDivergence`, legacy Euler/RK4 | explicit update consumption | residual read | legacy pressure path | update after assembled RHS |
| `DensityBasedTime` | explicit tableau update | residual read | unified density Euler/SSPRK3/RK4 | each stage, no allocation |
| `IEquationSystemCoupling` implementations | source/coupling RHS | source residual | common density coupling | one residual per same patch index |
| IBM | conservative correction only in current path | no Field flux/residual reader found | Ghost/forcing/constraint | IBM must not gain a second residual |
| tests | no direct Field flux/residual test caller found | N/A | N/A | add structural test |

Pressure-based Eulerian solver uses its own phase/matrix workspace. The legacy single-field
`stepPressure` route invokes the compressible explicit assembly, so it must pass the same
`PatchWorkspace`; no pressure corrector mathematics or ordering is to change.

## 4. Residual semantic classification

`Residual::{resX,resY,resZ,source,local,global,globalMask}` are all RHS-stage storage:

| Member | Semantics | Cross-timestep physical history? | Stable mesh identity? | Decision |
| --- | --- | --- | --- | --- |
| `resX/resY/resZ` | directional assembled face flux cache | No | No | workspace |
| `source` | explicit source accumulator | No | No | workspace |
| `local` | local strong residual before distributed assembly | No | No | workspace |
| `global` | SUM result replicated to participating local points | No | No | workspace |
| `globalMask` | marks validity of current `global` result | No | No; identity remains Field GlobalDof metadata | workspace |

`FluxField::values_` contains only candidate/canonicalized numerical `F*`. Canonical face
metric/mask/owner arrays remain in `Field::GeometryStorage` and are not moved.

## 5. Workspace-owner decision

| Candidate | Assessment | Decision |
| --- | --- | --- |
| A. `CompressibleAlgorithm` owns collection | correct lifecycle and patch association, but needs a per-patch aggregate | selected collection owner |
| B. equation-system owns | equation has no lifetime beyond one assembly and should not retain RK workspace | reject |
| C. existing solver workspace | no existing reusable per-patch Flux+Residual owner exists | unavailable |
| D. minimal `PatchWorkspace` aggregate | data-only, solver-side, maps exactly one patch to one flux/residual pair | selected representation |

`CompressibleAlgorithm` will own a contiguous `std::vector<PatchWorkspace>`. It is rebuilt only
when patch count changes and `ensureFor(Field)` allocates only when existing Flux/Residual shape
does not match current `TotalSize/NVar`. The vector index equals the participating patch index;
no `Field*` map or cell-loop lookup is introduced.

`Field`, `StateBundle`, `ExecutionRuntime`, and infrastructure do not own workspace storage.
`StateBundle::distributed` may retain non-owning Flux/Residual views only when runtime needs them.
Runtime consumes those views or explicit borrowed workspace parameters; it does not retain an owner.

## 6. Planned API migration

1. Remove `FluxField convectiveFlux_`, `Residual residual_`, all Field workspace setup/reset,
   and all Field flux/residual/source accessors.
2. Add `solver/algorithm/SF_patchWorkspace.h` with `ensureFor(const Field&)` and `clear()`;
   it stores only existing `FluxField` and `Residual` values.
3. Make `DensityBasedRHS::assembleAllPatches` accept an index-aligned workspace vector and
   validate count/shape before any operator runs.
4. Pass `FluxField&` / `Residual&` explicitly through `Equation::Compressible::System`,
   `RusanovEOS`, WENO/TENO dispatch, face assembly, viscous/laplacian, source assembly,
   residual assembly and explicit time functions. These types are core storage types, so equation
   and discretization do not include the algorithm aggregate.
5. Extend only the existing equation-coupling source-assembly call with aligned
   `std::vector<Residual*>`; migrate app coupling and source-producing models to write the same
   residual. This adds no core-to-solver dependency.
6. Change MPI halo / GlobalDof helpers to receive borrowed aligned Flux/Residual storage, or
   consume registered non-owning distributed views. Keep canonical metadata in Field and retain
   one all-patch barrier.
7. Migrate legacy `ddtDispatch`/Euler/RK4 signatures so their residual reader/reset is explicit;
   preserve all stage coefficients and call times.
8. Add a small structural test that creates N Fields and N workspaces, verifies layout association
   and verifies that Field no longer declares `FluxField`/`Residual`.

## 7. Invariants and non-goals

- allocation happens at solver composition/bind or explicit Field layout change, never per RHS
  stage; no Q, flux or residual copies;
- source clear remains at the same `System::begin` point; residual sign remains
  `Jac*divergence - Source`;
- all candidate fluxes complete before one canonical COPY; all local residuals complete before one
  GlobalDof SUM; no averaging;
- `thermodynamicCache_`, Field Q, geometry, boundary metadata, GlobalDof identity, IBM storage,
  `StateBundle` state ownership and pressure corrector semantics remain untouched;
- no architecture allowlist reduction is attempted; new `core -> solver` dependency must be zero.

## 8. Planned validation

After each migration slice: build and a targeted test. Final validation: architecture checker,
full build, CTest, four-rank Sod, WENO7 plus Euler/SSPRK3/RK4 representatives, serial ghost IBM,
serial forcing IBM, RPI wall-boiling smoke, and comparison of dt/time, conservation summaries,
state extrema and available norms. If output differs, investigate workspace association, reset
timing, source sign, face index, registration, canonical COPY and GlobalDof SUM before changing
tolerance.
