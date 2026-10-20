# SonicSolver Architecture Weight Audit

日期：2026-09-21  
状态：Architecture Contraction 修改前基线

## 1. 统计口径

统计范围为 `src/` 下的 `.h/.hpp/.cpp/.cc/.cxx`。LOC 是物理行数，只用于定位重量；最终分类同时考虑职责数量、依赖方向、public exposure、运行时 authority 和算法内聚性。

## 2. Directory LOC

| Directory | LOC |
|---|---:|
| `src/models` | 24,949 |
| `src/app` | 20,622 |
| `src/solver` | 18,368 |
| `src/infrastructure` | 15,945 |
| `src/methods` | 5,530 |
| `src/core` | 4,670 |

Public/source headers：371。

## 3. Top 30 production files by LOC

| LOC | File | Classification |
|---:|---|---|
| 2241 | `infrastructure/mesh/communication/SF_communicationPlan.cpp` | Algorithm-heavy |
| 2234 | `models/ibm/method/ghost/SF_ilwClosure.cpp` | Algorithm-heavy |
| 2053 | `app/gui/VTKView/SF_widget.cpp` | Compatibility/UI-heavy |
| 1529 | `app/gui/Controller/SF_controller.cpp` | Compatibility/UI-heavy |
| 1519 | `infrastructure/mpi/SF_haloExchange.cpp` | Algorithm-heavy |
| 1337 | `infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.cpp` | Algorithm-heavy |
| 1300 | `app/gui/IO/SF_project.cpp` | Compatibility/UI-heavy |
| 1221 | `app/gui/GUI/SF_mainWindowTree.cpp` | Compatibility/UI-heavy |
| 1213 | `solver/system/SF_systemBuilder.cpp` | Architecture-heavy |
| 954 | `infrastructure/mesh/createMesh/SF_meshGen.cpp` | Algorithm-heavy |
| 941 | `models/ibm/method/ghost/SF_weightBuilder.h` | Algorithm-heavy |
| 933 | `app/gui/Geometry/SF_editor.cpp` | Compatibility/UI-heavy |
| 910 | `app/gui/IO/SF_configDocument.cpp` | Compatibility/UI-heavy |
| 863 | `models/physics/multiphase/SF_multiphase.cpp` | Algorithm-heavy |
| 832 | `infrastructure/mesh/decompose/SF_meshDecompose.cpp` | Algorithm-heavy |
| 802 | `methods/numerics/scalar/SF_scalarTransport.h` | Algorithm-heavy / ownership debt |
| 771 | `app/application/model/SF_configParser.h` | Compatibility-heavy |
| 648 | `app/gui/GUI/SF_GUI.cpp` | Compatibility/UI-heavy |
| 621 | `solver/boundary/SF_thermodynamicClosure.cpp` | Algorithm-heavy |
| 619 | `solver/algorithm/pressureBased/eulerian/SF_pressureStepper.cpp` | Algorithm-heavy / transitional |
| 616 | `infrastructure/mesh/SF_mesh.cpp` | Algorithm-heavy |
| 581 | `models/physics/equationSet/SF_homogeneousMultiphaseEquationSet.h` | Algorithm-heavy |
| 578 | `solver/algorithm/pressureBased/SF_corrector.cpp` | Algorithm-heavy |
| 568 | `solver/system/SF_solvePlan.cpp` | Healthy, dense compiler logic |
| 563 | `solver/system/SF_systemPrinter.cpp` | Healthy, explain rendering |
| 547 | `solver/algorithm/SF_compressible.cpp` | Algorithm-heavy |
| 533 | `solver/boundary/reconstruction/ILW/SF_boundaryClosure.cpp` | Algorithm-heavy |
| 531 | `app/gui/ViewerControls/SF_toolbar.cpp` | Compatibility/UI-heavy |
| 521 | `models/ibm/SF_compositeIBM.cpp` | Algorithm-heavy |
| 510 | `app/cli/SF_cli.cpp` | Application-heavy |

`SF_communicationPlan.cpp`、`SF_ilwClosure.cpp`、`SF_haloExchange.cpp` 和 `SF_MultiBlockMesh.cpp` 很大，但内部复杂度来自 donor search、ILW、MPI ownership 和 mesh topology。本阶段不以 LOC 为由改变其 public architecture。

## 4. Key files

| File | LOC | Finding |
|---|---:|---|
| `solver/system/SF_systemBuilder.cpp` | 1213 | 同时拥有物理方程、model contribution、transformation、planning、runtime inference；Architecture-heavy |
| `solver/system/SF_resolvedSimulationSystem.h` | 459 | 将 raw/executable/numerical/plan/runtime/BuildRequest 暴露在同一 public header；Architecture-heavy |
| `app/application/model/SF_configParser.h` | 771 | time/term/pressure/IBM/boundary/source parser 集中；Compatibility-heavy |
| `app/application/execution/SF_singleFluid.cpp` | 468 | conservative composition 生命周期与 multi 重复；Architecture-heavy |
| `app/application/execution/SF_multiPatch.cpp` | 328 | conservative composition 生命周期与 single 重复；Architecture-heavy |
| `core/interfaces/SF_interfaces.h` | 222 | umbrella 聚合多个无关 port；Architecture-heavy |
| `core/interfaces/SF_executionRuntime.h` | 234 | contract 与 `LocalExecutionRuntime` concrete implementation 混合；Architecture-heavy |
| `app/application/model/compatibility/*` | 4501 | legacy reader 合理存在，但 native path 也回流此层；Compatibility-heavy |

## 5. Include fan-in

最高 fan-in：

| Includes | Header |
|---:|---|
| 32 | `core/config/SF_config.h` |
| 29 | `core/interfaces/SF_log.h` |
| 29 | `methods/numerics/structured/SF_structured.h` |
| 27 | `models/physics/phaseChange/SF_phaseChange.h` |
| 26 | `core/state/SF_valueTypes.h` |
| 26 | `core/mesh/SF_dimension.h` |
| 24 | `core/interfaces/SF_interfaces.h` |
| 19 | `core/field/SF_field.h` |
| 16 | `infrastructure/mesh/MultiBlockMesh/SF_MultiBlockMesh.h` |
| 16 | `infrastructure/io/private/SF_foamParser.h` |
| 15 | `solver/system/SF_resolvedSimulationSystem.h` |
| 13 | `core/state/SF_equationSet.h` |
| 12 | `core/config/SF_configTypes.h` |
| 9 | `core/interfaces/SF_executionRuntime.h` |

`SF_interfaces.h` 的 basename include 统计为 24；production 源中显式 `#include "SF_interfaces.h"` 为 22。它是当前最明显的非语义 compile coupling。

## 6. Include fan-out

除 GUI 外，最高 fan-out 的 production 文件为：

| Includes | File |
|---:|---|
| 30 | `app/application/model/compatibility/SF_compatibility.cpp` |
| 29 | `app/application/execution/SF_singleFluid.cpp` |
| 27 | `app/application/execution/SF_multiPatch.cpp` |
| 19 | `models/ibm/method/ghost/SF_ilwClosure.cpp` |
| 19 | `app/application/model/compatibility/SF_sourceReader.cpp` |
| 18 | `infrastructure/mpi/SF_haloExchange.cpp` |
| 18 | `infrastructure/mesh/communication/SF_communicationPlan.cpp` |
| 18 | `solver/algorithm/SF_compressible.h/.cpp` |
| 18 | `app/application/model/compatibility/SF_controlReader.cpp` |
| 18 | `app/application/model/compatibility/SF_fieldReader.cpp` |
| 18 | `app/application/model/compatibility/SF_thermoReader.cpp` |
| 16 | `app/application/execution/SF_eulerian.cpp` |

Single/multi conservative entry points 的高 fan-out 来自两份重复 composition，而不是数值算法本身。

## 7. Authority consumer counts

文本 consumer 基线：

| Authority/type | Production files |
|---|---:|
| `CaseConfig` | 35 |
| `SolverConfig` | 38 |
| `ResolvedSimulationSystem` | 21 |
| `SF_interfaces.h` explicit production includes | 22 |

### CaseConfig

消费分布于 application inspection/execution/output、native/legacy binding、mesh construction。主要问题是 application execution 收到完整 `CaseConfig`，即使只需要 runtime/output/field 的局部配置。

### SolverConfig

被 application、density/pressure algorithm、linear algebra 和 compatibility reader 同时消费。它仍是过宽的 aggregate configuration boundary。

### ResolvedSimulationSystem

被 inspection、single/multi/eulerian execution、density/pressure algorithms、state realizer、printer、validator 消费。它目前是 transitional aggregate，不应继续增加字段。

## 8. BuildRequest ownership baseline

| Field | Current owner/producer | Consumer | Future owner | Removal condition |
|---|---|---|---|---|
| `templateOrigin` | application inspection | builder/report | legacy provenance | all native instances selected by object type |
| `singleFluidPreset` | application inspection | single-fluid preset | EquationSystemInstance | instance reads/installs itself |
| `phaseNames` | application inspection | Eulerian pack | coupling | independent phase instances plus coupling resolution |
| `levelSet` | application inspection | builder | LevelSet model | model contribution installed directly |
| `homogeneousThermodynamics` | application inspection | builder/runtime inference | equation instance | homogeneous instance owns closure |
| `legacyMixture` | compatibility/application | builder | compatibility bridge | legacy mixture translated to typed contribution |
| `phaseChange` | application inspection | builder | phase-change model | model contributes requirement/equations |
| `transportedLegacyAlpha` | application inspection | builder | legacy mixture model | legacy alpha contribution object |
| `interfaceGhostFluid` | application inspection | builder | interface coupling | LevelSet/GFM contribution installed directly |
| `turbulence` | application inspection | builder | turbulence model | model contribution installed directly |
| `turbulenceModel` | application inspection | builder | turbulence model | typed model descriptor contributes closure/equations |
| `turbulencePhaseNames` | application inspection | builder | turbulence/coupling | per-instance turbulence binding |
| `parallel` | application inspection | runtime inference | runtime capability | infer from ExecutionRuntime resources |
| `constraintGlobalDofAvailable` | application inspection | validator | runtime capability | capability set supplied separately |
| `distributedLinearSystemAvailable` | application inspection | validator | runtime capability | capability set supplied separately |
| `composition` | native/legacy binding | builder | equation instances/coupling | objects read typed sections |
| `immersed` | application inspection | builder | IBM model | descriptor contributes directly |

## 9. Classification summary

### Architecture-heavy

- `SF_systemBuilder.cpp`
- `SF_resolvedSimulationSystem.h`
- `SF_singleFluid.cpp` / `SF_multiPatch.cpp`
- `SF_inspection.cpp`
- `SF_interfaces.h`
- `SF_executionRuntime.h`

### Compatibility-heavy

- `app/application/model/compatibility/*`
- `app/application/model/SF_configParser.h`
- native `Model::Description -> in-memory dictionary -> readCase()` bridge

### Algorithm-heavy

- communication plan, halo exchange, mesh topology
- IBM ILW/weight construction
- pressure corrector and compressible RHS
- scalar transport kernels

### Healthy or bounded

- `solver/run/SF_planExecutor.*`: consumes `CompiledSolvePlan + OpRegistry` and does not select physics
- `solver/system/SF_numericalCompiler.*`: lowers executable terms using typed recipes
- `solver/system/SF_solvePlan.*`: compiles structured control flow
- `core/state/SF_stateBundle.h`: physical state/storage binding authority

## 10. Modification priorities

1. Move source, turbulence, LevelSet and IBM mathematical metadata out of the giant builder to their real model owners.
2. Remove umbrella interface includes and move local runtime implementation to infrastructure.
3. Share conservative execution composition between one and many patches.
4. Let native typed objects bypass dictionary reserialization where their readers are implemented.
5. Reduce BuildRequest fields only after their new owners are active and tested.

`solver/run` and the large algorithm-heavy implementation files are frozen for this phase.
