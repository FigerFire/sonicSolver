# CaseConfig Ownership Audit (Gate B)

对应阶段的 §11。对 `CaseConfig` 的**每一个**字段给出：

```text
字段 / 当前 writer / 当前 readers / 未来 owner / 移除条件
```

并标注它属于哪一层：

* **Equation** — 描述"解什么方程"
* **Model** — 描述"哪些物理模型参与"
* **Coupling** — 描述"相互依赖的方程如何达到一致"
* **Numerics** — 描述"怎么离散"
* **Solve** — 描述"怎么推进/求解"
* **Runtime** — 描述"运行期并行与后端能力"
* **IO-Output** — 描述"写什么结果、写到哪"
* **Legacy** — 过渡期兼容结构，有明确移除条件

`CaseConfig` 的定义在 [SF_caseConfig.h](../src/infrastructure/io/SF_caseConfig.h)。

---

## 1. Writer 汇总

`CaseConfig` 只有两个 writer，都是 `CaseAdapter` 的 member function：

| writer | 位置 | 写入字段 |
|--------|------|----------|
| `CaseAdapter::build()`（经 `SF_nativeBinding.cpp`） | `compatibility/SF_nativeBinding.cpp:214-217` | `caseName`, `createMesh`, `meshParameterFile`, `modelDescription`, `composition.instances` |
| `CaseAdapter::decodeNativeCase()`（经 `SF_nativeDecode.cpp`） | `compatibility/SF_nativeDecode.cpp:920-1090` | `time`, `parallel`, `writeInitial` |
| `CaseAdapter::buildCaseConfig()` | `compatibility/SF_compatibility.cpp:73-111` | 其余全部字段 |
| `CaseAdapter::resultWriterConfig()` | `compatibility/SF_compatibility.cpp:113-125` | `output` 子对象 |

`CaseConfig` 一旦由 `CaseAdapter` 交出即视为只读；生产代码里不存在任何第二个写入点。

---

## 2. 字段审计

### 2.1 身份与磁盘位置（IO-Output）

| 字段 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|------|-------------|--------------|------------|----------|
| `modelDescription` | `SF_nativeBinding.cpp:216`、`SF_model.cpp:53` | `SF_output.cpp:27`、`SF_cli.cpp:320` | Model | 不需要：这是 `Model::Description` 的持有者，是模型层的真实输入 |
| `caseDir` | `buildCaseConfig():74` | `SF_environment.cpp:118,217`、`SF_mesh.cpp:109,133,170`、`SF_vtkOutput.cpp`(6)、`SF_compatibility.cpp:115` | IO-Output | 不需要：磁盘位置不是数值语义 |
| `caseName` | `nativeBinding:214`、`buildCaseConfig():75` | `SF_cli.cpp:319`、`SF_vtkOutput.cpp`(7)、`SF_vtkIndexWriter.cpp:300` | IO-Output | 不需要：输出命名 |
| `jobName` | `buildCaseConfig():76`、`resultWriterConfig():117` | `SF_vtkOutput.cpp`(6)、`SF_vtkIndexWriter.cpp:300` | IO-Output | 不需要：输出命名 |

### 2.2 网格与几何输入（IO-Output / Model）

| 字段 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|------|-------------|--------------|------------|----------|
| `meshParameterFile` | `nativeBinding:215`（仅 `createMesh` 为真时）、`buildCaseConfig():77` | `SF_mesh.cpp:110` | IO-Output（mesh generator 格式） | `createMesh` 路径被独立的 mesh 生成入口取代 |
| `meshFiles` | `buildCaseConfig():78` | `SF_mesh.cpp:128,130,170` | IO-Output | 不需要：既有网格文件清单 |
| `ibmGeometryFiles` | `buildCaseConfig():79` | `SF_environment.cpp:119,215,218` | Model（IBM 几何） | IBM 几何改由 `Model::ObjectDescriptor` 的 typed 参数携带 |

### 2.3 方程与物理系统（Equation / Model）

| 字段 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|------|-------------|--------------|------------|----------|
| `composition` | `nativeBinding:217`（`instances`）、`buildCaseConfig()` 间接 | `SF_systemBuilder.cpp:515-562`、`SF_inspection.cpp:116-142`、`SF_model.cpp:54` | Equation | `EquationSystemInstanceConfig` 直接由 `Model::Description` 驱动 `EquationSystem`，不再经 `CaseConfig` 中转 |
| `multiPhase` | `buildCaseConfig():80` | `SF_inspection.cpp:38-105`、`SF_runtimeConfig.cpp:122`、`SF_compatibility.cpp:80` | Model（Physical/phase system） | 相系统改由 typed `phaseSystem` section 直接构造，`MultiPhaseConfig` 不再作为 `CaseConfig` 成员 |
| `multiPhaseEnabled` | `buildCaseConfig():81` | `SF_inspection.cpp:38,87,97,100`、`SF_runtimeConfig.cpp:120` | Model | 与 `multiPhase` 同批移除。"启用与否"应表达为 block 是否存在（§19），而不是并列 bool |

### 2.4 数值策略（Numerics / Solve / Coupling）

`solver`（`FDM::SolverConfig`）是当前最大的复合字段，71 个 reader。它同时承载
Numerics / Solve / Coupling 三类语义，因此必须按子对象而不是整体来看：

| 子对象 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|--------|-------------|--------------|------------|----------|
| `solver.numerics.solver` | `buildCaseConfig():101`（仅 `solverPropertiesLoaded_`）、`Legacy::makeSolverConfig()` | `SF_environment.cpp:264-266`、`SF_inspection.cpp:108` | Coupling | 由 `RuntimeRequirements`/typed section 直接给出；`SF_environment.cpp` 已禁止再用它推断 runtime service |
| `solver.numerics.startTime` | `buildCaseConfig():86` | Numerics 消费者 | Numerics | 与 `RunControl` 统一为唯一时间 authority |
| `solver.pressure.workflow` / `.relaxation` / `.velocityRelaxation` | `buildCaseConfig():102-107` | `SF_pressureStepper` 路径 | Coupling | pressure-correction 关系由 equation/constraint 声明后移除 |
| `solver.ibm.method` / `.ibm.forcing` | `buildCaseConfig():87-88` | IBM 路径 | Model（IBM 语义） | IBM 方法按 §12 分类完后由 typed section 直接构造 |
| `solver.sources.wallHeat` / `.sources.enabled` | `buildCaseConfig():89-99` | source 装配 | Model（source 贡献） | wall-heat 作为 typed source 声明后移除 |
| `solver` 其余成员 | `Legacy::makeSolverConfig()` | 见 §4 | 各类 | `Legacy::makeSolverConfig()` 骨架被 typed 构造取代后整体移除 |

| 字段 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|------|-------------|--------------|------------|----------|
| `time`（`Time::RunControl`） | `SF_nativeDecode.cpp:920-955`（`decodeRuntimeSettings`） | `SF_application.cpp:62-63`、`SF_compatibility.cpp:84` | Solve | 不需要：`RunControl` 是时间推进的既有 typed 值对象 |
| `createMesh` | `nativeBinding:214`、`buildCaseConfig():109` | `SF_application.cpp:72`、`SF_compatibility.cpp:109` | IO-Output | 网格生成变为独立命令而不是 run 的一个 bool |

### 2.5 运行期并行（Runtime）

| 字段 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|------|-------------|--------------|------------|----------|
| `parallel`（`CaseParallelConfig`） | `SF_nativeDecode.cpp:1081-1090`（`decodeRuntimeSettings`） | `SF_environment.cpp:63-146`（51 处） | Runtime | 不需要：MPI 分区是运行期输入，且已有独立 typed 值对象。注意它描述的是**分区**，不是"是否并行"这一 solver 身份 |

### 2.6 输出（IO-Output）

| 字段 | 当前 writer | 当前 readers | 未来 owner | 移除条件 |
|------|-------------|--------------|------------|----------|
| `writeInitial` | `SF_nativeDecode.cpp:933-935` | `SF_report.cpp:121`、`SF_singleFluid.cpp:437`、`SF_eulerian.cpp:101`、`SF_multiPatch.cpp:273` | IO-Output | 不需要：纯输出策略 |
| `output`（`ResultWriterConfig`） | `resultWriterConfig():113-125` | `SF_application.cpp:67`、`SF_output.cpp:40-82` | IO-Output | 不需要：writer 的 typed 配置 |

---

## 3. 不得新增的字段形态（§12）

`CaseConfig` 与 `BuildRequest` **不得**再新增 `bool xxxEnabled` 形式。

现状核对：

* `CaseConfig` 现有两个此类 flag：`multiPhaseEnabled`（§2.3，已有移除条件）与
  `createMesh`（§2.4，描述命令而不是模型启用）。本次**没有**新增任何 flag。
* `BuildRequest` 侧在 §22/§23 中反而**删除了**两个 bool
  （`constraintGlobalDofAvailable`、`distributedLinearSystemAvailable`，见 §4）。
* 可选模型的启用语义统一由 §19 表达：

  ```text
  block absent                -> disabled
  block present               -> enabled by default
  block present + enabled:false -> configured but inactive
  ```

---

## 4. §22 / §23：BuildRequest 字段审计

`BuildRequest` 的唯一生产者是
[SF_inspection.cpp](../src/app/application/SF_inspection.cpp)（:83-143）。
逐字段审计结果：

| 字段 | 判定 | 处置 |
|------|------|------|
| `constraintGlobalDofAvailable` | **不是** case 属性：唯一生产者硬编码 `true`；语义是"后端是否具备该能力" | 删除；改为 `SF_systemBuilder.cpp` 内文件作用域 `constexpr bool kConstraintGlobalDofAvailable = true;` + WHY 注释 |
| `distributedLinearSystemAvailable` | 同上 | 删除；改为 `constexpr bool kDistributedLinearSystemAvailable = true;` |
| `templateOrigin` | 是 `CaseConfig`/`composition` 的投影，且是 builder 的真实输入 | 保留 |
| `singleFluidPreset` / `phaseNames` / `levelSet` / `homogeneousThermodynamics` / `legacyMixture` / `phaseChange` / `transportedLegacyAlpha` / `turbulence` / `parallel` / `composition` / `immersed` | 都是 builder 的真实判定输入，尚无 typed owner | 保留 |

删除依据可证：原表达式为 `!constraintDof || true` 与 `true`，替换为常量后求值结果逐位相同，
因此是行为保持的删除。

---

## 5. §21：explain-only 字段的分离

`ResolvedSimulationSystem` 上有三个字段**只被 printer / CLI / 测试读取**，从不参与数值传播：

| 字段 | 唯一读者 | 处置 |
|------|----------|------|
| `templateOrigin` | `SF_systemPrinter.cpp:487-489`、`SF_cli.cpp:306` | 移入 `CaseClassification` |
| `densityBehavior` | 同上、`test_pisoArchitecture.cpp:460-461` | 移入 `CaseClassification` |
| `thermodynamicCompressibility` | 同上、`test_pisoArchitecture.cpp:477-478` | 移入 `CaseClassification` |

分离后 `ResolvedSimulationSystem` 的 loose 字段数从 12 降到 10，新类型
[SF_caseClassification.h](../src/solver/system/SF_caseClassification.h) 用注释明确声明
"EXPLAIN-ONLY，不是数值 authority"。

**注意**：`ResolvedSimulationSystem::formulation` **没有**被移动。它是 builder 的真实
输入（`SF_systemBuilder.cpp:656` 把它传给 `NumericalCompiler::compile`），因此不属于
§21 的 explain-only 范畴。

---

## 6. 结论

* `CaseConfig` 的每个字段都有唯一 writer 和明确 owner 分类（§2）。
* 本次没有新增任何 `bool xxxEnabled`（§12）。
* `BuildRequest` 删除了 2 个可由后端能力常量替代的字段（§22/§23）。
* 3 个 explain-only 字段被隔离到专门的类型（§21）。
* 仍标记为 **Legacy** 的字段均有可执行的移除条件，且移除条件都取决于 Equation/Model
  typed 化的完成度，而不是取决于"看起来是否像 OpenFOAM"。
