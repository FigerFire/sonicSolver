# Native IO Contraction Audit (Gate A)

对应阶段的 §3 Gate A。目标不是"文件变短"，而是证明原生语义 YAML 到强类型配置之间
**不存在** 一层 OpenFOAM 形状的中间文档；并把每一个 OpenFOAM 遗留结构按性质分类，
给出删除顺序与删除依据。

本审计只覆盖 case 输入侧（case → `CaseConfig`）。数值、离散、边界、MPI 语义不在范围内。

---

## 1. 真实调用链（收缩后）

```text
case.yaml
  │
  ├─ test/<case>/mesh/mesh.yaml          -> Description::mesh
  ├─ test/<case>/fields/fields.yaml      -> Description::fields(注册)
  ├─ test/<case>/fields/internalField.yaml
  ├─ test/<case>/fields/boundaries.yaml
  ├─ test/<case>/solvers/solvers.yaml    -> Description::solver
  ├─ test/<case>/solvers/algorithm.yaml  -> Description::solver 内 algorithm 段
  ├─ test/<case>/solvers/numerics.yaml   -> Description::numerics
  ├─ test/<case>/solvers/runtime.yaml    -> Description::runtime
  ├─ test/<case>/models/models.yaml      -> Description::objects
  └─ test/<case>/equations/equations.yaml / algorithms/algorithms.yaml（可选）
        │
        ▼
SF::CaseIO::read(path)                      infra/io/case/SF_case.cpp
        │   syntax -> Model::Description（语义对象，尚未解释）
        ▼
SF::Application::ModelLoader::build(model)  app/application/model/SF_model.cpp:42
        │   分离 built-in object 与 plugin object；收集 EquationSystemInstanceConfig
        ▼
SF::CaseAdapter::build(description)         compatibility/SF_nativeBinding.cpp
        │   语义对象 -> typed sections_（NativeCaseSections）
        │   field 语义名 -> Model::FieldDescriptor，boundaries 按 mesh patch set 展开
        │   不读任何字典内容，不构造 document path
        ▼
SF::CaseAdapter::decodeNativeCase()         compatibility/SF_nativeDecode.cpp
        │   typed section -> 数值语义值（CaseConfig 的 typed 成员）
        │   绝不按文件名或目录名推断语义
        ▼
SF::CaseConfig                              infra/io/SF_caseConfig.h
```

关键不变量：

* `CaseAdapter` 内部 **没有** `documents_` 这类"伪路径 -> 字典文本"的映射。
  唯一入口是 `NativeCaseSections` 的 typed slot 和 `fields`。
* `decodeNativeCase()` 的每个分支由 **语义对象是否存在** 决定（`hasAlgorithm`、
  `hasThermoDynamics`、…），而不是由某个 OpenFOAM 文件名是否存在决定。
* `field()` 只按语义名查找；`boundaries` 已在 typed 阶段按 mesh patch set 展开
  `default`，因此解码阶段不做目录名/文件路径推断。
* `CaseIO::read()` 不再有任何 OpenFOAM 分支；`ModelLoader::read()` 对没有
  `case.yaml` 的目录直接 fail fast（`SF_model.cpp:57`），不再回落到旧解析器。

---

## 2. 分类清单

分类含义（§3）：

| 类 | 含义 | 处置 |
|----|------|------|
| 1 | 真正的原生语义对象 | 保留，且必须是唯一的解码入口 |
| 2 | OpenFOAM 形状的内部兼容结构 | 迁移完成后删除 |
| 3 | 死兼容（无 reader / 不可达） | 立即删除 |
| 4 | 纯解析 helper | 保留，仅命名遗留 |

### 2.1 类 1：原生语义对象（保留）

| 结构 | 位置 | 说明 |
|------|------|------|
| `NativeCaseSections` | `compatibility/SF_compatibility.h:31` | 14 个 typed slot + `fields`，一个 native 语义对象一个 slot |
| `Model::Description` | `core/model/SF_model.h:39` | 格式中立的语义描述，IO 与执行共用 |
| `Model::FieldDescriptor` | `core/model/SF_model.h:21` | field 元数据；`storage`/`location`/`type` 已是强类型 |
| `Description::runtime/numerics/solver/mesh` | `core/model/SF_model.h:43-46` | 语义对象的 typed 容器 |
| `Model::FactoryRegistry<CaseConfig>` | `compatibility/SF_compatibility.h:11` | plugin object 的强类型构造入口，IO 不枚举模型参数名 |
| `CaseConfig` | `infrastructure/io/SF_caseConfig.h:38` | 解码产物，solver/mesh/writer 的唯一配置值对象 |

### 2.2 类 2：OpenFOAM 形状的内部兼容结构（待迁移后删除）

| 结构 | 位置 | 现状 | 删除条件 |
|------|------|------|----------|
| `Legacy::makeSolverConfig()` | `compatibility/SF_legacyConfig.h/.cpp` | 被 `SF_compatibility.cpp:84` 与 `SF_nativeDecode.cpp:12` 调用，产出 `FDM::SolverConfig` 的 legacy 默认骨架 | `FDM::SolverConfig` 的每个成员都能由 typed section 直接构造 |
| `SF::` legacy globals + `ParserState` | `models/initial/SF_parameter.h/.cpp`、`compatibility/SF_parserState.h` | `ParserState` 在事务里 save/restore 约 82 个全局变量（静态 mutex 保护） | `decodeNativeCase()` 的每个写入点都改成 typed 局部量 |
| `SF_foamParser.h/.cpp` | `infrastructure/io/private/` | 被 `SF_nativeDecode.cpp:16`、`SF_casePath.cpp:6`、`SF_meshDictReader.cpp:10` 引用；后两者是 mesh generator 格式，属合法 | `SF_nativeDecode.cpp` 不再需要任何字典 token 时，前一条依赖消失 |

> 注意：`SF_parameter.*` 与 `ParserState` 仍然承载**真实**的解码中转语义，不能
> 因为"名字像 OpenFOAM"就删除。它们属于类 2，删除条件是迁移完成，而不是命名。

### 2.3 类 3：死兼容（已删除）

| 结构 | 位置（删除前） | 删除依据 |
|------|----------------|----------|
| `documents_` 伪路径 -> 文档文本映射 | `SF_compatibility.h` | 被 `sections_` + `fields` 取代；typed 之外无第二入口 |
| `parseOpenFOAM*` / `parse*File` 系列 | `SF_compatibility.h/.cpp`、`SF_nativeDecode.cpp` | 全部重命名为 `decode*`；不再有"先解析成文档再解释"的阶段 |
| `isFoamControlDictPath` | `infrastructure/io/private/SF_casePath.h:15`、`.cpp:50` | 唯一调用者是已删除的 ctor 分支 |
| `isOpenFOAMCase()` | `SF_compatibility.cpp` | 删除后无任何调用点（曾误加回一份，已再次删除并复核） |
| `addIBMGeometryFile` 的 `constant/triSurface/` 前缀分支 | `SF_nativeDecode.cpp:~2295` | 不可达：typed 路径不会给出该前缀；改为 §22 明确 fail fast |
| dynamic mesh / moving boundary 警告路径 | `SF_nativeDecode.cpp` | typed section 无法表达该模型，路径不可达 |
| 6 个只写不读的 legacy 全局 (`ibmFiles`, `gravityFiles`, `mrfFiles`, `wallHeatSourceFiles`, `dynamicMeshFiles`, `movingBoundaryFiles`) | `SF_parameter.h/.cpp`、`SF_parserState.h` | 全仓 grep 证明 0 个 reader |
| `meshParamFile_` | `SF_compatibility.h` | 语义就是 mesh generator 路径，改名 `meshGenerator_` |

删除这 8 项之后，`check_program_decomposition.py --mode scope` 会持续守住
`parseOpenFOAM` / `documents_` / `isFoamControlDictPath` 不再出现。

### 2.4 类 4：纯解析 helper（保留，仅命名遗留）

以下函数只做"外部语法 token -> 值"的转换，不含任何语义决策，因此保留；
所有 `#include` 与调用点都已在本次收缩中复核过。

| helper | 位置 | 用途 |
|--------|------|------|
| `foamIntValue` `foamDoubleValue` `foamBoolValue` `foamVectorValue` `foamScalarValue` | `SF_nativeDecode.cpp` | token -> 数值 |
| `foamUnquote` `foamWordList` `foamLower` `stripTokenQuotes` `upperToken` `normalize` | `SF_nativeDecode.cpp` | token 规范化 |
| `foamBCTypeValue` `foamThermalBCTypeValue` `foamWallHeatSourceType` | `SF_nativeDecode.cpp` | BC 名字 -> 强类型枚举 |
| `foamIntegerTriple` `foamNumbers` | `SF_nativeDecode.cpp` | token 列表 -> 整数三元组/数列表 |
| `sourceTokenKey` `sourceSchemeMentions` `appendSourceSchemeToken` | `SF_nativeDecode.cpp` | source term 名字规范化 |
| `readFoamWord` | `SF_nativeDecode.cpp` | 词读取 |
| `appendUnique` | `SF_nativeDecode.cpp` | 容器去重 |
| `isPlainRelativeFileName` | `SF_nativeDecode.cpp` | 路径合法性（唯一调用点：`addIBMGeometryFile` 的 fail fast） |

这些名字里的 `foam*` 只描述历史来源，不描述语义；重命名会带来无收益的全仓 diff，
因此本次保留。它们不是"中间文档翻译"，不违反 §5/§6/§7。

### 2.5 类 4 之外的 IO helper（合法保留）

| 结构 | 位置 | 性质 |
|------|------|------|
| `SF_documents.h/.cpp` | `compatibility/` | 23 LOC 纯磁盘代理（存在性、目录、文本读取），不解析语义 |
| `SF_sourceParserUtils.h/.cpp` | `compatibility/private/` | 20/54 LOC source token 工具 |
| `SF_foamParser.h/.cpp` | `infrastructure/io/private/` | 381 LOC mesh generator（`blockMeshDict`）词法读取器——`blockMeshDict` 是**网格格式**，不是 case 格式 |
| `SF_casePath.h/.cpp` | `infrastructure/io/private/` | 路径规范化 |

---

## 3. 前后对照

| 维度 | 收缩前 | 收缩后 |
|------|--------|--------|
| 语义入口 | 按 OpenFOAM 文件名/目录名推断 | `NativeCaseSections` typed slot + 语义 field 名 |
| 中间层 | 语义对象 -> OpenFOAM 形状文档 -> 再解析 | 语义对象 -> typed spec（一步） |
| 外部 OpenFOAM case | 由 `isOpenFOAMCase()` 分支支持 | 正式终止（§4）；`ModelLoader::read` 明确抛错 |
| `SF_nativeDecode.cpp` | 3412 LOC | 3232 LOC |
| 死兼容结构 | 8 项 | 0 项 |
| 回归保护 | 无 Gate A 守卫 | `resolvedSystemConsumerScope`（`--mode scope`）持续守卫 |

LOC 减少不是验收标准（§27）。验收标准是：类 2/类 3 结构按 §3 顺序被删除，且
删除后仍有测试覆盖替代路径：
[test_nativeIoNoFoamRoundTrip.cpp](../test/test_nativeIoNoFoamRoundTrip.cpp)
证明原生 case 目录不存在 `system/`、`constant/`、`0/`，typed section 直接被填充，
完整原生路径直接产出强类型 `CaseConfig`。

---

## 4. 仍未收缩的部分（诚实记录）

以下内容**没有**在本阶段删除，原因是它们仍承载真实语义，删除会改变数值行为或
触碰 §2 的禁区。

1. `Legacy::makeSolverConfig()` 骨架仍是 `CaseConfig::solver` 的构造起点（类 2）。
2. 约 82 个 legacy `SF::` 全局变量 + `ParserState` 事务（类 2，294 LOC）。
3. `SF_nativeDecode.cpp` 仍然是最大的单文件（3232 LOC）：它是 typed spec -> 强类型
   数值配置的唯一执行体，尚未按 equation/model 边界拆分。
4. `EulerianEquationSystem` 尚未迁移到 `EquationSystemInstance`（§25/§26 明确不在本
   阶段范围）。
5. `type: thermoDynamics` 语义对象路径目前没有任何 case 使用，未经端到端验证。
6. 少量 OpenFOAM 措辞残留在非 IO 层诊断字符串中，例如
   `parsing/SF_pressureConfigParser.h:63`、`execution/SF_singleFluid.cpp:99`、
   `execution/SF_multiPatch.cpp:155`、`execution/SF_eulerian.cpp:56`。它们是
   用户可见文案而非结构，不构成本阶段的 Gate A 阻塞项。
7. `tools/migrateRegistryCases.py` + `tools/compactRegistryCases.py` 是一次性迁移脚本，
   不属于 runtime IO 路径；按 §36 需要单独的调用点证明后才能删除。

---

## 5. 守卫

`python3 tools/check_program_decomposition.py --mode scope` 断言：

* 8 个已收窄的 consumer 不再包含 `SF_resolvedSimulationSystem.h`；
* `SF_compatibility.h` 仍声明 `NativeCaseSections`、`sections_`、按语义名的 `field()`、
  `decodeNativeCase`，且 14 个 typed slot 全部存在；
* 原生 IO 路径中不再出现 `parseOpenFOAM` / `documents_` / `isFoamControlDictPath` /
  `documentFile(` / `dictionaryPath(`；
* `SF_environment.cpp` 不再用 `solverConfig.numerics.solver` 推断 runtime service。

该守卫注册为 ctest 用例 `resolvedSystemConsumerScope`。
