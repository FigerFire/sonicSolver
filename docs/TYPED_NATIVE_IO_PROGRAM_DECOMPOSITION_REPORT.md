# Phase Report — Typed Native IO + Program Decomposition

**Branch**: `refactor/explicit-time-plan-authority-20260918`
**Base commit**: `acdd66a refactor: move explicit timestep control to compiled plan`
**Constraint**: 按提示词修改代码，不修改任何数学算法。

本阶段有两个子阶段：

* **A — Typed Native IO**：删除原生 IO 路径中的 OpenFOAM 形状中间文档翻译，
  让原生语义 YAML 直接产出强类型 spec。
* **B — Program Decomposition**：拆解 `ResolvedSimulationSystem` 的 consumer，
  使每个 consumer 只依赖它真正需要的编译产物。

---

## A. 交付物

### A.1 代码

| 文件 | 变更 | 性质 |
|------|------|------|
| [SF_compatibility.h](../src/app/application/model/compatibility/SF_compatibility.h) | 完全重写（161 LOC） | `documents_` 伪路径映射 → `NativeCaseSections` typed slot + `fields`；`parseOpenFOAM*`/`parse*File` → `decode*`；新增只读 `sections()` |
| [SF_nativeBinding.cpp](../src/app/application/model/SF_nativeBinding.cpp) | 重写（220 LOC） | 直接填充 typed section；`meshParamFile_` → `meshGenerator_`；`thermoDynamics` 工厂显式 emit |
| [SF_nativeDecode.cpp](../src/app/application/model/compatibility/SF_nativeDecode.cpp) | 机械重写 3412 → 3232 LOC | 19 × `Model::Field` → `Model::FieldDescriptor`；`constant/triSurface/` 死分支 → §22 fail fast；约 16 条用户可见文案去 OpenFOAM 措辞 |
| [SF_documents.cpp](../src/app/application/model/compatibility/SF_documents.cpp) | 重写为 23 LOC | 纯磁盘代理 |
| [SF_compatibility.cpp](../src/app/application/model/compatibility/SF_compatibility.cpp) | 127 LOC | 删除 legacy ctor 分支 |
| `infrastructure/io/private/SF_casePath.h/.cpp` | 删除死函数 | `isFoamControlDictPath` |
| `models/initial/SF_parameter.h/.cpp`、`compatibility/SF_parserState.h` | 删除 8 个死全局 | 全仓 grep 证明 0 reader |
| 4 × `test/<case>/solvers/algorithm.yaml` | corrector key 重命名 | 值与原 C++ 默认值相同，行为不变 |

### A.2 文档

| 文档 | 对应要求 |
|------|----------|
| [docs/native-io-contraction-audit.md](../docs/native-io-contraction-audit.md) | §3 Gate A：真实调用链 + 4 类结构分类 + 删除依据 |
| [docs/case-config-ownership.md](../docs/case-config-ownership.md) | §11 + §21 + §22/§23 |
| [docs/resolved-system-consumer-audit.md](../docs/resolved-system-consumer-audit.md) | §13 Gate B |
| 本文 | §38 A-I |

### A.3 sub-phase B 代码

| 文件 | 变更 |
|------|------|
| [SF_runtimeRequirements.h](../src/solver/system/SF_runtimeRequirements.h)（103 LOC） | 新 owner：`RuntimeRequirements`（6 成员）、`Program`（4 const 引用）、3 个 `requires*` 声明 |
| [SF_caseClassification.h](../src/solver/system/SF_caseClassification.h)（27 LOC） | **新增**：3 个 explain-only 事实 |
| [SF_resolvedSimulationSystem.h](../src/solver/system/SF_resolvedSimulationSystem.h)（76 LOC） | loose 字段 17 → 10；新增 `CompilationResult` + `programOf()` / `compilationOf()` |
| [SF_buildRequest.h](../src/solver/system/SF_buildRequest.h)（41 LOC） | 删除 2 个后端能力 flag |
| [SF_systemBuilder.cpp](../src/solver/system/SF_systemBuilder.cpp) | `requires*` 定义 + 委托聚合重载；`hasEquationPrefix` 收窄；classification 写入；2 个 capability 常量 |
| [SF_compressible.h/.cpp](../src/solver/algorithm/SF_compressible.h) | §16 收窄 |
| [SF_pressureStepper.h/.cpp](../src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.h) | §17 收窄（PISO 数学未改） |
| [SF_stateRealizer.h/.cpp](../src/solver/system/SF_stateRealizer.h) | 收窄 |
| [SF_environment.h/.cpp](../src/app/application/SF_environment.h) | §15 收窄 |
| [SF_resolvedEquationSystem.h](../src/solver/system/SF_resolvedEquationSystem.h) | 新增 `hasEquationPrefix(ExecutableEquationSystem&, string_view)` |
| [SF_solvePlan.cpp](../src/solver/system/SF_solvePlan.cpp) | 删除重复的局部 `hasEquationPrefix` |

### A.4 测试与守卫

| 名称 | 类型 | 对应要求 |
|------|------|----------|
| `nativeIoNoFoamRoundTrip` | C++（[test_nativeIoNoFoamRoundTrip.cpp](../test/test_nativeIoNoFoamRoundTrip.cpp)） | §10 / §35 |
| `resolvedSystemConsumerScope` | Python（[check_program_decomposition.py](../tools/check_program_decomposition.py) `--mode scope`） | §15/§16/§17 / §35 |
| `resolvedSystemFieldCount` | Python（同工具 `--mode fieldcount`） | §20 / §35 |
| `programDoesNotOwnRuntimeState` | C++（[test_programDecomposition.cpp](../test/test_programDecomposition.cpp)） | §18 / §35 |
| `check_architecture.py` | 已修改 | §34：RUNTIME authority token 迁到 `SF_runtimeRequirements.h`；新增"聚合只有 `RuntimeRequirements runtime`"守卫 |

---

## B. 删除的东西（§27 的验收口径）

**本阶段的验收标准是删除了什么，不是文件变短。**

1. **删除了中间翻译层**：`documents_`（伪路径 → 字典文本）被 `NativeCaseSections` +
   `fields` 取代。解码不再有"先合成文档再解释"的阶段。
2. **删除了重复 authority**：
   * 6 个 runtime 需求向量 → 1 个 `RuntimeRequirements`；
   * 3 个 explain-only 事实 → 1 个 `CaseClassification`；
   * `BuildRequest` 中 2 个可由后端能力常量替代的 flag → 文件作用域 `constexpr`。
3. **删除了 8 个死兼容结构**（清单见
   [native-io-contraction-audit.md §2.3](../docs/native-io-contraction-audit.md)）。
4. **删除了死代码**：`isFoamControlDictPath`（声明 + 定义）、`isOpenFOAMCase()`、
   6 个只写不读的 legacy 全局、dynamic-mesh/moving-boundary 警告路径、
   `SF_solvePlan.cpp` 中的重复 `hasEquationPrefix`。
5. **删除了 consumer 的过度依赖**：5 个只需单一编译产物的 consumer 不再 include 聚合头。

---

## C. `Program` / `CompilationResult` 设计（§18/§19）

```cpp
struct Program {
    const ExecutableEquationSystem& equations;   // WHAT
    const CompiledNumericalSystem&  numerics;    // HOW
    const CompiledSolvePlan&        solve;       // ORDER
    const RuntimeRequirements&      runtime;     // RUNTIME
};

struct CompilationResult {
    const RawEquationSystem& raw;                // 编译输入
    const std::vector<TransformationRecord>& transformations;  // INSPECTION
    Program program;
};
```

**偏离记录**：§18 的示意使用值成员；实现改用 const 引用。

理由是**数值不变量**，不是风格：`CompressibleAlgorithm::bindSolvePlan` 会做身份检查
`&plan != &solve_`，若传入的 plan 不是它自己的就抛错。复制 `CompiledSolvePlan`
会破坏这个不变量。因此 `Program` / `CompilationResult` 是**视图**，
`sizeof(Program) == 4 * sizeof(const void*)` 由 `programDoesNotOwnRuntimeState` 固化。

§18 关于"不得携带 CaseConfig / mesh / Field / MPI communicator / IBM runtime state /
output settings / raw system / transformation history / printer strings"的约束逐条
核对结果见
[resolved-system-consumer-audit.md §5](../docs/resolved-system-consumer-audit.md)。

---

## D. §21 explain-only 字段

`templateOrigin`、`densityBehavior`、`thermodynamicCompressibility` 经审计确认只被
printer / CLI / 测试读取（`SF_systemPrinter.cpp:487-489`、`SF_cli.cpp:306`、
`test_pisoArchitecture.cpp:460-461,477-478`），从不参与数值传播，因此被归入
`CaseClassification`。

**偏离记录**：§21 的字面表述是"移出运行期传播路径"。实现选择**分组**（三个字段进入一个
明确标注 `EXPLAIN-ONLY` 的独立类型），而不是把它们完全从聚合中移走。理由是它们确实是
`ResolvedSimulationSystem` 的解码产物，且 explain 需要它们；分组已经达成"不再作为
RUNTIME 传播"的语义隔离目标。

`ResolvedSimulationSystem::formulation` **未被移动**：它是 builder 的真实输入
（`SF_systemBuilder.cpp:656` 传给 `NumericalCompiler::compile`），不属于 explain-only。

---

## E. §22/§23 BuildRequest

| 字段 | 依据 | 处置 |
|------|------|------|
| `constraintGlobalDofAvailable` | 唯一生产者 `SF_inspection.cpp:133` 硬编码 `true`；语义是后端能力声明，不是 case 属性 | 删除 → `constexpr bool kConstraintGlobalDofAvailable = true;` |
| `distributedLinearSystemAvailable` | 同上（`:134`） | 删除 → `constexpr bool kDistributedLinearSystemAvailable = true;` |

行为保持性可证：原表达式求值为 `!constraintDof || true` 与 `true`，替换为常量后结果
逐位相同。其余 12 个字段都是 builder 的真实判定输入，尚无 typed owner，因此保留。

---

## F. 冻结数值 baseline 验证（§23/§33）

| 项 | 结果 |
|----|------|
| Sod WENO7 + Rusanov, forwardEuler, t = 0.2 | **bit-identical**，`check_sod_rusanov_regression.py` 通过 |
| 步数 | 454（与 baseline 相同） |
| VTS SHA-256 | `8bf55a1c4973f3a2322aeb8785305bf0589cc9eb4cf89505203309f0ccba0768`（未变） |
| PISO (`pressureConstraintPiso`) | **byte-identical to frozen legacy result** |
| PISO VTS SHA-256 | `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511` |

**数值算法零改动**：本阶段没有触碰 WENO/TENO、Rusanov/Steger-Warming/Roe/HLLC、
central diffusion、RK 系数、stage time、CFL、EOS 数学、positivity、boundary/halo 顺序、
canonical face COPY、GlobalDof SUM、MPI ownership。

唯一的数值相邻改动是 `SF_pressureStepper` 的**参数形状**（从接收聚合改为接收
3 个显式引用）与该类 `equationIds()` helper 的入参类型；PISO 的求解数学、迭代顺序、
松弛因子均未改变，已由 byte-identical PISO 回归证实。

---

## G. §32 回归门

每个子阶段结束后都跑完整门；最终一次全部 GREEN：

| 命令 | 结果 |
|------|------|
| `cmake --build build --parallel 6` | exit 0，无 error/warning |
| `cmake --build build-tests --parallel 6` | exit 0 |
| `cmake --build build-ibm --target sonicSolver --parallel 10`（Release `-O3`，刻意未定义 `NDEBUG`，assert 保留） | exit 0 |
| `python3 tools/check_architecture.py` | exit 0；allowlist 仍为 **10** 条，未增长 |
| `python3 tools/check_program_decomposition.py --mode scope` / `--mode fieldcount` | exit 0 |
| `ctest --test-dir build-tests --output-on-failure` | **12/12 passed**（原 8 项 + 新增 4 项） |
| `check_sod_rusanov_regression.py` | passed（bit-identical）；Debug `build/` 与新 `build-ibm/`（`-O3`）两套二进制均通过 |
| `check_piso_regression.py` | passed（byte-identical）；两套二进制得到同一个 VTS SHA-256 |
| `sonicSolver check` 冒烟（`find test -name case.yaml`） | **30/30 case OK** |
| WENO7 + StegerWarming 单步（`weno7Steger` 覆盖，临时 case） | rc=0，`reconstruction=WENO7, numericalFlux=StegerWarming` |
| `mpirun -np 4` Sod（`test/Sod/sodCase`）`--steps 1` / `--steps 3` | rc=0，`Solver: densityBase`、`Mesh ready: 4 partition(s)` |
| 执行冒烟：Sod 系列 + `createMesh`/`thermalCase`/`pressureConstraintPiso`/`eulerianEulerianCase`/`rpiWallBoilingCase`/`surfaceKKTIBMCase`/`variationalIBMCase`/`guiCase/sodGuiFromSfm`，`--steps 3` | 全部 rc=0 |
| 10 个串行 IBM `cylinderFlow*`，`--steps 3` | **全部 rc=0**（Debug 与 `-O3` 两套二进制；两者 `time`/`dt` 逐位相同） |
| `test/IBM/IBMCase` 串行 | 预期 §22 fail-fast：rc=255，`runtime MPI ranks=1 must equal [parallel].split product=4` |
| `test/mixtureCase` | 预期 §22 fail-fast：`weno3Steger` 需要 five-variable PerfectGas（见 [executable-system-audit-2026-09-16.md §93](../docs/executable-system-audit-2026-09-16.md)） |

新增 4 个用例的注册位置：根 [CMakeLists.txt](../CMakeLists.txt) 的测试段。

> 注：`build-ibm/` 是本阶段为缩短 IBM 冒烟耗时临时创建的 Release 构建目录
> （`-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O3" -DBUILD_TESTS=OFF`，
> 刻意不定义 `NDEBUG` 以保留 assert）。它只用于**交叉验证 `-O3` 与 Debug 数值等价**
> （Sod bit-identical、PISO 同 SHA、10 个 IBM case 的 `time`/`dt` 逐位相同），
> 验证完成后已按 `.gitignore` 覆盖范围删除；交付物只依赖 `build/` 与 `build-tests/`。
> 上表中所有 Debug 构建 / `ctest` / `check_architecture.py` 的结果都在清理后重跑确认。

### G.1 有界冒烟口径（诚实说明）

单个 `cylinderFlow*` case 跑到 `endTime = 5.0` 需要约 11400 步。实测单核吞吐约
2.9 step/s（`--steps` 计数 / 计时），即**单 case 约 1 小时 CPU**；本机 10 核，
10 个 case 串行约 10 小时。因此最终门采用**有界执行冒烟**：

* 10 个 `cylinderFlow*` 统一 `--steps 3`：验证 typed native IO → executable system →
  compile → Plan → execute 的**全部执行路径**（RK stage、IBM 几何/forcing、
  state closure、IBM 各方法分支）都已跑通且无 fatal。
* 额外证据：早期批次中 `cylinderFlowDFMAugmentedLagrangian` 与
  `cylinderFlowFullyImplicitDLM` **完整跑到 `time = 5.000000e+00`**（各自 91907 行日志，
  末行 `[SF] Step time: time=5.000000e+00, dt=3.896002e-04`，随后写出
  `*_t5.vts` 与 `.pvd`）。这证明长时程推进本身没有回归。
* 本文件不再声称"10 个 case 全部跑到 endTime"；跑满 endTime 的长时程覆盖只有上述 2 个 case。
* 这个口径与本 repository 既有阶段的 IBM 冒烟口径一致：见
  [refactor-phase1a-report.md:155](../docs/refactor-phase1a-report.md)
  （`./build/sonicSolver --steps 1 CASE` 覆盖 every non-MPI `cylinderFlow*`）与
  [compiled-plan-completeness-migration.md:381](../docs/compiled-plan-completeness-migration.md)
  （one-step runtime diagnostic）。本阶段取 `--steps 3`，略严于既有口径。
* 与 [refactor-phase1a-report.md](../docs/refactor-phase1a-report.md) 记录的
  "3 个 KKT case 停在 `ConstraintGlobalDof -> owned HYPRE row -> lambda COPY` 未启用"
  相比，本阶段这 3 个 case（`DFMAugmentedLagrangian`、`DFMImplicitSelfPropelled`、
  `FullyImplicitDLM`）均 rc=0，不再停在能力缺口上。

---

## H. 预期取消但未做的清单

* 未引入 `AlgorithmContext` / `ProgramContext` / `CaseManager` / `ConfigManager` /
  `InputContext` / `CaseService` / `ConfigurationFacade`（§8/§16 禁止）。
* 未恢复 `TimeScheme` / `explicitStageCount` / `ddtDispatch`，未创建
  `ExplicitSolver` / `ImplicitSolver` / `RK4Solver`（§30 禁止）。
* 未修改 `SF_communicationPlan.cpp` / `SF_ilwClosure.cpp` / `SF_haloExchange.cpp` /
  `SF_MultiBlockMesh.cpp`（§31）。
* 未重组 `core/`（§30）。
* 未迁移 `EulerianEquationSystem`，未改动 `create("eulerian","liquid"|"gas"|"particles")`（§25/§26）。
* 未把 God object 改名为 `ProgramContext`。

---

## I. 剩余 legacy 与已推迟问题（诚实记录）

### I.1 native IO 侧

| 剩余项 | 规模 | 性质 | 移除条件 |
|--------|------|------|----------|
| `Legacy::makeSolverConfig()` 骨架 | `SF_legacyConfig.h/.cpp` | §3 类 2 | `FDM::SolverConfig` 每个成员都能由 typed section 直接构造 |
| legacy `SF::` 全局 + `ParserState` 事务 | 约 82 个全局，294 LOC，静态 mutex 保护 | §3 类 2 | `decodeNativeCase()` 每个写入点改为 typed 局部量 |
| `SF_nativeDecode.cpp` 单文件 | 3232 LOC | 尚未按 equation/model 边界拆分 | 下一阶段 |
| `foam*` 值转换 helper | 约 21 个函数 | §3 类 4 | 保留：只描述历史来源，重命名无收益 |
| `SF_foamParser.h/.cpp` | 381 LOC | 合法：`blockMeshDict` 是**网格格式**而非 case 格式 | 不需要 |
| `tools/migrateRegistryCases.py` + `compactRegistryCases.py` | 一次性迁移脚本 | 非 runtime 路径 | 需要独立的调用点证明（§36） |
| 少量 OpenFOAM 措辞残留 | `parsing/SF_pressureConfigParser.h:63`、`execution/SF_singleFluid.cpp:99`、`execution/SF_multiPatch.cpp:155`、`execution/SF_eulerian.cpp:56` | 用户可见文案，非结构 | 文案清理，不阻塞 Gate A |
| `type: thermoDynamics` 语义对象路径 | 无 case 使用 | **未端到端验证** | 需要新增测试 case |

### I.2 架构侧

| 剩余项 | 说明 |
|--------|------|
| composition 层仍持有聚合 | `SF_execution` / `SF_flowLoop` / `SF_singleFluid` / `SF_eulerian` / `SF_multiPatch` / `SF_application` / `SF_inspection`。这是装配点的正常状态 |
| `SF_systemValidator` / `SF_systemPrinter` 仍接收整个聚合 | 应收窄为各自需要的编译产物（下一阶段） |
| `SF_systemBuilder.cpp` ~995 LOC | 同时负责 `BuildRequest` 判定、capability 需求、编译调用、聚合装配。本阶段只删除冗余逻辑（§24） |
| `transformation history` 仍在聚合上 | 只被 explain 读取；去留取决于后续 explain 统一入口 |
| 3 个编译层内部头文件仍 include 聚合头（类型可见性） | `SF_equationContribution.h`、`SF_numericalCompiler.h`、`SF_solvePlan.h`；消除需要重新安排声明位置。另外 `SF_systemBuilder.h` / `SF_systemPrinter.h` / `SF_systemValidator.h` 也 include 聚合头，但它们是生产者 / 真实 consumer，属正常状态 |

### I.3 与本阶段无关的既存失败（有据可查，不在任何回归门内）

| 问题 | 证据 |
|------|------|
| `test/IBM/IBMCase` 在 `mpirun -np 4` 下报 `Resolved unknown 'rho' expects storage 'conservative' ... expected=2, actual=10` | 本阶段开始前即存在，见 [refactor-run-composition-root-report.md §9-10](../docs/refactor-run-composition-root-report.md) |
| `cylinderFlowPeskinMPI2` 在 `ImmersedForcingSystem::reduceSurfaceVectors` 的 `MPI_Allreduce` 死锁 | 已验证是 Peskin 路径特有：`cylinderFlowFictitiousDomainMPI2`（body dofs=138，`max|Ju-Us|=0`）与 `cylinderFlowVelocityForcingMPI2` 均正常退出 |
| `test/mixerVessel2DLevelSet` 在 `mpirun -np 2` 下报 `Resolved unknown 'rho' expects storage 'conservative' on every participating patch; expected=3, actual=5.` | 与第 1 行的 `IBMCase` 故障**同一签名、同一层**（distributed multi-patch storage 解析）。该故障已由 [refactor-run-composition-root-report.md §10](../docs/refactor-run-composition-root-report.md) 用**重构前二进制复现**证明为既存问题；本阶段未触碰 storage 解析代码。串行下该 case 按 `split: 2` 正确 fail-fast（rc=255），属 §22 预期行为 |

以上都不属于本阶段交付范围，且被本阶段之前的文档记录为既存状态。

---

## J. §26 变更汇报（10 点）

1. **影响哪一层**
   IO 层（case 输入侧）+ solver/system 层的所有权与接口形状。数值层未改。

2. **属于哪一类**
   sub-phase A：既不是"加项"也不是"加方程"也不是"改怎么解方程"，而是**删除中间翻译层**
   （§5/§6/§7）。sub-phase B：**改怎么解方程**中的"调用形状"部分——不改变任何离散或
   时间推进数学，只收窄参数契约。

3. **ownership before/after**

   | 对象 | before | after |
   |------|--------|-------|
   | 原生 case 语义入口 | `documents_` 伪路径 map | `NativeCaseSections` + `fields` |
   | RUNTIME 需求 | 聚合上 6 个 loose 向量 | `RuntimeRequirements` |
   | explain-only 事实 | 聚合上 3 个 loose 字段 | `CaseClassification` |
   | 后端能力 flag | `BuildRequest` 2 个 bool | `SF_systemBuilder.cpp` 内 `constexpr` |
   | `CompiledSolvePlan` identity | 隐式 | `Program::solve` 必须引用聚合自身（有意固化） |

4. **execution order before/after**
   `CaseIO::read → ModelLoader::build → CaseAdapter::build → decodeNativeCase → CaseConfig`
   顺序未变；变的是中间不再出现"合成 OpenFOAM 文档 → 再解析"这一步。
   solver 执行顺序完全未变（Sod 454 步 bit-identical、PISO byte-identical）。

5. **public API before/after**

   | 变更 | 类型 |
   |------|------|
   | `CaseAdapter::parseOpenFOAM*` / `parse*File` → `decode*` | rename |
   | `CaseAdapter::sections()` 新增 | 新增只读观察入口 |
   | `realizeState(const ResolvedSimulationSystem&, …)` → `(const ExecutableEquationSystem&, const RuntimeRequirements&, …)` | 收窄 |
   | `CompressibleAlgorithm` ctor：聚合 → 4 个显式引用 | 收窄 |
   | `PressureStepper` ctor：聚合 → 3 个显式引用 | 收窄 |
   | `Environment::build(…, const ResolvedSimulationSystem&, …)` → `(…, const RuntimeRequirements&, …)` | 收窄 |
   | `programOf()` / `compilationOf()` 新增 | 新增 |
   | `BuildRequest` 删除 2 个字段 | 删除 |
   | `hasEquationPrefix(ExecutableEquationSystem&, string_view)` 新增 | 新增 |

6. **dependency before/after**
   5 个 consumer（10 个文件：`SF_planExecutor`、`SF_stateRealizer`、`SF_environment`、
   `SF_compressible`、`SF_pressureStepper`）的 `SF_resolvedSimulationSystem.h` include 被移除，
   改为各自的编译产物头。当前直接 include 聚合头的文件为 10 个。依赖方向未反转
   （§16 的目标依赖主干保持）。

6b. **诚实说明**：本分支上新文件均未纳入 git 跟踪（`git status` 显示为 untracked），
   因此无法给出重构前 include 集合的精确 diff。上面记录的是**当前实测值**，
   以及由 `tools/check_program_decomposition.py --mode scope` 固化的收窄清单。

7. **MPI/parallel semantics 是否变化**
   **未变**。`ResolvedSimulationSystem` 从不携带 MPI communicator；`RuntimeRequirements`
   只把既有的需求向量聚合在一起，SUM/COPY 语义、canonical face、owner 规则都未触碰。
   `mpirun -np 4` Sod 单步通过。

8. **numerical semantics 是否变化**
   **未变**。Sod WENO7/Rusanov t=0.2 bit-identical（454 步，VTS SHA-256 未变）；
   PISO byte-identical（VTS SHA-256 `4bc3e0b…72511` 未变）。

9. **regression tests**
   新增 4 个（`nativeIoNoFoamRoundTrip`、`resolvedSystemConsumerScope`、
   `resolvedSystemFieldCount`、`programDoesNotOwnRuntimeState`），并修改
   `check_architecture.py` 的 RUNTIME authority 断言。全套 12/12 通过。

10. **deferred issues**
    见 §I。最重要三项：`SF_nativeDecode.cpp` 尚未按 equation/model 边界拆分；
    `SF_systemValidator` / `SF_systemPrinter` 仍接收聚合；
    `type: thermoDynamics` 路径未经端到端验证。

---

## K. 下一阶段（§39）

**Eulerian EquationSystem Instance Migration。**

本阶段已为该迁移铺好前置条件：

* Equation 层已有 typed 入口（`NativeCaseSections` + `Model::Description`），
  不再需要"文档 → 语义"的翻译；
* `ExecutableEquationSystem` 已是可独立传递的编译产物，consumer 已按 WHAT 收窄；
* 聚合字段数与 RUNTIME authority 已有守卫，迁移过程中不会悄悄重新长出重复 authority。
