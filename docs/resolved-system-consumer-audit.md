# Resolved System Consumer Audit (Gate B)

对应阶段的 §13 Gate B / §14-§19。审计 `ResolvedSimulationSystem` 的每一个 consumer：
它读了哪些字段、这些字段属于哪一类语义（WHAT / HOW / ORDER / RUNTIME / INSPECTION），
以及收缩之后它应该只依赖哪个编译产物。

`ResolvedSimulationSystem` 定义在
[SF_resolvedSimulationSystem.h](../src/solver/system/SF_resolvedSimulationSystem.h)。

---

## 1. 分类口径

| 类 | 含义 | 对应编译产物 |
|----|------|--------------|
| **WHAT** | 正在解哪些方程/约束 | `RawEquationSystem` → `ExecutableEquationSystem` |
| **HOW** | 每个算子怎么离散 | `CompiledNumericalSystem` |
| **ORDER** | 时间怎么推进、stage 怎么组织 | `TimeRecipe` → `CompiledSolvePlan` |
| **RUNTIME** | 运行期能力/后端/workspace 需求 | `RuntimeRequirements` |
| **INSPECTION** | 只被 explain/printer/CLI/测试读取，不参与数值传播 | `CaseClassification` + `TransformationRecord` |

收缩原则（§14/§15/§16/§17）：一个 consumer 如果只需要一类语义，它就该只接收那一类
编译产物，而不是聚合对象。

---

## 2. 聚合字段与 owner

收缩前 17 个 loose 字段 → 现在 10 个：

| # | 字段 | 类 | owner |
|---|------|----|-------|
| 1 | `formulation` | WHAT | 主未知量表述；**是 builder 的真实输入**（`SF_systemBuilder.cpp:656` 传给 `NumericalCompiler::compile`），因此不属 explain-only |
| 2 | `timeRecipe` | ORDER | 时间方法与 stage 数的唯一 authority（§30） |
| 3 | `classification` | INSPECTION | `CaseClassification`（§21，explain-only） |
| 4 | `rawSystem` | WHAT | 编译输入，只由 `CompilationResult` 暴露 |
| 5 | `executableSystem` | WHAT | `Program::equations` |
| 6 | `numericalSystem` | HOW | `Program::numerics` |
| 7 | `transformations` | INSPECTION | `TransformationRecord` 历史，只由 `CompilationResult` 暴露 |
| 8 | `executionPolicies` | ORDER/RUNTIME | solve plan 生成输入 |
| 9 | `solvePlan` | ORDER | `Program::solve` |
| 10 | `runtime` | RUNTIME | `RuntimeRequirements` |

被合并掉的 7 个字段（17 → 10）：

```text
runtime + executionCapabilities + requirements
  + workspaceRequirements + providerRequirements + runtimeServiceRequirements
        -> RuntimeRequirements runtime            (6 -> 1)

templateOrigin + densityBehavior + thermodynamicCompressibility
        -> CaseClassification classification      (3 -> 1)
```

---

## 3. Consumer 清单

### 3.1 已收窄（不再包含聚合头）

| consumer | 收缩前读到的 | 现在依赖 | 类 |
|----------|--------------|----------|-----|
| `src/solver/run/SF_planExecutor.h/.cpp` | 整个聚合 | `CompiledSolvePlan` | ORDER |
| `src/solver/system/SF_stateRealizer.h/.cpp` | 聚合（只为 `executableSystem`） | `ExecutableEquationSystem` + `RuntimeRequirements` | WHAT + RUNTIME |
| `src/solver/algorithm/SF_compressible.h/.cpp` | 聚合 | `ExecutableEquationSystem`, `CompiledNumericalSystem`, `CompiledSolvePlan`, `RuntimeRequirements` | WHAT + HOW + ORDER + RUNTIME |
| `src/solver/algorithm/pressureBased/eulerian/SF_pressureStepper.h/.cpp` | 聚合 | `ExecutableEquationSystem`, `CompiledSolvePlan`, `RuntimeRequirements` | WHAT + ORDER + RUNTIME |
| `src/app/application/SF_environment.h/.cpp` | 聚合（含 `solverConfig.numerics.solver` 推断 runtime service） | `RuntimeRequirements` + mesh/IBM runtime 输入 | RUNTIME |

`SF_compressible.h:62` 的注释显式记录了收缩结论：
"不接收 `ResolvedSimulationSystem`，也不构造 `AlgorithmContext`。"

### 3.2 合法持有聚合（composition 层）

这些文件位于"装配层"，它们的职责就是读取多类语义并把它们分派给上表的 consumer。
按 §14 它们可以持有聚合，但**不得**把聚合继续向下传递。

| consumer | 读到的字段 | 类 |
|----------|-----------|-----|
| `src/app/application/execution/SF_execution.h/.cpp` | `executableSystem`, `numericalSystem`, `solvePlan`, `runtime` | 全类（分派） |
| `src/app/application/execution/SF_flowLoop.h/.cpp` | 同上 | 全类（分派） |
| `src/app/application/execution/SF_singleFluid.cpp` | 同上 | 全类（分派） |
| `src/app/application/execution/SF_eulerian.cpp` | 同上 | 全类（分派） |
| `src/app/application/execution/SF_multiPatch.cpp` | 同上 | 全类（分派） |
| `src/app/application/SF_application.cpp` | `runtime`, `timeRecipe`, `classification` | 全类（顶层编排） |
| `src/app/application/SF_inspection.h/.cpp` | 全部（生成 `BuildRequest`） | 全类（检查） |

### 3.3 生产者与 explain 层

| consumer | 职责 | 类 |
|----------|------|-----|
| `src/solver/system/SF_systemBuilder.h/.cpp` | 生产聚合 | 全类 |
| `src/solver/system/SF_systemValidator.h/.cpp` | 校验 | WHAT + ORDER + RUNTIME |
| `src/solver/system/SF_systemPrinter.h/.cpp` | explain | INSPECTION |
| `src/models/ibm/SF_ibmExplain.h` | explain（仅注释引用） | INSPECTION |
| `test/test_pisoArchitecture.cpp` | 断言 `classification` | INSPECTION |
| `test/test_programDecomposition.cpp` | 断言 `Program` 视图语义 | INSPECTION |

### 3.4 间接消费者（只经编译产物）

| consumer | 依赖 | 说明 |
|----------|------|------|
| `src/solver/system/SF_equationContribution.h` | aggregate header 只为类型可见性 | 不读聚合数据成员 |
| `src/solver/system/SF_numericalCompiler.h` | 同上 | 产出 `CompiledNumericalSystem` |
| `src/solver/system/SF_solvePlan.h` | 同上 | 产出 `CompiledSolvePlan` |

---

## 4. 前后对照

| 维度 | 收缩前 | 收缩后 |
|------|--------|--------|
| 聚合 loose 字段 | 17 | 10 |
| 只需一类语义的 consumer | 5 个直接接收聚合 | 5 个只接收各自编译产物 |
| `AlgorithmContext` / `ProgramContext` | 风险：为压缩参数而引入 | 未引入（§16 明确禁止） |
| 聚合头的 include 者 | 未记录（historically 无法从 git 还原：相关文件在本分支上未受版本控制） | **10** 个文件（`SF_inspection.h`、`SF_execution.h`、`SF_flowLoop.h`、`SF_equationContribution.h`、`SF_numericalCompiler.h`、`SF_solvePlan.h`、`SF_systemBuilder.h`、`SF_systemPrinter.h`、`SF_systemValidator.h`、`test_programDecomposition.cpp`）。另有 5 个 consumer（10 个文件）已退出，见下 |
| `Program` 携带状态 | 未定义 | 4 个 const 引用的纯视图，`sizeof == 4 * sizeof(const void*)` |
| `compileSystemProducts` 的返回 | 三个并列 loose 字段 | `CompilationResult{ raw, transformations, program }` |

---

## 5. Program 视图的边界（§18）

`Program` 只包含四类编译产物：

```cpp
struct Program {
    const ExecutableEquationSystem& equations;
    const CompiledNumericalSystem& numerics;
    const CompiledSolvePlan& solve;
    const RuntimeRequirements& runtime;
};
```

**明确不含**（§18 禁止项，逐条核对）：

| 禁止项 | 是否包含 | 依据 |
|--------|----------|------|
| `CaseConfig` | 否 | `SF_runtimeRequirements.h` 不 include `SF_caseConfig.h`；测试刻意不 include 它 |
| mesh | 否 | 成员只有 4 个 |
| `Field` | 否 | 同上 |
| MPI communicator | 否 | 同上 |
| IBM runtime state | 否 | 同上 |
| output settings | 否 | 同上 |
| raw system | 否 | raw system 在 `CompilationResult`，不在 `Program` |
| transformation history | 否 | 同上 |
| printer strings | 否 | 同上 |

`CompilationResult`（§19）承担"编译过程记录"：

```cpp
struct CompilationResult {
    const RawEquationSystem& raw;
    const std::vector<TransformationRecord>& transformations;
    Program program;
};
```

### 5.1 与 §18 写法的偏离（必须记录）

§18 的示意使用值成员；实现改用 **const 引用**。理由：

`CompressibleAlgorithm::bindSolvePlan` 会做身份检查 `&plan != &solve_`，若 plan 不是
它自己的就抛错。复制 `CompiledSolvePlan` 会破坏这个不变量，从而改变行为。因此
`Program`/`CompilationResult` 是**视图**而不是容器，§18 关于"不得携带运行期状态"的
约束由"只有引用、`sizeof` 等于四个指针"来保证，而不是由值语义来保证。

该偏离已由
[test_programDecomposition.cpp](../test/test_programDecomposition.cpp) 固化：
成员数量、成员顺序/类型、`sizeof`、以及每个引用都指向聚合自身（不是副本）。

---

## 6. 守卫

| 守卫 | 断言 |
|------|------|
| `resolvedSystemConsumerScope`（`tools/check_program_decomposition.py --mode scope`） | §3.1 的 10 个文件不再 include `SF_resolvedSimulationSystem.h`；`SF_environment.cpp` 不再用 `solverConfig.numerics.solver` |
| `resolvedSystemFieldCount`（`--mode fieldcount`） | 聚合 loose 字段数 ≤ 10，只允许下降 |
| `programDoesNotOwnRuntimeState`（`test/test_programDecomposition.cpp`） | `Program` 恰好 4 个成员、纯引用视图、无 owned state |
| `check_architecture.py` | `RuntimeRequirements` 是 RUNTIME 的唯一 authority；聚合只有 `RuntimeRequirements runtime`，不得重复声明 loose 需求向量 |

---

## 7. 仍未收缩的部分（诚实记录）

1. **composition 层仍持有聚合**（§3.2 的文件）。这是当前架构的正常状态：它们是
   装配点。真正需要继续收缩的是 §3.3 的 `SF_systemValidator` / `SF_systemPrinter`，
   它们目前仍接收整个聚合。
2. `SF_systemBuilder.cpp` 约 995 LOC，同时负责 `BuildRequest` 判定、capability 需求
   生成、编译调用与聚合装配。它的拆分属于下一阶段，本阶段只删除了冗余逻辑（§24）。
3. `transformation history` 目前是 `std::vector<TransformationRecord>`，只被 explain
   读取；它是否应该继续存在于聚合中，取决于后续 explain 的统一入口设计。
4. §3.4 记录的 3 个**编译层内部头文件**仍 include 聚合头，只为类型可见性
   （`SF_equationContribution.h`、`SF_numericalCompiler.h`、`SF_solvePlan.h`）。
   消除它们需要重新安排 `RawEquationSystem` / `CompiledNumericalSystem` 的声明位置，
   属于下一阶段。此外 `SF_systemBuilder.h` / `SF_systemPrinter.h` / `SF_systemValidator.h`
   也 include 聚合头，但它们是聚合的真实生产者 / consumer，属于 §3.2 / §3.3 的正常状态。
