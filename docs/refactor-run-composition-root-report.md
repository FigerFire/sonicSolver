# Run Composition Root 拆分报告

> 将 `SF_run.cpp` 的 `runConfiguredCase()` 从 300+ 行巨型装配流程压缩为
> `parse → read → compile → validate → build → validate → execute`
> 六步 composition root，并把 CLI 校验、program 校验、environment 构建、
> environment 校验分别落到独立文件。

## 0. 变更文件

| 文件 | 变更 |
|---|---|
| [SF_run.cpp](../../src/app/application/run/SF_run.cpp) | 重写为 6 步编排 + 匿名 `execute()` 路由 |
| [SF_preflight.h](../../src/app/application/SF_preflight.h) / [SF_preflight.cpp](../../src/app/application/SF_preflight.cpp) | 新增：CLI 解析 + Phase-1 program 校验 |
| [SF_environment.h](../../src/app/application/SF_environment.h) / [SF_environment.cpp](../../src/app/application/SF_environment.cpp) | 新增：DomainKind + ExecutionEnvironment + 环境构建/校验 + DomainBuilder |
| [CMakeLists.txt](../../src/app/application/CMakeLists.txt) | 注册 `SF_preflight.cpp`、`SF_environment.cpp` 到 `SF_application` |

## 1. 影响哪一层

Application Composition Root（`SF::Application::runConfiguredCase`）及其两个
新增的 application 层模块：

- `SF_preflight.*` —— Phase-1 前置校验（CLI + program 能力）
- `SF_environment.*` —— Phase-2 执行环境构建与 runtime 校验

不触碰 Equation / Discretization / Coupling / Time / Linear Algebra /
Infrastructure backend 的任何数值实现。

## 2. 属于“加项 / 加方程 / 改求解方式”中的哪一类

**三者都不是。** 这是一次纯 composition-root 的职责拆分：

- 没有给任何方程增加 term；
- 没有新增 equation / constraint / unknown；
- 没有改变任何算子离散或时间推进策略。

它只是把“装配阶段”的职责边界切开：数学编译（compile）与执行拓扑
（MPI / single / multi-patch）不再混在同一个函数里。

## 3. Ownership before / after

| 对象 | Before | After |
|---|---|---|
| CLI 解析 | `runConfiguredCase` 内联 | `Preflight::parseCommandLine` |
| usage 输出 | `runConfiguredCase` 内联 | `Preflight::reportUsage` |
| program 能力 fail-fast | `runConfiguredCase` 内联 | `Preflight::validateProgram` |
| MPI / mesh / storage / IBM runtime 构建 | `runConfiguredCase` 内联 | `Environment::build` |
| runtime topology 校验 | `runConfiguredCase` 内联 | `Environment::validate` |
| createMesh 路径 | `runConfiguredCase` 内联 | `Environment::buildMeshOnly` |
| single/multi 布局判定 | `runConfiguredCase` 内联 `needsMultiFieldMPI` | `Environment::resolveSingleFieldLayout`（DomainBuilder） |
| 最终 runner 分发 | `runConfiguredCase` 内联 `if/else` | 匿名 `execute()`（只看 `env.domain`） |

`ExecutionEnvironment`（栈对象）是本次运行的运行时资源 owner：
`ParallelContext`（move-only，`unique_ptr`）、`Field` / `MultiBlockMesh`、
`CompositeIB` / `IB`。`ResultWriter` 仍由 `runConfiguredCase` 栈持有，
`Environment` 只持非拥有指针。

## 4. Execution order before / after

Before（原 `runConfiguredCase`）：

```
parse CLI → read case → configure output → linearMPI 推断 → init MPI
→ create mesh / inspectCase → Ghost 预处理 → single/multi 判定
→ runner 路由 → boundary empty-dims → IBM setup → distributed IBM 校验
→ execute
```

After：

```
runConfiguredCase
├── Preflight::parseCommandLine     (Phase-1)
├── ModelLoader::read
├── ModelLoader::configureOutput
├── inspectCase  →  ResolvedSimulationSystem   (compile)
├── Preflight::validateProgram      (Phase-1 fail-fast)
├── Environment::build              (Phase-2: mesh/MPI/storage/IBM runtime)
├── Environment::validate           (Phase-2: collective 校验)
└── execute                         (只跑 CompiledSolvePlan)
```

分支语义保持逐条一致（见 §7 与回归测试）：

- `createMesh` 在 `inspectCase` 之前短路到 `buildMeshOnly`；
- multi-patch 路由仍发生在 empty-dims / IBM collective 之前；
- Ghost IBM 只在 MPI 分支作为 source-zone preprocessor 注入；
- variational forcing 多 patch 仍在 build 阶段 fail-fast 拒绝；
- IBM `setup` / `usesForcing` 仍先本地记录，再由 `validate()` 做
  `allRanksAgree` 归约。

## 5. Public API before / after

Before：`SF_run.cpp` 只暴露 `SF::Application::runConfiguredCase(int,char**)`，
其余逻辑全部内联为局部变量与匿名 lambda。

After：

```cpp
// SF_preflight.h
namespace SF::Application::Preflight {
    struct RunRequest { ok, usageRequested, error, exitCode,
                        casePath, initialOutputOnly, stepLimit };
    RunRequest parseCommandLine(int argc, char* argv[]);
    int  reportUsage(const RunRequest&, const char* programName);
    bool validateProgram(const System::ResolvedSimulationSystem&);
}

// SF_environment.h
namespace SF::Application {
    enum class DomainKind { Single, DistributedMultiPatch };
    struct ExecutionEnvironment { ... };

    namespace Environment {
        bool build(ExecutionEnvironment&, const CaseConfig&,
                   const CaseInspection&, int&, char**&, ResultWriter&);
        bool validate(ExecutionEnvironment&);
        int  buildMeshOnly(const CaseConfig&, int&, char**&, ResultWriter&);
    }
}
```

`runConfiguredCase` 对外签名不变；新增符号全部位于
`SF::Application::{Preflight,Environment}` 命名空间，不与既有符号冲突。

## 6. Dependency before / after

Before：`SF_run.cpp` 直接 include mesh / MPI / IBM / boundary / runner /
inspection / output 等十余个头文件，单文件承担全部依赖。

After：

```
SF_run.cpp        → SF_preflight.h, SF_environment.h, SF_runners.h,
                    SF_systemPrinter.h, SF_inspection.h, SF_model.h
SF_preflight.cpp  → SF_resolvedSimulationSystem.h, SF_log.h
SF_environment.cpp → mesh / MPI / IBM / boundary / inspection / output
```

依赖方向仍符合主干：

```
Application/Composition → Resolved System → ... → Infrastructure backend
```

`SF_preflight` 不依赖 mesh / MPI / IBM / boundary；`SF_environment` 不依赖
`SF_run`。没有新增对 solver 数值实现的依赖。

## 7. MPI / parallel semantics 是否变化

**未变化。** 关键点逐条核对：

- `ParallelContext` 的 distributed-context 请求条件由
  `linearMPI = solver == PressureBased` 改为
  `needsDistributedContext = parallel.enabled || solver == PressureBased`。
  这是同一布尔值的身份重命名：HYPRE backend 直接使用 `MPI_COMM_WORLD` 并
  调用 `HYPRE_Initialize()`，pressure-based 串行 case 仍需先 `MPI_Init`。
  请求条件数值等价，且 `&& !createMesh` 的语义由 `buildMeshOnly` 在 build
  之前短路、并以 `false` 构造 `ParallelContext` 来保持。
- `if (parallel.enabled && parallel.active())` 分支条件不变，因此
  pressure-based 串行 case 仍进入串行 `setupComplexMesh` 分支，与旧行为一致。
- `resolveSingleFieldLayout` 的 `parallel.allRanksAgree(...)` 仍是所有 rank
  按相同顺序执行的 collective，避免 single/multi 路由发散造成的 MPI 死锁。
- multi-patch 早返回仍跳过 single-field 的 empty-dims 与 IBM collective。
- `validate()` 的 IBM `allRanksAgree` 归约顺序与旧实现一致。
- canonical face / GlobalDof 的 COPY / SUM 语义未触碰。

## 8. Numerical semantics 是否变化

**未变化。** 本次未修改任何数值路径：dt、stage times、mass/momentum/energy
装配、rho/p extrema、norms、canonical face、GlobalDof 均未触及。
Sod 短步进（3 步）输出的 dt 序列与预期一致。

## 9. Regression tests

| # | 命令 | 结果 |
|---|---|---|
| 1 | `mpiexec -np 4 ./build/sonicSolver --steps 3 test/Sod/sodCase` | ✅ exit 0，3 步 dt 正常 |
| 2 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowDFMAugmentedLagrangian` | ✅ exit 0 |
| 3 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowDFMExplicitSelfPropelled` | ✅ exit 0 |
| 4 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowDFMFractionalStepSelfPropelled` | ✅ exit 0 |
| 5 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowDFMImplicitSelfPropelled` | ✅ exit 0 |
| 6 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowFictitiousDomain` | ✅ exit 0 |
| 7 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowFullyImplicitDLM` | ✅ exit 0 |
| 8 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowGhost` | ✅ exit 0 |
| 9 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowPeskin` | ✅ exit 0 |
| 10 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowVelocityForcing` | ✅ exit 0 |
| 11 | `./build/sonicSolver --steps 2 test/IBM/cylinderFlowVelocityForcingBP` | ✅ exit 0 |
| 12 | `./build/sonicSolver --initial-output test/IBM/cylinderFlowGhost` | ✅ exit 0 |
| 13 | `./build/sonicSolver test/createMesh` | ✅ exit 0 |
| 14 | `mpiexec -np 4 ./build/sonicSolver test/IBM/IBMCase` | ❌ exit 1（**pre-existing**，见 §10） |

## 10. Deferred issues

1. **IBMCase（PIMPLE + HYPRE 多 patch）失败为既有问题，非本次引入。**
   用重构前的原始 `SF_run.cpp` 重建二进制后复测，得到完全相同的错误：

   ```
   Fatal: Resolved unknown 'rho' expects storage 'conservative' on every
   participating patch; expected=2, actual=10.
   ```

   报错位于 solver/system 层的 storage 解析（`packed-distributed:conservative`
   vs 期望 `conservative`），与 composition-root 无关，不属于本次任务范围。

2. **`ExecutionDomain` 尚未实体化为独立类。** 本次用 `DomainKind` 枚举 +
   `ExecutionEnvironment` 结构体表达 topology；`SerialDomain /
   DistributedSinglePatchDomain / DistributedMultiPatchDomain` 的类化拆分
   （storage/execution topology 接口）留待后续，当前 runner 仍由
   `execute()` 按 `env.domain` 二路分发到 `runSingleField / runMultiPatch`。

3. **`runSingleField / runMultiPatch` 尚未合并为单一 `Run::execute`。**
   二者仍是 storage-layout 命名的入口；长期应下沉为
   `Program + Environment → Run::execute` 的单一执行器，single/multi 只存在于
   infrastructure/runtime adapter。本次为保持行为等价，只把分发收敛到一处。

4. **`SF_preflight::validateProgram` 目前只检查 `RuntimeStatus::Unsupported`。**
   Plan / Transformer / OpId provider / IBM+Eulerian 组合 / pressure-turbulence
   policy 等数学级校验仍在 `inspectCase`（`ResolvedSimulationSystem` 构建期）
   抛出，未逐一迁入 Phase-1 显式 report。迁移需先厘清 `System::RuntimeReport`
   可承载的字段，避免制造第二套 authority。

5. **`linearMPI` 身份判断虽已改名，但 pressure-based 仍以
   `solver == PressureBased` 作为 HYPRE 需要 MPI context 的代理。** 长期应改为
   `program.requirements.distributedLinearAlgebra` 这类 capability 声明，
   由 Environment 校验，而不是由 solver 枚举推断。
