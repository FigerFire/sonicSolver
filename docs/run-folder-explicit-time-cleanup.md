# Run Folder Responsibility Cleanup 与 Explicit Time 解耦报告

日期：2026-09-18

## 1. Scope

本轮只清理 application run 目录职责，并把显式时间积分从 density formulation 命名中解耦。未迁移 density lifecycle 到 `CompiledSolvePlan`，未修改 CFD 公式、stage 系数、stage time、边界/halo 顺序或 MPI 语义。

基线分支为 `equation-runtime-authority-20260916`，基线提交为 `5c915181e64ee4024dd68c75f3160a2fd6eb7292`。开始时工作树已有大量未提交的架构迁移文件；本轮保留这些内容，没有 reset、clean 或回滚用户修改。本轮分支为 `refactor/run-explicit-time-cleanup-20260918`。

## 2. Before directory responsibilities

`app/application/run/` 同时包含顶层运行入口、physics execution composition、equation/interface coupling、IBM interface adapter、runtime config、日志格式化和输出字段构造。`solver/algorithm/SF_densityBasedTime.*` 实际实现 Euler/SSPRK3/RK4 stage orchestration，却以 primary formulation 命名。

## 3. After directory responsibilities

| Before | After | Responsibility |
| --- | --- | --- |
| `solver/algorithm/SF_densityBasedTime.*` | `solver/algorithm/time/SF_explicit.*` | generic explicit integration |
| `app/application/run/SF_equationCoupling.*` | `solver/equation/coupling/SF_equationCoupling.*` | equation coupling |
| `app/application/run/SF_interfaceCoupling.*` | `solver/equation/coupling/SF_interfaceCoupling.*` | interface equation coupling |
| `app/application/run/SF_services.*` | `app/application/adapters/SF_ibmAdapters.*` | IBM adapters |
| `app/application/run/SF_runtime.*` | `app/application/model/SF_runtimeConfig.*` + `app/application/output/SF_report.*` | runtime configuration 与报告格式分离 |
| `app/application/run/SF_output.*` | `app/application/output/SF_fields.*` | output field composition |
| `app/application/run/SF_singleFluid.cpp` | `app/application/execution/SF_singleFluid.cpp` | execution composition |
| `app/application/run/SF_eulerian.cpp` | `app/application/execution/SF_eulerian.cpp` | execution composition |
| `app/application/run/SF_multiPatch.cpp` | `app/application/execution/SF_multiPatch.cpp` | execution-domain composition |
| `app/application/run/SF_executionBuilder.cpp` | `app/application/execution/SF_executionBuilder.cpp` | compiled execution binding |
| `app/application/run/SF_executionAssemblers.h` | `app/application/execution/SF_executionAssemblers.h` | execution assembly declarations |
| `app/application/run/SF_runners.h` | `app/application/execution/SF_runners.h` | execution composition entry declarations |

没有增加 Manager、Coordinator、ServiceLocator 或 forwarding wrapper。

## 4. DensityBasedTime -> Explicit Time

显式积分实现位于 `src/solver/algorithm/time/SF_explicit.{h,cpp}`，namespace 为 `SF::Time::Explicit`。`CompressibleAlgorithm` 通过该路径调用 `advance`。迁移只改变文件位置、include、namespace 和调用名；对基线源码与新实现进行逐行差异核对后，Euler、SSPRK3、RK4 的系数、state snapshot、stage combination、publish 和 validation 代码均未改变。

`methods/numerics/time/SF_Euler.*`、`SF_RK4.*`、`SF_time.*` 仍存在。它们包含单 Field 更新 kernel 及早期 callback-based stage driver；`solver/algorithm/time/SF_explicit.*` 则负责当前多 patch/registered-variable 显式 stage orchestration。因此目前仍有两层时间推进 authority 重叠。本轮按要求只记录，不合并。后续应在 density lifecycle 迁入 compiled plan 后，让 plan 拥有 stage control，`Time::Explicit` 与低层 kernels 只保留数值更新职责。

## 5. Equation coupling relocation

equation coupling 和 level-set/interface coupling 已移入 `solver/equation/coupling/`，namespace 为 `SF::Equation::Coupling`。输入仍是 `StateBundle`、phase/interface fields 与现有 model；输出仍是原 RHS/source/closure contribution。没有改变 source、halo、commit、reinitialization 或 interface 数值公式。

## 6. IBM adapter relocation

原 generic `SF_services` 中的 IBM adapters 已移入 `app/application/adapters/SF_ibmAdapters.*`，namespace 为 `SF::Application::Adapters`。它们仍只把 concrete IBM implementation 适配为 solver immersed interfaces，不选择 runner、不拥有 timestep、不创建 communicator，也不改变 IBM 方法。

## 7. Runtime config/report split

`SF_runtime.*` 被拆为：

- `model/SF_runtimeConfig.*`：`CaseConfig -> Mesh/IBM runtime config`。
- `output/SF_report.*`：time/run-control/solver configuration 的人类可读格式化。

report 代码不写 solver state，也不参与 runtime branch selection。

## 8. Output relocation

VTK scalar field view 构造迁移到 `output/SF_fields.*`。`ResultWriter` 调用、字段选择和数据内容未改变。`app/application/model/SF_output.*` 是 case model 中既有 output configuration，和已删除的 `run/SF_output.*` 不是同一职责，因此保留。

## 9. Execution directory

single-fluid、multi-patch、Eulerian 及 execution builder 迁移到 `app/application/execution/`。这些文件继续组装 physics object、`StateBundle`、providers、runtime services 和 compiled program；本轮没有拆分大文件，也没有改写 numerical lifecycle。为避免与现有 `SF::Execution` runtime namespace 冲突，现有 `SF::Application::Runners` namespace 保持不变。

## 10. Deleted obsolete files

旧位置已删除，且没有 forwarding wrapper：

```text
solver/algorithm/SF_densityBasedTime.{h,cpp}
app/application/run/SF_equationCoupling.{h,cpp}
app/application/run/SF_interfaceCoupling.{h,cpp}
app/application/run/SF_services.{h,cpp}
app/application/run/SF_runtime.{h,cpp}
app/application/run/SF_output.{h,cpp}
app/application/run/SF_singleFluid.cpp
app/application/run/SF_eulerian.cpp
app/application/run/SF_multiPatch.cpp
app/application/run/SF_executionBuilder.cpp
app/application/run/SF_executionAssemblers.h
app/application/run/SF_runners.h
```

## 11. Include/build updates

algorithm、equation 和 application CMake source lists 已指向新路径；application/environment、inspection、CLI 和 execution includes 已更新。架构守卫新增以下约束：`run/` 只允许四个 orchestration 文件；旧路径必须不存在；`SF_explicit` 不得依赖 application/run；equation coupling 不得依赖 application execution/run；源码中禁止旧 `DensityBasedTime` 标识。

stale-reference 搜索结果：旧 density time、run coupling、services、run runtime 与 run output 引用均为 0。唯一 `#include "SF_output.h"` 来自 `app/application/model/SF_output.cpp`，它引用同目录的 case output configuration，不是旧 run output builder。

## 12. Numerical invariants

以下内容保持不变：Euler/SSPRK3/RK4 tableau、stage time、stage state combination、`DensityBasedRHS`、WENO/TENO、flux split、pressure correction、IBM/ILW/KKT/HYPRE、boundary/halo order、canonical face COPY、GlobalDof SUM 和 `StateBundle` semantics。

ownership before/after：显式 RK orchestration 从 density-named algorithm file 迁入通用 time algorithm namespace；state 和 workspace owner 未变。execution order before/after 完全相同。public numerical API 只发生 namespace/name 迁移。MPI/parallel semantics 与 numerical semantics 均未改变。

## 13. Build/test results

### Build

- `cmake --build build --parallel 4`：通过，`sonicSolver` 完整链接。
- `cmake -S . -B build-tests -DBUILD_TESTS=ON -DBUILD_GUI=OFF`：通过。
- `cmake --build build-tests --parallel 4`：通过，测试与 production targets 均完成，无 compile/link failure。

### Tests

- `pisoArchitecture`：通过。
- `pisoNumericalRegression`：通过。
- `numericalFluxContract`：失败，仍是迁移前已知的 reconstructed Rusanov constant-stencil 偏差（例如 component 0 为 36.0551，对比 physical 36）。本轮未修改该通量公式或 tolerance。
- architecture checker：通过；618 个 source files，allowlisted dependency edges 为 10，没有新增 dependency violation。
- serial velocity-forcing IBM one-step smoke：通过；`max|Ju-Us|=1.9613e-12`，Schur CG 9 iterations，relative residual 报告为 0。

### Pre/post numerical regression

所有 case 使用相同 configuration，各推进一步：

| Case | Scheme/path | dt/final time | Final VTS SHA-256 before/after | Result |
| --- | --- | --- | --- | --- |
| single Sod WENO7 | Euler + WENO7 + Steger-Warming | `6.681531e-04` | `1ab9ccf4cfb3083459de82fa2312e9eaafc0dbe4be896cdf354d9aafd9d55d07` | bitwise identical |
| single Ghost IBM | RK4 + WENO5 + Lax-Friedrichs | `4.392680e-04` | `b9d70a4855121102dfb351841332986ae12080862f0facc2aa696d194f329b0e` | bitwise identical |
| pressure PISO | compiled pressure plan | `1.659315e-06` | `11cba96d719f1b00ced3a01237ba489c636ffec13de3c0c516144af47ca88040` | bitwise identical |

Ghost RK4 final closure 仍为 `min(rho)=1.22383`、`min(p)=101208`。诊断文本的迁移前后比较无差异。该证据覆盖 Euler 与 RK4；SSPRK3 实现经源码 diff 证明公式体未修改，仓库当前没有独立 SSPRK3 case。

## 14. Remaining run/ files and why they belong there

`app/application/run/` 最终只包含：

```text
SF_run.h
SF_run.cpp
SF_runFlow.h
SF_runFlow.cpp
```

`SF_run` 提供 top-level command orchestration；`SF_runFlow` 连接统一 time driver 与已构造的 execution callback。二者不实现 equation coupling、IBM adapter、output field composition 或 explicit integrator。

## 15. Remaining large files

- `execution/SF_singleFluid.cpp`：462 行，仍包含 single-fluid composition 和 legacy density execution wiring。
- `execution/SF_multiPatch.cpp`：326 行，仍包含 multi-patch execution composition。
- `execution/SF_eulerian.cpp`：171 行，保留 Eulerian composition。
- `algorithm/time/SF_explicit.cpp`：494 行，包含 Euler、SSPRK3、RK4 的冻结 stage orchestration。
- `methods/numerics/time/SF_time.h`：402 行，保留旧 callback-based time update helpers。

本轮没有为拆文件新增 Manager/Adapter/Interface。

## 16. Known numerical issues untouched

`numericalFluxContract` 的 reconstructed Rusanov constant-stencil mismatch 仍存在；它属于既有 high-order numerical issue。本轮没有调整公式、fallback 或 tolerance。density lifecycle、`IFlowAlgorithm`、`makeFlowAlgorithm`、BindingPass、AssemblyPlanRegistry、Ghost IBM MPI 和新 physics capability 也保持原状态。

## 17. Next task

下一轮建议为 **Explicit Time Plan Authority Migration**：把 `CompressibleAlgorithm::stepDensity` 的 legacy global stage lifecycle lowering 到 `CompiledSolvePlan -> StageLoop -> OpRegistry`，并让 `Time::Explicit` 只保留数值 stage helper。届时再评估删除 `IFlowAlgorithm` / `makeFlowAlgorithm`。本轮没有提前实施这些内容。
