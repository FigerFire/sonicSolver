# Phase 28A 修订执行计划

日期：2026-09-26。采用用户同意的评审修订，原计划的旁支清理不作为分布式压力验收前置条件。

## 数值与 authority 审计

本轮属于 execution/linear-algebra/runtime 层的“改变怎么解”，不新增物理 term。调用链保持 ExecutableSystem → CompiledSolvePlan → PlanExecutor → OpRegistry → PressureOperators。StateBundle/已实现 U、p view 拥有 physical state；PressureOperators 拥有局部数值 workspace；runtime/backend 拥有通信和编号映射。

当前 `assembleMomentum` 的 outer>1 分支采用 U_n + dt R(U_k)，不是严格 forward Euler 的 U_n + dt R(U_n)。保留现有算术，明确记录 Legacy fixed-time spatial-residual reevaluation；不冻结 RHS，不宣称已经实现通用 backward Euler。首先以线性衰减 oracle 固定这个事实。TimeRecipe 继续决定物理 step topology，既有 provider 的时间层局限必须在 explain 可见。

## 基线

HEAD 537889b05234c34217aa1a8eff7ab7844b69b4b9；它不是完整 Phase27B 源码。完整 dirty source、diff、manifest、CMakeCache 与 binary 已保存到 `/private/tmp/sonic-phase28a-baseline-20260926-154421`。修改前 build 成功、architecture check 2 条既有 allowlist；完整 CTest 和 production 对照单独记录结果。

## 顺序

1. 源码/可执行基线、时间层 oracle、MPI substrate 预检。
2. 将 pressure fragment 的 physical prepare/commit 移到 generic planner composition；单次 time.commit，保持 state/diagnostics 相对顺序。
3. 必要的 capability/explain 修正；保留未实现组合 Unsupported。
4. global dt 在创建 rAU 前归约；U/p/pPrime/rAU halo；canonical face response 与 flux COPY。
5. GlobalDofSystem → DistributedLinearSystem → HYPRE；稳定 reference entity；owner-only rows。
6. owner-only correction/global diagnostics；分别归约 delta 与 scale 后归一化。
7. matrix/halo/face/reference 测试，PISO MPI，SIMPLE/PIMPLE MPI，再完整 production regression。
8. 中文报告及两个设计文档/src/ARCHITECTURE.md 同步更新。

## 边界和验证

- processor interface 不是 physical BC；不平均 shared state/face flux。
- serial 同执行配置保持原记录精度；MPI 按 physical DOF 对齐，容差在执行前设定。
- referenceCell 表达原单块 interior-cell 顺序下的实体，而非 HYPRE row。
- 预通信错误经一致检查；通信中不可恢复错误由 infrastructure fail-stop，不在 numerical code 调用 MPI。
- 首差异检查 dt → halo state → matrix/RHS/face response → correction → flux。
- 保留旧 pressure gauge-dependent convergence normalization，明确为 Legacy；不偷偷改变 baseline。
- 空 owner-row rank/subcommunicator/真实多块能力须明确 Supported/Unsupported。
- 不实现 implicit/IMEX/新 source/IBM pressure。ProviderCatalog 全面迁移、model cycle、Field/IO/ILW 大改不作为前置；无实际消费者时不添加 IMEX 空类型。
