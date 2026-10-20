# Phase 28A — 分布式压力算子与架构收紧

日期：2026-09-27。中文报告；文档仅保存在本地，本轮未提交或推送 GitHub。

## 1. 结论与范围

本轮实现了**单逻辑源块分解、每 rank 一个参与 patch**的常密度 pressure operators，复用原 PISO/SIMPLE/PIMPLE Plan 和 numerical providers。已执行 serial、2-rank x/y、4-rank 对照；不是任意多块压力系统的完整实现。

保留原 production 数学。fixed-time predictor 的时间层审计发现旧 outer 更新为 U_n + dt R(U_k)，本轮将其明确标为 Legacy，而非宣称它是严格 forwardEuler。多连通二维压力零空间也保留为 Legacy；没有偷偷新增 gauge row 改变串行结果。

原工作区包含 Phase27B 未提交修改，不能用 HEAD 代替数值基线。原 HEAD 为 `537889b05234c34217aa1a8eff7ab7844b69b4b9`；完整 dirty source、diff、manifest、CMakeCache、可执行文件已冻结于 `/private/tmp/sonic-phase28a-baseline-20260926-154421`。当前 diff 同时包含这些先前修改；本报告只声明本轮增量。

## 2. 状态分类

| 状态 | 项目 |
|---|---|
| Implemented | global dt；pressure state/workspace halo；owner-only pressure rows；分布式 HYPRE solve；reference entity；canonical face response/flux COPY；global convergence；通用 physical prepare/commit；显式不适用 preset 报错 |
| Implemented、限定范围 | 单源正交块分解为 2/4 ranks，PISO、SIMPLE/PIMPLE 复用原 operations；MPI 不引入第二套 Plan |
| Interface-only | 本轮未为未来 IMEX/ProviderCatalog 添加无消费者的空接口 |
| Legacy | R(U_k) 的 fixed-time predictor；pressure norm 对 gauge 敏感；二维 inactive-z 多连通层只有单个 reference row；scalar halo 临时打包；已有模型/pressure 特化 |
| Unsupported | 多本地 patch/真实多块压力装配、空 owner-row rank、subcommunicator、通用 implicit/IMEX、pressure RK stage treatment、新常密度 IBM/model term 和非正交数值能力 |

## 3. 层级、ownership 和依赖

属于“改变怎么求方程”，不是新增物理 term 或 equation。调用链仍为：

```text
ExecutableEquationSystem / CompiledNumericalSystem
  → CompiledSolvePlan → PlanExecutor → OpRegistry
  → PressureOperators → GlobalDofSystem → DistributedLinearSystem
  → HYPRE backend
```

| 语义 | Before | After |
|---|---|---|
| 物理状态 | 已实现 U/p views、StateBundle | 不变，不复制第二份 Q |
| 物理 clock | pressure commit callback 写 state.time/step | generic `time.commit` 唯一提交；pressure callback 只提交数值状态 |
| physical prepare/commit | pressure preset fragment 内嵌 | fragment 标记 `physicalStepBody`；通用 planner 注入 |
| numerical workspace | PressureOperators 本地向量 | 不变；增加一次构造的数学实体/后端编号映射 |
| pressure matrix | 局部 SparseSystem/SolverSession | GlobalDofSystem；backend row 仅在 DistributedLinearSystem 映射 |
| reference | 局部 row 语义 | 全局稳定物理实体；非 rank-0 owner 也有效 |
| 通信 | 仅串行 pressure assumption | 稳定注入 IExecutionRuntime；MPI 仍在 infrastructure/backend |
| canonical view | 接口中的整块临时 face buffer | borrowed DistributedFieldView；仅打包通信 payload |

关键 API：`PressureOperators::bind(state, runtime)`；`HaloExchange::synchronizeCanonicalFaceFlux(views)` 替代旧 Field/vector overload。未引入 Manager、ServiceLocator 或 MPI 专用压力 solver。

关键文件：`src/solver/algorithm/pressure/SF_pressureOperators.{h,cpp}`、`src/solver/system/SF_pressureCoupling.cpp`、`SF_solvePlan.cpp`、`SF_systemBuilder.cpp`、`src/solver/algorithm/SF_singleFluidStepper.cpp`、`src/infrastructure/mpi/SF_parallelCoordinator.cpp`、`SF_haloExchange.{h,cpp}`、`backend/SF_mpiBackend.cpp`、`src/solver/linearAlgebra/backend/hypre/SF_hypreBackend.cpp`。

Architecture checker：628 source files，既有 allowlisted dependency edges **2 → 2**；没有新增 core→solver 或 numerical→MPI include。重复 header basename 仍按既有规则报告，未借本轮搬文件。

## 4. 时间语义与 execution order

`test/test_fixedTimePressure.cpp` 增加 scalar decay oracle：重复 U_n − dt a U_k 收敛到 U_n/(1+dt a)，而严格 FE 为 U_n(1−dt a)。因此不能将多次 outer 称为“同一个 FE RHS 的一致性修正”。本轮不改变这一旧数学；`SF_systemPrinter.cpp` explain 明示。

prepare → fixed-time begin → 原 momentum/correction/relaxation/flux/convergence 顺序 → pressure state commit → 唯一 time.commit。PISO 不新增 outer；SIMPLE/PIMPLE 仍由原 Plan 控制固定时间的 outer loop。未修改 RK 系数、stage time、显式 density lifecycle 或 pressure numerical formula。

MPI 中先算 local dt，再 global MIN，最后构造 dt/rho 的 rAU。所有 ranks 的 pressure mobility 使用同一个物理 dt。

## 5. 分布式装配与 reference

仅 canonical owner 写压力行；邻接 remote halo 通过稳定 GlobalDofId 参与列引用。rank 的连续行范围由 runtime 分配，数学 ID 与 backend row 不混用。局部 boundary/processor corner 依据实际 physical side masks 选择 inward donor，不能将处理器边缘当成物理墙。

referenceCell 按原单块 interior entity 顺序映射，通过全局计数/极值查找稳定实体；不是“每 rank 固定一个 row”，也不是 HYPRE row 0。测试 referenceCell=150 在 y 分区下由非零 owner 持有，并检查唯一 gauge row、参考压力和参考层绝对压力。

初始 matrix/RHS 以物理坐标逐行对齐，比较 row set、column set、coefficients 和 RHS；固定容差 1e-12。后续 HYPRE 分区相关迭代差异通过最终 U/p 和 continuity 对照，不要求跨分区浮点 bitwise。

## 6. Halo / canonical / reduction

U、p、pPrime、rAU 使用 owner COPY；pressure predictor 的 face response 和 face flux 在候选就绪后 canonical COPY，再装配。corrected face flux 保留持久 workspace。GlobalDof contribution SUM 没有改成 COPY 或 average。

max divergence 全局 MAX。outer convergence 分别 global MAX(deltaU)、MAX(scaleU)、MAX(deltaP)、MAX(scaleP)，之后求 ratio；不取局部 ratio 的 MAX。flux delta 也全局 MAX。测试 oracle 用 numerator/denominator 极值位于不同 ranks 的情况防止错误恢复。

一个 rank 的 zero RHS 不能令其提前离开全局线性求解。HYPRE exact-zero 和 zero-initial-guess near-zero 判断使用全局范数及原配置容差；没有调大 tolerance 或改预条件器。非有限 workspace 在通信前一致拒绝；异常展开由 infrastructure MPI_Abort fail-stop，避免 peers 永久等待。Application 在 environment 析构前向 stderr 打印原始异常和 rank，避免 fail-stop 隐藏原因。expected-failure 测试与正常成功测试分开。

## 7. First-difference 调查及修复

1. 初次分布式 velocity 差异来自 physical boundary 与 processor edge 相交时选错 inward donor。修复 physical-side 判定后，初始矩阵/RHS 一致，U 误差降至 1e-10 量级。
2. HYPRE near-zero RHS 在 MPI 下未复用 serial absolute tolerance 初值收敛判定。补齐全局范数，保留原 tolerance。
3. 去掉 remote scratch 后 Identifier vector API 从本地 view 推导消息尺寸，导致 MPI_ERR_TRUNCATE。单 patch/rank 使用已有直接 scalar owner-COPY 路径；最终 MPI 测试通过。
4. 长时间 PISO pressure 常数偏移来自旧 inactive-z 多连通 nullspace，不是 workspace 索引或 assembly sign 错误。未加额外 pin 掩盖它。
5. 非平凡 SIMPLE/PIMPLE 测例耗尽 outer budget 不代表控制流出错。另加静止平衡 case 检验每 timestep 一次 outer 后所有 ranks 同时 early exit；原非平凡对照保留。

## 8. 显式 preset 与 native 输入

Native reader 不再把缺省 algorithm 视作显式注册 PIMPLE。`output/SF_report.cpp` 删除从 SolverConfig 默认值重复打印的 Coupling algorithm/Pressure solve 块，实际 coupling/binding 由已存在的 resolved system report 解释。显式注册但不存在数学 constraint 的 preset 为 Invalid；可选且未声明的 request 可以 Inactive。未知 algorithm 不能因 density 标签而忽略。

Ghost 与12个同类 IBM fixture 旧 `algorithm: PIMPLE` 实际不执行；逐个运行 check 确认缺少数学 constraint 后，移除这项无效 preset，保留原 linear-solver controls。没有按 density 标签泛化拒绝合法组合。用 frozen binary 和当前 binary 读取同一归一化 fixture，Ghost20 输出 hash 一致。未用 solver-family flag 恢复 silent ignore。`tools/check_cli_init_recipe.py` 同时检查 CLI 默认 recipe 可读和显式不适用 coupling 被拒绝。

## 9. 数值基线

修改前：完整 Debug build 成功，19/19 CTest 通过。修改后：正常 build 和独立 clean Debug build 成功；22/22 CTest 通过。Clean 日志 `/private/tmp/phase28a-clean-build.log`、`/private/tmp/phase28a-clean-ctest.log`，最终 clean CTest 耗时281.78秒（与另一轮 IBM smoke 并行运行）。链接器有既有 duplicate static-library warning，无 compile/link failure。

| Production 对照 | 结果 |
|---|---|
| 4-rank Sod，20 steps | diagnostic sequence、VTS hashes 与冻结 binary 相同 |
| serial Ghost，20 steps | diagnostic sequence、VTS hashes 相同 |
| Eulerian，2 steps | VTS hashes 相同 |
| legacy PISO，1 step | diagnostic sequence、VTS hashes 相同 |
| constant-density PISO，2 steps | diagnostic sequence、VTS hashes 相同 |

旧 PISO final SHA-256：`4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`。详见 `/private/tmp/phase28a-production-resume.log`。哈希覆盖输出字段的记录精度，不宣称所有内部浮点中间量逐 bit 相同。

## 10. 新增 MPI pressure 回归

`tools/check_distributed_pressure.py` 固定 U absolute tolerance=1e-8、p shape/参考层 absolute tolerance=1e-6、matrix/RHS tolerance=1e-12。没有为失败调整比较 tolerance。16×16×1 cell 网格；channel 2048 steps 到 t=10，并检查 Poiseuille analytical L2≤1.2e-3。

| Case / partition | U max delta | p shape delta |
|---|---:|---:|
| cavity 2×1 | 6.728e-11 | 7.146e-8 |
| cavity 1×2 | 7.072e-11 | 8.396e-8 |
| cavity 2×2 | 3.543e-10 | 1.419e-7 |
| channel 2×1 | 4.912e-10 | 1.999e-8 |
| channel 1×2 | 2.639e-10 | 1.893e-8 |
| channel 2×2 | 4.797e-10 | 1.609e-8 |
| SIMPLE 2×1，8 steps | 1.421e-11 | 1.164e-10 |
| PIMPLE 2×1，8 steps | 1.421e-11 | 1.164e-10 |
| 静止平衡 SIMPLE/PIMPLE 2×1 | 0 | 0 |

所有 case 的 dt/physical clock 一致；初次矩阵/RHS 通过；全局 outer 次数一致。MPI 在 host 环境运行，不是 sandbox launch failure。测试输出保留在 `/var/folders/nv/1ngf5p455x5gy26g38l7dw7w0000gn/T/sonic-distributed-pressure-jq2k74v2`。

`test/test_distributedPressure.cpp` 另验证非零 rank 的 dt 限制、global ratio、某 rank 零 RHS/另一 rank 非零、全局零/近零 RHS 和单 rank 异常 fail-stop。

## 11. IBM 覆盖边界

本轮没有改 IBM 数学或 pressure IBM capability。`sonicsolver-ibm-validation` 方法矩阵用于防止以 Ghost 代替 forcing/KKT 验收。

三个 implicit/augmented fixture 在 frozen/current 两端均被拒绝：配置了 diffusion recipe，但 equation system 没有 consuming term/operation。这是既存 Invalid fixture，不是 KKT numerical pass，也不能由此推断其底层算法已验证。其余七个串行方法在去掉既存 inactive density PIMPLE 后配置检查一致通过。

完整5秒 self-propelled 扩展试跑已停止，保留输出；它不是本轮接受证据。短程 smoke 的最终结果如下，明确不替代完整时间验收：serial `DFMExplicitSelfPropelled`、`DFMFractionalStepSelfPropelled`、`FictitiousDomain`、`Peskin`、`VelocityForcing`、`VelocityForcingBP`，以及 `FictitiousDomain`/`VelocityForcing` 的 MPI-2 和 MPI-4，共10组 frozen/current 输出 VTS 哈希相同。Ghost 已由20步 production 对照覆盖。记录：`/private/tmp/phase28a-ibm-final.log`。仓库当前缺少 CMake 可选引用的 `test_distributedKKT.cpp`、`test_sparseCanonicalCopy.cpp`、`test_distributedSurfaceSchur.cpp`，所以22项 CTest不包含这些测试。

## 12. 性能与分配

编号和 owner row 集合在 bind 时建立，不每 stage 重建。canonical face 直接借用 workspace view，不复制第二份完整 FluxField。一个 patch/rank 的 transient scalar exchange 不再为每个远端 patch 分配全尺寸 scratch。仍有原 scalar component packing、linear algebra buffers 和旧 flux iterate snapshot；不声称本轮消除了所有 per-stage allocation。没有增加第二个物理 state、clock、registry 或 MPI pressure lifecycle。

## 13. Architecture guards 与 deferred

守卫禁止 PressureOperators 直接 MPI、SparseSystem/backend row authority、preset-name runtime branch，并要求 GlobalDof、global dt、全局分母归约和 canonical COPY contract。已有 planner 不得持有 pressure OpId/preset 分支守卫仍通过。

后续独立工作：多本地 patch/真实多块；空 rank/subcommunicator；多连通 pressure gauge；gauge-invariant convergence normalization；fixed-time 时间 recipe 数值定义；通用 implicit/IMEX；provider catalog全面收敛；未消费模型依赖；pressure IBM、非正交实现。没有进入 Phase28B，没有修改 density RK、native IO 架构或 Field storage。

## 14. Self review

- physical state/clock/term authority 未新增第二份；Plan 和 runtime provider 一致。
- numerical formula、RK/stage timing、已有串行 boundary/halo order、residual sign 未改。
- COPY/SUM/canonical owner 语义未改；新增的是之前不具备的压力分布式执行能力。
- 未新增 clamp、降阶、solver/EOS fallback；显式无效设置变为 fail-fast。
- 没有把零空间常数偏移隐藏为“绝对压力一致”；没有把短程 IBM smoke 写成完整验收。

## 15. IBM 追加诊断证据

MPI Peskin-2/4 的原冻结 executable 在180秒内未退出；当前 executable 立即 fail-stop，不能计作 numerical pass。为区分旧故障与本轮引入，复制完整冻结源码，仅修改 `SF_application.cpp` 的 stderr 异常输出及 `SF_mpiBackend.cpp` 的异常 fail-stop（文件差异清单 `/private/tmp/phase28a-old-diagnostic-diff.log`），独立重建。旧源码 MPI-2 和当前版本均明确报告：

```text
Peskin interpolation produced a singular mass diagonal.
```

证据：`/private/tmp/phase28a-old-peskin-diagnostic.log`、`/private/tmp/phase28a-peskin-diagnostic.log`。旧数值代码未变，未加入 mass clamp 或跳过 marker 的 fallback。本轮修复了异常被 MPI_Finalize 等待吞掉的诊断/退出行为，**没有修复 Peskin 分布式数学**。MPI-4 观察到原版超时/当前失败，但旧源码诊断重建只额外复现 MPI-2，不能夸大覆盖。
