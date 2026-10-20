# 高阶数值回归隔离追踪

日期：2026-09-14。此阶段冻结 Phase 4C、IO consistency cleanup、Field 后续拆分和 IBM 重构；没有修改 CFL 默认值、WENO/TENO 权重、Steger–Warming、physical-state validation、canonical COPY、GlobalDof SUM 或 Phase 4B workspace ownership。

## 1. 对照版本与可复现性

目标是以相同 case、mesh、CFL、Euler、初值和 rank 比较旧版和当前版本的第一物理步。

| 项目 | old | current |
| --- | --- | --- |
| 可执行文件 | `build-phase1a-tests/sonicSolver`，mtime 2026-09-11 10:42 | `build/sonicSolver` |
| Phase 4B workspace | Field-owned flux/residual | `PatchWorkspace` owns flux/residual |
| single WENO7 contract | `EquationSetRusanov` | `HighOrderPerfectGas` |
| multi TENO5 contract | 历史报告记录为 `HighOrderPerfectGas` | `HighOrderPerfectGas` |
| Git revision | 不可用；`.git` 指向不存在的 `/Users/lpf/.sonicSolver-git/sonicSolver.git` | 当前工作目录 |

没有可用的 pre-Phase3/Phase4 Git object、worktree 或完整旧源码快照。`/private/tmp/sonic-io-old` 只含空目录；不能作为旧版构建。9 月 11 日的旧二进制可用于运行，但无法加逐阶段 trace。

因此，**single WENO7 不能构成 old/current 同方法比较**：旧版本名义上配置 WENO7 + Steger–Warming，实际构造 `CompressibleAlgorithm` 时选择 `EquationSetRusanov`。这与 Phase 3C 文档中的旧 composition 一致。其“通过到 t=5”的历史基线是 EOS-Rusanov 基线，不是已覆盖真实 WENO7 的基线。

4-rank TENO5 的旧版和当前版都是实际 high-order contract，才是有效的 workspace/lifecycle 对照；但本环境的 `mpirun -np 4` 需要主机 PRTE socket。主机运行申请被自动审批以账户用量限制拒绝，沙箱则报 PRTE bind/socket failure。该项没有被记作数值失败，也没有以不同 rank 替代。

## 2. 新增的只读单步追踪

新增 [SF_highOrderTrace.h](/Volumes/SSD_LPF/sonicSolver/sonicSolver/src/solver/algorithm/SF_highOrderTrace.h:1)，默认关闭；仅当 `SF_HIGH_ORDER_TRACE=1` 且 `state.step == 0` 才输出。它只读取现有 `Field`、`FluxField`、`Residual`，不会写数值 state。

- [SF_compressible.cpp](/Volumes/SSD_LPF/sonicSolver/sonicSolver/src/solver/algorithm/SF_compressible.cpp:271) 记录初始 Q、boundary-prepared Q、legacy/EOS CFL 和 gamma。
- [SF_densityBasedRHS.cpp](/Volumes/SSD_LPF/sonicSolver/sonicSolver/src/solver/algorithm/SF_densityBasedRHS.cpp:115) 记录 candidate flux、方向 residual、canonical COPY 前后、local residual 和 GlobalDof SUM 后 residual。
- [SF_densityBasedTime.cpp](/Volumes/SSD_LPF/sonicSolver/sonicSolver/src/solver/algorithm/SF_densityBasedTime.cpp:309) 记录 Euler/SSPRK3/RK4 每次 Q 更新后的状态。
- [SF_parallelCoordinator.cpp](/Volumes/SSD_LPF/sonicSolver/sonicSolver/src/infrastructure/mpi/SF_parallelCoordinator.cpp:25) 仅在 trace 开启时验证每个本 rank block 恰有一个 patch/workspace 映射，并打印 block、owner rank、patch、`Field*`、`FluxField*` 和 `Residual*`；canonical 与 GlobalDof barrier 都检查。

使用同一个临时 `endStep: 1` WENO7 case 验证 trace 不干扰数值：开启 trace 与关闭 trace 的一步 VTK 文件逐字节相同。

## 3. 第一物理步：serial WENO7 Sod

case 是 `test/Sod/sodCase_weno7` 的无结果目录副本，只把 runtime 改为 `endStep: 1`、`writeInterval: 1`。mesh、初值、CFL=0.5、Euler 和 WENO7/Steger–Warming 均未改。

当前 trace 的关键值如下（interior 统计遵循求解器的 `forFluidInterior`；boundary checkpoint 统计完整 storage/ghost）。

| Checkpoint | old | current | Same? |
| --- | --- | --- | --- |
| initial Q / initial VTK | 旧/当前初始 `.vts` 字节完全相同 | `sum(rho)=2367.75`，`sum(rhoE)=5.785500000000001e8`，`p=[10000,100000]` | 是 |
| CFL dt | `6.681531047810609e-4` | legacy=`6.681531047809716e-4`；EOS=`6.681531047810609e-4`；active=`6.681531047810609e-4` | 是，roundoff 量级 |
| boundary Q / ghost | 无旧版内部 trace | 完整 storage `sum(rho)=17653.75`；ghost-only `sum(rho)=15286`、`sum(rhoE)=3.7352e9`，左右 ghost 与对应物理端一致 | 未采集 |
| gamma | 无旧版内部 trace | config=1.4，EquationSet=1.4，传入 `divDispatch`=1.4 | 未采集 |
| candidate face flux | 不可比；旧版实际 Rusanov | XI discontinuity 左侧 `F=(11.869151331370059,5500,0,0,4258741.7198523311)`；右侧 `F≈(-1.0057e-10,999.999999965916,0,0,2.6438e-5)` | 不适用 |
| directional residual | 不可比；旧版实际 Rusanov | XI mass `L1=451.02775062204245`，XI energy `L1=1.6183218536205542e8`；discontinuity-left cell 的 XI `F=(11.869151331370059,5500,0,0,4258741.7198523311)` | 不适用 |
| canonical flux | single patch 无共享 face；前后统计相同 | candidate 与 COPY 后 FluxField 统计逐项相同 | 是（single 无 canonical mutation） |
| local residual | 无旧版内部 trace | mass `L1=9020.555012440855`，energy `L1=3.2366437072210212e9` | 未采集 |
| global residual | single patch 无 GlobalDof shared contribution | `globalMask=0` | 不适用 |
| Q after Euler update | old 是 EOS-Rusanov；新是 high-order，必然不同 | `rho=[0.1249999999953,1.0000000000004]`，`p=[9999.99999956,100000.00000003]` | 否，方法不同 |

旧/当前一步后 VTK 的最大绝对差也证明这个差异发生在对流方法选择，而非输入：`rho=3.0070896868e-2`、`rhoU=6.6785435157e-10`、`rhoE=3.2991502581e2`、`p=1.2533557114e2`。该差异不能用来判定 Phase 4B regression。

## 4. 4-rank TENO5 Sod

历史 Phase 3A 记录的 first `dt` 是 `6.681531e-4`，并记录 500 步到 `t=2.435715e-1` 的正常 high-order multi-patch 结果。当前真实 TENO5 + Steger–Warming + HighOrderPerfectGas 在约 `t=0.0479` 出现负压，但此环境无法启动主机 MPI 以运行新增 trace 或生成 old/current 的同一步 VTK。

| Checkpoint | old | current | Same? |
| --- | --- | --- | --- |
| initial Q | 待主机 MPI trace | 待主机 MPI trace | 未判定 |
| dt | 历史首步 `6.681531e-4` | 待主机 MPI trace | 未判定 |
| boundary Q | 待主机 MPI trace | 待主机 MPI trace | 未判定 |
| candidate flux | 待主机 MPI trace | 待主机 MPI trace | 未判定 |
| directional residual | 待主机 MPI trace | 待主机 MPI trace | 未判定 |
| canonical flux / COPY | 待主机 MPI trace | 新断言已就位，待执行 | 未判定 |
| local residual | 待主机 MPI trace | 待主机 MPI trace | 未判定 |
| GlobalDof SUM | 待主机 MPI trace | 新断言已就位，待执行 | 未判定 |
| Q after stage | 待主机 MPI trace | 待主机 MPI trace | 未判定 |

主机恢复后应执行：

```sh
SF_HIGH_ORDER_TRACE=1 mpirun -np 4 ./build-phase1a-tests/sonicSolver run /private/tmp/high-order-old-teno5
SF_HIGH_ORDER_TRACE=1 mpirun -np 4 ./build/sonicSolver run /private/tmp/high-order-new-teno5
```

注意旧二进制没有新 trace；要完成逐 checkpoint 对照，需要恢复对应旧 revision 后以同一 trace 编译，或先以旧二进制的 step-1 VTK 对当前 step-1 VTK 做输出级对照，再回到旧 revision 加 trace。不能以 serial multi-patch 或不同 YAML 替代。

## 5. 高阶稳定性诊断（不是修复）

所有 run 从同一 serial WENO7 Sod 副本开始；只改了本次临时 case 的指定变量，未改任何仓库默认配置。

| 临时实验 | 结果 | 含义 |
| --- | --- | --- |
| WENO7 + Steger–Warming，CFL=0.5 | `t=0.1297871` 后 `splitStegerWarming` 对 `rho=0.0754065, p=-300.189` fail-fast | 已复现问题；没有 fallback、clamp 或自动降 dt |
| WENO7 + Lax–Friedrichs，CFL=0.5 | 到 `t=0.15`，无 failure | 空间 flux splitting/耗散是重要嫌疑，但不是结论 |
| WENO7 + Steger–Warming，CFL=0.125 | 到 `t=0.15`，无 failure | CFL/temporal positivity robustness 是重要嫌疑，但不是结论 |

## 6. 当前结论与分类

不能诚实地在 A（workspace）、B（lifecycle）或 C（high-order numerical-method defect）之间做最终分类：

1. single WENO7 的历史“正确”基线没有执行 WENO7，因此它排除了将该基线视作 old high-order evidence；它不能证明 A 或 B。
2. 当前 single trace 证明输入、gamma、CFL、workspace clear 后的 candidate flux、残差与 Euler update 已可逐阶段观察；未发现这一步的 workspace 映射或 reset 失败。
3. 真正可裁定 A/B/C 的 multi TENO5 同 rank old/current trace 仍未完成，原因是主机 MPI 被环境限制拒绝，而不是求解器 assertion。

当前仅能给出 **C 候选**：真实 WENO7 + Steger–Warming 在 CFL=0.5 失稳，而 LF 或 0.25 CFL 均通过到 `t=0.15`。这不是对高阶方法的修复建议，也不是将问题定性为 C；在完成 4-rank TENO5 逐阶段对照前，不应修改高阶格式或归责 Phase 4B。

## 7. 未改变的语义与验证

- `ninja -C build sonicSolver` 通过；构建时仍有既有的 `ninja: premature end of file; recovering` warning，但 link 成功。
- `python3 tools/check_architecture.py` 通过：allowlisted dependency edges = 10，未新增 core -> solver edge。
- trace 关闭与开启的一步 WENO7 VTK 逐字节一致。
- 未改变 canonical COPY、GlobalDof SUM、residual 符号、source clear timing、workspace allocation 条件、geometry 或 IBM storage。
- 未运行完整 CTest 或 4-rank host regression；它们不能据此标为通过。

完成 host MPI 对照前，保持 Phase 4C、IO consistency cleanup、Field 和 IBM 重构暂停。
