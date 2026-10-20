# 高阶数值回归的首差异隔离追踪

日期：2026-09-15  
范围：仅为定位 4-rank TENO5 Sod 的 first divergence 加入运行时追踪、严格调试断言和离线 VTK 对比工具。本轮没有修改 CFL、WENO/TENO、Steger-Warming、通量公式、边界/halo 顺序、canonical COPY、GlobalDof SUM、Field/PatchWorkspace ownership 或任何物理有效性规则。

> **2026-09-15 后续黑盒结论。** 旧 4-rank executable 的名义
> `TENO5 + StegerWarming` 一步输出对 `LaxFriedrichs` 替换完全不敏感，
> 而 current 明确的 `HighOrderPerfectGas` contract 对同一替换敏感。因此旧
> 4-rank TENO5 已不具备“已验证真实高阶 baseline”的资格。本报告原先提出的
> old-source step-0 checkpoint 对照暂停；详见
> [旧版 TENO5 Contract 黑盒验证与首差异前置隔离](old-teno5-contract-blackbox.md)。
> 随后的 4-rank LLDB 取证进一步确认 old active Euler RHS 直接执行
> `RusanovEOS::div`，high-order dispatch 与 splitter 均未执行；详见
> [旧版 executable 数值执行路径取证](old-executable-execution-path-forensics.md)。

## 1. 修改内容与边界

| 路径 | 修改 | 运行时影响 |
| --- | --- | --- |
| `src/solver/algorithm/SF_highOrderTrace.h` | 新增 `SF_HIGH_ORDER_TRACE_STEP`（默认物理步 `0`）的严格解析；追踪包含 `StateBundle::step`、stage 和 stage time；加入 `workspaceAfterClear()`。 | 仅在 `SF_HIGH_ORDER_TRACE=1` 且选中该物理步时读取并输出。`workspaceAfterClear()` 逐项检查 `FluxField` 与完整 `Residual` 存储均为零，发现非零即 fail fast；不写入任何数值数组。 |
| `src/solver/algorithm/SF_densityBasedRHS.cpp` | 在原有 `equations.begin(workspace.convectiveFlux, workspace.residual)` 返回后调用上述检查；建立 stage trace context。 | 清零调用、RHS 顺序和 workspace 生命周期不变。 |
| `src/solver/algorithm/SF_compressible.cpp` | 在 density 物理步起点建立 trace context。 | 只读日志状态。 |
| `src/infrastructure/mpi/SF_parallelCoordinator.cpp` | canonical-face barrier 与 GlobalDof barrier 的已有 workspace mapping trace 使用同一 step selector，并打印 rank、global block、owner rank、patch、Field/Flux/Residual 地址。 | MPI 调用、COPY/SUM 类型及其顺序不变。追踪仍位于 infrastructure，solver 层未引入 MPI。 |
| `tools/compare_high_order_steps.py` | 新增独立的 ASCII VTK/PVD 比较工具，比较 `rho`、`rhoU/V/W`、`rhoE`、`p` 的 count、L1、L2、L∞、相对差和 machine-roundoff scale。 | 不启动、不修改求解器。只有显式给出 `--abs-tol` 与 `--rel-tol` 才声明 `FIRST_DIVERGENCE_STEP`。 |

`SF_HIGH_ORDER_TRACE_STEP=N` 的 `N` 是 `StateBundle::step`，即物理步编号；输出文件的 ordinal `1` 是初始输出 `0` 后的第一个 commit，对应物理步 `0`，不能混为同一个编号。

## 2. 非侵入性验证

### 构建与架构

```text
ninja -C build sonicSolver
```

完整链接成功。Ninja 仍打印既有的 `premature end of file; recovering` warning，但没有编译或链接失败。

```text
python3 tools/check_architecture.py
```

检查通过：580 个源文件，`allowlisted dependency edges = 10`，没有新增 `core -> solver` edge。追踪代码没有改变 allowlist，也没有新增 solver 对 MPI implementation 的依赖。

### trace 开关不改变 serial 数值输出

以临时副本运行 `test/Sod/sodCase_weno7` 两个物理步：一份关闭 trace，另一份设置 `SF_HIGH_ORDER_TRACE=1 SF_HIGH_ORDER_TRACE_STEP=1`。该 case 的临时写出控制只改为逐步输出，求解配置未改。

```text
VTK_OUTPUTS_BYTE_IDENTICAL=yes
```

离线比较的三个输出 frame（0、1、2）中，六个比较量均为零差异。并且关闭/选择第 1 步的 run 不会输出 `step=0` trace。选中的第 1 步中，workspace-after-clear 对 FluxField 和 Residual 都报告 `l1=0`、`linf=0`、`non-zero-count=0`。

这证明 trace selector 和其只读检查没有改变该 serial case 的输出；它不构成 old/current 高阶路径等价性的证明。

## 3. 4-rank TENO5 Sod：首个输出差异

比较对象均来自同一临时 case 内容（4 rank、densityBase、singleFluid、Euler、TENO5、characteristic、StegerWarming、CFL 0.5、shared interface flux）；唯一临时 case 改动是 `writeInterval=1`，用于逐步输出。

```text
old binary:     ./build-phase1a-tests/sonicSolver
current binary: ./build/sonicSolver
comparison:     python3 tools/compare_high_order_steps.py OLD/result CURRENT/result \
                  --abs-tol 1e-12 --rel-tol 1e-10
```

| 输出 ordinal | `rho` L∞ | `rhoE` L∞ | `p` L∞ | 最大相对差 | 判定 |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0 | 0 | 0 | 0 | within tolerance |
| 1 | 3.0070896851028928e-02 | 3.2991502812525141e+02 | 1.2533557205884426e+02 | 1 | divergent |
| 2 | 4.8100329023607613e-02 | 6.3724397923249926e+03 | 2.5585549623650586e+03 | 1 | divergent |
| 110 | 1.9260366784397739e-01 | 7.3077156061762449e+04 | 2.3917750525128304e+04 | 1.996149076502381 | divergent |

比较工具的结果为：

```text
FIRST_DIFFERENT_STEP = 1
FIRST_DIVERGENCE_STEP = 1
SERIES_LENGTH_MISMATCH = old=501, current=111, compared=0..110
```

PVD 的 `timestep` 字段在这个 case 中是输出 ordinal，不是物理时间，因此 `FIRST_TIME_DIFFERENT_STEP = none` 仅表示 PVD ordinal 相同，不能用来说明物理时间相同。求解器日志给出的物理时间为：

| 已提交物理步 | old `time, dt` | current `time, dt` |
| ---: | --- | --- |
| 0 | `6.681531e-04, 6.681531e-04` | `6.681531e-04, 6.681531e-04` |
| 1 | `1.272174e-03, 6.040207e-04` | `1.225742e-03, 5.575894e-04` |

因此初始 state 与首个 dt 相同，差异在**物理步 0 的 assembly / commit 内**已经出现；它不是多步 roundoff 持续放大后的首个可见差异。之后 CFL 从不同 state 计算出不同 dt。

old binary 完成 500 步至 `time=2.435715e-01`；current binary 写出 `Sod_000110` 后，下一次 assembly 在 `splitStegerWarming` 因 `p=-872.514` fail fast。该 fail fast 未被绕过或替换为 fallback。

## 4. current 物理步 0 的检查点

以 host MPI 运行：

```text
SF_HIGH_ORDER_TRACE=1 SF_HIGH_ORDER_TRACE_STEP=0 \
mpirun -np 4 ./build/sonicSolver run --steps 1 <临时 TENO5 Sod case>
```

current run 的 composition 输出为：

```text
Convection dispatch : TENO5 + StegerWarming + HighOrderPerfectGas
```

四个 rank 的 local patch 都在 RHS 入口得到独立的 solver-owned workspace。每个 workspace 在原有 `equations.begin(...)` 之后满足：

| 项目 | 四个 local patch 的观测 |
| --- | --- |
| FluxField storage | `count=116280, l1=0, linf=0, non-zero-count=0` |
| Residual storage | `count=193800, l1=0, linf=0, non-zero-count=0` |
| canonical barrier mapping | 每个 local owner block 的 `matching-workspaces=1`；非本 rank block 为 0 |
| GlobalDof barrier mapping | 同样每个 local owner block 的 `matching-workspaces=1`；Residual 地址指向该 local PatchWorkspace |

追踪还记录了 current path 的候选面通量、canonical COPY 后面通量、方向 residual、GlobalDof SUM 前后 residual、RHS 和 final Q。这些数据已存在于本次 host run log；它们是为恢复 old source 后的同点比较准备的检查点，不能单独证明 old/current 的任一语义相同。

## 5. A/B/C 分类

| 假设 | 结论 | 证据与限制 |
| --- | --- | --- |
| A. Phase 4B workspace ownership / clear / patch association | 未裁决 | current 的四个 workspace 在 clear 后严格为零，且本地 block/workspace 映射一对一。这排除了当前 run 的“clear 后残留值”这一种表现，但旧二进制没有相同 trace，尚无法比较其 clear 时点、地址关联和每个 checkpoint。 |
| B. lifecycle / canonical COPY / GlobalDof SUM | 未裁决 | current trace 显示 canonical 与 GlobalDof barrier 对 local workspace 的一对一选择；不改变 COPY 或 SUM。没有 old 的相同 checkpoint，不能据此证明两个 revision 的 barrier 前后值和调用时点相同。 |
| C. high-order reconstruction / flux path | 候选，未证明 | 差异在首个物理步内出现，且 current 走 `TENO5 + StegerWarming + HighOrderPerfectGas`。这使高阶路径成为应继续比对的候选，但还没有 old 的 stage/RHS/face-flux checkpoint，不能宣称 first difference 已定位到 reconstruction、flux splitting 或任何具体公式。 |

结论：本轮证据排除“连续 roundoff amplification 才首次可见”的解释；尚不足以在 A/B/C 中做根因归属。没有把 workspace extraction 回滚，也没有据此修改任何数值算法。

## 6. 下一项唯一有效的隔离动作

在旧 binary 的实际 dispatch / thermodynamic contract 被重新确认前，不能取得
old source 就直接把它当作 real-high-order 对照。后续首先需要调查 old algorithm
composition / dispatch；只有该 contract 被确认后，才应取得对应 source revision 并在其上应用同一组只读 trace，以相同 rank、case、物理步 `0` 比较以下序列：

```text
initial Q
→ workspace-after-clear
→ boundary-ready Q
→ candidate face flux
→ canonical COPY 后 face flux
→ directional residual
→ GlobalDof SUM 前/后 residual
→ RHS
→ final Q
```

首个不一致 checkpoint 才决定下一步是只检查 A、只检查 B，还是进入 C 的高阶格式诊断。在该 checkpoint 出现前，Phase 4C、IO consistency cleanup、Field/IBM 重构和数值公式修复继续冻结。

## 7. 不变量复核

| 项目 | 结果 |
| --- | --- |
| CFL 公式 / dt 选择 | 未修改 |
| WENO/TENO、Steger-Warming、Lax-Friedrichs、Roe、Rusanov | 未修改 |
| stage coefficient、stage time、RK/Euler 生命周期 | 未修改 |
| physical boundary、halo、IBM 顺序 | 未修改 |
| canonical face `COPY` 语义 | 未修改 |
| GlobalDof `SUM` 语义 | 未修改 |
| `FluxField` / `Residual` ownership 与 layout | 未修改 |
| Field / geometry / IBM storage | 未修改 |
| hidden fallback / clamp | 未新增 |
