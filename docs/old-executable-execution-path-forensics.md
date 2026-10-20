# 旧版 executable 数值执行路径取证

日期：2026-09-15  
目标 executable：`./build-phase1a-tests/sonicSolver`  
control executable：`./build/sonicSolver`

本阶段只检查既有 binary、旧 build artifacts，并用 LLDB 对临时 4-rank case 做函数命中取证。没有修改 production solver source、executable、IO schema、数值配置默认值或任何求解/MPI 语义。

## 1. 最终分类

```text
E1 — Rusanov bypass confirmed

OLD_EXECUTION_PATH = RusanovEOS::div
                     (EquationSetRusanov execution path, CONFIRMED)
OLD_HIGH_ORDER_DISPATCH = BYPASSED / NOT EXECUTED
OLD_SPLITTER_DISPATCH = NOT EXECUTED
HISTORICAL_TENO5_BASELINE = INVALID AS A HIGH-ORDER BASELINE
```

在 old Steger 与 old LF 两个 4-rank case 中，每个 rank 都从 active Euler RHS 进入 `SF::Numerics::RusanovEOS::div`。`SF::divDispatch`、`splitStegerWarming`、`splitLaxFriedrichs` 和 `TENO5::teno5_core` 的命中数全部为零。

因此，old/current physical step 0 差异首先是 numerical composition mismatch：old 执行 EOS-aware Rusanov，current 执行 `TENO5 + selected splitter + HighOrderPerfectGas`。该对照不能用于裁定 Phase 4B workspace、canonical COPY 或 GlobalDof SUM regression。

## 2. Binary 与 debug-symbol 状态

| 项目 | old | current |
| --- | --- | --- |
| `file` | Mach-O 64-bit executable arm64 | Mach-O 64-bit executable arm64 |
| 大小/时间 | 4.4 MiB，2026-09-11 10:42 | 4.4 MiB，2026-09-15 09:55 |
| Mach-O UUID | `EC19757A-6632-3572-869C-CFFDD2E07CD2` | `6EC098A0-D01E-391C-BE0B-4741525F1B9A` |
| stripped | No；`nm` 可见 external 与 local C++ symbols | No |
| executable 内嵌 `__DWARF` | 无有效 `.debug_info` 内容 | 无有效 `.debug_info` 内容 |
| `.dSYM` | 未发现 | 未发现 |
| object debug info | old build tree 保留 DWARF v5 object files | 本轮未依赖 |
| LLDB source/function resolution | 可用；由 executable debug map 与保留的 old object files 提供 | 可用 |

old 的 `SF_compressible.cpp.o` 与 `SF_rusanovEOS.cpp.o` 都是 arm64 Mach-O object，包含 DWARF v5。其 compile units 分别记录旧构建时的源路径和 `build-phase1a-tests` 工作目录。object 标记为 optimized，因此部分变量/栈帧不可用，但函数入口、调用栈和 breakpoint hit count 可可靠解析。

LLDB 显示的行号来自 old object 的 debug mapping；同一路径上的当前源文件内容可能已经变化。本报告只使用运行时函数名与 caller 关系，不把当前磁盘上的对应行文本当成 old source 证据。

## 3. old/current symbol inventory

| Symbol / keyword | old | current | Notes |
| --- | --- | --- | --- |
| `EquationSetRusanov` 字符串 | absent | present | old 中没有这个后期 enum 字符串；不能用字符串缺失否定 Rusanov 函数执行。 |
| `HighOrderPerfectGas` 字符串 | absent | present | current composition contract 可见；old 无该 contract 名称。 |
| `RusanovEOS::div` | `0x1001dbde8` | `0x1001e24d4` | 两者都链接；old runtime 已确认命中。 |
| `Compressible::System::convection` | `0x100072778` | `0x1000780c4` | 签名因 workspace extraction 不同。 |
| `divDispatch` | `0x1000727e4` | `0x1000781bc` | 两者都链接；old runtime 命中 0。 |
| `splitStegerWarming` | `0x1001dccf4` | `0x1001e340c` | 两者都链接；old Steger/LF runtime 均命中 0。 |
| `splitLaxFriedrichs` | `0x1001dd23c` | `0x1001e3954` | 两者都链接；old Steger/LF runtime 均命中 0。 |
| `TENO5::teno5_core` | `0x100090228` | `0x100095ff0` | 两者都链接；old Steger/LF runtime 均命中 0。 |

该 inventory 证明 high-order 与 Rusanov code 都被链接进 old executable，不能单独说明实际执行了哪一个。本报告的路径判定来自第 5～7 节的 runtime breakpoint evidence。

## 4. 旧 build artifact 证据

`build-phase1a-tests/` 仍保留：

```text
CMakeCache.txt
build.ninja
compile_commands.json
.ninja_log
old object files
old static libraries
```

关键 artifact 包括：

```text
src/solver/algorithm/.../SF_compressible.cpp.o
src/solver/algorithm/.../SF_multiPatch.cpp.o
src/solver/algorithm/.../SF_multiPatchTime.cpp.o
src/solver/equation/.../compressible/SF_compressible.cpp.o
src/methods/numerics/.../Flux/SF_rusanovEOS.cpp.o
src/methods/numerics/.../Flux/SF_sw.cpp.o
src/methods/numerics/.../Flux/SF_lxF.cpp.o
src/methods/numerics/.../convection/SF_TENO5.cpp.o
```

`build.ninja` 的 old link command 同时链接 `SF_equation`、`SF_flux`、`SF_teno5`、`SF_weno3/5/7`、`SF_solverAlgorithm`。这解释了 old binary 为什么同时拥有 Rusanov 与 high-order symbols；它仍然只是 code-presence evidence。

## 5. 4-rank LLDB 取证方法

所有 runtime hit 实验保持原始 4-rank execution。`mpirun` 启动四个独立 LLDB，每个 rank 写独立日志，避免日志互相覆盖。两组 case 都只运行 physical step 0：

```text
old Steger: TENO5 + characteristic + StegerWarming
old LF:     TENO5 + characteristic + LaxFriedrichs
```

每个 rank 对下列实际存在的 symbol 设置 auto-continue breakpoint：

```text
SF::Numerics::RusanovEOS::div
SF::divDispatch
SF::Riemann::splitStegerWarming
SF::Riemann::splitLaxFriedrichs
SF::TENO5::teno5_core
```

命中计数 run 均正常完成一步并以 process status 0 退出。另做一次 first-hit backtrace run，在 Rusanov 首次命中后主动结束 debugger；该 run 的 PRTE nonzero 是 debugger 在 MPI_Finalize 前终止造成的预期取证行为，不是 solver assertion 或数值失败。

## 6. old Steger runtime function hits

四个 rank 的结果完全一致：

| Breakpoint | rank 0 | rank 1 | rank 2 | rank 3 |
| --- | ---: | ---: | ---: | ---: |
| `RusanovEOS::div` aggregate | 2369 | 2369 | 2369 | 2369 |
| └ main `RusanovEOS::div` entry | 1 | 1 | 1 | 1 |
| └ internal cell lambda | 2368 | 2368 | 2368 | 2368 |
| `divDispatch` | 0 | 0 | 0 | 0 |
| `splitStegerWarming` | 0 | 0 | 0 | 0 |
| `splitLaxFriedrichs` | 0 | 0 | 0 | 0 |
| `TENO5::teno5_core` | 0 | 0 | 0 | 0 |

名义 Steger case 没有进入 Steger splitter，也没有进入任何 high-order dispatch/reconstruction 入口。

## 7. old LF runtime function hits

四个 rank 的结果同样完全一致，并与 old Steger 的命中矩阵相同：

| Breakpoint | rank 0 | rank 1 | rank 2 | rank 3 |
| --- | ---: | ---: | ---: | ---: |
| `RusanovEOS::div` aggregate | 2369 | 2369 | 2369 | 2369 |
| └ main `RusanovEOS::div` entry | 1 | 1 | 1 | 1 |
| └ internal cell lambda | 2368 | 2368 | 2368 | 2368 |
| `divDispatch` | 0 | 0 | 0 | 0 |
| `splitStegerWarming` | 0 | 0 | 0 | 0 |
| `splitLaxFriedrichs` | 0 | 0 | 0 | 0 |
| `TENO5::teno5_core` | 0 | 0 | 0 | 0 |

名义 LF 配置也没有进入 LF splitter。YAML splitter 在 startup log 中被读取和回显，但没有影响 active convection RHS。

## 8. Relevant backtraces

old Steger 和 old LF 的四个 rank 都得到同一条 active RHS caller 链。省略标准库转发帧后为：

```text
ddtDispatch(... Euler ...) RHS lambda
  ↓
CompressibleAlgorithm::step(Field&, double) RHS lambda
  ↓
SF::Numerics::RusanovEOS::div(Field&, EquationSet::Model const&)
```

代表性 old Steger rank 0 原始关键帧：

```text
frame #1  SF::Numerics::RusanovEOS::div(...)
          old SF_rusanovEOS.cpp:65
frame #2  CompressibleAlgorithm::step(Field&, double)::$_0::operator()(...)
          old SF_compressible.cpp:187
frame #9  ddtDispatch(...)::lambda(Field&)::operator()(...)
          old SF_time.h:89
```

old LF rank 0 的函数链、函数地址和 old debug-map 行号相同。其他三个 rank 也都在同一 Rusanov entry 停止。

这排除了“Rusanov 只用于 validation、construction 或未使用 precomputation”的解释：它由 Euler 的 active RHS callback 直接调用，并完成了整个一步运行。

## 9. 实际 equation/convection execution path

old binary 的 observable path 是：

```text
YAML parser
  ├─ reads convection = TENO5
  └─ reads flux = StegerWarming / LaxFriedrichs
          ↓
startup logger
  ├─ prints TENO5
  └─ prints selected splitter
          ↓
old CompressibleAlgorithm Euler RHS callback
          ↓
RusanovEOS::div(Field&, EquationSet::Model const&)
          ↓
high-order divDispatch / splitter / TENO5 are bypassed
```

因此 parsed `FluxSplitter` 在 old production density RHS construction 之前失去作用。可执行文件证据与高度可疑假设完全一致：配置被解析和打印，但 numerical algorithm 实际绑定到 EOS-aware Rusanov path。

由于 old binary 没有后期 `ConvectionThermodynamicContract` enum 字符串，本报告不声称 old source 中存在一个字面名为 `EquationSetRusanov` 的 enum 值；确认的是与该 contract 数值语义对应的 `RusanovEOS::div` 被 active RHS 执行。

## 10. First-step signature 交叉验证

| Comparison | rho L∞ | rhoE L∞ | p L∞ |
| --- | ---: | ---: | ---: |
| known single old-Rusanov vs current-high-order | 3.0070896868e-02 | 3.2991502581e+02 | 1.2533557114e+02 |
| 4-rank old nominal-TENO5 vs current-TENO5 | 3.0070896851028928e-02 | 3.2991502812525141e+02 | 1.2533557205884426e+02 |
| absolute signature difference | 1.6971073851790308e-11 | 2.3152514359026100e-06 | 9.1884426467458979e-07 |
| relative signature difference | 5.6436872921639083e-10 | 7.0177204386810490e-09 | 7.3310732905355730e-09 |

此前该接近性只能作为 corroborating evidence。现在 debugger 独立证明 old nominal TENO5 执行 Rusanov、且 high-order splitter 从未执行，两类证据已经闭环。

## 11. Historical baseline 资格

historical old 4-rank “TENO5 + StegerWarming” baseline 的稳定运行只能证明其实际 Rusanov composition 能运行到历史终点。它没有验证：

```text
TENO5 reconstruction
Steger-Warming splitting
HighOrderPerfectGas production dispatch
```

因此：

```text
HISTORICAL_TENO5_BASELINE = INVALID AS A HIGH-ORDER REGRESSION BASELINE
PHASE4B_OLD_CURRENT_COMPARISON = NOT ADMISSIBLE
```

不应恢复 old/current identical step-0 workspace trace，因为两者在 candidate flux 之前已执行不同的数学方法。

## 12. 下一阶段唯一允许的动作

调查路线改为：

```text
定位 first build/revision：
TENO5 + StegerWarming + HighOrderPerfectGas
在 production multi-patch RHS 中首次真实执行
          ↓
用 runtime function hit 确认：
divDispatch > 0
splitStegerWarming > 0
TENO5::teno5_core > 0
RusanovEOS::div = 0
          ↓
只有该版本才可候选为 current high-order regression baseline
```

在找到这个 baseline 前，继续冻结 Phase 4C、IO、Field、IBM、PatchWorkspace、lifecycle、canonical/GlobalDof 和高阶数值修复。

## 13. 安全与变更复核

| 项目 | 结果 |
| --- | --- |
| production source modified | No |
| old/current executable modified | No |
| PatchWorkspace / Field modified | No |
| CFL / WENO / TENO / SW / LF modified | No |
| fallback / limiter / clamp added | No |
| canonical COPY / GlobalDof SUM modified | No |
| build required | No；取证为 read-only，未触发 rebuild |
| numerical regression rerun | 仅运行已授权的一步 4-rank function-hit forensic；正常 hit-count runs 退出 0 |

