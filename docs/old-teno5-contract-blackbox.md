# 旧版 TENO5 Contract 黑盒验证与首差异前置隔离

日期：2026-09-15  
范围：只运行既有 old/current executable 与临时 case 副本。没有修改生产源码、PatchWorkspace、lifecycle、canonical COPY、GlobalDof SUM、CFL、TENO/WENO、Steger-Warming、physical-state validation 或 IO。

> **后续 executable 取证已闭环。** 4-rank LLDB 证明 old Steger 与 old LF
> 都由 Euler RHS 直接调用 `RusanovEOS::div`，而 `divDispatch`、两种 splitter
> 和 `TENO5::teno5_core` 均未执行。最终 E1 结论见
> [旧版 executable 数值执行路径取证](old-executable-execution-path-forensics.md)。

## 1. 要回答的问题与判别方法

先前的 old/current 4-rank Sod 比较显示：初始输出 ordinal 0 完全相同，ordinal 1 已显著分叉。因此差异已收缩到 physical step 0 的 assembly / commit 内。

在把这个差异归因于 workspace、lifecycle 或高阶格式以前，先验证 old executable 是否真的响应其名义的 `TENO5 + StegerWarming + HighOrderPerfectGas` contract。判别器必须是 executable 的一步输出，而不是 case 文件名、历史文档或 current composition log。

所有临时 case 都复制自 `test/Sod/sodCase`，并共同作了唯一的非数值变动 `writeInterval: 1`，以写出 ordinal 0 和 1。mesh、初始场、4 rank、Euler、CFL=0.5、boundary、geometry、runtime 及 `interfaceFlux: shared` 均相同。

| Case | 对流配置 | 结果 |
| --- | --- | --- |
| old A | `TENO5 + characteristic + StegerWarming` | 成功运行 1 个物理步。 |
| old B | `WENO7 + characteristic + StegerWarming` | 合法配置名，但 case 网格只有 3 层 ghost；旧 binary 正确 fail fast，因为 WENO7 需要至少 4 层。未修改 mesh 绕过。 |
| old C | explicit Rusanov / EquationSetRusanov | 当前 YAML 的 `ConvectionScheme` 入口只接受 WENO3/WENO5/TENO5/WENO7；没有合法 direct Rusanov contract 配置入口。未修改源码制造此 case。 |
| old D | `TENO5 + characteristic + LaxFriedrichs` | 成功运行 1 个物理步。 |
| current control A/D | 与 old A/D 完全相同的临时 case，使用 `build/sonicSolver` | 两个 case 均成功运行 1 个物理步。 |

## 2. old binary 自身的 composition 输出

命令：

```text
mpirun -np 4 ./build-phase1a-tests/sonicSolver run --steps 1 <temporary-case>
```

old A 实际打印：

```text
EquationSet       : singleFluid/perfectGas cold start
Formulation       : conservativeFluxDifference
Convection scheme : TENO5
Flux method       : StegerWarming
Interface flux    : sharedInterfaceFlux
```

old D 将 `Flux method` 打印为 `LaxFriedrichs`，其余同上。这说明旧 executable 会读取并回显 YAML 值。

但 old stdout 没有打印下列任何实际 dispatch / thermodynamic contract 信息：

```text
Convection dispatch
HighOrderPerfectGas
EquationSetRusanov
Rusanov
dispatch
```

因此：

```text
OLD_BINARY_CONTRACT_LOGGING = unavailable
```

日志只能说明它读取了配置，不能证明该配置进入了数值路径。

## 3. 黑盒输出比较

比较使用：

```text
python3 tools/compare_high_order_steps.py OLD_RESULT CURRENT_RESULT \
  --abs-tol 1e-12 --rel-tol 1e-10
```

工具对每一个 frame 的 `rho`、`rhoU`、`rhoV`、`rhoW`、`rhoE`、`p` 都比较 count、L1、L2、L∞ 与最大相对差。

### old A（TENO5/Steger）与 old D（TENO5/LF）

| 输出 ordinal | 六个字段的 L∞ | L1 / L2 | VTK SHA-256 | 判定 |
| ---: | --- | --- | --- | --- |
| 0 | 全部 0 | 全部 0 | 完全相同 | identical |
| 1 | 全部 0 | 全部 0 | `Sod.pvd`、`.vtm`、4 个 `.vts` 均完全相同 | identical |

比较工具输出：

```text
FIRST_DIFFERENT_STEP = none
FIRST_DIVERGENCE_STEP = none
```

这不是只比较最大误差：old A 与 old D 的所有结果文件逐文件 SHA-256 清单无差异。

### current control A（TENO5/Steger）与 current control D（TENO5/LF）

current binary 明确打印：

```text
Convection dispatch : TENO5 + StegerWarming + HighOrderPerfectGas
Convection dispatch : TENO5 + LaxFriedrichs + HighOrderPerfectGas
```

ordinal 1 的黑盒差异为：

| 字段 | L∞ | L1 | L2 | 最大相对差 |
| --- | ---: | ---: | ---: | ---: |
| rho | 3.1720471991655175e-02 | 4.1871023028984702 | 3.6444047706763544e-01 | 1.3439478483439987e-01 |
| rhoU | 2.1316282072803006e-14 | 1.5681948032242490e-12 | 1.3787144962267104e-13 | 5.0000000000000000e-01 |
| rhoV | 7.5960213596438312e-16 | 6.0768170877150651e-14 | 6.7940880474815989e-15 | 5.0000000000000000e-01 |
| rhoW | 0 | 0 | 0 | 0 |
| rhoE | 2.9103830456733704e-11 | 1.4551915228366852e-09 | 1.8521564287013979e-10 | 2.7222782452671327e-16 |
| p | 1.1893572720719021e+02 | 1.0774658074786916e+04 | 1.1166837661237225e+03 | 5.7691180195443975e-03 |

四个 ordinal-1 `.vts` 文件的 SHA-256 都不同。由此可知，在相同 Sod 首步、相同 YAML 入口和相同 4-rank runtime 中，这个 splitter 替换是能够被真实 high-order contract 观测到的有效 discriminator；old 的零差异不能归咎于 Sod 首步天然无法区分 Steger 与 LF。

## 4. WENO7 与 explicit Rusanov 的处理

old WENO7 case 在 mesh validation 阶段 fail fast：

```text
mesh:canonicalMetric declares nGhost=3, but convection=WENO7 and ILW=0
require at least 4 ghost layers.
```

这符合实现能力检查。没有修改 mesh、ILW 或 solver 来使 WENO7 强行运行。

Rusanov 不是现有 `convection.default` 可选 scheme，当前合法 YAML 入口也没有单独的 `ConvectionThermodynamicContract` 键。因此没有制造一个伪造的 explicit-Rusanov case。

不需要为 reconstruction 再构造 smooth case：本轮更强的 flux-splitting discriminator 已经在 old binary 中完全失效、而在 current 真实 high-order contract 中明确有效。按本阶段的停止规则，应先暂停 A/B/C 与旧源码 trace 恢复，检查 old algorithm composition / dispatch；平滑场不能使旧 TENO5/Steger 基线重新获得已证实的 high-order 资格。

## 5. First-step signature comparison

先前已确认的 single nominal WENO7 比较中，old 是 `EquationSetRusanov`、current 是 `HighOrderPerfectGas`。将该记录的首步 signature 与本轮 4-rank old-TENO5/current-TENO5 signature 并列：

| Comparison | rho L∞ | rhoE L∞ | p L∞ |
| --- | ---: | ---: | ---: |
| previous single old-Rusanov vs current-high-order | 3.0070896868e-02 | 3.2991502581e+02 | 1.2533557114e+02 |
| 4-rank old nominal TENO5 vs current TENO5 | 3.0070896851028928e-02 | 3.2991502812525141e+02 | 1.2533557205884426e+02 |
| absolute signature difference | 1.6971073851790308e-11 | 2.3152514359026100e-06 | 9.1884426467458979e-07 |
| relative signature difference | 5.6436872921639083e-10 | 7.0177204386810490e-09 | 7.3310732905355730e-09 |

这组接近性是强线索：4-rank old nominal TENO5 的首步输出很可能属于与先前 old single Rusanov 相同的数值 contract。它不是单独的证明；本报告的 contract 判定依据是 old splitter 不敏感、且 current real-high-order control 对同一 splitter 替换敏感。

## 6. 最终判定与停止点

```text
CONTRACT CONTRADICTED

The historical 4-rank TENO5 baseline did not demonstrably exercise
the assumed high-order TENO5/Steger-Warming contract.
```

更精确地说，old executable 回显 `TENO5` 与 `StegerWarming`，但其 ordinal-1 数值输出对 `StegerWarming -> LaxFriedrichs` 完全不敏感；该替换在 current 明确的 `HighOrderPerfectGas` contract 中产生可观测差异。因此 historical 4-rank TENO5 baseline 不再具备 Phase 4B high-order numerical regression baseline 的资格。

后续 LLDB function-hit 取证已经证明 old 的 active Euler RHS 直接调用
`RusanovEOS::div`，且 high-order dispatch、两种 splitter 与 TENO5 reconstruction
均未命中。因此内部执行路径已经从黑盒推断提升为 E1 runtime-confirmed。

后续顺序必须是：

```text
old algorithm composition / dispatch investigation
→ old contract 的可观测确认
→ 再决定是否恢复 old source 并加同一套 step-0 trace
```

在此前继续冻结 Phase 4C、IO cleanup、Field/IBM refactor、PatchWorkspace ownership、lifecycle、canonical COPY、GlobalDof SUM、CFL 与高阶数值修改。

## 7. 本阶段未改变的语义

| 项目 | 结果 |
| --- | --- |
| 生产 solver 源码 | 未修改 |
| temporary case 中的数值配置 | 仅用于独立 black-box 对照；未写回仓库 case |
| physical boundary / halo / IBM order | 未修改 |
| canonical COPY / GlobalDof SUM | 未修改 |
| CFL、WENO/TENO、splitter 算法、positivity handling | 未修改 |
| hidden fallback / limiter / clamp | 未新增 |
| architecture checker / build | 本轮无生产代码变更，未重新构建 |
