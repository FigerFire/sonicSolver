# 真实 High-Order Production Baseline 历史构建定位

日期：2026-09-15

## 1. 调查边界与判定方法

本次调查只读取现存 executable、build artifact、符号表与运行时输出，没有修改 production source、数值格式、workspace ownership、MPI 语义或 case 的 mesh/初始场。

搜索范围限定为当前 repository 及其现存 build tree。共找到三个真实 `sonicSolver` executable。`build/._sonicSolver`、`build-phase1a-tests/._sonicSolver` 和 `build-tests/._sonicSolver` 是 4096-byte AppleDouble 元数据文件，不是 Mach-O executable，未列入候选。

运行时使用相同的 4-rank、density-based、Euler、TENO5、characteristic、CFL 0.5、shared-interface-flux、one-step Sod discriminator：

- SW case：`TENO5 + StegerWarming`
- LF case：`TENO5 + LaxFriedrichs`

LLDB 分别对以下精确函数设置 one-shot、auto-continue breakpoint：

```text
SF::Numerics::RusanovEOS::div(
SF::divDispatch(
SF::Riemann::splitStegerWarming(
SF::Riemann::splitLaxFriedrichs(
SF::TENO5::teno5_core(
```

one-shot breakpoint 命中后会从最终 `breakpoint list` 消失；未命中的 breakpoint 保留并显示 `hit count = 0`。四个 rank 都以 status 0 退出。这里匹配的是 `RusanovEOS::div`，不会把 CFL 计算使用的 `RusanovEOS::deltaT` 计为 Rusanov convection。

## 2. HISTORICAL_BINARY_INVENTORY

三个 executable 的 SHA-256、Mach-O UUID 均不同，因此不存在需要按内容合并的重复 binary。

| Build | executable mtime | Size (bytes) | SHA-256 | Mach-O UUID | File type | Stripped |
| --- | --- | ---: | --- | --- | --- | --- |
| `build-tests/sonicSolver` | 2026-09-10 14:41:43 +0800 | 4,120,040 | `3e13d4e629d931d1d514ce164b9c99330f87ebdc1295185f41a69dd7ab2355bb` | `023DDA80-1CF8-354B-85A2-11DA9B91A172` | Mach-O 64-bit arm64 executable | No；`nm` 可见 external、non-external 与 cold symbols |
| `build-phase1a-tests/sonicSolver` | 2026-09-11 10:42:59 +0800 | 4,585,400 | `423abc82179880ed180fa8f2005db057025e0c046f565f5966ba5ee99f2fbfaa` | `EC19757A-6632-3572-869C-CFFDD2E07CD2` | Mach-O 64-bit arm64 executable | No；`nm` 可见 external、non-external 与 cold symbols |
| `build/sonicSolver` | 2026-09-15 09:55:06 +0800 | 4,665,208 | `d22b67411700b5a8c992a419d8f1bdfce1148da615e1a272e6476cc9c7a1da39` | `6EC098A0-D01E-391C-BE0B-4741525F1B9A` | Mach-O 64-bit arm64 executable | No；`nm` 可见 external、non-external 与 cold symbols |

mtime 只用于描述 artifact chronology，不用于推断执行路径。下面的分类以 ABI/build artifact 和运行时函数命中为依据。

## 3. Build chronology 与静态 composition-era 分类

| Build | 相关 artifact | Constructor ABI | Density lifecycle / workspace ABI | Static era |
| --- | --- | --- | --- | --- |
| `build-tests` | `SF_compressible.cpp.o` 2026-09-10 14:15:01；`SF_multiPatch.cpp.o` 11:18:55；`SF_singleFluid.cpp.o` 14:41:42 | `CompressibleAlgorithm(SolverConfig)` | 无 `stepDensity`、`DensityBasedRHS::assembleAllPatches`、`PatchWorkspace` 或 thermodynamic-contract symbol；`System::convection(Field&, AssemblyContext const&)` | S0 — pre-contract composition |
| `build-phase1a-tests` | `SF_compressible.cpp.o` 2026-09-11 10:42:36；`SF_multiPatch.cpp.o` 10:42:35；`SF_singleFluid.cpp.o` 10:42:57 | `CompressibleAlgorithm(SolverConfig)` | 无 `stepDensity`、`DensityBasedRHS::assembleAllPatches`、`PatchWorkspace` 或 thermodynamic-contract symbol；`System::convection(Field&, AssemblyContext const&)` | S0 — pre-contract composition |
| `build` | `SF_densityBasedRHS.cpp.o` 2026-09-15 09:54:52；`SF_compressible.cpp.o` 09:54:53；`SF_multiPatch.cpp.o` 09:55:05；`SF_singleFluid.cpp.o` 09:55:06 | `CompressibleAlgorithm(SolverConfig, ConvectionThermodynamicContract)` | 有 `stepDensity`、`DensityBasedRHS::assembleAllPatches(... vector<PatchWorkspace>& ..., ConvectionThermodynamicContract)`；`System::convection(Field&, FluxField&, Residual&, AssemblyContext const&)` | S2 — explicit contract composition present |

`strings` 给出同样的分界：两个 S0 binary 都没有 `EquationSetRusanov`、`HighOrderPerfectGas`、`PatchWorkspace` 或 `Convection dispatch` 字符串；当前 S2 binary 包含这些可见 composition/diagnostic 字符串。三个 binary 都链接了 `RusanovEOS::div`、`divDispatch`、两个 splitter 和 `TENO5::teno5_core`，所以仅凭这些数值 symbol 的存在无法判定 active path。

历史 build artifact 没有显示 S1 中间态：现存序列从两个 S0 executable 直接跳到一个 S2 executable。这个事实只说明现存 artifact 有缺口，不证明中间 build 从未存在。

## 4. Runtime function-hit matrix

下表中的 `>0` 表示四个 rank 各自均至少命中一次，`0` 表示四个 rank 各自最终 breakpoint hit count 均为零。`build-phase1a-tests` 的非 one-shot 取证还记录到每 rank `RusanovEOS::div` aggregate hit count 2369（入口 1 次、内部 cell lambda 2368 次）；其余四个目标函数均为 0。

| Build / case | `RusanovEOS::div` | `divDispatch` | SW | LF | TENO5 | Rank exits | Runtime result |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `build-tests` / SW | >0 | 0 | 0 | 0 | 0 | 4/4 status 0 | Rusanov density convection |
| `build-tests` / LF | >0 | 0 | 0 | 0 | 0 | 4/4 status 0 | Rusanov density convection |
| `build-phase1a-tests` / SW | >0 | 0 | 0 | 0 | 0 | 4/4 status 0 | Rusanov density convection |
| `build-phase1a-tests` / LF | >0 | 0 | 0 | 0 | 0 | 4/4 status 0 | Rusanov density convection |
| `build` / SW | 0 | >0 | >0 | 0 | >0 | 4/4 status 0 | High-order TENO5/SW |
| `build` / LF | 0 | >0 | 0 | >0 | >0 | 4/4 status 0 | High-order TENO5/LF |

临时运行时证据位于：

```text
/private/tmp/arch-build-tests-sw-rank-{0,1,2,3}.log
/private/tmp/arch-build-tests-lf-rank-{0,1,2,3}.log
/private/tmp/old-sw-lldb-rank-{0,1,2,3}.log
/private/tmp/old-lf-lldb-rank-{0,1,2,3}.log
/private/tmp/arch-current-sw-rank-{0,1,2,3}.log
/private/tmp/arch-current-lf-rank-{0,1,2,3}.log
```

## 5. SW/LF numerical discriminator

### S0 builds

`build-tests` 的 SW 与 LF ordinal-1 四个 block 文件逐 rank SHA-256 完全相同。六个输出字段 `rho`、`rhoU`、`rhoV`、`rhoW`、`rhoE`、`p` 的 L1、L2、L∞ 差异全部为零。`build-phase1a-tests` 的既有黑盒结果相同。因此配置中的 splitter 替换没有进入两个 S0 executable 的 active convection。

### 当前 S2 build

当前 build 的 SW 与 LF ordinal-1 四个 block 文件逐 rank SHA-256 全部不同。关键差异为：

| Field | L1 | L2 | L∞ | 最大相对差 |
| --- | ---: | ---: | ---: | ---: |
| `rho` | 4.18710230289847019e+00 | 3.64440477067635438e-01 | 3.17204719916551747e-02 | 1.34394784834399872e-01 |
| `rhoE` | 1.45519152283668518e-09 | 1.85215642870139785e-10 | 2.91038304567337036e-11 | 2.72227824526713271e-16 |
| `p` | 1.07746580747869157e+04 | 1.11668376612372253e+03 | 1.18935727207190212e+02 | 5.76911801954439751e-03 |

`rhoU` 与 `rhoV` 只有 roundoff-scale 差异，`rhoW` 相同；`rho` 与 `p` 已给出明确、可观察的 splitter sensitivity。它与严格函数命中矩阵共同满足 H2 判据。

## 6. Historical binary classification

| Build | mtime | Static era | Rusanov div | divDispatch | SW | LF | TENO5 | Numerical SW/LF sensitivity | Classification |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `build-tests/sonicSolver` | 2026-09-10 14:41:43 +0800 | S0 | >0 in SW/LF | 0 | 0 | 0 | 0 | ordinal 1 identical | R0 — confirmed Rusanov production path |
| `build-phase1a-tests/sonicSolver` | 2026-09-11 10:42:59 +0800 | S0 | >0 in SW/LF | 0 | 0 | 0 | 0 | ordinal 1 identical | R0 — confirmed Rusanov production path |
| `build/sonicSolver` | 2026-09-15 09:55:06 +0800 | S2 | 0 | SW/LF >0 | SW case >0 | LF case >0 | SW/LF >0 | ordinal 1 differs | H2 — confirmed real HighOrderPerfectGas production path |

## 7. FIRST_H2_BUILD 与 workspace ownership

现存 build 中第一个、也是唯一通过 runtime contract 的 H2 是：

```text
build/sonicSolver
```

它的 workspace model 是 **W1 — PatchWorkspace-owned flux/residual**。这一判断不依赖目录名或 mtime，binary ABI 已直接保留以下证据：

```text
DensityBasedRHS::assembleAllPatches(
    ...,
    std::vector<PatchWorkspace>&,
    ...,
    ConvectionThermodynamicContract)

System::convection(
    Field&,
    FluxField&,
    Residual&,
    AssemblyContext const&)
```

当前对应源码也显示 `PatchWorkspace` 直接拥有一个 `FluxField convectiveFlux` 和一个 `Residual residual`，并由 `DensityBasedRHS` 按 participating patch 接收 workspace collection。这与 binary 的 W1 ABI 一致。

两个更早的 S0 binary 使用 `System::convection(Field&, AssemblyContext const&)`，没有 `PatchWorkspace` symbol，属于 W0-era evidence；但它们的 runtime classification 都是 R0，不是 H2。因此现存 artifact 中没有 `H2 + W0` executable。

## 8. Phase 4B 与 high-order activation 的相对关系

从现存、可运行 artifact 能证明的顺序是：

```text
S0 + W0 + R0
    build-tests
    build-phase1a-tests

        [没有现存 executable 覆盖这个 transition window]

S2 + W1 + H2
    build
```

因此不能把 high-order activation 与 Phase 4B workspace extraction 分离为历史 A/B 对照。现存的稳定历史输出证明的是 Rusanov convection 的稳定性，不是 TENO5/Steger-Warming 的稳定性。它既不能证明 Phase 4B 造成当前 high-order instability，也不能证明 Phase 4B 与其无关。

“旧版稳定 TENO5”这一历史 baseline 前提已被 runtime evidence 否定；更窄的“Phase 4B 是否改变真实 high-order 数值”假设则因为缺少 `H2 + W0` baseline 而无法用现存 executable 检验。

如果未来恢复旧 source 并人为接入 `HighOrderPerfectGas`，该结果只能称为 controlled reconstruction experiment，不能称为 historical production baseline。

## 9. Rusanov 语义区分

当前 H2 density production convection 不调用 `RusanovEOS::div`，但 timestep estimation 仍调用 `RusanovEOS::deltaT` 计算 wave-speed-based CFL。前者是本次分类的 convection path，后者是独立的 timestep calculation；本报告没有修改或否定后者。

## 10. 结论与下一步边界

```text
FIRST_CONFIRMED_HIGH_ORDER_BUILD = build/sonicSolver
FIRST_CONFIRMED_HIGH_ORDER_TIME = 2026-09-15 09:55:06 +0800 (executable artifact mtime)
FIRST_CONFIRMED_HIGH_ORDER_WORKSPACE_MODEL = W1 — PatchWorkspace-owned FluxField/Residual

PRE_PHASE4B_HIGH_ORDER_BASELINE = NO

PHASE4B_HIGH_ORDER_REGRESSION_HYPOTHESIS =
    NOT_TESTABLE_WITH_EXISTING_BASELINES

NEXT_VALID_INVESTIGATION =
    HIGH-ORDER NUMERICAL VALIDATION：独立验证 TENO5 reconstruction、
    Steger-Warming splitting、LF control、smooth convergence、shock-tube behavior、
    CFL stability envelope、positivity robustness、multi-patch/single-patch equivalence。
```

本次只定位 historical baseline，没有开始上述 numerical validation，也没有重新开启 Phase 4B architecture regression isolation。
