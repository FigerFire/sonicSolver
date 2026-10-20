# Phase 4B：从 Field 分离 Residual / Face-Flux Workspace

## 1. 完整调用矩阵

| 写入者 | 读取者 | 显式 storage |
|---|---|---|
| 对流重构、Rusanov、低阶通量 | canonical-face COPY、残差装配 | `FluxField&` |
| 对流/粘性/源项 | RK/Euler、GlobalDof SUM | `Residual&` |
| 算法 | 全部 workspace consumer | `PatchWorkspace` 按 patch 索引 |

## 2. Ownership Before

`Field` 同时拥有守恒量、几何、边界/IBM metadata、`convectiveFlux_` 和
`residual_`；所有离散、时间积分和 MPI 都通过 Field accessor 读写后两项。

## 3. Chosen Workspace Owner

`CompressibleAlgorithm::workspaces_` 是唯一 owner。它与
`StateBundle::patches` 保持一对一索引，`PatchWorkspace` 只含已有的
`FluxField` 与 `Residual`，没有指向 Field 的缓存指针。

## 4. Why Field / StateBundle / Runtime Do Not Own It

Field 只保留 physical state、几何、边界和 IBM metadata；StateBundle 只组合
physical state；Runtime 在 canonical COPY 和 GlobalDof SUM 屏障临时借用显式
传入的 storage，既不持有也不注册为第二个 owner。

## 5. Representation / Allocation / Resize / Clear

`PatchWorkspace::ensureFor(field)` 仅在 `MX/MY/MZ/NVar` 或 face 数改变时
重分配。每个 RHS stage 由 `Equation::Compressible::System::begin` 按原顺序
清空 flux 和 residual；RK4 原有 stage 间残差清空保留为对 `Residual` 的清空。
不存在每 step 临时分配，也不复制 `Field` 或 `Q`。

## 6. Field API Removed

删除 `Field` 的 `FluxField`/`Residual` 成员、setup/resize allocation、
`clearResidual`、`clearConvectiveFlux`、`clearSource`、`ConvectiveFlux`、
`ResX/Y/Z`、`Source`、local/global residual accessor。保留 canonical face
几何索引；它不是数值 flux workspace。

## 7. RHS / Discretization / Sources

`DensityBasedRHS` 对每个 patch 显式传入 workspace。WENO/TENO、Rusanov、
低阶 flux、粘性、重力、MRF、WallHeat、界面张力、相变耦合均改为接收
`FluxField&` 或 `Residual&`。`Field` 没有 solver include，core→solver 依赖为零。

## 8. Pressure / IBM

压力分支的 legacy `ddtDispatch` 也改为显式消费 algorithm workspace residual。
IBM ghost、forcing、KKT 和 boundary 调用顺序未改；本阶段没有把 IBM workspace
加入 `PatchWorkspace`。

## 9. Runtime / MPI Semantics

canonical-face 阶段仍在所有 patch 产生候选 F* 后执行一次：Runtime 只借用
`Field* + FluxField* + Residual*`，MPI COPY 后重装配结构 residual。GlobalDof
仍是 local contributor → SUM → owner → broadcast，storage 改为借用 residual。
没有平均 face flux，没有改变 COPY/SUM 语义。

## 10. Numerical Lifecycle

Euler、SSPRK3、RK4/RK4 stage 节点、CFL、boundary/halo、IBM、源项、pressure
corrector 顺序均未改变。此变更只改变 workspace 的 owner 和参数传递；预期不改变
数值结果。

## 11. Quantitative Result

| 指标 | Before | After |
|---|---:|---:|
| Field workspace members | 2 | 0 |
| Field workspace public API families | 5 | 0 |
| workspace authoritative owners | Field per patch | Algorithm vector |
| canonical barrier locations | 1 | 1 |
| GlobalDof SUM barrier locations | 1 | 1 |
| core → solver includes | 0 | 0 |
| Field direct semantic data members | 16 | 14 |
| Field public callable declarations | 约 118 | 约 100 |
| Field writable API families | 51 | 46 |
| Field responsibility categories | 6 | 5 |
| Field-owned FluxField | 1 | 0 |
| Field-owned Residual | 1 | 0 |
| solver-owned FluxField per participating patch | 0 | 1 |
| solver-owned Residual per participating patch | 0 | 1 |
| additional numerical storage copies | 0 | 0 |

## 12. Validation

架构检查通过：allowlisted dependency edges 为 10，新增 core → solver 为 0。
静态搜索确认 Field 不再拥有 `FluxField`、`Residual`，且不存在旧 forwarding
accessor。按 object 的增量编译已覆盖 field、residual、flux、equation、density
RHS/time、application coupling、execution runtime 和 MPI workspace 接口。

本 closeout 环境中，`cmake --build build --parallel 4` 在重新生成 Ninja metadata
后仍显示 `ninja: premature end of file; recovering`，并由执行通道在约 30 秒后中断，
没有报告 compile failure 或 link failure。因此不能将完整 build、完整 CTest、host
MPI CTest、4-rank Sod 或冻结数值回归标为通过。状态为 **implementation complete,
validation pending**。

## 13. Deferred Issues

pressure-based time ownership、scalarTransport ownership、multiphase temperature
advance、model derivative / SF_viscous、Field remaining responsibilities、IBM
architecture、distributed KKT、HYPRE ownership 和 Equation DSL/AssemblyPlan 均未处理。

## 14. Closeout Self Review

| 项目 | 结果 |
|---|---|
| Full build passed | 未确认：host build 被执行通道中断 |
| Full CTest including host MPI passed | 未确认 |
| 4-rank Sod regression passed | 未确认 |
| Field owns FluxField / Residual | No / No |
| Second FluxField / Residual exists | No / No |
| Per-stage allocation introduced | No |
| Core → solver dependency introduced | No |
| Canonical COPY / GlobalDof SUM semantics changed | No / No |
| Residual sign / source clear timing changed | No / No |
| Geometry / IBM storage modified | No |
