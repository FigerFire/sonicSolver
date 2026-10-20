# IBM 正式算例运行说明（2026-10-10）

本轮整理本目录直接包含 `case.yaml` 的 27 个正式算例。历史诊断与冻结目录保持原状。24 个已支持算例完成实际 10 步启动验证（12 个串行、12 个 MPI）；另外 3 个 KKT 算例仍是 Unsupported，输入检查通过不代表它们能够执行。

## 直接运行

在仓库根目录：

```sh
./build/sonicSolver run test/IBM/cylinderFlowGhost
./build/sonicSolver run test/IBM/cylinderFlowViscousGhost
./build/sonicSolver run test/IBM/cylinderFlowSSTGhost
mpirun -np 2 ./build/sonicSolver run test/IBM/cylinderFlowPeskinLinearMPI2
mpirun -np 4 ./build/sonicSolver run test/IBM/cylinderFlowVelocityForcingLinearMPI4
```

其它 MPI 算例按下表 rank 数启动；普通串行算例直接 `run`。数值结果写入各算例的 `result/`。本轮验证在 `test/caseReadiness-20261010/` 副本中执行，没有清理或覆盖正式算例的旧结果。若手动运行希望清除旧帧，请自行使用 `sonicSolver cleanResult <case>`。

## 配置目录

所有 solver 配置位于 `solvers/`，文件头统一为 `object: solver`。STATE 选择为 `solvers/state.yaml`，方程选择为 `solvers/equationRegistry.yaml`，执行选择为 `solvers/execution.yaml`；它们由 `solvers/solvers.yaml` 显式注册。热物性及 IBM 仍位于 `models/`，使用 `object: models`。

## 共同设置

- WHAT/STATE：质量、动量、总能量守恒方程，主变量 `rho/rhoU/rhoE`；只有两个 SST 算例另外推进 `k/omega` 并计算 `mu_t`。其余保留的 k/epsilon 初值不代表启用 RAS 输运。
- 物性：`models/thermoDynamics.yaml` 显式选择 PerfectGas+hConst+const，gamma=1.4，R=287.05 J/(kg K)，动力黏度 mu=1.7894e-5 Pa s，Pr=0.72。空间配置不再保管 mu/Pr。
- 网格不变：实际 `mesh.sfm` 为 x=[0,100]、y=[0,20]、z=[0,0.1]，100×20×1 单元，101×21×2 节点。不是扩大后的 100×40×1；STL 没有移动。按 SI 物性解释时坐标为 m，求解器不自动把 m 换成 mm。
- 圆柱截面名义中心 (20,10)，半径 5；原始圆柱 STL 轴向范围 [-2.5,2.5]，穿过薄层域。surface 方法使用各自 cylinderSurface.stl；以对应 IBM.geometryFiles 为准。
- 初始 U=(2.36643,0,0) m/s，rho=1.225 kg/m³，p=101325 Pa。普通算例 k=0.005、epsilon=0.0005，但 DNS 配置不执行湍流 PDE。SST 初始 k=1e-6、omega=1，epsilon 保留但不参与 SST 输运。
- 外边界保持原配置：Left 的 U/rho/p 为上述 fixedValue；Right 与 sideWalls 均为 zeroGradient；frontBack 为 empty。这里没有改成右侧固定压力或上下 symmetry。SST 的 k/omega 入流分别 fixedValue 1e-6/1，其余出口/侧边 zeroGradient。
- 时间与空间：classicalRK4，WENO5 characteristic + LaxFriedrichs，shared canonical flux，CENTRAL2 黏性/导热。CFL=0.3。普通 forcing 类 maxDeltaT=0.01 s；三个有界黏性/SST Ghost 案例为 1e-5 s。
- Ghost 的圆柱壁面为 stationaryNoSlipAdiabatic：静止无滑移、绝热，Euler ILW 关闭（value=0），使用既有黏性 Ghost 实现。forcing/DLM 用自己方法的壁面速度目标，不套用 Ghost wallClosure。
- forcing/DLM 中的 ILW=3 是已有外部物理边界处理，仍保留；不能把它解释成 Ghost Euler/slip 圆柱壁面。prescribed 刚体目标速度为零；selfPropelled 算例则匹配随耦合更新的刚体速度。BP 是有限 penalty 逼近，不是精确乘子约束。

## 算例、方法、rank 和终止条件

| 算例 | IBM 数值方法 | rank / split | 配置 endTime / maxDeltaT | 本轮验证 |
|---|---|---|---|---|
| cylinderFlowDFMAugmentedLagrangian | dfmAugmentedLagrangian / partitionOfUnity | 1 / [1, 1, 1] | 5.0 / 0.01 | Unsupported KKT，10步 t=未推进 |
| cylinderFlowDFMExplicitSelfPropelled | dfmExplicitSelfPropelled / 不适用 | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392368e-03 |
| cylinderFlowDFMFractionalStepSelfPropelled | dfmFractionalStepSelfPropelled / 不适用 | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392217e-03 |
| cylinderFlowDFMImplicitSelfPropelled | dfmImplicitSelfPropelled / partitionOfUnity | 1 / [1, 1, 1] | 5.0 / 0.01 | Unsupported KKT，10步 t=未推进 |
| cylinderFlowFictitiousDomain | dfmFractionalStepPrescribed / 不适用 | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392217e-03 |
| cylinderFlowFictitiousDomainCrossingMPI2 | dfmFractionalStepPrescribed / 不适用 | 2 / [1, 2, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392217e-03 |
| cylinderFlowFictitiousDomainCrossingMPI4 | dfmFractionalStepPrescribed / 不适用 | 4 / [1, 4, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392217e-03 |
| cylinderFlowFictitiousDomainMPI2 | dfmFractionalStepPrescribed / 不适用 | 2 / [2, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392217e-03 |
| cylinderFlowFictitiousDomainMPI4 | dfmFractionalStepPrescribed / 不适用 | 4 / [4, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392217e-03 |
| cylinderFlowFullyImplicitDLM | dfmImplicitPrescribed / partitionOfUnity | 1 / [1, 1, 1] | 5.0 / 0.01 | Unsupported KKT，10步 t=未推进 |
| cylinderFlowGhost | stationaryNoSlipAdiabatic / 不适用 | 1 / [1, 1, 1] | 1.0 / 0.01 | 启动通过，10步 t=4.392446e-03 |
| cylinderFlowPeskin | peskinOriginal / partitionOfUnity | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392520e-03 |
| cylinderFlowPeskinLinear | peskinOriginal / linearReproducing | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392502e-03 |
| cylinderFlowPeskinLinearMPI2 | peskinOriginal / linearReproducing | 2 / [1, 2, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392502e-03 |
| cylinderFlowPeskinLinearMPI4 | peskinOriginal / linearReproducing | 4 / [1, 4, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392502e-03 |
| cylinderFlowPeskinMPI2 | peskinOriginal / partitionOfUnity | 2 / [2, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392520e-03 |
| cylinderFlowPeskinMPI4 | peskinOriginal / partitionOfUnity | 4 / [4, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392520e-03 |
| cylinderFlowSSTGhost | stationaryNoSlipAdiabatic / 不适用 | 1 / [1, 1, 1] | 0.02 / 0.00001 | 启动通过，10步 t=1.000000e-04 |
| cylinderFlowSSTGhostAuthored | stationaryNoSlipAdiabatic / 不适用 | 1 / [1, 1, 1] | 0.02 / 0.00001 | 启动通过，10步 t=1.000000e-04 |
| cylinderFlowVelocityForcing | velocityForcingFTS / partitionOfUnity | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392429e-03 |
| cylinderFlowVelocityForcingBP | velocityForcingBP / 不适用 | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392367e-03 |
| cylinderFlowVelocityForcingLinear | velocityForcingFTS / linearReproducing | 1 / [1, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392398e-03 |
| cylinderFlowVelocityForcingLinearMPI2 | velocityForcingFTS / linearReproducing | 2 / [1, 2, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392398e-03 |
| cylinderFlowVelocityForcingLinearMPI4 | velocityForcingFTS / linearReproducing | 4 / [1, 4, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392398e-03 |
| cylinderFlowVelocityForcingMPI2 | velocityForcingFTS / partitionOfUnity | 2 / [2, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392429e-03 |
| cylinderFlowVelocityForcingMPI4 | velocityForcingFTS / partitionOfUnity | 4 / [4, 1, 1] | 5.0 / 0.01 | 启动通过，10步 t=4.392429e-03 |
| cylinderFlowViscousGhost | stationaryNoSlipAdiabatic / 不适用 | 1 / [1, 1, 1] | 0.02 / 0.00001 | 启动通过，10步 t=1.000000e-04 |
