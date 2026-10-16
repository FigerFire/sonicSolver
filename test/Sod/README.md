# Sod 正式算例运行说明（2026-10-10）

三个算例均已完整运行至 0.2 s。算例内仍使用各自显式配置的空间格式，不存在运行时自动格式回退。

```sh
./build/sonicSolver run test/Sod/sodCase
./build/sonicSolver run test/Sod/sodCase_weno7
./build/sonicSolver run test/Sod/sodCase_weno7_t0p2
```

实际网格 x=[0,100]、y=[0,20]、z=[0,0.1]，100×20×1 单元，101×21×2 节点；不自动缩放长度单位。间断位于 x=50：左 rho=1、p=100000 Pa，右 rho=0.125、p=10000 Pa，两侧 U=0。主变量 rho/rhoU/rhoE，质量/动量/总能量方程，PerfectGas+hConst+const，gamma=1.4、R=287.05，mu=0、Pr=0.72，无黏性项，无 IBM/RAS/MRF 执行。保留的 phi/k/epsilon 配置没有新增这些输运方程；explain 中 momentum RHS=0，保留的 gravity 文件也没有形成源项。

边界文件及其顺序保持原样：U/p/rho 的 Right 和 sideWalls 为 zeroGradient，frontBack 为 empty；Left 未显式声明这三项，继承既有处理。未用新增 Left 项替换旧配置，因此冻结 WENO7 结果保持。

三个算例 startTime=0、endTime=0.2、无 endStep 截断，CFL=0.5，串行 nProcs=1/split=[1,1,1]，writeInitial=true，每 0.05 s 输出。原 sodCase 的 endTime=5/endStep=500 和按50步输出改为同样的 0.2 s 时间窗口；原 sodCase_weno7 的 endTime=5 改为0.2。

| 算例 | 时间 recipe | 空间格式 | 结果 |
|---|---|---|---|
| sodCase | SSPRK3 | TENO5 characteristic + StegerWarming | 0.2 s 完成 |
| sodCase_weno7 | forwardEuler | weno7Rusanov | 0.2 s 完成 |
| sodCase_weno7_t0p2 | forwardEuler | weno7Rusanov | 0.2 s 完成；冻结 dt/场哈希回归通过 |

原 sodCase 的 forwardEuler + TENO5/StegerWarming 在约0.076 s出现负压力。因此这里只显式修改该算例的时间 recipe 为 SSPRK3，保留空间离散、CFL 和原始初值，未修改源码算法，也未给 solver 增加 fallback。该算例的新结果不作为旧 forwardEuler 结果的等价回归。

热物性修改 `models/thermoDynamics.yaml`，时间/离散修改 `solvers/numerics.yaml`，终止与输出修改 `solvers/runtime.yaml`。旧 `models/equations.yaml` 已移入 `solvers/equations.yaml`，IBM 的 inactive 配置由 models 注册。完整原生选取在 solvers/state.yaml、solvers/equationRegistry.yaml、solvers/execution.yaml。

本轮完整输出保留在 `test/caseReadiness-20261010/full/sodCaseSSPRK3/result`、`full/sodCase_weno7/result`、`full/sodCase_weno7_t0p2/result`，原目录旧结果没有覆盖。
