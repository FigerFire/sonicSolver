# variationalIBMCase

该算例验证密度基 FDM 上的第一阶段 variational/DLM 体约束 IBM。计算域、字段和
刚体运动均使用 SI；圆柱 STL 是 `t=0` 参考构型。

```text
densityBased explicit predictor
  -> volume constraint C u = U_s
  -> local KKT elimination
  -> lambda_b, body force/torque, mechanical-work energy update
```

当前设置为均匀 `20 m/s` 空气流与静止刚体，推进 50 步、每 5 步输出一次。
背景网格采用与 `IBMCase` 相同的约 `50` 点/物体直径分辨率，但保持单 block，
因为当前体约束乘子的多 patch/MPI 所有权尚未实现。正确运行应报告非零
`constrainedCells`、满足所设约束容差，并在 `result/` 输出
`IBMConstraintMask`、`IBMMultiplierX/Y/Z`，可直接观察圆柱内部约束和近场响应。

`geometryFiles cylinder.STL;` 使用统一路径规则，自动读取
`constant/triSurface/cylinder.STL`。显式地址必须加引号，例如
`geometryFiles "constant/triSurface/cylinder.STL";`。

运行：

```sh
./build/sonicSolver test/variationalIBMCase
```

该算例不验证 surface multiplier、压力—乘子单体 KKT、MPI 或两向 FSI；这些组合
当前均应明确失败。
