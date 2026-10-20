# Surface 一阶矩／力矩回归

与 canonical case 使用相同冻结 mesh、STL、flow numerics 和 endTime=5，仅显式选择 `linearReproducing`。MPI variants 使用 y 切分，使 surface 支撑跨越真实分区界面。`mesh.sfm` 为输入 fixture；VTK preview / result 为生成文件。

```sh
./build/sonicSolver run test/IBM/cylinderFlowPeskinLinear
```

这是数值方法变化；旧 `partitionOfUnity` baseline 继续保留。新方法要求完整 3D affine 支撑，不满秩会明确失败；权重可以有符号，不能 clipping 或静默退回旧归一化。
