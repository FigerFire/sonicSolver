# Fractional DLM 跨分区支撑回归（MPI-2）

复制 canonical `cylinderFlowFictitiousDomain` 的 source mesh、STL、数值与物理配置；仅沿 y 切分 `[1,2,1]`，使 body constraints 跨越 y=10 分区界面。保留 endTime=5、writeInterval=0.1。

```sh
mpirun -np 2 ./build/sonicSolver run test/IBM/cylinderFlowFictitiousDomainCrossingMPI2
```

138 个 canonical owner body points 在两个 rank 上按 78/60 分布，物理总体积 6.9；空 body rank 合法。与串行 reference 比较物理体积加权 Q、constraint residual 和 force/power；不能仅按普通 x 切分 case 通过就宣称 body 支撑跨分区已经覆盖。
