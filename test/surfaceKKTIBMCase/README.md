# surfaceKKTIBMCase

该算例验证压力基表面变分 IBM 的单体
`(p', deltaU, Lambda_s, q)` KKT。默认 `motionMode motivation`，刚体质量、
主转动惯量、外力和外力矩全部从 `constant/IBMProperties` 读取。

表面使用 STL 三角形重心积分，`J` 使用显式半径的 Wendland C2 紧支撑核并做
partition-of-unity 归一化；传播严格采用带表面积权的离散伴随 `S=J^T`。

运行：

```sh
./build/sonicSolver test/surfaceKKTIBMCase
```

将 `motionMode` 改为 `rotate`、`externalForce` 置零并设置
`externalTorque`，即可验证转动惯量方程。
