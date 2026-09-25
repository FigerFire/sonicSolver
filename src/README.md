# src 源码导航

`src/` 只保留六个目录域。目录域用于归类，实际模块仍由各自的调度头文件和
CMake target 独立维护，不能因为物理位置相邻就相互读取内部实现。

```text
src/
├── app/             程序装配与 CLI/GUI 前端
├── solver/          四层求解架构
├── core/            公共状态、存储、配置和接口契约
├── methods/         可复用纯数学工具与数值方法
├── models/          物理、湍流、IBM 和初值模型
└── infrastructure/  IO、网格与并行通信
```

## 依赖方向

```text
app
  -> solver/algorithm
       -> solver/equation
            -> solver/discretization + solver/boundary
                 -> solver/linearAlgebra（隐式项需要时）

solver 按需调用 methods、models、core 和 infrastructure；
methods 不依赖求解流程，core 不拥有具体物理模型。
```

## 目录入口

- `app/application/SF_application.h`：case 装配和运行入口。
- `solver/algorithm/SF_solverAlgorithm.h`：第一层，求解流程。
- `solver/equation/SF_equation.h`：第二层，方程表达与组装。
- `solver/discretization/SF_discretization.h`：第三层，方程项离散调度。
- `solver/boundary/SF_boundary.h`：第三层，边界离散调度。
- `solver/linearAlgebra/SF_linearAlgebra.h`：第四层，线性系统求解。
- `methods/numerics/SF_numerics.h`：数值方法入口。
- `methods/math/SF_math.h`：纯数学方法入口。
- `core/config/SF_config.h`：强类型配置入口。

新增代码时，先根据职责选择目录域，再通过所属模块的调度入口公开；不要新增
跨目录域的大型 umbrella header，也不要把方程组装、数值核或 IO 混入同一个文件。
