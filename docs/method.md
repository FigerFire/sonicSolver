# R4 方程方法的数学入口

2026-11-19。起点是 kernel 轮交付后的脏工作树，HEAD 为 `b2866848360e4fcb9b39c0efdf46ce69bd918c78`。本轮属于“怎么解方程”的编译、资源描述与运行绑定重构。方程和项仍来自控制组合，已有数值核继续执行。

## 从三个调用读源码

三个调用都由 native 输入或内置预制展开成四轴。[SF_lowering.cpp](../src/solver/system/SF_lowering.cpp) 的 `compileProgram` 解析实际 HOW occurrence、目标和 STATE 输入。[SF_builtinNumerics.cpp](../src/solver/system/SF_builtinNumerics.cpp) 的 `resolveNumerics` 决定有效 WHICH，[SF_builtinProviders.cpp](../src/solver/system/SF_builtinProviders.cpp) 的 `compileDirectEquation` 调用对应数学入口。预制 HOW 和用户显式 HOW 都沿这条链路。

| HOW 调用 | 声明、方法编译与输入版本 | 实际执行、所有者和读者 |
| --- | --- | --- |
| `sst.effectiveViscosity: SST.muEff -> SST.muEff` | [SF_turbulenceSystemContribution.cpp](../src/models/turbulence/SF_turbulenceSystemContribution.cpp) 声明层流黏度加所选涡黏度的关系及 `Current` 读取。WHICH 选 LocalAffine。[SF_affine.cpp](../src/solver/system/SF_affine.cpp) 的 `compileLocalAffine` 从等式生成仿射系数和余项，声明发布与空间需求 | [SF_scalarBinding.cpp](../src/solver/algorithm/scalar/SF_scalarBinding.cpp) 的 `bindAffine` 建立每次调用的局部求值实例；`AffineOperations::solve/write` 写入物化工作区、做 NearestInterior 延拓并发布。Cell 动量读同一 `SST.muEff: Current` |
| `sst.viscosity: mu_t -> mu_t` | SST 声明明确列出密度、速度或守恒速度关系、k/omega 和参数。WHICH 选 SSTCentral2。[SF_sstMethod.cpp](../src/solver/system/SF_sstMethod.cpp) 的 `compileSstClosure` 绑定操作数和输出角色。fixedVolume Step 读取声明的 Current 输入；守恒 Stage 沿各次声明读取 Stage 视图，终端接受后从 Current 重算 | `bindSst` 连接版次读取和 `SstClosureOperations::refresh`。[SF_sstClosure.h](../src/solver/discretization/SF_sstClosure.h) 提供 `sstEddyViscosity`。模型数组由 `Turbulence::Manager` 拥有，由 [SF_singleFluid.cpp](../src/app/application/execution/SF_singleFluid.cpp) 按编译发布和实际切片提供。非局部闭合、两条 `SST.nu.*` 和 `sst.effectiveViscosity` 读取该次发布；守恒黏性仍走原 Manager 端口 |
| `predictor: E_U -> U`，或固定外迭代的 Working U | [SF_model.cpp](../src/app/application/model/SF_model.cpp) 解析 `momentumStress`，[SF_equationContribution.cpp](../src/solver/system/SF_equationContribution.cpp) 的 `finalizeMathematics` 替换 `E_U/viscU`。WHICH 选 CellExplicit、算术面应力和源项实现。[SF_cell.cpp](../src/solver/system/SF_cell.cpp) 的 `compileCellEquation / lowerCellEquation` 同处，解析压力、速度、`SST.muEff: Current` 和额外源的真实读写 | [SF_cellBinding.cpp](../src/solver/algorithm/scalar/SF_cellBinding.cpp) 的 `bindCell` 建立每 occurrence 的 `CellEquationExecution`。`evaluate/advance` 装配 RHS 和发布预测；物理 U 仍由 STATE 的字段所有者提供。候选值、Working 和压力响应属于该预测实例。`FaceFluxExecution` 按 `predictor` 借用响应 |

默认 PISO 的 `momentum/momentum.diffusion` 与 [bind.md](bind.md) 的显式 `E_U/viscU` 都经过同一组合。显式输入实际保留 `appliedForce`，以 `source.field` 读取 `aForce: Current` 并按原权重进入 RHS。应力 occurrence 保持原位置、符号、外层权重和周围项。WHICH 的 `stress.central2.arithmetic-faces` 经 [SF_fusedSpatial.cpp](../src/solver/system/SF_fusedSpatial.cpp) 的目录绑定，最终由 [SF_cellEquation.h](../src/solver/discretization/SF_cellEquation.h) 的 `stressArithmeticCentral2/stressFace` 使用实际系数视图，每张定向面计算一次牵引，并以相反符号装配相邻单元。随后 `CellEquationExecution` 发布预测速度及原 `dt/rho` 响应。[SF_faceBinding.cpp](../src/solver/algorithm/SF_faceBinding.cpp) 的 `FaceFluxExecution::predict` 借该实例，[SF_pressureOperators.cpp](../src/solver/algorithm/pressure/SF_pressureOperators.cpp) 的 `correctFlux` 发布原校正面数组 `phi`，固定时间 k/omega 继续读取该版本。PISO 的 Physical 与 SIMPLE/PIMPLE 的 Working 预测各沿已有流程。这里的证据来自源码追踪。

## 方法入口迁移

| 起点位置或对象 | 交付位置与职责 |
| --- | --- |
| `SF_scalarMethod.cpp::compileLocalAffine`、`SF_methodObjects.h` 的前置声明、标量头中的私有契约 | [SF_affine.cpp](../src/solver/system/SF_affine.cpp) 与 [SF_affine.h](../src/solver/system/SF_affine.h)。等式分解、版本能力和本地发布要求同处；原公共前置声明已删 |
| `SF_scalarMethod.cpp::compileCellEquation / lowerCellEquation`、标量头中的 Cell 契约 | [SF_cell.cpp](../src/solver/system/SF_cell.cpp) 与 [SF_cell.h](../src/solver/system/SF_cell.h)。资源收集和逐项 lowering 共用应力结构判定。删除压力响应的重复 gradient 指针，核验原单项梯度及已有势变量、容量条件 |
| `SF_scalarMethod.cpp::LinearBlockMethod / addScalarMethod` | [SF_block.cpp](../src/solver/system/SF_block.cpp) 与 [SF_block.h](../src/solver/system/SF_block.h)。`linearBlockMethod()` 在原内置目录注册，真实块装配、线性求解与发布操作继续沿原接口；旧注册包装已删 |
| `SF_scalarMethod.cpp::validateScalarProgram` | [SF_executionContract.cpp](../src/solver/system/SF_executionContract.cpp) 的 `validateScalarExecution`。保留原冻结位置，读取实际发布者、别名写入、区间、通量生产和消费以及唯一提交 |
| `fieldResource`、两个 `equivalentMath`、标量方法头中的 `ScalarOps` | `equationFieldResource` 与 STATE 资源实现相邻，保留原字段请求和 owner。标量/SST 共用 [SF_formula.cpp](../src/core/system/SF_formula.cpp) 的 `sameFormulaExpression` 及原规范化编码。操作 ID 集中在 [SF_operationIds.h](../src/core/system/SF_operationIds.h)，原标识字符串保持一致 |
| `bindSingleFluidOperations` 中二次 Cell 类型比较、Stage 标量数统计和两份区间列表 | 一次方法选择后消费 `BoundEquation` 中的实际借用端口和区间回调。两类区间动作进入一份有序列表；二次类型判断与旧分组已删 |

`SF_scalarMethod.cpp` 从 873 行缩至 414 行，保留普通标量编译中共享的表达式、RHS、边界资源和时间能力处理。节点与 metric 输运的几何和边界、容量及共享通量、固定时间缺陷求值仍由各真实实现处理。线性块具有不同的未知量集合和装配生命周期。SST 非局部核需要梯度及壁距，继续保留实际操作。已有执行端口仍在原 binder，避免为文件切分增加转发链。带权源项、重复 occurrence 及其累加顺序沿用原代码。系统 CMake 接入三份新实现；控制层、说明和 binder 使用各数学契约的头部。公开 case 输入、WHICH 名称和操作 ID 沿用原接口。

## 绑定、时间与所有权

- **单实例端口。** 原统计通过“Stage 参与者且类型非 Cell”推断标量单实例端口的使用。现在 `bindScalar` 仅在借用 `SolverServices::scalarEquation` 时交回该端口。公共循环核验该实际端口恰有一个借用者，并继续拒绝同时提供实例表。普通实例表仍检查键、资源匹配和生命周期。
- **共同区间。** 原 Cell 回调先准备速度与边界、提出步长，再处理纯步长提议的固定时间标量。现在按是否声明 `prepare` 在启动时稳定排序，组内保留 HOW 顺序。仍逐调用完成准备及提议，再设置共同 dt，保存 predictor 的固定基值和标量 OldTime，随后发布 OldTime。区间接受沿同一顺序。改变的是启动回调表达，dt 提议及快照的算术和调用位置沿用原实现。
- **共享 Stage。** 实例化时挂接输出，Step 的面资源也在此时取得指定 predictor 实例。随后解析输入，再安装操作。Stage 中全部参与者完成入口快照和 Ready 后才计算 RHS，全部 RHS 完成后统一推进。最终物理接受后重算守恒 SST 的 Current 黏度与后续关系，最后执行唯一 TimeCommit。Stage 模型数组和原 once-step RAS 扫描仍保有各自数学及频率。
- **静态与动态检查。** `validateScalarExecution` 是原编译收尾的跨调用检查，未增加运行检查。单流体操作表封存路径、独立执行入口和可重绑 Eulerian 的绑定检查保持原设计。执行期仍由消费点校验 revision、边界及区间。公开 OldTime、物理字段、通量、方法工作区和 clock 的 owner 保持原分工。没有改变串行牵引装配或现役 MPI 的 canonical COPY、残差 SUM。

## 三种维护动作

| 目标 | 维护者要改的声明 |
| --- | --- |
| 给仿射关系加一项 | 在 [bind.md](bind.md) 的显式 `E_mu` RHS 上加一个动力黏度常量，例如 `1.0e-6`、量纲 `[1,-1,-1,0,0,0,0]`。用字段代替常量时，在该 occurrence 声明受支持的读取。LocalAffine 编译现有加法，EddyMu 和应力消费者沿原发布 |
| 调整合法顺序 | 在 [bind.md](bind.md) 的 fixedVolume authored HOW 中交换 `sst.sink.omega` 和 `sst.effectiveViscosity`。两次调用写不同结果，无互相读取，仍位于所需发布之后、predictor 之前。压力流程及 k/omega Frozen 系数生命周期保持。对于阶段预制，应按其真实调度依赖核对顺序 |
| 换已支持的方法 | 对常系数纯扩散且满足两者共同 STATE 和时间要求的普通标量，把 WHICH 从 `ScalarDiffusionCentral2` 改为 `ExplicitCentral2Arithmetic`，选择既有算术面实现。若配置显式算子绑定，应一并使用相容选择。对于 fixedVolume 黏度也可按 [bind.md](bind.md) 选择 LocalAffine 及与之相容的 `E_mu` 公式，这会改变物理闭合 |

上述操作只改相关声明，沿同一 per-call 编译和绑定。新数值算子应先补数学签名与实际实现，再修改所属方法及选择入口。R4 是本轮主体，R1/R2 只支撑源码导航和入口核对，R3/R5/R7/R8 限于本次触及的依赖。R6 的 once-step RAS 原位扫描及其他固定核仍按现役数值差异保留。若替换原位 RAS，应保留并比较 k/omega 或 k/epsilon 更新顺序、限制器、边界、分项 RHS 和黏度发布。涉及融合或并行应力的替代应补齐边界就绪状态、dt、阶段时刻、定向面预算与 canonical COPY/SUM 对照。

## 构建与证据

按本轮源码快照计算，修改 25 个现有生产文件、新增 6 个生产文件，删除整文件 0 个。新增 637 行、删除 549 行，净增 88 行，包含迁移行和 CMake，排除文档、测试源码及起点脏树。新增文件为三组 `SF_affine`、`SF_cell`、`SF_block` 的实现与私有契约。逐文件明细和基线记录在 [delta.json](../test/method/evidence/delta.json)、[source.diff](../test/method/evidence/source.diff) 和 [start-status.txt](../test/method/evidence/start-status.txt)。两个现有测试源码同步采用仿射契约所在头部。

```sh
DEVELOPER_DIR=/Library/Developer/CommandLineTools cmake --build build --target sonicSolver -j 4
```

生产目标编译与链接完成，退出码 0。构建仍报告既有枚举 switch 分支和 JSON pointer 弃用警告。实际命令、日志与产物标识见 [build.log](../test/method/evidence/build.log) 和 [finish.json](../test/method/evidence/finish.json)。

本轮只做源码追踪、差异审阅和生产构建。未运行测试、算例、回归、诊断、check/explain 或架构脚本，未数值验收。待数值对照包括默认与显式 PISO 首步 RHS、牵引预算、预测与校正通量，SIMPLE/PIMPLE 的固定基值及 OldTime，守恒 Stage 的最终 Current 重算与唯一提交。既有 MPI 参考及历史未通过项保留。
