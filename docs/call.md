# 单次方程调用的编译边界

2026-11-19。以 R4 交付时的工作树为起点。本轮调整“怎么解方程”的编译分工，涉及 HOW、WHICH 与空间编译接线。物理、时间与空间实现继续使用原数值代码。

## 生产入口

原 `Lowering::compileProgram` 一并索引树、匹配 STATE/WHICH、调方法并校验返回契约。随后空间编译用 `equationMethod` 再查直接方法，另以方法名识别 Rhie–Chow。当前生产链的职责如下。

| 阶段 | 输入与决定 | 交付给下一阶段 |
| --- | --- | --- |
| [SF_lowering.cpp](../src/solver/system/SF_lowering.cpp) 的 `compileProgram / compileCalls` | 规范化 HOW，确定叶节点身份与作用域。建立方程／目标引用，核对空、重复和无主绑定。由 [SF_stateRealization.cpp](../src/solver/system/SF_stateRealization.cpp) 的 `resolveCallInputs` 将 STATE 读版本绑定到对应调用 | 一个已定位的叶节点。块叶同时携带原有的方程集和未知量集 |
| [SF_numericalCompiler.cpp](../src/solver/system/SF_numericalCompiler.cpp) 的 `compileCall` | [SF_builtinNumerics.cpp](../src/solver/system/SF_builtinNumerics.cpp) 的 `resolveNumerics` 依通配、方程和步骤覆盖合成 WHICH。每份选择先规范化算子 occurrence，再合成并核对地址。方法编译检查符号、显式 STATE 输入、所选时间能力和操作契约 | 已有的 `CompiledEquationCall`，包括实际操作或片段、读写、状态效果和数学数据。需要后续空间编译的实现交出函数入口 |
| `compileProgram` 与已有空间编译 | 合并目标与附加视图、资源、融合成员、残差和发布。展开 recipe 前释放原索引，展开后按 occurrence 顺序核对并重建父作用域引用。此后 `NumericalCompiler::compile` 用完整程序所确定的 recipe 和算子目录，调用逐次选定的空间入口 | 在原时点完成空间绑定，并由 `SolvePlanner::compile` 冻结计划。Explain 和运行绑定消费该结果 |

同一 `E` 的两个 HOW 步骤各获得一个 `EquationCall` 和 `CompiledEquationCall`。步骤引用保留独立的读版本。按方程绑定的 WHICH 可以同时应用于两次调用，按步骤绑定仍能覆盖各自选择。有效算子以方程内 occurrence 寻址，空间结果额外标记所属调用 occurrence。没有按方程共享可变工作区。

## 三条生产链

- **Cell 与压力响应。** [SF_equationContribution.cpp](../src/solver/system/SF_equationContribution.cpp) 的 `finalizeMathematics` 按控制层引用替换 `E_U/viscU`，并保留 `appliedForce` 和外层权重。单次编译解析 `viscU` 的 WHICH、`SST.muEff: Current`、预测目标与工作视图。 [SF_cell.cpp](../src/solver/system/SF_cell.cpp) 的 `compileCellEquation` 处理调用契约，已选的 `lowerCellEquation` 在空间阶段绑定 `stress.central2.arithmetic-faces`、压力梯度和源项。`bindCell` 实现实际输入视图，`CellEquationExecution::evaluate` 消费编译项，由 [SF_cellEquation.h](../src/solver/discretization/SF_cellEquation.h) 计算定向面牵引并装至相邻单元。预测实例拥有候选值、Working 输出和 `dt/rho` 响应。`bindSingleFluidOperations` 按 predictor occurrence 将这一实例交给 [FaceFluxExecution](../src/solver/algorithm/SF_faceBinding.cpp)，压力沿该面资源完成校正。
- **普通标量与 Stage SST。** [SF_scalarMethod.cpp](../src/solver/system/SF_scalarMethod.cpp) 的 `compileScalarEquation` 消费目标、所选标量形式、源项、边界和真实通量依赖，生成本调用数据。LocalAffine 的 `compileLocalAffine` 与 SSTCentral2 的 `compileSstClosure` 同经直接方法入口，分别提供局部表达式和非局部闭合。Stage SST 的黏度闭合按 STATE 读取该阶段状态并发布，容量输运用自己的标量调用和共享质量通量。HOW 保持闭合、关系和各 RHS 的相对位置。共同 recipe 和 `StageGroup` 仍管入口快照、所有 RHS 就绪屏障、统一推进及最后 Current 刷新。
- **线性块。** HOW 的 `BlockSolve` 将原有方程集和未知量集交给 [SF_block.cpp](../src/solver/system/SF_block.cpp) 的 `LinearBlockMethod::compileBlock`，必须有步骤级 WHICH，例如 `LinearBlockCentral2Hypre`。第一个方程和目标仍提供兼容的引用锚点。方法保留方程行顺序、未知量布局、原恒定 Central2 扩散和 HYPRE 选择，并生成单个 assemble／solve／publish 片段，不产生物理时钟操作。

## 删减与 R5 接口

- 删除 `lowerDirectSpatial`、第二次直接方法查找、空间阶段的 Rhie–Chow 名称分支和 provider 目录参数。直接方法映射只在本地编译时选择。Rhie–Chow 随已有面契约交出空间入口，`ConservativeResidual` 与旧守恒压力实现共用既有融合 lowering。遗留 `IProvider` 的空间虚调用及重复 `usesSpatialRecipes` 改为交出可空的静态函数入口；入口无需借用 provider 对象。
- `compileDirectEquation` 的空结果表示该入口无此方法，然后按同一 WHICH 名称查仍在用的目录；所选实现的异常附带调用上下文并传播。保留普通标量形式查找、压力各关系、块接口及 Eulerian、IBM/KKT、once-step RAS 等真实专用入口。它们仍有自身布局、生命周期或完整数学消费条件。
- 删除整树编译中的第二次算子 AST 地址遍历、调用指针到块指针的映射，以及展开后重建的步骤指针列表。合成后的 occurrence 仍作唯一性检查。原始索引限于 `compileCalls`，融合使用展开后建立的父索引。
- R5 可复用每次调用的 `providerContract`、操作、读写效果、输出视图、字段需求、残差和响应引用。运行绑定仍使用私有契约类型选择 binder，按 predictor occurrence 接线后再绑定读取，各方法仍有自己的端口解析和回调。旧流体、模型和约束操作所有者仍从冻结计划查询操作需求。这些属于运行接线，空间方法已在编译中选定。

## 所有权、顺序与证据

公开四轴语法、内置方法名和操作 ID 沿用原定义。接口变化限于内部单次编译和空间接线，四个引用内部接口的测试源码已同步。物理由原数组所有者持有，数值工作区依 occurrence 建立，时钟仍由物理步终端提交。predictor 基值和标量 OldTime 的保存顺序、显式 RHS 屏障、面牵引符号与 canonical COPY／残差 SUM 均沿原实现。没有新增 fallback，实际数值等价仍待对照。

相对本轮起点有 8 个生产文件变更，新增 245 行、删除 203 行，净增 42 行，无新增或删除生产文件。统计排除文档、测试源码和已有脏树。详见 [逐文件增删](../test/call/evidence/delta.json) 和 [源码差异](../test/call/evidence/source.diff)。

`DEVELOPER_DIR=/Library/Developer/CommandLineTools cmake --build build --target sonicSolver -j 4` 生产编译和链接成功，退出码 0。日志记录已有的枚举分支与 JSON 指针弃用警告。 [构建日志](../test/call/evidence/build.log) · [产物标识](../test/call/evidence/finish.json)

源码核对覆盖同方程多步骤、`E_U/viscU`、标量和 Stage SST，以及线性块的输入和顺序。本轮未运行算例、测试、回归、diagnostics、check/explain 或架构脚本。待指定算例验证预测与校正通量、变系数牵引预算、固定时间基值、阶段时刻和最终发布，以及并行 COPY/SUM 参考。
