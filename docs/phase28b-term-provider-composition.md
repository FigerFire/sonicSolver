# Phase 28B：结构收敛与 Equation Term Provider 组合

状态用语：**Implemented** 表示代码和相应回归已运行；**Interface-only** 表示有注册/编译入口，但没有完整生产执行入口；**Legacy** 表示仍由已验证的既有数值路径执行；**Unsupported** 表示明确拒绝执行。

## 1. Scope

本轮保持 `Equation Composition → Transformation → Numerical Compilation → SolvePlan → Provider Resolution → Runtime` 主干，不修改 PISO/SIMPLE/PIMPLE 控制流、压力公式、WENO/TENO、IBM 几何或 MPI COPY/SUM 语义。新增 Momentum 物理项应通过 equation term 和 term provider 进入现有 pressure predictor。

## 2. Frozen baseline

开始时 HEAD 为 `537889b05234c34217aa1a8eff7ab7844b69b4b9`，工作树已有上一阶段的未提交修改；可恢复快照位于 `/private/tmp/sonic-phase28b-baseline-20260927`。修改前完整 build 成功，架构 allowlist 为 2 条 `methods → solver` 已知边，host CTest 24/24；sandbox 内 3 个 MPI 用例因 PRTE socket 无法启动，host 重跑均通过。基线 solver 与 pressure VTS 的 SHA 已存于该快照。

## 3. Removed dead/legacy files

**Implemented**：删除无 production include、实例或 CMake 编译入口的 `src/app/application/model/legacy/` 五个文件；清理 `src/` 中的 AppleDouble `._*`，`.gitignore` 已忽略 `._*` 和 `.DS_Store`。没有删除其 Git 历史。

## 4. SolveStrategyKind cleanup

**Implemented**：`SolveBlock` 带 typed `ExecutionPolicyKind`，EulerianStepper 的执行分支消费它；删除无人使用的 `SolveStrategyDescriptor` 和配置 accessor。`SolveStrategyKind` 仍以 **Legacy** metadata 留在 composition、validation、explain；Planner 对旧 metadata 与 typed policy 冲突时抛错，不静默选择一方。

## 5. Equation coupling ownership cleanup

**Implemented**：与 MultiphaseModel/LevelSet/PhaseChange 直接耦合的 concrete implementation 从 `solver/equation/coupling/` 移到 `models/physics/equationRuntime/`；移除 generic equation umbrella 中的 concrete Eulerian 引用。该移动不改变任何数学调用。

## 6. Eulerian adapter ownership

**Implemented**：既有 `SF_eulerianEquationAdapter` 源码迁到 `solver/algorithm/eulerian/equation/` 并保持专用 target/公式。它仍是 **Legacy** specialized execution，不声称已泛化为 equation-term runtime。

## 7. Integer identifier migration

**Implemented**：pressure GlobalDof entity/backend row storage 使用 `int64_t`；runtime/coordinator/halo/MPI backend 增加窄整数 identifier COPY 与整数 collective。测试含大于 `2^53` 的 identity，禁止经过 `double`。单 rank 多 pressure block、插值式 identifier halo 仍 **Unsupported**，会 fail-fast；已有 2/4-rank pressure 回归通过。

## 8. Built-in OpId centralization

**Implemented**：`core/system/SF_operationIds.h` 集中 pressure、Momentum、velocity/flux correction、IBM constraint、flow/explicit time 与 Eulerian 内置 OpId，plan/operation binding 使用常量。开放的自定义 OpId 仍是字符串；provider ID、capability 名称与 trace label 不混入 OpId 类型。

## 9. EquationRole

**Implemented**：`EquationId` 保持开放身份，`EquationRole` 表达 Mass/Momentum/Energy 等数学角色。内置 composition 明确赋 role；gravity/MRF 的 target 查 role 和 solved unknown，不从 `E_MOMENTUM*` 前缀推断。自定义 ID 只要显式声明 Momentum role，也可匹配相同 source provider。

## 10. SystemContribution before/after

**Implemented**：从保存 `RawEquationSystem&` 的对象改为可拷贝的纯值贡献；contributor 显式接收 raw 输入，result 仅保存方程/term/constraint/closure/policy 值。测试让 raw 离开作用域后继续复制并检查 contribution。`BuildRequest::userContributions` 为 C++ typed 贡献入口，不建立新 context/manager。

## 11. CompiledTermBinding

**Implemented**：扩展既有 `BoundTerm`，保留原 equation ID、实际 transformed execution equation ID、EquationRole、左右侧、ordinal、recipe、provider 与窄 kernel callback。它与纯数学 `Equation::Term` 分离；原 Momentum term 仅映射一次到 `E_MOMENTUM_PREDICTOR`。自定义 primitive source 可以由 provider 定义执行而不冒充一个内置 `SourceKind` recipe；保守格式没有 recipe 的 source 会在编译时 fail-fast。

## 12. TermProviderCatalog

**Implemented**：独立于 operation provider catalog。描述符拥有 ID、typed match、可选 recipe 与窄数值 callback；NumericalCompiler 解析并冻结结果。内置注册包括 conservative flux、primitive Upwind1、Central2/4、conservative/primitive gravity、MRF 和 energy wallHeat。C++ 调用者可从 built-in catalog 增加 descriptor，并通过 BuildRequest 在 compilation 时注入。

## 13. Deterministic resolution

**Implemented**：0 匹配为 Unsupported、1 匹配为 Resolved、多匹配为 Invalid；冲突 ID 排序后报告，与注册次序无关。term 的浮点累加顺序按 equation declaration ordinal，而不取决于 catalog 容器迭代。测试覆盖空匹配、错误 role、缺少 U unknown、错误 recipe 和双注册歧义。

## 14. Transient term provider

**Legacy**：时间推进仍由 `CompiledTimeRecipe` 和现有 time operations 执行；pressure fixed-time predictor 仍是 `U_n + dt R(U_k)`。`ddt(U)` 还不是普通 `BoundTerm` 的 kernel callback，不把现有 fixed-time 数学伪装为通用 transient provider。BackwardEuler/DIRK/IMEX 继续 **Unsupported**。

## 15. Convection provider

**Implemented**：primitive Upwind1 term 解析为 `convection.primitiveUpwind1`，其已冻结 kernel 在 pressure Momentum assembly 中执行；`oldFlux` workspace 仍由 pressure operation 准备。原 face flux、upwind 选择与求和顺序保持。保守 WENO/TENO flux 仍为 **Legacy** fused kernel，仅由绑定选择 recipe，不重写。

## 16. Diffusion provider

**Implemented**：Central2 term 解析为 `diffusion.central2`，pressure predictor 调用其已冻结 kernel。Central4 descriptor 可用于现有支持它的编译路径，但 pressure constant-density 对 Central4 仍 **Unsupported**；没有静默降级。原 `nu * laplacian` 公式保持。

## 17. Gravity provider

**Implemented**：pressure primitive Momentum 通过 `source.gravity.primitive` 读取源配置、区域与 geometry，返回加速度并只加一次 `dt*a`。conservative Momentum 继续使用原 `rho*a` 和能量功率实现；编译端用同一数学 term 的不同 storage contract 选择 provider。PISO/SIMPLE/PIMPLE 均无 gravity 特例。

## 18. MRF provider

**Implemented**：pressure primitive Momentum 通过 `source.mrf.primitive` 复用既有旋转参考系公式的纯加速度入口；测试验证 Coriolis 符号/量值。conservative MRF source math 未改，pressure operation 不检查 MRF model name。

## 19. Turbulence stress provider

**Unsupported**：当前 constant-density pressure runtime 没有完整的 turbulent effective stress 与 k/epsilon/omega transported equation execution。启用 RAS 时 explain 明确列 `closure.turbulence` 与 `equation.turbulence-transport` 不可用，run fail-fast。**Interface-only** 的 term catalog 可容纳将来注册的 stress provider，但本轮没有假造“仅有黏度项就可运行”的 turbulent PISO。

## 20. Custom source provider

**Implemented**（C++ composition）：测试给 raw Momentum 增加 `source(customForce)`，注册不占用内置 SourceKind 的 primitive callback，确认 compilation Runnable、目标为 transformed predictor、SolvePlan operations 不变。native YAML 的任意新 source schema 仍 **Interface-only**；没有 expression interpreter。

## 21. Conservative-path migration

**Legacy** fused conservative RHS 保持原 source dispatcher 与 WENO/TENO workspace/顺序；NumericalCompiler 已不按 gravity/MRF/wallHeat 名称作中央选择。先收敛 binding authority，不在这一轮重写已验证的 conservative source formulas。

## 22. Pressure-constrained Momentum migration

**Implemented**：`PressureOperators::bind` 只消费编译后的 Momentum terms 与已解析 callback，验证唯一 Upwind1、至多一个 Central2、source 必须在 RHS 且具有 primitive provider。assembly 按固定次序计算 convection、diffusion 和 source；source 在 fixed-time outer iteration 中按当前 `U_k` 重评估，但始终基于 `U_n` 加一次 `dt`，不重复累计物理推进。

## 23. Runtime requirement derivation

**Implemented**：compiled term provider ID 进入 `providerRequirements`；操作 provider 与 term provider 共同出现在 explain/validation。普通 term provider 没有 MPI include，已有 halo/canonical face/reduction 仍由 ExecutionRuntime 完成。若以后 term 需要额外 halo 或 reduction，应扩展 descriptor 的明确 requirement，本轮 gravity/MRF 没有新增通信需求。

## 24. Serial tests

**Implemented**：`pressureTermProviders` 新 CTest 在 PISO/SIMPLE/PIMPLE 分别比较无源、gravity、MRF：输出有限且速度差非零，物理时钟和 CompiledSolvePlan 相同；turbulence case 验证 Unsupported 与 run 拒绝。`providerCatalogResolution` 覆盖 role、unknown、recipe、歧义、贡献生命周期及 MRF 公式。

## 25. MPI tests

**Implemented**：gravity 2-rank 与串行采用相同 16×16 case 和 `g=0.01`，首压力矩阵/RHS 按物理 row 一致；最终 `U` 最大差 `5.71e-11`，pressure shape 最大差 `6.81e-08`，在既有 `1e-8/1e-6` 容差内。`distributedPressureRegression` 包含 2×1、1×2、2×2 及 SIMPLE/PIMPLE 组合，在 host 通过。

## 26. Conservative regressions

**Implemented**：使用保存于 `/private/tmp/sonic-phase28a1-baseline-20260927/sonicSolver` 的更早独立 binary 与当前文件运行相同 WENO7/Rusanov 两步案例；gravity、MRF、wallHeat 各自的 dt 序列和最终 VTS SHA 完全一致。它不是紧邻 Phase 28B 的二进制快照，因此该证据还需与 Phase 28A.1 已通过的 frozen regressions 一起解读。对应最终 SHA 为 `553c00f3…`, `fb13e31f…`, `884c278d…`。wallHeat explain 显示 `source.wallHeat.energy`，没有接入 pressure-only Momentum。

## 27. Pressure regressions

**Implemented**：既有 `pisoNumericalRegression`、`constantDensityPressureRegression`、Poiseuille、fixed-time pressure、distributed pressure 全部通过。primitive Upwind1/Central2 kernel 提取前后的 PISO control、PISO gravity/MRF、SIMPLE gravity、PIMPLE MRF 两步 VTS SHA 逐一一致。legacy PISO frozen SHA 仍为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`。

## 28. Production regressions

**Implemented**：4-rank Sod 20 步、serial Ghost IBM 20 步、Eulerian 2 步与上述独立 binary 输出 VTS SHA 一致；Sod/Ghost 的 dt 序列一致。`sodRusanovNumericalRegression` 冻结 WENO7/Rusanov 全程通过。RPI wall-boiling 在该旧 binary 与当前版均因未声明 `S_EE_PIMPLE` execution policy 而 fail-fast，是既有 **Unsupported** case，不能当作通过的 wall-boiling 回归。

## 29. Architecture guards

**Implemented**：checker 禁止 central NumericalCompiler、provider resolver、SingleFluidStepper、PressureOperators 恢复 gravity/MRF/wallHeat model-name dispatch 或 `E_MOMENTUM` prefix role 推断，并检查 SystemContribution 不再持有 raw 引用。架构 checker 通过，allowlist 仍为 2，无新增反向依赖。

## 30. Remaining Legacy

- **Legacy** conservative fused convection/source executor、Eulerian specialized equation adapter、fixed-time pressure transient closure。
- **Legacy** 用于 explain/一致性校验的 SolveStrategyKind metadata；它不再分派 pressure/Eulerian runtime。
- **Legacy** Field 仍同时保存 state、geometry、boundary metadata、GlobalDof/IBM 分类和 thermo cache；本轮未拆。

## 31. Remaining Unsupported

- **Unsupported** pressure turbulent transported equations及其完整 stress coupling；arbitrary multi-block pressure、empty-row MPI rank、非正交 pressure correction。
- **Unsupported** pressure + IBM 双 constraint、通用隐式时间积分、任意新 native IO source schema。
- **Unsupported** 当前 RPI wall-boiling fixture 的 `S_EE_PIMPLE` policy 缺失；本轮不修改 Eulerian execution semantics。

## 32. Recommended next step

优先补齐独立的 turbulent Momentum stress 与其 transport execution contract，然后考虑 `Momentum + D U=0 + J U=U_b` 的 Pressure/IBM constraint composition。不要建立 PISO/SIMPLE/PIMPLE 各自的物理源项分支。

最终核验：clean build、host CTest、架构 checker、`git diff --check`；数值/并行数学与 Phase 28B 前基线一致。未开始 Pressure+IBM、implicit time、VOF 或 Field 重设计。
