# Phase 28C：Model-Owned Provider 与 Contribution 收口

状态用语：**Implemented** 表示已有生产执行和回归证据；**Interface-only** 表示只有组合/编译接口；**Legacy** 表示仍由既有专用数值路径执行；**Unsupported** 表示明确拒绝。

## 1. Scope

**Implemented**：把 gravity、MRF、wallHeat 的 term provider 注册和参数编译移到 model 侧；把 model contribution 的选择移到 application；pressure Momentum 与 conservative source 执行均消费冻结的 term binding。没有新增物理公式、改变 PISO/SIMPLE/PIMPLE Plan、压力矩阵、Rhie–Chow、stage time、MPI COPY/SUM 或 Field 状态布局。

## 2. Frozen baseline

修改前 HEAD `537889b05234c34217aa1a8eff7ab7844b69b4b9`，且工作树已有 Phase 28A/28B 未提交修改。快照、dirty patch、旧二进制、构建/CTest/checker 日志保存在 `/private/tmp/sonic-phase28c-baseline-20260928/`。修改前完整 build 成功，host CTest **25/25**；sandbox CTest 的 3 个 PRTE socket 失败是运行环境限制，host 重跑通过。架构 allowlist 为 2 条既有 `methods → solver` 边。

## 3. solver/system model dependencies before

**Legacy（修改前）**：`SF_systemBuilder.cpp` include 并直接调用 physics/IBM contributor；`SF_buildRequest.h` include LevelSet/Turbulence 具体 spec；`SF_termProviderCatalog.cpp` include body-acceleration model 并列出 gravity/MRF/wallHeat；`SF_solverSystem` CMake 直接链接 `SF_physics`、`SF_ibm`、`SF_turbulence`。共 5 条直接 model include、3 条直接 model target link。

## 4. solver/system dependencies after

**Implemented**：`solver/system` 中 `#include "models/..."` 为 0，`SF_solverSystem` 不直接链接上述三个 model target。`BuildRequest` 消费中立 `SystemContribution` 和 `SourceTermProviderDescriptor` 值。model header 只依赖 `core/system`/`core/interfaces` 中立类型，不 include `solver/system` 实现。旧 Eulerian/legacy mixture template lowering 仍属于 **Legacy**，本轮未重写。

## 5. Provider ownership before

**Legacy（修改前）**：通用 `TermProviderCatalog::builtIn()` 同时列数值格式和三个具体 model，pressure source callback 每次接收整份 `SourceConfig`；conservative `SourceTerm` 另有三个 model 的中央 dispatcher。

## 6. Provider ownership after

**Implemented**：通用 catalog 只登记 primitive Upwind、conservative convection 与 Central diffusion。各 model 目录分别提供 source descriptor：`gravity/SF_accelerationProvider.*`、`mrf/SF_frameProvider.*`、`heat/SF_wallFluxProvider.*`。application 把它们作为值放进 `BuildRequest`，NumericalCompiler 解析唯一匹配并冻结 kernel。零匹配和多匹配继续显式 Unsupported/Invalid。

## 7. Gravity provider migration

**Implemented**：`source.gravity.primitive` 和 `source.gravity.conservative` ID 不变。zone/acceleration 以值捕获，编译后原 config 可销毁。primitive 返回原加速度；conservative 仍调用原 `Source::Gravity::addTranslation`，保持 `rho*a` 与能量功率数学。没有 runtime `SourceConfig` 查找。

## 8. MRF provider migration

**Implemented**：两个原 ID 不变。编译时冻结 zone、center、反向 frame `omega`、frame velocity；primitive 沿用 `rotatingFrameAcceleration`，conservative 沿用 `addRotatingFrame`。Coriolis 符号与数值单测、生产 A/B 已覆盖。

## 9. WallHeat provider migration

**Implemented**：`source.wallHeat.energy` 保留原 ID，仅匹配 Energy/Enthalpy role；conservative kernel 沿用 `Source::WallHeat::addSource`。pressure-only 系统没有相应方程时，model term 明确失败，不能静默忽略。新的压力能量方程数值执行仍 **Unsupported**。

## 10. Provider compiled data contract

**Implemented**：中立 descriptor 提供 matcher、recipe、owner 和按值捕获的 primitive/conservative kernel factory；`ResolvedTermProvider`/`BoundTerm` 保存编译后 callback、provider owner 与 data availability。没有 `shared_ptr<void>`、服务定位器或虚类 payload 层级。测试在原 gravity/MRF config 及自定义 model 对象退出作用域后解析并执行 kernel。

## 11. SourceConfig removal

**Implemented**：`CompiledPressureOperatorConfig`、`PressureOperators`、`Compressible::AssemblyContext` 不再持有/传递整份 `SourceConfig`。保守格式的 source 清空时机仍在 diffusion/source 装配点；通用 `SourceTerm::Sp` 只依声明顺序执行编译 kernel，不含 SourceKind/model 清单。

## 12. PressureOperators dependency before/after

**Implemented**：原先的 `momentumSources_`、独立 convection/diffusion callback 与 `config_.sources` 改为一份按声明 ordinal 保存的 `BoundTerm` 序列。assembly 仍按原浮点顺序先累计 source acceleration，再求每分量的 convection/diffusion，最后加 `dt*a`；矩阵和 Rhie–Chow 代码未改。

## 13. BuildRequest before/after

**Implemented**：去掉 LevelSet/Turbulence 具体 spec 和 IBM descriptor 指针；增加中立 model/user contribution 列表与 source descriptor 列表。旧物理配置只在 application/model composition 边界解释。尚存的 homogeneous/legacy mixture 配置与 capability 检查是 **Legacy**，没有作为新 model 扩展接口。

## 14. Model contribution lowering

**Implemented**：application 生成 Source、LevelSet、Turbulence、IBM 的纯值 contribution；SystemBuilder 统一 apply。model contribution 自带 provider requirement（例如 `equation.level-set`、`closure.turbulence`、`ibm.constraint`），SystemBuilder 不调用具体 contributor。role-targeted source 只在实际 equation role 上扩展，缺目标时 fail-fast。

## 15. Custom model extensibility proof

**Implemented**：`test_pisoArchitecture.cpp` 中 `TestBodyForceModel` 独立贡献 Momentum term、编译数据为 0.01 的 provider，并在 model 对象销毁后通过既有 compiler/plan 路径解析和调用 kernel；无需修改 NumericalCompiler、SystemBuilder、PressureOperators、PISO/SIMPLE/PIMPLE 或 PlanExecutor 的 model 分派。native YAML 自定义 schema 仍是 **Interface-only**。

## 16. Serial gravity/MRF regression

**Implemented**：PISO gravity、SIMPLE gravity、PISO MRF、PIMPLE MRF 两步 A/B 的全部 VTS、dt 序列与首压力矩阵文件 SHA 一致。`pressureTermProviders` 定向测试通过。

## 17. MPI gravity/MRF regression

**Implemented**：host 上 gravity 的 serial/MPI-2 对照首压力矩阵与 RHS 等价；最终 U 最大差 `5.70890329e-11`，pressure shape 差 `6.80884114e-08`，均在冻结容差内。完整 distributed pressure CTest 的 cavity/channel 2×1、1×2、2×2 与 SIMPLE/PIMPLE 组合通过。MRF 本轮执行了 serial PISO/PIMPLE A/B；未单独声称 MRF MPI 数值基线。IBM 的 FictitiousDomainMPI2 与 VelocityForcingMPI2 短副本在新旧 binary 中均运行 46 步，VTS/dt 字节相同；PeskinMPI2 两版均在 rank 1 因 `Peskin interpolation produced a singular mass diagonal` fail-stop，这是实际数值失败而非 PRTE 启动故障。

## 18. Conservative regression

**Implemented**：两步 conservative gravity/MRF/wallHeat 的输出 VTS 与 dt 序列均与本轮冻结 binary 字节一致；冻结的 WENO7/Rusanov CTest 通过。保守格式 WENO/TENO convection 本身仍为 **Legacy** fused kernel，本轮只迁移源项绑定。

## 19. Production regressions

**Implemented**：最终 binary 与冻结二进制对 10 组生产输入重新运行 A/B：PISO/SIMPLE gravity、PISO/PIMPLE MRF、conservative gravity/MRF/wallHeat、4-rank Sod 20 步、serial Ghost IBM 20 步、Eulerian 2 步，全部 VTS、可用 dt/clock 序列及首压力矩阵 SHA 字节一致。model provider 文件定名后重新清洁构建并再次执行这 10 组 A/B，结果仍全部一致。对 `test/IBM/cylinderFlow*` 中全部 10 个串行内置方法制作 `endTime=0.02` 短副本，7 个可运行方法的 VTS/dt 与旧 binary 相同；DFM AugmentedLagrangian、DFM Implicit SelfPropelled、FullyImplicitDLM 三个 fixture 在新旧 binary 中均因已选择但未被方程消费的 diffusion recipe 同样 fail-fast，未误报为本轮可用。最终完整 `cmake --build build --clean-first --parallel 4` 成功；重命名后 host CTest 再次 **25/25**，其中 MPI/HYPRE 测试实际运行；架构 checker 与 `git diff --check` 通过。legacy `pressureConstraintPiso` 最终 VTS SHA-256 为 `4bc3e0bef4f252614aa4e68651c38c634baeec5fa28ca6ffae5ff1215cf72511`，与冻结值完全一致。

## 20. Architecture guards

**Implemented**：checker 禁止 `solver/system → models/*` include、`SF_solverSystem` 直接 model target link、generic catalog/SourceTerm 中的具体 model 清单、PressureOperators 中 raw model config、SystemBuilder 中具体 contributor 调用。旧的 2 条依赖 allowlist 没有为本轮扩张。

## 21. Remaining Legacy

**Legacy**：homogeneous/legacy mixture 的旧 template 和 capability lowering、Eulerian specialized equation adapter、保守格式 fused WENO/TENO convection、fixed-time pressure predictor 的专用数学、`SolveStrategyKind` explain metadata。本文不把这些称为通用 equation runtime。

## 22. Remaining Unsupported

**Unsupported**：pressure turbulent stress 与 k/epsilon/omega transport、pressure+IBM 双 constraint、通用隐式时间、非正交压力、多本地 patch/任意多块 pressure；pressure-only wallHeat 缺 Energy/Enthalpy 方程时明确拒绝。上述 3 个 IBM fixture 的 diffusion recipe 不匹配、PeskinMPI2 的奇异质量对角线和 RPI wall-boiling fixture 原有 `S_EE_PIMPLE` policy 缺失均不在本轮数值能力扩张范围；没有为让它们运行而加 fallback。

## 23. Readiness for turbulence phase

**Interface-only**：中立 model contribution/provider descriptor 可承载未来 Momentum stress term 及 transport equations，但尚无 turbulent stress provider 数值实现。后续须先定义张量应力数学 term，再接 k/epsilon/omega transport；不能把 scalar `nuEff` 作为唯一 turbulence contract。
