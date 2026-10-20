# Executable Equation Authority 审计

日期：2026-09-15

> 本文记录修改前基线；完成后的 ownership 与验证结果见
> `docs/executable-equation-authority-refactor-report.md`。

## 当前 authority map

| 数学对象 | Description authority | Runtime authority | 重复/缺失 |
|---|---|---|---|
| density mass/momentum/energy | `System::ResolvedSimulationSystem::equations` | `Equation::Compressible::System::definition_` | 三个方程被独立创建两次；runtime 只按 `TermKind` 判断整组算子是否执行 |
| Eulerian phase continuity/momentum/enthalpy | `System::build()` 按 phase 名创建 descriptor | `PhaseEquationAssembler` 的固定成员函数及 `PressureStepper::stepImpl()` | description 与 assembly sequence 独立维护 |
| pressure/shared-pressure/volume-fraction constraints | resolved equations/constraints/`S_EE_PIMPLE` | `PressureStepper` correction 循环与 `PhaseEquationAssembler` | solve block 当前不驱动 runtime |
| turbulence transport | 仅 closure 文本 | `Turbulence::EquationSystem` 根据 config/model 隐式创建 k/epsilon/omega storage 并求解 | k/epsilon/omega unknown、equation、solve block 未注册 |
| density gravity/MRF/wall heat | resolved system 未描述具体 terms | `SourceTerm::Sp` 中 `switch(SourceKind)` | central execution switch 是唯一 runtime selection authority |
| Eulerian gravity/MRF/wall heat | resolved system 未描述具体 terms | `PressureStepper` constructor switch 创建 `PhaseSourceRegistry` | 有 registry，但 composition 仍由 central switch 决定 |
| phase change | capability/closure 级描述 | homogeneous coupling 或 Eulerian phase/source 实现 | contribution 没有关联到目标 equation definition |
| IBM | resolved unknown/equation/constraint/solve block | specialized IBM ports/coupling | 数学 identity 已注册，assembly binding 仍由 adapter 完成 |

## 问题回答

1. **只存在于 ResolvedSimulationSystem 的 equation**：`E_LEVEL_SET`、`E_PHASE_MASS`、`E_LEGACY_ALPHA` 以及 IBM descriptor 导入的部分 equation identity；它们的实际执行由 specialized coupling/model path 完成。
2. **只存在于 runtime specialized system 的 equation**：具体 turbulence transport equations；`Turbulence::EquationSystem` 按模型隐式决定 k/epsilon/omega。
3. **重复描述的 equation**：density mass/momentum/energy；Eulerian phase continuity/momentum/enthalpy；pressure correction 与共享压力关系。
4. **使用 registry 的 source**：Eulerian `PhaseSourceRegistry` 按注册顺序调用 `PhaseEquationSource::add`。
5. **使用 central switch 的 source**：density `SourceTerm::Sp`；Eulerian `PressureStepper` constructor 也用 `SourceKind` switch 构造 registry。
6. **隐式 transport equation**：RAS k-epsilon、RAS k-omega SST；homogeneous phase mass、legacy alpha、level-set transport 只有粗粒度 descriptor，尚未与 executable definition 绑定。
7. **Workflow::Plan runtime consumer**：application 日志、ModuleGraph/mesh capability validation 和 execution builder preparation；`runFlow` 与 stepper 尚不遍历 stages。
8. **solveBlocks runtime consumer**：Workflow builder 与 application executor selection；stepper 内部仍固定实现真正的 block sequence。

## 调用链

```text
ResolvedSimulationSystem
  -> Workflow::Plan / check / explain / application executor selection

CompressibleAlgorithm
  -> owned Equation::Compressible::System
  -> DensityBasedRHS
  -> System::convection / diffusionAndSources
  -> SourceTerm::Sp switch

PressureStepper
  -> owned PhaseEquationAssembler
  -> owned Turbulence::EquationSystem
  -> owned PhaseSourceRegistry
  -> hard-coded PIMPLE loop
```

## 迁移顺序

1. 在现有 `EquationDescriptor` 上关联 `Equation::Definition`，不新建平行 IR。
2. 让 compressible transitional adapter 读取 resolved definitions。
3. 引入窄的 contribution registration，先迁移 source selection 并保持注册顺序。
4. 由具体 turbulence model 注册 unknown/equation/closure/solve block。
5. 从 definition 和 numerics lower 出 AssemblyPlan；specialized kernel 保留。
6. 让 stepper 校验并消费 Workflow solve-block stages，再执行原有数学顺序。

所有阶段保持 state、workspace、boundary/halo、stage time、canonical COPY、
GlobalDof SUM 和 source application order 不变。
