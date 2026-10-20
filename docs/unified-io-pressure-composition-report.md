# 统一 IO 与压力约束 Composition 迁移报告

日期：2026-09-17

本报告的状态词严格限定为：**Implemented**、**Interface-only**、**Legacy**、**Unsupported**。

## 1. Previous IO

**Legacy**。原生 case 主要由 `solvers/solvers.yaml`、`models/models.yaml` 和兼容字典构造 `SolverConfig`；`type densityBase/pressureBase` 仍是旧 runtime 的兼容输入。

## 2. New user-facing IO

**Implemented**。原生 loader 现在可选读取：

```text
equations/equations.yaml
algorithms/algorithms.yaml
models/thermoDynamics.yaml
```

它们被绑定到中立的 `EquationCompositionConfig`，而不是向用户暴露 Raw/Executable system、transformer 或 plan。旧目录仍可共存，未移动 production case。

## 3. Built-in field catalog

**Implemented**。`rho`、`U`、`p`、`T`、`h`、`e`、`alpha`、`k`、`omega`、`epsilon` 可以只声明 `output`、`outputName`、`restart` 的 IO override。未知 symbol 必须有完整 field 声明；`new: false` 显式报错。完整旧式声明继续按 **Legacy** 兼容。

Field 元数据未加入 `conserved`、`compressible` 或 `pressureBased` 标记。

## 4. Equation registry

**Implemented**（builtin），**Interface-only**（custom）。`use: [Momentum, Continuity, Energy]` 被解析为 equation composition；未知 builtin 和文件形式的 custom equation 都 fail-fast。`add/extend/replace/disable` 的输入位置已保留，但尚无 custom equation lowering，不能静默忽略。

## 5. thermoDynamics registry

**Implemented**。`models/thermoDynamics.yaml` 是单流体物性 composition 的唯一新入口，包含 `equationOfState`、可选 `thermo`、可选 `transport` 与对应 properties。

## 6. equationOfState providers

**Implemented**（composition/validation）。支持 `perfectGas` 与 `rhoConst`。`rhoConst` 要求正的 `properties.equationOfState.rho`，并解析为 constant density 与 zero thermodynamic compressibility。

## 7. caloric thermo providers

**Implemented**（selection/validation），**Legacy**（runtime material model）。支持选择 `hConst` 与 `janaf`；Energy equation 缺少 caloric provider 时立即报：`Energy equation requires caloric thermodynamics...`。现有 SingleFluid EquationSet 仍是旧数值物性实现。

## 8. transport providers

**Implemented**（selection/validation），**Legacy**（runtime material model）。支持 `const` 与 `sutherland` selection。`rhoConst` Momentum 要求显式 transport；没有 transport 不会自动假定零黏度。

## 9. removal of energy key

**Implemented**。`thermoDynamics.energy` 被明确拒绝；能量、焓或内能的求解变量属于 equations registry。

## 10. Algorithm registry

**Implemented**。`Explicit`、`SIMPLE`、`PISO`、`PIMPLE` 被解析为 user-facing preset。内部仅向 composition 提供 transformation request 和 execution policy；输入不出现 transformer 名称。

## 11. rhoConst lowering

**Implemented**（数学 composition），**Unsupported**（数值执行）。

```text
Momentum + Continuity + rhoConst
  -> raw U, p, rho closure
  -> div(U)=0 constraint
  -> PressureConstraint transformation
```

该路径不会生成 `rhoE`，也不会用 PerfectGas 伪造密度可压缩性。它尚未绑定常密度 U/p predictor 和 pressure matrix，故 backend 为 `Unsupported`，application 在 explain 后 fail-fast。

## 12. PressureConstraint transformation

**Implemented**。一个 shared `PressureConstraintTransformer` 由 Momentum、`C_INCOMPRESSIBILITY`、pressure multiplier 的数学 contract 触发；不检查 density/pressure solver identity。`pPrime` 只在 executable system 中生成，storage 标识为 `TransientWorkspace`，不在 raw system 中。

## 13. Lid cavity

**Unsupported**。常密度 cavity 的 composition、constraint 和 PISO plan 已可解释，但没有 U/p momentum predictor、压力边界和常密度 pressure matrix 的数值算子，因此没有伪造的 cavity 数值回归、divergence 监控或中心线速度结果。

## 14. Sod

**Implemented**（composition regression），**Legacy**（runtime execution）。

```text
Continuity + Momentum + Energy + perfectGas + hConst
  -> variable rho, rhoU, rhoE
  -> no incompressibility constraint
  -> explicit conservative route
```

`test_pisoArchitecture` 验证这条 composition 不受 pressure transformation 污染。现有 Sod 的已验证高阶数值执行链未被修改，也未在本轮重新作为错误 baseline 比较。

## 15. PISO generalization

**Implemented**（PISO N、N>0 且 non-orthogonal correctors 为 0）。`ExecutionCapabilitySignature` 代替 exact equation-ID set，包含 pressure constraint、constant density、conservative state、predictor/correction binding、auxiliary schedule 与 corrector counts。PISO(2) 现在保持 generic plan 路由；PISO(1) 的既有 operation order 不变。

**Interface-only**：PISO non-orthogonal loop。现有 corrector workspace 不支持在不改变 matrix/velocity-correction 语义的情况下重复 pressure solve，因而不自动展开为错误循环。

## 16. SIMPLE/PIMPLE

**Legacy**。它们仍共享 PressureConstraint transformation 与 structured execution policy，但 numerical backend 仍由 LegacyPressureExecutionAdapter 提供。outer-loop、relaxation 和 convergence 的 generic execution 尚未实现。

## 17. workspace resources

**Implemented**。pressure correction unknown `pPrime` 从 `SpecializedExecutor:pressureCorrection` 迁为 `TransientWorkspace:pressureCorrection`；它不可 restart、不可 output，且不是 physical state authority。

## 18. architecture guards

**Implemented**。`tools/check_architecture.py` 新增检查：generated pressure correction 必须是 transient workspace；generic pressure routing 必须拥有 `ExecutionCapabilitySignature`；旧 `isMinimalGenericPiso` exact-ID routing 被禁止。现有 guard 继续禁止 transformer 直接 MPI 和根据 density/pressure solver identity 判断 transformer applicability。

## 19. legacy compatibility

**Legacy**。`SolverConfig.numerics.solver` 和旧 pressure/density adapters 暂存以维持 production case。它们不覆盖已声明的 semantic composition；constant-density composition 不会被重定向到 legacy pressure/PerfectGas route。

## 20. interface-only features

- custom equation parser/lowering；
- explicit `add/extend/replace/disable` lowering；
- PISO non-orthogonal loop；
- generic SIMPLE/PIMPLE executor；
- generic auxiliary scalar scheduling；
- caloric/transport provider 到现有 numerical EOS 的完全解绑。

## 21. unsupported features

- serial或MPI constant-density U/p predictor、pressure boundary、Poisson matrix 与 cavity runtime；
- distributed generic pressure plan；
- generic IBM/KKT lowering；
- 自定义 equation 的数值执行。

所有这些组合均应 fail-fast；没有 Rusanov、PerfectGas、legacy pressure 或其它算法替代 fallback。

## 22. next work

先实现并独立验证 constant-density momentum predictor 与 U/p pressure-correction operator，再接入 lid-driven cavity：边界闭合、`max/L2 div(U)`、momentum/pressure residual、中心线速度。随后将 PISO non-orthogonal loop、SIMPLE/PIMPLE outer loop 逐项绑定到相同 generated equation/operator，而不是复制 pressure mathematics。

## 验证

- `cmake --build build-piso-tests --target SF_caseAdapter --parallel 4`：通过。
- `cmake --build build-piso-tests --target SF_application --parallel 4`：通过。
- `ctest --test-dir build-piso-tests -R pisoArchitecture --output-on-failure`：通过。
- `python3 tools/check_architecture.py --quiet`：通过；allowlisted dependency edges = 10。

未运行 lid cavity、完整 Sod 数值、MPI 或 IBM smoke，因为本轮没有完成其所需的常密度 runtime，且不应把 skeleton composition 当作数值验证。

## 最终职责表

| 对象 | authoritative responsibility |
| --- | --- |
| Field | symbol、dimensions、IO、initial/boundary metadata |
| Equation | mathematics、transported/conservative/constraint role、solved variable |
| equationOfState | p-rho-T relation 与 compressibility relation |
| thermo | Cp/Cv/h/e caloric relations |
| transport | mu/nu/kappa 等输运关系 |
| Algorithm | transformation request 与 execution policy |
| Numerics | operator discretization、time integration、linear backend |
| SolvePlanner | execution order、loops、block solve、subcycle |
| Runtime | execute compiled plan；不重定义数学系统 |
