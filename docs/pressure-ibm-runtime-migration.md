# Pressure / IBM Runtime Migration

## 1. Scope

本次只迁移 Eulerian pressure runtime control authority；未改 pressure、interphase、IBM、HYPRE、boundary 或 MPI 数学。

## 2. Source audit

真实三层循环原位于 `EulerianEulerian::PressureStepper::stepImpl()`：outer、pressureCorrector、`nonOrthogonalCorrectors + 1`。

## 3. Before dependency graph

`Plan -> PressureStepper::advance -> stepImpl loops -> equation helpers`。

## 4. After dependency graph

`CompiledSolvePlan -> Run::PlanExecutor -> Run::OpRegistry -> PressureStepper numerical helper callbacks`。

## 5. Old PressureStepper lifecycle

已删除 `stepImpl()`。其原有操作被逐项保留为 callbacks：interphase、source、turbulence、momentum diagonal/flux/predictor、pressure correction/publish/sync、phase/flux correction、energy、commit 和 diagnostics。

## 6. Compiled PIMPLE Plan

`S_EE_PIMPLE` lowering 为 `EE.step`，其中 `EE.outer`、`EE.pressure`、`EE.nonOrthogonal` 是嵌套 `Loop` nodes。

## 7. Loop mapping

| Loop | Source | Plan |
|---|---|---|
| outer | `outerCorrectors` | `EE.outer.repetitions` |
| pressure corrector | `pressureCorrectors` | `EE.pressure.repetitions` |
| nonOrthogonal | `0..nonOrthogonalCorrectors` | `EE.nonOrthogonal.repetitions = N + 1` |

## 8. Operation table

`ee.step.begin`、`ee.interphase.compute`、`ee.sources.assemble`、`ee.momentum.*`、`ee.pressure.*`、`ee.phase.correct`、`ee.faceFlux.*`、`ee.energy.solve`、`ee.turbulence.*`、`ee.step.commit`。

## 9. Pressure operation providers

**Implemented**：provider callbacks 位于 `PressureStepper::registerOperations()`，调用既有 `PhaseEquationAssembler`、boundary、turbulence 和 workspace helper。

## 10. IBM classification

Ghost/ILW 是 **boundary closure**；Peskin、DFM、VelocityForcing、Brinkman 经现有 `projectPredictedState()` 是 **sequential constraint**；monolithic DLM/KKT 是 **coupled block**。

## 11. Built-in IBM integration table

| IBM method | Mathematical type | Runtime insertion | Status |
|---|---|---|---|
| Ghost + ILW | boundary closure | existing boundary pipeline | Legacy |
| Peskin | diffuse constraint | `projectPredictedState()` | Legacy |
| DFM explicit/fractional | forcing/constraint | `projectPredictedState()` | Legacy |
| FTS/BP | velocity/penalty forcing | `projectPredictedState()` | Legacy |
| monolithic DLM/KKT | coupled HYPRE block | existing monolithic API | Legacy |
| augmented Lagrangian | coupled HYPRE block | existing monolithic API | Legacy |

IBM was not inserted into Eulerian PIMPLE because the current Eulerian runner does not compose IBM ports into `PressureStepper`; claiming integration would be false. Existing density/pressure adapters remain the IBM runtime path.

## 12. Ghost IBM execution path

**Legacy** boundary pipeline; unchanged physical boundary -> ghost/ILW -> halo order.

## 13. Peskin/DFM/velocity-forcing execution path

**Legacy** existing `projectPredictedState()` dispatch; formulas unchanged.

## 14. KKT execution path

**Legacy** existing monolithic API/HYPRE implementation; distributed capability remains unchanged.

## 15. MPI ownership invariants

**Implemented**: no callback accesses MPI directly. Existing `ExecutionRuntime` continues canonical face synchronization, owner COPY and residual/constraint SUM.

## 16. Numerical invariants preserved

**Implemented**: callbacks reuse existing numerical helpers and preserve source order, boundary positions and loop counts.

## 17. Deleted authority

**Implemented**: `PressureStepper::stepImpl()` and all three solver-level loops are removed. Plan controls those loops.

## 18. Remaining Legacy

IBM integration into Eulerian PIMPLE, density RK, `IFlowAlgorithm`, BindingPass and AssemblyPlanRegistry.

## 19. Serial regression

**Unsupported**: full production baseline could not complete in this execution window; only compile and architecture/unit test evidence exists.

## 20. MPI regression

**Unsupported**: not executed in this execution window.

## 21. IBM regression

**Unsupported**: not executed in this execution window.

## 22. First-difference analysis

**Unsupported**: no production comparison exists yet.

## 23. Unsupported paths

No new IBM capability was claimed. Custom IBM providers remain Unsupported; existing distributed KKT support/fail-fast behavior is unchanged.

## 24. Next step

Run Eulerian serial/MPI baseline and migration regression, then compose existing IBM ports into Eulerian operation registration without changing their algorithms.
