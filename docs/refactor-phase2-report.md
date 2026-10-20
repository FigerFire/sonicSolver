# Phase 2: Solver State and Services Consolidation

## 1. State Ownership Before

`SolverState` duplicated bundle, field/fields and clock members. `StateBundle`
also held canonical, patches, clock, equation binding and registry. Steppers
selected and rewrote those entries through fallbacks.

## 2. State Ownership After

`StateBundle` is the physical-state contract: patches, transported variables,
equation binding, distributed views, time, dt and step. `SolverState` now only
holds the bundle pointer and the driver's maximum timestep. No Field is copied.

## 3. Patch State Representation

| Before | After |
| --- | --- |
| `field`, `fields`, `canonical`, `patches` selected state | `patches` is the only enumeration |
| `canonical` single-patch compatibility pointer | removed |
| empty patches fell back | validation fails fast; `singlePatch()` checks size 1 |

`canonical` here was never MPI canonical-entity ownership. MPI COPY/SUM and
canonical-face semantics remain in ExecutionRuntime and ParallelCoordinator.

## 4. Time Ownership

Before, solver state, bundle, and all three algorithms had writable clocks;
the first advance copied state into algorithm counters and later copied it
back. After, steppers read and commit `StateBundle::{time,dt,step}` directly.
`StepResult` is only a post-commit snapshot. Restart comes from the bundle;
stage, IBM target-time, boundary, and observer expressions retain existing
`time`, `time + dt`, and stage-fraction timing.

## 5. VariableRegistry Ownership

Single-fluid now registers directly into `StateBundle::transported`; coupling,
integration, output, and distributed registration share it. The copied
`multiPhaseVariables -> bundle.transported` path is removed. Multi-patch keeps
one registry per incompatible local patch, with no bundle duplicate.

## 6. EquationSet / EOS Binding

`CompressibleAlgorithm::equationSet_` and its setter are removed. The bundle
bears the density EquationSet; `validatePatches()` checks its identity against
participating Fields. CFL, Rusanov flux and EOS validation consume it.

## 7. SolverServices Binding

Applications now bind once before `Time::Driver`:

```cpp
algorithm.bindServices(services);
algorithm.advance(state);
```

`advance` no longer receives a service table or overrides/restores pointers.

## 8. Removed Setter APIs

Algorithm runtime, IBM, equation-system, observer, boundary-pipeline,
transport, transported-variable, and EquationSet setters were removed.
`Execution::Runtime::setParallelCoordinator` was unused and removed.

## 9. Removed Legacy APIs

Removed `setMultiPhaseAuxFields`, `legacyTransportedVariables`, and unused
`ITimeStepReducer`/`LocalTimeStepReducer`. Runtime `globalMinimum` is the sole
global-CFL reduction contract.

## 10. step / advance Public API

`advance(SolverState&)` is the public timestep entry. Former public `step`
overloads became private `stepImpl` bodies without changing their operations.

## 11. Compatibility Branch Count

Active-field fallback branches: **5 -> 0**. Temporary service
override/restore paths: **3 -> 0**. Boundary-pipeline/manual ordering remains
unchanged for Phase 3.

## 12. State Authoritative Source Count

| Metric | Before | After |
| --- | ---: | ---: |
| SolverState state-entry members | 7 | 2 |
| writable clock representations | 6 | 1 |
| single-fluid registry sources | 2 | 1 |
| density EOS binding sources | 3 | 1 logical binding |
| patch collection forms | 4 | 1 |

## 13. Service Authoritative Path Count

| Metric | Before | After |
| --- | ---: | ---: |
| algorithm service setters | 16 | 0 |
| temporary overrides | 3 | 0 |
| timestep service paths | 2 | 1 |
| dt global-reduction paths | 2 | 1 |

## 14. Runtime / MPI Semantics

Unchanged: active bundle attachment, halo freshness, COPY, residual SUM, and
canonical-face assembly retain their existing contracts and order.

## 15. Numerical Lifecycle

Unchanged: no numerical formula, RK coefficient/stage time, boundary/halo/IBM
order, pressure corrector, KKT structure, data layout, or diagnostics timing
changed.

## 16. Build/Test Before/After

| Check | Before | After |
| --- | --- | --- |
| architecture checker | 10 allowlisted edges | 10 allowlisted edges |
| Debug build | passed | passed |
| CTest | 6/6 | 6/6 |
| 4-rank Sod | passed | passed |
| RPI wall boiling, one step | passed | passed; same `dt=1e-5` summary |
| serial IBM smoke | 7 pass; 3 known KKT stops | 7/7 passed |

## 17. Deferred Issues

Single/multi lifecycle merge, methods/numerics time relocation, scalar
transport equation ownership, model derivative ownership, viscous relocation,
phase boundary ownership, level-set HJ-WENO ownership, Field split, IBM
architecture, distributed monolithic KKT, HYPRE backend ownership, and
Equation DSL/AssemblyPlan remain deferred. Phase 3 has not started.
