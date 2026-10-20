# Phase 2 State and Service Consolidation Plan

## Scope and invariants

This phase changes ownership, non-owning references, and public timestep
entrypoints.  It does not change field storage, numerical formulas, RK
coefficients or stage times, physical-boundary/halo/IBM order, canonical-face
assembly, global-DOF SUM, pressure correction, or KKT layout.  `StateBundle`
is the existing container closest to an authoritative simulation state, so no
`SimulationState`, context, adapter, locator, or new service abstraction will
be introduced.

The pre-change architecture checker reports 10 allowlisted dependency edges.
They are unrelated to this phase and are not a target of the change.

## 1. Solver-state truth table (before)

| Representation | Creator | Writer | Reader | Authoritative? | Alias? | Can diverge? |
| --- | --- | --- | --- | --- | --- | --- |
| `SolverState::bundle` | application runners | application only | all three steppers | Intended, but optional | no | Yes: all other `SolverState` members can select a different field or clock. |
| `SolverState::field` | single-fluid/Eulerian runners; steppers write it back | application and steppers | compressible, pressure, multi-patch fallback | No | Intended single-patch alias | Yes. |
| `SolverState::fields` | multi-patch runner; steppers write it back | application and multi-patch stepper | multi-patch and compressible fallback | No | Intended patch alias | Yes. |
| `StateBundle::canonical` | all runners | applications and all three steppers | compressible/pressure and bundle fallback | No | Single-patch compatibility pointer | Yes: it is independent from `patches.front()`. It has no MPI-canonical-entity meaning. |
| `StateBundle::patches` | all runners | application and multi-patch stepper | multi-patch; distributed registration | Intended patch collection | no | Yes: empty collection falls back to `canonical`; multi-patch writes it back. |
| `StateBundle::transported` | single-fluid runner copies `multiPhaseVariables` | local registry and copied bundle registry can both be mutated | compressible stepper | Intended only after copy | No, it is a copy | Yes. |
| `VariableRegistry` / per-patch `variables` | single/multi application composition | registration code | equation coupling, VTK, distributed registration | Per workflow/patched state | no | Single-fluid has a copied second registry; multi-patch has one registry per incompatible patch field. |

Current lifetimes are application-owned `Field`, model, registry, runtime, and
services. `StateBundle` is non-owning for field storage and holds its own
registry descriptions.  It has no general consistency check.  `Field` storage
is never copied by the current state path.

### Planned state ownership

`StateBundle::patches` will be the only patch enumeration, including
single-patch runs (`size()==1`).  `canonical` will be removed because its only
observed role is the single-patch compatibility pointer, not geometry or MPI
ownership.  `StateBundle::singlePatch()` will fail fast unless exactly one
non-null patch is present.  `SolverState` will become a thin timestep request:
the bundle pointer and the driver-provided maximum timestep only.  It will no
longer carry field, fields, time, dt, or step.

For the single-fluid workflow the existing `StateBundle::transported` registry
will be populated directly and passed by reference to coupling, integration,
distributed registration, and output.  It will no longer be copied from
`multiPhaseVariables`.  Multi-patch retains one registry per local patch until
the separately deferred scalar-transport ownership work; there is no duplicate
bundle registry for those entries today.

## 2. Time truth table (before)

| Representation | Initialization / writer | Readers | Issue |
| --- | --- | --- | --- |
| `StateBundle::{time,dt,step}` | application initializes time; steppers write all three after success | output through `StepResult`; runtime has no clock use | Intended state clock but only optional. |
| `SolverState::{time,dt,step}` | application initializes time; steppers write all three | steppers initialise their internal counters from it | A second writable clock. |
| `CompressibleAlgorithm::{physicalTime_,dt_,stepCounter_}` | config start time; first `advance` overwrites time/step; `step` advances it | stage time, IBM target time, observer | A third writable clock. |
| `MultiPatchAlgorithm::{physicalTime_,dt_,stepCounter_}` | same pattern | stage/IBM time and observer | A third writable clock. |
| `PressureStepper::{physicalTime_,stepCounter_}` and local `dt` | first `advance` copies solver state; then advances | pressure step and result | A third writable clock. |
| `Time::Driver` | starts from run control and consumes `AdvanceResult` | output/termination | Driver snapshot is not state storage. |
| `StepResult` | constructed after commit | `Time::Driver` | Correctly a result snapshot, but mirrors the above sources. |

Restart starts from `caseConfig.time.startTime`, which flows into
`SolverState::time` only on the first advance, then gets copied into each
algorithm counter and written back.  IBM target time currently uses
`physicalTime_ + dt_`; observers use the same counters.  Output consumes the
post-commit `StepResult` through `Time::Driver`.

### Planned time ownership

`StateBundle::{time,dt,step}` will be the only writable physical clock.
Algorithms will bind and validate the same bundle on their first advance,
perform their existing stage sequence using those members, and commit directly
to them.  `StepResult` remains a snapshot.  The first-step restart value is
read directly from the bundle, so it is not reset to the numerical-config
start time.  The expressions used for stage, moving IBM, and boundary target
times remain `time`, `time + dt`, and the same existing values at each call
site.

## 3. Equation/EOS binding truth table (before)

| Binding | Creator / writer | Consumers | Divergence risk |
| --- | --- | --- | --- |
| `Field::equationSet()` | equation-set factory binds it before stepping | boundary, EOS state, flux helpers, validation | It is the actual field-storage binding. |
| `StateBundle::equations` | single-fluid runner copies `field.equationSet()` | compressible `advance` | May differ from Field. |
| `CompressibleAlgorithm::equationSet_` | application setter, then per-call bundle override/restore | CFL, flux, thermo validation | Can differ from both. |
| multi-patch `Equation::Compressible::System` | algorithm member | Euler flux assembly | No separately supplied bundle binding in current multi-patch application. |
| pressure-based PhaseSystem equation state | PhaseSystem construction / Field binding | pressure equations, closure | Distinct physical formulation, not a duplicate density-based EOS object. |

The plan treats `StateBundle::equations` as the simulation binding for a
density-based bundle and verifies that each participating patch has the same
shared EquationSet object when one is present.  The algorithm setter/cache will
be removed; CFL, numerical flux, and validation read the bundle binding.  The
pressure-based PhaseSystem remains outside this density-EOS unification.

## 4. Service truth table (before)

| Dependency | Constructor injection | Setter injection | Per-call service | Fallback / temporary override |
| --- | --- | --- | --- | --- |
| `SolverServices` | no | individual algorithm setters | yes | every `advance` saves, replaces non-null entries, then restores. |
| `ExecutionRuntime` | `Runtime(parallel)` | unused `Runtime::setParallelCoordinator`; algorithm setters | all steppers | algorithm pointer is overridden per step; attached bundle can change. |
| boundary applicator/pipeline | applicator built from config | none / algorithm setter | compressible and multi-patch | pipeline optional; manual order remains when absent. |
| IBM ports | adapters built by application | several algorithm setters | compressible and multi-patch | per-call non-null port set replaces prior ports. |
| transport / equation coupling / observer | application-owned | algorithm setters | compressible/multi-patch | per-call override/restore. |
| `TimeStepReducer` | no | compressible setter | compressible and pressure | competes with `ExecutionRuntime::globalMinimum`. |

### Planned service ownership

Each `INavierStokesStepper` receives exactly one explicit `SolverServices`
value through `bindServices(...)` before the time driver begins.  `advance`
will no longer accept services, override members, or restore old members.  The
struct remains a plain explicit dependency collection and does not perform
lookup, routing, or fallback selection.  Runtime remains application-owned;
the algorithm attaches the first active bundle and rejects a different bundle
after stepping starts.  The unused `Runtime::setParallelCoordinator` and the
unused `ITimeStepReducer` path will be removed; global dt reduction uses the
already-required `ExecutionRuntime::globalMinimum` when a runtime is bound.
The physical boundary pipeline and its manual fallback are intentionally left
unchanged, including their order.

## 5. API migration sequence

1. Add `StateBundle` validation and single-patch access; migrate applications
   to `patches`, then remove `canonical`.
2. Remove the single-fluid registry copy and make the bundle registry the
   shared input to integration, MPI registration, coupling, and output.
3. Make the bundle equation binding authoritative with fail-fast consistency
   checks; remove the algorithm EquationSet setter/cache.
4. Replace setter/per-call service paths with one `bindServices` path; delete
   `ITimeStepReducer` and the unused runtime coordinator setter.
5. Make `advance(SolverState&)` the only public timestep entry and move the
   existing `step` implementations private without changing their bodies or
   order.
6. Remove the unused `setMultiPhaseAuxFields` and
   `legacyTransportedVariables` after confirming no source/test call sites.
7. Remove duplicated `SolverState` clock/field members and use the bundle
   directly for restart, stage, IBM, observer, and result snapshots.

After each substep: rebuild and run focused tests.  Full post-change validation
will repeat the architecture check, build, six CTests, 4-rank Sod, one-step RPI
wall boiling, and the seven serial IBM smoke cases.  The three known monolithic
KKT IBM stops are baseline deferred failures.

## Deferred work

This plan intentionally leaves single/multi lifecycle merging, RK/Euler/SSPRK
unification, methods/numerics time relocation, scalar-transport equation
ownership, model derivative ownership, viscous relocation, phase-boundary and
level-set HJ-WENO ownership, Field responsibility split, IBM redesign,
distributed monolithic KKT, HYPRE backend ownership, and Equation DSL /
AssemblyPlan untouched.
