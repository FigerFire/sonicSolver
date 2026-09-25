# Phase 3B Plan — Density-Based Explicit Time Integration

## Scope

Phase 3A is present: `StateBundle::equations` is the density thermodynamic
binding; `DensityBasedRHS` is the sole shared spatial RHS; canonical-face and
GlobalDof barriers each have one implementation; and stage/correction
publication are explicit.  Phase 3B replaces only duplicated explicit time
integration.  It does not move final boundary/commit tails, modify spatial RHS,
or alter pressure-based stepping.

## RK truth table before modification

`R(Q)` below denotes the existing stored spatial residual with the established
sign convention, so an explicit update is `Q - dt R(Q)`.  Transported
variables use the matching `q + dt rhs` convention.

| Scheme | Single implementation | Multi implementation | State formula | RHS times | Final formula | Difference |
| --- | --- | --- | --- | --- | --- | --- |
| Euler | `methods/numerics/time` through `ddtDispatch` | `SF_multiPatchTime.cpp::stepEuler` | `Q1=Q0-dt R(Q0)` | `t` | `Q1` | same algebra; separate field/patch loops |
| SSPRK3 | `SSPRK3::ddt` through `ddtDispatch` | `stepSSPRK3` | `Q1=Q0-dt R(Q0)`; `Q2=3/4 Q0+1/4(Q1-dt R(Q1))`; `Q3=1/3 Q0+2/3(Q2-dt R(Q2))` | `t`, `t+dt`, `t+dt/2` | `Q3` | same Shu–Osher coefficients and nodes; separate storage/loops |
| RK4 | `Time::rk4` plus registry wrapper | `stepRK4` | `Q2=Q0-dt/2 R1`; `Q3=Q0-dt/2 R2`; `Q4=Q0-dt R3` | `t`, `t+dt/2`, `t+dt/2`, `t+dt` | `Q0-dt(R1+2R2+2R3+R4)/6` | same coefficients/nodes; separate storage/loops |

The SSPRK3 node order is intentionally retained although it differs from the
more common presentation of an autonomous RK tableau.  It is the existing
non-autonomous stage-time contract, not a Phase 3B correction target.

## State and registry truth table

| Representation | Meaning | Owner during a stage | Must remain temporary? |
| --- | --- | --- | --- |
| Field conservative values at timestep entry | committed `Qn` | `StateBundle::patches` | no |
| `q0`, `k1`…`k4` | RK snapshots/RHS workspaces | explicit integrator | yes |
| Field conservative values between stages | `Qstage` | explicit integrator | yes |
| `StateBundle::transported` or per-patch coupling registry | integrated variable storage/RHS | existing registry owner | values are stage state; snapshots are temporary |
| `StateBundle::time` | committed physical clock | bundle | yes: never write at stage time |

For a single patch, the bundle registry is used.  For multi-patch coupling,
the existing per-patch registry is used.  The unified integrator will select
the existing registry for each participating patch; it will not copy, merge,
or create a second registry.

## Chosen implementation

Extract the proven vector-of-patches implementation from
`SF_multiPatchTime.cpp` into one internal `solver/algorithm` function.  Its
inputs will be the existing `StateBundle`, stable `SolverServices`, selected
scheme, and explicit RHS/publish/closure callbacks.  This is ordinary function
composition, not a manager, adapter, context, or public API.

The one loop is stage outermost:

```text
snapshot Qn
for each existing stage:
  DensityBasedRHS for all patches at current stage time
  update all patches and their existing registry variables
  publish all stage state
  validate all stage state
```

`DensityBasedRHS` continues to own boundary/halo/IBM preparation and spatial
synchronization.  The integrator only applies existing tableau algebra.

## Frozen contracts

- Euler, SSPRK3, and RK4 coefficients and RHS time nodes above.
- `StateBundle::time` remains committed time during every stage.
- `beginStep` remains once per physical timestep; flow correction remains after
  all stages; final boundary/commit tails remain in their current algorithms.
- Boundary, halo, IBM, canonical COPY, GlobalDof SUM, and stage publication
  order remain unchanged.

## Planned validation

Before and after: architecture checker, build, 6 CTest, 4-rank Sod, serial
WENO7 Sod, RPI wall boiling, and IBM smoke.  Add a small deterministic stage
algebra regression that records Euler/SSPRK3/RK4 RHS nodes and verifies one and
multiple independent patches receive the same tableau updates.
