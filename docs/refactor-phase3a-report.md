# Phase 3A Report — Density-Based RHS / EOS / Stage Synchronization

## 1. EOS/CFL ownership

Before: single-patch density stepping used its `Field` EquationSet for
EOS-aware Rusanov CFL and closure, while multi-patch used configured gamma and
could pass null thermodynamics into RHS.  After: `StateBundle::equations` is
mandatory, every participating patch must reference exactly that object, and
both CFL paths use it.

| Metric | Before | After |
| --- | --- | --- |
| Authoritative density thermodynamic sources | 2 | 1: `StateBundle::equations` |
| Multi-patch CFL source | configured gamma | bound EquationSet |
| Patch/bundle binding invariant | none | mandatory before timestep |
| Unsupported high-order EOS | implicit fixed-gamma use | fail fast |

The new `Model::perfectGasGamma()` is a capability query.  For high-order
WENO/TENO/Steger-Warming, the timestep check requires that capability and
requires its gamma to equal `idealGasGamma` within floating-point roundoff.
Its error identifies the selected scheme, selected EquationSet capability, and
the currently required PerfectGas thermodynamic family.

## 2. Shared RHS contract

`src/solver/algorithm/SF_densityBasedRHS.{h,cpp}` is the common sequence for
both density algorithms: stage boundary preparation; `prepareRHS`; convection
candidates for all patches; one canonical-face barrier; diffusion/sources and
coupling; then local residual staging and one GlobalDof SUM barrier.  Existing
RK dispatch, commits, and observer timing remain owned by each algorithm.

## 3. Single/multi RHS duplication

| Implementation | Before | After |
| --- | --- | --- |
| RHS sequencing/barriers | 2 copies | 1 shared implementation |
| explicit-stage publication | 2 copies | 1 shared implementation |
| EOS closure traversal | 2 copies | 1 shared implementation |

The public single and multi algorithm classes remain separate; this is not a
lifecycle merge.

## 4. Canonical-face barrier

Every patch still first produces candidate convective flux.  Exactly one
`ExecutionRuntime::finalize` then selects/copies canonical `F*` and assembles
the existing `+F*`/`-F*` terms.  The barrier stays outside the patch loop.

## 5. GlobalDof SUM barrier

Each patch still stages its local residual before one
`accumulateGlobalDof("conservativeResidual")` finalization.  These are SUM
contributions, never averages of replicated state.

## 6. Stage publication

The shared publisher retains the former write set: `conservative`, integrated
`StateBundle::transported` variables, and equation-system RHS variables.  It
deduplicates declared names, publishes after explicit integration, and does
not copy Field storage or introduce another state container.

## 7. Correction write set

`FlowAlgorithmResult` now distinguishes a correction that ran from one that
wrote conservative state.

| Result | Density runtime publication |
| --- | --- |
| no correction | none |
| correction without conservative write | none |
| correction writes conservative state | `writeOwned("conservative")` |

## 8. Correction publication

`DensityBasedRHS::publishCorrectedConservativeState()` is the one density
publication site and runs only after a confirmed conservative write.  No
temporary service override, state copy, or new service path was added.

## 9. Timestep-begin contract

The density contract checks a non-empty patch collection, a bundle EquationSet,
identity of every patch binding, PerfectGas capability for the high-order path,
and gamma agreement.  IBM/ILW single-fluid composition attaches that existing
PerfectGas binding without reinitializing conservative state.  The one-density
PerfectGas boundary dispatch reuses the established ILW kernel and preserves
its stage/order; non-PerfectGas ILW remains rejected.

## 10. Commit contract

No commit ordering changed.  Single patch keeps its existing commit/final
boundary placement; multi patch keeps its existing final boundary/commit
placement.  Diagnostics and output retain their original committed-state path.

## 11. Final closure/commit investigation

The pre-existing final-closure placement difference remains deferred to Phase
3.  The shared helper handles stage closure only and does not move a final
boundary, final closure, or commit call.

## 12. Numerical regression

The frozen data is in [baseline-phase3a-numerics.md](baseline-phase3a-numerics.md).
The final 4-rank Sod run matched its frozen final line exactly:

```text
time=2.435715e-01, dt=5.352234e-04
```

Its initial step sequence and final conservative/primitive metrics match the
recorded precision.  The serial WENO7 EquationSet Sod baseline completed after
the shared RHS/EOS implementation with frozen final `t=5.0`,
`dt=3.206358e-04`; the final subsequent ILW applicability adjustment is not
reachable in that no-ILW case.  RPI wall boiling completed successfully.

The serial ghost IBM smoke entered the new binding contract and completed
repeated RK4 closure stages without capability or closure failure before the
bounded interactive runner was stopped.  The short 4-rank IBM case completed
mesh, IBM, and canonical-face setup before the same runner limit.  These are
smoke observations, not substitutes for a full long-duration IBM regression.

## 13. Architecture check

`python3 tools/check_architecture.py` passes with **10 allowlisted dependency
edges**, unchanged.  It retains the pre-existing 17 duplicate header basename
diagnostics.  No dependency edge was added.

An obsolete CMake source entry for nonexistent
`structured/SF_dimension.h` was removed from
`src/solver/discretization/CMakeLists.txt`; this changes no executable
numerical behavior.

## 14. Deferred RK/lifecycle work

Single/multi lifecycle merge; RK/Euler/SSPRK3/RK4 unification; final
closure/commit convergence; and time-integration ownership remain deferred.

## 15. Remaining lifecycle differences

Single and multi now consume the same patch collection, EquationSet binding,
RHS orchestration, canonical-face barrier, GlobalDof barrier, and services.
Their pre-existing time dispatch and final commit/boundary placement still
differ and remain Phase 3 work.

## Deferred discretization capability

High-order WENO/TENO plus Steger-Warming currently consumes a PerfectGas-
specific thermodynamic formulation rather than generic `EquationSet` data.  It
was deliberately not replaced with `RusanovEOS::div`, which would alter spatial
discretization and dissipation.  A separate **Discretization / Thermodynamics
Interface Consolidation** phase should establish:

```text
EquationSet / ThermodynamicModel
        -> thermodynamic quantities required by numerical flux
        -> high-order reconstruction / flux splitting
```

## Validation record

| Check | Result |
| --- | --- |
| `cmake --build build --parallel 4` | passed |
| architecture checker | passed; allowlist = 10 |
| Phase 1A CTest | 6/6 passed |
| 4-rank Sod | passed; frozen final timestep matched |
| serial WENO7 EquationSet Sod | completed after shared RHS/EOS change; frozen metrics matched |
| RPI wall boiling | passed |
| serial ghost IBM | closure smoke reached; bounded before full case end |

## Explicit non-changes

- No Rusanov replacement for multi-patch WENO/TENO/Steger-Warming.
- No WENO/TENO, Steger-Warming, eigensystem, flux-splitting, or reconstruction
  formula change.
- No RK coefficient/stage-time, boundary/halo order, canonical COPY, GlobalDof
  SUM, pressure corrector, KKT, HYPRE, Field layout, or MPI ownership change.
