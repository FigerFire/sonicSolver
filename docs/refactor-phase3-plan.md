# Phase 3 Density-Based Lifecycle Plan

## Preconditions

The checked source has one `StateBundle` state/clock, `patches` is the only
patch enumeration, services bind once, and `advance(SolverState&)` is the only
public timestep entry.  Baseline architecture debt is 10 allowlisted edges.

## Actual stage comparison

| Stage | CompressibleAlgorithm | MultiPatchAlgorithm | Same math? | Preserve |
| --- | --- | --- | --- | --- |
| validate state | requires exactly one patch | validates nonempty patch collection | yes | fail-fast state checks |
| boundary/halo/ghost | pipeline, otherwise BC -> halo -> ghost IBM -> halo | same pipeline/manual sequence | yes | exact order |
| CFL | EquationSet-aware Rusanov when bound | per patch gamma `deltaT` | no | multi must use bundle EquationSet or fail explicitly |
| global dt | runtime min then driver/config cap | same | yes | reduction before cap |
| begin step | equation coupling once for field | coupling once for fields | yes | exactly once |
| stage boundary | before every RHS | before every RHS | yes | stage time |
| convection | one patch then canonical finalize | all patches then one finalize | mathematically same | all-patch barrier |
| diffusion/source/coupling | after finalize | after finalize | yes | barrier order |
| GlobalDof residual | finalize after RHS | finalize after RHS | yes | SUM not COPY |
| RK update | `ddtDispatch` Euler/SSPRK3/RK4 | `SF_multiPatchTime.cpp` | duplicated | select multi-patch stage loop |
| transported variables | bundle registry passed to dispatch | per-patch equation registry | same intent | update alongside Q stage |
| publish/validate | WriteOwned conservative/registered then validate | same, deduplicated registrations | same | current freshness contracts |
| correction | may WriteOwned after forcing correction | no corresponding publish branch | differs | preserve only correction write-set |
| commit/final closure | commit, then closure except forcing branch | closure, then commit | differs | do not merge until contract proven |
| clock/observer | commit clock then observer | same | yes | observer after successful commit |

## Authoritative implementation

`MultiPatchAlgorithm` and `SF_multiPatchTime.cpp` are the implementation to
retain because their stage loop already enforces all-patch convection,
canonical-face finalization, remaining RHS, GlobalDof accumulation, and update
barriers.  `CompressibleAlgorithm` will be migrated away after its EquationSet
aware CFL/flux/validation and forcing-correction write contract are represented
there without changing coefficients or call timing.

## Required decisions before merge

1. The multi-patch path must pass `StateBundle::equations` to CFL, convection,
   diffusion/source, and closure validation.  If a participating patch lacks
   the same binding, fail before the timestep; do not use gamma fallback.
2. `flowAlgorithm.correct()` must expose or retain its existing write-set.
   Only a correction that writes conservative owner state may issue the existing
   WriteOwned publication.  No additional halo or closure may be added merely
   to make paths look alike.
3. `commitStep` must stay after correction and occur once.  The final boundary
   order is not merged until the coupling contract proves whether it needs
   pre-commit or post-commit state.  Add a focused invariant/regression first
   if that proof is absent.

## Target lifecycle

```
validate bundle -> compute every patch CFL -> runtime global min -> beginStep
for each existing RK stage:
  boundary/halo/IBM at stage time
  all-patch convection -> canonical-face finalize
  all-patch diffusion/source/coupling -> GlobalDof SUM
  existing stage update -> publish -> validate
flow/IBM correction -> publish only its actual write-set
commitStep -> existing required final closure -> commit bundle clock -> observer
```

No MPI conditional, new manager/context/adapter, pressure-based modification,
or change to RK coefficients, EOS/flux/viscous/source/IBM formula is allowed.
