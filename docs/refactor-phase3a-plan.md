# Phase 3A: Density-Based RHS / EOS / Stage Synchronization Plan

## Confirmed thermodynamic/discretization decision

`StateBundle::equations` is the authoritative density-based thermodynamic
binding. CFL, state closure, thermodynamic validation, and capability
validation use that binding. Multi-patch WENO/TENO with Steger-Warming keeps
its existing PerfectGas-specific implementation; it is not replaced with
`RusanovEOS::div`. Before a density timestep, every participating field must
reference the bundle binding and the high-order path must prove that it is a
PerfectGas binding with the configured gamma. Unsupported bindings fail before
the first stage. This phase does not add generic-EOS WENO/TENO,
characteristic, flux-splitting, or reconstruction support.

## Scope and preserved lifecycle

This phase changes only the shared density-based RHS, EOS, and state-publication
contract. `CompressibleAlgorithm`, `MultiPatchAlgorithm`, their existing RK
implementations, `SF_multiPatchTime.cpp`, and the existing final-boundary /
`commitStep` ordering remain in place. No flux formula, stage coefficient or
stage time, boundary/halo/IBM ordering, canonical-face ownership, or GlobalDof
SUM operation changes.

## Baseline findings

`StateBundle::equations` is already the single-patch EquationSet binding. The
multi-patch application does not populate it, and multi-patch CFL, convection,
and closure instead use the fixed-gamma path. Its RHS context passes null
thermodynamics. This creates two different closures for the same density-based
equation family.

| Contract | Single before | Multi before | Phase 3A action |
| --- | --- | --- | --- |
| EOS binding | bundle EquationSet when present | no bundle EquationSet | bind one perfect-gas EquationSet for multi-patch single-fluid cases; reject an unbound density state |
| CFL | `RusanovEOS::deltaT` | fixed-gamma `deltaT` | every patch uses the bundle EquationSet |
| convection | EquationSet in `AssemblyContext` | null context member | pass the same bundle EquationSet |
| state closure | EquationSet-aware | fixed-gamma physical validation | use EquationSet-aware validation |
| RHS barriers | one-patch instance | all-patch instance | one shared vector-of-patches implementation, retaining the all-patch barrier |

## Timestep-start binding invariant

Density-based algorithms will validate before CFL that the bundle has an
EquationSet and every participating patch points to that same object. The
diagnostic will include the patch index and expected/actual binding addresses.
It will not alter a Field binding. Pressure-based workflows are outside this
contract.

## Shared RHS implementation shape

Use ordinary internal functions in `solver/algorithm`, with explicit arguments
instead of a new context, manager, adapter, or interface. Both algorithms will
pass a vector of participating patches; a single-patch vector has size one.
The shared sequence is fixed as:

1. boundary/halo/IBM preparation at the existing stage time;
2. `prepareRHS` for every patch;
3. `begin` and convection candidate computation for every patch;
4. exactly one canonical-face finalization barrier;
5. diffusion/source and equation coupling for every patch; and
6. local residual staging for every patch followed by exactly one GlobalDof
   `SUM` barrier.

The helper selects the pre-existing service transport model first and the
per-patch equation-coupling transport model only when the former is absent.
That preserves the prior single and multi service selections.

## Correction write set

The current `FlowAlgorithmResult` reports only `performedCorrection`, so it
cannot state whether conservative owner state was written. Its actual writers
are:

| Algorithm / correction | Conservative write | Existing publication | Required Phase 3A change |
| --- | --- | --- | --- |
| density based, no IBM constraint | no | none | none |
| density based, constraint projection | yes | single only, conditional on constraint | report conservative write and publish in either density path |
| pressure corrector | yes | no density-path publication rule | leave pressure lifecycle untouched |

Add the smallest result extension, `wroteConservativeState`, and set it in the
existing correction implementations. Density-based algorithms will publish the
existing `WriteOwned("conservative")` contract only when that flag is true.
No halo exchange is added after the publication: the next existing boundary
preparation remains the reader.

## Commit / final closure investigation

`CompressibleAlgorithm` currently performs `commitStep` then final closure,
except after immersed forcing where it defers the next read to the next step.
`MultiPatchAlgorithm` performs final closure then `commitStep`. The coupling
interfaces receive the same one `beginStep`, per-stage `prepareRHS` /
`assembleRHS`, and one post-correction `commitStep`; their required ordering is
not proven identical. Phase 3A preserves both orderings and only makes
correction publication explicit.

## Validation

The frozen baseline is recorded in `docs/baseline-phase3a-numerics.md`. After
each source-edit group: build and run focused density tests. Final validation
will repeat architecture checking, build, CTest, four-rank Sod, serial
density-based Sod, RPI wall boiling, and serial IBM smoke as applicable.
