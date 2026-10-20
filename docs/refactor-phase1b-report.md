# Phase 1B CFD numerics ownership report

Date: 2026-09-11

## Status

Phase 1B made the two ownership moves that preserve all stated dependency
rules.  It does **not** satisfy the requested end-state of moving every CFD
spatial implementation out of `methods/numerics`: the remaining candidates
are directly consumed by `models` or `infrastructure`, and moving them would
create new forbidden direct dependencies.  This report records those concrete
blockers instead of adding allowlist entries or forwarding wrappers.

## 1. Numerics Ownership Before

`src/methods/numerics` contained 4,910 lines across Field-independent kernels,
Field-aware face/structured/viscous/scalar operators, state validation, IBM
kernels, and time integration.  The Phase 1A architecture check had 23 exact
allowlist entries, including four `methods -> solver` entries:

1. `structured/SF_iteration.h -> solver/discretization/structured/SF_dimension.h`;
2. `viscous/SF_viscous.h -> solver/discretization/structured/SF_vectorCalculus.h`;
3. `scalar/SF_scalarTransport.h -> solver/boundary/SF_boundaryGeometry.h`;
4. `scalar/SF_scalarTransport.h -> solver/boundary/reconstruction/ILW/SF_boundaryClosure.h`.

## 2. File Classification

The complete pre-change classification is in
[refactor-phase1b-plan.md](refactor-phase1b-plan.md).  It covers every source
and header under `methods/numerics`, including the required convection, Flux,
structured, viscous, scalar, and time areas.

The essential result is:

- scalar WENO/TENO kernels, Riemann array kernels, Newtonian tensor algebra,
  and generic immersed row kernels are Field-independent;
- WENO field/stencil/IBM assembly, numerical face fluxes, structured traversal,
  viscous RHS assembly, scalar RHS/boundary handling, and Field validation are
  CFD discretization or CFD state logic;
- Euler/RK4/time helpers are lifecycle logic and remain deferred.

## 3. Files Moved

| Before | After | Change |
| --- | --- | --- |
| `methods/numerics/convection/SF_WENO.h` | `solver/discretization/convection/SF_WENO.h` | Field-aware WENO direction/stencil/IBM/face-flux assembly. |
| `methods/numerics/viscous/SF_newtonian.h` | `methods/math/SF_newtonian.h` | Stateless Newtonian stress and traction algebra. |

`SF_div.h` now includes the WENO assembly header through its explicit solver
path.  `SF_numerics.h` no longer presents that solver-owned assembly entry.
There is no forwarding header at the old WENO path.

## 4. Pure Math Kernels Retained in methods

The WENO3, WENO5, WENO7, and TENO5 implementations remain in
`methods/numerics/convection`: each consumes a fixed scalar stencil and emits a
scalar reconstruction using unchanged candidate polynomials, smoothness
indicators, epsilon, and nonlinear weights.  `SF_riemann.{h,cpp}` remains there
as the array-only Euler eigensystem/matrix kernel.  The pure Newtonian tensor
kernel now lives under `methods/math`.

Generic immersed interpolation, spreading, regularized-kernel, and local
constraint algebra headers also remain in `methods`; IBM ownership itself was
not changed.

## 5. CFD Discretization Moved to solver

`solver/discretization/convection/SF_WENO.h` now owns the code that:

- traverses Field faces and selects active directions;
- acquires conservative stencils and face metrics;
- evaluates IBM/non-fluid and ILW-closed stencil conditions;
- selects the configured low-order IBM method when applicable;
- writes canonical face fluxes.

It depends downward on retained numerical kernels in `methods`, which conforms
to `solver/discretization -> methods`.  It does not create an upward solver
dependency.

## 6. Deferred Time Integration

`methods/numerics/time/SF_Euler.*`, `SF_RK4.*`, and `SF_time.h` are unchanged.
They retain Euler updates, RK stage handling, callbacks, and state commits.
`solver/algorithm/SF_multiPatchTime.cpp` is also unchanged.

Time integration debt remains intentionally deferred.

## 7. CMake Dependency Before/After

Before, `SF_discretization` consumed WENO assembly indirectly through the
methods/numerics public include set.  After, its `target_sources` explicitly
lists `convection/SF_WENO.h`; the header itself uses explicit paths for core
data and retained methods kernels.  The existing target direction remains:

```text
SF_discretization -> SF_riemann / SF_flux / SF_weno* -> core
```

No `methods` target was linked to a solver target.  No global include directory
was added.

## 8. Architecture Allowlist Before/After

```text
23 -> 23
```

No entry was added, broadened, or removed.  The WENO assembly header was not
one of the four allowlisted `methods -> solver` edges, so its ownership move
does not change that count.

The count cannot safely shrink within this phase because moving `SF_iteration`,
`SF_viscous`, or `SF_scalarTransport` would create direct dependencies that the
same guard correctly rejects:

| Candidate move | Existing external consumers that would become forbidden direct consumers |
| --- | --- |
| structured traversal/canonical-face headers | `models/physics`, `models/turbulence`, `models/ibm`, `infrastructure/mesh`, and `infrastructure/mpi` |
| `SF_viscous.h` | `models/turbulence/SF_turbulence.cpp` and `models/turbulence/RAS/SF_kOmegaSST.cpp` |
| `SF_scalarTransport.h` | `models/physics/multiphase/SF_multiphase.cpp` |
| Flux headers | `infrastructure/mpi/SF_haloExchange.cpp` uses `SF_lxF.h` and canonical-face utilities |

Replacing those calls or moving their ownership is a model/infrastructure/IBM
architecture change, prohibited by the Phase 1B scope.  Allowlisting the new
direct edges would hide the issue; a compatibility header at the old methods
path would preserve the reverse dependency.  Neither was done.

## 9. Public API Before/After

No namespace, function signature, class, or runtime API changed.  The only
include-path changes are the intentional ownership moves listed above:

- consumers of Field-aware WENO assembly include
  `solver/discretization/convection/SF_WENO.h`;
- consumers of pure Newtonian tensor algebra include
  `methods/math/SF_newtonian.h`.

The obsolete `methods/numerics/convection/SF_WENO.h` entry was removed rather
than kept as a compatibility wrapper.  No Manager, Adapter, Context, Service,
Registry, or new interface was introduced.

## 10. Numerical Invariants Checked

The move changed declarations, include paths, and CMake ownership only.  The
following WENO assembly and kernel semantics were not edited:

- WENO/TENO candidate polynomials, smoothness indicators, epsilon, nonlinear
  weights, and stencil offsets;
- face metric normalization and orientation;
- Riemann/numerical flux selection and formulae;
- IBM/ILW branch order, ghost interpretation, and fail-fast behavior;
- canonical-face owner calculation, residual sign, and MPI COPY/SUM semantics.

`SF_newtonian.h` changed only location; its stress and traction expressions are
byte-for-byte unchanged.

No separate numerical-equivalence kernel test was added.  The existing CTest,
four-rank Sod, and IBM smoke cases exercise the moved WENO public path without
requiring a new test framework.

## 11. Build/Test Before/After

| Check | Before | After |
| --- | --- | --- |
| Full Debug CLI/GUI build | pass | pass |
| `BUILD_TESTS=ON` build | pass | pass |
| CTest (including MPI/HYPRE paths) | 6/6 pass | 6/6 pass |
| `mpirun -np 4 ./build/sonicSolver run test/sodCase` | pass | pass, 500 steps |
| `python3 tools/check_architecture.py` | pass, 23 allowlisted | pass, 23 allowlisted |

Post-change serial IBM one-step smoke results match the Phase 1A baseline:
DFM explicit self-propelled, DFM fractional-step self-propelled, fictitious
domain, ghost, Peskin, velocity forcing, and velocity-forcing BP pass.  DFM
augmented Lagrangian, DFM implicit self-propelled, and fully implicit DLM stop
at the pre-existing unimplemented distributed monolithic KKT path after RK4.

## 12. Remaining methods -> solver edges

All four exact Phase 1A entries remain intentionally:

| Source | Target | Why it cannot safely be removed in Phase 1B |
| --- | --- | --- |
| `structured/SF_iteration.h` | structured dimension state | Its Field traversal is used directly by models and infrastructure. |
| `viscous/SF_viscous.h` | vector calculus | Turbulence closure code directly calls its velocity-gradient API. |
| `scalar/SF_scalarTransport.h` | boundary geometry | Multiphase model owns this scalar API and invokes its boundary work. |
| `scalar/SF_scalarTransport.h` | ILW boundary closure | The same model-owned API invokes scalar ILW closure. |

The scalar header also contains `advanceForwardEuler`; separating spatial RHS
from it without moving its callers would require a forwarding facade, while
moving callers would modify model architecture.  Both are outside this phase.

## 13. Deferred model dependencies

The nineteen model-to-solver/boundary allowlist entries remain unchanged.  They
cover turbulence, phase-system, RPI, level-set, and IBM consumers of structured
dimension/vector-calculus or boundary geometry.  They require a separate
decision about closure inputs, gradient timing, and boundary ownership.

HYPRE backend ownership, SolverState, SolverServices, Field responsibilities,
single/multi-patch lifecycle, IBM architecture, and Equation DSL/AssemblyPlan
are also unchanged.

## Self review

- Compatibility wrapper added: no.
- Adapter/Manager/Service/Context/Registry added: no.
- Time integration or stage order changed: no.
- Boundary invocation order changed: no.
- MPI COPY/SUM/canonical semantics changed: no.
- Numerical formulas changed: no.
