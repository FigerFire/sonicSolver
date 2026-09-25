# Phase 1C Model / Equation / Discretization responsibility report

## 1. Model Dependency Classification

The pre-change allowlist had 23 entries, including seventeen model-to-solver
edges. The classification was written before source changes in
[refactor-phase1c-plan.md](refactor-phase1c-plan.md).

| Category | Pre-change entries | Result |
| --- | ---: | --- |
| A — model consumes a Field-aware derivative | 5 | Deferred. The model still invokes a discrete gradient or curl, so moving it would require the equation/discretization caller to preserve stage, boundary, and halo timing. |
| B — model consumes boundary geometry/topology | 7 | Resolved. The queries now come from core/mesh/SF_meshBoundaryGeometry.h. |
| C — model also owns PDE/boundary behavior | 2 | Phase boundary remains a deferred ownership issue; level-set HJ-WENO remains a deferred model-owned scheme. Its dimension dependency was removed. |
| D — dimension/geometry metadata | 5 | Resolved by moving SF_dimension.h to core/mesh. |

The two infrastructure and four methods baseline entries were not redesigned.
One methods entry naturally disappeared because SF_iteration.h now consumes
neutral dimension metadata.

## 2. SF_dimension Responsibility Before/After

Before, src/solver/discretization/structured/SF_dimension.h held only
active-direction state, inactive-plane normal metadata, dimension count, and
text formatting. It contains no stencil, derivative, face assembly, or
discretization operator.

It now lives at src/core/mesh/SF_dimension.h, unchanged in behavior. All
callers include that explicit path. This removes direct solver/discretization
dependencies from level-set HJ-WENO, IBM ghost ILW/weight construction,
turbulence metadata users, and methods/numerics/structured/SF_iteration.h.

EMPTY parsing and the decision to update this metadata remain in
src/solver/boundary/SF_boundaryGeometry.h, because that action consumes BC
settings and is boundary configuration. No active-direction value, reset point,
or update order changed.

## 3. SF_vectorCalculus Responsibility Before/After

src/solver/discretization/structured/SF_vectorCalculus.h remains in
discretization. Its gradSampled, scalar/vector grad, curl, metric conversion,
and offsets read Field metrics or Field samples. They are Field-aware discrete
derivatives, not pure math.

Pure tensor operations used by these call sites (symm, curl, and doubleDot
after a gradient is available) already remain in methods/math/SF_tensor.h. No
vector-calculus code was copied into models or moved to methods.

## 4. Boundary Geometry Responsibility Before/After

Before, src/solver/boundary/SF_boundaryGeometry.h mixed two concerns:

- neutral structured mesh/topology queries: physical point, set axis, point
  coordinates, normal estimate, and boundary-set analysis;
- EMPTY BC configuration that resets/deactivates active directions.

The neutral queries now have one implementation in
src/core/mesh/SF_meshBoundaryGeometry.h, under
SF::StructuredMesh::BoundaryGeometry. Model callers now use that header
directly. src/solver/boundary/SF_boundaryGeometry.h retains EMPTY
configuration and imports the neutral functions for its existing boundary and
frozen scalar users; no formula or error path is duplicated.

This is a header-only move inside the source root exported by SF_headers;
there was no compiled source or CMake target ownership to change.

## 5. Turbulence Changes

The following turbulence files now consume neutral geometry/metadata rather
than solver boundary or discretization headers:

- src/models/turbulence/SF_turbulence.cpp
- src/models/turbulence/SF_equationModelOps.cpp
- src/models/turbulence/SF_equationSystem.cpp

Their boundary application, wall-distance construction, and model/source call
order are unchanged. SF_equationModelOps.cpp, SF_kOmegaSSTEquation.cpp, and
RAS/SF_kOmegaSST.cpp still directly use SF_vectorCalculus.h for same-call
Field-aware derivatives. Moving those derivatives into an equation caller
would require an explicit proof that the gradient is computed after the same
halo and boundary updates and against the same stage state. That proof is not
available in this phase, so the paths remain unchanged.

## 6. Phase/RPI/Wall Model Changes

SF_phaseBoundary.cpp, RPI/SF_wallMapping.cpp,
interphase/SF_wallLubrication.cpp, and heat/SF_wallHeatSource.h now use
neutral mesh topology/geometry queries. Wall normal estimation, wall-cell
mapping, source formulas, phase boundary ghost writes, and their invocation
order are unchanged.

SF_phaseBoundary.cpp still implements phase-field boundary behavior. It is a
model-owned PDE/boundary path and is deferred to equation/boundary ownership
work. src/models/physics/multiphase/SF_multiphase.cpp, advanceTemperature,
and advanceForwardEuler were not modified.

## 7. Level-Set Changes

src/models/physics/interfaceModel/levelSet/SF_hjWeno.cpp now includes
core/mesh/SF_dimension.h for active-direction metadata. Its HJ-WENO
derivative, upwinding, reconstruction, Field sampling, and timing were not
changed.

HJ-WENO remains a model-owned PDE-specific numerical scheme. Its eventual
equation/discretization ownership decision is deferred.

## 8. IBM Ghost Dependency Changes

src/models/ibm/method/ghost/SF_ilwClosure.cpp and SF_weightBuilder.h now
consume core/mesh/SF_dimension.h. The moved API is only the 2D active-plane
metadata read by projected geometry and ILW local frames.

No ILW polynomial, sample selection, classification, donor weights, IBM
runtime state, or ghost reconstruction behavior changed.

## 9. SF_viscous Ownership Re-evaluation

src/methods/numerics/viscous/SF_viscous.h remains in place and remains
allowlisted for its dependency on discrete vector calculus. The legacy
RAS/SF_kOmegaSST.cpp and turbulence code still invoke its velocityGradientAt
path. Moving it now would either keep a forwarding path or change derivative
ownership/timing, so neither is safe in this phase.

## 10. Architecture Allowlist

23 -> 10

Removed entries:

- one methods -> solver/discretization entry for SF_iteration.h;
- five models -> solver/discretization dimension-metadata entries;
- seven models -> solver/boundary geometry/topology entries.

| Remaining source | Target | Why retained | Future decision | Could change numerical results? |
| --- | --- | --- | --- | --- |
| models/turbulence/SF_equationModelOps.cpp | solver/discretization/structured/SF_vectorCalculus.h | Same-call velocity gradient for closure quantities. | Equation/discretization must provide the same-stage gradient. | Yes, if stage/halo/BC timing changes. |
| models/turbulence/SF_kOmegaSSTEquation.cpp | solver/discretization/structured/SF_vectorCalculus.h | k/omega gradients feed SST blending. | Provide derivative inputs through the turbulence equation path. | Yes. |
| models/turbulence/RAS/SF_kOmegaSST.cpp | solver/discretization/structured/SF_vectorCalculus.h | Legacy scalar gradients and vorticity. | Consolidate legacy and registered turbulence equation paths. | Yes. |
| models/physics/phaseSystem/interphase/SF_lift.cpp | solver/discretization/structured/SF_vectorCalculus.h | Continuous-phase curl for lift closure. | Equation computes/provides same-state curl. | Yes. |
| models/physics/phaseSystem/interphase/SF_turbulentDispersion.cpp | solver/discretization/structured/SF_vectorCalculus.h | Volume-fraction gradient for dispersion closure. | Equation computes/provides same-state gradient. | Yes. |
| methods/numerics/scalar/SF_scalarTransport.h | solver/boundary/SF_boundaryGeometry.h | Frozen model-owned scalar PDE boundary path. | Equation/lifecycle consolidation. | Yes. |
| methods/numerics/scalar/SF_scalarTransport.h | solver/boundary/reconstruction/ILW/SF_boundaryClosure.h | Frozen scalar ILW closure path. | Equation/lifecycle consolidation. | Yes. |
| methods/numerics/viscous/SF_viscous.h | solver/discretization/structured/SF_vectorCalculus.h | Existing Field-aware viscous derivative helper. | Migrate callers after gradient timing is explicit. | Yes. |
| infrastructure/mpi/backend/SF_hypreBackend.cpp | solver/linearAlgebra/hypre/SF_hypre.h | Existing HYPRE backend coupling. | Infrastructure/linear-algebra ownership phase. | Potentially, if distributed assembly changes. |
| infrastructure/mpi/backend/SF_schurPreconditioner.h | solver/linearAlgebra/SF_blockSchur.h | Existing Schur backend coupling. | Infrastructure/linear-algebra ownership phase. | Potentially. |

## 11. Remaining models -> solver edges

The five remaining model edges are the first five rows in section 10. They
remain only for Field-aware discrete derivatives. No derivative was moved,
cached, recomputed, or passed through a new context. This keeps the original
evaluation timing and numerical results.

## 12. Remaining methods -> solver edges

The three remaining methods edges are the scalar transport and viscous rows in
section 10. SF_scalarTransport.h is explicitly frozen because it is coupled
to model-owned scalar time integration. SF_viscous.h remains because its
turbulence callers still need the original direct derivative path.

## 13. Numerical Timing Invariants

- Active-direction state is reset and configured by the same boundary setup
  calls as before.
- Geometry/axis/normal query function bodies and failure behavior are
  unchanged; only their owner and qualified call path changed.
- Gradient, curl, halo synchronization, boundary application, wall source
  evaluation, turbulence correction, and equation assembly order are
  unchanged.
- No time scheme, phase temperature advance, pressure corrector, multi-patch
  lifecycle, MPI COPY/SUM/canonical ownership, or numerical formula changed.

## 14. Build/Test Before/After

Both baseline and post-change checks passed:

| Check | Before | After |
| --- | --- | --- |
| python3 tools/check_architecture.py | 23 allowlisted entries | 10 allowlisted entries, no new violation |
| cmake --build build --parallel 4 | Passed | Passed |
| ctest --test-dir build-phase1a-tests --output-on-failure --parallel 1 | 6/6 passed | 6/6 passed |
| mpirun -np 4 ./build/sonicSolver run test/sodCase | Passed | Passed |
| ./build/sonicSolver --steps 1 test/rpiWallBoilingCase | Not required baseline case | Passed; WallHeat, RPI, k-epsilon, and phase pressure coupling ran |
| Serial test/IBM/cylinderFlow methods | 7 passed; 3 KKT stops | 7 passed; the same 3 KKT stops |

The three serial IBM stops are cylinderFlowDFMAugmentedLagrangian,
cylinderFlowDFMImplicitSelfPropelled, and cylinderFlowFullyImplicitDLM. Each
stops because distributed monolithic IBM KKT is not enabled until its owned
HYPRE row/lambda path is connected. This was present before the Phase 1C
changes.

## 15. Deferred Changes

- Provide Field-aware gradients/curl from equation/discretization to
  turbulence, lift, and turbulent-dispersion closure evaluation while proving
  identical stage, halo, and boundary timing.
- Define equation ownership for phase boundary behavior and model-owned
  HJ-WENO.
- Consolidate SF_scalarTransport.h with equation/lifecycle ownership.
- Revisit SF_viscous.h only after the legacy turbulence derivative callers
  have a safe replacement.
- Keep HYPRE/Schur backend ownership for its dedicated infrastructure phase.

## Self Review

- Solver numerical code copied into a model: no.
- Field-aware derivative moved into methods: no.
- Context, Adapter, Manager, Service, or compatibility class added: no.
- Gradient, boundary, halo, or model-source timing changed: no.
- Numerical formulas, MPI ownership, scalar time integration, and
  single/multi-patch lifecycle changed: no.

