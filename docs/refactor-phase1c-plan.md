# Phase 1C model, equation, and discretization responsibility plan

## Baseline

`python3 tools/check_architecture.py` reports 23 allowlisted dependency
entries: four `methods -> solver`, seventeen `models -> solver`, and two
`infrastructure -> solver`. Before-source baselines passed the Debug build,
six CTest targets, and the four-rank Sod case. The serial IBM smoke baseline
has seven passing methods and three pre-existing distributed-KKT stops.

## Dependency classification before code changes

| Current source | Current target | Type | Current use | Phase 1C decision |
| --- | --- | --- | --- | --- |
| `models/turbulence/SF_turbulence.cpp` | boundary geometry | B | Structured patch axis for scalar turbulence BC. | Move neutral topology lookup to `core/mesh`; preserve BC timing. |
| `models/turbulence/SF_turbulence.cpp` | dimension | D | Active-space metadata. | Move metadata to `core/mesh`. |
| `models/turbulence/SF_equationModelOps.cpp` | boundary geometry | B | Physical-point and active-patch lookup. | Move neutral geometry functions to `core/mesh`. |
| `models/turbulence/SF_equationModelOps.cpp` | dimension | D | Filter-width active dimensions. | Move metadata to `core/mesh`. |
| `models/turbulence/SF_equationModelOps.cpp` | vector calculus | A | Field-aware velocity gradient for strain and vorticity. | Defer: timing and halo/BC readiness belong to the equation/discretization caller. |
| `models/turbulence/SF_equationSystem.cpp` | boundary geometry | B | Patch/physical-cell topology for transported turbulence variables. | Move neutral geometry functions to `core/mesh`. |
| `models/turbulence/SF_kOmegaSSTEquation.cpp` | vector calculus | A | Field-aware gradients of k and omega. | Defer: changing derivative acquisition would change model evaluation timing. |
| `models/turbulence/RAS/SF_kOmegaSST.cpp` | vector calculus | A | Field-aware gradients and vorticity. | Defer with the legacy turbulence path. |
| `models/physics/phaseSystem/SF_phaseBoundary.cpp` | boundary geometry | C | Phase model currently applies PDE boundary values and ghosts. | Move only neutral topology queries; defer equation/boundary ownership. |
| `models/physics/heat/SF_wallHeatSource.h` | boundary geometry | B | Wall patch axis and mesh coordinates. | Move neutral geometry functions to `core/mesh`. |
| `models/physics/phaseChange/RPI/SF_wallMapping.cpp` | boundary geometry | B | Wall-cell/patch topology. | Move neutral geometry functions to `core/mesh`. |
| `models/physics/phaseSystem/interphase/SF_wallLubrication.cpp` | boundary geometry | B | Wall points and patch normal. | Move neutral geometry functions to `core/mesh`. |
| `models/physics/interfaceModel/levelSet/SF_hjWeno.cpp` | dimension | C/D | Active-direction metadata used by model-owned HJ-WENO. | Move metadata only; retain HJ-WENO and record its PDE ownership as deferred. |
| `models/ibm/method/ghost/SF_ilwClosure.cpp` | dimension | D | 2D projection metadata. | Move metadata only; retain ILW algorithm. |
| `models/ibm/method/ghost/SF_weightBuilder.h` | dimension | D | 2D projected IBM geometry. | Move metadata only; retain IBM classification and weights. |
| `models/physics/phaseSystem/interphase/SF_lift.cpp` | vector calculus | A | Field-aware curl of continuous-phase velocity. | Defer: caller must supply same-stage curl before closure evaluation. |
| `models/physics/phaseSystem/interphase/SF_turbulentDispersion.cpp` | vector calculus | A | Field-aware volume-fraction gradient. | Defer for the same timing reason. |

The four `methods -> solver` and two `infrastructure -> solver` entries are
outside this phase. `SF_scalarTransport.h` is explicitly frozen. The
`SF_viscous.h` move is reconsidered only after the turbulence derivative
dependencies are understood; no move is planned because legacy turbulence
still invokes it.

## Safe ownership change

`SF_dimension.h` contains process-wide active-direction and inactive-plane
metadata, with no stencil or derivative operation. It will move unchanged to
`core/mesh/SF_dimension.h`.

`SF_boundaryGeometry.h` contains two responsibilities. The Field-based
structured topology and coordinate helpers (`BoundarySetInfo`, axis/point
queries, physical-point tests, normal estimation, and set analysis) will move
unchanged to `core/mesh/SF_meshBoundaryGeometry.h`. EMPTY configuration remains
in the solver boundary header because it consumes BC settings and changes the
active-space state. Model call sites will use the neutral mesh header.

No vector-calculus function moves: every function in `SF_vectorCalculus.h`
either reads `Field` metrics or samples Field data. Tensor algebra remains in
`methods/math/SF_tensor.h`.

## Invariants for implementation

- Copy functions and state operations without changing bodies, formulae, or
  failure behavior.
- Retain all callers and invocation order; only their include/namespace path
  changes for neutral metadata and geometry.
- Do not change time integration, scalar transport, turbulence source timing,
  halo synchronization, IBM ILW, RPI/heat formulae, MPI ownership, or solver
  lifecycle.
- Do not introduce a Context, Adapter, Manager, Service, or compatibility
  forwarding class.
