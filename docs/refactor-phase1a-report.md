# Phase 1A include ambiguity and architecture guard report

Date: 2026-09-11

## 1. Baseline build/test

The baseline is recorded in [refactor-baseline.md](refactor-baseline.md).
Before the refactor, the CLI/GUI build and all six CTest targets passed.  The
MPI targets require the host process runtime; a sandbox-only execution cannot
bind PRTE sockets.

## 2. Include changes

This phase changes include spelling only.  It does not move a numerical
implementation, change a call sequence, or alter a solver algorithm.

The following cross-module ambiguous includes now name their owning module:

| Call-site area | Explicit header owner selected |
| --- | --- |
| `infrastructure/mesh`, `infrastructure/mpi`, `infrastructure/io` writers | `core/field/SF_field.h` |
| `core/interfaces`, `core/state`, `solver/linearAlgebra` | `core/field/SF_field.h` and `core/interfaces/SF_executionRuntime.h` |
| `methods/math`, IBM ILW closure, boundary ILW closure | `methods/math/discrete/SF_polynomial.h` |
| application runner sources | `app/application/run/SF_output.h` or `app/application/model/SF_output.h`, according to the symbol namespace |
| mesh loader and generated-mesh writer | `infrastructure/io/mesh/SF_sfmMesh.h` and `infrastructure/io/mesh/SF_generatedMeshWriter.h` |

The runner output split was validated by compilation: `Application::Output`
functions belong to `run/SF_output.h`, while `ModelLoader::configureOutput`
belongs to `model/SF_output.h`.  The explicit paths revealed and corrected an
otherwise hidden wrong-header selection in three runner sources.

## 3. Duplicate header status

The repository still has 17 duplicate basenames.  They remain source-layout
debt, but this phase prevents new groups and reports all bare call sites:

`SF_IO`, `SF_application`, `SF_compressible`, `SF_empty`, `SF_equationCoupling`,
`SF_eulerian`, `SF_executionRuntime`, `SF_field`, `SF_fixedValue`, `SF_model`,
`SF_output`, `SF_polynomial`, `SF_state`, `SF_structured`, `SF_symmetry`,
`SF_time`, and `SF_zeroGradient`.

The baseline is stored in `tools/architecture_allowlist.json`.  The guard fails
if a new duplicate basename appears or an existing group gains another path.
Existing bare includes are printed as debt rather than treated as a compiler
resolution failure, because their actual resolution still depends on each
target's include search order.

No rename is recommended in Phase 1A.  Qualifying cross-module includes gives
an immediately reviewable safety improvement without expanding the scope into
class or API renaming.

## 4. CMake include exposure changes

`src/infrastructure/io/CMakeLists.txt` now exposes only its own directory as
`SF_io` public include state.  Its implementation-only directories
`io/mesh`, `infrastructure/mesh`, and `infrastructure/mesh/createMesh` are
private.  This reduces the chance that an unrelated consumer resolves a bare
header through `SF_io`'s transitive include search path.

No global CMake include directory was added.  `SF_headers` still exports the
source root, which is why the guard continues to report duplicate-basename
debt rather than declaring it resolved.

## 5. `SF_caseConfig` / `RunControl` dependency fix

**Previous dependency:**

```text
infrastructure/io/SF_caseConfig.h
  -> solver/algorithm/time/SF_time.h
  -> Time::RunControl
```

This made infrastructure configuration depend on an algorithm-layer time
driver, contrary to the `infrastructure -> solver` prohibition in
`ARCHITECTURE.md`.

`Time::RunControl` is a plain configuration value, not an algorithm.  It now
lives in `src/core/config/SF_runControl.h`; both the time driver and
`SF_caseConfig.h` include that core header.  The `SF::Time::RunControl` name,
six members, and all default values are unchanged.  `SF_timeDriver` declares
its public dependency on `SF_solverConfig` in CMake.

**Numerical-result impact:** none expected.  This is a declaration relocation
and an include-graph correction; the run-control data, time-driver code, and
runtime call sites are unchanged.

## 6. Architecture guard design

`tools/check_architecture.py` scans quoted includes in C/C++ headers and source
files under `src/`.  It resolves explicit paths, same-directory includes, and
unique basename includes, then enforces these direction rules:

- `core` cannot include `solver` or `models`;
- `methods` cannot include `solver`;
- `infrastructure` cannot include `solver`;
- `solver/linearAlgebra` cannot include `equation` or `algorithm`;
- `solver/discretization` and `solver/boundary` cannot include `equation` or
  `algorithm`;
- `models` cannot include `solver/boundary` or `solver/discretization` unless
  the exact edge is a recorded baseline exception.

Exact existing exceptions and duplicate-header groups are versioned in
`tools/architecture_allowlist.json`.  The checker exits nonzero for a new
forbidden edge, a new duplicate basename, or an added path in a duplicate
group.  It does not use a broad module-level exemption.

Run it with:

```sh
python3 tools/check_architecture.py
```

After Phase 1A it scanned 614 source files, found 23 explicitly allowlisted
existing dependency edges, and found no new violations.

## 7. Existing dependency allowlist

The 23 exact entries are retained as visible debt, grouped here by reason:

| Existing direction | Exact source areas | Deferred ownership issue |
| --- | --- | --- |
| `methods -> solver` | `SF_iteration.h`, `SF_scalarTransport.h` (two boundary edges), `SF_viscous.h` | CFD discretization helpers remain under `methods/numerics` |
| `models -> solver/discretization` | level-set HJ-WENO, IBM ghost weight/ILW closure, lift, turbulent dispersion, turbulence equation helpers and k-omega SST | Model code still consumes solver structured dimension/vector-calculus operators |
| `models -> solver/boundary` | phase boundary, wall heat, RPI wall mapping, wall lubrication, turbulence equation helpers | Model closures still consume boundary geometry |
| `infrastructure -> solver/linearAlgebra` | MPI HYPRE backend and Schur preconditioner | HYPRE backend ownership remains mixed with infrastructure |

Each record contains the exact source path, target header, rule, and reason;
see `tools/architecture_allowlist.json` for the authoritative list.

## 8. New dependency violations: before/after

| Measure | Before | After |
| --- | ---: | ---: |
| New forbidden dependencies introduced by this phase | 0 | 0 |
| Existing direct reverse/debt edges observed in the architecture audit | 24 | 23 |
| `infrastructure/io/SF_caseConfig.h -> solver/algorithm/time/SF_time.h` | present | removed |
| New duplicate-basename groups | 0 | 0 |

The 24-to-23 reduction is the `RunControl` dependency correction.  The other
23 pre-existing edges are guarded as exact allowlist entries, not normalized or
hidden by broad exceptions.

## 9. Tests before/after

| Check | Before | After |
| --- | --- | --- |
| `cmake --build build --parallel 4` | pass | pass |
| CTest build with `BUILD_TESTS=ON` | pass | pass |
| CTest, six targets including MPI/HYPRE paths | 6/6 pass in host MPI runtime | 6/6 pass in host MPI runtime |
| `mpirun -np 4 ./build/sonicSolver run test/sodCase` | pass | pass |
| `python3 tools/check_architecture.py` | not yet installed | pass; 23 allowlisted, no new violation |

Supplementary post-change serial IBM one-step smoke coverage used
`./build/sonicSolver --steps 1 CASE` for every non-MPI `cylinderFlow*` case:

- Passed: DFM explicit self-propelled, DFM fractional-step self-propelled,
  fictitious domain, ghost, Peskin, velocity forcing, and velocity-forcing BP.
- Reached the existing monolithic KKT feature stop: DFM augmented Lagrangian,
  DFM implicit self-propelled, and fully implicit DLM.  Each completed its RK4
  stages and then reported that the `ConstraintGlobalDof -> owned HYPRE row ->
  lambda COPY` path is not enabled.  No IBM algorithm source was changed in
  this phase, so these post-change-only results are recorded as a known
  limitation rather than attributed to the include refactor.

The Sod case completed all 500 configured steps and emitted outputs through
step 500.  A pre-change result hash was not captured, so this validates runtime
completion and unchanged configuration rather than bitwise equivalence.

## 10. Deferred issues

These were intentionally not changed:

- ownership of CFD numerical helpers presently under `methods/numerics`;
- `models -> solver/discretization` and `models -> solver/boundary` edges;
- HYPRE/Schur backend ownership between `infrastructure` and `linearAlgebra`;
- `SolverState`, `SolverServices`, and any duplicate state/service entry;
- single-patch versus multi-patch lifecycle consolidation;
- the scope and representations owned by `Field`;
- IBM architecture and its distributed monolithic KKT limitation;
- Equation DSL / AssemblyPlan evolution.

No `Manager`, `Adapter`, `Service`, `Context`, `Registry`, or other new
architecture abstraction was introduced.  No public runtime API, solver
lifecycle, MPI ownership rule, numerical operator, field representation, or
IBM method implementation changed.  Apart from the moved plain configuration
declaration, the source changes are include qualification and CMake include
visibility only.
