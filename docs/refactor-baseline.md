# Phase 1A refactor baseline

Date: 2026-09-11

## Build configuration

The existing `build/` tree was configured with Ninja, Debug, CLI and GUI enabled,
MPI/HYPRE enabled, and `BUILD_TESTS=OFF`.  It uses AppleClang 21 and Open MPI
5.0.9.  OpenMP was not found, so this configuration is single-core apart from
MPI execution.

The source baseline built successfully before Phase 1A with:

```sh
cmake --build build --parallel 4
```

Because the existing tree did not enable CTest targets, a separate,
non-source build tree was configured before the change:

```sh
cmake -S . -B build-phase1a-tests -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug -DBUILD_CLI=ON -DBUILD_GUI=ON \
  -DBUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-phase1a-tests --parallel 4
ctest --test-dir build-phase1a-tests --output-on-failure --parallel 1
```

All six tests passed in the host MPI runtime:

| Test | Baseline result |
| --- | --- |
| `registryIO` | pass |
| `surfaceSchurProjection` | pass |
| `distributedConstraintLayout` | pass |
| `sparseCanonicalCopy` | pass |
| `distributedMinimalKKT` | pass |
| `distributedSurfaceSchur` | pass |

The same CTest invocation in the filesystem sandbox passed the three non-MPI
tests but could not start the three MPI tests because PRTE was denied socket/
binding permission.  This was an execution-environment restriction; the
identical suite passed when run with the host MPI runtime.

The serial command `./build/sonicSolver run test/sodCase` is not a valid
single-rank regression for this case: its configured partition product is four.
The established baseline used:

```sh
mpirun -np 4 ./build/sonicSolver run test/sodCase
```

It completed successfully.  No result hash was captured before the run, so
this baseline establishes completion and runtime configuration, not a bitwise
field comparison.

The CTest suite covered IBM surface Schur code before the change.  The ten
serial `cylinderFlow*` runtime smoke cases were not run before the source edit;
they are therefore reported as post-change-only supplementary coverage in the
Phase 1A report rather than presented as a before/after comparison.
