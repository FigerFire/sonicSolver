# Dynamic equation feedback acceptance

The user authorized these two serial groups and their metrics/trace comparisons on 2026-11-01. All four configurations completed, and the scoped checks in `metrics.py` passed. This is not a full regression/MPI/model-combination claim.

- `density/baseline`: uniform flow on the Sod-derived mesh, classical RK4, 0–1 s / 40 steps, no source.
- `density/feedback`: same setup plus Stage axial force-density and mechanical-power relations feeding actual RHS.
- `pressure/baseline`: constant-density transient PISO Poiseuille, 0–10 s / 2048 steps, original pressure drop and no-slip plate boundaries.
- `pressure/feedback`: same setup plus a Step Current-velocity relation feeding predictor acceleration.

Registered SFM input meshes are versioned. Current runtime files use `createMesh:false`: true is a mesh-generation-only mode that exits without CFD. `prepare.py` is a historical input/bootstrap generator and never invokes a solver; it refuses existing cases. Its freshly created inputs need one explicit mesh-generation pass, then `createMesh:false` before physical execution. Do not overwrite these frozen inputs or existing results.

All YAML follows production native IO: one registry per object directory and solver documents directly in `solvers/`. Output is under each case's `result/`, with final trace runs in `result/trace-04/`. `run-01.log` preserves mesh-only/entry failures; `run-02.log` is the first complete CFD run; `run-03.log` preserves the first traced run; `run-04.log` is the final traced run after explicitly restricting the new verification snapshots to serial execution. Local `evidence/` holds authorization, exact configurations, binary hashes, build logs and metric results. These generated outputs remain unversioned.

`metrics.py` only reads those four cases; it does not launch a solver. Run it from this directory after authorized CFD execution. Integer GlobalDof identities remain integers when reading matrix CSV, never floating-point keys.

Precise WHAT/STATE/HOW/WHICH, geometry, IC/BC, units, output settings, commands, actual findings and limitations: local `docs/loop.md`. Original Sod/Poiseuille inputs and results remain untouched.
