# Dynamic equation feedback inputs

Prepared only. No solver, check/explain, regression or numerical validation was run for this task. Run only after the user specifies these cases (AGENTS.md §35).

- `density/baseline`: uniform flow copied from Sod's mesh input, classical RK4, 0–1 s, no source.
- `density/feedback`: same input plus Stage relations for axial force density and mechanical power; both feed actual conservative RHS.
- `pressure/baseline`: constant-density transient PISO Poiseuille, 0–10 s (endStep cap 2200), original pressure drop and no-slip plate boundaries.
- `pressure/feedback`: same input plus a Step Current-velocity relation feeding actual predictor acceleration.

All YAML is production native IO: one registry in each object directory, solver documents directly in `solvers/`. `prepare.py` creates isolated inputs and never invokes a binary. It refuses to overwrite a case. Generated mesh and numerical output remain in each case's `mesh/` and `result/`; logs/checkpoints must stay here. Original cases/results remain untouched.

Precise settings, equations, ownership, source units, proposed comparisons and unresolved limits: local `docs/loop.md`. Start snapshot/build logs/hashes: local `evidence/`. No history evidence/results are uploaded with this fixture.
