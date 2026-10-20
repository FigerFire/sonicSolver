# Phase 3A Numerical Baseline

Captured before Phase 3A source edits on 2026-09-11.

The totals and norms below are calculated from the final VTK point fields using
tensor-product trapezoidal integration over each written structured piece.
Adjacent pieces are each integrated over their own extent, so their shared
boundary points have the usual half weight on both sides. This is a lightweight
post-processing record, not a new solver regression framework.

## Four-rank Sod

Command: `mpirun -np 4 ./build/sonicSolver test/sodCase`

The case stops at `endStep=500`, physical time `2.435715e-01`.

| quantity | value |
| --- | ---: |
| first six `(time, dt)` | `(6.681531e-04, 6.681531e-04)`, `(1.272174e-03, 6.040207e-04)`, `(1.850019e-03, 5.778452e-04)`, `(2.416037e-03, 5.660178e-04)`, `(2.975014e-03, 5.589775e-04)`, `(3.529599e-03, 5.545843e-04)` |
| last `(time, dt)` | `(2.435715e-01, 5.352234e-04)` |
| total mass | `9.779741906218e+01` |
| total momentum `(x, y, z)` | `(2.361223645979e+04, 6.996363528082e-12, 0)` |
| total energy | `2.169411924704e+07` |
| min / max rho | `4.132308946792e-01` / `7.507655531855e-01` |
| min / max p | `3.026904794569e+04` / `6.793422191548e+04` |
| L2 rho / p | `7.008343666007e+00` / `5.431969111109e+05` |

## Single-patch density-based Sod

Command: `./build/sonicSolver test/sodCase_weno7Check`

The serial case uses the density-based EquationSet path and reaches physical
time `5.0` after 9438 steps.

| quantity | value |
| --- | ---: |
| first six `(time, dt)` | `(6.681531e-04, 6.681531e-04)`, `(1.272174e-03, 6.040207e-04)`, `(1.850019e-03, 5.778452e-04)`, `(2.416037e-03, 5.660178e-04)`, `(2.975014e-03, 5.589775e-04)`, `(3.529599e-03, 5.545843e-04)` |
| last `(time, dt)` | `(5.000000e+00, 3.206358e-04)` |
| total mass | `8.392751270150e+01` |
| total momentum `(x, y, z)` | `(2.372785305663e+04, 3.522265384126e-12, 0)` |
| total energy | `1.909181973501e+07` |
| min / max rho | `4.086256292029e-01` / `6.810781004577e-01` |
| min / max p | `3.024601352767e+04` / `6.060132650364e+04` |
| L2 rho / p | `5.951833779191e+00` / `4.471500421890e+05` |
