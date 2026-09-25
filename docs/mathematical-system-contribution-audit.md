# Mathematical System Contribution Audit

| Module | Current entry | Mathematical contribution | Runtime path | Status |
|---|---|---|---|---|
| EOS/caloric/transport | composition + EquationSet | Closure | property/thermo refresh | Implemented |
| Gravity/MRF/wall heat | source contribution | Equation term | existing phase/source assembly | Legacy |
| turbulence | `addTurbulenceEquations` | closure + optional equations | Eulerian turbulence provider | Legacy |
| Eulerian multiphase | phase template | phase equations + shared constraints + interphase closures/terms | PhaseSystem / PressureStepper callbacks | Legacy |
| Ghost/ILW IBM | IBM descriptor + boundary pipeline | boundary/interface closure | physical BC -> ghost/ILW -> halo | Legacy |
| Peskin | immersed descriptor | constraint | `projectPredictedState()` | Legacy |
| DFM explicit/fractional | immersed descriptor | constraint | `projectPredictedState()` | Legacy |
| FTS/BP/Brinkman | immersed descriptor | constraint/penalty term | `projectPredictedState()` | Legacy |
| monolithic DLM/KKT | immersed descriptor + transform | multiplier unknown + constraint + coupled block | prepare/solve/accept KKT | Legacy |
| augmented Lagrangian | immersed descriptor + transform | multiplier unknown + constraint + coupled block | existing KKT path | Legacy |

`System::addImmersed` contributes mathematical descriptors; it does not construct MPI/HYPRE runtime objects. `projectPredictedState()` selects an existing built-in local formula and does not own outer PIMPLE/RK/timestep control flow. Ghost IBM remains boundary closure.
