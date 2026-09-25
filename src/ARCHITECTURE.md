# SonicSolver Architecture — Equation Composition, System Formulation, Solve Planning

**Revision:** 2026-09-24 v2  
**Purpose:** converge the solver architecture after the pressure–velocity coupling design review.

---

## 1. One sentence architecture

SonicSolver is not a collection of solver families.

It is a system that:

```text
composes equations
    ↓
formulates/transforms constrained equation systems
    ↓
derives executable equations/operators
    ↓
plans their execution
    ↓
binds numerical/runtime providers
    ↓
executes
```

The stable design rule is:

> Presets and models never own a hidden solver path. They only contribute fields, equations, constraints, transformations, numerical requirements, and solve-plan fragments to the same compiled system that an advanced user could construct manually.

---

## 2. The three fundamental questions

### 2.1 Equation Composition — WHAT physics/mathematics exists?

This layer answers:

```text
What fields exist?
What equations exist?
What closures exist?
What constraints exist?
```

Examples:

```text
Mass
Momentum
Energy
Alpha
k
omega
solid momentum
induction equation
IBM constraint
incompressibility constraint
```

Sources may be:

```text
built-in presets
registered models
user additions
user replacements/overrides
```

All sources must lower to the same equation IR.

---

### 2.2 System Formulation / Transformation — how is the mathematical system made solvable?

This layer is between raw equations and execution.

It may:

```text
introduce a multiplier
derive an algorithmic equation
eliminate a variable
construct a Schur complement
create a pressure-correction equation
construct a projection
form a KKT system
construct an augmented-Lagrangian system
split a constrained solve into predictor / multiplier / correction
```

This is where pressure–velocity coupling belongs conceptually.

For constant-density flow:

```text
Momentum
+
div(U) = 0
    ↓
PressureConstraintFormulation
    ↓
MomentumPredictor
PressureCorrection/Poisson
VelocityCorrection
FluxCorrection
```

For IBM:

```text
Momentum
+
J(U) = Ub
    ↓
ImmersedConstraintFormulation
    ↓
Monolithic KKT
or
Predictor + multiplier solve + correction
```

Pressure coupling and IBM share the architecture of:

```text
base equation
+
constraint
+
multiplier
+
formulation
```

but they do not have to share one numerical implementation.

---

### 2.3 Equation Execution — in what order/how are the executable operations solved?

This layer answers:

```text
explicit or implicit?
segregated or block coupled?
A → B → C?
repeat A/B to convergence?
subcycle?
RK stages?
outer correctors?
pressure correctors?
non-orthogonal correctors?
```

The output is the `CompiledSolvePlan`.

Typical control-flow IR:

```text
Sequence
Loop
StageLoop
Assemble
Solve
Correct
Update
Commit
```

`PlanExecutor` interprets this IR.

It does not infer physics and does not derive equations.

---

## 3. Correct position of SIMPLE / PISO / PIMPLE

SIMPLE, PISO and PIMPLE are not only solve-order labels.

They have two responsibilities:

```text
A. formulation / factorization
   derive the executable pressure-constraint operators/equations

B. execution pattern
   define predictor/corrector/outer-loop ordering
```

Therefore a coupling model should contribute both:

```text
EquationSystemTransformer
+
SolvePlanContribution
```

Conceptually:

```text
SIMPLE
    inspect Momentum + IncompressibilityConstraint
    derive pressure-correction formulation
    add predictor / pressure / velocity / flux correction operations
    contribute outer iteration + relaxation schedule

PISO
    inspect the same constrained base system
    derive compatible projection/correction operators
    contribute predictor + N corrector schedule

PIMPLE
    combine an outer SIMPLE-like loop with inner PISO-like corrections
    without becoming a top-level solver identity
```

Do not implement:

```cpp
if (simple) runSimpleSolver();
if (piso)   runPisoSolver();
```

Do not reduce them to a single enum that is interpreted only inside a stepper.

---

## 4. Compressible / incompressible are presets, not solver identities

`compressible` and `incompressible` may exist as convenient built-in equation presets.

Example:

```text
compressible preset
    contributes:
        Mass
        Momentum
        Energy
        thermodynamic closure

incompressible preset
    contributes:
        Momentum
        rho = rho0 closure
        div(U) = 0 constraint
```

These presets must not directly choose:

```text
a stepper class
SIMPLE/PISO/PIMPLE
RK/Euler
WENO/TENO
MPI
HYPRE
IBM method
```

After composition, preset identity remains provenance only.

An advanced user may construct an equivalent system manually.

Equivalent resolved equations/constraints/models must produce equivalent compilation regardless of whether they came from a preset.

---

## 5. Density-based / pressure-based are not top-level configuration identities

Do not keep `DensityBasedSolver` / `PressureBasedSolver` as architecture categories.

Do not make a user-facing enum the authority that decides the physical equation pack.

Instead introduce an internal compiled object such as:

```cpp
struct StateRealization {
    // primary unknown roles
    // transported vs derived variables
    // multiplier roles
    // conservative groups
    // storage bindings
};
```

It is derived from:

```text
ExecutableEquationSystem
+
active formulations
+
unknown/constraint roles
```

Examples:

```text
Mass + Momentum + Energy
rho/rhoU/rhoE transported
no pressure multiplier constraint
    -> conservative transported-state realization

Momentum + rhoConst + div(U)=0
pressure multiplier
pressure-correction formulation
    -> U/p pressure-constraint realization
```

Thus “density-based” and “pressure-based” may remain useful `explain` descriptions, but they are compiled consequences, not top-level solver selectors.

Temporary compatibility input such as `densityBase` / `pressureBase` may be decoded during migration, but it must not survive as runtime dispatch authority.

---

## 6. Models are system contributors

A model is not primarily “an object called once per timestep”.

A registered model contributes to the mathematical/execution system.

A model may contribute any subset of:

```text
fields
equations
terms
closures
constraints
transformations/formulations
numerical requirements
solve-plan fragments
runtime requirements
```

Examples:

```text
kOmegaSST
    fields: k, omega, nut
    equations: E_k, E_omega
    closure: nut(...)
    terms: turbulence contribution to flow equations

IBM
    constraint: J(U)-Ub=0
    multiplier: lambda
    term: J^T lambda
    formulation: KKT / projection / forcing
    runtime requirement: geometry / marker ownership

SIMPLE
    requirement: Momentum + incompressibility constraint
    formulation: pressure correction
    operations: predictor / pressure / velocity / flux correction
    plan: outer iteration + relaxation

PISO
    requirement: Momentum + incompressibility constraint
    formulation: pressure projection/correction
    plan: predictor + N correctors
```

No model has a privileged hidden runtime path.

---

## 7. Registered-but-inapplicable models

Never silently ignore a model.

A registered model/formulation has explicit status:

```text
not registered
active
inactive
invalid
unsupported
```

Example:

```text
User registers SIMPLE.

Raw system:
    transported rho
    conservative mass equation
    momentum
    energy
    no incompressibility constraint

Result:
    SIMPLE inactive/incompatible

Reason:
    no supported pressure-multiplier constraint system was found.
```

Do not implement:

```cpp
if (densityBased) ignoreSimple = true;
```

The decision comes from equation/constraint requirements, not a solver-family label.

In strict mode, an explicitly registered but incompatible coupling model should be a configuration error.

---

## 8. Presets and user freedom

All usage levels share one compilation pipeline.

### Level 0 — fully default

```text
default fields
default equation preset
default numerics
default execution policies
```

### Level 1 — model/preset composition

```text
register SIMPLE/PISO/PIMPLE
register IBM
register turbulence
register phase change
register MRF
```

### Level 2 — reuse fields, customize equations

```text
extend momentum with Lorentz force
disable energy
replace one transport equation
```

### Level 3 — full custom system

```text
register fields
define equations
define constraints
define formulations
define solve plan
```

All four lower to:

```text
ResolvedSimulationSystem
```

There is only one runtime authority chain.

---

## 9. Override semantics

Do not rely on same-name “last definition wins”.

Use explicit actions:

```text
add
extend
replace
disable
```

Suggested precedence:

```text
built-in defaults
    ↓
registered model contributions
    ↓
user additions/extensions
    ↓
explicit user replacements/disables
    ↓
validation
```

Examples:

```text
extend Momentum:
    + LorentzForce

replace Momentum:
    MyMomentumEquation

disable Energy
```

This is important for reproducible research cases.

---

## 10. Compilation pipeline

Target chain:

```text
User Input
    │
    ├── Built-in Presets
    ├── Registered Models
    └── User Definitions
            ↓
      System Composition
            ↓
      Raw Equation System
            ↓
      System Formulations / Transformers
            │
            ├── Pressure constraint formulation
            ├── IBM constraint formulation
            ├── shared-pressure formulation
            └── custom
            ↓
      Executable Equation System                    WHAT
            ↓
      State Realization
            ↓
      Numerical Compiler                            HOW
            ↓
      Compiled Numerical System
            ↓
      Solve Planner                                 ORDER
            ↓
      Compiled Solve Plan
            ↓
      Provider + Runtime Requirements               RUNTIME
            ↓
      validate()
            ↓
      Runtime provider binding
            ↓
      OpRegistry
            ↓
      PlanExecutor
            ↓
      committed state
```

The layers have distinct authority.

No later stage may reconstruct an earlier decision from raw configuration.

---

## 11. Runtime naming

The current `SingleFluidStepper` name is wrong because it describes a physics label that no longer owns the runtime path.

Do **not** replace it with another top-level state-family solver if the object can be made generic.

Preferred target:

```text
SingleFluidStepper -> SingleFluidStepper
```

or simply:

```text
FlowStepper
```

Its responsibility:

```text
own runtime operation registry/workspaces
bind compiled operations
execute CompiledSolvePlan
commit state
```

Subordinate providers may remain specialized:

```text
ConservativeRHS
PressurePredictor
PressureCorrection
IBMConstraintProvider
TurbulenceProvider
```

Therefore:

```text
ConservativeRHS -> ConservativeRHS
```

is still an appropriate rename.

Likewise:

```text
EulerianEulerian::EulerianStepper -> EulerianStepper
```

because Eulerian phase execution is its real domain; PIMPLE is a formulation/plan contribution, not its identity.

Avoid:

```text
CompressibleSolver
IncompressibleSolver
DensitySolver
PressureSolver
PimpleSolver
```

as top-level runtime identities.

---

## 12. Pressure coupling architecture

The base mathematical system should contain the constraint, not a user-written pressure Poisson equation:

```text
Momentum
+
IncompressibilityConstraint
```

The coupling formulation derives algorithmic equations/operators:

```text
MomentumPredictor
PressureCorrection/Poisson
VelocityCorrection
FluxCorrection
```

The solve planner then lowers the registered coupling preset to control flow.

Example PIMPLE:

```text
Sequence step
    prepare
    dt.compute
    step.begin

    Loop outerCorrectors
        MomentumPredictor

        Loop pressureCorrectors
            Loop nonOrthogonalCorrectors
                PressureAssemble
                PressureSolve

            VelocityCorrection
            FluxCorrection

        convergence / relaxation policy

    commit
    time.commit
```

The exact mathematical operators remain owned by the formulation/provider.

The plan owns only execution ordering and loop structure.

---

## 13. Relationship to KKT / constrained systems

Pressure coupling and IBM can share generic metadata/contracts such as:

```text
ConstraintDescriptor
MultiplierDescriptor
ConstraintOperator role
Adjoint/spreading role
Formulation requirement
```

Potential common formulation categories:

```text
Monolithic
Projection
Fractional
ApproximateSchur
AugmentedLagrangian
```

However, do not force SIMPLE/PISO and IBM KKT into one numerical implementation merely because both originate from constrained saddle-point systems.

Architecture sharing is encouraged.

Mathematical over-generalization is not.

---

## 14. SystemBuilder responsibilities

The current `SF_systemBuilder.cpp` must converge to a compiler orchestrator, not a God file.

Logical stages:

```text
compose defaults/presets/models/user definitions
apply explicit overrides
validate raw system
select/apply formulations
build executable system
realize state
compile numerics
compile solve plan
derive provider/runtime requirements
final validation
```

Split implementation by real stage while keeping one small public build entry.

Do not introduce a second workflow framework.

---

## 15. P0 priorities

### P0-A — remove solver-family authority

- remove runtime dispatch from `compressible/incompressible`
- remove `DensityBased/PressureBased` as physical equation selectors
- stop using a top-level `PrimaryForm` to choose the equation pack
- introduce compiled `StateRealization`
- keep old input names only as migration aliases if required

### P0-B — reposition pressure coupling

- SIMPLE/PISO/PIMPLE become coupling-model presets
- each contributes a system formulation + solve-plan fragment
- planner must not derive the pressure equation itself
- stepper must not own the SIMPLE/PISO/PIMPLE lifecycle

### P0-C — rename runtime classes

- `SingleFluidStepper -> SingleFluidStepper` (preferred)
- `ConservativeRHS -> ConservativeRHS`
- `Eulerian EulerianStepper -> EulerianStepper`
- remove physical-family wording from comments and explain output where it implies runtime identity

### P0-D — compiled authority

- every stepper consumes compiled numerical/time/dt authority
- no Eulerian/raw-config reread of CFL/maxDeltaT when compiled values exist
- provider selection comes from compiled operation/provider requirements

### P0-E — truthful rhoConst + PIMPLE path

- constant-density composition creates/retains incompressibility constraint
- PIMPLE formulation derives required algorithmic pressure operations
- if numerical provider exists, run it
- if not, fail specifically as `Unsupported: missing <operation/provider>`
- never fail merely because “this is not a pressure solver”
- never fake PIMPLE by repeating a full physical timestep

---

## 16. P1 priorities

- split `SF_systemBuilder.cpp` by compilation stage
- return contribution `.cpp` files to their owning model targets
- return `PlanExecutor` to the run target
- reduce application execution files to runtime composition
- make `OpRegistry`/workspace lifetime consistent across steppers
- remove ambiguous duplicate public header names or use path-qualified includes
- add architecture guards for forbidden solver-family dispatch

---

## 17. P2 priorities

- split `SF_nativeDecode.cpp` by input domain
- extract Eulerian diagnostics from commit
- classify remaining raw config reads into model / closure / numerical authority
- delete legacy parser/config files already out of production targets
- remove obsolete include directories
- implement/document explicit add/extend/replace/disable equation semantics
- update `explain` to show contribution provenance and active/inactive/unsupported models

---

## 18. Explain output target

`explain` should expose the compiled reasoning without hidden dispatch.

Example:

```text
EQUATION SOURCES
  builtin preset : incompressible
  model          : kOmegaSST
  coupling model : PIMPLE
  model          : IBM ghost

RAW EQUATIONS
  E_MOMENTUM
  C_INCOMPRESSIBILITY
  E_K
  E_OMEGA

FORMULATIONS
  pressure.constraint.pimple       active
  ibm.ghost                        active

DERIVED OPERATIONS
  momentum.predict
  pressure.assemble
  pressure.solve
  velocity.correct
  flux.correct

STATE REALIZATION
  transported : U, k, omega
  multiplier  : p
  derived     : rho=rho0

NUMERICAL HOW
  convection  : ...
  diffusion   : ...
  dt policy   : ...

SOLVE PLAN
  outer loop ...
  pressure corrector loop ...

PROVIDERS
  ...

RUNTIME
  MPI ...
  HYPRE ...
```

A registered but irrelevant coupling model must be visible:

```text
coupling SIMPLE : inactive
reason: incompressibility/pressure-multiplier constraint not present
```

No silent ignore.

---

## 19. Acceptance criteria

Architecture is considered converged when:

```text
[ ] no runtime solver family is selected from compressible/incompressible
[ ] density-based/pressure-based are descriptions, not top-level dispatch authorities
[ ] state realization is compiled from the executable system/formulations
[ ] SIMPLE/PISO/PIMPLE contribute both formulation and solve-plan semantics
[ ] no `runSimpleSolver()` / `runPisoSolver()` hidden path exists
[ ] stepper classes are named by runtime domain, not compressibility/coupling strategy
[ ] a coupling preset is activated by equation/constraint requirements
[ ] explicit incompatible model registration is never silently ignored
[ ] presets and equivalent manual definitions lower to the same IR
[ ] builtin equations and user equations use the same Equation IR
[ ] compiled HOW/ORDER/RUNTIME remain authoritative
[ ] CMake ownership matches module ownership
[ ] no CFD formula is changed by architecture cleanup
```

Once these conditions hold, freeze the architecture and return focus to numerical/physical capability.
