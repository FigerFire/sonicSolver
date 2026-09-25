/// @file SF_pressureStepper.cpp
/// @brief 双欧拉各相预测、共享压力校正和能量更新的数学顺序。

#include "solver/algorithm/eulerian/SF_eulerianStepper.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/run/SF_planExecutor.h"

#include "SF_phaseGravity.h"
#include "SF_phaseMRF.h"
#include "SF_phaseWallHeat.h"
#include "SF_phaseChange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace SF::EulerianEulerian {

namespace {
std::vector<std::string> equationIds(
        const System::ExecutableEquationSystem& equations) {
    std::vector<std::string> ids;
    ids.reserve(equations.equations.size());
    for (const auto& equation : equations.equations) {
        ids.push_back(equation.id);
    }
    return ids;
}

struct PhaseSourceRegistration {
    FDM::SourceKind kind;
    bool wallHeat;
    std::unique_ptr<Physics::PhaseSystems::PhaseEquationSource> (*make)(
        const FDM::SourceConfig&);
};

const std::array<PhaseSourceRegistration,3>& phaseSourceContributions() {
    static const std::array<PhaseSourceRegistration,3> registry{{
        {FDM::SourceKind::Gravity,false,
         [](const FDM::SourceConfig& config) {
             return Physics::PhaseSystems::makePhaseGravitySource(
                 config.gravity);
         }},
        {FDM::SourceKind::MRF,false,
         [](const FDM::SourceConfig& config) {
             return Physics::PhaseSystems::makePhaseMRFSource(
                 config.rotating);
         }},
        {FDM::SourceKind::WallHeat,true,
         [](const FDM::SourceConfig& config) {
             return Physics::PhaseSystems::makePhaseWallHeatSource(
                 config.wallHeat);
         }}
    }};
    return registry;
}

State::VariableDescriptor distributedScalar(
        const std::string& name, int depth,
        State::UpdatePolicy policy = State::UpdatePolicy::LocalImplicit,
        State::FieldLocation location = State::FieldLocation::Cell) {
    State::VariableDescriptor descriptor;
    descriptor.name = name;
    descriptor.location = location;
    descriptor.updatePolicy = policy;
    descriptor.halo.depth = depth;
    descriptor.halo.stages = State::HaloSyncStage::BeforeRHS
        | State::HaloSyncStage::AfterUpdate
        | State::HaloSyncStage::BeforeOutput;
    return descriptor;
}
}

EulerianStepper::EulerianStepper(
        Physics::PhaseSystems::PhaseSystem& system,
        FDM::SolverConfig config,
        const System::ExecutableEquationSystem& equations,
        const System::CompiledNumericalSystem& numerics,
        const System::CompiledSolvePlan& solvePlan,
        const System::RuntimeRequirements& requirements)
    : system_(system), executable_(equations), numerics_(numerics),
      solve_(solvePlan),
      runtime_(requirements), config_(std::move(config)),
      turbulence_(config_.turbulence),
      assemblyPlans_(equations.equationDefinitions, equationIds(equations)),
      equations_(system_, config_.pressure, workspace_) {
    FDM::validatePressureCorrectionConfig(config_.pressure);
    // compiled HOW 与 raw config 必须一致（runtime 只执行 compiled policy）。
    if (numerics_.dt.cfl != config_.numerics.cfl
        || numerics_.dt.maxDeltaT != config_.numerics.maxDeltaT) {
        throw std::runtime_error(
            "EulerianStepper received a compiled dt policy that disagrees "
            "with its solver configuration.");
    }
    if (numerics_.phaseTransport.convection
        != FDM::PhaseConvectionScheme::Upwind) {
        throw std::runtime_error(
            "EulerianStepper has no provider for the compiled phase "
            "convection scheme.");
    }
    // Eulerian runtime 的 domain 是 phase state + shared-pressure constraint；
    // 判据来自 executable system，而不是 density/pressure 标签。
    if (!System::hasConstraint(equations,"C_SHARED_PRESSURE")) {
        throw std::runtime_error(
            "EulerianStepper requires the shared-pressure constraint "
            "formulation in its executable equation system.");
    }
    bool wallHeatEnabled = false;
    for (FDM::SourceKind kind : config_.sources.enabled) {
        const auto found = std::find_if(
            phaseSourceContributions().begin(),
            phaseSourceContributions().end(),
            [kind](const PhaseSourceRegistration& item) {
                return item.kind == kind;
            });
        if (found == phaseSourceContributions().end()) {
            throw std::runtime_error(
                "No Eulerian equation contribution is registered for SourceKind.");
        }
        wallHeatEnabled = wallHeatEnabled || found->wallHeat;
        sourceRegistry_.add(found->make(config_.sources));
    }
    if (Physics::PhaseChange::normalizeModel(
            system_.config().phaseChange.model) == "rpi"
        && !wallHeatEnabled) {
        throw std::runtime_error(
            "Eulerian RPI requires a wallHeatFlux boundary so the unified "
            "wall-boiling ledger can assemble mass, momentum and enthalpy.");
    }
    workspace_.setupLike(system_);
    turbulence_.initialize(system_);
    if (turbulence_.hasTransportEquations()
        != System::hasEquationPrefix(executable_,"E_TURB_")) {
        throw std::runtime_error(
            "Resolved turbulence equations do not match the active "
            "Eulerian turbulence model.");
    }
    equations_.attachTurbulence(&turbulence_);
    equations_.bindAssemblyPlans(assemblyPlans_);
}

void EulerianStepper::registerState(State::StateBundle& state) {
    Field& geometry = const_cast<Field&>(system_.geometry());
    const int depth = geometry.NG();
    state.registerAuxiliary(
        distributedScalar("pressure", depth,
                          State::UpdatePolicy::PressureCorrection),
        0, geometry, system_.sharedPressure());
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        auto& phaseState = system_.phases()[phase];
        auto& primary = phaseState.primary;
        auto& primitive = phaseState.primitive;
        const std::string prefix = "phase" + std::to_string(phase) + ".";
        auto addLocalScalar = [&](const std::string& suffix,
                                  ScalarField& value) {
            state.distributed.add(State::scalarView(
                prefix + suffix, 0, geometry, value, 0,
                State::HaloSyncStage::None, State::ExchangeKind::None));
        };
        auto addLocalVector = [&](const std::string& suffix,
                                  Physics::PhaseSystems::PhaseVectorField& value) {
            state.distributed.add(State::scalarComponentsView(
                prefix + suffix, 0, geometry,
                {&value[0], &value[1], &value[2]}, 0,
                State::HaloSyncStage::None, State::ExchangeKind::None));
        };
        state.registerAuxiliary(
            distributedScalar(prefix + "mass", depth),
            0, geometry, primary.phaseMass);
        state.distributed.add(State::scalarComponentsView(
            prefix + "momentum", 0, geometry,
            {&primary.momentum[0], &primary.momentum[1],
             &primary.momentum[2]}, depth,
            State::HaloSyncStage::BeforeRHS
                | State::HaloSyncStage::AfterUpdate
                | State::HaloSyncStage::BeforeOutput));
        state.registerAuxiliary(
            distributedScalar(prefix + "enthalpy", depth),
            0, geometry, primary.phaseEnthalpy);
        addLocalScalar("alpha", primitive.alpha);
        addLocalScalar("density", primitive.density);
        addLocalVector("velocity", primitive.velocity);
        addLocalScalar("temperature", primitive.temperature);
        addLocalScalar("primitiveEnthalpy", primitive.enthalpy);

        auto& diagonal = workspace_.momentumDiagonal.at(phase);
        state.registerAuxiliary(
            distributedScalar(prefix + "momentumDiagonal", depth,
                              State::UpdatePolicy::DerivedOnly),
            0, geometry, diagonal);
        addLocalVector("predictedVelocity",
                       workspace_.predictedVelocity.at(phase));
        addLocalVector("previousVelocity",
                       workspace_.previousVelocity.at(phase));
        addLocalScalar("previousMass",
                       workspace_.previousPhaseMass.at(phase));
        addLocalVector("previousMomentum",
                       workspace_.previousMomentum.at(phase));
        addLocalScalar("previousEnthalpy",
                       workspace_.previousPhaseEnthalpy.at(phase));
        state.distributed.add(State::workspaceView(
            prefix + "volumeFaceFlux", 0, geometry,
            workspace_.faceFlux.at(phase).volume, 3, depth,
            State::HaloSyncStage::None,
            State::ExchangeKind::CanonicalFaceFlux,
            State::FieldLocation::Face));
        state.distributed.add(State::workspaceView(
            prefix + "massFaceFlux", 0, geometry,
            workspace_.faceFlux.at(phase).mass, 3, depth,
            State::HaloSyncStage::None,
            State::ExchangeKind::CanonicalFaceFlux,
            State::FieldLocation::Face));

        if (phase < turbulence_.states().size()) {
            auto& turbulence = turbulence_.states()[phase];
            auto addTurbulence = [&](const std::string& suffix,
                                     ScalarField& value) {
                if (value.empty()) return;
                state.registerAuxiliary(
                    distributedScalar(prefix + suffix, depth),
                    0, geometry, value);
            };
            addTurbulence("k", turbulence.kineticEnergy.variable);
            addTurbulence("epsilon", turbulence.dissipation.variable);
            addTurbulence("omega", turbulence.specificDissipation.variable);
            addTurbulence("mu_t", turbulence.eddyViscosity);
        }
    }
    state.registerAuxiliary(
        distributedScalar("pressureCorrection", depth,
                          State::UpdatePolicy::PressureCorrection),
        0, geometry, workspace_.pressureCorrection);
    state.registerAuxiliary(
        distributedScalar("previousPressure", depth,
                          State::UpdatePolicy::DerivedOnly),
        0, geometry, workspace_.previousPressure);
}

void EulerianStepper::applyBoundaryAndSynchronize() {
    boundary_.applySharedPressure(
        system_.geometry(), system_.sharedPressure(),
        config_.boundaries.energyFromPressure);
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        boundary_.apply(system_.geometry(), system_.phaseProperties(phase),
                        system_.phases()[phase]);
    }
    synchronizePrimaryState();
    turbulence_.applyBoundary(system_);
    synchronizeTurbulenceState();
}

void EulerianStepper::synchronizeTurbulenceState() {
    if (!services_.executionRuntime || !turbulence_.active()) return;
    const int depth = system_.geometry().NG();
    std::vector<Execution::FieldAccess> writes;
    std::vector<Execution::FieldAccess> reads;
    for (size_t phase = 0; phase < turbulence_.states().size(); ++phase) {
        const auto& state = turbulence_.states()[phase];
        const std::string prefix = "phase" + std::to_string(phase) + ".";
        const std::vector<std::pair<std::string, const ScalarField*>> fields{
            {prefix + "k", &state.kineticEnergy.variable},
            {prefix + "epsilon", &state.dissipation.variable},
            {prefix + "omega", &state.specificDissipation.variable},
            {prefix + "mu_t", &state.eddyViscosity}};
        for (const auto& item : fields) {
            if (item.second->empty()) continue;
            writes.push_back(Execution::writeOwned(item.first));
            reads.push_back(Execution::readHalo(item.first,depth));
        }
    }
    if (!writes.empty()) {
        services_.executionRuntime->finalize({"Eulerian turbulence update",writes});
        services_.executionRuntime->prepare({"Eulerian turbulence stencil",reads});
    }
}

void EulerianStepper::synchronizePrimaryState() {
    if (services_.executionRuntime) {
        const int depth = system_.geometry().NG();
        std::vector<Execution::FieldAccess> writes{
            Execution::writeOwned("pressure")};
        std::vector<Execution::FieldAccess> reads{
            Execution::readHalo("pressure",depth)};
        for (size_t phase=0; phase<system_.phases().size(); ++phase) {
            const std::string prefix="phase"+std::to_string(phase)+".";
            for (const auto* suffix : {"mass","momentum","enthalpy"}) {
                writes.push_back(Execution::writeOwned(prefix+suffix));
                reads.push_back(Execution::readHalo(prefix+suffix,depth));
            }
        }
        services_.executionRuntime->finalize({"Eulerian primary-state update",writes});
        services_.executionRuntime->prepare({"Eulerian primary-state stencil",reads});
    }
    system_.recoverPrimitiveState();
}

void EulerianStepper::synchronizePressureCorrection() {
    if (!services_.executionRuntime) return;
    const int depth=system_.geometry().NG();
    services_.executionRuntime->finalize({
        "Eulerian pressure-correction update",
        {Execution::writeOwned("pressureCorrection")}});
    services_.executionRuntime->prepare({
        "Eulerian pressure-correction stencil",
        {Execution::readHalo("pressureCorrection",depth)}});
}

void EulerianStepper::synchronizeMomentumDiagonal() {
    if (!services_.executionRuntime) return;
    const int depth=system_.geometry().NG();
    std::vector<Execution::FieldAccess> writes;
    std::vector<Execution::FieldAccess> reads;
    for (size_t phase=0; phase<workspace_.momentumDiagonal.size(); ++phase) {
        const std::string name="phase"+std::to_string(phase)
            +".momentumDiagonal";
        writes.push_back(Execution::writeOwned(name));
        reads.push_back(Execution::readHalo(name,depth));
    }
    services_.executionRuntime->finalize({"Eulerian momentum-diagonal update",writes});
    services_.executionRuntime->prepare({"Eulerian momentum interpolation stencil",reads});
}

void EulerianStepper::assembleCanonicalPhaseFlux() {
    if (!services_.executionRuntime) return;
    if (!state_) {
        throw std::runtime_error(
            "Eulerian canonical face synchronization requires bound state.");
    }
    for (size_t phase=0; phase<workspace_.faceFlux.size(); ++phase) {
        const std::string prefix="phase"+std::to_string(phase)+".";
        for (const std::string& name : {
                 prefix+"volumeFaceFlux", prefix+"massFaceFlux"}) {
            const auto registered = state_->distributed.select(
                name, State::HaloSyncStage::None);
            if (registered.empty()) {
                throw std::runtime_error(
                    "Eulerian canonical face workspace '"+name
                    +"' was not bound before timestep execution.");
            }
            std::vector<State::DistributedFieldView> borrowed;
            borrowed.reserve(registered.size());
            for (const auto* view : registered) borrowed.push_back(*view);
            // Face workspaces remain solver-owned.  The runtime only borrows
            // their registered views to apply the canonical owner -> COPY
            // synchronization; no second face-flux storage is created.
            services_.executionRuntime->synchronizeTransient(borrowed);
        }
    }
}

double EulerianStepper::stableTimeStep(double cfl) {
    const double transportLimit = equations_.stableTimeStep(cfl);
    // compiled dt policy 是 dt 上限的唯一 authority。
    const double probeDt = std::min(transportLimit, numerics_.dt.maxDeltaT);
    if (!std::isfinite(probeDt) || probeDt <= 0.0) {
        throw std::runtime_error(
            "Eulerian stable-time-step probe is invalid.");
    }
    applyBoundaryAndSynchronize();
    system_.computeInterphase(probeDt, workspace_.previousVelocity);
    sourceRegistry_.assemble(system_, probeDt);
    turbulence_.prepare(system_);
    synchronizeTurbulenceState();
    return std::min(
        std::min(
            transportLimit,
            equations_.sourceTimeStep(
                numerics_.phaseTransport.sourceCfl)),
        turbulence_.sourceTimeStep(
            system_, numerics_.phaseTransport.sourceCfl));
}

void EulerianStepper::bindServices(FDM::SolverServices services) {
    if (servicesBound_) {
        throw std::runtime_error(
            "EulerianStepper services are already bound for this lifecycle.");
    }
    services_ = services;
    servicesBound_ = true;
}

void EulerianStepper::bindSolvePlan(
        const System::CompiledSolvePlan& plan) {
    if (solveStagesBound_) {
        throw std::runtime_error(
            "EulerianStepper solve plan is already bound.");
    }
    if (&plan != &solve_) {
        throw std::runtime_error(
            "EulerianStepper received a plan other than its resolved plan.");
    }
    for (const auto& block : plan.blocks) {
        pimpleStageActive_ = pimpleStageActive_
            || block.strategyKind == FDM::SolveStrategyKind::PressureVelocityCoupling;
    }
    if (!pimpleStageActive_) {
        throw std::runtime_error(
            "Eulerian execution requires a pressure-velocity coupling strategy.");
    }
    const auto required = System::SolvePlanner::requiredOperations(plan);
    const bool planSolvesTurbulence = std::find(
        required.begin(),required.end(),"ee.turbulence.solve")
        != required.end();
    if (planSolvesTurbulence != turbulence_.hasTransportEquations()) {
        throw std::runtime_error(
            "Compiled Plan turbulence operations do not match the resolved "
            "Eulerian turbulence equations.");
    }
    solveStagesBound_ = true;
}

void EulerianStepper::bindState(State::StateBundle& state) {
    state.validatePatches();
    if (&state.singlePatch() != &system_.geometry()) {
        throw std::runtime_error(
            "Eulerian EulerianStepper requires its PhaseSystem geometry "
            "as the single StateBundle patch.");
    }
    if (state_ && state_ != &state) {
        throw std::runtime_error(
            "EulerianStepper cannot switch StateBundle after stepping begins.");
    }
    state_ = &state;
    equations_.setExecutionRuntime(services_.executionRuntime);
    if (state.distributed.empty()) {
        state.registerConservativeState();
        registerState(state);
    }
    if (services_.executionRuntime) services_.executionRuntime->attachState(state);
}

void EulerianStepper::prepare(FDM::SolverState& state) {
    if (!servicesBound_ || !solveStagesBound_) {
        throw std::runtime_error(
            "EulerianStepper requires services and a solve plan before state preparation.");
    }
    if (!state.bundle) {
        throw std::runtime_error(
            "EulerianStepper::prepare requires a StateBundle.");
    }
    bindState(*state.bundle);
    realizedState_ = System::realizeState(executable_,runtime_,*state.bundle);
}

FDM::StepResult EulerianStepper::advance(FDM::SolverState& state) {
    if (!servicesBound_) {
        throw std::runtime_error("EulerianStepper requires bindServices before advance.");
    }
    if (!solveStagesBound_) {
        throw std::runtime_error(
            "EulerianStepper requires bindSolvePlan before advance.");
    }
    if (!state.bundle) {
        throw std::runtime_error("EulerianStepper::advance requires a StateBundle.");
    }
    bindState(*state.bundle);
    double dt = 0.0;
    // solver-owned registry 跨 timestep 复用；绑定顺序不变。
    operations_.clear();
    registerOperations(operations_,state,dt);
    operations_.retain(System::assignedOperations(
        runtime_,{"flow.eulerian-pressure"}));
    const Run::PlanTraceContext trace{state_->step,state_->time,&dt};
    Run::PlanExecutor::execute(solve_,operations_,&trace);
    return {true,dt,state_->time,state_->step,false,false,{}};
}

void EulerianStepper::registerOperations(
        Run::OpRegistry& operations,
        FDM::SolverState& solverState,
        double& dt) {
    lastSummary_ = {};
    operations.bind("ee.dt.compute",[&] {
        dt=std::min({
            solverState.maximumTimeStep,
            numerics_.dt.maxDeltaT,
            stableTimeStep(numerics_.dt.cfl)});
        if (services_.executionRuntime) {
            dt=services_.executionRuntime->globalMinimum(dt);
        }
        if(!std::isfinite(dt)||dt<=0.0
            ||dt>solverState.maximumTimeStep) {
            throw std::runtime_error(
                "Eulerian EulerianStepper produced an invalid time step.");
        }
    });
    operations.bind("ee.step.begin",[&] {
        if (!std::isfinite(dt) || dt <= 0.0) {
            throw std::runtime_error(
                "Eulerian Plan began a step before computing a valid dt.");
        }
        if (pendingPressureCorrection_) {
            throw std::runtime_error(
                "Eulerian Plan retained an unpublished pressure correction.");
        }
        equations_.resetTurbulenceDiagnostics(); equations_.refreshRowMap();
        applyBoundaryAndSynchronize(); system_.validateState("pre-interphase");
    });
    operations.bind("ee.interphase.compute",[&] { system_.computeInterphase(dt,workspace_.previousVelocity); });
    operations.bind("ee.sources.assemble",[&] { sourceRegistry_.assemble(system_,dt); });
    if (turbulence_.active()) {
        operations.bind("ee.turbulence.prepare",[&] {
            turbulence_.prepare(system_);
            synchronizeTurbulenceState();
        });
    }
    operations.bind("ee.sources.validate",[&] {
        equations_.validateSourceStep(dt);
    });
    operations.bind("ee.momentum.diagonal",[&] { equations_.initializeMomentumDiagonal(dt); });
    operations.bind("ee.momentum.flux",[&] { equations_.buildMomentumInterpolatedFlux(); });
    operations.bind("ee.faceFlux.canonical",[&] { assembleCanonicalPhaseFlux(); });
    operations.bind("ee.continuity.assemble",[&] { equations_.assembleContinuity(dt,lastSummary_); });
    operations.bind("ee.boundary.prepare",[&] { applyBoundaryAndSynchronize(); });
    operations.bind("ee.momentum.solve",[&] { equations_.solveMomentumPredictors(dt); });
    operations.bind("ee.interphase.correct",[&] { equations_.applySemiImplicitInterphase(dt); });
    operations.bind("ee.boundary.afterMomentum",[&] { applyBoundaryAndSynchronize(); });
    operations.bind("ee.diagonal.sync",[&] { synchronizeMomentumDiagonal(); });
    operations.bind("ee.momentum.flux.after",[&] { equations_.buildMomentumInterpolatedFlux(); });
    operations.bind("ee.faceFlux.canonical.after",[&] { assembleCanonicalPhaseFlux(); });
    operations.bind("ee.pressure.solve",[&] {
        if (pendingPressureCorrection_) {
            throw std::runtime_error(
                "Pressure solve started before the previous result was published.");
        }
        pendingPressureCorrection_=equations_.solvePressureCorrection();
    });
    operations.bind("ee.pressure.publish",[&] {
        if (!pendingPressureCorrection_) {
            throw std::runtime_error(
                "Pressure publish has no pending linear-solve result.");
        }
        equations_.setPressureCorrection(
            pendingPressureCorrection_->solution);
        lastSummary_.pressureIterations+=
            pendingPressureCorrection_->iterations;
        lastSummary_.pressureResidual=
            pendingPressureCorrection_->relativeResidual;
        ++lastSummary_.pressureCorrections;
        pendingPressureCorrection_.reset();
    });
    operations.bind("ee.pressure.sync",[&] { synchronizePressureCorrection(); });
    operations.bind("ee.phase.correct",[&] { equations_.correctAllPhases(); });
    operations.bind("ee.faceFlux.correct",[&] { equations_.correctCanonicalFaceFlux(); });
    operations.bind("ee.boundary.afterPressure",[&] { applyBoundaryAndSynchronize(); });
    operations.bind("ee.faceFlux.canonical.pressure",[&] { assembleCanonicalPhaseFlux(); });
    operations.bind("ee.energy.solve",[&] { equations_.solvePhaseEnergy(dt); });
    if (turbulence_.hasTransportEquations()) {
        operations.bind("ee.turbulence.solve",[&] {
            equations_.solveTurbulence(dt);
        });
    }
    operations.bind("ee.boundary.final",[&] { applyBoundaryAndSynchronize(); });
    operations.bind("ee.outer.validate",[&] { system_.validateState("PIMPLE outer corrector"); ++lastSummary_.outerCorrectors; });
    operations.bind("ee.step.commit",[&] {
    workspace_.commitTimeLevel(system_);
    turbulence_.commit(system_);
    lastSummary_.pressureStructureRebuilds =
        equations_.pressureReuseStatistics().structureRebuilds;
    lastSummary_.pressureLinearSolves =
        equations_.pressureReuseStatistics().solves;
    const auto turbulenceReuse =
        equations_.turbulenceReuseStatistics();
    lastSummary_.turbulenceIterations =
        equations_.turbulenceIterations();
    lastSummary_.turbulenceResidual =
        equations_.turbulenceResidual();
    lastSummary_.turbulenceStructureRebuilds =
        turbulenceReuse.structureRebuilds;
    lastSummary_.turbulenceLinearSolves = turbulenceReuse.solves;
    for (int cell = 0; cell < system_.geometry().TotalSize(); ++cell) {
        int i = 0, j = 0, k = 0;
        system_.geometry().getIJK(cell, i, j, k);
        const auto& geometry = system_.config().eulerianEulerian;
        if (!Ops::solved(
                system_.geometry(), i, j, k,
                geometry.axisymmetric, geometry.radialCoordinate)) continue;
        double emptyThickness = 1.0;
        bool representativePlane = true;
        const std::array<int,3> coordinates{i,j,k};
        for (int axis = 0; axis < 3; ++axis) {
            if (Math::isDirectionActiveIndex(axis)) continue;
            if (coordinates[(size_t)axis] != system_.geometry().NG()) {
                representativePlane = false;
                break;
            }
            const int last = system_.geometry().NG()
                + (axis == 0 ? system_.geometry().NX()
                   : axis == 1 ? system_.geometry().NY()
                               : system_.geometry().NZ()) - 1;
            const auto coordinate = [&](int index) {
                return axis == 0 ? system_.geometry().X(index,j,k)
                    : axis == 1 ? system_.geometry().Y(i,index,k)
                                : system_.geometry().Z(i,j,index);
            };
            emptyThickness = std::abs(
                coordinate(last)-coordinate(system_.geometry().NG()));
            if (!std::isfinite(emptyThickness)
                || emptyThickness <= 0.0) {
                throw std::runtime_error(
                    "Eulerian diagnostics found invalid empty thickness.");
            }
        }
        if (!representativePlane) continue;
        const double volume = 1.0 / Ops::inverseCellVolume(
            system_.geometry(), i, j, k,
            geometry.axisymmetric, geometry.radialCoordinate)
            / emptyThickness;
        double sum = 0.0;
        for (const auto& phase : system_.phases()) {
            sum += phase.primitive.alpha.values()[(size_t)cell];
            lastSummary_.totalPhaseMass += volume
                * phase.primary.phaseMass.values()[(size_t)cell];
            lastSummary_.totalPhaseEnthalpy += volume
                * phase.primary.phaseEnthalpy.values()[(size_t)cell];
        }
        lastSummary_.totalWallHeat += volume * (
            system_.sources().wallBoilingConvectiveHeat
                .values()[(size_t)cell]
            + system_.sources().wallBoilingQuenchingHeat
                .values()[(size_t)cell]
            + system_.sources().wallBoilingEvaporativeHeat
                .values()[(size_t)cell]);
        lastSummary_.maxAlphaSumError = std::max(
            lastSummary_.maxAlphaSumError, std::abs(sum - 1.0));
    }
    });
    operations.bind("ee.time.commit",[&] {
        if (pendingPressureCorrection_) {
            throw std::runtime_error(
                "Eulerian timestep cannot commit with an unpublished "
                "pressure correction.");
        }
        state_->dt=dt;
        state_->time+=dt;
        ++state_->step;
    });
}

} // namespace SF::EulerianEulerian
