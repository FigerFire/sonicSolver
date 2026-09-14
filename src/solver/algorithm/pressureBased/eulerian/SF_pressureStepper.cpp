/// @file SF_pressureStepper.cpp
/// @brief 双欧拉各相预测、共享压力校正和能量更新的数学顺序。

#include "solver/algorithm/pressureBased/eulerian/SF_pressureStepper.h"

#include "SF_phaseGravity.h"
#include "SF_phaseMRF.h"
#include "SF_phaseWallHeat.h"
#include "SF_phaseChange.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::EulerianEulerian {

namespace {
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

PressureStepper::PressureStepper(
        Physics::PhaseSystems::PhaseSystem& system,
        FDM::SolverConfig config)
    : system_(system), config_(std::move(config)),
      turbulence_(config_.turbulence),
      equations_(system_, config_.pressure.workflow, workspace_) {
    FDM::validateSolverProperties(config_.pressure.workflow);
    if (config_.pressure.workflow.type
        != FDM::SolverAlgorithm::PressureBased) {
        throw std::runtime_error(
            "EulerianEulerianPressureStepper requires type pressureBase.");
    }
    bool wallHeatEnabled = false;
    for (FDM::SourceKind kind : config_.sources.enabled) {
        switch (kind) {
            case FDM::SourceKind::Gravity:
                sourceRegistry_.add(
                    Physics::PhaseSystems::makePhaseGravitySource(
                        config_.sources.gravity));
                break;
            case FDM::SourceKind::MRF:
                sourceRegistry_.add(
                    Physics::PhaseSystems::makePhaseMRFSource(
                        config_.sources.rotating));
                break;
            case FDM::SourceKind::WallHeat:
                wallHeatEnabled = true;
                sourceRegistry_.add(
                    Physics::PhaseSystems::makePhaseWallHeatSource(
                        config_.sources.wallHeat));
                break;
        }
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
    equations_.attachTurbulence(&turbulence_);
}

void PressureStepper::registerState(State::StateBundle& state) {
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

void PressureStepper::applyBoundaryAndSynchronize() {
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

void PressureStepper::synchronizeTurbulenceState() {
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

void PressureStepper::synchronizePrimaryState() {
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

void PressureStepper::synchronizePressureCorrection() {
    if (!services_.executionRuntime) return;
    const int depth=system_.geometry().NG();
    services_.executionRuntime->finalize({
        "Eulerian pressure-correction update",
        {Execution::writeOwned("pressureCorrection")}});
    services_.executionRuntime->prepare({
        "Eulerian pressure-correction stencil",
        {Execution::readHalo("pressureCorrection",depth)}});
}

void PressureStepper::synchronizeMomentumDiagonal() {
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

void PressureStepper::assembleCanonicalPhaseFlux() {
    if (!services_.executionRuntime) return;
    std::vector<Execution::FieldAccess> fluxes;
    for (size_t phase=0; phase<workspace_.faceFlux.size(); ++phase) {
        const std::string prefix="phase"+std::to_string(phase)+".";
        fluxes.push_back(
            Execution::writeCanonicalFace(prefix+"volumeFaceFlux"));
        fluxes.push_back(
            Execution::writeCanonicalFace(prefix+"massFaceFlux"));
    }
    services_.executionRuntime->finalize({"Eulerian canonical phase flux",fluxes});
}

double PressureStepper::stableTimeStep(double cfl) {
    const double transportLimit = equations_.stableTimeStep(cfl);
    const double probeDt = std::min(
        transportLimit, config_.numerics.maxDeltaT);
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
                config_.pressure.workflow.phaseSourceCfl)),
        turbulence_.sourceTimeStep(
            system_, config_.pressure.workflow.phaseSourceCfl));
}

void PressureStepper::bindServices(FDM::SolverServices services) {
    if (servicesBound_) {
        throw std::runtime_error(
            "PressureStepper services are already bound for this lifecycle.");
    }
    services_ = services;
    servicesBound_ = true;
}

void PressureStepper::bindState(State::StateBundle& state) {
    state.validatePatches();
    if (&state.singlePatch() != &system_.geometry()) {
        throw std::runtime_error(
            "Eulerian PressureStepper requires its PhaseSystem geometry "
            "as the single StateBundle patch.");
    }
    if (state_ && state_ != &state) {
        throw std::runtime_error(
            "PressureStepper cannot switch StateBundle after stepping begins.");
    }
    state_ = &state;
    equations_.setExecutionRuntime(services_.executionRuntime);
    if (state.distributed.empty()) {
        state.registerConservativeState();
        registerState(state);
    }
    if (services_.executionRuntime) services_.executionRuntime->attachState(state);
}

FDM::StepResult PressureStepper::advance(FDM::SolverState& state) {
    if (!servicesBound_) {
        throw std::runtime_error("PressureStepper requires bindServices before advance.");
    }
    if (!state.bundle) {
        throw std::runtime_error("PressureStepper::advance requires a StateBundle.");
    }
    bindState(*state.bundle);
    double dt=std::min({
        state.maximumTimeStep,
        config_.numerics.maxDeltaT,
        stableTimeStep(config_.numerics.cfl)});
    if (services_.executionRuntime) {
        dt=services_.executionRuntime->globalMinimum(dt);
    }
    if(!std::isfinite(dt)||dt<=0.0 ||dt>state.maximumTimeStep) {
        throw std::runtime_error(
            "Eulerian PressureStepper produced an invalid time step.");
    }
    lastSummary_=stepImpl(dt);
    state_->dt=dt;
    state_->time+=dt;
    ++state_->step;
    return {true,dt,state_->time,state_->step,false,false,{}};
}

StepSummary PressureStepper::stepImpl(double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "EulerianEulerian step requires finite positive dt.");
    }
    StepSummary summary;
    equations_.resetTurbulenceDiagnostics();
    equations_.refreshRowMap();
    applyBoundaryAndSynchronize();
    system_.validateState("pre-interphase");

    for (int outer = 0;
         outer < config_.pressure.workflow.outerCorrectors; ++outer) {
        system_.computeInterphase(dt, workspace_.previousVelocity);
        sourceRegistry_.assemble(system_, dt);
        turbulence_.prepare(system_);
        synchronizeTurbulenceState();
        equations_.validateSourceStep(dt);
        equations_.initializeMomentumDiagonal(dt);
        equations_.buildMomentumInterpolatedFlux();
        assembleCanonicalPhaseFlux();
        equations_.assembleContinuity(dt, summary);
        applyBoundaryAndSynchronize();

        equations_.solveMomentumPredictors(dt);
        equations_.applySemiImplicitInterphase(dt);
        applyBoundaryAndSynchronize();
        synchronizeMomentumDiagonal();
        equations_.buildMomentumInterpolatedFlux();
        assembleCanonicalPhaseFlux();

        for (int correction = 0;
             correction < config_.pressure.workflow.pressureCorrectors;
             ++correction) {
            for (int nonOrthogonal = 0;
                 nonOrthogonal
                     <= config_.pressure.workflow.nonOrthogonalCorrectors;
                 ++nonOrthogonal) {
                const auto result = equations_.solvePressureCorrection();
                equations_.setPressureCorrection(result.solution);
                synchronizePressureCorrection();
                equations_.correctAllPhases();
                equations_.correctCanonicalFaceFlux();
                applyBoundaryAndSynchronize();
                assembleCanonicalPhaseFlux();
                summary.pressureIterations += result.iterations;
                summary.pressureResidual = result.relativeResidual;
                ++summary.pressureCorrections;
            }
        }
        equations_.solvePhaseEnergy(dt);
        equations_.solveTurbulence(dt);
        applyBoundaryAndSynchronize();
        system_.validateState("PIMPLE outer corrector");
        ++summary.outerCorrectors;
    }
    workspace_.commitTimeLevel(system_);
    turbulence_.commit(system_);
    summary.pressureStructureRebuilds =
        equations_.pressureReuseStatistics().structureRebuilds;
    summary.pressureLinearSolves =
        equations_.pressureReuseStatistics().solves;
    const auto turbulenceReuse =
        equations_.turbulenceReuseStatistics();
    summary.turbulenceIterations =
        equations_.turbulenceIterations();
    summary.turbulenceResidual =
        equations_.turbulenceResidual();
    summary.turbulenceStructureRebuilds =
        turbulenceReuse.structureRebuilds;
    summary.turbulenceLinearSolves = turbulenceReuse.solves;
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
            summary.totalPhaseMass += volume
                * phase.primary.phaseMass.values()[(size_t)cell];
            summary.totalPhaseEnthalpy += volume
                * phase.primary.phaseEnthalpy.values()[(size_t)cell];
        }
        summary.totalWallHeat += volume * (
            system_.sources().wallBoilingConvectiveHeat
                .values()[(size_t)cell]
            + system_.sources().wallBoilingQuenchingHeat
                .values()[(size_t)cell]
            + system_.sources().wallBoilingEvaporativeHeat
                .values()[(size_t)cell]);
        summary.maxAlphaSumError = std::max(
            summary.maxAlphaSumError, std::abs(sum - 1.0));
    }
    return summary;
}

} // namespace SF::EulerianEulerian
