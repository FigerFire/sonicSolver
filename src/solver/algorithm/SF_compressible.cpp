/// @file SF_compressible.cpp
/// @brief 单流体可压缩算法的边界、RHS、时间步与校正流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "solver/algorithm/SF_compressible.h"
#include "solver/algorithm/SF_densityBasedRHS.h"
#include "solver/algorithm/time/SF_explicit.h"
#include "solver/algorithm/SF_highOrderTrace.h"
#include "SF_numericsPolicy.h"
#include "SF_physicalState.h"
#include "SF_rusanovEOS.h"
#include "solver/algorithm/SF_solverAlgorithm.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace SF::SolverAlgorithm {
CompressibleAlgorithm::CompressibleAlgorithm(
        FDM::SolverConfig config,
        const System::ResolvedSimulationSystem& system)
    : config_(std::move(config))
    , resolved_(system)
    , boundaryApplicator_(config_.boundaries)
    , equations_(system.equationDefinitions)
    , flowAlgorithm_(FDM::makeFlowAlgorithm(config_)) {
    FDM::validateNumericsConfig(config_.numerics);
}

void CompressibleAlgorithm::bindSolveStages(
        const std::vector<FDM::SolveStage>& stages) {
    if (solveStagesBound_) {
        throw std::runtime_error(
            "CompressibleAlgorithm solve stages are already bound.");
    }
    if (stages.size() != resolved_.solveBlocks.size()) {
        throw std::runtime_error(
            "CompressibleAlgorithm workflow does not match the resolved solve blocks.");
    }
    for (size_t index = 0; index < stages.size(); ++index) {
        const auto& stage = stages[index];
        const auto& block = resolved_.solveBlocks[index];
        if (stage.id != block.id || stage.equations != block.equations
            || stage.constraints != block.constraints
            || stage.strategy != block.strategyKind) {
            throw std::runtime_error(
                "CompressibleAlgorithm workflow stage '" + stage.id
                + "' differs from resolved solve block '" + block.id + "'.");
        }
        densitySolveBlockActive_ = densitySolveBlockActive_
            || stage.strategy == FDM::SolveStrategyKind::ExplicitTimeIntegration;
        pressurePredictorActive_ = pressurePredictorActive_
            || stage.strategy == FDM::SolveStrategyKind::SegregatedPredictor;
        pressureCorrectionActive_ = pressureCorrectionActive_
            || stage.strategy == FDM::SolveStrategyKind::PressureCorrection;
    }
    if (config_.numerics.solver == FDM::SolverAlgorithm::DensityBased) {
        if (!densitySolveBlockActive_) {
            throw std::runtime_error(
                "Density execution requires an explicit-time solve strategy.");
        }
    } else if (!pressurePredictorActive_ || !pressureCorrectionActive_) {
        throw std::runtime_error(
            "Pressure execution requires typed predictor and correction strategies.");
    }
    solveStagesBound_ = true;
}

void CompressibleAlgorithm::bindServices(FDM::SolverServices services) {
    if (servicesBound_) {
        throw std::runtime_error(
            "CompressibleAlgorithm services are already bound for this lifecycle.");
    }
    services_ = services;
    servicesBound_ = true;
}

void CompressibleAlgorithm::prepare(FDM::SolverState& state) {
    if (!servicesBound_ || !solveStagesBound_) {
        throw std::runtime_error(
            "CompressibleAlgorithm requires services and solve stages before state preparation.");
    }
    if (!state.bundle) {
        throw std::runtime_error(
            "CompressibleAlgorithm::prepare requires a StateBundle.");
    }
    bindState(*state.bundle);
    realizedState_ = System::realizeState(resolved_,*state.bundle);
}

void CompressibleAlgorithm::bindState(State::StateBundle& state) {
    const bool firstBinding = state_ == nullptr;
    state.validatePatches();
    if (config_.numerics.solver == FDM::SolverAlgorithm::DensityBased) {
        DensityBasedRHS::validateTimestepState(state, config_);
    }
    if (state_ && state_ != &state) {
        throw std::runtime_error(
            "CompressibleAlgorithm cannot switch StateBundle after stepping begins.");
    }
    state_ = &state;
    ensureWorkspaces(state.patches);
    if (services_.executionRuntime) services_.executionRuntime->attachState(state);
    if (firstBinding && config_.numerics.solver == FDM::SolverAlgorithm::DensityBased) {
        emitConvectionContract();
    }
}

void CompressibleAlgorithm::ensureWorkspaces(const std::vector<Field*>& fields) {
    if (fields.empty()) {
        throw std::runtime_error("CompressibleAlgorithm requires non-empty patches.");
    }
    workspaces_.resize(fields.size());
    for (size_t index = 0; index < fields.size(); ++index) {
        if (!fields[index]) {
            throw std::runtime_error("CompressibleAlgorithm received null patch.");
        }
        workspaces_[index].ensureFor(*fields[index]);
    }
}

double CompressibleAlgorithm::physicalTime() const {
    if (!state_) throw std::runtime_error("CompressibleAlgorithm has no bound StateBundle.");
    return state_->time;
}

double CompressibleAlgorithm::timeStep() const {
    if (!state_) throw std::runtime_error("CompressibleAlgorithm has no bound StateBundle.");
    return state_->dt;
}

void CompressibleAlgorithm::prepareBoundaryState(
        const std::vector<Field*>& fields, double time, double dt) {
    DensityBasedRHS::prepareBoundaryState(
        fields, config_, boundaryApplicator_, services_, time, dt);
}

void CompressibleAlgorithm::emitTimeStep() {
    if (!services_.observer) return;

    std::ostringstream detail;
    detail << "step=" << state_->step
           << ", time=" << state_->time
           << ", dt=" << state_->dt;

    FDM::SolverMessage message;
    message.kind = FDM::SolverMessageKind::TimeStep;
    message.topic = "Time step";
    message.detail = detail.str();
    message.step = state_->step;
    message.time = state_->time;
    message.dt = state_->dt;
    services_.observer->onSolverMessage(message);
}

void CompressibleAlgorithm::emitConvectionContract() {
    if (!services_.observer || !state_ || !state_->equations) return;
    const auto gamma = state_->equations->perfectGasGamma();
    std::ostringstream detail;
    detail << "reconstruction="
           << FDM::toString(config_.numerics.convection)
           << ", numericalFlux=" << FDM::toString(config_.numerics.flux)
           << ", thermodynamics=PerfectGas(gamma=" << *gamma << ")"
           << ", execution=reconstructed-face flux";
    services_.observer->onSolverMessage({
        FDM::SolverMessageKind::Setup,
        "Convection contract", detail.str(),
        state_->step, state_->time, state_->dt});
}

void CompressibleAlgorithm::validateStateClosure(Field& field, const char* stage) {
    const Numerics::PhysicalStateSummary summary =
        Numerics::validatePhysicalState(field, stage);
    if (!services_.observer) return;
    services_.observer->onSolverMessage({FDM::SolverMessageKind::StateClosure,
        "State closure", Numerics::formatPhysicalStateSummary(stage, summary),
        state_->step, state_->time, state_->dt});
}

void CompressibleAlgorithm::validateDensityStateClosure(
        const std::vector<Field*>& fields, const char* stage) {
    DensityBasedRHS::validateStateClosure(
        fields, *state_, services_, stage, false);
}

void CompressibleAlgorithm::correctTransportModel(
        const std::vector<Field*>& fields) {
    if (!services_.transportModel) return;
    const auto reads = services_.transportModel->distributedReadFields();
    const auto writes = services_.transportModel->distributedWriteFields();
    const int haloDepth = services_.transportModel->distributedHaloDepth();
    auto writeContract = [&](const char* name) {
        std::vector<Execution::FieldAccess> accesses;
        for (const auto& fieldName : writes) {
            accesses.push_back(Execution::writeOwned(fieldName));
        }
        if (services_.executionRuntime && !accesses.empty()) {
            services_.executionRuntime->finalize({name, accesses});
        }
    };
    auto readContract = [&](const char* name) {
        if (reads.empty()) return;
        if (haloDepth <= 0) {
            throw std::runtime_error(
                "Transport model declared distributed reads without "
                "a positive stencil halo depth.");
        }
        std::vector<Execution::FieldAccess> accesses;
        for (const auto& fieldName : reads) {
            accesses.push_back(Execution::readHalo(fieldName, haloDepth));
        }
        if (services_.executionRuntime) {
            services_.executionRuntime->prepare({name, accesses});
        }
    };
    for (Field* field : fields) {
        if (!field) continue;
        services_.transportModel->applyBoundary(*field);
        writeContract("transport boundary update");
        readContract("transport correction stencil");
        services_.transportModel->correct(*field, state_->dt);
        services_.transportModel->applyBoundary(*field);
        writeContract("transport corrected state");
        readContract("transport diffusion stencil");
    }
}

void CompressibleAlgorithm::stepPressure(Field& field, double maximumTimeStep) {
    if (!std::isfinite(maximumTimeStep) || maximumTimeStep <= 0.0) {
        throw std::runtime_error(
            "SolverAlgorithm::CompressibleAlgorithm requires a finite positive time-step limit.");
    }
    prepareBoundaryState({&field}, state_->time, 0.0);
    ensureWorkspaces({&field});

    state_->dt = state_->equations
        ? Numerics::RusanovEOS::deltaT(field, *state_->equations, config_.numerics.cfl)
        : deltaT(field, config_.numerics.cfl,
                 config_.numerics.idealGasGamma);
    if (services_.executionRuntime) state_->dt = services_.executionRuntime->globalMinimum(state_->dt);
    state_->dt = std::min({state_->dt, config_.numerics.maxDeltaT, maximumTimeStep});

    if (services_.equationSystem) services_.equationSystem->beginStep(field, state_->dt);

    correctTransportModel({&field});

    auto assembleRHS = [&](Field& f, double stageTime) {
        DensityBasedRHS::assembleAllPatches(
            {&f}, workspaces_, config_, equations_, *state_, services_, boundaryApplicator_,
            stageTime);
    };

    auto postStage = [&](Field& f, const char* stage) {
        DensityBasedRHS::publishIntegratedState({&f}, *state_, services_);
        validateStateClosure(f, stage);
    };
    ddtDispatch(field, workspaces_.front().residual,
                state_->time, state_->dt, config_.numerics.time,
                assembleRHS, postStage, &state_->transported);

    FDM::FlowAlgorithmContext flowContext;
    flowContext.fields = {&field};
    flowContext.targetTime = state_->time + state_->dt;
    flowContext.executionRuntime = services_.executionRuntime;
    flowContext.equationSystem = services_.equationSystem;
    flowContext.immersed = services_.immersed;
    flowContext.prepareBoundaryState = [this, &field]() {
        prepareBoundaryState({&field}, state_->time + state_->dt, state_->dt);
    };
    const FDM::FlowAlgorithmResult flowResult =
        flowAlgorithm_->correct(flowContext, state_->dt);
    const bool immersedForcingCorrected =
        flowResult.performedCorrection && services_.immersed.constraint != nullptr;
    if (immersedForcingCorrected && services_.executionRuntime) {
        services_.executionRuntime->finalize({
            "immersed forcing corrected state",
            {Execution::writeOwned("conservative")}});
    }
    if (flowResult.performedCorrection) {
        validateStateClosure(field, "flow correction");
        if (services_.observer) {
            services_.observer->onSolverMessage({FDM::SolverMessageKind::StateClosure,
                                        "Flow correction",
                                        flowResult.detail,
                                        state_->step, state_->time, state_->dt});
        }
    }

    if (services_.equationSystem) services_.equationSystem->commitStep(field, state_->dt);

    // 对没有约束 IBM 校正的普通路径，保持既有的时间步尾边界闭合。
    // 强制 IBM 路径的下一个消费者是下一时间步开头的 boundary state，
    // 因而在那里进行唯一的 ReadHalo；不得在此再隐式重复同步。
    if (!immersedForcingCorrected) {
        prepareBoundaryState({&field}, state_->time + state_->dt, state_->dt);
    }
    state_->time += state_->dt;
    ++state_->step;
    emitTimeStep();
}

void CompressibleAlgorithm::finishDensityStep(
        const std::vector<Field*>& fields,
        const FDM::FlowAlgorithmResult& flowResult) {
    // All patches share one mathematical lifecycle.  DensityBasedRHS/Time
    // already completed spatial stage barriers; Runtime remains responsible
    // for distributed publication.  Commit model history before exposing the
    // now read-ready timestep to observers and output.
    DensityBasedRHS::publishCorrectedConservativeState(
        flowResult.wroteConservativeState, services_);
    if (flowResult.performedCorrection) {
        validateDensityStateClosure(fields, "flow correction");
        if (services_.observer) {
            services_.observer->onSolverMessage({
                FDM::SolverMessageKind::StateClosure,
                "Flow correction", flowResult.detail,
                state_->step, state_->time, state_->dt});
        }
    }
    if (services_.equationSystem) {
        services_.equationSystem->commitStep(fields, state_->dt);
    }
    prepareBoundaryState(fields, state_->time + state_->dt, state_->dt);
    state_->time += state_->dt;
    ++state_->step;
    emitTimeStep();
}

void CompressibleAlgorithm::stepDensity(
        const std::vector<Field*>& fields, double maximumTimeStep) {
    if (!std::isfinite(maximumTimeStep) || maximumTimeStep <= 0.0) {
        throw std::runtime_error(
            "SolverAlgorithm::CompressibleAlgorithm requires a finite positive time-step limit.");
    }
    HighOrderTrace::beginPhysicalStep(state_->step);
    if (HighOrderTrace::activeFor(state_->step)) {
        HighOrderTrace::conservative("initial Q", fields);
        HighOrderTrace::gamma(config_, *state_->equations);
    }
    prepareBoundaryState(fields, state_->time, 0.0);
    if (HighOrderTrace::activeFor(state_->step)) {
        HighOrderTrace::conservative(
            "boundary-prepared Q", fields, true);
    }

    double localDt = std::numeric_limits<double>::max();
    for (Field* field : fields) {
        if (!field) continue;
        localDt = std::min(localDt, Numerics::RusanovEOS::deltaT(
            *field, *state_->equations, config_.numerics.cfl));
    }
    state_->dt = services_.executionRuntime
        ? services_.executionRuntime->globalMinimum(localDt) : localDt;
    state_->dt = std::min({state_->dt, config_.numerics.maxDeltaT,
                           maximumTimeStep});
    if (HighOrderTrace::activeFor(state_->step)) {
        for (std::size_t patch = 0; patch < fields.size(); ++patch) {
            if (!fields[patch]) continue;
            const double legacyDt = deltaT(
                *fields[patch], config_.numerics.cfl,
                config_.numerics.idealGasGamma);
            const double equationDt = Numerics::RusanovEOS::deltaT(
                *fields[patch], *state_->equations, config_.numerics.cfl);
            std::cout << std::setprecision(17)
                      << "[SF TRACE] dt patch=" << patch
                      << " legacy=" << legacyDt
                      << " equationSetRusanov=" << equationDt
                      << " active=" << state_->dt << '\n';
        }
    }

    if (services_.equationSystem) {
        services_.equationSystem->beginStep(fields, state_->dt);
    }
    correctTransportModel(fields);

    ensureWorkspaces(fields);
    Time::Explicit::advance(
        fields, workspaces_, *state_, config_.numerics.time, services_.equationSystem,
        [this](const std::vector<Field*>& patches,
               std::vector<PatchWorkspace>& workspaces, double stageTime) {
            DensityBasedRHS::assembleAllPatches(
                patches, workspaces, config_, equations_, *state_, services_,
                boundaryApplicator_, stageTime);
        },
        [this](const std::vector<Field*>& patches) {
            DensityBasedRHS::publishIntegratedState(
                patches, *state_, services_);
        },
        [this](const std::vector<Field*>& patches, const char* stage) {
            validateDensityStateClosure(patches, stage);
        });

    FDM::FlowAlgorithmContext flowContext;
    flowContext.fields = fields;
    flowContext.targetTime = state_->time + state_->dt;
    flowContext.executionRuntime = services_.executionRuntime;
    flowContext.equationSystem = services_.equationSystem;
    flowContext.immersed = services_.immersed;
    flowContext.prepareBoundaryState = [this, &fields]() {
        prepareBoundaryState(fields, state_->time + state_->dt, state_->dt);
    };
    finishDensityStep(fields, flowAlgorithm_->correct(flowContext, state_->dt));
}

FDM::StepResult CompressibleAlgorithm::advance(FDM::SolverState& state) {
    FDM::StepResult result;
    if (!servicesBound_) {
        result.message = "CompressibleAlgorithm requires bindServices before advance.";
        return result;
    }
    if (!solveStagesBound_) {
        result.message =
            "CompressibleAlgorithm requires bindSolveStages before advance.";
        return result;
    }
    if (!state.bundle) {
        result.message = "CompressibleAlgorithm::advance requires a StateBundle.";
        return result;
    }
    bindState(*state.bundle);
    if (densitySolveBlockActive_) {
        stepDensity(state.bundle->patches, state.maximumTimeStep);
    } else if (pressurePredictorActive_ && pressureCorrectionActive_) {
        stepPressure(state.bundle->singlePatch(), state.maximumTimeStep);
    } else {
        result.message =
            "CompressibleAlgorithm has no executable bound solve block.";
        return result;
    }
    result.accepted = true;
    result.dt = state.bundle->dt;
    result.time = state.bundle->time;
    result.step = state.bundle->step;
    return result;
}

} // namespace SF::SolverAlgorithm
