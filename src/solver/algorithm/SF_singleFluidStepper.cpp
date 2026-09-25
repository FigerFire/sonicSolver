/// @file SF_compressible.cpp
/// @brief 注册 numerical callbacks，由 CompiledSolvePlan 驱动 timestep。
///
/// Data flow:
///   StateBundle + numerical providers -> OpRegistry
///       -> CompiledSolvePlan -> PlanExecutor -> committed state/time
///
/// 本文件不拥有 RK stage loop；stage 数学由 Time::Explicit 实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "solver/algorithm/SF_singleFluidStepper.h"
#include "solver/algorithm/SF_conservativeRHS.h"
#include "solver/algorithm/time/SF_explicit.h"
#include "solver/algorithm/SF_highOrderTrace.h"
#include "solver/algorithm/immersed/SF_constraintOps.h"
#include "solver/algorithm/immersed/SF_immersedStrategy.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/system/SF_numericalCompiler.h"
#include "SF_numericsPolicy.h"
#include "SF_physicalState.h"
#include "SF_rusanovEOS.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace SF::SolverAlgorithm {
SingleFluidStepper::SingleFluidStepper(
        FDM::SolverConfig config,
        const System::ExecutableEquationSystem& equations,
        const System::CompiledNumericalSystem& numerics,
        const System::CompiledSolvePlan& solvePlan,
        const System::RuntimeRequirements& requirements)
    : config_(std::move(config))
    , executable_(equations)
    , numerics_(numerics)
    , solve_(solvePlan)
    , runtime_(requirements)
    , boundaryApplicator_(config_.boundaries)
    , equations_(equations.equationDefinitions) {
    FDM::validateNumericsConfig(config_.numerics);
    // compiled authority 与 raw config 必须一致：runtime 只执行 compiled HOW，
    // 因此这里拒绝"input 说 A、compiled 说 B"的装配。
    if (numerics_.time.recipe.id() != config_.numerics.timeRecipe.id()
        || numerics_.dt.cfl != config_.numerics.cfl
        || numerics_.dt.maxDeltaT != config_.numerics.maxDeltaT) {
        throw std::runtime_error(
            "SingleFluidStepper received a compiled numerical system whose "
            "time/step policy disagrees with its solver configuration.");
    }
}

void SingleFluidStepper::bindSolvePlan(
        const System::CompiledSolvePlan& plan) {
    if (solvePlanBound_) {
        throw std::runtime_error(
            "SingleFluidStepper solve plan is already bound.");
    }
    if (&plan != &solve_) {
        throw std::runtime_error(
            "SingleFluidStepper received a plan other than its resolved plan.");
    }
    requiredOperations_ = System::SolvePlanner::requiredOperations(plan);
    if (requiresOperation("pressure.prepare")) {
        genericPisoCorrector_ = std::make_unique<PressureBased::Corrector>(
            config_.boundaries,config_.pressure,
            config_.numerics.idealGasGamma);
    }
    if (requiresOperation("ibm.kkt.solve")) {
        monolithicKkt_ = std::make_unique<PressureBased::MonolithicKKT>(
            config_.pressure,config_.numerics.idealGasGamma);
    }
    solvePlanBound_ = true;
}

void SingleFluidStepper::bindServices(FDM::SolverServices services) {
    if (servicesBound_) {
        throw std::runtime_error(
            "SingleFluidStepper services are already bound for this lifecycle.");
    }
    services_ = services;
    servicesBound_ = true;
}

void SingleFluidStepper::prepare(FDM::SolverState& state) {
    if (!servicesBound_ || !solvePlanBound_) {
        throw std::runtime_error(
            "SingleFluidStepper requires services and a solve plan before state preparation.");
    }
    if (!state.bundle) {
        throw std::runtime_error(
            "SingleFluidStepper::prepare requires a StateBundle.");
    }
    bindState(*state.bundle);
    realizedState_ = System::realizeState(executable_,runtime_,*state.bundle);
}

void SingleFluidStepper::bindState(State::StateBundle& state) {
    const bool firstBinding = state_ == nullptr;
    state.validatePatches();
    if (requiresOperation("explicit.stage.execute")) {
        ConservativeRHS::validateTimestepState(
            state,config_,numerics_);
    }
    if (state_ && state_ != &state) {
        throw std::runtime_error(
            "SingleFluidStepper cannot switch StateBundle after stepping begins.");
    }
    state_ = &state;
    ensureWorkspaces(state.patches);
    if (services_.executionRuntime) services_.executionRuntime->attachState(state);
    if (firstBinding && requiresOperation("explicit.stage.execute")) {
        emitConvectionContract();
    }
}

bool SingleFluidStepper::requiresOperation(
        std::string_view operation) const {
    return std::any_of(
        requiredOperations_.begin(),requiredOperations_.end(),
        [&](const System::OpId& value) { return value == operation; });
}

void SingleFluidStepper::ensureWorkspaces(const std::vector<Field*>& fields) {
    if (fields.empty()) {
        throw std::runtime_error("SingleFluidStepper requires non-empty patches.");
    }
    workspaces_.resize(fields.size());
    for (size_t index = 0; index < fields.size(); ++index) {
        if (!fields[index]) {
            throw std::runtime_error("SingleFluidStepper received null patch.");
        }
        workspaces_[index].ensureFor(*fields[index]);
    }
}

double SingleFluidStepper::physicalTime() const {
    if (!state_) throw std::runtime_error("SingleFluidStepper has no bound StateBundle.");
    return state_->time;
}

double SingleFluidStepper::timeStep() const {
    if (!state_) throw std::runtime_error("SingleFluidStepper has no bound StateBundle.");
    return state_->dt;
}

void SingleFluidStepper::prepareBoundaryState(
        const std::vector<Field*>& fields, double time, double dt) {
    ConservativeRHS::prepareBoundaryState(
        fields,numerics_.requiredHaloWidth,
        boundaryApplicator_,services_,time,dt);
}

void SingleFluidStepper::emitTimeStep() {
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

void SingleFluidStepper::emitConvectionContract() {
    if (!services_.observer || !state_ || !state_->stateModel) return;
    const auto gamma = state_->stateModel->perfectGasGamma();
    std::ostringstream detail;
    const auto& convection = System::NumericalCompiler::requireUniqueRecipe(
        numerics_,FDM::TermRole::Convection);
    detail << "reconstruction="
           << FDM::toString(convection.convection())
           << ", numericalFlux=" << FDM::toString(convection.flux())
           << ", recipe=" << FDM::toString(convection.id())
           << ", thermodynamics=PerfectGas(gamma=" << *gamma << ")"
           << ", execution=reconstructed-face flux";
    services_.observer->onSolverMessage({
        FDM::SolverMessageKind::Setup,
        "Convection contract", detail.str(),
        state_->step, state_->time, state_->dt});
}

void SingleFluidStepper::validateStateClosure(Field& field, const char* stage) {
    const Numerics::PhysicalStateSummary summary =
        Numerics::validatePhysicalState(field, stage);
    if (!services_.observer) return;
    services_.observer->onSolverMessage({FDM::SolverMessageKind::StateClosure,
        "State closure", Numerics::formatPhysicalStateSummary(stage, summary),
        state_->step, state_->time, state_->dt});
}

void SingleFluidStepper::validateDensityStateClosure(
        const std::vector<Field*>& fields, const char* stage) {
    ConservativeRHS::validateStateClosure(
        fields, *state_, services_, stage, false);
}

void SingleFluidStepper::correctTransportModel(
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

void SingleFluidStepper::bindPressureOps(
        Run::OpRegistry& operations,
        Field& field,
        double maximumTimeStep,
        PressureBased::CorrectionSummary& pressureSummary) {
    if (!genericPisoCorrector_) {
        throw std::runtime_error(
            "Generic PISO plan has no pressure operator provider.");
    }
    if (state_->patches.size() != 1 || state_->patches.front() != &field) {
        throw std::runtime_error(
            "Generic PISO capability currently requires exactly one patch.");
    }
    if (!state_->transported.empty() || services_.transportModel
        || services_.equationSystem || services_.immersed.boundary
        || services_.immersed.constraint || services_.immersed.system) {
        throw std::runtime_error(
            "Generic PISO capability supports laminar single-fluid state "
            "without transported, phase, turbulence, or IBM services.");
    }

    operations.bind("pressure.prepare",[this,&field,maximumTimeStep] {
        if (!std::isfinite(maximumTimeStep) || maximumTimeStep <= 0.0) {
            throw std::runtime_error(
                "Generic PISO requires a finite positive time-step limit.");
        }
        prepareBoundaryState({&field},state_->time,0.0);
        ensureWorkspaces({&field});
        state_->dt = state_->stateModel
            ? Numerics::RusanovEOS::deltaT(
                field,*state_->stateModel,numerics_.dt.cfl)
            : deltaT(field,numerics_.dt.cfl,
                     config_.numerics.idealGasGamma);
        if (services_.executionRuntime) {
            state_->dt = services_.executionRuntime->globalMinimum(state_->dt);
        }
        state_->dt = std::min({state_->dt,numerics_.dt.maxDeltaT,
                               maximumTimeStep});
    });
    operations.bind("momentum.assemble",[this,&field] {
        ConservativeRHS::assembleAllPatches(
            {&field},workspaces_,config_,numerics_,
            equations_,*state_,services_,
            boundaryApplicator_,state_->time);
    });
    operations.bind("momentum.solve",[this,&field] {
        Time::Explicit::forwardEuler(
            field,workspaces_.front().residual,state_->dt);
        ConservativeRHS::publishIntegratedState({&field},*state_,services_);
        validateStateClosure(field,"Euler final");
    });
    operations.bind("pressure.boundary.prepare",[this,&field] {
        prepareBoundaryState(
            {&field},state_->time+state_->dt,state_->dt);
        genericPisoCorrector_->setExecutionRuntime(
            services_.executionRuntime);
        genericPisoCorrector_->setInterfaceJumpProvider({});
    });
    operations.bind("pressure.assemble",[this,&field] {
        genericPisoCorrector_->assemble({&field},state_->dt);
    });
    operations.bind("pressure.solve",[this] {
        genericPisoCorrector_->solve();
    });
    operations.bind("pressure.update.prepare",[this] {
        genericPisoCorrector_->preparePressureUpdate();
    });
    operations.bind("velocity.correct",[this] {
        genericPisoCorrector_->correctVelocity();
    });
    operations.bind("flux.correct",[this] {
        genericPisoCorrector_->correctFlux();
    });
    operations.bind("pressure.correction.commit",[this,&pressureSummary] {
        pressureSummary = genericPisoCorrector_->commitPressureUpdate();
    });
    operations.bind("pressure.step.commit",[this,&field,&pressureSummary] {
        validateStateClosure(field,"flow correction");
        if (services_.observer) {
            std::ostringstream detail;
            detail << "iterations=" << pressureSummary.iterations
                   << ", residual=" << pressureSummary.finalResidual
                   << ", maxDiv(before)="
                   << pressureSummary.maxDivergenceBefore
                   << ", maxDiv(after)="
                   << pressureSummary.maxDivergenceAfter
                   << ", interfaceFaces=" << pressureSummary.interfaceFaces
                   << ", maxJump(target/correction)="
                   << pressureSummary.maxTargetPressureJump << "/"
                   << pressureSummary.maxCorrectionPressureJump
                   << ", HYPRE(rebuilds/solves)="
                   << pressureSummary.structureRebuilds << "/"
                   << pressureSummary.linearSolves;
            services_.observer->onSolverMessage({
                FDM::SolverMessageKind::StateClosure,
                "Flow correction",detail.str(),
                state_->step,state_->time,state_->dt});
        }
        prepareBoundaryState(
            {&field},state_->time+state_->dt,state_->dt);
        state_->time += state_->dt;
        ++state_->step;
        emitTimeStep();
    });
}

void SingleFluidStepper::bindExplicitOps(
        Run::OpRegistry& operations,
        const std::vector<Field*>& fields,
        double maximumTimeStep,
        Time::Explicit::Workspace& explicitWorkspace) {
    operations.bind("flow.step.prepare",[this,&fields,maximumTimeStep] {
        if (!std::isfinite(maximumTimeStep) || maximumTimeStep <= 0.0) {
            throw std::runtime_error(
                "SolverAlgorithm::SingleFluidStepper requires a finite positive time-step limit.");
        }
        HighOrderTrace::beginPhysicalStep(state_->step);
        if (HighOrderTrace::activeFor(state_->step)) {
            HighOrderTrace::conservative("initial Q", fields);
            HighOrderTrace::gamma(config_, *state_->stateModel);
        }
        prepareBoundaryState(fields, state_->time, 0.0);
        if (HighOrderTrace::activeFor(state_->step)) {
            HighOrderTrace::conservative("boundary-prepared Q", fields, true);
        }
    });
    operations.bind("flow.dt.compute",[this,&fields,maximumTimeStep] {
        double localDt = std::numeric_limits<double>::max();
        for (Field* field : fields) {
            if (!field) continue;
            localDt = std::min(localDt, Numerics::RusanovEOS::deltaT(
                *field, *state_->stateModel, numerics_.dt.cfl));
        }
        state_->dt = services_.executionRuntime
            ? services_.executionRuntime->globalMinimum(localDt) : localDt;
        state_->dt = std::min({state_->dt, numerics_.dt.maxDeltaT,
                               maximumTimeStep});
        if (HighOrderTrace::activeFor(state_->step)) {
            for (std::size_t patch = 0; patch < fields.size(); ++patch) {
                if (!fields[patch]) continue;
                const double legacyDt = deltaT(
                    *fields[patch], numerics_.dt.cfl,
                    config_.numerics.idealGasGamma);
                const double equationDt = Numerics::RusanovEOS::deltaT(
                    *fields[patch], *state_->stateModel, numerics_.dt.cfl);
                std::cout << std::setprecision(17)
                          << "[SF TRACE] dt patch=" << patch
                          << " legacy=" << legacyDt
                          << " stateModelRusanov=" << equationDt
                          << " active=" << state_->dt << '\n';
            }
        }
    });
    operations.bind("flow.step.begin",[this,&fields,&explicitWorkspace] {
        if (services_.equationSystem) {
            services_.equationSystem->beginStep(fields, state_->dt);
        }
        correctTransportModel(fields);
        ensureWorkspaces(fields);
        Time::Explicit::begin(
            explicitWorkspace,fields,*state_,numerics_.time.recipe,
            services_.equationSystem);
    });
    operations.bind("explicit.stage.execute",
        [this,&fields,&explicitWorkspace](
                const Run::ExecutionContext& context) {
            Time::Explicit::executeStage(
                explicitWorkspace,context.stageIndex,fields,workspaces_,*state_,
                services_.equationSystem,
                [this](const std::vector<Field*>& patches,
                       std::vector<PatchWorkspace>& workspaces,
                       double stageTime) {
                    ConservativeRHS::assembleAllPatches(
                        patches,workspaces,config_,numerics_,
                        equations_,*state_,services_,
                        boundaryApplicator_,stageTime);
                },
                [this](const std::vector<Field*>& patches) {
                    ConservativeRHS::publishIntegratedState(
                        patches,*state_,services_);
                },
                [this](const std::vector<Field*>& patches,const char* stage) {
                    validateDensityStateClosure(patches,stage);
                });
        });
}

FDM::StepResult SingleFluidStepper::advance(FDM::SolverState& state) {
    FDM::StepResult result;
    if (!servicesBound_) {
        result.message = "SingleFluidStepper requires bindServices before advance.";
        return result;
    }
    if (!solvePlanBound_) {
        result.message =
            "SingleFluidStepper requires bindSolvePlan before advance.";
        return result;
    }
    if (!state.bundle) {
        result.message = "SingleFluidStepper::advance requires a StateBundle.";
        return result;
    }
    bindState(*state.bundle);

    const auto& fields = state.bundle->patches;
    // solver-owned registry / stage storage / summary 跨 timestep 复用：
    // 普通 timestep 只 overwrite，不重新分配。
    operations_.clear();
    pressureSummary_ = PressureBased::CorrectionSummary{};
    correctionPerformed_ = false;
    correctionWroteConservative_ = false;
    correctionDetail_.clear();

    if (requiresOperation("explicit.stage.execute")) {
        bindExplicitOps(
            operations_,fields,state.maximumTimeStep,explicitWorkspace_);
    }
    if (requiresOperation("pressure.prepare")) {
        bindPressureOps(
            operations_,state.bundle->singlePatch(),state.maximumTimeStep,
            pressureSummary_);
    }
    if (requiresOperation("ibm.constraint.project")
        && services_.immersed.constraint && services_.immersed.system) {
        operations_.bind("ibm.constraint.project",[&] {
            const auto immersed = ImmersedAlgorithm::projectConstraint(
                fields,state_->time+state_->dt,state_->dt,
                services_.executionRuntime,services_.immersed);
            correctionPerformed_ = immersed.performed;
            correctionWroteConservative_ = immersed.performed;
            correctionDetail_ = immersed.detail;
        });
    }
    if (requiresOperation("ibm.kkt.solve") && monolithicKkt_
        && services_.immersed.constraint && services_.immersed.system) {
        operations_.bind("ibm.kkt.solve",[&] {
            if (fields.size() != 1 || !fields.front()) {
                throw std::runtime_error(
                    "ibm.kkt.solve currently requires one local Field.");
            }
            if (services_.executionRuntime
                && services_.executionRuntime->distributed()) {
                throw std::runtime_error(
                    "ibm.kkt.solve has no distributed ConstraintGlobalDof "
                    "provider for the pressure block.");
            }
            ImmersedAlgorithm::validateMonolithicProvider(
                *services_.immersed.system);
            services_.immersed.constraint->setExecutionRuntime(
                services_.executionRuntime);
            prepareBoundaryState(
                fields,state_->time+state_->dt,state_->dt);
            if (services_.equationSystem) {
                services_.equationSystem->preparePressureCorrection(fields);
            }
            const auto kkt = monolithicKkt_->correct(
                *fields.front(),state_->time+state_->dt,state_->dt,
                *services_.immersed.constraint);
            std::ostringstream detail;
            detail << "KKT iterations=" << kkt.iterations
                   << ", residual=" << kkt.relativeResidual
                   << ", maxDiv(before/after)="
                   << kkt.maxDivergenceBefore << "/"
                   << kkt.maxDivergenceAfter
                   << ", HYPRE(rebuilds/solves)="
                   << kkt.structureRebuilds << "/" << kkt.linearSolves
                   << ", " << kkt.immersed.detail;
            correctionPerformed_ = true;
            correctionWroteConservative_ = true;
            correctionDetail_ = detail.str();
        });
    }
    if (requiresOperation("flow.step.commit")) {
        operations_.bind("flow.step.commit",[&] {
            ConservativeRHS::publishCorrectedConservativeState(
                correctionWroteConservative_,services_);
            if (correctionPerformed_) {
                validateDensityStateClosure(fields,"flow correction");
                if (services_.observer) {
                    services_.observer->onSolverMessage({
                        FDM::SolverMessageKind::StateClosure,"Flow correction",
                        correctionDetail_,state_->step,state_->time,state_->dt});
                }
            }
            if (services_.equationSystem) {
                services_.equationSystem->commitStep(fields,state_->dt);
            }
            prepareBoundaryState(fields,state_->time+state_->dt,state_->dt);
        });
    }
    if (requiresOperation("time.commit")) {
        operations_.bind("time.commit",[this] {
            state_->time += state_->dt;
            ++state_->step;
            emitTimeStep();
        });
    }

    operations_.retain(System::assignedOperations(
        runtime_,{"flow.conservative","ibm.constraint"}));
    Run::PlanExecutor::validateBindings(solve_,operations_);
    const Run::PlanTraceContext trace{
        state_->step,state_->time,&state_->dt};
    Run::PlanExecutor::execute(solve_,operations_,&trace);
    result.accepted = true;
    result.dt = state.bundle->dt;
    result.time = state.bundle->time;
    result.step = state.bundle->step;
    return result;
}

} // namespace SF::SolverAlgorithm
