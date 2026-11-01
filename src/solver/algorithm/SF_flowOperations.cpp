#include "core/system/SF_operationIds.h"
/// @file SF_singleFluidStepper.cpp
/// @brief 原 flow/RK/pressure workspace 的稳定 owner；只绑定冻结的单操作内核。
///
/// Data flow:
///   StateBundle + numerical providers -> OpRegistry
///       -> CompiledSolvePlan -> PlanExecutor -> committed state/time
///
/// 本文件不拥有 RK stage loop；stage 数学由 Time::Explicit 实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "solver/algorithm/SF_flowOperations.h"
#include "solver/algorithm/SF_pressureProviderBinding.h"
#include "solver/algorithm/SF_conservativeRHS.h"
#include "solver/algorithm/time/SF_explicit.h"
#include "solver/algorithm/SF_highOrderTrace.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_stateRealization.h"
#include "core/system/SF_formula.h"
#include "core/system/SF_scheduleIds.h"
#include "solver/system/SF_numericalCompiler.h"
#include "SF_numericsPolicy.h"
#include "solver/system/SF_scalarMethod.h"
#include "SF_physicalState.h"
#include "SF_rusanovEOS.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace SF::SolverAlgorithm {
FlowOperations::FlowOperations(
        FDM::BoundaryConfig boundaries,FDM::PressureCorrectionConfig pressure,double pressureGamma,
        const System::ExecutableEquationSystem& equations,
        const System::CompiledNumericalSystem& numerics,
        const System::CompiledSolvePlan& solvePlan,
        const System::RuntimeRequirements& requirements,
        State::StateBundle& state, const FDM::SolverServices& services,
        System::StateRealization& realized, const double& maximumTimeStep)
    : pressureConfig_(std::move(pressure)),pressureGamma_(pressureGamma)
    , executable_(equations)
    , numerics_(numerics)
    , solve_(solvePlan)
    , runtime_(requirements)
    , boundaryApplicator_(std::move(boundaries))
    , services_(services), state_(&state), realizedState_(realized)
    , maximumTimeStep_(maximumTimeStep) {
    if(!std::isfinite(numerics_.dt.cfl) || numerics_.dt.cfl<=0 || !std::isfinite(numerics_.dt.maxDeltaT) || numerics_.dt.maxDeltaT<=0)
        throw std::runtime_error("Flow execution requires a finite positive compiled dt policy");
    if (System::requiredPackedComponents(executable_.state,solve_.compiledProgram,"conservative")>0) {
        equations_=std::make_unique<SF::Equation::Compressible::System>();
    }
    if (!equations_ && (System::requiresOperation(runtime_,System::OpIds::FlowStepPrepare)
        || System::requiresOperation(runtime_,"flow.stage.rhs")
        || System::requiresOperation(runtime_,System::OpIds::ExplicitStageExecute)))
        throw std::runtime_error("Requested conservative operations require declared packed conservative STATE backing.");

}

void FlowOperations::bindSolvePlan(
        const System::CompiledSolvePlan& plan) {
    if (solvePlanBound_) {
        throw std::runtime_error(
            "FlowOperations solve plan is already bound.");
    }
    if (&plan != &solve_) {
        throw std::runtime_error(
            "FlowOperations received a plan other than its resolved plan.");
    }
    requiredOperations_ = System::SolvePlanner::requiredOperations(plan);
    if ((requiresOperation(System::OpIds::ExplicitStageExecute) || requiresOperation("flow.stage.rhs"))) {
        const auto reason=System::validateCompiledConservativeStage(executable_,plan);
        if (!reason.empty()) throw std::runtime_error(reason);
    }
    if (requiresOperation(System::OpIds::PressurePrepare)) {
        auto providers=bindPressureProvider({executable_,numerics_,solve_,runtime_},{
            boundaryApplicator_.definitions(),pressureConfig_,pressureGamma_});
        pressureOperators_=std::move(providers.constant);
        conservativePressureCorrector_=std::move(providers.conservative);
    }
    solvePlanBound_ = true;
}

void FlowOperations::prepare() {
    bindSolvePlan(solve_);
    bindState(*state_);
    if (conservativePressureCorrector_) {
        if (state_->patches.size()!=1)
            throw std::runtime_error("Conservative pressure provider requires one patch.");
        conservativePressureCorrector_->bindStateViews(realizedState_,*state_->patches.front());
    }
    if (pressureOperators_) {
        if (!services_.executionRuntime)
            throw std::runtime_error("Pressure operations require a bound execution runtime.");
        pressureOperators_->bind(realizedState_,*services_.executionRuntime);
    }
}

void FlowOperations::bindTemporalStateViews() {
    std::size_t index=0;
    for (const auto& demand:solve_.compiledProgram.stateViews) {
        if (demand.kind!=System::StateViewKind::OldTime
            && demand.kind!=System::StateViewKind::Stage) continue;
        if (demand.producer) continue; // Generic producer binding owns derived Stage views.
        if (std::none_of(solve_.compiledProgram.steps.begin(),solve_.compiledProgram.steps.end(),[&](const auto& call) {
            return call.temporalRhs=="flow.stage.rhs" && call.temporalResidual && call.target.symbol==demand.symbol;
        })) continue;
        for (std::size_t patch=0;patch<state_->patches.size();++patch) {
            if (index==temporalViews_.size())
                temporalViews_.push_back(std::make_unique<State::DistributedFieldView>());
            auto& view=*temporalViews_[index++];
            Field* geometry=state_->patches[patch];
            view.name=demand.storage;view.blockId=(int)patch;view.geometry=geometry;
            view.components=demand.components;view.haloDepth=geometry->NG();
            view.exchange=State::ExchangeKind::None;
            // This view is owned by this implementation; an owning capture would form a cycle.
            view.read=[this,geometry,patch,demand,symbolPointer=&executable_.state.at(demand.symbol)](int cell,int component) {
                if (component<0 || component>=demand.components)
                    throw std::out_of_range("Temporal STATE view component is out of range.");
                const auto& symbol=*symbolPointer;
                if (symbol.storageBinding==System::StorageBinding::NamedDistributed) {
                    auto* registry=services_.equationSystem ? services_.equationSystem->variables(*geometry) : &state_->transported;
                    if (!registry) throw std::runtime_error("Temporal scalar has no active registry.");
                    const auto& variables=registry->variables();
                    const auto found=std::find_if(variables.begin(),variables.end(),[&](const auto& value) {return value.descriptor.name==symbol.storageKey;});
                    if (found==variables.end()) throw std::runtime_error("Temporal scalar missing registered storage: "+symbol.storageKey);
                    if (demand.kind==System::StateViewKind::OldTime)
                        return explicitWorkspace_.patches.at(patch).registered.at(found-variables.begin()).q0.at((size_t)cell);
                    if (!explicitWorkspace_.active || explicitWorkspace_.nextStage!=demand.stage)
                        throw std::runtime_error("Scalar Stage STATE view is outside its active stage lifetime.");
                    return found->value->values().at((size_t)cell);
                }
                if (demand.kind==System::StateViewKind::OldTime)
                    return explicitWorkspace_.patches.at(patch).q0.at(
                        (size_t)cell*geometry->NVar()+demand.componentOffset+component);
                if (!explicitWorkspace_.active || explicitWorkspace_.nextStage!=demand.stage)
                    throw std::runtime_error("Stage STATE view is outside its active stage lifetime.");
                int i=0,j=0,k=0;geometry->getIJK(cell,i,j,k);
                return (*geometry)(i,j,k,demand.componentOffset+component);
            };
            view.write=[](int,int,double) {
                throw std::runtime_error("Temporal STATE view is published by the time provider.");
            };
            realizedState_.bindView(demand.symbol,demand.kind,view,demand.stage);
        }
    }
}

void FlowOperations::bindState(State::StateBundle& state) {
    const bool firstBinding = state_ == nullptr;
    state.validatePatches();
    if ((requiresOperation(System::OpIds::ExplicitStageExecute) || requiresOperation("flow.stage.rhs"))) {
        ConservativeRHS::validateTimestepState(
            state,numerics_);
    }
    if (state_ && state_ != &state) {
        throw std::runtime_error(
            "FlowOperations cannot switch StateBundle after stepping begins.");
    }
    state_ = &state;
    if (!pressureOperators_) ensureWorkspaces(state.patches);
    if (firstBinding && (requiresOperation(System::OpIds::ExplicitStageExecute) || requiresOperation("flow.stage.rhs"))) {
        emitConvectionContract();
    }
}

bool FlowOperations::requiresOperation(
        std::string_view operation) const {
    return std::any_of(
        requiredOperations_.begin(),requiredOperations_.end(),
        [&](const System::OpId& value) { return value == operation; });
}

void FlowOperations::ensureWorkspaces(const std::vector<Field*>& fields) {
    if (fields.empty()) {
        throw std::runtime_error("FlowOperations requires non-empty patches.");
    }
    workspaces_.resize(fields.size());
    for (size_t index = 0; index < fields.size(); ++index) {
        if (!fields[index]) {
            throw std::runtime_error("FlowOperations received null patch.");
        }
        workspaces_[index].ensureFor(*fields[index]);
    }
}

void FlowOperations::prepareBoundaryState(
        const std::vector<Field*>& fields, double time, double dt) {
    ConservativeRHS::prepareBoundaryState(
        fields,numerics_.requiredHaloWidth,
        boundaryApplicator_,services_,time,dt);
}

void FlowOperations::emitConvectionContract() {
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

void FlowOperations::validateStateClosure(Field& field, const char* stage) {
    const Numerics::PhysicalStateSummary summary =
        Numerics::validatePhysicalState(field, stage);
    if (!services_.observer) return;
    services_.observer->onSolverMessage({FDM::SolverMessageKind::StateClosure,
        "State closure", Numerics::formatPhysicalStateSummary(stage, summary),
        state_->step, state_->time, state_->dt});
}

void FlowOperations::validateDensityStateClosure(
        const std::vector<Field*>& fields, const char* stage) {
    ConservativeRHS::validateStateClosure(
        fields, *state_, services_, stage, false);
}

void FlowOperations::bindPressureOps(
        Run::OpRegistry& operations,
        Field& field,
        const double& maximumTimeStep,
        PressureBased::CorrectionSummary& pressureSummary) {
    if (!conservativePressureCorrector_) {
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

    operations.bind(System::OpIds::PressurePrepare,"flow.conservative",[owner=shared_from_this(),this,&field,&maximumTimeStep] {
        resetStep();
        if (!std::isfinite(maximumTimeStep) || maximumTimeStep <= 0.0) {
            throw std::runtime_error(
                "Generic PISO requires a finite positive time-step limit.");
        }
        prepareBoundaryState({&field},state_->time,0.0);
        ensureWorkspaces({&field});
        state_->validateThermodynamicBinding();
        state_->dt=Numerics::RusanovEOS::deltaT(field,*state_->stateModel,numerics_.dt.cfl);
        if (services_.executionRuntime) {
            state_->dt = services_.executionRuntime->globalMinimum(state_->dt);
        }
        state_->dt = std::min({state_->dt,numerics_.dt.maxDeltaT,
                               maximumTimeStep});
    });
    operations.bind(System::OpIds::MomentumAssemble,"flow.conservative",[owner=shared_from_this(),this,&field] {
        ConservativeRHS::assembleAllPatches(
            {&field},workspaces_,numerics_,
            *equations_,*state_,services_,
            boundaryApplicator_,state_->time);
    });
    operations.bind(System::OpIds::MomentumSolve,"flow.conservative",[owner=shared_from_this(),this,&field] {
        Time::Explicit::forwardEuler(
            field,workspaces_.front().residual,state_->dt);
        ConservativeRHS::publishIntegratedState({&field},*state_,services_);
        validateStateClosure(field,"Euler final");
    });
    operations.bind(System::OpIds::PressureBoundaryPrepare,"flow.conservative",[owner=shared_from_this(),this,&field] {
        prepareBoundaryState(
            {&field},state_->time+state_->dt,state_->dt);
        conservativePressureCorrector_->setExecutionRuntime(
            services_.executionRuntime);
        conservativePressureCorrector_->setInterfaceJumpProvider({});
    });
    operations.bind(System::OpIds::PressureAssemble,"flow.conservative",[owner=shared_from_this(),this,&field] {
        conservativePressureCorrector_->assemble({&field},state_->dt);
    });
    operations.bind(System::OpIds::PressureSolve,"flow.conservative",[owner=shared_from_this(),this] {
        conservativePressureCorrector_->solve();
    });
    operations.bind(System::OpIds::PressureUpdatePrepare,"flow.conservative",[owner=shared_from_this(),this] {
        conservativePressureCorrector_->preparePressureUpdate();
    });
    operations.bind(System::OpIds::VelocityCorrect,"flow.conservative",[owner=shared_from_this(),this] {
        conservativePressureCorrector_->correctVelocity();
    });
    operations.bind(System::OpIds::FluxCorrect,"flow.conservative",[owner=shared_from_this(),this] {
        conservativePressureCorrector_->correctFlux();
    });
    operations.bind(System::OpIds::PressureCorrectionCommit,"flow.conservative",[owner=shared_from_this(),this,&pressureSummary] {
        pressureSummary = conservativePressureCorrector_->commitPressureUpdate();
    });
    operations.bind(System::OpIds::PressureStepCommit,"flow.conservative",[owner=shared_from_this(),this,&field,&pressureSummary] {
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
    });
}

void FlowOperations::bindConstantPressureOps(
        Run::OpRegistry& operations,const double& maximumTimeStep) {
    if (!pressureOperators_ || state_->patches.size()!=1) {
        throw std::runtime_error(
            "Constant-density pressure operations require one local patch.");
    }
    operations.bind(System::OpIds::PressurePrepare,"flow.pressure-operators",[owner=shared_from_this(),this,&maximumTimeStep] {
        resetStep();
        state_->dt=pressureOperators_->prepare(maximumTimeStep);
    });
    operations.bind(System::OpIds::PressureStepBegin,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->beginFixedTimeStep();
    });
    operations.bind(System::OpIds::PressureIterationBegin,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->beginIteration();
    });
    operations.bind(System::OpIds::MomentumAssemble,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->assembleMomentum();
        pressureOperators_->traceState("momentum",state_->step,state_->time);
    });
    operations.bind(System::OpIds::MomentumSolve,"flow.pressure-operators",[owner=shared_from_this(),this](const Run::ExecutionContext& context) {
        if (!context.target)
            throw std::runtime_error("Momentum numerical leaf lacks typed target storage.");
        pressureOperators_->solveMomentum(context.target->kind);
        if (context.target->kind==System::TargetKind::Physical)
            state_->distributed.markModified("velocity");
    });
    operations.bind(System::OpIds::PressureBoundaryPrepare,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->preparePressureBoundary();
    });
    operations.bind(System::OpIds::PressureAssemble,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->assemblePressure();
    });
    operations.bind(System::OpIds::PressureSolve,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->solvePressure();
    });
    operations.bind(System::OpIds::PressureUpdatePrepare,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->preparePressureUpdate();
        state_->distributed.markModified("pressure");
    });
    operations.bind(System::OpIds::VelocityCorrect,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->correctVelocity();
        state_->distributed.markModified("velocity");
    });
    operations.bind(System::OpIds::FluxCorrect,"flow.rhie-chow",[owner=shared_from_this(),this] {
        pressureOperators_->traceState("flux-before",state_->step,state_->time+state_->dt);
        pressureOperators_->correctFlux();
        pressureOperators_->traceState("flux-after",state_->step,state_->time+state_->dt);
    });
    operations.bind(System::OpIds::PressureCorrectionCommit,"flow.pressure-operators",[owner=shared_from_this(),this] {
        const auto summary=pressureOperators_->commitCorrection();
        if (services_.observer && !pressureOperators_->fixedTimeActive()) {
            std::ostringstream detail;
            detail << "iterations=" << summary.iterations
                   << ", relativeResidual=" << summary.relativeResidual
                   << ", maxDiv(before/after)="
                   << summary.maxDivergenceBefore << "/"
                   << summary.maxDivergenceAfter;
            services_.observer->onSolverMessage({
                FDM::SolverMessageKind::StateClosure,
                "Pressure continuity",detail.str(),
                state_->step,state_->time,state_->dt});
        }
    });
    operations.bind(System::OpIds::PressureRelaxationApply,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->applyRelaxation();
        state_->distributed.markModified("pressure");
        state_->distributed.markModified("velocity");
    });
    operations.bind(System::OpIds::PressureFluxConsistencyRestore,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->restoreFluxConsistency();
    });
    operations.bind(System::OpIds::PressureConvergenceEvaluate,"flow.pressure-operators",
                    [owner=shared_from_this(),this](const Run::ExecutionContext& context) {
        pressureOperators_->evaluateConvergence();
        context.signals->publish(System::kPressureOuterConvergedSignal,
                                 pressureOperators_->converged());
        if (services_.observer) {
            const auto summary=pressureOperators_->commitCorrection();
            std::ostringstream continuity;
            continuity << "iterations=" << summary.iterations
                       << ", relativeResidual=" << summary.relativeResidual
                       << ", maxDiv(before/after)="
                       << summary.maxDivergenceBefore << "/"
                       << summary.maxDivergenceAfter;
            services_.observer->onSolverMessage({
                FDM::SolverMessageKind::StateClosure,
                "Pressure continuity",continuity.str(),
                state_->step,state_->time,state_->dt});
            std::ostringstream detail;
            detail << "maxDelta=" << pressureOperators_->iterationDelta()
                   << ", velocityDelta="
                   << pressureOperators_->velocityIterationDelta()
                   << ", pressureDelta="
                   << pressureOperators_->pressureIterationDelta()
                   << ", maxDiv=" << pressureOperators_->continuityDefect()
                   << ", candidateDiv="
                   << pressureOperators_->candidateContinuityDefect()
                   << ", fluxDelta="
                   << pressureOperators_->fluxIterationDelta()
                   << ", converged="
                   << (pressureOperators_->converged() ? "true" : "false");
            services_.observer->onSolverMessage({
                FDM::SolverMessageKind::StateClosure,
                "Fixed-time iteration",detail.str(),
                state_->step,state_->time,state_->dt});
        }
    });
    operations.bind(System::OpIds::PressureIterationEnd,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->endIteration();
    });
    operations.bind(System::OpIds::PressureStepCommit,"flow.pressure-operators",[owner=shared_from_this(),this] {
        pressureOperators_->commitStep();
        pressureOperators_->traceState("accepted",state_->step,state_->time+state_->dt);
    });
}

void FlowOperations::bindExplicitOps(
        Run::OpRegistry& operations,
        const std::vector<Field*>& fields,
        const double& maximumTimeStep,
        Time::Explicit::Workspace& explicitWorkspace) {
    operations.bind(System::OpIds::FlowStepPrepare,"flow.conservative",[owner=shared_from_this(),this,&fields,&maximumTimeStep] {
        resetStep();
        if (!std::isfinite(maximumTimeStep) || maximumTimeStep <= 0.0) {
            throw std::runtime_error(
                "SolverAlgorithm::FlowOperations requires a finite positive time-step limit.");
        }
        HighOrderTrace::beginPhysicalStep(state_->step);
        if (HighOrderTrace::activeFor(state_->step)) {
            HighOrderTrace::conservative("initial Q", fields);
            HighOrderTrace::gamma(*state_->stateModel);
        }
        prepareBoundaryState(fields, state_->time, 0.0);
        if (HighOrderTrace::activeFor(state_->step)) {
            HighOrderTrace::conservative("boundary-prepared Q", fields, true);
        }
    });
    operations.bind(System::OpIds::FlowDtCompute,"flow.conservative",[owner=shared_from_this(),this,&fields,&maximumTimeStep] {
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
                    *state_->stateModel->perfectGasGamma());
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
    operations.bind(System::OpIds::FlowStepBegin,"flow.conservative",[owner=shared_from_this(),this,&fields,&explicitWorkspace] {
        beginExplicitStep();
    });
    operations.bind(System::OpIds::ExplicitStageExecute,"flow.conservative",
        [owner=shared_from_this(),this,&fields,&explicitWorkspace](
                const Run::ExecutionContext& context) {
            executeExplicitStage(context.stageIndex);
        });
}

double FlowOperations::proposeTimeStep() const {
    double dt=std::numeric_limits<double>::max();
    for (auto* field:state_->patches) dt=std::min(dt,Numerics::RusanovEOS::deltaT(*field,*state_->stateModel,numerics_.dt.cfl));
    return services_.executionRuntime?services_.executionRuntime->globalMinimum(dt):dt;
}
void FlowOperations::beginExplicitStep() {
        if (services_.equationSystem) {
            services_.equationSystem->beginStep(state_->patches, state_->dt);
        }
        for (Field* field:state_->patches) {
            auto* registry=services_.equationSystem ? services_.equationSystem->variables(*field) : &state_->transported;
            if (!registry) continue;
            for (const auto& variable:registry->variables()) {
                const auto policy=variable.descriptor.updatePolicy;
                if (policy!=State::UpdatePolicy::Explicit && policy!=State::UpdatePolicy::BoundedExplicit
                    && policy!=State::UpdatePolicy::HamiltonJacobi) continue;
                if (std::none_of(solve_.compiledProgram.steps.begin(),solve_.compiledProgram.steps.end(),[&](const auto& call) {
                    return call.temporalResidual && call.source.target.symbol==variable.descriptor.name;
                })) throw std::runtime_error("Registered scalar has no compiled temporal HOW occurrence: "+variable.descriptor.name);
            }
        }
        ensureWorkspaces(state_->patches);
        Time::Explicit::begin(
            explicitWorkspace_,state_->patches,*state_,numerics_.time.recipe,
            services_.equationSystem);
        if(solve_.compiledProgram.temporalParticipants.empty())bindTemporalStateViews();
}
void FlowOperations::executeExplicitStage(int stage) {
            Time::Explicit::executeStage(
                explicitWorkspace_,stage,state_->patches,workspaces_,*state_,
                services_.equationSystem,
                [owner=shared_from_this(),this](const std::vector<Field*>& patches,
                       std::vector<PatchWorkspace>& workspaces,
                       double stageTime) {
                    ConservativeRHS::assembleAllPatches(
                        patches,workspaces,rhsNumerics(),
                        *equations_,*state_,services_,
                        boundaryApplicator_,stageTime);
                },
                [owner=shared_from_this(),this](const std::vector<Field*>& patches) {
                    ConservativeRHS::publishIntegratedState(
                        patches,*state_,services_);
                },
                [owner=shared_from_this(),this](const std::vector<Field*>& patches,const char* stage) {
                    validateDensityStateClosure(patches,stage);
                });
}
void FlowOperations::bindFieldSources() {
    if(pressureOperators_)pressureOperators_->bindFieldSources();
    if(solve_.compiledProgram.temporalParticipants.empty())return;
    if (std::none_of(numerics_.operators.begin(),numerics_.operators.end(),[](const auto& term){return bool(term.stageSource.evaluate) || bool(term.fieldSource.bind);})) return;
    stageNumerics_=std::make_unique<System::CompiledNumericalSystem>(numerics_);
    for(auto& term:stageNumerics_->operators)if(term.fieldSource.bind) {
        std::vector<System::FieldSource::Reader> reads;
        for(const auto& use:term.fieldSource.inputs)reads.push_back(realizedState_.reader(use,numerics_.time.recipe.stageCount()));
        const int offset=term.output=="rhoU"?1:4,count=term.output=="rhoU"?3:1;
        term.conservativeSource=[value=term.fieldSource.bind(std::move(reads)),offset,count](Field& field,Residual& residual) {
            Math::forFluidInterior(field,[&](int i,int j,int k) {
                const int cell=field.getIdx(i,j,k);
                for(int c=0;c<count;++c)residual.source(i,j,k,offset+c)+=value(cell,c);
            });
        };
    }
    for (auto& term:stageNumerics_->operators) if(term.stageSource.evaluate) {
        std::vector<System::StageConservativeSource::Reader> reads;
        for(const auto& symbol:term.stageSource.reads) reads.push_back(realizedState_.stageReader(symbol,numerics_.time.recipe.stageCount()));
        term.conservativeSource=[kernel=term.stageSource.evaluate,reads=std::move(reads)](Field& field,Residual& residual){kernel(field,residual,reads);};
    }
}
Time::TemporalCallbacks FlowOperations::temporalParticipant() {
    if (state_->patches.size()!=1 || (services_.executionRuntime && services_.executionRuntime->distributed())
        || services_.immersed.boundary || services_.immersed.constraint || services_.immersed.system || services_.equationSystem
        || (services_.transportModel && !std::any_of(solve_.compiledProgram.steps.begin(),solve_.compiledProgram.steps.end(),[](const auto& call){return call.stageLocal && std::find(call.coefficientStates.begin(),call.coefficientStates.end(),"mu_t")!=call.coefficientStates.end();})))
        throw std::runtime_error("Unsupported shared flow stage group: serial single patch without auxiliary/immersed stage treatment required.");
    bindTemporalStateViews();
    auto owner=shared_from_this();Time::TemporalCallbacks result;result.identity="flow.conservative";
    result.stepSize=[owner]{return owner->proposeTimeStep();};
    result.snapshot=[owner]{owner->beginExplicitStep();};
    result.prepare=[owner](const auto& c) {
        const double time=owner->state_->time+owner->numerics_.time.recipe.stage(c.stageIndex).abscissa*owner->state_->dt;
        owner->prepareBoundaryState(owner->state_->patches,time,owner->state_->dt);
    };
    result.rhs=[owner](const auto& c) {
        Time::Explicit::evaluateRHS(owner->explicitWorkspace_,c.stageIndex,owner->state_->patches,owner->workspaces_,*owner->state_,
            [owner](const auto& patches,auto& workspace,double time) {
                ConservativeRHS::assembleAllPatches(patches,workspace,owner->rhsNumerics(),*owner->equations_,
                    *owner->state_,owner->services_,owner->boundaryApplicator_,time);
            });
        const auto& recipe=System::NumericalCompiler::requireUniqueRecipe(owner->rhsNumerics(),FDM::TermRole::Convection);
        auto& workspace=owner->workspaces_.front();
        workspace.massFlux->publish(owner->state_->singlePatch(),workspace.convectiveFlux,
            {owner->state_->step,c.stageIndex,owner->state_->time+owner->numerics_.time.recipe.stage(c.stageIndex).abscissa*owner->state_->dt},FDM::toString(recipe.id()));
    };
    result.advance=[owner](const auto& c){owner->executeExplicitStage(c.stageIndex);
        for(auto& workspace:owner->workspaces_)workspace.massFlux->invalidate();};
    result.validatePublish=[owner] {
        if (owner->explicitWorkspace_.active) throw std::runtime_error("Flow publication before complete stages.");
        owner->validateDensityStateClosure(owner->state_->patches,"common stage final");
    };
    result.publish=[owner]{owner->prepareBoundaryState(owner->state_->patches,owner->state_->time+owner->state_->dt,owner->state_->dt);};
    return result;
}

StageMassFluxView FlowOperations::stageMassFlux(StageFluxIdentity stamp) const {
    if(workspaces_.size()!=1)throw std::runtime_error("Shared canonical mass flux requires one participating patch.");
    const auto& recipe=System::NumericalCompiler::requireUniqueRecipe(rhsNumerics(),FDM::TermRole::Convection);
    return {workspaces_.front().massFlux,state_->singlePatch(),stamp,FDM::toString(recipe.id())};
}

void FlowOperations::resetStep() {
    if ((requiresOperation(System::OpIds::ExplicitStageExecute) || requiresOperation("flow.stage.rhs")))
        ConservativeRHS::validateTimestepState(*state_,numerics_);
    pressureSummary_={};
    immersedCorrection_={};
}

void FlowOperations::bindOperations(Run::OpRegistry& operations) {
    if ((requiresOperation(System::OpIds::ExplicitStageExecute) || requiresOperation("flow.stage.rhs")))
        bindExplicitOps(operations,state_->patches,maximumTimeStep_,explicitWorkspace_);
    if (requiresOperation(System::OpIds::PressurePrepare)) {
        if (pressureOperators_)
            bindConstantPressureOps(operations,maximumTimeStep_);
        else if (conservativePressureCorrector_)
            bindPressureOps(operations,state_->singlePatch(),maximumTimeStep_,pressureSummary_);
        else throw std::runtime_error("Requested pressure operations have no bound pressure workspace.");
    }
    auto finalizeFlow=[owner=shared_from_this(),this] {
        ConservativeRHS::publishCorrectedConservativeState(immersedCorrection_.performed,services_);
        if (immersedCorrection_.performed) {
            validateDensityStateClosure(state_->patches,"flow correction");
            if (services_.observer)
                services_.observer->onSolverMessage({FDM::SolverMessageKind::StateClosure,
                    "Flow correction",immersedCorrection_.detail,state_->step,state_->time,state_->dt});
        }
    };
    if (requiresOperation(System::OpIds::FlowFinalizeBegin))
        operations.bind(System::OpIds::FlowFinalizeBegin,"flow.conservative",finalizeFlow);
    if (requiresOperation(System::OpIds::FlowStepCommit))
        operations.bind(System::OpIds::FlowStepCommit,"flow.conservative",[owner=shared_from_this(),this,finalizeFlow] {
            if (!requiresOperation(System::OpIds::FlowFinalizeBegin)) finalizeFlow();
            if (services_.equationSystem) services_.equationSystem->commitStep(state_->patches,state_->dt);
            prepareBoundaryState(state_->patches,state_->time+state_->dt,state_->dt);
        });
    auto assigned=requiredOperations_;
    // StageGroup owns temporal callbacks. Keep only requested, non-temporal
    // operations here; a provider label never acquires an entire operation set.
    for (const auto& p:solve_.compiledProgram.temporalParticipants)
        assigned.erase(std::remove_if(assigned.begin(),assigned.end(),[&](const auto& op) {
            return op==p.snapshot || op==p.prepareStage || op==p.rhs || op==p.advance || op==p.publish;
        }),assigned.end());
    operations.retain(assigned);
}
} // namespace SF::SolverAlgorithm
