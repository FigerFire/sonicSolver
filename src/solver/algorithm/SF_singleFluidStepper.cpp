#include "SF_singleFluidStepper.h"
#include "SF_singleFluidBinding.h"
#include "core/system/SF_operationIds.h"
#include "core/config/SF_numericsPolicy.h"
#include <sstream>
#include <stdexcept>
#include <utility>
namespace SF::SolverAlgorithm {
SingleFluidStepper::SingleFluidStepper(FDM::SolverConfig config,
    const System::ExecutableEquationSystem& equations,const System::CompiledNumericalSystem& numerics,
    const System::CompiledSolvePlan& plan,const System::RuntimeRequirements& runtime)
    :config_(std::move(config)),executable_(equations),numerics_(numerics),solve_(plan),runtime_(runtime) {
    if (System::requiresProvider(runtime_,"flow.conservative")
        || System::requiresProvider(runtime_,"flow.pressure-operators"))
        FDM::validateNumericsConfig(config_.numerics);
    if (numerics_.time.recipe.id()!=config_.numerics.timeRecipe.id()
        || numerics_.dt.cfl!=config_.numerics.cfl
        || numerics_.dt.maxDeltaT!=config_.numerics.maxDeltaT)
        throw std::runtime_error("SingleFluidStepper received a compiled numerical system whose time/step policy disagrees with its solver configuration.");
}
void SingleFluidStepper::bindServices(FDM::SolverServices services) {
    if (servicesBound_) throw std::runtime_error("SingleFluidStepper services are already bound for this lifecycle.");
    services_=services;servicesBound_=true;
}
void SingleFluidStepper::bindSolvePlan(const System::CompiledSolvePlan& plan) {
    if (solvePlanBound_ || &plan!=&solve_) throw std::runtime_error("SingleFluidStepper requires its unbound resolved plan.");
    solvePlanBound_=true;
}
void SingleFluidStepper::prepare(FDM::SolverState& state) {
    if (!servicesBound_ || !solvePlanBound_ || !state.bundle)
        throw std::runtime_error("SingleFluidStepper requires services, solve plan and StateBundle before preparation.");
    if (prepared_) throw std::runtime_error("SingleFluidStepper runtime operations are already bound.");
    state.bundle->validatePatches();state_=state.bundle;
    boundPatchIdentities_=state_->patches;
    maximumTimeStep_=state.maximumTimeStep;
    if (services_.executionRuntime) services_.executionRuntime->attachState(*state_);
    realizedState_=System::realizeState(executable_.state,runtime_,*state_,solve_.compiledProgram.stateViews);
    bindSingleFluidOperations(operations_,config_,executable_,numerics_,solve_,runtime_,
                             *state_,services_,realizedState_,maximumTimeStep_);
    for (const auto& binding:runtime_.operationBindings) {
        if (binding.operation!=System::OpIds::TimeCommit) continue;
        if (binding.status!=System::BindingStatus::Resolved)
            throw std::runtime_error("Physical clock commit has no frozen provider binding.");
        operations_.bind(binding.operation,binding.provider,[this] {
            state_->time+=state_->dt;++state_->step;emitTimeStep();
        });
    }
    Run::PlanExecutor::validateBindings(solve_,operations_);
    prepared_=true;
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

double SingleFluidStepper::physicalTime() const {
    if (!state_) throw std::runtime_error("SingleFluidStepper has no bound StateBundle.");
    return state_->time;
}
double SingleFluidStepper::timeStep() const {
    if (!state_) throw std::runtime_error("SingleFluidStepper has no bound StateBundle.");
    return state_->dt;
}
FDM::StepResult SingleFluidStepper::advance(FDM::SolverState& state) {
    FDM::StepResult result;
    if (!prepared_) { result.message="SingleFluidStepper requires prepare before advance.";return result; }
    if (state.bundle!=state_) throw std::runtime_error("SingleFluidStepper cannot switch its authoritative StateBundle after binding.");
    state_->validatePatches();
    if (state_->patches!=boundPatchIdentities_)
        throw std::runtime_error("Participating patch identities changed after immutable runtime binding.");
    maximumTimeStep_=state.maximumTimeStep;
    const Run::PlanTraceContext trace{state_->step,state_->time,&state_->dt};
    Run::PlanExecutor::execute(solve_,operations_,&trace);
    result.accepted=true;result.dt=state_->dt;result.time=state_->time;result.step=state_->step;
    return result;
}
}
