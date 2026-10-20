#include "SF_immersedBinding.h"
#include "solver/algorithm/immersed/SF_constraintOps.h"
#include "solver/algorithm/immersed/SF_immersedStrategy.h"
#include "solver/algorithm/pressureBased/SF_kkt.h"
#include "solver/system/SF_immersedMethods.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
namespace SF::SolverAlgorithm {
void bindImmersedOperations(Run::OpRegistry& operations,
    const System::RuntimeRequirements& requirements, const System::CompiledSolvePlan& plan,
    const FDM::SolverConfig& config, State::StateBundle& state,
    const FDM::SolverServices& services, FDM::ImmersedConstraintResult& correction,
    std::function<void(const std::vector<Field*>&,double,double)> prepareBoundaryState) {
    const auto required=[&](std::string_view id) {
        return std::any_of(requirements.operationBindings.begin(),requirements.operationBindings.end(),
            [&](const auto& binding) { return binding.status==System::BindingStatus::Resolved
                && binding.provider=="ibm.constraint" && binding.operation==id; });
    };
    for (const auto& binding:requirements.operationBindings) {
        if (binding.status!=System::BindingStatus::Resolved || binding.provider!="ibm.constraint") continue;
        if (binding.operation!=System::OpIds::IbmConstraintProject && binding.operation!=System::OpIds::IbmKktSolve)
            throw std::runtime_error("Unimplemented immersed operation '"+binding.operation+"' for provider ibm.constraint.");
        if (!services.immersed.constraint || !services.immersed.system)
            throw std::runtime_error("Missing immersed constraint/system ports for provider ibm.constraint, operation '"+binding.operation+"'.");
    }
    std::shared_ptr<PressureBased::MonolithicKKT> kktProvider;
    if (required(System::OpIds::IbmKktSolve))
        kktProvider=std::make_shared<PressureBased::MonolithicKKT>(config.pressure,config.numerics.idealGasGamma);
    if (required(System::OpIds::IbmConstraintProject)
        && services.immersed.constraint && services.immersed.system) {
        operations.bind(System::OpIds::IbmConstraintProject,"ibm.constraint",[&state,&services,&correction,prepareBoundaryState,kktProvider] {
            const auto immersed = ImmersedAlgorithm::projectConstraint(
                state.patches,state.time+state.dt,state.dt,
                services.executionRuntime,services.immersed);
            correction.performed = immersed.performed;
            correction.detail = immersed.detail;
        });
    }
    if (required(System::OpIds::IbmKktSolve) && kktProvider
        && services.immersed.constraint && services.immersed.system) {
        operations.bind(System::OpIds::IbmKktSolve,"ibm.constraint",[&state,&services,&correction,prepareBoundaryState,kktProvider] {
            if (state.patches.size() != 1 || !state.patches.front()) {
                throw std::runtime_error(
                    "ibm.kkt.solve currently requires one local Field.");
            }
            if (services.executionRuntime
                && services.executionRuntime->distributed()) {
                throw std::runtime_error(
                    "ibm.kkt.solve has no distributed ConstraintGlobalDof "
                    "provider for the pressure block.");
            }
            ImmersedAlgorithm::validateMonolithicProvider(
                *services.immersed.system);
            services.immersed.constraint->setExecutionRuntime(
                services.executionRuntime);
            prepareBoundaryState(
                state.patches,state.time+state.dt,state.dt);
            if (services.equationSystem) {
                services.equationSystem->preparePressureCorrection(state.patches);
            }
            const auto kkt = kktProvider->correct(
                *state.patches.front(),state.time+state.dt,state.dt,
                *services.immersed.constraint);
            std::ostringstream detail;
            detail << "KKT iterations=" << kkt.iterations
                   << ", residual=" << kkt.relativeResidual
                   << ", maxDiv(before/after)="
                   << kkt.maxDivergenceBefore << "/"
                   << kkt.maxDivergenceAfter
                   << ", HYPRE(rebuilds/solves)="
                   << kkt.structureRebuilds << "/" << kkt.linearSolves
                   << ", " << kkt.immersed.detail;
            correction.performed = true;
            correction.detail = detail.str();
        });
    }

    if (System::requiresProvider(requirements,"ibm.boundary") && !services.immersed.boundary)
        throw std::runtime_error("Compiled Ghost boundary has no bound closure port.");
    System::validateImmersedBindings(plan,services.immersed.system,
        System::requiresProvider(requirements,"ibm.boundary"));
}
}
