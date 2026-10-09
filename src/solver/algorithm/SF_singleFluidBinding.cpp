#include "SF_singleFluidBinding.h"
#include "SF_flowOperations.h"
#include "scalar/SF_scalarBinding.h"
#include "solver/system/SF_scalarMethod.h"
#include "transport/SF_transportBinding.h"
#include "interface/SF_levelSetBinding.h"
#include "immersed/SF_immersedBinding.h"
#include <algorithm>
#include <stdexcept>
namespace SF::SolverAlgorithm {
void bindSingleFluidOperations(Run::OpRegistry& operations,const FDM::SolverConfig& config,
    const System::ExecutableEquationSystem& equations,const System::CompiledNumericalSystem& numerics,
    const System::CompiledSolvePlan& plan,const System::RuntimeRequirements& requirements,
    State::StateBundle& state,const FDM::SolverServices& services,
    System::StateRealization& realized,const double& maximumTimeStep) {
    std::shared_ptr<FlowOperations> flow;
    if (System::requiresProvider(requirements,"flow.conservative") || System::requiresProvider(requirements,"flow.pressure-operators")) {
        flow=std::make_shared<FlowOperations>(config,equations,numerics,plan,requirements,state,services,realized,maximumTimeStep);
        flow->prepare();flow->bindOperations(operations);
    }
    if (!plan.compiledProgram.temporalParticipants.empty()) {
        std::vector<Time::TemporalCallbacks> callbacks;
        // Bind all stage views before constructing cross-provider readers.
        if (flow) callbacks.push_back(flow->temporalParticipant());
        auto scalars=scalarParticipants(numerics,plan,state,services,realized,maximumTimeStep);
        if (flow) flow->bindStageSources();
        callbacks.insert(callbacks.end(),scalars.begin(),scalars.end());
        std::vector<Time::TemporalCallbacks> ordered;
        for (const auto& p:plan.compiledProgram.temporalParticipants) {
            const auto found=std::find_if(callbacks.begin(),callbacks.end(),[&](const auto& c){return c.identity==p.identity;});
            if (found==callbacks.end()) throw std::runtime_error("Missing temporal numerical owner: "+p.identity);
            ordered.push_back(*found);
        }
        Time::bindStageGroup(operations,plan,numerics,state,realized,maximumTimeStep,std::move(ordered));
        return;
    }
    bindTransportOperations(operations,requirements,state,services);
    bindLevelSetOperations(operations,plan,state,services,realized);
    bindImmersedOperations(operations,requirements,plan,config,state,services,
        flow->correctionResult(),[flow](const auto& fields,double time,double dt) {
            flow->prepareBoundaryState(fields,time,dt);
        });
}
}
