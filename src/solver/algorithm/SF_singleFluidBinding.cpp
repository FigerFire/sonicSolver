#include "SF_singleFluidBinding.h"
#include "SF_flowOperations.h"
#include "scalar/SF_scalarBinding.h"
#include "scalar/SF_thermalOperations.h"
#include "solver/system/SF_scalarMethod.h"
#include "transport/SF_transportBinding.h"
#include "interface/SF_levelSetBinding.h"
#include "immersed/SF_immersedBinding.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
namespace SF::SolverAlgorithm {
void bindSingleFluidOperations(Run::OpRegistry& operations,const FDM::BoundaryConfig& boundaries,const FDM::PressureCorrectionConfig& pressure,double pressureGamma,
    const System::ExecutableEquationSystem& equations,const System::CompiledNumericalSystem& numerics,
    const System::CompiledSolvePlan& plan,const System::RuntimeRequirements& requirements,
    State::StateBundle& state,const FDM::SolverServices& services,
    System::StateRealization& realized,const double& maximumTimeStep) {
    std::shared_ptr<FlowOperations> flow;
    if (System::requiresOperation(requirements,System::OpIds::FlowStepPrepare) || System::requiresOperation(requirements,System::OpIds::PressurePrepare)) {
        flow=std::make_shared<FlowOperations>(boundaries,pressure,pressureGamma,equations,numerics,plan,requirements,state,services,realized,maximumTimeStep);
        flow->prepare();flow->bindOperations(operations);
    }
    if (!plan.compiledProgram.temporalParticipants.empty()) {
        std::vector<Time::TemporalCallbacks> callbacks;
        // Bind all stage views before constructing cross-provider readers.
        if (flow) callbacks.push_back(flow->temporalParticipant());
        auto scalars=scalarParticipants(numerics,plan,state,services,realized,maximumTimeStep,
            flow?std::function<StageMassFluxView(StageFluxIdentity)>([flow](auto stamp){return flow->stageMassFlux(stamp);}):std::function<StageMassFluxView(StageFluxIdentity)>{});
        auto thermal=thermalParticipants(numerics,plan,state,services,realized);
        callbacks.insert(callbacks.end(),thermal.begin(),thermal.end());
        callbacks.insert(callbacks.end(),scalars.begin(),scalars.end());
        std::vector<Time::TemporalCallbacks> ordered;
        for (const auto& p:plan.compiledProgram.temporalParticipants) {
            const auto found=std::find_if(callbacks.begin(),callbacks.end(),[&](const auto& c){return c.identity==p.identity;});
            if (found==callbacks.end()) throw std::runtime_error("Missing temporal numerical owner: "+p.identity);
            ordered.push_back(*found);
        }
        auto readers=bindAlgebraicOperations(operations,numerics,plan,state,services,realized);
        Time::bindStageGroup(operations,plan,numerics,state,realized,maximumTimeStep,std::move(ordered),std::move(readers));
        if (flow) flow->bindFieldSources();
        return;
    }
    auto readers=bindAlgebraicOperations(operations,numerics,plan,state,services,realized);
    if(!readers.empty()){realized.bindDerivedViews();for(auto& bind:readers)bind();}
    if(!flow)return;
    flow->bindFieldSources();
    bindTransportOperations(operations,requirements,state,services);
    bindLevelSetOperations(operations,plan,state,services,realized);
    bindImmersedOperations(operations,requirements,plan,pressure,pressureGamma,state,services,
        flow->correctionResult(),[flow](const auto& fields,double time,double dt) {
            flow->prepareBoundaryState(fields,time,dt);
        });
}
}
