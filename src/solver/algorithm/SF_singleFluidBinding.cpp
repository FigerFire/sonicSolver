#include "SF_singleFluidBinding.h"
#include "SF_flowOperations.h"
#include "scalar/SF_scalarBinding.h"
#include "scalar/SF_cellBinding.h"
#include "SF_faceBinding.h"
#include "solver/system/SF_affine.h"
#include "solver/system/SF_block.h"
#include "solver/system/SF_scalarMethod.h"
#include "solver/system/SF_sstMethod.h"
#include "transport/SF_transportBinding.h"
#include "interface/SF_levelSetBinding.h"
#include "immersed/SF_immersedBinding.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <map>
#include <typeindex>
namespace SF::SolverAlgorithm {
namespace {
void bindInterval(Run::OpRegistry& ops,const System::RuntimeRequirements& requirements,
        State::StateBundle& state,System::StateRealization& views,const double& limit,
        std::vector<IntervalBinding> calls) {
    const auto action=[&](std::string_view id)->const System::ResolvedOperationBinding& {
        const auto found=std::find_if(requirements.operationBindings.begin(),requirements.operationBindings.end(),
            [&](const auto& binding){return binding.operation==id && binding.status==System::BindingStatus::Resolved;});
        if(found==requirements.operationBindings.end())
            throw std::runtime_error("Missing interval operation binding: "+std::string(id));
        return *found;
    };
    const auto& prepare=action(System::OpIds::CellPrepare);
    const auto& accept=action(System::OpIds::CellAccept);
    // Prepared calls establish the interval's geometry/velocity state before
    // proposal-only clients. Keep that order for the private base snapshot,
    // OldTime publication and acceptance, and preserve HOW order in each set.
    std::stable_partition(calls.begin(),calls.end(),[](const auto& call){return bool(call.prepare);});
    const bool publishesOldTime=std::any_of(calls.begin(),calls.end(),[](const auto& call){return call.publishesOldTime;});
    ops.bind(prepare.operation,prepare.provider,[calls,&state,&limit,&views,publishesOldTime] {
        double dt=limit;
        for(const auto& call:calls) {
            if(call.prepare)call.prepare();
            dt=std::min(dt,call.stepSize());
        }
        if(!std::isfinite(dt) || dt<=0)throw std::runtime_error("Invalid common equation interval");
        state.dt=dt;
        if(publishesOldTime)views.beginOldTime();
        for(const auto& call:calls)if(call.snapshot)call.snapshot();
        if(publishesOldTime)views.publishOldTime();
    });
    ops.bind(accept.operation,accept.provider,[calls=std::move(calls)] {
        for(const auto& call:calls)if(call.accept)call.accept();
    });
}
}
void bindSingleFluidOperations(Run::OpRegistry& operations,const FDM::BoundaryConfig& boundaries,const FDM::PressureCorrectionConfig& pressure,double pressureGamma,
    const System::ExecutableEquationSystem& equations,const System::CompiledNumericalSystem& numerics,
    const System::CompiledSolvePlan& plan,const System::RuntimeRequirements& requirements,
    State::StateBundle& state,const FDM::SolverServices& services,
    System::StateRealization& realized,const double& maximumTimeStep) {
    const auto& program=plan.compiledProgram;
    const bool stages=!program.temporalParticipants.empty();
    std::map<std::string,std::vector<const System::ResidualDefinition*>> residuals;
    for(const auto& definition:program.residuals)residuals[definition.occurrence].push_back(&definition);

    auto publications=std::make_shared<std::set<std::string>>();
    std::shared_ptr<CellEquationExecution> response;
    const System::CompiledEquationCall* faceCall=nullptr;
    const bool fixedResponse=System::requiresOperation(requirements,System::OpIds::PressureStepBegin);
    using Bind=std::function<BoundEquation(const EquationBindingInput&)>;
    const std::map<std::type_index,Bind> methods{
        {typeid(System::CompiledAffineEquation),[&](const auto& input){return bindAffine(input,publications);}},
        {typeid(System::CompiledSstClosure),bindSst},
        {typeid(System::CompiledScalarEquation),bindScalar},
        {typeid(System::CompiledLinearBlock),bindBlock},
        {typeid(System::CompiledCellEquation),[&](const auto& input) {
            const bool responds=numerics.facePrediction && numerics.facePrediction->response==input.call.source.occurrence;
            std::shared_ptr<CellEquationExecution> instance;
            auto result=bindCell(input,responds,responds && fixedResponse,instance);
            if(responds)response=std::move(instance);
            return result;
        }}
    };
    std::set<std::string> scalarPorts;
    std::vector<BoundEquation> bindings;
    std::vector<IntervalBinding> interval;
    std::map<std::string,Time::TemporalCallbacks> participants;
    std::size_t singleScalarUsers=0;
    // Attach every output before any consumer resolves a reader.
    for(const auto& call:program.steps) {
        for(const auto& resource:call.physicalResources)
            if(resource.fields && (resource.owner=="native.scalar" || resource.owner=="model.scalar"))
                scalarPorts.insert(resource.id);
        if(call.faceFlux && !call.faceFlux->corrected) {
            if(faceCall)throw std::runtime_error("Multiple independent face layouts/owners are not implemented: "+call.source.occurrence);
            faceCall=&call;
        }
        const auto type=std::type_index(call.providerContract.type());
        const auto method=methods.find(type);
        if(method==methods.end())continue; // Fused and specialized operation owners use their existing binding.
        auto bound=method->second({call,numerics,state,services,realized,maximumTimeStep,residuals[call.source.occurrence]});
        if(bound.singleScalarPort)++singleScalarUsers;
        if(!bound.stage.identity.empty()) {
            auto identity=bound.stage.identity;
            if(!participants.emplace(std::move(identity),std::move(bound.stage)).second)
                throw std::runtime_error("Repeated temporal equation binding: "+call.source.occurrence);
        }
        if(bound.interval.stepSize)interval.push_back(std::move(bound.interval));
        bindings.push_back(std::move(bound));
    }
    if(services.scalarEquation && (!services.scalarInstances.empty() || singleScalarUsers!=1))
        throw std::runtime_error("Single scalar data port cannot bind multiple occurrences or coexist with instance table.");
    for(const auto& item:services.scalarInstances)
        if(item.data.expired() || !scalarPorts.erase(item.occurrence))
            throw std::runtime_error("Invalid, duplicate or expired scalar instance data key: "+item.occurrence);

    std::shared_ptr<FaceFluxExecution> faces;
    if(faceCall) {
        if(!response || !services.executionRuntime)
            throw std::runtime_error("Missing face response occurrence/runtime: "+faceCall->faceFlux->response);
        faces=std::make_shared<FaceFluxExecution>(*faceCall,response,state,*services.executionRuntime,realized);
    }
    std::shared_ptr<FlowOperations> flow;
    if (System::requiresOperation(requirements,System::OpIds::FlowStepPrepare) || System::requiresOperation(requirements,System::OpIds::PressurePrepare)) {
        flow=std::make_shared<FlowOperations>(boundaries,pressure,pressureGamma,equations,numerics,plan,requirements,state,services,realized,maximumTimeStep);
        flow->prepare(faces);
        if(stages) {
            auto callbacks=flow->temporalParticipant();
            auto identity=callbacks.identity;
            if(!participants.emplace(std::move(identity),std::move(callbacks)).second)
                throw std::runtime_error("Repeated fused temporal owner binding");
        }
    }

    FluxReaders fluxes;
    if(stages && flow)fluxes.mass=[flow](auto stamp){return flow->stageMassFlux(stamp);};
    if(faces)fluxes.volume=[faces](std::string_view symbol,int axis,int lower){return faces->published(symbol,axis,lower);};
    if(stages || !bindings.empty())realized.bindDerivedViews();
    for(auto& bound:bindings)if(bound.connect)bound.connect(fluxes);
    if(flow)flow->bindFieldSources();

    for(auto& bound:bindings)if(bound.install)bound.install(operations);
    bindVariableCommits(operations,program.root,state,realized);
    if(faces)operations.bind(faceCall->faceFlux->publication,faceCall->backendProvider,[faces]{faces->predict();});
    if(flow)flow->bindOperations(operations);
    if(stages) {
        std::vector<Time::TemporalCallbacks> ordered;
        for(const auto& participant:program.temporalParticipants) {
            const auto found=participants.find(participant.identity);
            if(found==participants.end())throw std::runtime_error("Missing temporal numerical owner: "+participant.identity);
            ordered.push_back(std::move(found->second));participants.erase(found);
        }
        if(!participants.empty())throw std::runtime_error("Unconsumed temporal numerical owner: "+participants.begin()->first);
        Time::bindStageGroup(operations,plan,numerics,state,realized,maximumTimeStep,std::move(ordered));
    }
    if(System::requiresOperation(requirements,System::OpIds::CellPrepare))
        bindInterval(operations,requirements,state,realized,maximumTimeStep,std::move(interval));

    // Increment samplers borrow the declared output. Equation-defect/linear
    // samplers remain numerical-method responsibilities, never guessed here.
    for(const auto& residual:program.residuals)if(!stages && !operations.containsResidual(residual.id,residual.provider)
            && residual.kind==System::ResidualKind::StateChange) {
        System::StateUse use{residual.target,System::StateVersion::Current};use.kind=residual.view;
        auto read=realized.reader(use,numerics.time.recipe.stageCount());
        auto* mesh=&state.singlePatch();const int components=equations.state.at(residual.target).components;
        operations.bindResidual(residual.id,residual.provider,[read=std::move(read),mesh,components](const auto&) {
            std::vector<double> result;const int g=mesh->NG();
            for(int k=g;k<g+mesh->NZ();++k)for(int j=g;j<g+mesh->NY();++j)for(int i=g;i<g+mesh->NX();++i)
                for(int c=0;c<components;++c)result.push_back(read(mesh->getIdx(i,j,k),c));
            return result;
        });
    }
    if(!flow || stages)return;
    bindTransportOperations(operations,requirements,state,services);
    bindLevelSetOperations(operations,plan,state,services,realized);
    bindImmersedOperations(operations,requirements,plan,pressure,pressureGamma,state,services,
        flow->correctionResult(),[flow](const auto& fields,double time,double dt) {
            flow->prepareBoundaryState(fields,time,dt);
        });
}
}
