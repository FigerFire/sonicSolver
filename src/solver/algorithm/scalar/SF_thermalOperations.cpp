#include "SF_thermalOperations.h"
#include "solver/system/SF_thermalMethod.h"
#include "methods/numerics/scalar/SF_scalarTransport.h"
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
namespace SF::SolverAlgorithm {
namespace {
class ThermalOperations {
    State::StateBundle& state_;
    Field& field_;
    const FDM::ThermalEquationData& data_;
    FDM::Scalar::TransportConfig config_;
    ScalarField stage_;
    std::vector<double> old_,rhs_;
    State::DistributedFieldView oldView_,stageView_,rhsView_;
    bool ready_=false;
public:
    ThermalOperations(const System::CompiledEquationCall& call,State::StateBundle& state,
            const FDM::ThermalEquationData& data,System::StateRealization& realized)
        :state_(state),field_(state.singlePatch()),data_(data),stage_(*data.temperature) {
        const auto& p=std::any_cast<const std::map<std::string,double>&>(call.providerContract);
        config_.name=data.temperature->name();config_.convectionEnabled=p.at("convection");
        config_.diffusionEnabled=p.at("diffusion");config_.diffusivity=p.at("D");
        config_.ilwBoundaryEnabled=p.at("ilw");config_.ilwAccuracyOrder=(int)p.at("ilwOrder");
        if(config_.ilwBoundaryEnabled && Math::activeDimensionCount()<2)
            throw std::runtime_error("Unsupported thermal ILW: requires an active tangential direction.");
        if(config_.convectionEnabled && (field_.NVar()!=5 || !field_.hasStateModel()
            || !field_.stateModel()->perfectGasGamma()))
            throw std::runtime_error("Unsupported thermal velocity port: bound five-variable PerfectGas state required by the current scalar convection kernel.");
        config_.bounds={p.at("bounded")!=0,p.at("lower"),p.at("upper"),(FDM::Scalar::BoundMode)(int)p.at("boundMode")};
        auto* physical=realized.at(call.target.symbol).fields.at(0);
        if (physical->geometry!=&field_ || physical->name!=config_.name || !data.temperature->isCompatibleWith(field_)
            || physical->read(0,0)!=data.temperature->values().at(0))
            throw std::runtime_error("Thermal STATE and physical port disagree.");
        // Pointer equality of the registered scalar, not a value comparison, proves the alias.
        const auto& vars=state.transported.variables();
        if(std::count_if(vars.begin(),vars.end(),[&](const auto& v){return v.descriptor.name==physical->name && v.value==data.temperature;})!=1)
            throw std::runtime_error("Thermal port must alias the authoritative registered physical scalar.");
        old_.resize(field_.TotalSize());rhs_.resize(field_.TotalSize());
        oldView_=State::workspaceView(call.oldTimeWorkspace,0,field_,old_,1,0,State::HaloSyncStage::None);
        rhsView_=State::workspaceView(call.residualWorkspace,0,field_,rhs_,1,0,State::HaloSyncStage::None);
        stageView_=State::scalarView(call.stageWorkspace,0,field_,stage_);
        realized.bindView(call.target.symbol,System::StateViewKind::OldTime,oldView_);
        realized.bindView(call.target.symbol,System::StateViewKind::Stage,stageView_,0);
        realized.bindView(call.target.symbol,System::StateViewKind::Workspace,rhsView_);
    }
    void snapshot(){stage_.values()=data_.temperature->values();old_=stage_.values();ready_=false;}
    void boundary(){
        FDM::Scalar::applyBoundaryConditions(field_,stage_,data_.boundaries,config_.ilwBoundaryEnabled,config_.ilwAccuracyOrder);
        FDM::Scalar::applyBounds(field_,stage_,config_.bounds);
        old_=stage_.values();
    }
    void rhs(){FDM::Scalar::assembleRHS(field_,stage_,config_,rhs_);}
    void advance(){
        if(!std::isfinite(state_.dt) || state_.dt<=0)throw std::runtime_error("Thermal dt must be finite and positive.");
        Math::forInterior(field_,[&](int i,int j,int k){
            if(field_.CellFlag(i,j,k)!=FLUID_CELL)return;
            int n=stage_.getIdx(i,j,k);const double value=old_[n]+state_.dt*rhs_[n];
            if(!std::isfinite(value))throw std::runtime_error("Non-finite thermal update.");
            stage_.values()[n]=value;
        });
        FDM::Scalar::applyBoundaryConditions(field_,stage_,data_.boundaries,config_.ilwBoundaryEnabled,config_.ilwAccuracyOrder);
        const auto report=FDM::Scalar::applyBounds(field_,stage_,config_.bounds);
        if(report.violated() && config_.bounds.mode==FDM::Scalar::BoundMode::Project)
            FDM::Scalar::applyBoundaryConditions(field_,stage_,data_.boundaries,config_.ilwBoundaryEnabled,config_.ilwAccuracyOrder);
        ready_=true;
    }
    void validate()const{if(!ready_)throw std::runtime_error("Thermal publication before numerical advance.");}
    void publish(){validate();data_.temperature->values()=stage_.values();state_.distributed.markModified(config_.name);ready_=false;}
};
}
std::vector<Time::TemporalCallbacks> thermalParticipants(const System::CompiledNumericalSystem& numerics,
        const System::CompiledSolvePlan& plan,State::StateBundle& state,const FDM::SolverServices& services,
        System::StateRealization& realized) {
    std::vector<Time::TemporalCallbacks> result;
    for(const auto& call:plan.compiledProgram.steps)if(call.backendProvider==System::ThermalProvider){
        if(!services.thermalEquation || !services.thermalEquation->temperature || state.patches.size()!=1
            || (services.executionRuntime && services.executionRuntime->distributed())
            || numerics.time.recipe.id()!=FDM::TimeRecipeId::ForwardEuler || !result.empty())
            throw std::runtime_error("Unsupported thermal binding: one independent scalar, serial single patch and Forward Euler required.");
        auto op=std::make_shared<ThermalOperations>(call,state,*services.thermalEquation,realized);
        Time::TemporalCallbacks cb;cb.identity=call.source.occurrence;
        cb.stepSize=[]{return std::numeric_limits<double>::max();};
        cb.snapshot=[op]{op->snapshot();};cb.prepare=[op](const auto&){op->boundary();};
        cb.rhs=[op](const auto&){op->rhs();};cb.advance=[op](const auto&){op->advance();};
        cb.validatePublish=[op]{op->validate();};cb.publish=[op]{op->publish();};result.push_back(std::move(cb));
    }
    return result;
}
}
