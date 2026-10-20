/// @file SF_scalarResources.cpp
/// @brief 将已编译的命名标量需求绑定到唯一物理数组；时间循环仍由 executeEquations 执行。
#include "app/application/execution/SF_execution.h"
#include "app/application/output/SF_fields.h"
#include "core/field/SF_scalarField.h"
#include "core/state/SF_state.h"
#include "core/interfaces/SF_log.h"
#include "infrastructure/execution/SF_localExecutionRuntime.h"
#include "infrastructure/io/SF_caseConfig.h"
#include "solver/system/SF_scalarMethod.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
namespace SF::Application::Execution::Detail {
namespace {
double finiteValue(const Model::Parameters& value,const std::string& source) {
    if(!value.is_number() || !std::isfinite(value.get<double>()))throw std::runtime_error(source+": finite scalar value required.");
    return value.get<double>();
}
FDM::ScalarEquationData boundaryData(const Field& field,const Model::FieldDescriptor& descriptor) {
    FDM::ScalarEquationData data;data.target=descriptor.symbol;
    const int g=field.NG(),shape[3]={field.NX(),field.NY(),field.NZ()};
    for(const auto& item:descriptor.boundaries.items()) {
        const auto source="fields/"+descriptor.name+"/boundaries/"+item.key();
        Model::Schema{{{"type","string",true},{"value","number"}},false}.validate(item.value(),source);
        const auto type=item.value().at("type").get<std::string>();
        if(type!="fixedValue" && type!="zeroGradient")throw std::runtime_error(source+": native scalar provider supports fixedValue/zeroGradient only.");
        const auto kind=type=="fixedValue"?FDM::ScalarBoundaryKind::FixedValue:FDM::ScalarBoundaryKind::ZeroGradient;
        if(type=="fixedValue" && !item.value().contains("value"))throw std::runtime_error(source+": fixedValue requires value.");
        const double value=type=="fixedValue"?finiteValue(item.value().at("value"),source):0;
        const auto found=field.getAllSets().find(item.key());
        if(found==field.getAllSets().end())throw std::runtime_error(source+": mesh boundary set absent.");
        std::set<int> physical,covered;
        for(int cell:found->second) {
            int i,j,k;field.getIJK(cell,i,j,k);
            if(i>=g && i<g+shape[0] && j>=g && j<g+shape[1] && k>=g && k<g+shape[2])physical.insert(cell);
        }
        for(int axis=0;axis<3;++axis)if(shape[axis]>1)for(int side=0;side<2;++side) {
            std::set<int> face;
            for(int cell:physical) {int i,j,k;field.getIJK(cell,i,j,k);const int at[3]={i,j,k};if(at[axis]==g+(side?shape[axis]-1:0))face.insert(cell);}
            if(face.size()!=static_cast<std::size_t>(shape[(axis+1)%3])*shape[(axis+2)%3])continue;
            auto& bc=data.boundary[2*axis+side];
            if(bc.kind!=FDM::ScalarBoundaryKind::Unspecified)throw std::runtime_error(source+": overlapping scalar boundary face.");
            bc.kind=kind;if(kind==FDM::ScalarBoundaryKind::FixedValue)bc.value=[value](double,double,double,double){return value;};
            covered.insert(face.begin(),face.end());
        }
        if(physical.empty() || covered!=physical)throw std::runtime_error(source+": complete axis-aligned boundary face(s) required.");
    }
    for(int axis=0;axis<3;++axis)if(shape[axis]>1)for(int side=0;side<2;++side)
        if(data.boundary[2*axis+side].kind==FDM::ScalarBoundaryKind::Unspecified)
            throw std::runtime_error("Missing scalar boundary backing: "+descriptor.symbol+" axis="+std::to_string(axis)+" side="+std::to_string(side));
    return data;
}
}
int executeScalarEquations(Field& geometry,ResultWriter& writer,const FDM::SolverConfig& config,
        const System::ResolvedSimulationSystem& system,const CaseConfig& input,bool initialOnly) {
    if(geometry.NVar()!=0 || geometry.stateModel())throw std::runtime_error("Native scalar resources require geometry-only Field without Q/EOS.");
    if(!input.modelDescription)throw std::runtime_error("Native scalar resources require registered field declarations.");
    std::vector<std::unique_ptr<ScalarField>> arrays;
    std::vector<std::shared_ptr<const FDM::ScalarEquationData>> ports;
    State::StateBundle state;state.patches={&geometry};state.time=input.time.startTime;
    FDM::SolverServices services;::SF::Execution::LocalRuntime runtime;services.executionRuntime=&runtime;
    for(const auto& symbol:system.executableSystem.state.symbols()) {
        if(symbol.storageBinding!=System::StorageBinding::NamedDistributed || symbol.components!=1)
            throw std::runtime_error("Unsupported native physical storage binding: "+symbol.id);
        const auto& fields=input.modelDescription->fields;
        const auto found=std::find_if(fields.begin(),fields.end(),[&](const auto& f){return f.symbol==symbol.id;});
        if(found==fields.end())throw std::runtime_error("Missing registered storage backing: "+symbol.id);
        auto physical=std::make_unique<ScalarField>();physical->setupLike(geometry,symbol.storageKey,finiteValue(found->initial,"fields/"+found->name+"/initial"));
        for(const auto& item:found->sets.items()) {
            const auto source="fields/"+found->name+"/initial/"+item.key();
            Model::Schema{{{"type","string",true},{"value","number",true}},false}.validate(item.value(),source);
            if(item.value().at("type")!="fixedValue")throw std::runtime_error(source+": unsupported scalar initialization.");
            const double value=finiteValue(item.value().at("value"),source);
            const auto set=geometry.getAllSets().find(item.key());
            if(set==geometry.getAllSets().end() || set->second.empty())throw std::runtime_error(source+": initialization mesh set absent.");
            for(int cell:set->second)physical->values().at(cell)=value;
        }
        state.distributed.add(State::scalarView(symbol.storageKey,0,geometry,*physical));arrays.push_back(std::move(physical));
    }
    for(const auto& call:system.solvePlan.compiledProgram.steps) {
        if(call.backendProvider!=System::ScalarOps::Provider)throw std::runtime_error("Unsupported native resource provider: "+call.backendProvider);
        const auto& contract=std::any_cast<const System::CompiledScalarEquation&>(call.providerContract);
        if(!contract.sources.empty())throw std::runtime_error("Native scalar known-source callback binding is Unsupported: "+call.source.occurrence);
        const auto& fields=input.modelDescription->fields;
        const auto descriptorAt=std::find_if(fields.begin(),fields.end(),[&](const auto& f){return f.symbol==call.target.symbol;});
        if(descriptorAt==fields.end())throw std::runtime_error("Missing registered scalar port backing: "+call.target.symbol);
        const auto& descriptor=*descriptorAt;
        auto port=std::make_shared<FDM::ScalarEquationData>(boundaryData(geometry,descriptor));
        services.scalarInstances.push_back({call.source.occurrence,port});ports.push_back(std::move(port));
    }
    const auto output=Output::registeredStateFields(system.executableSystem.state,state,geometry);
    auto saveStep=[&](int step){writer.saveStep(geometry,step,state.time,output);};
    auto saveTime=[&](double time){writer.save(geometry,time,output);};
    SF::broadcast("Native resources: ","geometry-only, Q=0, EOS=none, scalar arrays="+std::to_string(arrays.size()));
    if(input.writeInitial || initialOnly) {if(input.time.writeByStep)saveStep(0);else saveTime(state.time);}
    if(initialOnly)return 0;
    Time::DriverCallbacks callbacks;callbacks.saveStep=saveStep;callbacks.saveTime=saveTime;
    callbacks.report=[&](int step,const Time::AdvanceResult& result){SF::broadcast("Scalar step: ","step="+std::to_string(step)+", time="+std::to_string(result.time)+", dt="+std::to_string(result.dt));};
    return executeEquations(config,system,state,services,input.time,callbacks);
}
}
