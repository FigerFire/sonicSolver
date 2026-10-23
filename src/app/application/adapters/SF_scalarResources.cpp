/// @file SF_scalarResources.cpp
/// @brief 将已编译的命名标量需求绑定到唯一物理数组；时间循环仍由 executeEquations 执行。
#include "SF_scalarResources.h"
#include "solver/system/SF_resolvedSimulationSystem.h"
#include "solver/system/SF_stateRealizer.h"
#include "core/state/SF_stateBundle.h"
#include "core/field/SF_scalarField.h"
#include "core/state/SF_state.h"
#include "core/interfaces/SF_log.h"
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
std::vector<double> fieldValue(const Model::Parameters& value,int components,const std::string& source) {
    if(components==1)return {finiteValue(value,source)};
    if(!value.is_array() || value.size()!=static_cast<std::size_t>(components))throw std::runtime_error(source+": explicit vector components required.");
    std::vector<double> result;for(const auto& c:value)result.push_back(finiteValue(c,source));return result;
}
FDM::ScalarEquationData boundaryData(const Field& field,const Model::FieldDescriptor& descriptor,int components) {
    FDM::ScalarEquationData data;data.target=descriptor.symbol;
    const int g=field.NG(),shape[3]={field.NX(),field.NY(),field.NZ()};
    for(const auto& item:descriptor.boundaries.items()) {
        const auto source="fields/"+descriptor.name+"/boundaries/"+item.key();
        Model::Schema{{{"type","string",true},{"value","any"}},false}.validate(item.value(),source);
        const auto type=item.value().at("type").get<std::string>();
        if(type!="fixedValue" && type!="zeroGradient")throw std::runtime_error(source+": native scalar provider supports fixedValue/zeroGradient only.");
        const auto kind=type=="fixedValue"?FDM::ScalarBoundaryKind::FixedValue:FDM::ScalarBoundaryKind::ZeroGradient;
        if(type=="fixedValue" && !item.value().contains("value"))throw std::runtime_error(source+": fixedValue requires value.");
        const auto value=type=="fixedValue"?fieldValue(item.value().at("value"),components,source):std::vector<double>(components,0);
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
            bc.kind=kind;if(kind==FDM::ScalarBoundaryKind::FixedValue) {
                if(components==1)bc.value=[value](double,double,double,double){return value[0];};
                else bc.componentValue=[value](double,double,double,double,int c){return value.at(c);};
            }
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
ScalarResources bindScalarResources(Field& geometry,State::StateBundle& state,
        const System::ResolvedSimulationSystem& system,const CaseConfig& input,int block,
        const std::vector<System::ResourceView>& alreadyBound) {
    ScalarResources resources;
    const auto descriptor=[&](const std::string& id)->const Model::FieldDescriptor& {
        if (!input.modelDescription) throw std::runtime_error("Native scalar resources require registered field declarations.");
        const auto& fields=input.modelDescription->fields;
        const auto found=std::find_if(fields.begin(),fields.end(),[&](const auto& f){return f.symbol==id;});
        if(found==fields.end())throw std::runtime_error("Missing registered storage backing: "+id);
        return *found;
    };
    std::map<std::string,std::string> allocated;
    std::vector<State::DistributedFieldView> pending;
    std::set<std::string> occurrences;
    const auto allocate=[&](const std::string& id) {
        const auto& symbol=system.executableSystem.state.at(id);
        if (symbol.storageBinding!=System::StorageBinding::NamedDistributed || (symbol.components!=1 && symbol.components!=3)
            || symbol.componentOffset!=0)
            throw std::runtime_error("Unsupported native physical storage binding: "+symbol.id);
        const auto [previous,inserted]=allocated.emplace(symbol.storageKey,id);
        if (!inserted) {
            if (previous->second!=id) throw std::runtime_error("Conflicting native initialization owners for storage: "+symbol.storageKey);
            return;
        }
        if (!state.distributed.select(symbol.storageKey,State::HaloSyncStage::None).empty())
            throw std::runtime_error("Native resource already has a physical owner: "+symbol.storageKey);
        const auto& found=descriptor(id);
        const auto initial=fieldValue(found.initial,symbol.components,"fields/"+found.name+"/initial");
        std::vector<ScalarField*> components;
        for(int c=0;c<symbol.components;++c) {
            auto physical=std::make_unique<ScalarField>();
            physical->setupLike(geometry,symbol.storageKey,initial[c]);
            for(const auto& item:found.sets.items()) {
                const auto source="fields/"+found.name+"/initial/"+item.key();
                Model::Schema{{{"type","string",true},{"value","any",true}},false}.validate(item.value(),source);
                if(item.value().at("type")!="fixedValue")throw std::runtime_error(source+": unsupported field initialization.");
                const auto value=fieldValue(item.value().at("value"),symbol.components,source);
                const auto set=geometry.getAllSets().find(item.key());
                if(set==geometry.getAllSets().end() || set->second.empty())throw std::runtime_error(source+": initialization mesh set absent.");
                for(int cell:set->second)physical->values().at(cell)=value[c];
            }
            components.push_back(physical.get());resources.arrays.push_back(std::move(physical));
        }
        pending.push_back(State::scalarComponentsView(symbol.storageKey,block,geometry,components,geometry.NG(),State::HaloSyncStage::None));
    };
    for(const auto& call:system.solvePlan.compiledProgram.steps) for(const auto& required:call.physicalResources) {
        if (required.owner!="native.scalar") continue;
        if(required.owner!="native.scalar" || required.type!=typeid(FDM::ScalarEquationData)
            || required.writable || required.id!=call.source.occurrence || !occurrences.insert(required.id).second)
            throw std::runtime_error("Unsupported/duplicate physical resource contract: "+required.id+" owner="+required.owner);
        const auto* scalarContract=std::any_cast<System::CompiledScalarEquation>(&call.providerContract);
        if (!scalarContract) throw std::runtime_error("Physical scalar resource lacks its compiled implementation contract: "+required.id);
        const auto& contract=*scalarContract;
        if(!contract.sources.empty())throw std::runtime_error("Native scalar known-source callback binding is Unsupported: "+required.id);
        allocate(call.target.symbol);
        for (const auto& id:contract.stageReads)
            if (system.executableSystem.state.at(id).storageBinding==System::StorageBinding::NamedDistributed) allocate(id);
        const bool stencil=std::any_of(system.solvePlan.compiledProgram.steps.begin(),system.solvePlan.compiledProgram.steps.end(),[&](const auto& consumer){
            const auto* c=std::any_cast<System::CompiledScalarEquation>(&consumer.providerContract);
            return c && ((c->conservativeTransport && c->target==call.target.symbol)
                || std::any_of(c->diffusions.begin(),c->diffusions.end(),[&](const auto& d){return d.input==call.target.symbol;}));
        });
        auto port=std::make_shared<FDM::ScalarEquationData>();port->target=call.target.symbol;
        if(stencil)*port=boundaryData(geometry,descriptor(call.target.symbol),contract.components);
        resources.views.push_back(System::ResourceView::bind(required,port,*port,required.owner));
        resources.instances.push_back({required.id,port});resources.ports.push_back(std::move(port));
    }
    auto completeBindings=alreadyBound;
    completeBindings.insert(completeBindings.end(),resources.views.begin(),resources.views.end());
    System::validatePhysicalResources(system.solvePlan.compiledProgram,completeBindings);
    // Publish borrowed storage only after every requested array/port validates.
    // A failed bind cannot leave dangling registrations in the active state.
    for(auto& view:pending)state.distributed.add(std::move(view));
    for (const auto& instance:resources.instances)
        SF::broadcast("Bound physical resource: ",instance.occurrence+" owner=native.scalar target="+instance.data.lock()->target);
    return resources;
}
}
