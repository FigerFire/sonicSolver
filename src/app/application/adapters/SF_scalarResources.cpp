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
        const auto existing=state.distributed.select(symbol.storageKey,State::HaloSyncStage::None);
        State::DistributedFieldView* backing=nullptr;
        for(auto* view:existing)if(view->blockId==block) {
            if(backing || view->geometry!=&geometry || view->components!=symbol.components || !view->read || !view->write)
                throw std::runtime_error("Conflicting existing STATE backing shape/layout/owner: "+id);
            view->validate();if(view->requireReadable)view->requireReadable();backing=view;
        }
        if(backing)return;
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
        if (required.owner!="native.scalar" && required.owner!="model.scalar") continue;
        if(required.type!=typeid(FDM::ScalarEquationData)
            || required.writable || !occurrences.insert(required.id).second)
            throw std::runtime_error("Unsupported/duplicate physical resource contract: "+required.id+" owner="+required.owner);
        if(!required.fields)throw std::runtime_error("Missing compiled field resource requirements: "+required.id);
        const auto& need=*required.fields;
        for(const auto& id:need.states) {
            const auto& symbol=system.executableSystem.state.at(id);
            if(!symbol.producer && symbol.storageBinding==System::StorageBinding::NamedDistributed)allocate(id);
        }
        std::vector<const System::ResourceView*> provided;
        for(const auto& view:alreadyBound)if(view.requirement().id==required.id)provided.push_back(&view);
        if(provided.size()>1)throw std::runtime_error("Duplicate field resource backing: "+required.id);
        std::shared_ptr<const FDM::ScalarEquationData> port;
        if(!provided.empty()) {
            const auto& actual=provided.front()->requirement();
            if(actual.owner!=required.owner || actual.type!=required.type || actual.writable!=required.writable)
                throw std::runtime_error("Field port type/owner/access conflict: "+required.id);
            provided.front()->requireAvailable();port=provided.front()->read<FDM::ScalarEquationData>();
        } else {
            if(required.owner=="model.scalar")throw std::runtime_error("Missing model scalar boundary resource: "+required.id+" target="+need.target);
            auto created=std::make_shared<FDM::ScalarEquationData>();created->target=need.target;
            if(need.localOnly && !descriptor(need.target).boundaries.empty())
                throw std::runtime_error("Unsupported local boundary/halo resource: "+required.id);
            if(need.boundaryRequired)*created=boundaryData(geometry,descriptor(need.target),system.executableSystem.state.at(need.target).components);
            resources.views.push_back(System::ResourceView::bind(required,created,*created,required.owner));
            resources.ports.push_back(created);port=std::move(created);
        }
        if(port->target!=need.target)throw std::runtime_error("Field port target binding mismatch: "+required.id+" requires "+need.target);
        for(const auto& source:need.sources)if(!port->sources.count(source) || !port->sources.at(source))
            throw std::runtime_error("Missing known-source port binding: occurrence="+required.id+" target="+need.target+" source="+source);
        if(need.boundaryRequired)for(int axis=0;axis<3;++axis) {
            const int size=axis==0?geometry.NX():axis==1?geometry.NY():geometry.NZ();
            if(size<=1)continue;
            for(int side=0;side<2;++side) {
                const auto& bc=port->boundary[2*axis+side];
                if(bc.kind!=FDM::ScalarBoundaryKind::FixedValue && bc.kind!=FDM::ScalarBoundaryKind::ZeroGradient)
                    throw std::runtime_error("Missing scalar boundary backing: "+need.target);
                if(bc.kind==FDM::ScalarBoundaryKind::FixedValue && !(system.executableSystem.state.at(need.target).components==1?bool(bc.value):bool(bc.componentValue)))
                    throw std::runtime_error("Missing component-correct fixedValue callback: "+need.target);
            }
        }
        resources.instances.push_back({required.id,port});
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
