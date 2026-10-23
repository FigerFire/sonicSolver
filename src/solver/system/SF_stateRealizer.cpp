/// @file SF_stateRealizer.cpp
/// @brief Resolved unknown runtime binding implementation.

#include "SF_stateRealizer.h"
#include "SF_stateRealization.h"
#include <cmath>

#include <algorithm>
#include <stdexcept>

namespace SF::System {
void StateRealization::beginStage(int index,int count,double time) {
    if (epoch_->phase!=StagePhase::Closed || index<0 || index>=count || !std::isfinite(time))
        throw std::runtime_error("Stage epoch already active or invalid index/time.");
    epoch_->index=index;epoch_->time=time;epoch_->phase=StagePhase::Preparing;
    for (auto& view:derived_) if (view->descriptor.kind==StateViewKind::Stage) view->ready=false;
}
void StateRealization::stageReady() {
    if (epoch_->phase!=StagePhase::Preparing) throw std::runtime_error("Stage preparation barrier out of order.");
    epoch_->phase=StagePhase::Reading;
}
void StateRealization::rhsReady() {
    if (epoch_->phase!=StagePhase::Reading) throw std::runtime_error("Stage RHS barrier out of order.");
    epoch_->phase=StagePhase::Advancing;
}
void StateRealization::finishStage() {
    if (epoch_->phase!=StagePhase::Advancing) throw std::runtime_error("Stage close before RHS/advance completion.");
    epoch_->phase=StagePhase::Closed;epoch_->index=-1;
}
void StateRealization::requireStageRead(int index,double time) const {
    if (epoch_->phase!=StagePhase::Reading || epoch_->index!=index || epoch_->time!=time)
        throw std::runtime_error("Stage STATE read outside its readable index/time epoch.");
}
void StateRealization::requireStageAdvance(int index,double time) const {
    if (epoch_->phase!=StagePhase::Advancing || epoch_->index!=index || epoch_->time!=time)
        throw std::runtime_error("Stage advance before complete RHS barrier or wrong index/time.");
}
std::function<double(int,int)> StateRealization::stageReader(std::string_view symbol,int stages) {
    std::vector<State::DistributedFieldView*> fields;
    for (int i=0;i<stages;++i) {
        const auto& item=view(symbol,StateViewKind::Stage,i);
        if (item.fields.size()!=1) throw std::runtime_error("Stage reader requires one bound patch: "+std::string(symbol));
        fields.push_back(item.fields.front());
    }
    return [epoch=std::weak_ptr<StageEpoch>(epoch_),fields=std::move(fields)](int cell,int component) {
        const auto scope=epoch.lock();
        if (!scope || scope->phase!=StagePhase::Reading || scope->index<0)
            throw std::runtime_error("Expired or unprepared Stage STATE view.");
        return fields.at(scope->index)->read(cell,component);
    };
}

const RealizedUnknown& StateRealization::at(std::string_view id) const {
    const auto found = std::find_if(
        unknowns_.begin(),unknowns_.end(),[id](const RealizedUnknown& value) {
            return value.descriptor && value.descriptor->id == id;
        });
    if (found == unknowns_.end()) {
        throw std::runtime_error(
            "Runtime state has no realized unknown '"+std::string(id)+"'.");
    }
    return *found;
}

namespace {
std::vector<State::DistributedFieldView*> localFields(State::StateBundle& state,const std::string& storage) {
    const auto candidates=state.distributed.select(storage,State::HaloSyncStage::None);
    std::vector<State::DistributedFieldView*> result;
    for (const auto* geometry:state.patches) {
        State::DistributedFieldView* selected=nullptr;
        for (auto* candidate:candidates) if (candidate && candidate->geometry==geometry) {
            if (selected) throw std::runtime_error("Duplicate STATE storage on active patch: "+storage);
            selected=candidate;
        }
        if (selected) result.push_back(selected);
    }
    return result;
}
}

StateRealization realizeState(
        const StateRegistry& symbols,
        const RuntimeRequirements& requirements, State::StateBundle& state,
        const std::vector<CompiledStateView>& views) {
    state.validatePatches();
    StateRealization result;result.state_=&state;result.symbols_=&symbols;
    auto demands=views;compileDerivedViews(symbols,demands);
    result.unknowns_.reserve(symbols.symbols().size());
    for (const StateSymbol& unknown : symbols.symbols()) {
        RealizedUnknown realized;
        realized.descriptor = &unknown;
        if (unknown.producer) { result.unknowns_.push_back(std::move(realized));continue; }
        if (!unknown.runtimeStorageRequired
            || unknown.storageBinding == StorageBinding::TransientWorkspace
            || unknown.storageBinding == StorageBinding::SpecializedExecutor) {
            result.unknowns_.push_back(std::move(realized));
            continue;
        }
        if (unknown.storageKey.empty()) {
            throw std::runtime_error(
                "Resolved unknown '"+unknown.id
                +"' requires runtime storage but has no storage binding key.");
        }
        realized.fields = localFields(state,unknown.storageKey);
        if (realized.fields.size() != state.patches.size()) {
            throw std::runtime_error(
                "Resolved unknown '"+unknown.id+"' expects storage '"
                +unknown.storageKey+"' on every participating patch; expected="
                +std::to_string(state.patches.size())+", actual="
                +std::to_string(realized.fields.size())+".");
        }
        for (const auto* field : realized.fields) {
            if (!field || unknown.componentOffset < 0
                || field->components < unknown.componentOffset+unknown.components) {
                throw std::runtime_error(
                    "Runtime storage '"+unknown.storageKey
                    +"' cannot represent resolved unknown '"+unknown.id+"'.");
            }
        }
        result.unknowns_.push_back(std::move(realized));
    }
    for (const auto& descriptor:demands) if (descriptor.producer && descriptor.kind==StateViewKind::Physical)
        result.bindDerived(descriptor);
    for (const auto& descriptor:demands) {
        if (descriptor.producer && descriptor.kind==StateViewKind::Physical) continue;
        RealizedStateView realized{descriptor,{}};
        if (descriptor.producer) {
            result.views_.push_back(std::move(realized));continue;
        }
        if (descriptor.owner==StateViewOwner::PhysicalState
            || (descriptor.owner==StateViewOwner::NumericalProvider
                && descriptor.kind==StateViewKind::Physical
                && symbols.at(descriptor.symbol).storageBinding==StorageBinding::ProviderDistributed)) {
            auto fields=localFields(state,descriptor.storage);
            if (fields.size()!=state.patches.size())
                throw std::runtime_error("Missing physical STATE storage: "+descriptor.storage);
            for (auto* field:fields) {
                if (!field || field->components<descriptor.componentOffset+descriptor.components)
                    throw std::runtime_error("Physical STATE view layout mismatch: "+descriptor.symbol);
            }
            realized.fields=std::move(fields);
            auto found=std::find_if(result.unknowns_.begin(),result.unknowns_.end(),
                [&](const auto& base) { return base.descriptor->id==descriptor.symbol; });
            if (found!=result.unknowns_.end() && found->fields.empty()) found->fields=realized.fields;
        } else if (descriptor.owner==StateViewOwner::CompilerWorkspace) {
            if (descriptor.location!=VariableLocation::EulerianCell)
                throw std::runtime_error("Generic STATE workspace currently requires cell location.");
            const auto& base=result.at(descriptor.symbol);
            for (std::size_t patch=0;patch<state.patches.size();++patch) {
                auto storage=std::make_unique<StateRealization::ViewStorage>();
                auto& geometry=*state.patches[patch];
                const int total=geometry.TotalSize();
                storage->values.assign((std::size_t)total*descriptor.components,0.0);
                // Working starts as a view-specific snapshot. Correction is a delta.
                if (descriptor.kind==StateViewKind::Working) {
                    if (base.fields.size()!=state.patches.size() && !base.descriptor->constantValue)
                        throw std::runtime_error("Working STATE snapshot has no physical backing: "+descriptor.symbol);
                    for (int c=0;c<descriptor.components;++c)
                        for (int cell=0;cell<total;++cell)
                            storage->values[(std::size_t)c*total+cell]=
                                base.descriptor->constantValue ? *base.descriptor->constantValue
                                    : base.fields[patch]->read(cell,base.descriptor->componentOffset+c);
                }
                storage->field=State::workspaceView(descriptor.storage,(int)patch,
                    geometry,storage->values,descriptor.components);
                realized.fields.push_back(&storage->field);result.owned_.push_back(std::move(storage));
            }
        }
        result.views_.push_back(std::move(realized));
    }
    for (const WorkspaceRequirement& requirement
         : requirements.workspaceRequirements) {
        const auto fields = state.distributed.select(
            requirement.id,State::HaloSyncStage::None);
        if (fields.size() != state.patches.size()) {
            throw std::runtime_error(
                "Required solver workspace '"+requirement.id
                +"' was not bound on every participating patch before run.");
        }
        for (const auto* field : fields) {
            if (!field || field->components != requirement.components) {
                throw std::runtime_error(
                    "Solver workspace '"+requirement.id
                    +"' has a component-layout mismatch.");
            }
        }
    }
    return result;
}

void StateRealization::bindDerived(const CompiledStateView& descriptor) {
    if (std::any_of(derived_.begin(),derived_.end(),[&](const auto& v) {
        return v->descriptor.symbol==descriptor.symbol && v->descriptor.kind==descriptor.kind && v->descriptor.stage==descriptor.stage;
    })) return;
    const auto& symbol=symbols_->at(descriptor.symbol);
    const auto& producer=*descriptor.producer;
    for (const auto& id:symbol.dependencies) if (symbols_->at(id).producer) {
        if (descriptor.kind==StateViewKind::Physical) {
            std::vector<CompiledStateView> inputs;compileDerivedViews(*symbols_,inputs);
            for (const auto& input:inputs) if (input.symbol==id) bindDerived(input);
        } else {
            const auto input=std::find_if(views_.begin(),views_.end(),[&](const auto& v) {
                return v.descriptor.symbol==id && v.descriptor.kind==descriptor.kind && v.descriptor.stage==descriptor.stage;
            });
            if(input==views_.end())throw std::runtime_error("Unsupported missing derived dependency version: "+symbol.id+" -> "+id);
            const auto requested=input->descriptor;bindDerived(requested);
        }
    }
    const bool materialized=symbol.evaluation==StateEvaluation::Materialized;
    std::vector<State::DistributedFieldView*> outputs;
    if (materialized) outputs=localFields(*state_,descriptor.storage);
    if (materialized && outputs.size()!=state_->patches.size())
        throw std::runtime_error("Unsupported missing derived output storage: "+symbol.id+"/"+toString(descriptor.kind)+" producer="+producer.id);
    std::vector<State::DistributedFieldView*> fields;
    for (std::size_t patch=0;patch<state_->patches.size();++patch) {
        std::vector<State::DistributedFieldView> inputs;
        for (const auto& id:symbol.dependencies) {
            const auto& dependency=symbols_->at(id);
            const auto& available=descriptor.kind==StateViewKind::Physical ? at(id).fields
                : view(id,descriptor.kind,descriptor.stage).fields;
            if (available.size()!=state_->patches.size())
                throw std::runtime_error("Unsupported missing dependency view: "+symbol.id+" -> "+id+"/"+toString(descriptor.kind)+" producer="+producer.id);
            auto input=*available.at(patch);
            const int offset=descriptor.kind==StateViewKind::Physical?dependency.componentOffset:0;
            // Provider Stage views are already sliced; Physical packed aliases are not.
            if (input.geometry!=state_->patches[patch] || input.location!=State::FieldLocation::Cell
                || input.components<offset+dependency.components)
                throw std::runtime_error("Invalid derived input layout: "+symbol.id+" -> "+id);
            input.read=[read=input.read,offset](int cell,int component){return read(cell,offset+component);};
            input.components=dependency.components;
            input.write=[id](int,int,double){throw std::runtime_error("Execution contract: derived input is read-only: "+id);};
            inputs.push_back(std::move(input));
        }
        auto bound=std::make_shared<DerivedView>();bound->descriptor=descriptor;
        if (materialized) bound->field=*outputs.at(patch);
        else {
            auto& v=bound->field;v.name=symbol.id;v.blockId=(int)patch;v.geometry=state_->patches[patch];
            v.components=symbol.components;v.haloDepth=v.geometry->NG();v.exchange=State::ExchangeKind::None;
        }
        if (bound->field.components!=symbol.components || bound->field.geometry!=state_->patches[patch]
            || bound->field.location!=State::FieldLocation::Cell)
            throw std::runtime_error("Invalid derived output layout: "+symbol.id);
        auto* distributed=&state_->distributed;
        // As with distributed field handles, registration is complete before binding.
        std::vector<const unsigned long long*> revisions;
        const unsigned long long* outputRevision=nullptr;
        if (materialized) {
            for (const auto& input:inputs) revisions.push_back(&distributed->generation(input.name,input.blockId));
            outputRevision=&distributed->generation(bound->field.name,bound->field.blockId);
            bound->versions.resize(revisions.size());
        }
        bound->update=producer.bind(*state_->patches[patch],inputs,bound->field);
        if (materialized!=static_cast<bool>(bound->update) || !bound->field.read)
            throw std::runtime_error("Invalid STATE producer update/access contract: "+symbol.id+" producer="+producer.id);
        const auto weak=std::weak_ptr<DerivedView>(bound);
        const auto epoch=std::weak_ptr<StageEpoch>(epoch_);
        const std::string context=symbol.id+"/"+toString(descriptor.kind)+" stage="+std::to_string(descriptor.stage)+" producer="+producer.id;
        bound->field.requireReadable=[weak,epoch,inputs,revisions,outputRevision,materialized,distributed,context] {
            const auto value=weak.lock();const auto scope=epoch.lock();
            if (!value || !scope) throw std::runtime_error("Execution contract: expired derived view "+context);
            if (value->descriptor.kind==StateViewKind::Stage && (scope->phase!=StagePhase::Reading || scope->index!=value->descriptor.stage))
                throw std::runtime_error("Execution contract: wrong Stage epoch "+context);
            if (materialized) {
                if (!value->ready) throw std::runtime_error("Execution contract: unready derived dependencies "+context);
                if (*outputRevision!=value->outputVersion)
                    throw std::runtime_error("Execution contract: overwritten derived output "+context);
                for (std::size_t i=0;i<inputs.size();++i) {
                    if (inputs[i].requireReadable) inputs[i].requireReadable();
                    if (*revisions[i]!=value->versions[i])
                        throw std::runtime_error("Execution contract: stale derived dependency "+inputs[i].name+" "+context);
                }
            }
            if (value->descriptor.haloRequired && distributed->needsExchange(value->field.name,1))
                throw std::runtime_error("Execution contract: unpublished derived halo "+context);
        };
        bound->field.read=[validate=bound->field.requireReadable,read=bound->field.read](int cell,int component) {
            validate();return read(cell,component);
        };
        if (bound->update) {
            bound->update=[weak,epoch,revisions,outputRevision,compute=bound->update,distributed,context] {
                const auto value=weak.lock();const auto scope=epoch.lock();
                if (!value || !scope) throw std::runtime_error("Execution contract: expired derived update "+context);
                if (value->descriptor.kind==StateViewKind::Stage && (scope->phase!=StagePhase::Preparing || scope->index!=value->descriptor.stage))
                    throw std::runtime_error("Execution contract: update outside Stage preparation "+context);
                value->ready=false;
                for(std::size_t i=0;i<revisions.size();++i)value->versions[i]=*revisions[i];
                compute();
                for(std::size_t i=0;i<revisions.size();++i)if(*revisions[i]!=value->versions[i])
                    throw std::runtime_error("Execution contract: producer wrote its dependencies "+context);
                distributed->markModified(value->field.name,value->field.blockId);
                value->outputVersion=*outputRevision;value->ready=true;
            };
        }
        bound->field.write=[context](int,int,double){throw std::runtime_error("Execution contract: read-only derived view "+context);};
        bound->field.validate();fields.push_back(&bound->field);derived_.push_back(std::move(bound));
    }
    if (descriptor.kind==StateViewKind::Physical) {
        for (auto& base:unknowns_) if (base.descriptor->id==descriptor.symbol) base.fields=fields;
    }
    const auto found=std::find_if(views_.begin(),views_.end(),[&](const auto& v){
        return v.descriptor.symbol==descriptor.symbol && v.descriptor.kind==descriptor.kind && v.descriptor.stage==descriptor.stage;
    });
    if (found==views_.end()) views_.push_back({descriptor,std::move(fields)});
    else found->fields=std::move(fields);
}

void StateRealization::bindDerivedViews() {
    // Binding resolves handles only. It never runs a closure or inserts an Update.
    std::vector<CompiledStateView> pending;
    for (const auto& view:views_) if (view.descriptor.producer && view.fields.empty()) pending.push_back(view.descriptor);
    for (const auto& descriptor:pending) bindDerived(descriptor);
}
void StateRealization::updateView(std::string_view symbol,StateViewKind kind,int stage) {
    bool found=false;
    for (const auto& item:derived_) if (item->descriptor.symbol==symbol && item->descriptor.kind==kind && item->descriptor.stage==stage) {
        if (!item->update) throw std::runtime_error("Execution contract: no explicit update for readonly closure "+std::string(symbol));
        found=true;item->update();
    }
    if (!found) throw std::runtime_error("Execution contract: unbound derived update "+std::string(symbol)+"/"+toString(kind));
}

const RealizedStateView& StateRealization::view(std::string_view symbol,
        StateViewKind kind,int stage) const {
    const auto found=std::find_if(views_.begin(),views_.end(),[&](const auto& value) {
        return value.descriptor.symbol==symbol && value.descriptor.kind==kind
            && value.descriptor.stage==stage;
    });
    if (found==views_.end()) throw std::runtime_error("Unrequested STATE view: "+std::string(symbol)+"/"+toString(kind));
    if (found->fields.empty()) throw std::runtime_error("Numerical provider has not bound STATE view: "+std::string(symbol)+"/"+toString(kind));
    return *found;
}

bool StateRealization::requestsView(std::string_view symbol,StateViewKind kind,int stage) const {
    return std::any_of(views_.begin(),views_.end(),[&](const auto& value) {
        return value.descriptor.symbol==symbol && value.descriptor.kind==kind
            && value.descriptor.stage==stage;
    });
}

void StateRealization::bindView(std::string_view symbol,StateViewKind kind,
        State::DistributedFieldView& field,int stage) {
    field.validate();
    const auto found=std::find_if(views_.begin(),views_.end(),[&](const auto& value) {
        return value.descriptor.symbol==symbol && value.descriptor.kind==kind
            && value.descriptor.stage==stage;
    });
    if (found==views_.end() || found->descriptor.owner!=StateViewOwner::NumericalProvider
        || found->descriptor.producer)
        throw std::runtime_error("Provider cannot bind unrequested/owned STATE view: "+std::string(symbol));
    if (found->descriptor.components!=field.components)
        throw std::runtime_error("Provider STATE view component mismatch: "+std::string(symbol));
    const auto same=std::find_if(found->fields.begin(),found->fields.end(),
        [&](auto* item) { return item->blockId==field.blockId; });
    if (same!=found->fields.end()) {
        if (*same!=&field) throw std::runtime_error("Duplicate provider authority for STATE view: "+std::string(symbol));
    } else found->fields.push_back(&field);
    if (kind==StateViewKind::Physical) {
        const auto base=std::find_if(unknowns_.begin(),unknowns_.end(),
            [&](const auto& value) { return value.descriptor->id==symbol; });
        if (base!=unknowns_.end()) base->fields=found->fields;
    }
}

void validatePhysicalResources(const CompiledExecutionProgram& program,const std::vector<ResourceView>& bound) {
    std::vector<const ResourceViewRequirement*> required;
    for (const auto& call:program.steps) for (const auto& item:call.physicalResources) {
        if (item.id.empty() || item.owner.empty() || item.type==typeid(void))
            throw std::runtime_error("Invalid physical resource declaration: "+item.id);
        if (std::any_of(required.begin(),required.end(),[&](const auto* previous){return previous->id==item.id;}))
            throw std::runtime_error("Duplicate physical resource identity: "+item.id);
        required.push_back(&item);
    }
    for (const auto* item:required) {
        int matches=0;
        for (const auto& view:bound) {
            const auto& actual=view.requirement();
            if (actual.id!=item->id) continue;
            if (actual.owner!=item->owner || actual.type!=item->type || actual.writable!=item->writable)
                throw std::runtime_error("Physical resource type/owner/access conflict: "+item->id+" owner="+item->owner);
            view.requireAvailable();++matches;
        }
        if (matches!=1) throw std::runtime_error("Missing/duplicate physical resource binding: "+item->id+" owner="+item->owner);
    }
    for (const auto& view:bound)
        if (std::none_of(required.begin(),required.end(),[&](const auto* item){return item->id==view.requirement().id;}))
            throw std::runtime_error("Unconsumed physical resource binding: "+view.requirement().id);
}

} // namespace SF::System
