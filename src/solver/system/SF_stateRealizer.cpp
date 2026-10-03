/// @file SF_stateRealizer.cpp
/// @brief Resolved unknown runtime binding implementation.

#include "SF_stateRealizer.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System {

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

StateRealization realizeState(
        const StateRegistry& symbols,
        const RuntimeRequirements& requirements, State::StateBundle& state,
        const std::vector<CompiledStateView>& views) {
    state.validatePatches();
    StateRealization result;
    result.unknowns_.reserve(symbols.symbols().size());
    for (const StateSymbol& unknown : symbols.symbols()) {
        RealizedUnknown realized;
        realized.descriptor = &unknown;
        if (unknown.derivation!=StateDerivation::None) {
            for (std::size_t patch=0;patch<state.patches.size();++patch) {
                auto storage=std::make_unique<StateRealization::ViewStorage>();
                auto& view=storage->field;
                view.name=unknown.id;view.blockId=(int)patch;view.geometry=state.patches[patch];
                view.components=unknown.components;view.haloDepth=view.geometry->NG();
                view.exchange=State::ExchangeKind::None;
                view.read=[field=view.geometry,kind=unknown.derivation](int cell,int component) {
                    int i=0,j=0,k=0;field->getIJK(cell,i,j,k);
                    const auto value=field->thermodynamicState(i,j,k);
                    switch (kind) {
                    case StateDerivation::Velocity:return value.velocity.at((std::size_t)component);
                    case StateDerivation::Pressure:return value.pressure;
                    case StateDerivation::Temperature:return value.temperature;
                    case StateDerivation::Enthalpy:return (value.internalEnergyDensity+value.pressure)/value.density;
                    case StateDerivation::DynamicViscosity:return value.dynamicViscosity;
                    case StateDerivation::KinematicViscosity:return value.dynamicViscosity/value.density;
                    case StateDerivation::ThermalConductivity:return value.thermalConductivity;
                    case StateDerivation::None:break;
                    }
                    throw std::runtime_error("Unsupported STATE derivation.");
                };
                view.write=[](int,int,double) { throw std::runtime_error("Cannot write a derived STATE view."); };
                realized.fields.push_back(&view);result.owned_.push_back(std::move(storage));
            }
            result.unknowns_.push_back(std::move(realized));continue;
        }
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
        realized.fields = state.distributed.select(
            unknown.storageKey,State::HaloSyncStage::None);
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
    for (const auto& descriptor:views) {
        RealizedStateView realized{descriptor,{}};
        if (descriptor.owner==StateViewOwner::PhysicalState
            || (descriptor.owner==StateViewOwner::NumericalProvider
                && descriptor.kind==StateViewKind::Physical
                && symbols.at(descriptor.symbol).storageBinding==StorageBinding::ProviderDistributed)) {
            auto fields=state.distributed.select(descriptor.storage,State::HaloSyncStage::None);
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
    if (found==views_.end() || found->descriptor.owner!=StateViewOwner::NumericalProvider)
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

} // namespace SF::System
