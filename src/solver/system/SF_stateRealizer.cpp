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
        const ResolvedSimulationSystem& system, State::StateBundle& state) {
    state.validatePatches();
    StateRealization result;
    result.unknowns_.reserve(system.unknowns.size());
    for (const UnknownDescriptor& unknown : system.unknowns) {
        RealizedUnknown realized;
        realized.descriptor = &unknown;
        if (!unknown.runtimeStorageRequired
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
    for (const WorkspaceRequirement& requirement
         : system.workspaceRequirements) {
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

} // namespace SF::System
