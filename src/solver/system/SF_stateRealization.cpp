/// @file SF_stateRealization.cpp
/// @brief STATE REALIZATION 编译实现（role 分组只读 executable system）。

#include "SF_stateRealization.h"

#include <algorithm>

namespace SF::System {
namespace {

RealizedRoleGroup roleGroup(const UnknownDescriptor& unknown) {
    RealizedRoleGroup group;
    group.id = unknown.id;
    group.role = unknown.role;
    group.components = unknown.components;
    group.storageKey = unknown.storageKey;
    group.componentOffset = unknown.componentOffset;
    group.nameSpace = unknown.nameSpace;
    return group;
}

bool isTransportBinding(StorageBinding binding) {
    return binding == StorageBinding::PackedDistributed
        || binding == StorageBinding::NamedDistributed;
}

void append(std::vector<RealizedRoleGroup>& target,
            const UnknownDescriptor& unknown) {
    target.push_back(roleGroup(unknown));
}

} // namespace

CompiledStateRealization compileStateRealization(
        const ExecutableEquationSystem& system) {
    CompiledStateRealization realization;
    for (const UnknownDescriptor& unknown : system.unknowns) {
        switch (unknown.role) {
            case UnknownRole::Primary:
            case UnknownRole::Transported:
                if (isTransportBinding(unknown.storageBinding)) {
                    append(realization.transported,unknown);
                } else {
                    append(realization.derived,unknown);
                }
                break;
            case UnknownRole::Multiplier:
                append(realization.multipliers,unknown);
                break;
            case UnknownRole::Algebraic:
            case UnknownRole::Derived:
                append(realization.derived,unknown);
                break;
        }
    }
    for (const ConstraintDescriptor& constraint : system.constraints) {
        realization.constraints.push_back(constraint.id);
        if (!constraint.multiplierUnknown.empty()) {
            realization.multiplierConstraints = true;
        }
    }
    const auto hasTransported = [&](const std::string& id) {
        return std::any_of(
            realization.transported.begin(),realization.transported.end(),
            [&](const RealizedRoleGroup& group) { return group.id == id; });
    };
    const auto hasRole = [&](const std::string& id, UnknownRole role) {
        return std::any_of(
            system.unknowns.begin(),system.unknowns.end(),
            [&](const UnknownDescriptor& unknown) {
                return unknown.id == id && unknown.role == role;
            });
    };
    realization.conservativeTransportedMass = hasTransported("rho");
    realization.conservativeMomentum =
        hasTransported("rhoU") || hasTransported("U");
    realization.pressureMultiplier = hasRole("p",UnknownRole::Multiplier);
    realization.thermodynamicPressure = hasRole("p",UnknownRole::Derived);
    for (const RealizedRoleGroup& group : realization.transported) {
        if (group.id.rfind("phaseMass",0) == 0) {
            realization.phaseTransportedState = true;
            break;
        }
    }
    return realization;
}

} // namespace SF::System
