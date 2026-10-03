/// @file SF_stateRealization.cpp
/// @brief STATE REALIZATION 编译实现（role 分组只读 executable system）。

#include "SF_stateRealization.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System {
namespace {

RealizedRoleGroup roleGroup(const StateSymbol& unknown) {
    RealizedRoleGroup group;
    group.id = unknown.id;
    group.role = unknown.role;
    group.components = unknown.components;
    group.storageKey = unknown.storageKey;
    group.componentOffset = unknown.componentOffset;
    group.constantValue = unknown.constantValue;
    group.nameSpace = unknown.nameSpace;
    return group;
}

bool isTransportBinding(StorageBinding binding) {
    return binding == StorageBinding::PackedDistributed
        || binding == StorageBinding::NamedDistributed
        || binding == StorageBinding::ProviderDistributed;
}

void append(std::vector<RealizedRoleGroup>& target,
            const StateSymbol& unknown) {
    target.push_back(roleGroup(unknown));
}

} // namespace

CompiledStateRealization compileStateRealization(
        const StateRegistry& state,const std::vector<ConstraintDescriptor>& constraints) {
    CompiledStateRealization realization;
    for (const StateSymbol& unknown : state.symbols()) {
        const bool constraintMultiplier=std::any_of(constraints.begin(),constraints.end(),
            [&](const auto& constraint) { return constraint.multiplierUnknown==unknown.id; });
        if (unknown.role==StateRole::Algebraic && constraintMultiplier) {
            append(realization.multipliers,unknown);
            continue;
        }
        switch (unknown.role) {
            case StateRole::Primary:
            case StateRole::Transported:
                if (isTransportBinding(unknown.storageBinding)) {
                    append(realization.transported,unknown);
                } else {
                    append(realization.derived,unknown);
                }
                break;
            case StateRole::Multiplier:
                append(realization.multipliers,unknown);
                break;
            case StateRole::Algebraic:
            case StateRole::Derived:
                append(realization.derived,unknown);
                break;
        }
    }
    for (const ConstraintDescriptor& constraint : constraints) {
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
    const auto hasRole = [&](const std::string& id, StateRole role) {
        return std::any_of(
            state.symbols().begin(),state.symbols().end(),
            [&](const StateSymbol& unknown) {
                return unknown.id == id && unknown.role == role;
            });
    };
    realization.conservativeTransportedMass = hasTransported("rho");
    realization.conservativeMomentum = hasTransported("rhoU");
    realization.pressureMultiplier = hasRole("p",StateRole::Multiplier)
        || (state.isSolution("p") && std::any_of(constraints.begin(),constraints.end(),
            [](const auto& constraint) { return constraint.multiplierUnknown=="p"; }));
    realization.thermodynamicPressure = hasRole("p",StateRole::Derived);
    for (const RealizedRoleGroup& group : realization.transported) {
        if (group.id.rfind("phaseMass",0) == 0) {
            realization.phaseTransportedState = true;
            break;
        }
    }
    return realization;
}

CompiledStateView realizeTarget(const StateRegistry& state,CompiledTarget& target) {
    CompiledStateView view;
    view.symbol=target.symbol;
    const bool workspace=target.kind==TargetKind::Workspace;
    if (workspace) {
        view.location=target.workspaceLocation;view.ownership=target.workspaceOwnership;
        if (!target.resources.empty()) view.components=target.resources.front().components;
    }
    if (!workspace) {
        const auto& symbol=state.at(target.symbol);
        view.components=symbol.components;
        view.location=symbol.location;
        view.ownership=symbol.ownership;
        if (target.kind==TargetKind::Physical) {
            if (symbol.constantValue || symbol.derivation!=StateDerivation::None)
                throw std::runtime_error("Physical STATE target is a read-only closure: "+target.symbol);
            if (symbol.storageKey.empty())
                throw std::runtime_error("Physical STATE target has no storage contract: "+target.symbol);
            view.storage=symbol.storageKey;
            view.componentOffset=symbol.componentOffset;
            view.owner=symbol.storageBinding==StorageBinding::SpecializedExecutor
                || symbol.storageBinding==StorageBinding::ProviderDistributed
                || symbol.location!=VariableLocation::EulerianCell
                ? StateViewOwner::NumericalProvider : StateViewOwner::PhysicalState;
        }
    }
    if (target.kind!=TargetKind::Physical) {
        view.kind=target.kind==TargetKind::Working ? StateViewKind::Working
            : target.kind==TargetKind::Correction ? StateViewKind::Correction : StateViewKind::Workspace;
        view.storage=target.workspace.empty()
            ? "state."+target.symbol+"."+toString(view.kind) : target.workspace;
        target.workspace=view.storage;
        view.owner=target.viewOwner;
    }
    target.resources={{target.symbol,view.storage,view.componentOffset,view.components,
        ResourceAccessMode::Write,false,SynchronizationRequirement::WriteOwned}};
    return view;
}

void realizeTemporalViews(const StateRegistry& state,CompiledExecutionProgram& program,
        int stages) {
    for (const auto& call:program.steps) {
        if (!call.temporalResidual) continue;
        const auto& symbol=state.at(call.source.target.symbol);
        CompiledStateView old;
        old.symbol=symbol.id;old.kind=StateViewKind::OldTime;
        if (stages>1 && call.oldTimeWorkspace.empty())
            throw std::runtime_error("Temporal provider did not declare OldTime storage: "+call.source.occurrence);
        if (!call.publishesStageToPhysicalTarget && call.stageWorkspace.empty())
            throw std::runtime_error("Temporal provider did not declare Stage publication: "+call.source.occurrence);
        old.storage=call.oldTimeWorkspace;old.componentOffset=symbol.componentOffset;
        old.components=symbol.components;old.location=symbol.location;
        old.ownership=symbol.ownership;old.owner=StateViewOwner::NumericalProvider;
        if (stages>1) program.stateViews.push_back(old);
        for (int stage=0;stage<stages;++stage) {
            auto view=old;view.kind=StateViewKind::Stage;view.stage=stage;
            // The explicit backend publishes each stage into the physical packed
            // authority. This is a time-qualified alias, not another copy of Q.
            view.storage=call.publishesStageToPhysicalTarget
                ? call.target.resources.front().storage : call.stageWorkspace+"."+std::to_string(stage);
            program.stateViews.push_back(std::move(view));
        }
    }
}

} // namespace SF::System
