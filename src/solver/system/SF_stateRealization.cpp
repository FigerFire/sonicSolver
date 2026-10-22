/// @file SF_stateRealization.cpp
/// @brief STATE REALIZATION 编译实现（role 分组只读 executable system）。

#include "SF_stateRealization.h"

#include <algorithm>
#include <stdexcept>
#include <map>

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

void compileDerivedViews(const StateRegistry& state,std::vector<CompiledStateView>& views) {
    std::map<std::string,int> color;
    const auto visit=[&](const auto& self,const StateSymbol& symbol)->void {
        if (color[symbol.id]==2) return;
        if (color[symbol.id]==1) throw std::runtime_error("Invalid derived STATE dependency cycle: "+symbol.id);
        color[symbol.id]=1;
        if (symbol.producer) {
            const auto& p=*symbol.producer;
            if (p.id.empty() || p.versions.empty() || symbol.constantValue
                || symbol.evaluation==StateEvaluation::Direct)
                throw std::runtime_error("Invalid STATE producer declaration: "+symbol.id+" producer="+p.id);
            if (!p.bind) throw std::runtime_error("Unsupported missing STATE producer implementation: "+symbol.id+" producer="+p.id);
            Math::requireCompatible(symbol.signature,p.expression.signature(),"derived output "+symbol.id+" producer="+p.id);
            if (symbol.componentOffset!=0 || (symbol.evaluation==StateEvaluation::Materialized && symbol.storageKey.empty()))
                throw std::runtime_error("Unsupported derived output storage/slice: "+symbol.id+" producer="+p.id);
            state.validateSymbols(p.expression);
            const auto validate=[&](const auto& walk,const FormulaExpr& e)->void {
                if(e.kind()==FormulaExpr::Kind::Symbol && std::find(symbol.dependencies.begin(),symbol.dependencies.end(),e.name())==symbol.dependencies.end())
                    throw std::runtime_error("Invalid undeclared derived dependency: "+symbol.id+" -> "+e.name()+" producer="+p.id);
                for(const auto& a:e.arguments())walk(walk,a);
            };
            validate(validate,p.expression);
            for (const auto& dependency:symbol.dependencies) {
                if (!state.contains(dependency)) throw std::runtime_error("Invalid missing derived dependency: "+symbol.id+" -> "+dependency+" producer="+p.id);
                self(self,state.at(dependency));
            }
        }
        color[symbol.id]=2;
    };
    for (const auto& symbol:state.symbols()) visit(visit,symbol);
    for (const auto& symbol:state.symbols()) if (symbol.producer) {
        if (std::none_of(views.begin(),views.end(),[&](const auto& v){return v.symbol==symbol.id && v.kind==StateViewKind::Physical;})) {
            CompiledStateView view;view.symbol=symbol.id;view.storage=symbol.storageKey;
            view.components=symbol.components;view.componentOffset=symbol.componentOffset;
            view.location=symbol.location;view.ownership=symbol.ownership;view.owner=StateViewOwner::NumericalProvider;
            views.push_back(std::move(view));
        }
    }
    for (auto& view:views) if (state.contains(view.symbol) && state.at(view.symbol).producer) {
        const auto& symbol=state.at(view.symbol);const auto& p=*symbol.producer;
        if ((view.kind==StateViewKind::Stage ? view.stage<0 : view.stage!=-1))
            throw std::runtime_error("Invalid STATE view stage: "+view.symbol+"/"+toString(view.kind));
        // A Working/Correction workspace is not a projection of Physical STATE.
        // Its declared compiler/provider owner retains the existing storage authority.
        if (!view.producer && (view.owner==StateViewOwner::CompilerWorkspace
            || view.kind==StateViewKind::Working || view.kind==StateViewKind::Correction
            || view.kind==StateViewKind::Workspace)) continue;
        if (view.producer && view.producer!=symbol.producer)
            throw std::runtime_error("Invalid conflicting STATE producer: "+view.symbol+"/"+toString(view.kind));
        if (std::find(p.versions.begin(),p.versions.end(),view.kind)==p.versions.end()
            || symbol.location!=VariableLocation::EulerianCell || view.location!=symbol.location)
            throw std::runtime_error("Unsupported STATE producer view: "+view.symbol+"/"+toString(view.kind)+" producer="+p.id);
        if (view.components!=symbol.components)
            throw std::runtime_error("Invalid derived STATE view layout: "+view.symbol);
        view.producer=symbol.producer;
    }
}

CompiledStateRealization compileStateRealization(
        const StateRegistry& state,const std::vector<ConstraintDescriptor>& constraints) {
    std::vector<CompiledStateView> checked;compileDerivedViews(state,checked);
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
            const bool update=symbol.producer && symbol.evaluation==StateEvaluation::Materialized
                && target.producer==symbol.producer;
            if (symbol.readOnly() && !update)
                throw std::runtime_error("Physical STATE target is a read-only closure: "+target.symbol);
            if (symbol.storageKey.empty())
                throw std::runtime_error("Physical STATE target has no storage contract: "+target.symbol);
            view.storage=symbol.storageKey;
            view.componentOffset=symbol.componentOffset;
            view.owner=update || symbol.storageBinding==StorageBinding::SpecializedExecutor
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
    for(const auto& participant:program.temporalParticipants)for(const auto& name:participant.coefficientStates) {
        const auto& symbol=state.at(name);
        for(int stage=0;stage<stages;++stage) {
            CompiledStateView view;view.symbol=name;view.kind=StateViewKind::Stage;view.stage=stage;
            view.storage=symbol.storageKey+".stage-coefficient";view.components=symbol.components;view.componentOffset=symbol.componentOffset;
            view.location=symbol.location;view.ownership=symbol.ownership;view.owner=StateViewOwner::NumericalProvider;program.stateViews.push_back(std::move(view));
        }
    }
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
        if (stages>1 || call.synchronousStages) program.stateViews.push_back(old);
        if (!call.residualWorkspace.empty()) {
            auto residual=old;residual.kind=StateViewKind::Workspace;
            residual.storage=call.residualWorkspace;residual.componentOffset=0;
            program.stateViews.push_back(std::move(residual));
        }
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
