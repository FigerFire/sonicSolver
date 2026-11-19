/// @file SF_stateRealization.cpp
/// @brief STATE REALIZATION 编译实现（role 分组只读 executable system）。

#include "SF_stateRealization.h"
#include "core/interfaces/SF_scalarEquationData.h"

#include <algorithm>
#include <stdexcept>
#include <map>
#include <set>

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

int requiredPackedComponents(const StateRegistry& state,const CompiledExecutionProgram& program,std::string_view storage) {
    int components=0;
    std::set<std::string> visited;
    const auto inspect=[&](const auto& self,const std::string& id)->void {
        if(!visited.insert(id).second)return;
        const auto& symbol=state.at(id);
        if(symbol.storageBinding==StorageBinding::PackedDistributed && symbol.storageKey==storage) {
            if(symbol.componentOffset<0 || symbol.components<=0)
                throw std::runtime_error("Unsupported packed STATE backing: "+id+" storage="+symbol.storageKey);
            components=std::max(components,symbol.componentOffset+symbol.components);
        }
        if(symbol.producer)for(const auto& dependency:symbol.dependencies)self(self,dependency);
    };
    for(const auto& view:program.stateViews) {
        if(state.contains(view.symbol))inspect(inspect,view.symbol);
        else if(view.kind!=StateViewKind::Workspace || view.owner==StateViewOwner::PhysicalState)
            throw std::runtime_error("Unregistered physical STATE backing: "+view.symbol);
        // Numerical face/block workspace targets need no physical-state symbol.
    }
    return components;
}

void resolveCallInputs(const StateRegistry& state,const std::vector<EquationCall*>& steps) {
    std::set<std::string> configured;
    for(const auto& binding:state.callInputs()) {
        std::vector<EquationCall*> matches;
        for(auto* step:steps)if((binding.reference==ReferenceKind::Step?step->occurrence:step->equation)==binding.call)matches.push_back(step);
        if(matches.size()!=1)throw std::runtime_error("STATE call input binding is missing or ambiguous: "+binding.call+"; use one occurrence ID.");
        auto* step=matches.front();
        if(!configured.insert(step->occurrence).second || !step->inputs.empty())
            throw std::runtime_error("Duplicate STATE input binding for occurrence: "+step->occurrence);
        step->inputs=binding.reads;
    }
}

ResourceViewRequirement equationFieldResource(std::string id,const std::string& target,
        const std::vector<StateUse>& reads,std::vector<std::string> boundaries,
        std::vector<std::string> sources,bool local) {
    FieldDataRequirement fields;fields.target=target;fields.states={target};fields.localOnly=local;
    for(const auto& read:reads)fields.states.push_back(read.symbol);
    fields.boundaries=std::move(boundaries);fields.sources=std::move(sources);
    return {std::move(id),"native.scalar",typeid(FDM::ScalarEquationData),false,std::move(fields)};
}

void compileFieldResources(const StateRegistry& state,CompiledExecutionProgram& program) {
    std::set<std::string> boundaries;
    for(const auto& call:program.steps)for(const auto& resource:call.physicalResources)if(resource.fields) {
        const auto& field=*resource.fields;state.at(field.target);
        for(const auto& id:field.states)state.at(id);
        for(const auto& id:field.boundaries){state.at(id);boundaries.insert(id);}
    }
    for(auto& call:program.steps)for(auto& resource:call.physicalResources)if(resource.fields) {
        auto& field=*resource.fields;field.boundaryRequired=boundaries.count(field.target)!=0;
        if(field.localOnly && field.boundaryRequired)
            throw std::runtime_error("Unsupported local boundary/halo resource: occurrence="+call.source.occurrence+" equation="+call.source.equation+" method="+call.equationMethod+" target="+field.target);
    }
}

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
                || symbol.storageBinding==StorageBinding::TransientWorkspace
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

std::vector<FieldPublication> fieldPublications(const StateRegistry& state,const CompiledExecutionProgram& program) {
    std::vector<FieldPublication> publications;
    for(const auto& call:program.steps) {
        // Legacy whole-model operations publish their own storage. They do not
        // transfer authority for the arrays to independent numerical calls.
        if(call.legacyProvider)continue;
        for(const auto& effect:call.stateEffects) {
            if(!effect.publish || !state.contains(effect.symbol))continue;
            const auto& symbol=state.at(effect.symbol);
            if(symbol.storageBinding!=StorageBinding::ProviderDistributed
                || symbol.evaluation!=StateEvaluation::Materialized)continue;
            const auto physical=std::find_if(program.stateViews.begin(),program.stateViews.end(),[&](const auto& view) {
                return view.symbol==symbol.id && view.kind==StateViewKind::Physical;
            });
            CompiledStateView backing;
            if(physical!=program.stateViews.end())backing=*physical;
            else {
                // A grouped operation can publish more than its primary target.
                CompiledTarget target;
                target.symbol=symbol.id;
                target.producer=symbol.producer;
                backing=realizeTarget(state,target);
            }
            const auto collision=std::find_if(publications.begin(),publications.end(),[&](const auto& item) {
                const auto& previous=item.backing;
                return previous.storage==backing.storage
                    && previous.componentOffset<backing.componentOffset+backing.components
                    && backing.componentOffset<previous.componentOffset+previous.components;
            });
            if(collision!=publications.end())
                throw std::runtime_error("Overlapping field publications for storage="+backing.storage+
                    ": "+collision->occurrence+" / "+call.source.occurrence);
            publications.push_back({std::move(backing),call.source.occurrence,
                call.stageLocal?StateVersion::Stage:StateVersion::Current});
        }
    }
    return publications;
}

void realizeTemporalViews(const StateRegistry& state,CompiledExecutionProgram& program,
        int stages) {
    for(const auto& participant:program.steps)if(participant.coefficientVersion==StateVersion::Stage || participant.coefficientVersion==StateVersion::Current)for(const auto& name:participant.coefficientStates) {
        const auto& symbol=state.at(name);
        if(symbol.storageBinding==StorageBinding::TransientWorkspace && std::none_of(program.stateViews.begin(),program.stateViews.end(),[&](const auto& v){return v.symbol==name && v.kind==StateViewKind::Physical;})) {
            CompiledStateView physical;physical.symbol=name;physical.kind=StateViewKind::Physical;
            physical.storage=symbol.storageKey;physical.components=symbol.components;physical.componentOffset=symbol.componentOffset;
            physical.location=symbol.location;physical.ownership=symbol.ownership;physical.owner=StateViewOwner::NumericalProvider;
            program.stateViews.push_back(std::move(physical));
        }
        for(int stage=0;participant.coefficientVersion==StateVersion::Stage && stage<stages;++stage) {
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
    for(const auto& call:program.steps)if(call.stageLocal && call.coefficientProvides.empty())for(int stage=0;stage<stages;++stage){const auto& s=state.at(call.target.symbol);CompiledStateView v;v.symbol=s.id;v.kind=StateViewKind::Stage;v.stage=stage;v.storage=call.stageWorkspace;v.components=s.components;v.location=s.location;v.ownership=s.ownership;v.owner=StateViewOwner::NumericalProvider;program.stateViews.push_back(std::move(v));}
    for(const auto& call:program.steps)for(const auto& read:call.stateUses)
        if(read.version==StateVersion::Stage && !state.at(read.symbol).producer && !state.at(read.symbol).constantValue)
            for(int stage=0;stage<stages;++stage)
                if(std::none_of(program.stateViews.begin(),program.stateViews.end(),[&](const auto& v){return v.symbol==read.symbol && v.kind==StateViewKind::Stage && v.stage==stage;}))
                    throw std::runtime_error("Missing Stage publisher: step="+call.source.occurrence+" equation="+call.source.equation+" target="+call.target.symbol+" input="+read.symbol+" stage="+std::to_string(stage));
    for(const auto& call:program.steps)for(const auto& read:call.stateUses)if(read.version==StateVersion::Frozen) {
        const auto& symbol=state.at(read.symbol);
        const bool published=std::any_of(program.steps.begin(),program.steps.end(),[&](const auto& producer) {
            return producer.source.stepScoped && (std::find(producer.coefficientStates.begin(),producer.coefficientStates.end(),read.symbol)!=producer.coefficientStates.end()
                || (producer.equationMethod=="LocalAffine" && producer.target.symbol==read.symbol));
        });
        if(!published || !read.perIteration || !read.snapshot.empty() || read.halo || read.kind!=TargetKind::Physical
            || symbol.evaluation!=StateEvaluation::Materialized)
            throw std::runtime_error("Frozen coefficient needs an explicit Step materialized publisher and perIteration read: "+read.symbol);
        if(std::none_of(program.stateViews.begin(),program.stateViews.end(),[&](const auto& view){return view.symbol==read.symbol && view.kind==StateViewKind::Frozen;})) {
            CompiledStateView view;view.symbol=symbol.id;view.kind=StateViewKind::Frozen;view.storage=symbol.storageKey;
            view.components=symbol.components;view.location=symbol.location;view.ownership=symbol.ownership;view.owner=StateViewOwner::NumericalProvider;
            program.stateViews.push_back(std::move(view));
        }
    }
    const auto stageDependency=[&](const auto& self,const std::string& id,int stage)->void {
        if(std::any_of(program.stateViews.begin(),program.stateViews.end(),[&](const auto& v){return v.symbol==id && v.kind==StateViewKind::Stage && v.stage==stage;}))return;
        const auto& s=state.at(id);
        if(!s.producer)throw std::runtime_error("Missing Stage publisher for input: "+id);
        // Dependency cycles are validated by compileDerivedViews; this traversal has a local visiting guard.
        CompiledStateView v;v.symbol=id;v.kind=StateViewKind::Stage;v.stage=stage;v.storage=s.storageKey;
        v.components=s.components;v.owner=StateViewOwner::NumericalProvider;v.location=s.location;v.ownership=s.ownership;
        program.stateViews.push_back(std::move(v));
        for(const auto& dep:s.dependencies)self(self,dep,stage);
    };
    for(const auto& call:program.steps) {
        const auto inputs=[&](const auto& reads){for(const auto& read:reads)
            if(read.version==StateVersion::Stage && state.contains(read.symbol) && state.at(read.symbol).producer)
                for(int stage=0;stage<stages;++stage)stageDependency(stageDependency,read.symbol,stage);
        };
        inputs(call.stateUses);
    }

}

} // namespace SF::System
