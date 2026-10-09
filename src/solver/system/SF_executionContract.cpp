#include "SF_executionContract.h"
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
namespace SF::System {
namespace {
bool matches(const ExecutionScope& node,const std::string& selector) {
    if (selector.rfind("first-of:",0)==0) {
        const auto split=selector.find('|',9);
        return matches(node,selector.substr(9,split==std::string::npos?split:split-9))
            || (split!=std::string::npos && matches(node,selector.substr(split+1)));
    }
    if (selector=="kind:Commit") return node.kind==ExecutionKind::Commit;
    if (selector.rfind("equation:",0)==0) {
        if (node.kind==ExecutionKind::EquationCall) return node.step.equation==selector.substr(9);
        return std::any_of(node.children.begin(),node.children.end(),[&](const auto& child){return matches(child,selector);});
    }
    return node.id==selector || (node.kind==ExecutionKind::EquationCall && node.step.occurrence==selector);
}
struct Revision { int value=0,halo=0; std::string writer="initial state"; };
}
void composeDefaultPlacement(const ExecutionProgram& program,ExecutionScope& root) {
    const auto compose=[&](const auto& self,ExecutionScope& scope,bool top)->void {
        const auto index=[&](const std::string& selector) {
            std::vector<std::size_t> found;
            for (std::size_t i=0;i<scope.children.size();++i)
                if (matches(scope.children[i],selector)) found.push_back(i);
            if (found.empty() || (found.size()!=1 && selector.rfind("first-of:",0)!=0))
                throw std::runtime_error("HOW placement requires exactly one sibling '"+selector+"' in scope '"+scope.id+"'; no cross-scope ordering is inferred.");
            return found.front();
        };
        std::vector<std::pair<std::size_t,std::size_t>> edges;
        std::set<std::size_t> inserted;
        for (const auto& requirement:program.requirements) {
            if (requirement.scope.empty()? !top:requirement.scope!=scope.id) continue;
            const auto subject=index(requirement.occurrence);inserted.insert(subject);
            for (const auto& before:requirement.before) edges.push_back({subject,index(before)});
            for (const auto& after:requirement.after) edges.push_back({index(after),subject});
        }
        if (!inserted.empty()) {
            // Keep the selected recipe backbone, including whole nested control scopes.
            std::size_t previous=scope.children.size();
            for (std::size_t i=0;i<scope.children.size();++i) if (!inserted.count(i)) {
                if (previous<scope.children.size()) edges.push_back({previous,i});previous=i;
            }
            std::vector<std::size_t> order;std::set<std::size_t> done;
            while (order.size()<scope.children.size()) {
                std::vector<std::size_t> ready,placed;
                for (std::size_t i=0;i<scope.children.size();++i) {
                    if (done.count(i)) continue;
                    if (std::any_of(edges.begin(),edges.end(),[&](auto edge){return edge.second==i && !done.count(edge.first);})) continue;
                    ready.push_back(i);if (inserted.count(i)) placed.push_back(i);
                }
                if (ready.empty()) {
                    std::string unresolved;
                    for (std::size_t i=0;i<scope.children.size();++i) if (!done.count(i))
                        unresolved+=(scope.children[i].id.empty()?scope.children[i].step.equation:scope.children[i].id)+" ";
                    throw std::runtime_error("HOW dependency cycle in default placement scope '"+scope.id+"'; blocked occurrences: "+unresolved);
                }
                if (ready.size()>1) throw std::runtime_error("Ambiguous default HOW placement in scope '"+scope.id+"'; select an explicit HOW or declare the numerical placement relation.");
                const auto next=placed.empty()?ready.front():placed.front();done.insert(next);order.push_back(next);
            }
            auto original=std::move(scope.children);
            for (auto i:order) scope.children.push_back(std::move(original[i]));
        }
        for (auto& child:scope.children) self(self,child,false);
    };
    compose(compose,root,true);
}
void validatePlacement(const ExecutionProgram& program,const ExecutionScope& root) {
    std::map<std::string,const ExecutionScope*> scopes;
    const auto collect=[&](const auto& self,const ExecutionScope& scope)->void {
        if (!scope.id.empty()) {
            if (!scopes.emplace(scope.id,&scope).second) throw std::runtime_error("Duplicate HOW scope: "+scope.id);
        }
        for (const auto& child:scope.children) self(self,child);
    };
    collect(collect,root);
    std::map<const ExecutionScope*,std::vector<std::pair<std::size_t,std::size_t>>> edges;
    for (const auto& requirement:program.requirements) {
        const auto scope=requirement.scope.empty()?&root:
            scopes.count(requirement.scope)?scopes.at(requirement.scope):nullptr;
        if (!scope) throw std::runtime_error("Missing placement scope: "+requirement.scope);
        const auto index=[&](const std::string& selector) {
            std::vector<std::size_t> found;
            for (std::size_t i=0;i<scope->children.size();++i) if (matches(scope->children[i],selector)) found.push_back(i);
            if (found.empty() || (found.size()!=1 && selector.rfind("first-of:",0)!=0)) throw std::runtime_error("HOW placement requires exactly one sibling '"+selector+"' in scope '"+scope->id+"'; no cross-scope ordering is inferred.");
            return found.front();
        };
        const auto subject=index(requirement.occurrence);
        for (const auto& before:requirement.before) edges[scope].push_back({subject,index(before)});
        for (const auto& after:requirement.after) edges[scope].push_back({index(after),subject});
    }
    for (const auto& entry:edges) {
        const auto& children=entry.first->children;
        std::vector<int> color(children.size());std::vector<std::size_t> stack;
        const auto visit=[&](const auto& self,std::size_t vertex)->void {
            if (color[vertex]==2) return;
            if (color[vertex]==1) {
                std::string chain;
                for (auto i:stack) chain+=children[i].id+" -> ";
                throw std::runtime_error("HOW dependency cycle: "+chain+children[vertex].id);
            }
            color[vertex]=1;stack.push_back(vertex);
            for (auto edge:entry.second) if (edge.first==vertex) self(self,edge.second);
            stack.pop_back();color[vertex]=2;
        };
        for (std::size_t i=0;i<children.size();++i) visit(visit,i);
        for (auto edge:entry.second) if (edge.first>=edge.second)
            throw std::runtime_error("HOW placement violation in scope '"+entry.first->id+"': '"+children[edge.first].id+"' must precede '"+children[edge.second].id+"'; declared topology was not reordered.");
    }
}
std::string missingCapability(const CompiledExecutionProgram& program,
        const std::vector<std::string>& context,const CompiledEquationCall& consumer,
        const std::vector<CapabilityBinding>& boundContext) {
    std::set<std::string> available(context.begin(),context.end());
    auto bound=boundContext;
    for (const auto& call:program.steps) {
        available.insert(call.capabilities.begin(),call.capabilities.end());
        bound.insert(bound.end(),call.boundCapabilities.begin(),call.boundCapabilities.end());
    }
    for (const auto& item:bound) available.insert(item.name);
    for (const auto& requirement:consumer.capabilityRequirements) {
        if (!requirement.when.empty() && !available.count(requirement.when)) continue;
        const bool scoped=!requirement.equation.empty() || !requirement.fluidPort.empty() || !requirement.boundary.empty();
        const bool satisfied=scoped?std::any_of(bound.begin(),bound.end(),[&](const auto& item) {
            return item.name==requirement.required && item.equation==requirement.equation
                && item.fluidPort==requirement.fluidPort && item.boundary==requirement.boundary;
        }):available.count(requirement.required)>0;
        if (!satisfied) return "Occurrence '"+consumer.source.occurrence+"' ("+consumer.equationMethod+") requires capability '"+requirement.required+"'"
            +(requirement.when.empty()?"":" in context '"+requirement.when+"'")
            +(scoped?" bound to equation='"+requirement.equation+"', fluid-port='"+requirement.fluidPort+"', boundary='"+requirement.boundary+"'":"")+": "+requirement.reason;
    }
    return {};
}
void validateDataFlow(const StateRegistry& state,const CompiledExecutionProgram& program) {
    // Validate dependency graphs only for derived views, not coupled equation unknowns.
    std::map<std::string,int> color;std::vector<std::string> stack;
    const auto visit=[&](const auto& self,const std::string& id)->void {
        if (color[id]==2) return;
        if (color[id]==1) {
            std::string chain;for (const auto& item:stack) chain+=item+" -> ";
            throw std::runtime_error("Derived STATE dependency cycle: "+chain+id);
        }
        const auto& symbol=state.at(id);color[id]=1;stack.push_back(id);
        for (const auto& dependency:symbol.dependencies) if (state.contains(dependency)) self(self,dependency);
        stack.pop_back();color[id]=2;
    };
    for (const auto& symbol:state.symbols()) visit(visit,symbol.id);
    const auto overlaps=[&](const std::string& a,const std::string& b) {
        if (a==b) return true;
        if (!state.contains(a) || !state.contains(b)) return false;
        const auto& x=state.at(a);const auto& y=state.at(b);
        if (x.derivation!=StateDerivation::None || y.derivation!=StateDerivation::None
            || x.storageKey.empty() || x.storageKey!=y.storageKey
            || x.storageBinding!=y.storageBinding
            || x.location!=y.location) return false;
        return x.componentOffset<y.componentOffset+y.components && y.componentOffset<x.componentOffset+x.components;
    };
    std::map<std::string,Revision> revisions;
    std::map<std::string,std::map<std::string,int>> published;
    std::map<std::string,std::string> invalidLazy;
    std::set<std::string> snapshots;
    std::set<std::string> persistent;
    std::map<std::string,int> frozen;
    std::vector<std::string> loops;
    for (const auto& symbol:state.symbols()) {
        revisions[symbol.id]={};
        for (const auto& dependency:symbol.dependencies) published[symbol.id][dependency]=0;
    }
    int epoch=0;
    const auto version=[&](const auto& self,const std::string& id)->int {
        int result=revisions[id].value;
        if (state.contains(id)) for (const auto& dependency:state.at(id).dependencies) result+=self(self,dependency);
        return result;
    };
    const auto writer=[&](const auto& self,const std::string& id)->std::pair<int,std::string> {
        auto result=std::make_pair(revisions[id].value,revisions[id].writer);
        if (state.contains(id)) for (const auto& dependency:state.at(id).dependencies) {
            const auto upstream=self(self,dependency);if (upstream.first>result.first) result=upstream;
        }
        return result;
    };
    const auto check=[&](const auto& self,const std::string& id,const std::string& consumer,bool halo)->void {
        const auto& symbol=state.at(id);
        if (symbol.evaluation==StateEvaluation::Lazy && invalidLazy.count(id))
            throw std::runtime_error("Lazy STATE invalidation missing: consumer '"+consumer+"' reads '"+id+"'; writer '"+invalidLazy[id]+"' did not invalidate its versioned cache.");
        if (halo && revisions[id].halo!=revisions[id].value)
            throw std::runtime_error("STATE halo unavailable: consumer '"+consumer+"' reads '"+id+"'; writer '"+revisions[id].writer+"' has no halo publication.");
        for (const auto& dependency:symbol.dependencies) {
            self(self,dependency,consumer,halo);
            if (symbol.evaluation==StateEvaluation::Materialized && published[id][dependency]!=version(version,dependency))
                throw std::runtime_error("Stale STATE: consumer '"+consumer+"' requires current '"+id+"'; invalidating writer '"+writer(writer,dependency).second+"' changed '"+dependency+"'; missing publication of '"+id+"'.");
        }
    };
    const auto execute=[&](const auto& self,const ExecutionScope& node,bool stage)->void {
        if (node.kind==ExecutionKind::EquationCall) {
            const auto at=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& call){return call.source.occurrence==node.step.occurrence;});
            if (at==program.steps.end()) throw std::runtime_error("Missing compiled occurrence: "+node.step.occurrence);
            // Complete fused groups are one algebraic operation, not sequential partial writes.
            if (!at->fusionKey.empty() && at!=program.steps.begin() && (at-1)->fusionKey==at->fusionKey) return;
            std::vector<StateUse> uses;
            std::vector<StateEffect> effects;
            const auto members=at->fusionKey.empty()?1:at->fusionMembers.size();
            for (std::size_t i=0;i<members;++i) {
                const auto& member=*(at+i);
                uses.insert(uses.end(),member.stateUses.begin(),member.stateUses.end());
                for (const auto& effect:member.stateEffects) {
                    const auto previous=std::find_if(effects.begin(),effects.end(),[&](const auto& item){return item.symbol==effect.symbol;});
                    if (previous==effects.end()) effects.push_back(effect);
                    else if (previous->publish!=effect.publish || previous->halo!=effect.halo || previous->invalidatesLazy!=effect.invalidatesLazy)
                        throw std::runtime_error("Conflicting fused STATE effects: "+effect.symbol);
                }
            }
            for (const auto& workspace:at->workspaceRequires) if (!snapshots.count(workspace))
                throw std::runtime_error("EquationCall '"+node.step.equation+"' requires unavailable workspace '"+workspace+"' in this scope.");
            for (const auto& use:uses) {
                if (use.version==StateVersion::Frozen || use.version==StateVersion::OldTime || use.version==StateVersion::Lagged) {
                    const auto key=use.snapshot.empty()?use.symbol:use.snapshot;
                    if (!snapshots.count(key)) throw std::runtime_error("Unavailable "+std::string(toString(use.version))+" STATE '"+key+"' for consumer '"+node.step.occurrence+"'.");
                    if (use.perIteration && loops.empty()) throw std::runtime_error("Frozen STATE requires an enclosing iteration: "+key);
                    const auto binding=(use.perIteration ? loops.back()+"/iteration/" : std::string("physical/"))+node.id+":"+key;
                    if (use.version==StateVersion::Frozen && frozen.count(binding) && frozen[binding]!=revisions[key].value)
                        throw std::runtime_error("Frozen STATE overwritten within loop: "+key+" for "+node.step.occurrence);
                    frozen[binding]=revisions[key].value;
                } else {
                    if (!state.contains(use.symbol)) throw std::runtime_error("Unknown STATE use: "+use.symbol);
                    const bool commonStage=at->temporalResidual && !program.temporalParticipants.empty()
                        && std::any_of(program.temporalParticipants.begin(),program.temporalParticipants.end(),[&](const auto& p) {
                            return std::find(p.targets.begin(),p.targets.end(),use.symbol)!=p.targets.end();
                        });
                    if (use.version==StateVersion::Stage && !stage && !commonStage)
                        throw std::runtime_error("Unsupported stage STATE use outside declared StageLoop: "+node.step.occurrence);
                    check(check,use.symbol,node.step.occurrence,use.halo);
                }
            }
            for (const auto& effect:effects) {
                if (!state.contains(effect.symbol)) {
                    if (!effect.publish) throw std::runtime_error("Unknown STATE write: "+effect.symbol);
                    snapshots.insert(effect.symbol);revisions[effect.symbol].value=++epoch;continue;
                }
                const auto depends=[&](const auto& walk,const std::string& id)->bool {
                    if (overlaps(id,effect.symbol)) return true;
                    for (const auto& source:state.at(id).dependencies) if (walk(walk,source)) return true;
                    return false;
                };
                for (const auto& symbol:state.symbols()) if (symbol.evaluation==StateEvaluation::Lazy && depends(depends,symbol.id)) {
                    if (effect.invalidatesLazy) invalidLazy.erase(symbol.id);
                    else invalidLazy[symbol.id]=node.step.occurrence;
                }
                ++epoch;
                for (const auto& symbol:state.symbols()) if (overlaps(symbol.id,effect.symbol)) {
                    auto& revision=revisions[symbol.id];revision.value=epoch;revision.writer=node.step.occurrence;
                    if (effect.halo) revision.halo=revision.value;
                }
                if (effect.publish) for (const auto& dependency:state.at(effect.symbol).dependencies)
                    published[effect.symbol][dependency]=version(version,dependency);
            }
            for (std::size_t i=0;i<members;++i) {
                const auto& member=*(at+i);
                for (const auto& workspace:member.workspaceProvides) {
                    snapshots.insert(workspace);
                    if (std::none_of(effects.begin(),effects.end(),[&](const auto& effect){return effect.symbol==workspace;}))
                        revisions[workspace].value=++epoch;
                }
                persistent.insert(member.stepWorkspaces.begin(),member.stepWorkspaces.end());
            }
            return;
        }
        const auto savedSnapshots=snapshots;const auto savedFrozen=frozen;
        const int passes=(node.kind==ExecutionKind::Loop || node.kind==ExecutionKind::StageLoop) ? std::min(node.repetitions,2):1;
        const bool iterative=node.kind==ExecutionKind::Loop || node.kind==ExecutionKind::StageLoop;
        const auto label=node.id.empty()?std::string("<root-loop>"):node.id;
        if (iterative) loops.push_back(label);
        for (int i=0;i<passes;++i) {
            if (iterative) for (auto at=frozen.begin();at!=frozen.end();) {
                if (at->first.rfind(label+"/iteration/",0)==0) at=frozen.erase(at);else ++at;
            }
            for (const auto& child:node.children) self(self,child,stage || node.kind==ExecutionKind::StageLoop);
        }
        if (iterative) loops.pop_back();
        // Workspaces produced inside a control scope do not escape into siblings.
        if (node.kind==ExecutionKind::Loop || node.kind==ExecutionKind::StageLoop) {snapshots=savedSnapshots;snapshots.insert(persistent.begin(),persistent.end());frozen=savedFrozen;}
    };
    for (int physicalStep=0;physicalStep<2;++physicalStep) {
        snapshots.clear();persistent.clear();frozen.clear();
        // Lagged versions are immutable step-entry inputs, including legally zero initialized force.
        for (const auto& call:program.steps) for (const auto& use:call.stateUses)
            if (use.version==StateVersion::Lagged && state.contains(use.symbol)
                && state.at(use.symbol).availableVersion==StateVersion::Lagged && use.snapshot.empty())
                snapshots.insert(use.snapshot.empty()?use.symbol:use.snapshot);
        execute(execute,program.root,false);
    }
}
}
