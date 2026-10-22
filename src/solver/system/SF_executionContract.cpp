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
std::vector<const CompiledEquationCall*> equationConsumers(
        const ExecutableEquationSystem& system,const CompiledExecutionProgram& program,
        std::string_view equation) {
    const auto references=[](const CompiledEquationCall& call,std::string_view id) {
        return call.source.equation==id
            || std::find(call.sourceMathInputs.begin(),call.sourceMathInputs.end(),id)!=call.sourceMathInputs.end()
            || std::any_of(call.calls.begin(),call.calls.end(),[&](const auto& ref) {return ref.equation==id;});
    };
    const bool direct=std::any_of(program.steps.begin(),program.steps.end(),[&](const auto& call){return references(call,equation);});
    std::vector<const CompiledEquationCall*> result;
    for (const auto& call:program.steps) {
        bool uses=references(call,equation);
        if (!direct) for (const auto& constraint:system.constraints)
            if (std::find(constraint.equations.begin(),constraint.equations.end(),equation)!=constraint.equations.end())
                for (const auto& member:constraint.equations) uses=uses || references(call,member);
        if (uses) result.push_back(&call);
    }
    return result;
}
std::vector<const SolvePlanNode*> equationOperations(
        const ExecutableEquationSystem& system,const CompiledSolvePlan& plan,std::string_view equation) {
    const auto consumers=equationConsumers(system,plan.compiledProgram,equation);
    std::vector<const SolvePlanNode*> result;
    const auto visit=[&](const auto& self,const SolvePlanNode& node)->void {
        if (!node.operation.empty()) {
            bool consumes=false;
            for (const auto* call:consumers) {
                consumes=consumes || node.occurrence==call->source.occurrence;
                for (const auto& math:node.equationCalls)
                    consumes=consumes || math.equation==call->source.equation
                        || std::any_of(call->calls.begin(),call->calls.end(),[&](const auto& ref){return ref.equation==math.equation;});
            }
            if (consumes) result.push_back(&node);
        }
        for (const auto& child:node.children) self(self,child);
    };
    visit(visit,plan.root);
    return result;
}

void validateMathematicalCoverage(const ExecutableEquationSystem& system,
        const CompiledSolvePlan& plan,bool runnable) {
    for (const auto& call:plan.compiledProgram.steps) {
        (void)system.registry.at(call.source.equation);
        for (const auto& ref:call.calls) (void)system.registry.at(ref.equation);
        for (const auto& id:call.sourceMathInputs) (void)system.registry.at(id);
    }
    for (const auto& constraint:system.constraints) {
        if (constraint.equations.empty() && runnable)
            throw std::runtime_error("Constraint has no declared mathematical realization: "+constraint.id);
        for (const auto& id:constraint.equations) (void)system.registry.at(id);
    }
    for (const auto& closure:system.closures) {
        if (runnable && closure.equations.empty() && closure.states.empty())
            throw std::runtime_error("Closure has no mathematical/STATE reference: "+closure.description);
        for (const auto& id:closure.equations) (void)system.registry.at(id);
        for (const auto& id:closure.states) (void)system.state.at(id);
    }
    if (runnable) for (const auto& equation:system.registry.entries())
        if (equationConsumers(system,plan.compiledProgram,equation.id).empty())
            throw std::runtime_error("Active equation has no compiled consumer: "+equation.id);
    if (runnable) for (const auto& equation:system.registry.entries())
        if (equationOperations(system,plan,equation.id).empty())
            throw std::runtime_error("Active equation has no compiled operation: "+equation.id);
    const auto visit=[&](const auto& self,const SolvePlanNode& node)->void {
        for (const auto& call:node.equationCalls) (void)system.registry.at(call.equation);
        for (const auto& child:node.children) self(self,child);
    };
    visit(visit,plan.root);
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
namespace {
std::string coefficientOwner(const CompiledEquationCall& call) {
    return call.temporalOwner.empty()?call.source.occurrence:call.temporalOwner;
}
std::vector<const CompiledEquationCall*> coefficientPublishers(
        const StateRegistry* state,const CompiledExecutionProgram& program) {
    std::vector<const CompiledEquationCall*> publishers;
    for(const auto& call:program.steps) {
        if(call.coefficientProvides.empty()) {
            if(!call.coefficientStates.empty() || !call.coefficientPrepare.empty() || !call.coefficientPublish.empty())
                throw std::runtime_error("Coefficient publication lacks resource: "+call.source.occurrence);
            continue;
        }
        if(call.coefficientProvides.size()!=1 || call.coefficientStates.empty())
            throw std::runtime_error("Coefficient provider requires one explicit resource and its complete state set: "+call.source.occurrence);
        std::set<std::string> symbols;
        for(const auto& symbol:call.coefficientStates) {
            if(state)state->at(symbol);
            if(!symbols.insert(symbol).second)throw std::runtime_error("Duplicate coefficient state: "+symbol);
        }
        if(call.coefficientVersion==StateVersion::Stage) {
            if(!call.temporalResidual || call.coefficientPrepare.empty() || call.coefficientPublish.empty())
                throw std::runtime_error("Stage coefficient provider lacks preparation/publication treatment: "+call.source.occurrence);
        } else if(call.coefficientVersion==StateVersion::Frozen) {
            if(!call.coefficientPrepare.empty() || !call.coefficientPublish.empty()
                || std::none_of(call.stateEffects.begin(),call.stateEffects.end(),[&](const auto& e){return e.symbol==call.coefficientProvides.front() && e.publish;}))
                throw std::runtime_error("Frozen coefficient provider lacks physical-step publication: "+call.source.occurrence);
            for(const auto& symbol:call.coefficientStates)
                if(std::none_of(call.stateEffects.begin(),call.stateEffects.end(),[&](const auto& e){return e.symbol==symbol && e.publish;}))
                    throw std::runtime_error("Frozen coefficient state is not published: "+symbol);
        } else throw std::runtime_error("Unsupported coefficient provider treatment: "+std::string(toString(call.coefficientVersion)));
        const auto previous=std::find_if(publishers.begin(),publishers.end(),[&](const auto* p) {
            return coefficientOwner(*p)==coefficientOwner(call) || p->coefficientProvides==call.coefficientProvides
                || std::any_of(call.coefficientStates.begin(),call.coefficientStates.end(),[&](const auto& id) {
                    return std::find(p->coefficientStates.begin(),p->coefficientStates.end(),id)!=p->coefficientStates.end();
                });
        });
        if(previous!=publishers.end()) {
            const auto& p=**previous;
            if(!call.fusionKey.empty() && call.fusionKey==p.fusionKey
                && call.backendProvider==p.backendProvider && call.backendOperation==p.backendOperation
                && call.coefficientProvides==p.coefficientProvides && call.coefficientStates==p.coefficientStates
                && call.coefficientVersion==p.coefficientVersion)continue;
            throw std::runtime_error("Conflicting coefficient publishers: "+p.source.occurrence+" and "+call.source.occurrence);
        }
        publishers.push_back(&call);
    }
    return publishers;
}
}
void bindCoefficientContracts(const StateRegistry& state,CompiledExecutionProgram& program) {
    const auto publishers=coefficientPublishers(&state,program);
    for(auto& call:program.steps) {
        std::set<std::string> symbols;
        for(auto& use:call.coefficientUses) {
            state.at(use.symbol);
            if(!symbols.insert(use.symbol).second || use.acceptedVersions.empty()
                || !use.producer.empty() || !use.resource.empty())
                throw std::runtime_error("Invalid/unresolved coefficient request: "+call.source.occurrence+" / "+use.symbol);
            const CompiledEquationCall* selected=nullptr;
            for(const auto* p:publishers) {
                if(std::find(p->coefficientStates.begin(),p->coefficientStates.end(),use.symbol)==p->coefficientStates.end())continue;
                if(std::find(use.acceptedVersions.begin(),use.acceptedVersions.end(),p->coefficientVersion)==use.acceptedVersions.end())
                    throw std::runtime_error("Incompatible coefficient treatment: consumer "+call.source.occurrence+" / "+use.symbol+" cannot read "+toString(p->coefficientVersion));
                if(selected)throw std::runtime_error("Ambiguous coefficient publisher: "+use.symbol);
                selected=p;
            }
            if(!selected)throw std::runtime_error("Missing coefficient publisher: consumer "+call.source.occurrence+" / "+use.symbol);
            use.producer=coefficientOwner(*selected);use.resource=selected->coefficientProvides.front();use.version=selected->coefficientVersion;
            const auto read=std::find_if(call.stateUses.begin(),call.stateUses.end(),[&](const auto& r){return r.symbol==use.symbol;});
            if(read==call.stateUses.end())call.stateUses.push_back({use.symbol,use.version,false,use.version==StateVersion::Frozen?use.resource:""});
            else if(read->version!=use.version || read->halo)
                throw std::runtime_error("Coefficient read/version/halo incompatible with publisher: "+use.symbol);
            if(use.version==StateVersion::Stage && std::find(call.stageRequires.begin(),call.stageRequires.end(),use.resource)==call.stageRequires.end())
                call.stageRequires.push_back(use.resource);
        }
    }
    validateCoefficientContracts(state,program);
}
void validateCoefficientContracts(const StateRegistry& state,const CompiledExecutionProgram& program) {
    coefficientPublishers(&state,program);
    validateCoefficientContracts(program);
}
void validateCoefficientContracts(const CompiledExecutionProgram& program) {
    const auto publishers=coefficientPublishers(nullptr,program);
    for(const auto& call:program.steps) {
        for(const auto& read:call.stateUses) {
            const bool coefficient=std::any_of(publishers.begin(),publishers.end(),[&](const auto* p){return std::find(p->coefficientStates.begin(),p->coefficientStates.end(),read.symbol)!=p->coefficientStates.end();});
            if(coefficient && std::none_of(call.coefficientUses.begin(),call.coefficientUses.end(),[&](const auto& use){return use.symbol==read.symbol;}))
                throw std::runtime_error("Coefficient read lacks frozen producer binding: "+call.source.occurrence+" / "+read.symbol);
        }
        std::set<std::string> symbols;
        for(const auto& use:call.coefficientUses) {
            if(!symbols.insert(use.symbol).second)throw std::runtime_error("Duplicate coefficient consumer binding: "+use.symbol);
            const auto p=std::find_if(publishers.begin(),publishers.end(),[&](const auto* item) {
                return coefficientOwner(*item)==use.producer && item->coefficientProvides.front()==use.resource;
            });
            if(p==publishers.end() || (*p)->coefficientVersion!=use.version
                || std::find((*p)->coefficientStates.begin(),(*p)->coefficientStates.end(),use.symbol)==(*p)->coefficientStates.end()
                || std::find(use.acceptedVersions.begin(),use.acceptedVersions.end(),use.version)==use.acceptedVersions.end())
                throw std::runtime_error("Missing/incompatible frozen coefficient binding: "+call.source.occurrence+" / "+use.symbol);
            const auto count=std::count_if(call.stateUses.begin(),call.stateUses.end(),[&](const auto& read) {
                return read.symbol==use.symbol && read.version==use.version && !read.halo
                    && read.snapshot==(use.version==StateVersion::Frozen?use.resource:"");
            });
            if(count!=1)throw std::runtime_error("Coefficient stateUses/version/halo does not match publication: "+use.symbol);
        }
    }
    if(!program.temporalParticipants.empty())for(const auto* publisher:publishers)if(publisher->coefficientVersion==StateVersion::Stage
        && std::count_if(program.temporalParticipants.begin(),program.temporalParticipants.end(),[&](const auto& p){return p.identity==coefficientOwner(*publisher);})!=1)
        throw std::runtime_error("Stage coefficient publisher lacks a unique temporal owner: "+coefficientOwner(*publisher));
    for(const auto& participant:program.temporalParticipants) {
        std::vector<CoefficientUse> expected;
        for(const auto& call:program.steps)if(call.temporalResidual && coefficientOwner(call)==participant.identity)
            for(const auto& use:call.coefficientUses)if(std::none_of(expected.begin(),expected.end(),[&](const auto& x){return x.symbol==use.symbol;}))expected.push_back(use);
        if(expected.size()!=participant.coefficientUses.size())throw std::runtime_error("Temporal coefficient consumer set differs from compiled calls: "+participant.identity);
        for(std::size_t i=0;i<expected.size();++i) {
            const auto& a=expected[i];const auto& b=participant.coefficientUses[i];
            if(a.symbol!=b.symbol || a.resource!=b.resource || a.producer!=b.producer || a.version!=b.version || a.acceptedVersions!=b.acceptedVersions)
                throw std::runtime_error("Temporal coefficient consumer binding differs from compiled call: "+participant.identity);
        }
        const auto p=std::find_if(publishers.begin(),publishers.end(),[&](const auto* item){return coefficientOwner(*item)==participant.identity;});
        if(p==publishers.end()) {
            if(!participant.coefficientProvides.empty() || !participant.coefficientStates.empty()
                || !participant.coefficientPrepare.empty() || !participant.coefficientPublish.empty())
                throw std::runtime_error("Temporal coefficient publisher has no compiled owner: "+participant.identity);
            continue;
        }
        if(p==publishers.end() || participant.coefficientVersion!=StateVersion::Stage
            || participant.coefficientStates!=(*p)->coefficientStates || participant.coefficientProvides!=(*p)->coefficientProvides
            || participant.coefficientPrepare!=(*p)->coefficientPrepare || participant.coefficientPublish!=(*p)->coefficientPublish)
            throw std::runtime_error("Temporal coefficient publication differs from compiled provider: "+participant.identity);
    }
}
void validateDataFlow(const StateRegistry& state,const CompiledExecutionProgram& program) {
    validateCoefficientContracts(state,program);
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
        if (x.readOnly() || y.readOnly()
            || x.storageKey.empty() || x.storageKey!=y.storageKey
            || x.storageBinding!=y.storageBinding
            || x.location!=y.location) return false;
        return x.componentOffset<y.componentOffset+y.components && y.componentOffset<x.componentOffset+x.components;
    };
    std::map<std::string,Revision> revisions;
    std::map<std::string,std::map<std::string,int>> published;
    std::map<std::string,std::string> invalidLazy;
    std::set<std::string> snapshots;
    std::set<std::string> persistent,invalidated;
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
        if (node.kind==ExecutionKind::Commit) {
            for (const auto& call:program.steps) for (const auto& resource:call.commitRequires)
                if (!snapshots.count(resource)) throw std::runtime_error("Commit requires unpublished resource: "+resource);
        }
        if (node.kind==ExecutionKind::EquationCall) {
            const auto at=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& call){return call.source.occurrence==node.step.occurrence;});
            if (at==program.steps.end()) throw std::runtime_error("Missing compiled occurrence: "+node.step.occurrence);
            // Complete fused groups are one algebraic operation, not sequential partial writes.
            if (!at->fusionKey.empty() && at!=program.steps.begin() && (at-1)->fusionKey==at->fusionKey) return;
            std::vector<StateUse> uses;
            std::vector<StateEffect> effects;
            std::vector<CoefficientUse> coefficientUses;
            const auto members=at->fusionKey.empty()?1:at->fusionMembers.size();
            for (std::size_t i=0;i<members;++i) {
                const auto& member=*(at+i);
                uses.insert(uses.end(),member.stateUses.begin(),member.stateUses.end());
                coefficientUses.insert(coefficientUses.end(),member.coefficientUses.begin(),member.coefficientUses.end());
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
                            return std::find(p.targets.begin(),p.targets.end(),use.symbol)!=p.targets.end() || std::find(p.coefficientStates.begin(),p.coefficientStates.end(),use.symbol)!=p.coefficientStates.end();
                        });
                    if (use.version==StateVersion::Stage && !stage && !commonStage)
                        throw std::runtime_error("Unsupported stage STATE use outside declared StageLoop: "+node.step.occurrence);
                    // Common RHS calls read one prepared epoch. A coefficient
                    // is explicitly published before all RHS leaves, not between
                    // source HOW members whose writes commit after the barrier.
                    const auto coefficient=std::find_if(coefficientUses.begin(),coefficientUses.end(),[&](const auto& item){return item.symbol==use.symbol && item.version==StateVersion::Stage;});
                    if(use.version==StateVersion::Stage && commonStage && coefficient!=coefficientUses.end()) {
                        if(use.halo || std::none_of(program.temporalParticipants.begin(),program.temporalParticipants.end(),[&](const auto& p){
                            return p.identity==coefficient->producer && !p.coefficientPrepare.empty()
                                && std::find(p.coefficientProvides.begin(),p.coefficientProvides.end(),coefficient->resource)!=p.coefficientProvides.end()
                                && std::find(p.coefficientStates.begin(),p.coefficientStates.end(),use.symbol)!=p.coefficientStates.end();
                        }))throw std::runtime_error("Missing bound stage coefficient publication: "+use.symbol);
                    } else check(check,use.symbol,node.step.occurrence,use.halo);
                }
            }
            for (std::size_t i=0;i<members;++i) for (const auto& resource:(at+i)->workspaceInvalidates) {
                snapshots.erase(resource);persistent.erase(resource);invalidated.insert(resource);
                // An explicit replacement ends this frozen resource's lifetime.
                // Ordinary writes without invalidation still fail frozen checks.
                for (auto binding=frozen.begin();binding!=frozen.end();) {
                    const auto suffix=":"+resource;
                    if (binding->first.size()>=suffix.size() && binding->first.compare(binding->first.size()-suffix.size(),suffix.size(),suffix)==0)
                        binding=frozen.erase(binding);
                    else ++binding;
                }
            }
            for (const auto& effect:effects) {
                if (!state.contains(effect.symbol)) {
                    if (!effect.publish) throw std::runtime_error("Unknown STATE write: "+effect.symbol);
                    snapshots.insert(effect.symbol);invalidated.erase(effect.symbol);revisions[effect.symbol].value=++epoch;continue;
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
                    snapshots.insert(workspace);invalidated.erase(workspace);
                    if (std::none_of(effects.begin(),effects.end(),[&](const auto& effect){return effect.symbol==workspace;}))
                        revisions[workspace].value=++epoch;
                }
                for (const auto& workspace:member.stepWorkspaces)
                    if (snapshots.count(workspace)) persistent.insert(workspace);
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
        if (node.kind==ExecutionKind::Loop || node.kind==ExecutionKind::StageLoop) {snapshots=savedSnapshots;snapshots.insert(persistent.begin(),persistent.end());for (const auto& resource:invalidated) snapshots.erase(resource);frozen=savedFrozen;}
    };
    for (int physicalStep=0;physicalStep<2;++physicalStep) {
        snapshots.clear();persistent.clear();frozen.clear();invalidated.clear();
        // Lagged versions are immutable step-entry inputs, including legally zero initialized force.
        for (const auto& call:program.steps) for (const auto& use:call.stateUses)
            if (use.version==StateVersion::Lagged && state.contains(use.symbol)
                && state.at(use.symbol).availableVersion==StateVersion::Lagged && use.snapshot.empty())
                snapshots.insert(use.snapshot.empty()?use.symbol:use.snapshot);
        execute(execute,program.root,false);
    }
}
}
