/// @file SF_lowering.cpp
/// @brief 四轴引用解析、所选 provider 编译与控制流冻结；不访问运行时场数值。
#include "SF_lowering.h"
#include "SF_stateRealization.h"
#include "SF_executionContract.h"
#include "SF_formulaCompiler.h"
#include "SF_numericalCompiler.h"
#include "SF_solvePlan.h"
#include "SF_providerResolver.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <map>
#include <cmath>
#include <set>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <utility>
namespace SF::System::Lowering {
namespace {
SolvePlanNode lowerTree(const CompiledExecutionProgram&);

void expandRecipeScope(ExecutionScope& root,const std::vector<CompiledEquationCall>& calls) {
    const auto temporal=[&](const auto& self,const ExecutionScope& node)->bool {
        if(node.kind==ExecutionKind::StageLoop)throw std::runtime_error("Recipe-owned HOW cannot contain another authored StageLoop");
        if(node.kind==ExecutionKind::EquationCall || node.kind==ExecutionKind::BlockSolve) {
            const auto occurrence=node.kind==ExecutionKind::BlockSolve?node.id:node.step.occurrence;
            const auto call=std::find_if(calls.begin(),calls.end(),[&](const auto& c){return c.source.occurrence==occurrence;});
            return call!=calls.end() && call->temporalResidual;
        }
        return std::any_of(node.children.begin(),node.children.end(),[&](const auto& child){return self(self,child);});
    };
    const bool transient=temporal(temporal,root);
    bool stageActions=false;
    const auto classify=[&](const auto& self,ExecutionScope& node,ExecutionFrequency inherited)->ExecutionFrequency {
        auto frequency=node.frequency==ExecutionFrequency::Inherit?inherited:node.frequency;
        if(inherited==ExecutionFrequency::Stage && frequency==ExecutionFrequency::Step)
            throw std::runtime_error("HOW Step frequency cannot escape its Stage scope: "+node.id);
        if(node.kind==ExecutionKind::EquationCall || node.kind==ExecutionKind::BlockSolve) {
            const bool timed=temporal(temporal,node);
            if(timed) {
                if(node.kind!=ExecutionKind::EquationCall || frequency==ExecutionFrequency::Step)
                    throw std::runtime_error("Unsupported temporal HOW scope/frequency: "+node.id);
                frequency=ExecutionFrequency::Stage;
            } else if(transient && frequency==ExecutionFrequency::Inherit)
                throw std::runtime_error("Ambiguous HOW frequency at "+node.id+"; declare Step or Stage explicitly (no adjacent-stage inference)");
            if(!timed && frequency==ExecutionFrequency::Stage)stageActions=true;
        } else if(node.kind==ExecutionKind::Commit) {
            if(frequency==ExecutionFrequency::Stage)throw std::runtime_error("Commit cannot be repeated by a Stage frequency scope");
            return ExecutionFrequency::Step;
        } else {
            std::optional<ExecutionFrequency> children;
            for(auto& child:node.children) {
                const auto f=self(self,child,frequency);
                if(children && *children!=f) {
                    if(&node!=&root)throw std::runtime_error("HOW group crosses Step/Stage boundary: "+node.id);
                } else if(!children)children=f;
            }
            if(&node!=&root && children)frequency=*children;
            if(node.kind==ExecutionKind::Loop && temporal(temporal,node))
                throw std::runtime_error("Unsupported temporal HOW scope: temporal RHS cannot repeat in an ordinary Loop: "+node.id);
        }
        if(frequency==ExecutionFrequency::Inherit)frequency=ExecutionFrequency::Step;
        node.frequency=frequency;return frequency;
    };
    classify(classify,root,ExecutionFrequency::Inherit);
    std::size_t first=root.children.size(),last=0;
    for(std::size_t i=0;i<root.children.size();++i)if(root.children[i].frequency==ExecutionFrequency::Stage){first=std::min(first,i);last=i;}
    if(first==root.children.size() || !stageActions)return;
    if(!transient)throw std::runtime_error("Stage frequency requires an implemented temporal group");
    for(std::size_t i=first;i<=last;++i)if(root.children[i].frequency!=ExecutionFrequency::Stage)
        throw std::runtime_error("Non-contiguous Stage scope at "+root.children[i].id+"; source Step operations are never moved across stages");
    ExecutionScope stages;stages.kind=ExecutionKind::StageLoop;stages.id="recipe.stages";
    stages.origin={OriginKind::Generated,"HOW Stage frequency expanded by the selected time recipe"};
    stages.children.assign(root.children.begin()+first,root.children.begin()+last+1);
    root.children.erase(root.children.begin()+first,root.children.begin()+last+1);
    root.children.insert(root.children.begin()+first,std::move(stages));
}
void validateClockCommit(const SolvePlanNode& root) {
    int commits=0;
    const auto visit=[&](const auto& self,const SolvePlanNode& node,bool repeated,bool stage)->void {
        const bool container=node.kind==PlanNodeKind::Sequence || node.kind==PlanNodeKind::Loop || node.kind==PlanNodeKind::StageLoop;
        repeated=repeated || (container && node.repetitions!=1);
        stage=stage || node.kind==PlanNodeKind::StageLoop;
        if(node.operation==OpIds::TimeCommit)
            if(++commits>1 || repeated || stage)
                throw std::runtime_error("Physical clock Commit must occur once outside repeated or stage scopes.");
        for(const auto& child:node.children)self(self,child,repeated,stage);
    };
    visit(visit,root,false,false);
}
const ExecutionScope* stageScope(const CompiledExecutionProgram& program) {
    const ExecutionScope* found=nullptr;
    for(const auto& node:program.root.children)if(node.kind==ExecutionKind::StageLoop) {
        if(found)throw std::runtime_error("One implemented time group cannot own multiple Stage scopes");
        found=&node;
    }
    return found;
}
std::pair<std::vector<SolvePlanNode>,std::vector<SolvePlanNode>> stepPlacement(const CompiledExecutionProgram& program) {
    if(program.root.kind!=ExecutionKind::Sequence || program.root.children.empty()
        || program.root.children.back().kind!=ExecutionKind::Commit)
        throw std::runtime_error("Temporal placement requires ordered Step/Stage groups and one terminal Commit");
    const auto count=[](const auto& self,const ExecutionScope& node)->std::size_t {
        if(node.kind==ExecutionKind::EquationCall || node.kind==ExecutionKind::BlockSolve)return 1;
        std::size_t n=0;for(const auto& child:node.children)n+=self(self,child);return n;
    };
    std::vector<bool> stage;std::size_t next=0;
    for(const auto& node:program.root.children) {
        if(node.kind==ExecutionKind::Commit)continue;
        const auto n=count(count,node);if(!n || next+n>program.steps.size())throw std::runtime_error("Unimplemented temporal placement scope: "+node.id);
        const bool timed=program.steps[next].temporalResidual || program.steps[next].stageLocal;
        if(node.kind!=ExecutionKind::StageLoop)for(std::size_t i=0;i<n;++i)
            if((program.steps[next+i].temporalResidual || program.steps[next+i].stageLocal)!=timed)
                throw std::runtime_error("A source group crosses Step/Stage placement: "+node.id);
        const bool continuation=node.kind==ExecutionKind::EquationCall && next>0 && !program.steps[next].fusionKey.empty()
            && program.steps[next-1].fusionKey==program.steps[next].fusionKey;
        if(!continuation)stage.push_back(node.kind==ExecutionKind::StageLoop || timed);
        next+=n;
    }
    if(next!=program.steps.size() || program.loweredRoot.children.size()!=stage.size()+1)
        throw std::runtime_error("Temporal placement has unsupported root lifecycle decorations");
    const auto first=std::find(stage.begin(),stage.end(),true);
    const auto last=std::find(stage.rbegin(),stage.rend(),true);
    if(first==stage.end())throw std::runtime_error("Missing temporal placement group");
    const auto begin=static_cast<std::size_t>(first-stage.begin());
    const auto end=stage.size()-static_cast<std::size_t>(last-stage.rbegin());
    if(std::find(stage.begin()+begin,stage.begin()+end,false)!=stage.begin()+end)
        throw std::runtime_error("Unsupported mixed temporal topology: temporal calls must be contiguous; no Step call is moved");
    return {{program.loweredRoot.children.begin(),program.loweredRoot.children.begin()+begin},
            {program.loweredRoot.children.begin()+end,program.loweredRoot.children.end()-1}};
}
std::vector<SolvePlanNode> stageEvaluations(const CompiledExecutionProgram& program) {
    const auto at=std::find_if(program.loweredRoot.children.begin(),program.loweredRoot.children.end(),[](const auto& node){return node.kind==PlanNodeKind::StageLoop;});
    if(at==program.loweredRoot.children.end())throw std::runtime_error("Missing compiled Stage scope");
    std::set<std::string> evaluated;
    const auto callFor=[&](const SolvePlanNode& node)->const CompiledEquationCall& {
        const auto call=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& c){return c.source.occurrence==node.occurrence;});
        if(call==program.steps.end() || !node.children.empty() || node.operation!=call->backendOperation
            || (!call->stageLocal && !call->temporalResidual))
            throw std::runtime_error("Unsupported numerical leaf in authored Stage scope: "+node.id);
        return *call;
    };
    const auto lower=[&](const auto& self,const SolvePlanNode& scope)->SolvePlanNode {
        auto output=scope;output.children.clear();
        for(std::size_t i=0;i<scope.children.size();) {
            const auto& child=scope.children[i];
            if(child.operation.empty()) {
                if(child.kind!=PlanNodeKind::Sequence && child.kind!=PlanNodeKind::Loop)
                    throw std::runtime_error("Unsupported Stage evaluation scope: "+child.id);
                output.children.push_back(self(self,child));++i;continue;
            }
            const auto& call=callFor(child);auto leaf=child;
            if(call.stageLocal) {
                leaf.operation=call.backendOperation;++i;
            } else {
                const auto identity=call.temporalOwner.empty()?call.source.occurrence:call.temporalOwner;
                const auto participant=std::find_if(program.temporalParticipants.begin(),program.temporalParticipants.end(),[&](const auto& p){return p.identity==identity;});
                std::vector<const CompiledEquationCall*> members;
                for(const auto& step:program.steps)if(step.temporalResidual && (step.temporalOwner.empty()?step.source.occurrence:step.temporalOwner)==identity)members.push_back(&step);
                if(participant==program.temporalParticipants.end() || members.empty()
                    || !evaluated.insert(identity).second || i+members.size()>scope.children.size())
                    throw std::runtime_error("Stage numerical owner requires one complete adjacent HOW group: "+identity);
                for(std::size_t j=0;j<members.size();++j) {
                    const auto& member=callFor(scope.children[i+j]);
                    if(&member!=members[j] || member.temporalRhs!=participant->rhs)
                        throw std::runtime_error("Stage numerical owner crosses a HOW operation/scope or changes member order: "+identity);
                }
                // The original fused kernel consumes all adjacent mathematical members once.
                // An intervening algebraic call cannot be moved into or past this group.
                leaf.operation=participant->rhs;leaf.provider=participant->provider;
                leaf.occurrence=participant->identity;leaf.equationCalls=participant->calls;
                if(participant->calls.size()>1)leaf.target.reset();
                i+=members.size();
            }
            leaf.id=leaf.name=leaf.operation;output.children.push_back(std::move(leaf));
        }
        return output;
    };
    auto output=lower(lower,*at);
    if(evaluated.size()!=program.temporalParticipants.size())
        throw std::runtime_error("Stage scope does not evaluate every temporal numerical owner exactly once");
    return output.children;
}

void validateCoefficientInputs(const StateRegistry& state,const CompiledExecutionProgram& program) {
    std::vector<std::string> prepared;
    for(const auto& participant:program.temporalParticipants)
        prepared.insert(prepared.end(),participant.targets.begin(),participant.targets.end());
    for(const auto& participant:program.temporalParticipants)if(!participant.coefficientPrepare.empty()) {
        const auto ready=[&](const auto& self,const std::string& id,std::vector<std::string> path)->bool {
            if(std::find(path.begin(),path.end(),id)!=path.end())throw std::runtime_error("Invalid coefficient Stage dependency cycle: "+id);
            path.push_back(id);
            const auto& symbol=state.at(id);
            if(std::find(prepared.begin(),prepared.end(),id)!=prepared.end())return true;
            if(!symbol.producer || symbol.evaluation!=StateEvaluation::Lazy
                || std::find(symbol.producer->versions.begin(),symbol.producer->versions.end(),StateViewKind::Stage)==symbol.producer->versions.end())return false;
            return std::all_of(symbol.dependencies.begin(),symbol.dependencies.end(),[&](const auto& dependency){return self(self,dependency,path);});
        };
        for(const auto& input:participant.coefficientInputs) {
            if(input.version==StateVersion::Stage && !ready(ready,input.symbol,{}))
                throw std::runtime_error("Coefficient preparation precedes its required Stage publication: "+participant.identity+" / "+input.symbol+"; HOW is not reordered");
        }
        prepared.insert(prepared.end(),participant.coefficientStates.begin(),participant.coefficientStates.end());
    }
}
SolvePlanNode lowerTree(const CompiledExecutionProgram& program) {
    if (!program.loweredRoot.children.empty()) return program.loweredRoot;
    std::size_t next=0;
    std::set<std::string> iterations;
    std::map<std::string,const ExecutionScope*> iterationScopes;
    const auto scopeOf=[&](const auto& self,const ExecutionScope& node,const ExecutionScope* nearest,const ExecutionScope* signal)->void {
        if(node.kind==ExecutionKind::Loop){nearest=&node;if(!node.terminationSignal.empty())signal=&node;}
        if(node.kind==ExecutionKind::EquationCall || node.kind==ExecutionKind::BlockSolve) {
            const auto id=node.kind==ExecutionKind::BlockSolve?node.id:node.step.occurrence;
            const auto call=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& c){return c.source.occurrence==id;});
            if(call!=program.steps.end())iterationScopes[id]=call->iterationSignalRequired?signal:nearest;
        }
        for(const auto& child:node.children)self(self,child,nearest,signal);
    };
    scopeOf(scopeOf,program.root,nullptr,nullptr);
    const auto lower=[&](const auto& self,const ExecutionScope& source)->SolvePlanNode {
        SolvePlanNode result;
        switch (source.kind) {
            case ExecutionKind::Sequence: result.kind=PlanNodeKind::Sequence;break;
            case ExecutionKind::Loop: result.kind=PlanNodeKind::Loop;break;
            case ExecutionKind::StageLoop: result.kind=PlanNodeKind::StageLoop;break;
            case ExecutionKind::EquationCall: result.kind=PlanNodeKind::Update;break;
            case ExecutionKind::Commit: result.kind=PlanNodeKind::Sequence;break;
            case ExecutionKind::BlockSolve: result.kind=PlanNodeKind::BlockSolve;break;
        }
        result.repetitions=source.repetitions;
        result.terminationSignal=source.terminationSignal;
        result.minimumIterations=source.minimumIterations;
        result.residualAll=source.residualAll;result.failOnMaximum=source.failOnMaximum;
        for(const auto& id:source.residuals) {
            const auto at=std::find_if(program.residuals.begin(),program.residuals.end(),[&](const auto& r){return r.id==id;});
            if(at==program.residuals.end())throw std::runtime_error("Missing compiled residual: "+id);
            result.residuals.push_back(*at);
        }
        result.id=!source.id.empty() ? source.id
            : source.kind==ExecutionKind::EquationCall ? source.step.equation
            : std::string(toString(source.kind));
        result.name=result.id;
        if(source.kind==ExecutionKind::Commit && !source.commitSources.empty()) {
            result.kind=PlanNodeKind::Commit;result.operation=instanceOperation("state.commit",source.id);
            result.provider="execution.state";return result;
        }
        if (source.kind==ExecutionKind::EquationCall || source.kind==ExecutionKind::BlockSolve) {
            if (next>=program.steps.size())
                throw std::runtime_error("Compiled HOW has fewer method steps than source Steps.");
            const auto index=next++;
            const auto& method=program.steps[index];
            if (method.source.equation!=source.step.equation)
                throw std::runtime_error("Compiled HOW step disagrees with source HOW.");
            if (!method.fusionKey.empty() && index>0
                && program.steps[index-1].fusionKey==method.fusionKey) {
                result.kind=PlanNodeKind::Sequence;result.id.clear();
                return result;
            }
            if (!method.fragment.children.empty() || !method.fragment.operation.empty()) {
                result=method.fragment;
            } else {
                if (method.backendOperation.empty())
                    throw std::runtime_error("Compiled HOW step has no method fragment/backend.");
                result.operation=method.backendOperation;
                result.provider=method.backendProvider;
                result.equationCalls=method.calls;
            }
            if (!method.fusionKey.empty()) {
                const auto attach=[&](const auto& self,SolvePlanNode& node)->void {
                    if (node.children.empty()) for (std::size_t j=1;j<method.fusionMembers.size();++j) {
                        const auto& member=program.steps.at(index+j);
                        node.equationCalls.insert(node.equationCalls.end(),member.calls.begin(),member.calls.end());
                    }
                    for (auto& child:node.children) self(self,child);
                };
                attach(attach,result);
            }
            const auto bindTarget=[&](const auto& self,SolvePlanNode& node)->void {
                if (node.children.empty()) {
                    node.target=method.target;
                    node.occurrence=method.source.occurrence;
                }
                for (auto& child:node.children) self(self,child);
            };
            bindTarget(bindTarget,result);
        }
        for (const auto& child:source.children) {
            auto compiled=self(self,child);
            if (compiled.id.empty() && compiled.children.empty() && compiled.operation.empty()) continue;
            result.children.push_back(std::move(compiled));
        }
        const auto append=[&](std::vector<SolvePlanNode>& nodes,const std::vector<BoundOperation>& actions,bool commit=false) {
            for(const auto& action:actions) {
                if(action.operation.empty() || action.provider.empty() || action.operation==OpIds::TimeCommit)
                    throw std::runtime_error("Invalid numerical input/output action at HOW scope "+source.id);
                const auto old=std::find_if(nodes.begin(),nodes.end(),[&](const auto& n){return n.operation==action.operation;});
                if(old!=nodes.end()) {
                    if(old->provider!=action.provider)throw std::runtime_error("Conflicting numerical action owners: "+action.operation);
                    continue;
                }
                SolvePlanNode leaf;leaf.kind=commit?PlanNodeKind::Commit:PlanNodeKind::Update;
                leaf.id=leaf.name=leaf.operation=action.operation;leaf.provider=action.provider;
                nodes.push_back(std::move(leaf));
            }
        };
        if(&source==&program.root) {
            std::vector<SolvePlanNode> inputs;
            for(const auto& call:program.steps)append(inputs,call.stepInputs);
            if(!inputs.empty() && source.kind!=ExecutionKind::Sequence)
                throw std::runtime_error("Numerical step input preparation requires an ordered physical-step scope.");
            result.children.insert(result.children.begin(),inputs.begin(),inputs.end());
        }
        if(source.kind==ExecutionKind::Loop) {
            std::vector<SolvePlanNode> inputs,outputs;
            for(const auto& call:program.steps)if(iterationScopes[call.source.occurrence]==&source) {
                append(inputs,call.iterationInputs);append(outputs,call.iterationOutputs);
                iterations.insert(call.source.occurrence);
            }
            result.children.insert(result.children.begin(),inputs.begin(),inputs.end());
            result.children.insert(result.children.end(),outputs.begin(),outputs.end());
        }
        if(source.kind==ExecutionKind::Sequence && &source!=&program.root) {
            std::vector<SolvePlanNode> inputs;
            for(const auto& child:source.children)if(child.kind==ExecutionKind::EquationCall) {
                const auto call=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& c){return c.source.occurrence==child.step.occurrence;});
                if(call!=program.steps.end())append(inputs,call->groupInputs);
            }
            result.children.insert(result.children.begin(),inputs.begin(),inputs.end());
        }
        if(source.kind==ExecutionKind::Commit && source.commitSources.empty()) {
            const bool synchronous=stageScope(program) || std::any_of(program.steps.begin(),program.steps.end(),[](const auto& c){return c.temporalResidual && c.synchronousStages;});
            std::string clockOwner;
            for(std::size_t i=0;i<program.steps.size();++i) {
                const auto& call=program.steps[i];
                BoundOperation publication=call.physicalPublication;
                if(synchronous && call.temporalResidual)publication={call.temporalPublish,call.backendProvider};
                if(publication.operation.empty())continue;
                if(i>=next)throw std::runtime_error("HOW Commit precedes its required result: "+call.source.occurrence);
                append(result.children,{publication},true);
                if(!synchronous && !clockOwner.empty() && clockOwner!=publication.provider)
                    throw std::runtime_error("Physical Commit has incompatible publication/clock owners: "+source.id);
                clockOwner=publication.provider;
            }
            if(!result.children.empty()) {
                SolvePlanNode clock;clock.kind=PlanNodeKind::Commit;clock.id=clock.name=clock.operation=OpIds::TimeCommit;
                clock.provider=synchronous?TemporalOps::Provider:clockOwner;result.children.push_back(std::move(clock));
            }
        }
        if (source.kind==ExecutionKind::Commit && result.children.empty())
            throw std::runtime_error("Unsupported: HOW Commit has no declared physical result publication at scope "+source.id);
        return result;
    };
    auto result=lower(lower,program.root);
    for(const auto& call:program.steps)
        if((!call.iterationInputs.empty() || !call.iterationOutputs.empty()) && !iterations.count(call.source.occurrence))
            throw std::runtime_error("Numerical iteration snapshots/publication require an explicit HOW Loop: "+call.source.occurrence);
    if (next!=program.steps.size())
        throw std::runtime_error("Compiled HOW contains unconsumed method steps.");
    return result;
}

}
CompiledExecutionProgram compileProgram(
        const ExecutableEquationSystem& system,const StateRegistry& state,
        const ExecutionProgram& program,
        const std::vector<NumericalBinding>& bindings,
        const ProviderRegistry& registry,
        const CompiledTimeRecipe* time,
        const ITemporalMethod* temporalMethod) {
    CompiledExecutionProgram result;
    ExecutionScope root=program.root;
    if (!program.explicitOrder) {
        orderExecution(root);
        composeDefaultPlacement(program,root);
    }
    const auto address=[&](const auto& self,ExecutionScope& scope)->void {
        for (std::size_t i=0;i<scope.children.size();++i) {
            auto& child=scope.children[i];
            if (child.id.empty()) child.id=scope.id+"/"+std::to_string(i);
            if (child.kind==ExecutionKind::EquationCall && child.step.occurrence.empty()) child.step.occurrence=child.id;
            self(self,child);
        }
    };
    address(address,root);
    validatePlacement(program,root);
    const auto explicitFrequency=[&](const auto& self,const ExecutionScope& scope)->bool {
        return scope.frequency!=ExecutionFrequency::Inherit
            || std::any_of(scope.children.begin(),scope.children.end(),[&](const auto& child){return self(self,child);});
    };
    const bool recipeScope=program.recipeScopes || explicitFrequency(explicitFrequency,root);
    result.root=root;
    result.residuals=program.residuals;
    std::vector<EquationCall*> steps;
    std::vector<const ExecutionScope*> parents;
    std::map<const EquationCall*,const ExecutionScope*> blocks;
    const auto collect=[&](const auto& self,ExecutionScope& node,const ExecutionScope* parent,bool stage=false)->void {
        stage=stage || node.kind==ExecutionKind::StageLoop || node.frequency==ExecutionFrequency::Stage;
        if (node.kind==ExecutionKind::EquationCall) {
            if(!node.step.inputs.empty())throw std::runtime_error("HOW cannot declare input versions; bind the occurrence in STATE.");
            if (!node.children.empty())
                throw std::runtime_error("HOW Step cannot contain child nodes.");
            node.step.stageScoped=stage;
            steps.push_back(&node.step);
            parents.push_back(parent);
            return;
        }
        if (node.kind==ExecutionKind::Commit) {
            if(!node.children.empty())throw std::runtime_error("HOW Commit cannot contain child nodes.");
            if(node.commitSources.size()!=node.commitTargets.size())throw std::runtime_error("Variable Commit requires paired sources and targets.");
            std::set<std::string> targets;
            for(std::size_t i=0;i<node.commitSources.size();++i) {
                const auto& from=node.commitSources[i];const auto& to=node.commitTargets[i];
                if((from.kind!=TargetKind::Working && from.kind!=TargetKind::Correction) || to.kind!=TargetKind::Physical)
                    throw std::runtime_error("Variable Commit requires Working/Correction source and Physical target: "+node.id);
                const auto& a=state.at(from.symbol);const auto& b=state.at(to.symbol);
                if(a.storageBinding!=StorageBinding::NamedDistributed || b.storageBinding!=StorageBinding::NamedDistributed || a.componentOffset || b.componentOffset)
                    throw std::runtime_error("Unsupported variable Commit storage: named unsliced fields required");
                if(b.readOnly() || b.constantValue || !targets.insert(to.symbol).second)
                    throw std::runtime_error("Variable Commit target is read-only or duplicated: "+to.symbol);
                Math::requireCompatible(a.signature,b.signature,"Variable Commit");
                if(a.components!=b.components || a.location!=b.location || a.storageBinding!=b.storageBinding)
                    throw std::runtime_error("Variable Commit has incompatible STATE layout: "+node.id);
            }
            if(!node.commitSources.empty())result.operations.push_back({instanceOperation("state.commit",node.id),
                "explicit variable COPY; no physical time advance",OperationStage::StepCommit,
                {OperationCapability::StatePublication},{OriginKind::Generated,"HOW variable Commit"}});
            if(stage && !node.commitSources.empty())throw std::runtime_error("Variable Commit inside explicit StageLoop is not implemented.");
        }
        if (node.kind==ExecutionKind::Loop
            && (node.repetitions<1 || node.children.empty()
                || node.minimumIterations<1 || node.minimumIterations>node.repetitions))
            throw std::runtime_error("HOW Loop requires a body and valid iteration limits.");
        if(node.kind==ExecutionKind::BlockSolve) {
            if(node.equations.empty() || node.unknowns.empty())throw std::runtime_error("Implicit block requires equations and unknowns");
            std::vector<std::string> seen;
            for(const auto& id:node.equations){system.registry.at(id);if(std::find(seen.begin(),seen.end(),id)!=seen.end())throw std::runtime_error("Duplicate block equation: "+id);seen.push_back(id);}
            seen.clear();
            for(const auto& id:node.unknowns){const auto& s=state.at(id);if(s.readOnly() || s.constantValue)throw std::runtime_error("Block target is not writable: "+id);if(std::find(seen.begin(),seen.end(),id)!=seen.end())throw std::runtime_error("Duplicate block unknown: "+id);seen.push_back(id);bool present=false;for(const auto& eq:node.equations){const auto& e=system.registry.at(eq);for(const auto* side:{&e.lhs,&e.rhs}){const auto shape=inspectFormula(*side,id);present=present || std::find(shape.symbols.begin(),shape.symbols.end(),id)!=shape.symbols.end();}}if(!present)throw std::runtime_error("Block target absent from residual system: "+id);}
            if(!node.children.empty())throw std::runtime_error("BlockSolve cannot contain sequential children");
            node.step={node.equations.front(),{node.unknowns.front()},node.id};node.step.stageScoped=stage;
            blocks.emplace(&node.step,&node);steps.push_back(&node.step);parents.push_back(parent);return;
        }
        for (auto& child:node.children) self(self,child,&node,stage);
    };
    collect(collect,result.root,nullptr);
    std::vector<std::string> occurrences;
    for (const auto* step:steps) {
        (void)system.registry.at(step->equation);
        if(step->target.kind!=TargetKind::Workspace)(void)state.at(step->target.symbol);
        if (std::find(occurrences.begin(),occurrences.end(),step->occurrence)!=occurrences.end())
            throw std::runtime_error("Duplicate execution occurrence: "+step->occurrence);
        occurrences.push_back(step->occurrence);
    }
    resolveCallInputs(state,steps);
    for (const auto& binding:bindings) {
        if (binding.equation.empty() || (binding.method.empty() && binding.occurrence.empty()))
            throw std::runtime_error("Numerical binding requires equation and provider.");
        if (std::none_of(steps.begin(),steps.end(),
            [&](const EquationCall* call) {
                return (binding.equation=="*" || call->equation==binding.equation)
                    && (binding.occurrence.empty() || call->occurrence==binding.occurrence);
            }))
            throw std::runtime_error("Numerical binding addresses no execution occurrence: "+binding.equation);
    }
    for (const auto* step:steps) {
        auto effective=resolveNumerics(*step,bindings,blocks.count(step),&system);
        resolveOperatorReferences(system,*step,effective.value);
        const NumericalBinding* found=&effective.value;
        const auto& mathematical=system.registry.at(step->equation);
        if (step->target.kind!=TargetKind::Workspace) state.at(step->target.symbol);
        const IProvider* legacy=nullptr;
        CompiledEquationCall compiled;
        try {
            state.validateSymbols(mathematical.lhs);state.validateSymbols(mathematical.rhs);
            auto direct=blocks.count(step)?std::optional<CompiledEquationCall>{}:compileDirectEquation(system,*step,*found);
            if(direct)compiled=std::move(*direct);
            else {legacy=&registry.at(found->method);compiled=blocks.count(step)?legacy->compileBlock(system,*blocks.at(step),*found):legacy->compile(system,*step,*found);}
        }
        catch(const std::exception& error) {
            throw std::runtime_error(system.registry.at(step->equation).origin.source+": occurrence="+step->occurrence+" equation="+step->equation+" target="+step->target.symbol+" provider="+found->method+": "+error.what());
        }
        for(const auto& read:step->inputs)
            if(std::none_of(compiled.stateUses.begin(),compiled.stateUses.end(),[&](const auto& use){return use.symbol==read.symbol && use.version==read.version
                    && use.halo==read.halo && use.snapshot==read.snapshot && use.perIteration==read.perIteration && use.kind==read.kind;}))
                throw std::runtime_error("Unconsumed STATE input version: "+read.symbol+" at "+step->occurrence);
        if(!(legacy?legacy->usesSpatialRecipes():compiled.spatialTerms))for(const auto& op:found->operators)
            if(std::find(compiled.operatorBindings.begin(),compiled.operatorBindings.end(),op.provider)==compiled.operatorBindings.end())
                throw std::runtime_error("Unconsumed numerical operator choice: occurrence="+step->occurrence+" equation="+step->equation+" operator="+op.occurrence+" provider="+op.provider);
        compiled.numericalOperators=effective.value.operators;
        compiled.numericalOrigins=effective.origins;
        if(effective.value.time && (!compiled.temporalResidual || !time || *effective.value.time!=time->id()))
            throw std::runtime_error("Time recipe conflict/unused ddt choice: occurrence="+step->occurrence+" equation="+step->equation+" provider="+found->method+"; one resolved time recipe is required");
        compiled.legacyProvider=legacy!=nullptr;
        if(legacy)compiled.backendProvider=legacy->runtimeProvider();
        const auto ownFragment=[&](const auto& self,SolvePlanNode& node)->void {
            if (!node.operation.empty())
                if(legacy)node.provider=legacy->operationProvider(node.operation);
            for (auto& child:node.children) self(self,child);
        };
        ownFragment(ownFragment,compiled.fragment);
        if (compiled.source.equation!=step->equation || compiled.source.occurrence!=step->occurrence
            || compiled.source.target.symbol!=step->target.symbol || compiled.source.target.kind!=step->target.kind
            || compiled.target.symbol!=step->target.symbol || compiled.target.kind!=step->target.kind)
            throw std::runtime_error("Provider changed frozen equation occurrence/target semantics: "+step->occurrence);
        auto stateView=realizeTarget(state,compiled.target);
        const auto existing=std::find_if(result.stateViews.begin(),result.stateViews.end(),
            [&](const auto& view) { return view.symbol==stateView.symbol && view.kind==stateView.kind; });
        if (existing==result.stateViews.end()) result.stateViews.push_back(std::move(stateView));
        else if (existing->storage!=stateView.storage || existing->components!=stateView.components
                 || existing->owner!=stateView.owner)
            throw std::runtime_error("Conflicting storage authorities for STATE view: "+step->target.symbol);
        if (time && !compiled.temporalCapabilities.empty()
            && std::find(compiled.temporalCapabilities.begin(),compiled.temporalCapabilities.end(),time->id())
                ==compiled.temporalCapabilities.end())
            throw std::runtime_error("Unsupported: provider "+found->method+" cannot compile occurrence "
                +step->occurrence+" with selected temporal method "+FDM::toString(time->id()));
        if (!compiled.fragment.children.empty() || !compiled.fragment.operation.empty()) {
            if (!compiled.backendOperation.empty())
                throw std::runtime_error("EquationMethod has both fragment and backend operation authorities.");
            const auto validate=[&](const auto& self,const SolvePlanNode& node)
                -> void {
                if (node.children.empty() && node.operation.empty())
                    throw std::runtime_error("EquationMethod fragment has a leaf without OpId.");
                for (const auto& child:node.children) self(self,child);
            };
            validate(validate,compiled.fragment);
        } else if (compiled.backendOperation.empty()) {
            throw std::runtime_error("EquationMethod supplied neither a fragment nor backend OpId.");
        }
        if (compiled.stateEffects.empty() && compiled.temporalResidual
            && compiled.publishesStageToPhysicalTarget && compiled.target.kind==TargetKind::Physical) {
            for (const auto& written:compiled.writes) if (state.contains(written))
                compiled.stateEffects.push_back({written,false,false});
        }
        result.steps.push_back(std::move(compiled));
    }
    if(recipeScope) {
        expandRecipeScope(result.root,result.steps);
        // Expanding only the declared recipe group changes parent addresses,
        // not occurrence order or numerical selection. Rebind those references.
        parents.clear();steps.clear();std::size_t index=0;
        const auto reparent=[&](const auto& self,ExecutionScope& node,const ExecutionScope* parent,bool stage)->void {
            stage=stage || node.kind==ExecutionKind::StageLoop || node.frequency==ExecutionFrequency::Stage;
            if(node.kind==ExecutionKind::EquationCall || node.kind==ExecutionKind::BlockSolve) {
                node.step.stageScoped=stage;steps.push_back(&node.step);parents.push_back(parent);
                result.steps.at(index++).source=node.step;return;
            }
            for(auto& child:node.children)self(self,child,&node,stage);
        };
        reparent(reparent,result.root,nullptr,false);
    }
    // Fusion is explicitly authored by a provider contract and remains local
    // to adjacent sibling occurrences. Source WHAT/HOW and both targets survive.
    std::vector<std::string> fusedKeys;
    for (std::size_t i=0;i<result.steps.size();++i) {
        const auto& first=result.steps[i];
        if (first.fusionKey.empty()) continue;
        if (first.temporalResidual || first.fusionMembers.empty()
            || std::find(fusedKeys.begin(),fusedKeys.end(),first.fusionKey)!=fusedKeys.end())
            throw std::runtime_error("Unsupported duplicate/incompatible provider fusion: "+first.fusionKey);
        fusedKeys.push_back(first.fusionKey);
        const auto fragmentKey=[](const auto& self,const SolvePlanNode& node)->std::string {
            std::string key=std::to_string(static_cast<int>(node.kind))+":"+node.operation+":"+node.provider
                +":"+std::to_string(node.repetitions)+":"+std::to_string(node.minimumIterations)
                +":"+node.terminationSignal+":"+node.unsupportedReason;
            for (const auto& child:node.children) key+="["+self(self,child)+"]";
            return key;
        };
        const auto n=first.fusionMembers.size();
        for (std::size_t j=0;j<n;++j)
            if (std::any_of(first.fusionMembers.begin(),first.fusionMembers.begin()+j,
                    [&](const auto& previous) { return previous.equation==first.fusionMembers[j].equation
                        && previous.target==first.fusionMembers[j].target; }))
                throw std::runtime_error("Provider fusion declares a duplicate member: "+first.fusionKey);
        if (i+n>result.steps.size() || !parents[i] || parents[i]->kind!=ExecutionKind::Sequence)
            throw std::runtime_error("Provider fusion requires a complete adjacent pair/group: "+first.fusionKey);
        for (std::size_t j=0;j<n;++j) {
            const auto& member=result.steps[i+j];
            const auto sameMembers=member.fusionMembers.size()==n && std::equal(
                first.fusionMembers.begin(),first.fusionMembers.end(),member.fusionMembers.begin(),
                [](const auto& a,const auto& b) { return a.equation==b.equation && a.target==b.target; });
            if (fragmentKey(fragmentKey,member.fragment)!=fragmentKey(fragmentKey,first.fragment)
                || parents[i+j]!=parents[i] || !sameMembers || member.fusionKey!=first.fusionKey
                || member.equationMethod!=first.equationMethod || member.backendOperation!=first.backendOperation
                || member.backendProvider!=first.backendProvider
                || member.source.equation!=first.fusionMembers[j].equation
                || member.target.symbol!=first.fusionMembers[j].target)
                throw std::runtime_error("Provider fusion has missing, duplicate, nonadjacent or wrong-target member: "+first.fusionKey);
        }
        // A Commit or other source scope between calls cannot disappear in fusion.
        const auto& siblings=parents[i]->children;
        const auto at=std::find_if(siblings.begin(),siblings.end(),[&](const auto& node) {
            return node.kind==ExecutionKind::EquationCall && node.step.occurrence==first.source.occurrence;
        });
        if (at==siblings.end() || static_cast<std::size_t>(siblings.end()-at)<n)
            throw std::runtime_error("Provider fusion crosses an execution scope.");
        for (std::size_t j=0;j<n;++j)
            if (at[j].kind!=ExecutionKind::EquationCall || at[j].step.occurrence!=result.steps[i+j].source.occurrence)
                throw std::runtime_error("Provider fusion requires adjacent HOW siblings.");
        i+=n-1;
    }
    for(auto& call:result.steps)for(const auto& read:call.stateUses)if(read.kind!=TargetKind::Physical) {
        CompiledTarget target;target.symbol=read.symbol;target.kind=read.kind;
        target.workspace="state."+read.symbol+(read.kind==TargetKind::Working?".Working":".Correction");
        auto view=realizeTarget(state,target);
        if(view.owner!=StateViewOwner::CompilerWorkspace)throw std::runtime_error("Working STATE input requires compiler-owned workspace: "+read.symbol);
        const auto found=std::find_if(result.stateViews.begin(),result.stateViews.end(),[&](const auto& item){return item.symbol==view.symbol && item.kind==view.kind;});
        if(found==result.stateViews.end())result.stateViews.push_back(std::move(view));
        else if(found->storage!=view.storage || found->owner!=view.owner)throw std::runtime_error("Conflicting Working STATE resource: "+read.symbol);
    }
    for (auto& call:result.steps) for (auto& read:call.viewReads) {
        const auto bound=std::find_if(result.stateViews.begin(),result.stateViews.end(),[&](const auto& view) {
            return view.symbol==read.symbol && view.kind==read.kind && view.stage==read.stage;
        });
        if (bound==result.stateViews.end() || (bound->owner!=StateViewOwner::NumericalProvider && bound->owner!=StateViewOwner::CompilerWorkspace))
            throw std::runtime_error("Missing provider-owned input view: "+call.source.occurrence+" / "+read.symbol+" / "+toString(read.kind));
        read=*bound;
    }
    std::set<std::string> residualIds;
    for(auto& residual:result.residuals) {
        if(residual.id.empty() || !residualIds.insert(residual.id).second || !std::isfinite(residual.absoluteTolerance)
            || !std::isfinite(residual.relativeTolerance) || residual.absoluteTolerance<0 || residual.relativeTolerance<0
            || (residual.absoluteTolerance==0 && residual.relativeTolerance==0)
            || !std::isfinite(residual.scale) || residual.scale<=0)
            throw std::runtime_error("Invalid/duplicate residual definition: "+residual.id);
        std::vector<const CompiledEquationCall*> matches;
        for(const auto& call:result.steps)if(call.source.occurrence==residual.occurrence)matches.push_back(&call);
        if(matches.empty())for(const auto& call:result.steps)if(call.target.symbol==residual.occurrence || call.source.equation==residual.occurrence)matches.push_back(&call);
        if(matches.size()!=1)throw std::runtime_error("Residual evaluation object is missing or ambiguous: "+residual.id+" / "+residual.occurrence);
        const auto& call=*matches.front();
        if(std::find(call.residualCapabilities.begin(),call.residualCapabilities.end(),residual.kind)==call.residualCapabilities.end())
            throw std::runtime_error("Selected provider cannot compute requested residual kind: "+residual.id+" / "+call.equationMethod);
        if(call.backendOperation.empty())throw std::runtime_error("Selected residual has no bounded evaluation operation: "+residual.id);
        residual.occurrence=call.source.occurrence;residual.producerOperation=call.backendOperation;
        residual.provider=call.backendProvider;residual.target=call.target.symbol;
    }
    std::set<std::string> usedResiduals;
    const auto controls=[&](const auto& self,const ExecutionScope& scope,std::vector<std::string> activeSignals)->void {
        if(!scope.residuals.empty() && (scope.kind!=ExecutionKind::Loop || !scope.terminationSignal.empty()))
            throw std::runtime_error("Typed residual control requires an ordinary Loop without another termination signal: "+scope.id);
        if(!scope.terminationSignal.empty()) {
            if(std::find(activeSignals.begin(),activeSignals.end(),scope.terminationSignal)!=activeSignals.end())
                throw std::runtime_error("Nested loops cannot share a termination signal: "+scope.terminationSignal);
            activeSignals.push_back(scope.terminationSignal);
        }
        std::set<std::string> local;
        for(const auto& id:scope.residuals) {
            const auto at=std::find_if(result.residuals.begin(),result.residuals.end(),[&](const auto& r){return r.id==id;});
            if(at==result.residuals.end() || !local.insert(id).second)throw std::runtime_error("Missing/duplicate Loop residual: "+id);
            const auto owns=[&](const auto& walk,const ExecutionScope& n)->bool {
                if((n.kind==ExecutionKind::EquationCall || n.kind==ExecutionKind::BlockSolve) && n.step.occurrence==at->occurrence)return true;
                return std::any_of(n.children.begin(),n.children.end(),[&](const auto& child){return walk(walk,child);});
            };
            if(!owns(owns,scope))throw std::runtime_error("Residual producer is outside its Loop scope: "+id);
            usedResiduals.insert(id);
        }
        for(const auto& child:scope.children)self(self,child,activeSignals);
    };
    controls(controls,result.root,{});
    if(usedResiduals.size()!=result.residuals.size())throw std::runtime_error("Unconsumed residual definition; bind every criterion to a HOW Loop.");
    compileFieldResources(state,result);
    bindCoefficientContracts(state,result);
    for(const auto& call:result.steps) {
        const auto check=[&](const BoundOperation& action) {
            if(action.operation.empty() || action.provider.empty() || action.operation==OpIds::TimeCommit)
                throw std::runtime_error("Invalid numerical input/output action: "+call.source.occurrence);
        };
        for(const auto* actions:{&call.stepInputs,&call.iterationInputs,&call.iterationOutputs,&call.groupInputs})for(const auto& action:*actions)check(action);
        if(!call.physicalPublication.operation.empty())check(call.physicalPublication);
        for(const auto& operation:call.operations)if(operation.operation==OpIds::TimeCommit)
            throw std::runtime_error("Local numerical methods cannot declare physical clock operations: "+call.source.occurrence);
        const auto local=[&](const auto& self,const SolvePlanNode& node)->void {
            if(node.operation==OpIds::TimeCommit)throw std::runtime_error("Local numerical fragments cannot advance the physical clock: "+call.source.occurrence);
            for(const auto& child:node.children)self(self,child);
        };
        local(local,call.fragment);
        if(call.backendOperation==OpIds::TimeCommit)throw std::runtime_error("A numerical equation call cannot be the physical clock: "+call.source.occurrence);
    }
    bool needsCommit=false;for(const auto& call:result.steps)needsCommit=needsCommit || call.temporalResidual || !call.physicalPublication.operation.empty();
    const auto hasCommit=[&](const auto& self,const ExecutionScope& scope)->bool {
        if(scope.kind==ExecutionKind::Commit && scope.commitSources.empty())return true;
        return std::any_of(scope.children.begin(),scope.children.end(),[&](const auto& child){return self(self,child);});
    };
    if(needsCommit && !hasCommit(hasCommit,result.root))throw std::runtime_error("Physical-time HOW requires an explicit Commit.");
    result.loweredRoot=lowerTree(result);
    const auto overlap=[](const StateSymbol& a,const StateSymbol& b) {
        return a.id==b.id || (!a.storageKey.empty() && a.storageBinding==b.storageBinding && a.storageKey==b.storageKey
            && a.componentOffset<b.componentOffset+b.components && b.componentOffset<a.componentOffset+a.components);
    };
    std::vector<CompiledEquationCall*> transient;
    for (auto& call:result.steps) if (call.temporalResidual) transient.push_back(&call);
    const bool stageLocal=std::any_of(result.steps.begin(),result.steps.end(),[](const auto& c){return c.stageLocal;});
    const auto* stages=stageScope(result);
    const bool stageBody=stages!=nullptr;
    if(stageLocal && transient.empty())throw std::runtime_error("Stage algebraic evaluation requires an implemented temporal group");
    if (!transient.empty() && time) {
        if (!temporalMethod || temporalMethod->id()!=time->id())
            throw std::runtime_error("Explicit calls require a selected temporal recipe.");
        const auto compatible=[&](const auto& self,const ExecutionScope& node)->bool {
            if (node.kind==ExecutionKind::StageLoop) return stageBody;
            if (node.kind==ExecutionKind::Loop) {
                const auto temporalInside=[&](const auto& walk,const ExecutionScope& scope)->bool {
                    if (scope.kind==ExecutionKind::EquationCall)
                        return std::any_of(result.steps.begin(),result.steps.end(),[&](const auto& call) {
                            return call.source.occurrence==scope.step.occurrence && call.temporalResidual;
                        });
                    return std::any_of(scope.children.begin(),scope.children.end(),[&](const auto& child) { return walk(walk,child); });
                };
                if (temporalInside(temporalInside,node)) return false;
            }
            if (node.kind==ExecutionKind::Commit && std::none_of(result.root.children.begin(),result.root.children.end(),
                    [&](const auto& child) { return &child==&node; })) return false;
            for (const auto& child:node.children) if (!self(self,child)) return false;
            return true;
        };
        if (!compatible(compatible,result.root))
            throw std::runtime_error("Implementation capability missing: fused temporal provider cannot lower the authored HOW topology.");
        CompiledEquationCall fused;
        fused.temporalResidual=true;
        fused.backendOperation=transient.front()->backendOperation;
        fused.backendProvider=transient.front()->backendProvider;
        fused.temporalPreparation=transient.front()->temporalPreparation;
        fused.temporalStepSize=transient.front()->temporalStepSize;
        fused.temporalSnapshot=transient.front()->temporalSnapshot;
        const bool phased=stageBody || std::any_of(transient.begin(),transient.end(),[](const auto* call) {
            return call->synchronousStages;
        });
        for (auto* call:transient) {
            if (!phased && (call->backendOperation!=fused.backendOperation
                || call->backendProvider!=fused.backendProvider
                || call->temporalPreparation!=fused.temporalPreparation
                || call->temporalStepSize!=fused.temporalStepSize
                || call->temporalSnapshot!=fused.temporalSnapshot))
                throw std::runtime_error("Unsupported: temporal providers declare incompatible fusion backends.");
            fused.calls.insert(fused.calls.end(),call->calls.begin(),call->calls.end());
            call->temporalMethod=FDM::toString(time->id());
        }
        if (phased) {
            std::vector<std::string> written;
            for (const auto* call:transient) {
                if (call->temporalStagePrepare.empty() || call->temporalRhs.empty()
                    || call->temporalAdvance.empty() || call->temporalPublish.empty())
                    throw std::runtime_error("Unsupported temporal participant: split phase capability missing.");
                if (std::find(written.begin(),written.end(),call->target.symbol)!=written.end())
                    throw std::runtime_error("Duplicate temporal physical target: "+call->target.symbol);
                for (const auto& previous:written) {
                    const auto& a=state.at(previous);const auto& b=state.at(call->target.symbol);
                    if (overlap(a,b))
                        throw std::runtime_error("Overlapping temporal physical storage: "+previous+" / "+b.id);
                }
                written.push_back(call->target.symbol);
                const auto identity=call->temporalOwner.empty()?call->source.occurrence:call->temporalOwner;
                auto found=std::find_if(fused.temporalParticipants.begin(),fused.temporalParticipants.end(),
                    [&](const auto& p) {return p.identity==identity;});
                if (found==fused.temporalParticipants.end()) {
                    TemporalParticipant p;p.identity=identity;p.provider=call->backendProvider;
                    p.preparation=call->temporalPreparation;p.stepSize=call->temporalStepSize;
                    p.snapshot=call->temporalSnapshot;p.prepareStage=call->temporalStagePrepare;
                    p.coefficientPrepare=call->coefficientPrepare;p.coefficientPublish=call->coefficientPublish;
                    p.coefficientProvides=call->coefficientProvides;p.coefficientStates=call->coefficientStates;
                    p.coefficientVersion=call->coefficientVersion;p.coefficientInputs=call->coefficientInputs;
                    p.rhs=call->temporalRhs;p.advance=call->temporalAdvance;p.publish=call->temporalPublish;
                    fused.temporalParticipants.push_back(std::move(p));found=fused.temporalParticipants.end()-1;
                } else if (found->provider!=call->backendProvider || found->prepareStage!=call->temporalStagePrepare
                    || found->rhs!=call->temporalRhs || found->advance!=call->temporalAdvance || found->publish!=call->temporalPublish
                    || found->preparation!=call->temporalPreparation || found->stepSize!=call->temporalStepSize || found->snapshot!=call->temporalSnapshot
                    || found->coefficientPrepare!=call->coefficientPrepare || found->coefficientPublish!=call->coefficientPublish)
                    throw std::runtime_error("Conflicting temporal numerical owner: "+identity);
                for(const auto& use:call->coefficientUses)if(std::none_of(found->coefficientUses.begin(),found->coefficientUses.end(),[&](const auto& item){return item.symbol==use.symbol;}))found->coefficientUses.push_back(use);
                found->calls.insert(found->calls.end(),call->calls.begin(),call->calls.end());
                found->targets.push_back(call->target.symbol);
                found->stageRequires.insert(found->stageRequires.end(),call->stageRequires.begin(),call->stageRequires.end());
                found->stageProvides.insert(found->stageProvides.end(),call->stageProvides.begin(),call->stageProvides.end());
            }
            // A shared numerical owner is one RHS, so its mathematical calls must
            // be adjacent siblings even when no algebraic Stage call is present.
            for(const auto& participant:fused.temporalParticipants) {
                std::vector<std::size_t> members;
                for(std::size_t i=0;i<result.steps.size();++i)if(result.steps[i].temporalResidual
                    && (result.steps[i].temporalOwner.empty()?result.steps[i].source.occurrence:result.steps[i].temporalOwner)==participant.identity)members.push_back(i);
                const auto* parent=parents.at(members.front());
                if(!parent)throw std::runtime_error("Temporal numerical owner requires a sibling scope: "+participant.identity);
                const auto first=std::find_if(parent->children.begin(),parent->children.end(),[&](const auto& child){return child.kind==ExecutionKind::EquationCall && child.step.occurrence==result.steps[members.front()].source.occurrence;});
                if(first==parent->children.end() || static_cast<std::size_t>(parent->children.end()-first)<members.size())
                    throw std::runtime_error("Temporal numerical owner requires a complete adjacent HOW group: "+participant.identity);
                for(std::size_t i=0;i<members.size();++i)if(parents[members[i]]!=parent || first[i].kind!=ExecutionKind::EquationCall
                    || first[i].step.occurrence!=result.steps[members[i]].source.occurrence)
                    throw std::runtime_error("Temporal numerical owner crosses a HOW operation/scope: "+participant.identity+"; HOW is not reordered");
            }
            std::vector<std::string> available;
            for(const auto& p:result.steps)for(const auto& resource:p.coefficientProvides) {
                if(std::find(available.begin(),available.end(),resource)!=available.end())throw std::runtime_error("Duplicate same-stage coefficient producer: "+resource);
                available.push_back(resource);
            }
            for(const auto& participant:fused.temporalParticipants) {
                for(const auto& resource:participant.stageRequires)
                    if(std::find(available.begin(),available.end(),resource)==available.end())
                        throw std::runtime_error("Temporal occurrence '"+participant.identity+"' provider '"+participant.provider
                            +"' requires same-stage RHS resource '"+resource+"' before its consumer; available producer set does not satisfy authored HOW.");
                for(const auto& resource:participant.stageProvides) {
                    if(std::find(available.begin(),available.end(),resource)!=available.end())
                        throw std::runtime_error("Duplicated stage resource producer: "+resource+" at "+participant.identity);
                    available.push_back(resource);
                }
            }
            if (result.root.children.empty() || result.root.children.back().kind!=ExecutionKind::Commit
                || std::count_if(result.root.children.begin(),result.root.children.end(),[](const auto& n){return n.kind==ExecutionKind::Commit;})!=1)
                throw std::runtime_error("Common temporal group requires one terminal physical Commit.");
            result.temporalParticipants=fused.temporalParticipants;
            validateCoefficientInputs(state,result);
            for(const char* op:{TemporalOps::Dt,TemporalOps::Open,TemporalOps::Ready,TemporalOps::RhsReady,
                    TemporalOps::Close,TemporalOps::PublishReady})
                result.operations.push_back({op,op,std::string_view(op)==TemporalOps::PublishReady?OperationStage::StepCommit:OperationStage::Prepare,
                    {OperationCapability::TemporalSynchronization},{OriginKind::Generated,"common temporal recipe"}});
        }
        std::vector<SolvePlanNode> prefix,suffix;
        if(stageBody) {
            if(!phased || stages->repetitions!=1 || stages->minimumIterations!=1 || !stages->terminationSignal.empty())throw std::runtime_error("Stage algebraic HOW requires one recipe-owned Stage scope and terminal Commit");
            const auto stageShape=[&](const auto& self,const ExecutionScope& scope,bool loop)->void {
                if(scope.kind==ExecutionKind::EquationCall) {
                    const auto at=std::find_if(result.steps.begin(),result.steps.end(),[&](const auto& c){return c.source.occurrence==scope.step.occurrence;});
                    if(at==result.steps.end() || (loop && !at->stageLocal))
                        throw std::runtime_error("Explicit StageLoop iteration supports local algebraic calls only; temporal RHS cannot repeat: "+scope.step.occurrence);
                    return;
                }
                if(scope.kind!=ExecutionKind::Sequence && scope.kind!=ExecutionKind::Loop)
                    throw std::runtime_error("Unsupported control node inside explicit StageLoop: "+scope.id);
                for(const auto& child:scope.children)self(self,child,loop || scope.kind==ExecutionKind::Loop);
            };
            for(const auto& child:stages->children)stageShape(stageShape,child,false);
            std::vector<std::string> ready;
            for(const auto& call:result.steps) {
                if(!call.source.stageScoped)continue;
                if(!call.temporalResidual && !call.stageLocal)throw std::runtime_error("Selected implementation cannot execute in Stage scope: "+call.source.occurrence);
                for(const auto& use:call.stateUses)if(use.version==StateVersion::Stage && std::any_of(result.steps.begin(),result.steps.end(),[&](const auto& p){return p.stageLocal && p.target.symbol==use.symbol;}) && std::find(ready.begin(),ready.end(),use.symbol)==ready.end())throw std::runtime_error("Missing preceding Stage algebraic publication: "+use.symbol+" consumer="+call.source.occurrence);
                if(call.stageLocal && !call.coefficientProvides.empty())for(const auto& id:call.coefficientStates)ready.push_back(id);
                if(call.stageLocal && call.coefficientProvides.empty()) {
                    const auto& target=state.at(call.target.symbol);
                    if(std::any_of(transient.begin(),transient.end(),[&](const auto* c){return overlap(state.at(c->target.symbol),target);})
                        || std::any_of(ready.begin(),ready.end(),[&](const auto& id){return overlap(state.at(id),target);})
                        || std::any_of(result.temporalParticipants.begin(),result.temporalParticipants.end(),[&](const auto& p){return std::any_of(p.coefficientStates.begin(),p.coefficientStates.end(),[&](const auto& id){return overlap(state.at(id),target);});}))
                        throw std::runtime_error("Duplicate Stage state writer/storage: "+call.target.symbol);
                    ready.push_back(call.target.symbol);
                }
            }
        }
        if(stageBody || transient.size()!=result.steps.size()) {
            auto placement=stepPlacement(result);prefix=std::move(placement.first);suffix=std::move(placement.second);
        }
        result.temporalRoot=temporalMethod->compileFragment(*time,fused,prefix);
        if(stageBody) {
            auto at=std::find_if(result.temporalRoot.children.begin(),result.temporalRoot.children.end(),[](const auto& n){return n.kind==PlanNodeKind::StageLoop;});
            if(at==result.temporalRoot.children.end())throw std::runtime_error("Selected recipe omitted its StageLoop");
            auto first=std::find_if(at->children.begin(),at->children.end(),[&](const auto& n){return n.operation==fused.temporalParticipants.front().rhs;});
            auto last=std::find_if(first,at->children.end(),[](const auto& n){return n.operation==TemporalOps::RhsReady;});
            if(first==at->children.end() || last==at->children.end())throw std::runtime_error("Selected recipe omitted its RHS barriers");
            // Preparation/coefficient publication stays before this interval, unchanged.
            first=at->children.erase(first,last);auto evaluations=stageEvaluations(result);
            at->children.insert(first,evaluations.begin(),evaluations.end());
        }
        result.temporalRoot.children.insert(result.temporalRoot.children.end(),suffix.begin(),suffix.end());
        if (phased) {
            SolvePlanNode barrier;barrier.id=barrier.operation=TemporalOps::PublishReady;
            barrier.name=barrier.id;barrier.kind=PlanNodeKind::Commit;barrier.provider=TemporalOps::Provider;
            result.temporalRoot.children.push_back(barrier);
            for (const auto& p:fused.temporalParticipants) {
                SolvePlanNode publish;publish.kind=PlanNodeKind::Commit;publish.id=publish.name=publish.operation=p.publish;
                publish.provider=p.provider;publish.occurrence=p.identity;publish.equationCalls=p.calls;
                result.temporalRoot.children.push_back(std::move(publish));
            }
            for(const auto& p:fused.temporalParticipants)if(!p.coefficientPublish.empty()) {
                SolvePlanNode publish;publish.kind=PlanNodeKind::Commit;publish.id=publish.name=publish.operation=p.coefficientPublish;publish.provider=p.provider;
                result.temporalRoot.children.push_back(std::move(publish));
            }
            for(const auto& call:result.steps)if(call.stageLocal){SolvePlanNode publish;publish.kind=PlanNodeKind::Commit;publish.id=publish.operation=call.temporalPublish;publish.name="publish algebraic relation from committed physical inputs";publish.provider=call.backendProvider;publish.occurrence=call.source.occurrence;publish.equationCalls=call.calls;result.temporalRoot.children.push_back(std::move(publish));}
            barrier.id=barrier.name=barrier.operation=OpIds::TimeCommit;
            result.temporalRoot.children.push_back(barrier);
        }
        // Source Commit positions survive recipe expansion. Local publication
        // kernels implement outputs; only this process compiler appends the clock.
        for (std::size_t i=0;i<result.root.children.size();++i) {
            if (result.root.children[i].kind!=ExecutionKind::Commit) continue;
            if (i+1!=result.root.children.size())
                throw std::runtime_error("Unsupported: fused temporal provider requires a terminal Commit.");
            if (!phased) result.temporalRoot.children.push_back(result.loweredRoot.children.back());
        }
        result.hasTemporalRoot=true;
    }
    validateClockCommit(result.hasTemporalRoot?result.temporalRoot:result.loweredRoot);
    const auto& execution=result.hasTemporalRoot?result.temporalRoot:result.loweredRoot;
    const auto declaresClock=[&](const auto& self,const SolvePlanNode& node)->bool {
        return node.operation==OpIds::TimeCommit || std::any_of(node.children.begin(),node.children.end(),[&](const auto& c){return self(self,c);});
    };
    if(declaresClock(declaresClock,execution))result.operations.push_back({OpIds::TimeCommit,
        "advance the unique physical clock",OperationStage::StepCommit,{OperationCapability::PhysicalTimeCommit},
        {OriginKind::Generated,"HOW physical Commit"}});
    validateDataFlow(state,result);
    if (time) realizeTemporalViews(state,result,time->stageCount());
    compileDerivedViews(state,result.stateViews);
    return result;
}


std::string validateTemporalPlan(const CompiledSolvePlan& plan,int stages) {
    const auto& participants=plan.compiledProgram.temporalParticipants;
    if(participants.empty())return {};
    std::vector<std::pair<std::string,std::string>> expected,actual;
    const auto add=[&](const auto& id,const auto& provider){expected.emplace_back(id,provider);};
    const auto loopKey=[](const SolvePlanNode& n) {
        std::ostringstream key;key<<std::setprecision(17)<<n.id<<':'<<n.repetitions<<':'<<n.minimumIterations
            <<':'<<n.terminationSignal<<':'<<n.residualAll<<':'<<n.failOnMaximum;
        for(const auto& r:n.residuals)key<<':'<<r.id<<':'<<r.occurrence<<':'<<r.producerOperation<<':'<<r.provider
            <<':'<<static_cast<int>(r.kind)<<':'<<static_cast<int>(r.norm)<<':'<<static_cast<int>(r.normalization)
            <<':'<<r.absoluteTolerance<<':'<<r.relativeTolerance<<':'<<r.scale;
        return key.str();
    };
    const auto walkExpected=[&](const auto& self,const SolvePlanNode& n)->void {
        if(n.kind==PlanNodeKind::Loop)add("<Loop>",loopKey(n));
        if(!n.operation.empty())add(n.operation,n.provider);
        for(const auto& child:n.children)self(self,child);
        if(n.kind==PlanNodeKind::Loop)add("</Loop>",n.id);
    };
    std::pair<std::vector<SolvePlanNode>,std::vector<SolvePlanNode>> placement;
    const bool placed=stageScope(plan.compiledProgram) || std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[](const auto& c){return !c.temporalResidual;});
    try{if(placed)placement=stepPlacement(plan.compiledProgram);}catch(const std::exception& e){return e.what();}
    for(const auto& p:participants)for(const auto& op:p.preparation)add(op,p.provider);
    add(TemporalOps::Dt,TemporalOps::Provider);
    for(const auto& n:placement.first)walkExpected(walkExpected,n);
    for(const auto& p:participants)add(p.snapshot,p.provider);
    add("<StageLoop>",std::to_string(stages));add(TemporalOps::Open,TemporalOps::Provider);
    for(const auto& p:participants)add(p.prepareStage,p.provider);
    add(TemporalOps::Ready,TemporalOps::Provider);
    for(const auto& p:participants)if(!p.coefficientPrepare.empty())add(p.coefficientPrepare,p.provider);
    const bool actions=std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[](const auto& c){return c.stageLocal;});
    if(actions) {
        const auto walk=[&](const auto& self,const SolvePlanNode& n)->void {
            if(n.kind==PlanNodeKind::Loop)add("<Loop>",loopKey(n));
            if(!n.operation.empty())add(n.operation,n.provider);
            for(const auto& child:n.children)self(self,child);
            if(n.kind==PlanNodeKind::Loop)add("</Loop>",n.id);
        };
        try{for(const auto& n:stageEvaluations(plan.compiledProgram))walk(walk,n);}
        catch(const std::exception& e){return e.what();}
    }
    else for(const auto& p:participants)add(p.rhs,p.provider);
    add(TemporalOps::RhsReady,TemporalOps::Provider);
    for(const auto& p:participants)add(p.advance,p.provider);
    add(TemporalOps::Close,TemporalOps::Provider);add("</StageLoop>",std::string{});
    for(const auto& n:placement.second)walkExpected(walkExpected,n);
    add(TemporalOps::PublishReady,TemporalOps::Provider);
    for(const auto& p:participants)add(p.publish,p.provider);
    for(const auto& p:participants)if(!p.coefficientPublish.empty())add(p.coefficientPublish,p.provider);
    for(const auto& c:plan.compiledProgram.steps)if(c.stageLocal)add(c.temporalPublish,c.backendProvider);
    add(OpIds::TimeCommit,TemporalOps::Provider);
    bool valid=true;
    const auto walk=[&](const auto& self,const SolvePlanNode& n)->void {
        if(n.kind==PlanNodeKind::StageLoop) actual.emplace_back("<StageLoop>",std::to_string(n.repetitions));
        else if(n.kind==PlanNodeKind::Loop)actual.emplace_back("<Loop>",loopKey(n));
        else if(!n.children.empty() && n.kind!=PlanNodeKind::Sequence)valid=false;
        if(!n.operation.empty())actual.emplace_back(n.operation,n.provider);
        for(const auto& c:n.children)self(self,c);
        if(n.kind==PlanNodeKind::StageLoop)actual.emplace_back("</StageLoop>",std::string{});
        if(n.kind==PlanNodeKind::Loop)actual.emplace_back("</Loop>",n.id);
    };
    walk(walk,plan.root);
    if(!valid || actual!=expected)return "Invalid synchronous temporal plan: phase order, participant, StageLoop recipe or single terminal Commit differs from the compiled contract.";
    return {};
}

void compileSystem(ResolvedSimulationSystem& result,const ExecutionProgram& executionProgram,
        NumericalSelection& selected,
        bool additionalContributions) {
    result.realization=compileStateRealization(result.executableSystem.state,result.executableSystem.constraints);
    const auto temporalMethods=builtinTemporalMethods();
    const auto& temporalMethod=temporalMethods.at(selected.recipes.time.id());
    const auto compiledTime=temporalMethod.compile(selected.recipes.time);
    auto compiledExecution=compileProgram(result.executableSystem,result.executableSystem.state,executionProgram,
        selected.bindings,builtinProviders(),&compiledTime,&temporalMethod);
    const auto diffusion=[](const auto& self,const FormulaExpr& expression)->bool {
        if(expression.kind()==FormulaExpr::Kind::Operator && expression.operation()==Math::Op::Diffusion)return true;
        return std::any_of(expression.arguments().begin(),expression.arguments().end(),[&](const auto& child){return self(self,child);});
    };
    if(!selected.recipes.diffusion && selected.defaultDiffusion)for(const auto& call:compiledExecution.steps) {
        if(!call.spatialTerms && call.sourceMathInputs.empty())continue;
        const auto& equation=result.executableSystem.registry.at(call.source.equation);
        if(diffusion(diffusion,equation.lhs) || diffusion(diffusion,equation.rhs)) {selected.recipes.diffusion=selected.defaultDiffusion;break;}
    }
    for (const auto& id:result.executableSystem.state.solutionVariables()) {
        const auto& symbol=result.executableSystem.state.at(id);
        const bool owned=std::any_of(compiledExecution.steps.begin(),compiledExecution.steps.end(),[&](const auto& step) {
            if (std::find(step.writes.begin(),step.writes.end(),id)!=step.writes.end()) return true;
            return std::any_of(step.target.resources.begin(),step.target.resources.end(),[&](const auto& resource) {
                return resource.storage==symbol.storageKey && resource.access!=ResourceAccessMode::Read
                    && resource.componentOffset<=symbol.componentOffset
                    && resource.componentOffset+resource.components>=symbol.componentOffset+symbol.components;
            });
        });
        if (!owned) throw std::runtime_error("No compiled HOW output owns selected solution STATE '"+id+"'.");
    }
    TermProviderCatalog termProviders=TermProviderCatalog::builtIn();
    for (const auto& descriptor:selected.termProviders) termProviders.addSource(descriptor);
    result.numericalSystem=NumericalCompiler::compile(result.executableSystem,compiledExecution,
        selected.recipes,termProviders,compiledTime,builtinProviders(),selected.operatorBindings);
    result.numericalSystem.dt=selected.dt;
    result.numericalSystem.phaseTransport=selected.phaseTransport;
    result.numericalSystem.pressureFaceCoupling=selected.pressureFaceCoupling;
    result.numericalSystem.pressureOperator=selected.pressureOperator;
    result.numericalSystem.phaseOperator=selected.phaseOperator;
    result.executableSystem.operations.insert(result.executableSystem.operations.end(),
        compiledExecution.operations.begin(),compiledExecution.operations.end());
    // Compiled operation/resource declarations originate in selected providers or explicit HOW.
    for (const auto& call:compiledExecution.steps) {
        for (const auto& operation:call.operations) {
            const auto duplicate=std::find_if(result.executableSystem.operations.begin(),result.executableSystem.operations.end(),
                [&](const auto& existing) { return existing.operation==operation.operation; });
            if (duplicate==result.executableSystem.operations.end()) result.executableSystem.operations.push_back(operation);
        }
    }
    for (auto& step:compiledExecution.steps) {
        for (const auto& call:step.calls) {
            for (const auto& term:result.numericalSystem.operators) {
                if ((term.formulaId!=call.equation
                        && std::find(step.sourceMathInputs.begin(),
                            step.sourceMathInputs.end(),term.formulaId)
                           ==step.sourceMathInputs.end())
                    || term.provider.empty()) continue;
                if (std::find(step.operatorBindings.begin(),
                        step.operatorBindings.end(),term.provider)
                    ==step.operatorBindings.end())
                    step.operatorBindings.push_back(term.provider);
            }
        }
        if (step.temporalResidual && step.spatialTerms) {
            if (step.spatialTerms && step.operatorBindings.empty())
                throw std::runtime_error("ConservativeResidual has no bound spatial provider.");
            step.requirements.push_back("halo depth="
                +std::to_string(result.numericalSystem.requiredHaloWidth));
            step.requirements.push_back("canonical face COPY before residual assembly");
        }
    }
    result.solvePlan = SolvePlanner::compile(
        result.numericalSystem.time.recipe,
        compiledExecution);
    result.solvePlan.sourceProgram=executionProgram;
    result.solvePlan.sourceProgram.root=compiledExecution.root;
    result.solvePlan.compiledProgram=std::move(compiledExecution);
    result.runtime.capabilities=compileExecutionCapabilities(
        result.executableSystem,&result.solvePlan);
    if(!result.solvePlan.compiledProgram.temporalParticipants.empty()) {
        const auto reason=validateTemporalPlan(result.solvePlan,result.numericalSystem.time.recipe.stageCount());
        if(!reason.empty())throw std::runtime_error(reason);
    }
    validateDirectProgram(result.executableSystem,result.numericalSystem,result.solvePlan);
    result.runtime.operationBindings=compileOperationBindings(result.executableSystem,
        result.numericalSystem,result.solvePlan,additionalContributions);
    result.runtime.report=reportOperationBindings(result.solvePlan,result.runtime.operationBindings);
    validateMathematicalCoverage(result.executableSystem,result.solvePlan,result.runtime.report.status==RuntimeStatus::Runnable);
}

} // namespace SF::System::Lowering
