/// @file SF_methodObjects.cpp
/// @brief Built-in temporal mathematics and equation realization contracts.

#include "SF_methodObjects.h"
#include "SF_stateRealization.h"
#include "SF_executionContract.h"
#include "core/system/SF_operationIds.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::System {
void TemporalMethodRegistry::add(const ITemporalMethod& method) {
    if (std::any_of(methods_.begin(),methods_.end(),
        [&](const ITemporalMethod* existing) { return existing->id()==method.id(); }))
        throw std::runtime_error("Duplicate TemporalMethod registration.");
    methods_.push_back(&method);
}

const ITemporalMethod& TemporalMethodRegistry::at(FDM::TimeRecipeId id) const {
    for (const auto* method:methods_) if (method->id()==id) return *method;
    throw std::runtime_error("Selected TemporalMethod has no registered implementation.");
}

void ProviderRegistry::add(const IProvider& method) {
    if (std::any_of(methods_.begin(),methods_.end(),
        [&](const IProvider* existing) { return existing->id()==method.id(); }))
        throw std::runtime_error("Duplicate EquationMethod registration.");
    methods_.push_back(&method);
}

const IProvider& ProviderRegistry::at(std::string_view id) const {
    for (const auto* method:methods_) if (method->id()==id) return *method;
    throw std::runtime_error("EquationMethod '"+std::string(id)
        +"' has no registered implementation.");
}

CompiledExecutionProgram compileExecutionProgram(
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
    result.root=root;
    std::vector<const EquationCall*> steps;
    std::vector<const ExecutionScope*> parents;
    const auto collect=[&](const auto& self,const ExecutionScope& node,const ExecutionScope* parent)->void {
        if (node.kind==ExecutionKind::EquationCall) {
            if (!node.children.empty())
                throw std::runtime_error("HOW Step cannot contain child nodes.");
            steps.push_back(&node.step);
            parents.push_back(parent);
            return;
        }
        if (node.kind==ExecutionKind::Commit && !node.children.empty())
            throw std::runtime_error("HOW Commit cannot contain child nodes.");
        if (node.kind==ExecutionKind::Loop
            && (node.repetitions<1 || node.children.empty()
                || node.minimumIterations<1 || node.minimumIterations>node.repetitions))
            throw std::runtime_error("HOW Loop requires a body and valid iteration limits.");
        for (const auto& child:node.children) self(self,child,&node);
    };
    collect(collect,result.root,nullptr);
    std::vector<std::string> occurrences;
    for (const auto* step:steps) {
        if (std::find(occurrences.begin(),occurrences.end(),step->occurrence)!=occurrences.end())
            throw std::runtime_error("Duplicate execution occurrence: "+step->occurrence);
        occurrences.push_back(step->occurrence);
    }
    for (const auto& binding:bindings) {
        if (binding.equation.empty() || binding.method.empty())
            throw std::runtime_error("Numerical binding requires equation and provider.");
        if (std::none_of(steps.begin(),steps.end(),
            [&](const EquationCall* call) {
                return (binding.equation=="*" || call->equation==binding.equation)
                    && (binding.occurrence.empty() || call->occurrence==binding.occurrence);
            }))
            throw std::runtime_error("Numerical binding addresses no execution occurrence: "+binding.equation);
    }
    for (const auto* step:steps) {
        const NumericalBinding* found=nullptr;
        int precedence=-1;
        int countAtPriority=0;
        for (const auto& binding:bindings) {
            if (binding.equation!="*" && binding.equation!=step->equation) continue;
            if (!binding.occurrence.empty() && binding.occurrence!=step->occurrence) continue;
            const int rank=!binding.occurrence.empty() ? 2 : binding.equation!="*" ? 1 : 0;
            if (rank>precedence) { found=&binding;precedence=rank;countAtPriority=1; }
            else if (rank==precedence) ++countAtPriority;
        }
        if (countAtPriority>1)
            throw std::runtime_error("Conflicting highest-priority numerical bindings for occurrence: "+step->occurrence);
        if (!found)
            throw std::runtime_error("Implementation capability missing: no provider for equation occurrence "+step->equation);
        system.registry.at(step->equation);
        if (step->target.kind!=TargetKind::Workspace) state.at(step->target.symbol);
        const auto& selectedProvider=registry.at(found->method);
        auto compiled=selectedProvider.compile(system,*step,*found);
        compiled.backendProvider=selectedProvider.runtimeProvider();
        const auto ownFragment=[&](const auto& self,SolvePlanNode& node)->void {
            if (!node.operation.empty())
                node.provider=selectedProvider.operationProvider(node.operation);
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
                +":"+node.terminationSignal+":"+node.unsupportedReason+":"+std::to_string(node.legacyAdapter);
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
    result.loweredRoot=compileMethodProgram(result,&registry);
    std::vector<CompiledEquationCall*> transient;
    for (auto& call:result.steps) if (call.temporalResidual) transient.push_back(&call);
    if (!transient.empty() && time) {
        if (!temporalMethod || temporalMethod->id()!=time->id())
            throw std::runtime_error("Explicit calls require a selected temporal recipe.");
        const auto compatible=[&](const auto& self,const ExecutionScope& node)->bool {
            if (node.kind==ExecutionKind::StageLoop) return false;
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
        const bool phased=std::any_of(transient.begin(),transient.end(),[](const auto* call) {
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
                    if (a.storageBinding==b.storageBinding && a.storageKey==b.storageKey
                        && a.componentOffset<b.componentOffset+b.components && b.componentOffset<a.componentOffset+a.components)
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
                    p.rhs=call->temporalRhs;p.advance=call->temporalAdvance;p.publish=call->temporalPublish;
                    fused.temporalParticipants.push_back(std::move(p));found=fused.temporalParticipants.end()-1;
                } else if (found->provider!=call->backendProvider || found->prepareStage!=call->temporalStagePrepare
                    || found->rhs!=call->temporalRhs || found->advance!=call->temporalAdvance || found->publish!=call->temporalPublish)
                    throw std::runtime_error("Conflicting temporal numerical owner: "+identity);
                found->calls.insert(found->calls.end(),call->calls.begin(),call->calls.end());
                found->targets.push_back(call->target.symbol);
            }
            if (result.root.children.empty() || result.root.children.back().kind!=ExecutionKind::Commit
                || std::count_if(result.root.children.begin(),result.root.children.end(),[](const auto& n){return n.kind==ExecutionKind::Commit;})!=1)
                throw std::runtime_error("Common temporal group requires one terminal physical Commit.");
            result.temporalParticipants=fused.temporalParticipants;
            // The shared physical clock belongs to the temporal group even
            // when a participant still declares its standalone commit.
            for(auto& call:result.steps)
                call.operations.erase(std::remove_if(call.operations.begin(),call.operations.end(),
                    [](const auto& op){return op.operation==OpIds::TimeCommit;}),call.operations.end());
            for (const char* op:{TemporalOps::Dt,TemporalOps::Open,TemporalOps::Ready,TemporalOps::RhsReady,
                    TemporalOps::Close,TemporalOps::PublishReady,OpIds::TimeCommit})
                result.steps.front().operations.push_back({op,op,
                    (std::string_view(op)==OpIds::TimeCommit || std::string_view(op)==TemporalOps::PublishReady)
                        ?OperationStage::StepCommit:OperationStage::Prepare,
                    {OperationCapability::TemporalSynchronization},{OriginKind::Generated,"common temporal recipe"}});
        }
        std::vector<SolvePlanNode> prefix,suffix;
        const bool mixed=transient.size()!=result.steps.size();
        if (mixed) {
            // Ordered prefix -> contiguous temporal group -> ordered suffix -> Commit.
            // No domain identity enters the compiler and no authored call is moved.
            bool seenTemporal=false,seenSuffix=false;
            std::size_t prefixCount=0,suffixCount=0;
            for (const auto& call:result.steps) {
                if (call.temporalResidual) {
                    if (seenSuffix) throw std::runtime_error("Unsupported mixed temporal topology: temporal calls must be contiguous.");
                    seenTemporal=true;continue;
                }
                if (seenTemporal) seenSuffix=true;

            }
            const auto groupedCalls=[](const auto& self,const ExecutionScope& node)->bool {
                if (node.kind==ExecutionKind::EquationCall) return true;
                return (node.kind==ExecutionKind::Sequence || node.kind==ExecutionKind::Loop)
                    && !node.children.empty() && std::all_of(node.children.begin(),node.children.end(),
                        [&](const auto& child) { return self(self,child); });
            };
            const auto callCount=[](const auto& self,const ExecutionScope& node)->std::size_t {
                if (node.kind==ExecutionKind::EquationCall) return 1;
                std::size_t count=0;
                for (const auto& child:node.children) count+=self(self,child);
                return count;
            };
            if (result.root.kind!=ExecutionKind::Sequence || result.root.children.empty()
                || result.root.children.back().kind!=ExecutionKind::Commit
                || std::any_of(result.root.children.begin(),result.root.children.end()-1,[&](const auto& node) {
                    return !groupedCalls(groupedCalls,node);
                }))
                throw std::runtime_error("Unsupported mixed temporal topology: requires ordered calls/groups and terminal Commit.");
            // A source group remains one root node even if it contains several
            // local operations. Count scopes, not flattened equation leaves.
            std::size_t nextCall=0,temporalCount=0;
            for (const auto& node:result.root.children) {
                if (node.kind==ExecutionKind::Commit) continue;
                const auto count=callCount(callCount,node);
                const bool temporal=result.steps[nextCall].temporalResidual;
                for (std::size_t i=0;i<count;++i)
                    if (result.steps[nextCall+i].temporalResidual!=temporal)
                        throw std::runtime_error("Unsupported mixed temporal topology: a group crosses temporal placement.");
                const bool fusedContinuation=node.kind==ExecutionKind::EquationCall && nextCall>0
                    && !result.steps[nextCall].fusionKey.empty()
                    && result.steps[nextCall-1].fusionKey==result.steps[nextCall].fusionKey;
                if (!fusedContinuation) {
                    if (temporal) ++temporalCount;
                    else if (!temporalCount) ++prefixCount;
                    else ++suffixCount;
                }
                nextCall+=count;
            }
            const auto expectedChildren=prefixCount+temporalCount+suffixCount+1;
            if (result.loweredRoot.children.size()!=expectedChildren)
                throw std::runtime_error("Unsupported mixed temporal topology: root lifecycle decorations require explicit temporal placement.");
            prefix.assign(result.loweredRoot.children.begin(),result.loweredRoot.children.begin()+prefixCount);
            suffix.assign(result.loweredRoot.children.begin()+prefixCount+temporalCount,result.loweredRoot.children.end()-1);
        }
        result.temporalRoot=temporalMethod->compileFragment(*time,fused,prefix);
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
            barrier.id=barrier.name=barrier.operation=OpIds::TimeCommit;
            result.temporalRoot.children.push_back(barrier);
        }
        // Source commit positions survive temporal lowering. The provider owns
        // each commit's implementation; WHICH never invents a source commit.
        for (std::size_t i=0;i<result.root.children.size();++i) {
            if (result.root.children[i].kind!=ExecutionKind::Commit) continue;
            if (i+1!=result.root.children.size())
                throw std::runtime_error("Unsupported: fused temporal provider requires a terminal Commit.");
            if (!phased) result.temporalRoot.children.push_back(result.loweredRoot.children.back());
        }
        result.hasTemporalRoot=true;
    }
    validateDataFlow(state,result);
    if (time) realizeTemporalViews(state,result,time->stageCount());
    return result;
}

SolvePlanNode compileMethodProgram(const CompiledExecutionProgram& program, const ProviderRegistry* providers) {
    if (!program.loweredRoot.children.empty()) return program.loweredRoot;
    std::size_t next=0;
    const auto lower=[&](const auto& self,const ExecutionScope& source)->SolvePlanNode {
        SolvePlanNode result;
        switch (source.kind) {
            case ExecutionKind::Sequence: result.kind=PlanNodeKind::Sequence;break;
            case ExecutionKind::Loop: result.kind=PlanNodeKind::Loop;break;
            case ExecutionKind::StageLoop: result.kind=PlanNodeKind::StageLoop;break;
            case ExecutionKind::EquationCall: result.kind=PlanNodeKind::Update;break;
            case ExecutionKind::Commit: result.kind=PlanNodeKind::Sequence;break;
        }
        result.repetitions=source.repetitions;
        result.terminationSignal=source.terminationSignal;
        result.minimumIterations=source.minimumIterations;
        result.id=!source.id.empty() ? source.id
            : source.kind==ExecutionKind::EquationCall ? source.step.equation
            : std::string(toString(source.kind));
        result.name=result.id;
        if (source.kind==ExecutionKind::EquationCall) {
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
        if (providers) {
            std::vector<std::string> visited;
            for (const auto& call:program.steps) {
                if (std::find(visited.begin(),visited.end(),call.equationMethod)!=visited.end()) continue;
                visited.push_back(call.equationMethod);
                providers->at(call.equationMethod).lowerLifecycle(source,&source==&program.root,program.steps,result);
            }
        }
        if (source.kind==ExecutionKind::Commit && result.children.empty())
            throw std::runtime_error("Unsupported: no selected provider implements Commit scope "+source.id);
        return result;
    };
    auto result=lower(lower,program.root);
    if (next!=program.steps.size())
        throw std::runtime_error("Compiled HOW contains unconsumed method steps.");
    return result;
}

std::string validateTemporalPlan(const CompiledSolvePlan& plan,int stages) {
    const auto& participants=plan.compiledProgram.temporalParticipants;
    if(participants.empty())return {};
    std::vector<std::pair<std::string,std::string>> expected,actual;
    const auto add=[&](const auto& id,const auto& provider){expected.emplace_back(id,provider);};
    for(const auto& p:participants)for(const auto& op:p.preparation)add(op,p.provider);
    add(TemporalOps::Dt,TemporalOps::Provider);
    for(const auto& p:participants)add(p.snapshot,p.provider);
    add("<StageLoop>",std::to_string(stages));add(TemporalOps::Open,TemporalOps::Provider);
    for(const auto& p:participants)add(p.prepareStage,p.provider);
    add(TemporalOps::Ready,TemporalOps::Provider);
    for(const auto& p:participants)add(p.rhs,p.provider);
    add(TemporalOps::RhsReady,TemporalOps::Provider);
    for(const auto& p:participants)add(p.advance,p.provider);
    add(TemporalOps::Close,TemporalOps::Provider);add("</StageLoop>",std::string{});
    add(TemporalOps::PublishReady,TemporalOps::Provider);
    for(const auto& p:participants)add(p.publish,p.provider);
    add(OpIds::TimeCommit,TemporalOps::Provider);
    bool valid=true;
    const auto walk=[&](const auto& self,const SolvePlanNode& n)->void {
        if(n.kind==PlanNodeKind::StageLoop) actual.emplace_back("<StageLoop>",std::to_string(n.repetitions));
        else if(!n.children.empty() && n.kind!=PlanNodeKind::Sequence)valid=false;
        if(!n.operation.empty())actual.emplace_back(n.operation,n.provider);
        for(const auto& c:n.children)self(self,c);
        if(n.kind==PlanNodeKind::StageLoop)actual.emplace_back("</StageLoop>",std::string{});
    };
    walk(walk,plan.root);
    if(!valid || actual!=expected)return "Invalid synchronous temporal plan: phase order, participant, StageLoop recipe or single terminal Commit differs from the compiled contract.";
    return {};
}

} // namespace SF::System
