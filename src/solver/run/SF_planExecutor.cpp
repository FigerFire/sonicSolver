/// @file SF_planExecutor.cpp
/// @brief 解释 CompiledSolvePlan 的结构化 control flow。
///
/// Data flow:
///   CompiledSolvePlan -> Sequence / Loop / StageLoop
///       -> ExecutionContext -> OpRegistry callback
///
/// 本文件不解释 physics、time scheme 或 numerical formula。

#include "solver/run/SF_planExecutor.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <map>

namespace SF::Run {
namespace {
using LoopFrame = std::pair<std::string,std::pair<int,int>>;

bool traceEnabled() {
    const char* value = std::getenv("SF_PLAN_TRACE");
    return value && *value && std::string(value) != "0";
}

void traceOperation(
        const System::OpId& operation,
        const PlanTraceContext* context,
        const std::vector<LoopFrame>& loops,
        const ExecutionContext& execution) {
    if (!traceEnabled()) return;
    std::cout << std::setprecision(17) << "[SF PLAN]";
    if (context) {
        std::cout << " step=" << context->step
                  << " time=" << context->time
                  << " dt=" << (context->dt ? *context->dt : 0.0);
    }
    for (const auto& loop : loops) {
        std::cout << " " << loop.first << "="
                  << loop.second.first << "/" << loop.second.second;
    }
    if (execution.stageIndex >= 0) {
        std::cout << " stage=" << execution.stageIndex + 1
                  << "/" << execution.stageCount;
    }
    std::cout << " op=" << operation << '\n';
}

void executeNode(
        const System::SolvePlanNode& node,
        const OpRegistry& operations,
        const PlanTraceContext* trace,
        std::vector<LoopFrame>& loops,
        ExecutionContext& execution,
        std::map<System::OpId,std::size_t>& revisions) {
    if (node.repetitions <= 0) throw std::runtime_error("Compiled solve-plan node '"+node.id+"' has a non-positive repetition count.");
    if (node.kind == System::PlanNodeKind::Sequence || node.kind == System::PlanNodeKind::Loop) {
        if (node.kind == System::PlanNodeKind::Loop && node.children.empty()) throw std::runtime_error("Compiled solve-plan loop '"+node.id+"' is empty.");
        if (!node.terminationSignal.empty()) {
            if (node.kind != System::PlanNodeKind::Loop
                || node.minimumIterations < 1
                || node.minimumIterations > node.repetitions)
                throw std::runtime_error("Invalid loop termination contract for '"+node.id+"'.");
            execution.signals->reset(node.terminationSignal);
        }
        struct Metric { std::vector<double> previous; double initial=0; bool initialized=false; };
        std::vector<Metric> metrics(node.residuals.size());
        for(std::size_t i=0;i<node.residuals.size();++i)
            if(node.residuals[i].kind==System::ResidualKind::StateChange) {
                metrics[i].previous=operations.sampleResidual(node.residuals[i].id,node.residuals[i].provider,execution);
                if(metrics[i].previous.empty() || std::any_of(metrics[i].previous.begin(),metrics[i].previous.end(),[](double value){return !std::isfinite(value);}))
                    throw std::runtime_error("Invalid initial residual reference: "+node.residuals[i].id);
            }
        int completed=0;
        bool converged=node.residuals.empty();
        for (int repeat = 0; repeat < node.repetitions; ++repeat) {
            const auto before=revisions;
            if(!node.terminationSignal.empty()) execution.signals->reset(node.terminationSignal);
            if (node.kind == System::PlanNodeKind::Loop) {
                loops.push_back({node.id,{repeat+1,node.repetitions}});
            }
            for (const auto& child : node.children) {
                executeNode(child,operations,trace,loops,execution,revisions);
            }
            if (node.kind == System::PlanNodeKind::Loop) loops.pop_back();
            completed=repeat+1;
            if(!node.residuals.empty()) {
                converged=node.residualAll;
                for(std::size_t i=0;i<node.residuals.size();++i) {
                    const auto& definition=node.residuals[i];auto& metric=metrics[i];
                    const auto old=before.find(definition.producerOperation);
                    if(revisions[definition.producerOperation]==(old==before.end()?0:old->second))
                        throw std::runtime_error("Residual '"+definition.id+"' was not updated in Loop '"+node.id+"'.");
                    auto sample=operations.sampleResidual(definition.id,definition.provider,execution);
                    if(sample.empty() || (definition.kind==System::ResidualKind::StateChange && sample.size()!=metric.previous.size()))
                        throw std::runtime_error("Residual sampler changed layout or returned no values: "+definition.id);
                    double norm=0;
                    for(std::size_t j=0;j<sample.size();++j) {
                        const double value=definition.kind==System::ResidualKind::StateChange?sample[j]-metric.previous[j]:sample[j];
                        if(!std::isfinite(value) || !std::isfinite(sample[j]))throw std::runtime_error("Non-finite residual: "+definition.id);
                        if(definition.norm==System::ResidualNorm::Maximum) norm=std::max(norm,std::abs(value));
                        else norm=std::hypot(norm,value);
                    }
                    if(!metric.initialized){metric.initial=norm;metric.initialized=true;}
                    if(definition.kind==System::ResidualKind::StateChange)metric.previous=std::move(sample);
                    const double scale=definition.normalization==System::ResidualNormalization::Initial?metric.initial:definition.scale;
                    const bool relative=definition.normalization!=System::ResidualNormalization::None
                        && (scale==0?norm==0:norm/scale<=definition.relativeTolerance);
                    const bool accepted=norm<=definition.absoluteTolerance || relative;
                    converged=node.residualAll ? converged && accepted : converged || accepted;
                    if(traceEnabled())std::cout<<std::setprecision(17)<<"[SF RES] loop="<<node.id<<" iteration="<<repeat+1
                        <<" metric="<<definition.id<<" norm="<<norm<<" reference="<<scale<<" accepted="<<accepted<<'\n';
                }
                if(repeat+1>=node.minimumIterations && converged)break;
            }
            if (!node.terminationSignal.empty()
                && repeat+1 >= node.minimumIterations
                && execution.signals->value(node.terminationSignal)) break;
        }
        if(!node.residuals.empty())std::cout<<"[SF RES EXIT] loop="<<node.id<<" iterations="<<completed
            <<" reason="<<(converged?"converged":node.failOnMaximum?"maximum-failed":"maximum-continue-explicit")<<'\n';
        if(!node.residuals.empty() && !converged && node.failOnMaximum)
            throw std::runtime_error("Loop '"+node.id+"' reached maximum iterations without residual convergence.");
        return;
    }
    if (node.kind == System::PlanNodeKind::StageLoop) {
        if (node.children.empty()) {
            throw std::runtime_error(
                "Compiled solve-plan StageLoop '"+node.id+"' is empty.");
        }
        const ExecutionContext previous = execution;
        execution.stageCount = node.repetitions;
        for (int stage = 0; stage < node.repetitions; ++stage) {
            execution.stageIndex = stage;
            for (const auto& child : node.children) {
                executeNode(child,operations,trace,loops,execution,revisions);
            }
        }
        execution = previous;
        return;
    }
    if (node.kind == System::PlanNodeKind::Subcycle)
        throw std::runtime_error("PlanExecutor does not implement control node '"+std::string(System::toString(node.kind))+"'.");
    if (!node.children.empty()) throw std::runtime_error("Leaf solve-plan node '"+node.id+"' has child nodes.");
    if (node.operation.empty()) throw std::runtime_error("Executable solve-plan leaf '"+node.id+"' has no operation ID.");
    traceOperation(node.operation,trace,loops,execution);
    const auto* previousTarget=execution.target;
    execution.target=node.target ? &*node.target : nullptr;
    operations.invoke(node.operation,execution,node.provider);
    ++revisions[node.operation];
    execution.target=previousTarget;
}

void validateNode(
        const System::SolvePlanNode& node,
        const OpRegistry& operations) {
    if(!node.residuals.empty()) {
        if(node.kind!=System::PlanNodeKind::Loop || !node.terminationSignal.empty()
            || node.minimumIterations<1 || node.minimumIterations>node.repetitions)
            throw std::runtime_error("Typed residual control requires one ordinary Loop without a second termination signal: "+node.id);
        for(const auto& residual:node.residuals)
            if(!operations.containsResidual(residual.id,residual.provider))
                throw std::runtime_error("Residual '"+residual.id+"' has no sampler from its frozen numerical provider.");
    }
    if (node.children.empty()) {
        if (node.operation.empty()) {
            throw std::runtime_error(
                "CompiledSolvePlan leaf '"+node.id
                +"' has no runtime operation ID.");
        }
        if (!operations.contains(node.operation)) {
            throw std::runtime_error(
                "CompiledSolvePlan requires runtime operation '"
                +node.operation+"', but no numerical provider is bound.");
        }
        if (!operations.contains(node.operation,node.provider))
            throw std::runtime_error("CompiledSolvePlan operation '"+node.operation
                +"' requires frozen provider '"+node.provider+"', but its runtime owner differs.");
        return;
    }
    for (const auto& child : node.children) {
        validateNode(child,operations);
    }
}

}
void LoopSignalState::reset(const System::LoopSignalId& id) {
    if (id.empty()) throw std::runtime_error("Loop signal ID is empty.");
    values_.erase(id);
}
void LoopSignalState::publish(const System::LoopSignalId& id, bool value) {
    if (id.empty()) throw std::runtime_error("Loop signal ID is empty.");
    values_[id]=value;
}
bool LoopSignalState::value(const System::LoopSignalId& id) const {
    const auto found=values_.find(id);
    if (found==values_.end())
        throw std::runtime_error("Loop termination signal '"+id
            +"' was not published by its operation.");
    return found->second;
}
void OpRegistry::bind(System::OpId id, Operation operation) {
    bind(std::move(id),{},std::move(operation));
}
void OpRegistry::bind(System::OpId id, ContextOperation operation) {
    bind(std::move(id),{},std::move(operation));
}
void OpRegistry::bind(System::OpId id,std::string provider,Operation operation) {
    if (id.empty() || !operation) throw std::runtime_error("Cannot bind an empty plan operation.");
    if (contains(id)) throw std::runtime_error("Duplicate plan operation binding '"+id+"'.");
    bind(std::move(id),std::move(provider),
         [operation=std::move(operation)](const ExecutionContext&) {
             operation();
         });
}
void OpRegistry::bind(System::OpId id,std::string provider,ContextOperation operation) {
    if (id.empty() || !operation) throw std::runtime_error("Cannot bind an empty plan operation.");
    if (contains(id)) throw std::runtime_error("Duplicate plan operation binding '"+id+"'.");
    entries_.push_back({std::move(id),std::move(provider),std::move(operation)});
}
bool OpRegistry::contains(const System::OpId& id) const { return std::any_of(entries_.begin(),entries_.end(),[&](const Entry& entry) { return entry.id == id; }); }
bool OpRegistry::contains(const System::OpId& id,std::string_view provider) const {
    return std::any_of(entries_.begin(),entries_.end(),[&](const Entry& entry) {
        return entry.id==id && (provider.empty() || entry.provider==provider);
    });
}
void OpRegistry::bindResidual(std::string id,std::string provider,ResidualOperation operation) {
    if(id.empty() || !operation || std::any_of(residuals_.begin(),residuals_.end(),[&](const auto& r){return r.id==id;}))
        throw std::runtime_error("Empty or duplicate residual sampler: "+id);
    residuals_.push_back({std::move(id),std::move(provider),std::move(operation)});
}
bool OpRegistry::containsResidual(std::string_view id,std::string_view provider) const {
    return std::any_of(residuals_.begin(),residuals_.end(),[&](const auto& r){return r.id==id && r.provider==provider;});
}
std::vector<double> OpRegistry::sampleResidual(std::string_view id,std::string_view provider,const ExecutionContext& context) const {
    const auto at=std::find_if(residuals_.begin(),residuals_.end(),[&](const auto& r){return r.id==id && r.provider==provider;});
    if(at==residuals_.end())throw std::runtime_error("Missing frozen residual sampler: "+std::string(id));
    return at->operation(context);
}
void OpRegistry::retain(const std::vector<System::OpId>& assigned) {
    entries_.erase(std::remove_if(entries_.begin(),entries_.end(),
        [&](const Entry& entry) {
            return std::find(assigned.begin(),assigned.end(),entry.id)
                == assigned.end();
        }),entries_.end());
}
void OpRegistry::invoke(
        const System::OpId& id,
        const ExecutionContext& context,std::string_view provider) const {
    if (id.empty()) throw std::runtime_error("Executable solve-plan leaf has no operation binding.");
    const auto found = std::find_if(entries_.begin(),entries_.end(),[&](const Entry& entry) { return entry.id == id; });
    if (found == entries_.end()) throw std::runtime_error("No runtime implementation is bound for plan operation '"+id+"'.");
    if (!provider.empty() && found->provider!=provider)
        throw std::runtime_error("Runtime provider differs from frozen operation owner: "+id);
    found->operation(context);
}
void PlanExecutor::execute(
        const System::CompiledSolvePlan& plan,
        const OpRegistry& operations,
        const PlanTraceContext* trace) {
    validateBindings(plan,operations);
    std::vector<LoopFrame> loops;
    ExecutionContext execution;
    LoopSignalState signals;
    execution.signals=&signals;
    std::map<System::OpId,std::size_t> revisions;
    executeNode(plan.root,operations,trace,loops,execution,revisions);
}

void PlanExecutor::validateBindings(
        const System::CompiledSolvePlan& plan,
        const OpRegistry& operations) {
    validateNode(plan.root,operations);
}

} // namespace SF::Run
