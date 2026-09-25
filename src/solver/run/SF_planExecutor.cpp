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
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

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
        ExecutionContext& execution) {
    if (node.repetitions <= 0) throw std::runtime_error("Compiled solve-plan node '"+node.id+"' has a non-positive repetition count.");
    if (node.kind == System::PlanNodeKind::Sequence || node.kind == System::PlanNodeKind::Loop) {
        if (node.kind == System::PlanNodeKind::Loop && node.children.empty()) throw std::runtime_error("Compiled solve-plan loop '"+node.id+"' is empty.");
        for (int repeat = 0; repeat < node.repetitions; ++repeat) {
            if (node.kind == System::PlanNodeKind::Loop) {
                loops.push_back({node.id,{repeat+1,node.repetitions}});
            }
            for (const auto& child : node.children) {
                executeNode(child,operations,trace,loops,execution);
            }
            if (node.kind == System::PlanNodeKind::Loop) loops.pop_back();
        }
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
                executeNode(child,operations,trace,loops,execution);
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
    operations.invoke(node.operation,execution);
}

void validateNode(
        const System::SolvePlanNode& node,
        const OpRegistry& operations) {
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
        return;
    }
    for (const auto& child : node.children) {
        validateNode(child,operations);
    }
}

}
void OpRegistry::bind(System::OpId id, Operation operation) {
    if (id.empty() || !operation) throw std::runtime_error("Cannot bind an empty plan operation.");
    if (contains(id)) throw std::runtime_error("Duplicate plan operation binding '"+id+"'.");
    bind(std::move(id),
         [operation=std::move(operation)](const ExecutionContext&) {
             operation();
         });
}
void OpRegistry::bind(System::OpId id, ContextOperation operation) {
    if (id.empty() || !operation) throw std::runtime_error("Cannot bind an empty plan operation.");
    if (contains(id)) throw std::runtime_error("Duplicate plan operation binding '"+id+"'.");
    entries_.push_back({std::move(id),std::move(operation)});
}
bool OpRegistry::contains(const System::OpId& id) const { return std::any_of(entries_.begin(),entries_.end(),[&](const Entry& entry) { return entry.id == id; }); }
void OpRegistry::retain(const std::vector<System::OpId>& assigned) {
    entries_.erase(std::remove_if(entries_.begin(),entries_.end(),
        [&](const Entry& entry) {
            return std::find(assigned.begin(),assigned.end(),entry.id)
                == assigned.end();
        }),entries_.end());
}
void OpRegistry::invoke(
        const System::OpId& id,
        const ExecutionContext& context) const {
    if (id.empty()) throw std::runtime_error("Executable solve-plan leaf has no operation binding.");
    const auto found = std::find_if(entries_.begin(),entries_.end(),[&](const Entry& entry) { return entry.id == id; });
    if (found == entries_.end()) throw std::runtime_error("No runtime implementation is bound for plan operation '"+id+"'.");
    found->operation(context);
}
void PlanExecutor::execute(
        const System::CompiledSolvePlan& plan,
        const OpRegistry& operations,
        const PlanTraceContext* trace) {
    validateBindings(plan,operations);
    std::vector<LoopFrame> loops;
    ExecutionContext execution;
    executeNode(plan.root,operations,trace,loops,execution);
}

void PlanExecutor::validateBindings(
        const System::CompiledSolvePlan& plan,
        const OpRegistry& operations) {
    validateNode(plan.root,operations);
}

} // namespace SF::Run
