/// @file SF_solvePlan.cpp
/// @brief Structured solve-plan 编译及 transitional legacy block lowering。

#include "SF_solvePlan.h"
#include "SF_config.h"
#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>

namespace SF::System {
namespace {

SolvePlanNode leaf(
        PlanNodeKind kind,
        std::string id,
        std::string name,
        const ExecutionPolicy& policy,
        OpId operation = {}) {
    SolvePlanNode result;
    result.kind = kind;
    result.id = std::move(id);
    result.name = std::move(name);
    result.equations = policy.equations;
    result.constraints = policy.constraints;
    result.operation = operation;
    return result;
}

SolvePlanNode operationLeaf(
        PlanNodeKind kind,
        std::string id,
        std::string name,
        std::string equation,
        OpId operation) {
    SolvePlanNode result;
    result.kind = kind;
    result.id = std::move(id);
    result.name = std::move(name);
    if (!equation.empty()) result.equations.push_back(std::move(equation));
    result.operation = operation;
    return result;
}

ExecutionPolicy timePolicy(
        const ExecutableEquationSystem& system,
        const FDM::TimeRecipe& recipe) {
    ExecutionPolicy result;
    result.id = "TIME_RECIPE";
    result.name = std::string(FDM::toString(recipe.id()))+" time recipe";
    result.strategyName = FDM::toString(recipe.id());
    result.strategyKind = FDM::SolveStrategyKind::ExplicitTimeIntegration;
    for (const auto& equation : system.equations) {
        const auto& definition = system.equationDefinitions.at(equation.id);
        const auto hasTransient = [](const Equation::Expression& expression) {
            return std::any_of(
                expression.terms.begin(),expression.terms.end(),
                [](const Equation::Term& term) {
                    return term.kind == Equation::TermKind::Transient;
                });
        };
        if (equation.category != EquationCategory::PhysicalEquation
            && !hasTransient(definition.left) && !hasTransient(definition.right)) {
            continue;
        }
        result.equations.push_back(equation.id);
        result.unknowns.insert(result.unknowns.end(),
                               equation.solvedUnknowns.begin(),
                               equation.solvedUnknowns.end());
    }
    return result;
}

// Fragment lowering：planner 只把 stage 解析成 OpId、并按片段给出的顺序与
// 重复次数组装 control flow。它不知道 具体数学，也不知道 provider 名单。
SolvePlanNode lowerFragmentNode(
        const PlanFragmentNode& node,
        const ExecutableEquationSystem& system) {
    SolvePlanNode result;
    result.id = node.id;
    result.name = node.name;
    if (node.kind == PlanFragmentNode::Kind::Leaf) {
        result.kind = node.leafKind;
        result.unsupportedReason = node.unsupportedReason;
        if (!node.missingOperation.empty()) {
            // 具名缺失标记：该 schedule step 没有 formulation 声明。plan 仍然
            // 真实描述控制流，provider resolver 会把它报成 Unsupported。
            result.operation = node.missingOperation;
        } else if (!node.operation.empty()) {
            const auto found = std::find_if(
                system.operations.begin(),system.operations.end(),
                [&](const ExecutableOperation& item) {
                    return item.operation == node.operation;
                });
            if (found == system.operations.end()) {
                throw std::runtime_error(
                    "Plan fragment leaf '"+node.id
                    +"' references undeclared operation '"
                    +node.operation+"'.");
            }
            result.operation = found->operation;
        } else {
            const ExecutableOperation* operation =
                findExecutableOperation(system,node.stage);
            if (!operation) {
                throw std::runtime_error(
                    "Plan fragment leaf '"+node.id
                    +"' references stage '"+toString(node.stage)
                    +"' which the formulation never declared.");
            }
            result.operation = operation->operation;
        }
        if (!node.equation.empty()) result.equations.push_back(node.equation);
        return result;
    }
    result.kind = node.nodeKind;
    result.repetitions = node.repetitions;
    for (const PlanFragmentNode& child : node.children) {
        result.children.push_back(lowerFragmentNode(child,system));
    }
    return result;
}

SolvePlanNode compilePolicy(
        const ExecutionPolicy& policy) {
    if (policy.leafOperation.empty()) {
        throw std::runtime_error(
            "Execution policy '"+policy.id
            +"' has no declared operation for generic leaf lowering.");
    }
    return leaf(policy.leafKind,policy.id,policy.name,policy,
                policy.leafOperation);
}

SolvePlanNode compileExplicitStep(
        const ExecutionPolicy& time,
        const std::vector<ExecutionPolicy>& policies,
        const FDM::TimeRecipe& timeRecipe,
        std::set<std::string>& consumed) {
    const auto op = [](std::string id, std::string name, OpId operation) {
        return operationLeaf(PlanNodeKind::Update,std::move(id),
                             std::move(name),{},std::move(operation));
    };
    SolvePlanNode root;
    root.kind = PlanNodeKind::Sequence;
    root.id = "Explicit.step";
    root.name = "explicit time step";
    root.children.push_back(op(
        "Explicit.prepare","prepare physical step","flow.step.prepare"));
    root.children.push_back(op(
        "Explicit.dt","compute stable time step","flow.dt.compute"));
    root.children.push_back(op(
        "Explicit.begin","begin explicit step","flow.step.begin"));

    for (const auto& policy : policies) {
        if (policy.kind == ExecutionPolicyKind::BoundaryClosure) {
            // Boundary/interface closures are consumed by the fused explicit
            // RHS stage. ConservativeRHS preserves the validated
            // boundary/halo/IBM ordering.
            consumed.insert(policy.id);
        }
    }
    if (time.equations.empty()) {
        throw std::runtime_error(
            "Explicit time recipe has no transient equation to execute.");
    }
    SolvePlanNode stages = leaf(
        PlanNodeKind::StageLoop,"Explicit.stages",
        time.name+" fused explicit stages",time);
    stages.repetitions = timeRecipe.stageCount();
    stages.children.push_back(operationLeaf(
        PlanNodeKind::Update,"Explicit.stage.execute","execute explicit stage",
        {},"explicit.stage.execute"));
    root.children.push_back(std::move(stages));

    for (const auto& policy : policies) {
        if (!policy.leafOperation.empty()) {
            root.children.push_back(compilePolicy(policy));
            consumed.insert(policy.id);
        }
    }
    root.children.push_back(op(
        "Explicit.commit","commit physical state","flow.step.commit"));
    root.children.push_back(op(
        "Explicit.time.commit","commit time and step","time.commit"));
    return root;
}

void requireAllPoliciesConsumed(
        const std::vector<ExecutionPolicy>& policies,
        const std::set<std::string>& consumed) {
    for (const auto& policy : policies) {
        if (consumed.find(policy.id) == consumed.end()) {
            throw std::runtime_error(
                "SolvePlanner cannot lower active execution policy '"
                +policy.id+"' together with the selected structured plan.");
        }
    }
}

} // namespace

CompiledSolvePlan SolvePlanner::compile(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies,
        const FDM::TimeRecipe& timeRecipe,
        const std::vector<PlanFragment>& fragments) {
    CompiledSolvePlan result;
    result.root.kind = PlanNodeKind::Sequence;
    result.root.id = "PLAN_ROOT";
    result.root.name = "time-step sequence";

    auto ordered = policies;
    std::stable_sort(ordered.begin(),ordered.end(),
        [](const ExecutionPolicy& left, const ExecutionPolicy& right) {
            return left.priority < right.priority;
        });
    for (const auto& policy : ordered) {
        if (policy.strategyKind == FDM::SolveStrategyKind::AlgebraicUpdate) {
            throw std::runtime_error(
                "Execution policy '"+policy.id+"' has no typed strategy.");
        }
        result.blocks.push_back({
            policy.id,policy.name,policy.strategyName,policy.equations,
            policy.constraints,policy.unknowns,policy.strategyKind});
    }

    const auto activeFragment = std::find_if(
        fragments.begin(),fragments.end(),[](const PlanFragment& fragment) {
            return !fragment.nodes.empty();
        });
    const ExecutionPolicy resolvedTime = timePolicy(system,timeRecipe);
    if (activeFragment == fragments.end()
        || activeFragment->includeTimePolicy) {
        result.blocks.push_back({
            resolvedTime.id,resolvedTime.name,resolvedTime.strategyName,
            resolvedTime.equations,resolvedTime.constraints,resolvedTime.unknowns,
            resolvedTime.strategyKind});
    }
    if (activeFragment != fragments.end()) {
        // 顺序与重复来自 coupling preset 的片段；planner 只解析与组装。
        SolvePlanNode root;
        root.kind = PlanNodeKind::Sequence;
        root.id = activeFragment->rootId;
        root.name = activeFragment->name;
        for (const PlanFragmentNode& node : activeFragment->nodes) {
            root.children.push_back(lowerFragmentNode(node,system));
        }
        result.root = std::move(root);
        requireAllPoliciesConsumed(ordered,
            std::set<std::string>(activeFragment->consumedPolicies.begin(),
                                  activeFragment->consumedPolicies.end()));
        return result;
    }

    std::set<std::string> consumed;
    if (timeRecipe.topology() == FDM::TimeTopology::ExplicitStages) {
        const bool hasExplicitStage = std::any_of(
            system.operations.begin(),system.operations.end(),
            [](const ExecutableOperation& operation) {
                return operation.operation == "explicit.stage.execute";
            });
        if (!hasExplicitStage && !ordered.empty()) {
            std::string ids;
            for (const auto& policy : ordered) {
                if (!ids.empty()) ids += ", ";
                ids += policy.id;
            }
            throw std::runtime_error(
                "SolvePlanner has no PlanFragment or explicit stage operation "
                "for active execution policies: ["+ids+"].");
        }
        result.root = compileExplicitStep(
            resolvedTime,ordered,timeRecipe,consumed);
        requireAllPoliciesConsumed(ordered,consumed);
        return result;
    }
    for (const auto& policy : ordered) {
        result.root.children.push_back(
            compilePolicy(policy));
        consumed.insert(policy.id);
    }

    requireAllPoliciesConsumed(ordered,consumed);

    return result;
}

std::vector<OpId> SolvePlanner::requiredOperations(const CompiledSolvePlan& plan) {
    std::vector<OpId> result;
    const auto visit = [&](const auto& self, const SolvePlanNode& node) -> void {
        if (!node.operation.empty()
            && std::find(result.begin(),result.end(),node.operation) == result.end()) {
            result.push_back(node.operation);
        }
        for (const auto& child : node.children) self(self,child);
    };
    visit(visit,plan.root);
    return result;
}

} // namespace SF::System
