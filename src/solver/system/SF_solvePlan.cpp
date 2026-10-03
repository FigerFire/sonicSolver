#include "core/system/SF_operationIds.h"
/// @file SF_solvePlan.cpp
/// @brief Structured solve-plan 编译及 transitional legacy block lowering。

#include "SF_solvePlan.h"
#include "SF_methodObjects.h"
#include "core/system/SF_formula.h"
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
        const LegacyExecutionPolicy& policy,
        OpId operation = {}) {
    SolvePlanNode result;
    result.kind = kind;
    result.id = std::move(id);
    result.name = std::move(name);
    result.equations = policy.equations;
    result.constraints = policy.constraints;
    result.operation = operation;
    result.legacyAdapter = true;
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
    result.legacyAdapter = true;
    return result;
}

LegacyExecutionPolicy legacyTimeBlock(
        const ExecutableEquationSystem& system,
        const CompiledTimeRecipe& recipe,
        const CompiledExecutionProgram& compiledExecution) {
    LegacyExecutionPolicy result;
    result.id = "TIME_RECIPE";
    result.name = std::string(FDM::toString(recipe.id()))+" time recipe";
    result.strategyName = FDM::toString(recipe.id());
    result.strategyKind = FDM::SolveStrategyKind::ExplicitTimeIntegration;
    if (!compiledExecution.steps.empty()) {
        for (const auto& step:compiledExecution.steps) {
            for (const auto& input:step.sourceMathInputs)
                if (system.registry.contains(input))
                    result.equations.push_back(input);
            for (const auto& call:step.calls)
                if (system.registry.contains(call.equation))
                    result.equations.push_back(call.equation);
            for (const auto& symbol:std::vector<std::string>{step.source.target.symbol})
                if (std::any_of(system.state.symbols().begin(),system.state.symbols().end(),
                    [&](const StateSymbol& unknown) {
                        return unknown.id==symbol;
                    }))
                    result.unknowns.push_back(symbol);
        }
        return result;
    }
    for (const auto& equation : system.legacyEquations) {
        const auto& definition = system.legacyDefinitions.at(equation.id);
        const auto hasTransient = [](const SF::Equation::Expression& expression) {
            return std::any_of(
                expression.terms.begin(),expression.terms.end(),
                [](const SF::Equation::Term& term) {
                    return term.kind == SF::Equation::TermKind::Transient;
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
        const ExecutableEquationSystem& system,
        const CompiledExecutionProgram& compiledExecution) {
    SolvePlanNode result;
    result.id = node.id;
    result.name = node.name;
    if (node.kind == PlanFragmentNode::Kind::Leaf) {
        if (!node.step.empty()) {
            if (!node.operation.empty() || !node.missingOperation.empty()
                || !node.equationCalls.empty())
                throw std::runtime_error("MethodStepRef duplicates its method authority.");
            const auto matches=std::count_if(compiledExecution.steps.begin(),
                compiledExecution.steps.end(),[&](const CompiledEquationCall& step) {
                    return step.source.equation==node.step;
                });
            if (matches!=1)
                throw std::runtime_error("MethodStepRef '"+node.step
                    +"' requires exactly one compiled EquationMethod step.");
            const auto found=std::find_if(compiledExecution.steps.begin(),
                compiledExecution.steps.end(),[&](const CompiledEquationCall& step) {
                    return step.source.equation==node.step;
                });
            if (found->fragment.children.empty())
                throw std::runtime_error("MethodStepRef '"+node.step
                    +"' has no compiled numerical fragment.");
            const auto validateOperations=[&](const auto& self,
                const SolvePlanNode& methodNode)->void {
                if (methodNode.children.empty()
                    && std::none_of(system.operations.begin(),
                        system.operations.end(),
                        [&](const ExecutableOperation& operation) {
                            return operation.operation==methodNode.operation;
                        }))
                    throw std::runtime_error("MethodStepRef '"+node.step
                        +"' references undeclared operation '"
                        +methodNode.operation+"'.");
                for (const auto& child:methodNode.children)
                    self(self,child);
            };
            validateOperations(validateOperations,found->fragment);
            return found->fragment;
        }
        result.kind = node.leafKind;
        result.legacyAdapter = true;
        result.equationCalls = node.equationCalls;
        for (const CompiledMathRef& call:result.equationCalls) {
            (void)system.registry.at(call.equation);
            if (call.target.empty())
                throw std::runtime_error("Plan FormulaCall has no target: '"
                                         +call.equation+"'.");
        }
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
    result.terminationSignal = node.terminationSignal;
    result.minimumIterations = node.minimumIterations;
    if (!result.terminationSignal.empty()
        && (result.kind != PlanNodeKind::Loop
            || result.minimumIterations < 1
            || result.minimumIterations > result.repetitions)) {
        throw std::runtime_error(
            "Plan fragment '"+node.id
            +"' has an invalid loop-termination contract.");
    }
    for (const PlanFragmentNode& child : node.children) {
        result.children.push_back(lowerFragmentNode(child,system,compiledExecution));
    }
    return result;
}

SolvePlanNode compileLegacyPolicy(
        const LegacyExecutionPolicy& policy) {
    if (policy.leafOperation.empty()) {
        throw std::runtime_error(
            "Execution policy '"+policy.id
            +"' has no declared operation for generic leaf lowering.");
    }
    return leaf(policy.leafKind,policy.id,policy.name,policy,
                policy.leafOperation);
}

SolvePlanNode compileLegacyExplicitStep(
        const ExecutableEquationSystem& system,
        const LegacyExecutionPolicy& time,
        const std::vector<LegacyExecutionPolicy>& policies,
        const CompiledTimeRecipe& timeRecipe,
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
        "Explicit.prepare","prepare physical step",OpIds::FlowStepPrepare));
    root.children.push_back(op(
        "Explicit.dt","compute stable time step",OpIds::FlowDtCompute));
    root.children.push_back(op(
        "Explicit.begin","begin explicit step",OpIds::FlowStepBegin));

    for (const auto& policy : policies) {
        if (policy.kind == LegacyExecutionPolicyKind::BoundaryClosure) {
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
    auto stage=operationLeaf(
        PlanNodeKind::Update,"Explicit.stage.execute","execute explicit stage",
        {},OpIds::ExplicitStageExecute);
    for (const auto& id:time.equations) {
        const auto descriptor=std::find_if(system.legacyEquations.begin(),system.legacyEquations.end(),
            [&](const EquationDescriptor& equation) { return equation.id==id; });
        if (descriptor==system.legacyEquations.end())
            throw std::runtime_error("Explicit FormulaCall references missing equation '"
                                     +id+"'.");
        const auto& formula=system.registry.at(id);
        const auto hasDdt=[&](const auto& self,const FormulaExpr& expression,
                              const std::string& target) -> bool {
            if (expression.kind==FormulaExpr::Kind::Operator
                && expression.name=="ddt" && expression.arguments.size()==1
                && expression.arguments[0].kind==FormulaExpr::Kind::Symbol
                && expression.arguments[0].name==target) return true;
            return std::any_of(expression.arguments.begin(),expression.arguments.end(),
                [&](const FormulaExpr& child) { return self(self,child,target); });
        };
        for (const auto& target:descriptor->solvedUnknowns) {
            if (!hasDdt(hasDdt,formula.lhs,target))
                throw std::runtime_error("Explicit FormulaCall '"+id+"' targets '"
                    +target+"' without ddt(target) on its left side.");
            stage.equationCalls.push_back({id,target});
        }
    }
    if (stage.equationCalls.empty())
        throw std::runtime_error("Explicit stage has no FormulaCalls.");
    stages.children.push_back(std::move(stage));
    root.children.push_back(std::move(stages));

    for (const auto& policy : policies) {
        if (!policy.leafOperation.empty()) {
            root.children.push_back(compileLegacyPolicy(policy));
            consumed.insert(policy.id);
        }
    }
    root.children.push_back(op(
        "Explicit.commit","commit physical state",OpIds::FlowStepCommit));
    root.children.push_back(op(
        "Explicit.time.commit","commit time and step",OpIds::TimeCommit));
    return root;
}

SolvePlanNode compileExplicitStep(
        const ExecutableEquationSystem& system,
        const LegacyExecutionPolicy& time,
        const std::vector<LegacyExecutionPolicy>& policies,
        const CompiledTimeRecipe& timeRecipe,
        const CompiledExecutionProgram& compiledExecution,
        std::set<std::string>& consumed) {
    if (compiledExecution.hasTemporalRoot) {
        auto root=compiledExecution.temporalRoot;
        if (policies.empty()) return root;
        const auto hasCommit=[](const auto& self,const SolvePlanNode& node)->bool {
            if (node.kind==PlanNodeKind::Commit) return true;
            return std::any_of(node.children.begin(),node.children.end(),[&](const auto& child) { return self(self,child); });
        };
        const auto commit=std::find_if(root.children.begin(),root.children.end(),
            [&](const SolvePlanNode& node) { return hasCommit(hasCommit,node); });
        if (commit==root.children.end())
            throw std::runtime_error("Compiled temporal fragment has no state commit.");
        auto position=commit;
        for (const auto& policy:policies) {
            if (policy.kind==LegacyExecutionPolicyKind::BoundaryClosure)
                consumed.insert(policy.id);
            if (policy.leafOperation.empty()) continue;
            position=root.children.insert(position,compileLegacyPolicy(policy));
            ++position;
            consumed.insert(policy.id);
        }
        return root;
    }
    if (!compiledExecution.steps.empty())
        throw std::runtime_error("Migrated HOW has no compiled TemporalMethod fragment; "
                                 "legacy equation-target fallback is forbidden.");
    return compileLegacyExplicitStep(system,time,policies,timeRecipe,consumed);
}

void requireAllPoliciesConsumed(
        const std::vector<LegacyExecutionPolicy>& policies,
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
        const std::vector<LegacyExecutionPolicy>& policies,
        const CompiledTimeRecipe& timeRecipe,
        const std::vector<LegacyPlanFragment>& fragments,
        const CompiledExecutionProgram& compiledExecution) {
    CompiledSolvePlan result;
    result.root.kind = PlanNodeKind::Sequence;
    result.root.id = "PLAN_ROOT";
    result.root.name = "time-step sequence";

    auto ordered = policies;
    std::stable_sort(ordered.begin(),ordered.end(),
        [](const LegacyExecutionPolicy& left, const LegacyExecutionPolicy& right) {
            return left.priority < right.priority;
        });
    for (const auto& policy : ordered) {
        if ((policy.kind==LegacyExecutionPolicyKind::SegregatedPressureCorrection
                && policy.strategyKind!=FDM::SolveStrategyKind::PressureCorrection)
            || (policy.kind==LegacyExecutionPolicyKind::PressureVelocityFixedPoint
                && policy.strategyKind!=FDM::SolveStrategyKind::PressureVelocityCoupling)) {
            throw std::runtime_error(
                "Execution policy '"+policy.id
                +"': legacy strategy metadata disagrees with typed solve policy.");
        }
        if (policy.strategyKind == FDM::SolveStrategyKind::AlgebraicUpdate) {
            throw std::runtime_error(
                "Execution policy '"+policy.id+"' has no typed strategy.");
        }
        result.blocks.push_back({
            policy.id,policy.name,policy.strategyName,policy.equations,
            policy.constraints,policy.unknowns,policy.strategyKind,policy.kind});
    }

    const auto activeFragment = std::find_if(
        fragments.begin(),fragments.end(),[](const LegacyPlanFragment& fragment) {
            return !fragment.nodes.empty();
        });
    const LegacyExecutionPolicy resolvedTime = legacyTimeBlock(
        system,timeRecipe,compiledExecution);
    if (activeFragment == fragments.end()
        || activeFragment->includeTimePolicy) {
        result.blocks.push_back({
            resolvedTime.id,resolvedTime.name,resolvedTime.strategyName,
            resolvedTime.equations,resolvedTime.constraints,resolvedTime.unknowns,
            resolvedTime.strategyKind,resolvedTime.kind});
    }
    if (activeFragment != fragments.end()) {
        // 顺序与重复来自 coupling preset 的片段；planner 只解析与组装。
        SolvePlanNode root;
        root.kind = PlanNodeKind::Sequence;
        root.id = activeFragment->rootId;
        root.name = activeFragment->name;
        const auto appendStage = [&](OperationStage stage, const char* id) {
            PlanFragmentNode node;
            node.kind = PlanFragmentNode::Kind::Leaf;
            node.stage = stage;
            node.id = id;
            node.name = id;
            root.children.push_back(lowerFragmentNode(node,system,compiledExecution));
        };
        if (activeFragment->physicalStepBody)
            appendStage(OperationStage::Prepare,"PhysicalStep.prepare");
        for (const PlanFragmentNode& node : activeFragment->nodes) {
            root.children.push_back(lowerFragmentNode(node,system,compiledExecution));
        }
        if (activeFragment->physicalStepBody) {
            appendStage(OperationStage::StepCommit,"PhysicalStep.state.commit");
            root.children.push_back(operationLeaf(PlanNodeKind::Commit,
                "PhysicalStep.time.commit","commit physical clock",{},OpIds::TimeCommit));
        }
        result.root = std::move(root);
        requireAllPoliciesConsumed(ordered,
            std::set<std::string>(activeFragment->consumedPolicies.begin(),
                                  activeFragment->consumedPolicies.end()));
        return result;
    }

    if (!compiledExecution.steps.empty() && !compiledExecution.hasTemporalRoot) {
        if (compiledExecution.hasTemporalRoot)
            throw std::runtime_error("Structured HOW and temporal root both own the physical step.");
        result.root=compileMethodProgram(compiledExecution);
        requireAllPoliciesConsumed(ordered,
            std::set<std::string>(std::vector<std::string>{}.begin(),
                                  std::vector<std::string>{}.end()));
        return result;
    }

    std::set<std::string> consumed;
    if (timeRecipe.topology() == FDM::TimeTopology::ExplicitStages) {
        const bool hasExplicitStage = std::any_of(
            system.operations.begin(),system.operations.end(),
            [](const ExecutableOperation& operation) {
                return operation.operation == OpIds::ExplicitStageExecute;
            });
        if (!hasExplicitStage && !ordered.empty()) {
            std::string ids;
            for (const auto& policy : ordered) {
                if (!ids.empty()) ids += ", ";
                ids += policy.id;
            }
            throw std::runtime_error(
                "SolvePlanner has no LegacyPlanFragment or explicit stage operation "
                "for active execution policies: ["+ids+"].");
        }
        result.root = compileExplicitStep(
            system,resolvedTime,ordered,timeRecipe,compiledExecution,consumed);
        requireAllPoliciesConsumed(ordered,consumed);
        return result;
    }
    for (const auto& policy : ordered) {
        result.root.children.push_back(
            compileLegacyPolicy(policy));
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
