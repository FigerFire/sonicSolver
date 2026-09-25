/// @file SF_providerResolver.cpp
/// @brief Operation -> provider resolution uses typed needs and realized state.

#include "SF_providerResolver.h"

#include "SF_couplingStatus.h"
#include "SF_solvePlan.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace SF::System {
namespace {

bool hasCapability(const std::vector<OperationCapability>& values,
                   OperationCapability wanted) {
    return std::find(values.begin(),values.end(),wanted) != values.end();
}

bool supports(const std::vector<OperationCapability>& offered,
              const ExecutableOperation& operation) {
    return !operation.requirements.empty()
        && std::all_of(operation.requirements.begin(),operation.requirements.end(),
            [&](OperationCapability need) { return hasCapability(offered,need); });
}

const SolvePlanNode* findLeaf(const SolvePlanNode& node, const OpId& operation) {
    if (node.operation == operation) return &node;
    for (const SolvePlanNode& child : node.children) {
        if (const auto* found = findLeaf(child,operation)) return found;
    }
    return nullptr;
}

bool conservativePressureScheduleSupported(
        const ExecutionCapabilitySignature& signature,
        const std::vector<ExecutionPolicy>& policies,
        const FDM::TimeRecipe& timeRecipe) {
    const auto pressure = std::find_if(
        policies.begin(),policies.end(),[](const ExecutionPolicy& policy) {
            return policy.id == kPressureScheduleId;
        });
    return timeRecipe.topology() == FDM::TimeTopology::ExplicitStages
        && timeRecipe.stageCount() == 1 && policies.size() == 1
        && pressure != policies.end()
        && pressure->kind == ExecutionPolicyKind::SegregatedPressureCorrection
        && pressure->strategyKind == FDM::SolveStrategyKind::PressureCorrection
        && signature.pressureConstraint && signature.conservativeState
        && signature.momentumPredictor && signature.pressureCorrection
        && !signature.constantDensity && !signature.auxiliarySchedule
        && signature.outerCorrectors == 1
        && signature.pressureCorrectors > 0
        && signature.nonOrthogonalCorrectors == 0;
}

std::string realizationName(const CompiledStateRealization& realization) {
    if (realization.phaseTransportedState) return "Eulerian shared-pressure realization";
    if (realization.pressureMultiplier && !realization.conservativeTransportedMass) {
        return "single-fluid constant-density pressure-multiplier realization";
    }
    if (realization.conservativeTransportedMass) return "conservative transported-state realization";
    return "the compiled state realization";
}

} // namespace

ExecutionCapabilitySignature compileExecutionCapabilities(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies) {
    ExecutionCapabilitySignature result;
    result.pressureConstraint = hasConstraint(system,"C_INCOMPRESSIBILITY");
    result.constantDensity = std::any_of(
        system.unknowns.begin(),system.unknowns.end(),[](const UnknownDescriptor& item) {
            return item.id == "rho" && item.storageKey == "rhoConst";
        });
    result.conservativeState = std::any_of(
        system.unknowns.begin(),system.unknowns.end(),[](const UnknownDescriptor& item) {
            return item.storageBinding == StorageBinding::PackedDistributed
                && item.storageKey == "conservative";
        });
    result.momentumPredictor = std::any_of(
        system.compiledEquations.begin(),system.compiledEquations.end(),
        [](const CompiledEquation& item) {
            return item.operatorBinding == "momentum.predictor";
        });
    result.pressureCorrection = std::any_of(
        system.compiledEquations.begin(),system.compiledEquations.end(),
        [](const CompiledEquation& item) {
            return item.operatorBinding == "pressure.correction";
        });
    result.auxiliarySchedule = std::any_of(
        system.equations.begin(),system.equations.end(),
        [](const EquationDescriptor& item) {
            return item.id != "E_MASS" && item.id != "E_MOMENTUM"
                && item.id != "E_ENERGY" && item.id != "E_MOMENTUM_PREDICTOR"
                && item.category == EquationCategory::PhysicalEquation;
        });
    const auto pressure = std::find_if(
        policies.begin(),policies.end(),[](const ExecutionPolicy& item) {
            return item.id == kPressureScheduleId;
        });
    if (pressure != policies.end()) {
        result.outerCorrectors = pressure->repeatCount;
        result.pressureCorrectors = pressure->nestedRepeatCount;
        result.nonOrthogonalCorrectors = pressure->innerRepeatCount-1;
    }
    return result;
}

std::vector<ResolvedOperationBinding> resolveOperationBindings(
        const ExecutableEquationSystem& equations,
        const CompiledStateRealization& realization,
        const CompiledNumericalSystem& numerics,
        const CompiledSolvePlan& plan,
        const std::vector<ExecutionPolicy>& policies) {
    const auto required = SolvePlanner::requiredOperations(plan);
    const auto signature = compileExecutionCapabilities(equations,policies);
    const bool pressureSchedule = conservativePressureScheduleSupported(
        signature,policies,numerics.time.recipe);
    const std::vector<OperationCapability> conservativeCapabilities{
        OperationCapability::ConservativeExplicit,
        OperationCapability::PressureSchedule,
        OperationCapability::MomentumPredictor,
        OperationCapability::PressureCorrection,
        OperationCapability::PressureLinearSolve,
        OperationCapability::PressureBoundary,
        OperationCapability::VelocityCorrection,
        OperationCapability::FluxCorrection};
    const std::vector<OperationCapability> eulerianCapabilities{
        OperationCapability::EulerianPhaseExecution,
        OperationCapability::PressureLinearSolve};
    const std::vector<OperationCapability> immersedCapabilities{
        OperationCapability::ImmersedConstraint};
    std::vector<ResolvedOperationBinding> bindings;
    bindings.reserve(required.size());
    for (const OpId& id : required) {
        ResolvedOperationBinding binding;
        binding.operation = id;
        const auto found = std::find_if(
            equations.operations.begin(),equations.operations.end(),
            [&](const ExecutableOperation& item) { return item.operation == id; });
        if (found == equations.operations.end()) {
            const auto* leaf = findLeaf(plan.root,id);
            binding.reason = leaf && !leaf->unsupportedReason.empty()
                ? leaf->unsupportedReason
                : "no ExecutableOperation declares '"+id+"'";
        } else {
            const bool needsPressureSchedule = std::any_of(
                found->requirements.begin(),found->requirements.end(),
                [](OperationCapability need) {
                    return need != OperationCapability::ConservativeExplicit
                        && need != OperationCapability::EulerianPhaseExecution
                        && need != OperationCapability::ImmersedConstraint;
                });
            if (realization.conservativeTransportedMass
                && !realization.phaseTransportedState
                && supports(conservativeCapabilities,*found)
                && (!needsPressureSchedule || pressureSchedule)) {
                binding.provider = "flow.conservative";
            } else if (realization.phaseTransportedState
                       && supports(eulerianCapabilities,*found)) {
                binding.provider = "flow.eulerian-pressure";
            } else if (supports(immersedCapabilities,*found)) {
                binding.provider = "ibm.constraint";
            } else {
                binding.reason = "no provider implements '"+id+"' for "
                    +realizationName(realization);
            }
        }
        if (!binding.provider.empty()) binding.status = BindingStatus::Resolved;
        bindings.push_back(std::move(binding));
    }
    return bindings;
}

RuntimeReport reportOperationBindings(
        const CompiledSolvePlan& plan,
        const std::vector<ResolvedOperationBinding>& bindings) {
    RuntimeReport report;
    report.requiredOperations = SolvePlanner::requiredOperations(plan);
    for (const OpId& id : report.requiredOperations) {
        const auto found = std::find_if(
            bindings.begin(),bindings.end(),
            [&](const ResolvedOperationBinding& binding) {
                return binding.operation == id;
            });
        if (found == bindings.end()) {
            throw std::runtime_error(
                "Provider resolver omitted plan operation '"+id+"'.");
        }
        if (found->status == BindingStatus::Unsupported) {
            report.missingOperations.push_back(id);
            if (report.reason.empty()) report.reason = found->reason;
        }
    }
    report.status = report.missingOperations.empty()
        ? RuntimeStatus::Runnable : RuntimeStatus::Unsupported;
    if (!report.missingOperations.empty()) {
        report.reason += "; unresolved operations [";
        for (std::size_t i=0;i<report.missingOperations.size();++i) {
            if (i) report.reason += ", ";
            report.reason += report.missingOperations[i];
        }
        report.reason += "]";
    }
    return report;
}

} // namespace SF::System
