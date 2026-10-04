#include "SF_legacyNumerics.h"
#include "core/system/SF_operationIds.h"
#include "SF_stateRealization.h"
#include <algorithm>
namespace SF::System::Legacy {
std::string selectOperationProvider(const ExecutableOperation& operation,
        const ExecutableEquationSystem& equations,const CompiledNumericalSystem& numerics,
        const std::vector<LegacyExecutionPolicy>& policies) {
    const auto requires=[&](OperationCapability capability) {
        return std::find(operation.requirements.begin(),operation.requirements.end(),capability)!=operation.requirements.end();
    };
    const auto realization=compileStateRealization(equations.state,equations.constraints);
    const bool pressureOperation=std::any_of(operation.requirements.begin(),operation.requirements.end(),
        [](OperationCapability capability) {
            return capability!=OperationCapability::ConservativeExplicit
                && capability!=OperationCapability::EulerianPhaseExecution
                && capability!=OperationCapability::ImmersedConstraint;
        });
    if (realization.conservativeTransportedMass && !realization.phaseTransportedState
        && !pressureOperation)
        return "flow.conservative";
    return {};
}

bool supportedLegacyFluid(const EquationDescriptor& equation) {
    const auto solves=[&](const char* symbol) {
        return std::find(equation.solvedUnknowns.begin(),
                         equation.solvedUnknowns.end(),symbol)
            !=equation.solvedUnknowns.end();
    };
    return (equation.role==LegacyEquationRole::Mass && solves("rho"))
        || (equation.role==LegacyEquationRole::Momentum
            && (solves("rhoU") || solves("U")))
        || (equation.role==LegacyEquationRole::Energy && solves("rhoE"));
}


std::vector<CompiledMathRef> spatialInputs(const ExecutableEquationSystem& equations) {
    std::vector<CompiledMathRef> result;
    for (const auto& equation:equations.legacyEquations)
        if (supportedLegacyFluid(equation) && !equation.solvedUnknowns.empty())
            result.push_back({equation.id,equation.solvedUnknowns.front()});
    return result;
}

} // namespace SF::System::Legacy
