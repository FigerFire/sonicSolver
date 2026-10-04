#include "SF_executionComposition.h"
#include "SF_pressureCoupling.h"
#include "SF_eulerianCoupling.h"
#include "SF_legacyNumerics.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
namespace SF::System {
std::vector<LegacyPlanFragment> composeContributions(ResolvedSimulationSystem& result,const BuildRequest& request,
        ExecutionProgram& program,std::vector<NumericalBinding>& numerics) {
    // Module applicability follows contributed mathematics/default occurrences.
    // StateRealization only binds storage and never chooses this source topology.
    const auto predictor=std::find_if(program.root.children.begin(),program.root.children.end(),[](const auto& node) {
        return node.kind==ExecutionKind::EquationCall && node.step.equation=="momentum"
            && (node.step.target.symbol=="U" || node.step.target.symbol=="rhoU");
    });
    const bool native=request.coupling && isCouplingActive(result.coupling)
        && hasConstraint(result.rawSystem,"C_INCOMPRESSIBILITY")
        && predictor!=program.root.children.end();
    if (native) {
        const auto target=predictor->step.target.symbol;
        applyPressureExecution(program,*request.coupling,target);
        if (target=="rhoU") {
            numerics.erase(std::remove_if(numerics.begin(),numerics.end(),[](const auto& binding) {
                return binding.method=="ConservativeResidual" && binding.occurrence.empty()
                    && (binding.equation=="continuity" || binding.equation=="momentum" || binding.equation=="energy");
            }),numerics.end());
        }
        const auto pressure=pressureNumerics(*request.coupling,target);
        numerics.insert(numerics.end(),pressure.begin(),pressure.end());
    }
    if (hasConstraint(result.rawSystem,"C_SHARED_PRESSURE")) {
        if (!request.coupling || !isCouplingActive(result.coupling))
            throw std::runtime_error("Eulerian shared pressure requires an explicit supported native HOW coupling preset.");
        applyEulerianExecution(program,numerics,*request.coupling);
    } else if (request.coupling && isCouplingActive(result.coupling) && !native) {
        throw std::runtime_error("Single-fluid pressure requires a native default predictor occurrence; no legacy fallback is available.");
    }
    return {};
}
} // namespace SF::System
