#include "SF_executionComposition.h"
#include "SF_pressureCoupling.h"
#include "SF_eulerianCoupling.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
namespace SF::System {
void composeContributions(ResolvedSimulationSystem& result,const BuildRequest& request,
        ExecutionProgram& program,std::vector<NumericalBinding>& numerics,const PressureCouplingInputs* inputs) {
    const bool native=request.coupling && isCouplingActive(result.coupling) && inputs;
    if (native) {
        applyPressureExecution(program,*request.coupling,*inputs);
        numerics.erase(std::remove_if(numerics.begin(),numerics.end(),[&](const auto& binding) {
            if (binding.method!="ConservativeResidual") return false;
            const auto refers=[&](const EquationCall& call) {
                return binding.equation==call.equation && binding.occurrence==call.occurrence;
            };
            return refers(inputs->predictor) || std::any_of(inputs->fusedBalances.begin(),inputs->fusedBalances.end(),refers);
        }),numerics.end());
        const auto pressure=pressureNumerics(*request.coupling,*inputs);
        numerics.insert(numerics.end(),pressure.begin(),pressure.end());
    }
    if (hasConstraint(result.rawSystem,"C_SHARED_PRESSURE")) {
        if (!request.coupling || !isCouplingActive(result.coupling))
            throw std::runtime_error("Eulerian shared pressure requires an explicit supported native HOW coupling preset.");
        applyEulerianExecution(program,numerics,*request.coupling);
    } else if (request.coupling && isCouplingActive(result.coupling) && !native) {
        throw std::runtime_error("Single-fluid pressure requires explicit pressure method inputs; no legacy fallback is available.");
    }
}
} // namespace SF::System
