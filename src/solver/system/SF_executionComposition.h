#pragma once
#include "SF_resolvedSimulationSystem.h"
#include "SF_buildRequest.h"
namespace SF::System {
/// Compose pressure contributions before freezing source inputs.
void composeContributions(ResolvedSimulationSystem&,const BuildRequest&,
    ExecutionProgram&,std::vector<NumericalBinding>&);
}
