#pragma once
#include "SF_resolvedSimulationSystem.h"
#include "SF_buildRequest.h"
#include "core/system/SF_planFragment.h"
namespace SF::System {
/// Compose pressure and explicit legacy contributions before freezing source inputs.
std::vector<LegacyPlanFragment> composeContributions(ResolvedSimulationSystem&,const BuildRequest&,
    ExecutionProgram&,std::vector<NumericalBinding>&);
}
