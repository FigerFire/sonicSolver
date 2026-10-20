#pragma once
#include "solver/algorithm/SF_singleFluidBinding.h"
#include "solver/algorithm/time/SF_stageGroup.h"
#include "solver/algorithm/SF_stageFlux.h"
namespace SF::SolverAlgorithm {
std::vector<Time::TemporalCallbacks> scalarParticipants(const System::CompiledNumericalSystem&,
    const System::CompiledSolvePlan&,State::StateBundle&,const FDM::SolverServices&,
    System::StateRealization&,const double& maximumTimeStep,
    std::function<StageMassFluxView(StageFluxIdentity)> massFlux={});
}
