#pragma once
#include "solver/algorithm/SF_singleFluidBinding.h"
#include "solver/algorithm/time/SF_stageGroup.h"
#include "solver/algorithm/SF_stageFlux.h"
namespace SF::SolverAlgorithm {
/// @brief 仅绑定所选代数实现及其输入；不创建时间循环。
std::vector<std::function<void()>> bindAlgebraicOperations(Run::OpRegistry&,const System::CompiledNumericalSystem&,
    const System::CompiledSolvePlan&,State::StateBundle&,const FDM::SolverServices&,System::StateRealization&);
std::vector<Time::TemporalCallbacks> scalarParticipants(const System::CompiledNumericalSystem&,
    const System::CompiledSolvePlan&,State::StateBundle&,const FDM::SolverServices&,
    System::StateRealization&,const double& maximumTimeStep,
    std::function<StageMassFluxView(StageFluxIdentity)> massFlux={});
}
