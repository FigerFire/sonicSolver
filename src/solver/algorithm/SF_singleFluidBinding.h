#pragma once
#include "core/interfaces/SF_solverStepper.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_stateRealizer.h"
namespace SF::SolverAlgorithm {
void bindSingleFluidOperations(Run::OpRegistry&,const FDM::SolverConfig&,
    const System::ExecutableEquationSystem&,const System::CompiledNumericalSystem&,
    const System::CompiledSolvePlan&,const System::RuntimeRequirements&,
    State::StateBundle&,const FDM::SolverServices&,System::StateRealization&,const double& maximumTimeStep);
}
