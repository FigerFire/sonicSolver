#pragma once
#include "core/interfaces/SF_solverStepper.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_stateRealizer.h"
namespace SF::SolverAlgorithm {
void bindLevelSetOperations(Run::OpRegistry&, const System::CompiledSolvePlan&,
                           State::StateBundle&, const FDM::SolverServices&,
                           System::StateRealization&);
}
