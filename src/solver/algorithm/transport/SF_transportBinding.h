#pragma once
#include "core/interfaces/SF_solverStepper.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_runtimeRequirements.h"
namespace SF::SolverAlgorithm {
void bindTransportOperations(Run::OpRegistry&, const System::RuntimeRequirements&,
                             State::StateBundle&, const FDM::SolverServices&);
}
