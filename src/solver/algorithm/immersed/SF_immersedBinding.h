#pragma once
#include "core/interfaces/SF_solverStepper.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_runtimeRequirements.h"
namespace SF::SolverAlgorithm {
void bindImmersedOperations(Run::OpRegistry&, const System::RuntimeRequirements&,
    const System::CompiledSolvePlan&, const FDM::SolverConfig&, State::StateBundle&,
    const FDM::SolverServices&, FDM::ImmersedConstraintResult&,
    std::function<void(const std::vector<Field*>&,double,double)> prepareBoundary);
}
