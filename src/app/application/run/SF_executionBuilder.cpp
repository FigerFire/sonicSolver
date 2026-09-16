/// @file SF_executionBuilder.cpp
/// @brief 从 resolved solve blocks 选择 single-field executor composition。

#include "SF_runners.h"
#include "SF_executionAssemblers.h"

namespace SF::Application::Runners {

int runSingleField(
        Field& field, ResultWriter& writer,
        Parallel::ParallelContext& parallel, IBM::IB& ibm,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& system,
        const Workflow::Plan& workflow,
        const CaseConfig& caseConfig, bool ibmEnabled,
        bool initialOutputOnly, int localBlockId) {
    if (System::hasSolveStrategy(
            system,FDM::SolveStrategyKind::PressureVelocityCoupling)) {
        return Detail::executeEulerianEquations(
            field, writer, parallel, solverConfig, system, workflow,
            caseConfig,
            initialOutputOnly, localBlockId);
    }
    return Detail::executeConservativeEquations(
        field, writer, parallel, ibm, solverConfig, system, workflow,
        caseConfig,
        ibmEnabled, initialOutputOnly, localBlockId);
}

} // namespace SF::Application::Runners
