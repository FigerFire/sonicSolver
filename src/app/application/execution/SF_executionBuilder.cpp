/// @file SF_executionBuilder.cpp
/// @brief 根据 resolved mathematical strategy 组装 single-field execution provider。
///
/// Data flow:
///   ResolvedSimulationSystem + CompiledSolvePlan + field resources
///       -> conservative/Eulerian provider composition
///       -> existing execution assembler
///
/// 本文件不创建第二份 solve plan，也不实现 timestep 或数值公式。

#include "app/application/execution/SF_runners.h"
#include "app/application/execution/SF_executionAssemblers.h"

namespace SF::Application::Runners {

int runSingleField(
        Field& field, ResultWriter& writer,
        Parallel::ParallelContext& parallel, IBM::IB& ibm,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& system,
        const System::CompiledSolvePlan& plan,
        const CaseConfig& caseConfig, bool ibmEnabled,
        bool initialOutputOnly, int localBlockId) {
    if (System::hasSolveStrategy(
            system,FDM::SolveStrategyKind::PressureVelocityCoupling)) {
        return Detail::executeEulerianEquations(
            field, writer, parallel, solverConfig, system, plan,
            caseConfig,
            initialOutputOnly, localBlockId);
    }
    return Detail::executeConservativeEquations(
        field, writer, parallel, ibm, solverConfig, system, plan,
        caseConfig,
        ibmEnabled, initialOutputOnly, localBlockId);
}

} // namespace SF::Application::Runners
