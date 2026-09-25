/// @file SF_execution.cpp
/// @brief 根据 compiled execution requirements 组装 providers 并分发执行。
///
/// Data flow:
///   ResolvedSimulationSystem + CompiledSolvePlan + field resources
///       -> conservative/Eulerian provider composition
///       -> existing execution assembler
///
/// 本文件不创建第二份 solve plan，也不实现 timestep 或数值公式。

#include "app/application/execution/SF_execution.h"
#include "app/application/SF_application.h"
#include "SF_environment.h"

#include <stdexcept>

namespace SF::Application::Execution {

int executeSingle(
        Field& field, ResultWriter& writer,
        Parallel::ParallelContext& parallel, IBM::IB& ibm,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& system,
        const System::CompiledSolvePlan& plan,
        const CaseConfig& caseConfig, bool ibmEnabled,
        bool initialOutputOnly, int localBlockId) {
    if (System::requiresProvider(system,"flow.eulerian-pressure")) {
        return Detail::executeEulerianEquations(
            field, writer, parallel, solverConfig, system, plan,
            caseConfig,
            initialOutputOnly, localBlockId);
    }
    if (System::requiresProvider(system,"flow.conservative")) {
        return Detail::executeConservativeEquations(
            field, writer, parallel, ibm, solverConfig, system, plan,
            caseConfig,
            ibmEnabled, initialOutputOnly, localBlockId);
    }
    throw std::runtime_error(
        "Compiled system has no registered flow execution provider requirement.");
}

int execute(const System::ResolvedSimulationSystem& system,
            const System::CompiledSolvePlan& plan,
            const CaseConfig& caseConfig,
            ExecutionEnvironment& env,
            const RunRequest& request) {
    if (env.domain == DomainKind::DistributedMultiPatch) {
        return executeMulti(
            env.parallelMesh, *env.writer, *env.parallel, env.compositeIBM,
            caseConfig.solver, system, plan, caseConfig, env.ibmEnabled,
            request.initialOutputOnly);
    }
    return executeSingle(
        *env.activeField, *env.writer, *env.parallel, env.ibm,
        caseConfig.solver, system, plan, caseConfig, env.ibmEnabled,
        request.initialOutputOnly, env.localBlockId);
}

} // namespace SF::Application::Execution
