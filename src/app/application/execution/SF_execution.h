#pragma once

/// @file SF_execution.h
/// @brief 暴露 compiled program 到 concrete execution composition 的入口。
///
/// Data flow:
///   ResolvedSimulationSystem + CompiledSolvePlan + ExecutionEnvironment data
///       -> single-field or multi-patch runtime composition
///       -> existing numerical stepper
///
/// 该 contract 不解析命令行、不编译 mathematical system，也不定义数值公式。

#include "SF_config.h"
#include "SF_resolvedSimulationSystem.h"
#include "SF_flowLoop.h"

#include <string>

namespace SF {
class Field;
class ResultWriter;
class MultiBlockMesh;
struct CaseConfig;
namespace IBM {
class IB;
class CompositeIB;
}
namespace Parallel {
class ParallelContext;
}
namespace Application {

struct RunRequest;
struct ExecutionEnvironment;

namespace Execution {

/// @brief 执行 multi-patch 复合网格上的运行。
int executeMulti(
    MultiBlockMesh& mesh,
    ResultWriter& writer,
    Parallel::ParallelContext& parallel,
    IBM::CompositeIB& ibm,
    const FDM::SolverConfig& solverConfig,
    const System::ResolvedSimulationSystem& system,
    const System::CompiledSolvePlan& plan,
    const CaseConfig& caseConfig,
    bool ibmEnabled,
    bool initialOutputOnly);

/// @brief 执行 single-field 布局上的运行。
int executeSingle(
    Field& field,
    ResultWriter& writer,
    Parallel::ParallelContext& parallel,
    IBM::IB& ibm,
    const FDM::SolverConfig& solverConfig,
    const System::ResolvedSimulationSystem& system,
    const System::CompiledSolvePlan& plan,
    const CaseConfig& caseConfig,
    bool ibmEnabled,
    bool initialOutputOnly,
    int localBlockId);

/// @brief 只跑 CompiledSolvePlan；storage 布局由 env.domain 决定。
///
/// single / multi 只存在于 execution environment 的拓扑里，不再暴露成
/// application runner 的两个独立入口名。
int execute(const System::ResolvedSimulationSystem& system,
            const System::CompiledSolvePlan& plan,
            const CaseConfig& caseConfig,
            ExecutionEnvironment& env,
            const RunRequest& request);

namespace Detail {

/// @brief conservative single-field 方程装配入口。
int executeConservativeEquations(
    Field&, ResultWriter&, Parallel::ParallelContext&, IBM::IB&,
    const FDM::SolverConfig&, const System::ResolvedSimulationSystem&,
    const System::CompiledSolvePlan&,
    const CaseConfig&, bool ibmEnabled, bool initialOutputOnly,
    int localBlockId);

/// @brief Eulerian pressure 方程装配入口。
int executeEulerianEquations(
    Field&, ResultWriter&, Parallel::ParallelContext&,
    const FDM::SolverConfig&, const System::ResolvedSimulationSystem&,
    const System::CompiledSolvePlan&,
    const CaseConfig&, bool initialOutputOnly,
    int localBlockId);

} // namespace Detail
} // namespace Execution
} // namespace Application
} // namespace SF
