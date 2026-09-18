#pragma once

/// @file SF_runners.h
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
namespace Application::Runners {

int runMultiPatch(
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

int runSingleField(
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

} // namespace Application::Runners
} // namespace SF
