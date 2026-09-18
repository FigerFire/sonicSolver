#pragma once

/// @file SF_executionAssemblers.h
/// @brief single-field execution 的内部 state/provider 装配入口。
///
/// Data flow:
///   field + compiled plan + resolved providers
///       -> conservative or Eulerian execution composition
///       -> runFlow orchestration
///
/// 本文件不定义 timestep lifecycle、physics contribution 或 numerical kernel。

#include "SF_config.h"
#include "SF_resolvedSimulationSystem.h"

namespace SF {
class Field;
class ResultWriter;
struct CaseConfig;
namespace IBM { class IB; }
namespace Parallel { class ParallelContext; }
namespace Application::Runners::Detail {

int executeConservativeEquations(
    Field&, ResultWriter&, Parallel::ParallelContext&, IBM::IB&,
    const FDM::SolverConfig&, const System::ResolvedSimulationSystem&,
    const System::CompiledSolvePlan&,
    const CaseConfig&, bool ibmEnabled, bool initialOutputOnly,
    int localBlockId);

int executeEulerianEquations(
    Field&, ResultWriter&, Parallel::ParallelContext&,
    const FDM::SolverConfig&, const System::ResolvedSimulationSystem&,
    const System::CompiledSolvePlan&,
    const CaseConfig&, bool initialOutputOnly,
    int localBlockId);

} // namespace Application::Runners::Detail
} // namespace SF
