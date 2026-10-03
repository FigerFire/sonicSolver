#pragma once
#include "SF_transformation.h"
#include "core/system/SF_solveProgram.h"
namespace SF::System { struct CompiledNumericalSystem; }
namespace SF::System::Legacy {
std::vector<CompiledMathRef> spatialInputs(const ExecutableEquationSystem&);
/// Only explicit legacy plan leaves enter this compile-time compatibility adapter.
std::string selectOperationProvider(const ExecutableOperation&,const ExecutableEquationSystem&,
    const CompiledNumericalSystem&,const std::vector<LegacyExecutionPolicy>&);
}
