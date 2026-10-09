#pragma once
#include "core/system/SF_solveProgram.h"
namespace SF::System {
void composeDefaultPlacement(const ExecutionProgram& program,ExecutionScope& root);
void validatePlacement(const ExecutionProgram& program,const ExecutionScope& root);
void validateDataFlow(const StateRegistry& state,const CompiledExecutionProgram& program);
std::string missingCapability(const CompiledExecutionProgram& program,
    const std::vector<std::string>& context,const CompiledEquationCall& consumer,
    const std::vector<CapabilityBinding>& boundContext={});
}
