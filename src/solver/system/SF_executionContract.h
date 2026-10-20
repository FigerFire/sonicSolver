#pragma once
#include "core/system/SF_solveProgram.h"
namespace SF::System {
void composeDefaultPlacement(const ExecutionProgram& program,ExecutionScope& root);
void validatePlacement(const ExecutionProgram& program,const ExecutionScope& root);
/// @brief 将已选择 provider 的版本能力与唯一发布者冻结为读取契约。
void bindCoefficientContracts(const StateRegistry&,CompiledExecutionProgram&);
/// @brief Runtime 只核对冻结契约，不重选发布者或版本。
void validateCoefficientContracts(const CompiledExecutionProgram&);
/// @brief Compiler 同时核对物理符号身份与发布/消费完整性。
void validateCoefficientContracts(const StateRegistry&,const CompiledExecutionProgram&);
void validateDataFlow(const StateRegistry& state,const CompiledExecutionProgram& program);
std::string missingCapability(const CompiledExecutionProgram& program,
    const std::vector<std::string>& context,const CompiledEquationCall& consumer,
    const std::vector<CapabilityBinding>& boundContext={});
}
