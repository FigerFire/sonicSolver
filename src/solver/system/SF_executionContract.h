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
/// @brief 从冻结的调用与约束数学引用追溯消费者，不复制数学、不推断调度。
std::vector<const CompiledEquationCall*> equationConsumers(
    const ExecutableEquationSystem&,const CompiledExecutionProgram&,std::string_view equation);
/// @brief 只读查询实际 Plan leaves；保留融合/约束实现关系，不创建第二份计划。
std::vector<const SolvePlanNode*> equationOperations(
    const ExecutableEquationSystem&,const CompiledSolvePlan&,std::string_view equation);
/// @brief 拒绝活动数学无编译去向、悬空约束/融合引用；Unsupported 保留阻断说明。
void validateMathematicalCoverage(const ExecutableEquationSystem&,const CompiledSolvePlan&,bool runnable);
void validateDataFlow(const StateRegistry& state,const CompiledExecutionProgram& program);
std::string missingCapability(const CompiledExecutionProgram& program,
    const std::vector<std::string>& context,const CompiledEquationCall& consumer,
    const std::vector<CapabilityBinding>& boundContext={});
}
