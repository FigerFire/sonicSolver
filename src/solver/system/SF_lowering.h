#pragma once
/// @file SF_lowering.h
/// @brief 唯一四轴编译协调入口；普通函数，无场数组、clock 或运行时生命周期。
#include "SF_methodObjects.h"
#include "SF_numericalSelection.h"
#include "SF_resolvedSimulationSystem.h"
namespace SF::System::Lowering {
/// @brief 解析全部显式四轴输入，调用已选实现并冻结调用、资源需求和控制流。
CompiledExecutionProgram compileProgram(const ExecutableEquationSystem&,const StateRegistry&,
    const ExecutionProgram&,const std::vector<NumericalBinding>&,const ProviderRegistry&,
    const CompiledTimeRecipe* = nullptr,const ITemporalMethod* = nullptr);
/// @brief 生产组合的最终协调；空间 lowering 与 runtime 能力均消费同一编译结果。
void compileSystem(ResolvedSimulationSystem&,const ExecutionProgram&,NumericalSelection&,
    bool additionalContributions=false);
/// @brief 校验冻结的 temporal 拓扑，不执行 stage、不选择方法。
std::string validateTemporalPlan(const CompiledSolvePlan&,int stages);
}
