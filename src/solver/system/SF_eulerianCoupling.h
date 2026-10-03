#pragma once
#include "SF_couplingStatus.h"
#include "SF_methodObjects.h"
namespace SF::System {
/// @brief 显式 coupling preset -> 通用 EquationCall/Loop/Commit；无 runtime OpIds。
void applyEulerianExecution(ExecutionProgram& program,std::vector<NumericalBinding>& bindings,
    const CouplingPresetRequest& request);
/// @brief 注册逐相数学方法；运行实现共用现有专用 backend。
void addEulerianMethods(ProviderRegistry& providers);
}
