#pragma once
/// @file SF_eulerianTurbulence.h
/// @brief Eulerian 湍流编译绑定；不持有数组，不复制数学或运行时状态。
#include "SF_methodObjects.h"
#include <memory>
#include "core/config/types/SF_turbulenceConfigTypes.h"

namespace SF::System {
/// @brief 一项显式 model phase 在 PhaseSystem 中的真实槽位。
struct EulerianTurbulencePhaseBinding {
    std::string name;
    std::size_t slot;
};
/// @brief 所选 provider 验证 WHAT/STATE 后冻结的实现及 phase/model 分组。
struct CompiledEulerianTurbulenceContract {
    EquationRef equation;
    Target target;
    std::string method;
    FDM::TurbulenceModelKind model;
    bool transport=false;
    std::vector<EulerianTurbulencePhaseBinding> phases;
};
using EulerianTurbulenceContractPtr=std::shared_ptr<const CompiledEulerianTurbulenceContract>;
inline const CompiledEulerianTurbulenceContract* eulerianTurbulenceContract(const CompiledEquationCall& call) {
    const auto* value=std::any_cast<EulerianTurbulenceContractPtr>(&call.providerContract);
    return value ? value->get() : nullptr;
}
/// @brief 验证冻结契约与实际模型存储顺序，runtime 不重新解释 equation identity。
void validateEulerianTurbulenceBindings(const CompiledExecutionProgram& program,
    const std::vector<std::string>& phases,const FDM::TurbulenceConfig& config);
void addEulerianTurbulenceMethods(ProviderRegistry& providers);
}
