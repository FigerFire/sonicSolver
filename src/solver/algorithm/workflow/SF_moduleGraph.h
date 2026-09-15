#pragma once

/// @file SF_moduleGraph.h
/// @brief 运行装配阶段的模块能力图和缺失能力诊断。

#include "SF_capability.h"
#include "SF_config.h"
#include "SF_resolvedSimulationSystem.h"

#include <string>
#include <vector>

namespace SF::Workflow {

/// @brief 一个已选择模块提供和要求的能力声明。
struct ModuleDescriptor {
    std::string name;
    FDM::CapabilitySet provides;
    FDM::CapabilitySet requiresAll;
    std::vector<FDM::CapabilitySet> requiresOneOf;
};

/// @brief 经配置选择后的模块依赖图。
class ModuleGraph {
public:
    void add(ModuleDescriptor module);
    FDM::CapabilitySet providedCapabilities() const;
    void validate() const;
    std::string describe() const;

private:
    std::vector<ModuleDescriptor> modules_;
};

/// @brief 根据 workflow 请求和求解配置建立可组合模块图。
ModuleGraph makeModuleGraph(
    const System::ResolvedSimulationSystem& system,
    const FDM::SolverConfig& config);

} // namespace SF::Workflow
