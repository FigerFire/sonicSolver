#pragma once

/// @file SF_couplingStatus.h
/// @brief REGISTERED MODELS — coupling/model 注册状态（不含匹配逻辑）。
///
/// 该状态必须显式出现在 explain 与 validation 中，禁止"注册了但被静默忽略"。

#include "core/system/SF_equationIR.h"
#include "core/system/SF_scheduleIds.h"
#include "core/config/types/SF_pressureConfigTypes.h"

#include <string>
#include <vector>

namespace SF::System {

enum class CouplingStatus { Active, Inactive, Invalid, Unsupported };

const char* toString(CouplingStatus status);

/// @brief 稳定 schedule id：单压力/分离压力修正策略。
///
/// planner、model contribution、printer 必须引用同一常量，避免字符串 typo
/// 变成 runtime bug。
/// @brief 稳定 schedule id：Eulerian 相共享压力耦合策略。

/// @brief 一次 coupling preset 注册（来自 case 输入，不是运行时选择）。
struct CouplingPresetRequest {
    /// @brief 输入里的 preset 名字（SIMPLE/PISO/PIMPLE），只用于报告。
    std::string preset;
    FDM::PressureCouplingPreset presetKind = FDM::PressureCouplingPreset::PISO;
    int outerCorrectors = 1;
    int pressureCorrectors = 1;
    int nonOrthogonalCorrectors = 0;
    /// @brief 用户是否显式注册（case 声明），而不是框架默认。
    bool explicitlyRegistered = false;
    Provenance origin;
};

/// @brief coupling preset 的匹配结果与贡献清单。
struct CouplingReport {
    std::string id;
    std::string preset;
    CouplingStatus status = CouplingStatus::Inactive;
    std::string reason;
    std::vector<std::string> requirements;
    std::vector<std::string> derivedOperations;
    std::vector<std::string> requiredProviders;
    std::vector<std::string> derivedEquations;
};

inline bool isCouplingActive(const CouplingReport& report) {
    return report.status == CouplingStatus::Active;
}

} // namespace SF::System
