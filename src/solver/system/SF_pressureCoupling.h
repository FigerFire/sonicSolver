#pragma once

/// @file SF_pressureCoupling.h
/// Pressure coupling contributes mathematics, generic HOW and independent WHICH.
/// Single-fluid counts/order live only in HOW. Eulerian shared pressure retains
/// an explicit Legacy compatibility fragment until its own migration.

#include "SF_couplingStatus.h"
#include "SF_equationContribution.h"
#include "core/system/SF_planFragment.h"
#include "core/config/types/SF_pressureConfigTypes.h"

#include <string>
#include <vector>

namespace SF::System {

/// @brief 既有 conservative Corrector 的数学关系；同一描述用于注册与 provider 验证。
std::vector<Equation> conservativePressureRelations();

/// @brief 匹配 requirements，不做任何 contribution。
CouplingReport matchPressureCoupling(
    const CouplingPresetRequest& request,
    const RawEquationSystem& raw);

/// @brief 把 active preset 的 formulation/plan 贡献写入 composition。
///
/// 只有 status==Active 时才写入 transformation 请求；shared Eulerian 另有 Legacy policy；
/// inactive/invalid/unsupported 一律不写，并由调用方按状态报告或 fail。
/// @return 同一个 report（带上 derived operations/equations/providers）。
CouplingReport contributePressureCoupling(
    SystemCompositionBuilder& system,
    const CouplingPresetRequest& request,
    const RawEquationSystem& raw);

/// @brief 从 case 输入解析 coupling preset 请求。
CouplingPresetRequest couplingRequestFrom(
    const FDM::PressureCouplingConfig& coupling,
    bool explicitlyRegistered);

/// Single-fluid coupling contributes HOW directly. The
/// EquationMethods, not this tree, supply equation-level numerical leaves.
void applyPressureExecution(ExecutionProgram& program,
                            const CouplingPresetRequest& request,
                            std::string_view momentumTarget = "U");
/// Independent WHICH contribution; never mutates HOW.
std::vector<NumericalBinding> pressureNumerics(const CouplingPresetRequest& request,
    std::string_view momentumTarget = "U");

namespace Legacy {
/// Only unmigrated Eulerian shared pressure uses this compatibility fragment.

} // namespace Legacy

/// @brief 片段引用的 stage 在 executable system 中是否都有实现。
///
/// @param missing 输出缺失（未解析）的 operation id 列表。
bool couplingPlanResolved(const LegacyPlanFragment& fragment,
                          const ExecutableEquationSystem& executable,
                          std::vector<std::string>* missing);

} // namespace SF::System
