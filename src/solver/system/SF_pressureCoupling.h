#pragma once

/// @file SF_pressureCoupling.h
/// @brief FORMULATION — 压力–速度耦合 preset 的贡献与状态。
///
/// SIMPLE/PISO/PIMPLE 不是"求解器家族"，而是 pressure-constraint coupling
/// preset。注册一个 preset 只会贡献：
///
///   1. formulation/transformation 请求（pressure-constraint formulation）
///   2. derived algorithmic operations（predictor / assemble / solve /
///      velocity-correct / flux-correct 由 SolvePlanner 从 policy 展开）
///   3. solve-plan fragment（execution policy：外层/压力/非正交重复）
///
/// 它是否生效由 **resolved equation system 是否满足其数学要求** 决定，
/// 而不是由 density/pressure 标签决定。四种状态必须显式报告：
///
///   active       要求全部满足
///   inactive     要求不满足（例如守恒 transported rho 系统没有 div 约束）
///   invalid      注册本身与 resolved system 矛盾（必须 fail 配置校验）
///   unsupported  数学上成立，但缺少执行所需的 provider/operation

#include "SF_couplingStatus.h"
#include "SF_equationContribution.h"
#include "core/system/SF_planFragment.h"
#include "core/config/types/SF_pressureConfigTypes.h"

#include <string>
#include <vector>

namespace SF::System {

/// @brief 匹配 requirements，不做任何 contribution。
CouplingReport matchPressureCoupling(
    const CouplingPresetRequest& request,
    const RawEquationSystem& raw);

/// @brief 把 active preset 的 formulation/plan 贡献写入 composition。
///
/// 只有 status==Active 时才写入 transformation 请求与 execution policy；
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

/// @brief SolvePlanner 为某个 preset 展开的 operation id（稳定顺序）。
/// @brief 把 active preset 的执行顺序放入片段（formulation 之后调用）。
///
/// 片段只引用 `OperationStage`；具体 OpId 来自 executable operation
/// authority。formulation 未声明该 stage 时使用具名缺失标记，用于 Unsupported
/// 报告。
PlanFragment couplingPlanFragment(const CouplingPresetRequest& request);

/// @brief 相共享压力的执行片段；loop 与 leaf 顺序由 coupling contribution 持有。
PlanFragment sharedPressurePlanFragment(
    const ExecutionPolicy& policy,
    const ExecutableEquationSystem& executable);

/// @brief 片段引用的 operation 名字（顺序无关，仅用于报告与诊断）。
std::vector<std::string> couplingReportOperations(
    const CouplingPresetRequest& request);

/// @brief 片段引用的 stage 在 executable system 中是否都有实现。
///
/// @param missing 输出缺失（未解析）的 operation id 列表。
bool couplingPlanResolved(const PlanFragment& fragment,
                          const ExecutableEquationSystem& executable,
                          std::vector<std::string>* missing);

} // namespace SF::System
