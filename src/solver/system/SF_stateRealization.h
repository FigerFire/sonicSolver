#pragma once

/// @file SF_stateRealization.h
/// @brief STATE compiler output: base roles/storage and lazy semantic views.
/// Role summaries are reporting/backend capability data; they never schedule HOW.

#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"

#include <string>
#include <vector>

namespace SF::System {

/// @brief 一个 role group 的编译期视图。
struct RealizedRoleGroup {
    std::string id;
    StateRole role = StateRole::Primary;
    int components = 1;
    std::string storageKey;
    int componentOffset = 0;
    std::optional<double> constantValue;
    std::string nameSpace;
};

/// @brief STATE specification 导出的角色摘要。
struct CompiledStateRealization {
    /// @brief rho/rhoU/rhoE 之类的守恒量被 transported（density formulation）。
    bool conservativeTransportedMass = false;
    /// @brief 动量守恒量被 transported。
    bool conservativeMomentum = false;
    /// @brief 压力由代数约束的 multiplier 承担（div(U)=0 类约束）。
    bool pressureMultiplier = false;
    /// @brief 压力由 EOS closure 派生（p = p(rho,T)）。
    bool thermodynamicPressure = false;
    /// @brief 每相独立 transported state（Eulerian-Eulerian）。
    bool phaseTransportedState = false;
    /// @brief 存在 multiplier/constraint 对（DLM/KKT/约束代数）。
    bool multiplierConstraints = false;
    /// @brief 全部约束 id（稳定顺序）。
    std::vector<std::string> constraints;
    /// @brief transported / derived / multiplier 三组 role group。
    std::vector<RealizedRoleGroup> transported;
    std::vector<RealizedRoleGroup> derived;
    std::vector<RealizedRoleGroup> multipliers;

    /// @brief 该 realization 是否由守恒 transported state 驱动。
    bool conservativeState() const {
        return conservativeTransportedMass && conservativeMomentum;
    }
};

/// @brief 从 STATE specification 与已声明约束导出角色摘要。
CompiledStateRealization compileStateRealization(
    const StateRegistry& state,const std::vector<ConstraintDescriptor>& constraints);

/// @brief Generic state binding; providers may demand backend-owned storage only.
CompiledStateView realizeTarget(const StateRegistry& state,CompiledTarget& target);
/// @brief Temporal methods request versions without changing source HOW ordering.
void realizeTemporalViews(const StateRegistry& state,CompiledExecutionProgram& program,
    int stages);

} // namespace SF::System
