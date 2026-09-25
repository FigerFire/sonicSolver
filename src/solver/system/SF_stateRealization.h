#pragma once

/// @file SF_stateRealization.h
/// @brief STATE REALIZATION（编译期）— 由 executable equation system 的
///        unknown/constraint role 导出的求解状态语义。
///
/// 它回答：
///   - 哪些未知量是 transported conservative state？
///   - 哪些是 derived closure（例如 rho = rho0）？
///   - 哪些是 multiplier（例如 p 施加 div(U)=0）？
///   - 每个 role group 的 storage/index binding 是什么？
///
/// 它不回答 runtime topology（rank/patch/field pointer）；那些属于
/// `System::StateRealization`（realizeState）与 ExecutionEnvironment。
/// 关键不变量：这里没有任何 solver family 开关。同一条 executable system
/// 永远导出同一份 realization。

#include "core/system/SF_equationIR.h"

#include <string>
#include <vector>

namespace SF::System {

/// @brief 一个 role group 的编译期视图。
struct RealizedRoleGroup {
    std::string id;
    UnknownRole role = UnknownRole::Primary;
    int components = 1;
    std::string storageKey;
    int componentOffset = 0;
    std::string nameSpace;
};

/// @brief ExecutableEquationSystem 导出的 state realization。
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

/// @brief 从 executable equation system 导出 state realization。
CompiledStateRealization compileStateRealization(
    const ExecutableEquationSystem& system);

} // namespace SF::System
