#pragma once

/// @file SF_phaseState.h
/// @brief Eulerian 多流体每相独立状态；不写入主 `Field` 守恒量槽位。

#include "SF_scalarField.h"
#include "SF_phaseProperties.h"

#include <array>
#include <string>

namespace SF::Physics::PhaseSystems {

/// @brief Eulerian 相三分量场；分量保持全局 Cartesian 方向。
using PhaseVectorField = std::array<ScalarField, 3>;

/// @brief 每相唯一可推进主状态。
struct PhasePrimaryState {
    ScalarField phaseMass;       ///< alpha_k * rho_k。
    PhaseVectorField momentum;   ///< alpha_k * rho_k * U_k。
    ScalarField phaseEnthalpy;   ///< alpha_k * rho_k * h_k。

    /// @brief 按参考网格分配主状态。
    void setupLike(const Field& field, const std::string& phaseName);
};

/// @brief 由主状态和共享压力恢复的原始量缓存。
struct PhasePrimitiveCache {
    ScalarField alpha;
    ScalarField density;
    PhaseVectorField velocity;
    ScalarField temperature;
    ScalarField enthalpy;

    /// @brief 按参考网格分配派生缓存。
    void setupLike(const Field& field, const std::string& phaseName);
};

/// @brief 一个 Eulerian 相的主状态和非拥有语义的派生缓存。
struct PhaseState {
    std::string name;
    PhasePrimaryState primary;
    PhasePrimitiveCache primitive;

    /// @brief 按参考网格分配相状态。
    void setupLike(const Field& field, const std::string& phaseName);

    /// @brief 仅用于初值和边界提交，从原始量构造唯一主状态。
    void commitPrimitiveToPrimary(
        const Multiphase::PhaseProperties& properties);
};

} // namespace SF::Physics::PhaseSystems
