/// @file SF_LES.h
/// @brief LES 湍流模型选择与推进实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#pragma once

#include "SF_turbulence.h"

namespace SF {
namespace Turbulence {
namespace LES {

/// @brief Smagorinsky LES 模型。
class SmagorinskyModel final : public IModel {
public:
    /// @brief 返回所属大类名。
    /// @return `"LES"`。
    const char* familyName() const override { return "LES"; }

    /// @brief 返回模型名。
    /// @return `"Smagorinsky"`。
    const char* modelName() const override { return "Smagorinsky"; }

    /// @brief 初始化 LES 状态。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void initialize(const Field& flow,
                    ScalarFields& state,
                    const FDM::TurbulenceConfig& config) override;

    /// @brief LES 无独立标量边界，此处仅保留接口。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void applyBoundary(const Field& flow,
                       ScalarFields& state,
                       const FDM::TurbulenceConfig& config) override;

    /// @brief 基于局部应变率刷新 Smagorinsky 黏度。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    /// @param dt 当前时间步长（此模型不显式使用）。
    void correct(const Field& flow,
                 ScalarFields& state,
                 const FDM::TurbulenceConfig& config,
                 double dt) override;

    /// @brief 返回当前单元的 Smagorinsky 附加动力粘度。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return `mu_t = rho (Cs Delta)^2 |S|`。
    double eddyDynamicViscosity(const Field& flow,
                                const ScalarFields& state,
                                const FDM::TurbulenceConfig& config,
                                int i, int j, int k) const override;
};

} // namespace LES
} // namespace Turbulence
} // namespace SF
