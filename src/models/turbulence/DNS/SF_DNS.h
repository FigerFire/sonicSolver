/// @file SF_DNS.h
/// @brief DNS 无湍流闭式分支及能力声明。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#pragma once

#include "SF_turbulence.h"

namespace SF {
namespace Turbulence {
namespace DNS {

/// @brief DNS 模式，不添加任何附加湍流粘度。
class DirectNumericalSimulationModel final : public IModel {
public:
    /// @brief 返回所属大类名。
    /// @return `"DNS"`。
    const char* familyName() const override { return "DNS"; }

    /// @brief 返回模型名。
    /// @return `"DNS"`。
    const char* modelName() const override { return "DNS"; }

    /// @brief DNS 无额外标量状态，此处只清零缓存。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void initialize(const Field& flow,
                    ScalarFields& state,
                    const FDM::TurbulenceConfig& config) override;

    /// @brief DNS 无额外边界处理。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void applyBoundary(const Field& flow,
                       ScalarFields& state,
                       const FDM::TurbulenceConfig& config) override;

    /// @brief DNS 无附加湍流演化，仅保持附加粘度为0。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    /// @param dt 当前时间步长。
    void correct(const Field& flow,
                 ScalarFields& state,
                 const FDM::TurbulenceConfig& config,
                 double dt) override;

    /// @brief 返回0，表示不添加湍流附加粘度。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 0.0。
    double eddyDynamicViscosity(const Field& flow,
                                const ScalarFields& state,
                                const FDM::TurbulenceConfig& config,
                                int i, int j, int k) const override;
};

} // namespace DNS
} // namespace Turbulence
} // namespace SF
