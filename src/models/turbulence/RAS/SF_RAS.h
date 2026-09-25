/// @file SF_RAS.h
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#pragma once

#include "SF_turbulence.h"

namespace SF {
namespace Turbulence {
namespace RAS {

/// @brief RAS 两方程模型的公共基类。
class RASModelBase : public IModel {
protected:
    /// @brief 初始化单个湍流标量。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param settings 对应标量初值配置。
    /// @param slot 目标槽位。
    void initializeScalar(const Field& flow,
                          ScalarFields& state,
                          const std::vector<BCSetting<double>>& settings,
                          ScalarSlot slot) const;

    /// @brief 对单个湍流标量应用边界条件。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param settings 对应边界条件配置。
    /// @param slot 目标槽位。
    void applyScalarBC(const Field& flow,
                       ScalarFields& state,
                       const std::vector<BCSetting<double>>& settings,
                       ScalarSlot slot) const;

    /// @brief 将所有内部流体单元的附加动力粘度刷新到缓存。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void refreshEddyMu(const Field& flow,
                       ScalarFields& state,
                       const FDM::TurbulenceConfig& config) const;
};

} // namespace RAS
} // namespace Turbulence
} // namespace SF
