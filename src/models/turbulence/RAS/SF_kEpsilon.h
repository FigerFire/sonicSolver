/// @file SF_kEpsilon.h
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#pragma once

#include "SF_RAS.h"

namespace SF {
namespace Turbulence {
namespace RAS {

/// @brief 标准 k-epsilon RAS 模型（含湍流扩散项）。
class KEpsilonModel final : public RASModelBase {
public:
    /// @brief 返回所属大类名。
    /// @return `"RAS"`。
    const char* familyName() const override { return "RAS"; }

    /// @brief 返回模型名。
    /// @return `"kEpsilon"`。
    const char* modelName() const override { return "kEpsilon"; }

    /// @brief 初始化 k/epsilon 标量状态。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void initialize(const Field& flow,
                    ScalarFields& state,
                    const FDM::TurbulenceConfig& config) override;

    /// @brief 对 k/epsilon 应用模型 floors；不执行通用 BC。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void applyPhysicalBoundary(const Field& flow,
                       ScalarFields& state,
                       const FDM::TurbulenceConfig& config) override;

    /// @brief 返回只读cell物理闭合；solver保留原位数值sweep。
    FDM::TransportCellClosure onceStepClosure(const Field&,const ScalarFields&,
        const FDM::TurbulenceConfig&,int,int,int) const override;

    /// @brief 计算湍流附加动力粘度。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return `mu_t = Cmu * rho * k^2 / epsilon`。
    double eddyDynamicViscosity(const Field& flow,
                                const ScalarFields& state,
                                const FDM::TurbulenceConfig& config,
                                int i, int j, int k) const override;
};

} // namespace RAS
} // namespace Turbulence
} // namespace SF
