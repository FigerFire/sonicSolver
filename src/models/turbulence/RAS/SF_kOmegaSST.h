/// @file SF_kOmegaSST.h
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#pragma once

#include "SF_RAS.h"

namespace SF {
namespace Turbulence {
namespace RAS {

/// @brief Menter k-omega SST 模型。
///
/// 当前实现包含:
/// - Menter SST `F1/F2` 混合函数
/// - k 方程生产项 limiter 与 omega 方程交叉扩散项
/// - `mu_t = rho*a1*k/max(a1*omega, |Omega|*F2)` 湍黏度限幅
/// - 湍流标量扩散项 ∇·((μ + σ μ_t)∇φ)
class KOmegaSSTModel final : public RASModelBase {
public:
    /// @brief 壁面在 composition 时稳定注入；不改变时间循环。
    explicit KOmegaSSTModel(const FDM::IImmersedTurbulenceBoundary* wall=nullptr):wall_(wall) {}
    /// @brief 返回所属大类名。
    /// @return `"RAS"`。
    const char* familyName() const override { return "RAS"; }

    /// @brief 返回模型名。
    /// @return `"kOmegaSST"`。
    const char* modelName() const override { return "kOmegaSST"; }

    /// @brief 初始化 k/omega 标量状态。
    /// @param flow 主流场。
    /// @param state 湍流标量存储。
    /// @param config 强类型湍流配置。
    void initialize(const Field& flow,
                    ScalarFields& state,
                    const FDM::TurbulenceConfig& config) override;

    /// @brief 对 k/omega 应用 floors 和专属壁面关系；不执行通用 BC。
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
    /// @return 基于 Menter SST `F2` 限幅的 `mu_t`。
    double eddyDynamicViscosity(const Field& flow,
                                const ScalarFields& state,
                                const FDM::TurbulenceConfig& config,
                                int i, int j, int k) const override;
    FDM::TransportStageCoefficients stageCoefficients(const Field&,const FDM::TransportStageSamples&,
        const FDM::TurbulenceConfig&,int,int,int) const override;
private:
    const FDM::IImmersedTurbulenceBoundary* wall_;
};

} // namespace RAS
} // namespace Turbulence
} // namespace SF
