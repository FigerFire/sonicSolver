#pragma once

/// @file SF_correction.h
/// @brief 密度基 predictor 后的流动算法阶段。

#include "SF_flowAlgorithm.h"

namespace SF::DensityBased {

/// @brief 守恒密度基算法；主更新由通用显式 predictor 完成。
class Algorithm final : public FDM::IFlowAlgorithm {
public:
    const char* name() const override { return "densityBased"; }
    FDM::CapabilitySet capabilities() const override;
    FDM::FlowAlgorithmResult correct(
        FDM::FlowAlgorithmContext& context,
        double dt) override;
};

} // namespace SF::DensityBased
