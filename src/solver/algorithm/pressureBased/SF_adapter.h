#pragma once

/// @file SF_adapter.h
/// @brief 共享压力校正流动算法适配器。

#include "SF_flowAlgorithm.h"
#include "solver/algorithm/pressureBased/SF_corrector.h"
#include "solver/algorithm/pressureBased/SF_kkt.h"

namespace SF::PressureBased {

/// @brief 把 HYPRE 压力校正作为可注入的流动算法阶段。
class Algorithm final : public FDM::IFlowAlgorithm {
public:
    explicit Algorithm(const FDM::SolverConfig& config);

    const char* name() const override { return "pressureBased"; }
    FDM::CapabilitySet capabilities() const override;
    FDM::FlowAlgorithmResult correct(
        FDM::FlowAlgorithmContext& context,
        double dt) override;

private:
    Corrector corrector_;
    MonolithicKKT kkt_;
};

} // namespace SF::PressureBased
