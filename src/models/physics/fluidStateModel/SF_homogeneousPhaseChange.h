#pragma once

/// @file SF_homogeneousPhaseChange.h
/// @brief 无冗余 partialDensity FluidStateModel 的 Lee/RPI 相变源装配。

#include "SF_homogeneousMultiphaseStateModel.h"
#include "SF_phaseChange.h"
#include "SF_scalarField.h"
#include "SF_transferLedger.h"
#include "core/residual/SF_residual.h"

#include <memory>

namespace SF::Physics::FluidStateModel {

class HomogeneousPhaseChange {
public:
    HomogeneousPhaseChange(
        Multiphase::MultiPhaseConfig config,
        std::shared_ptr<const HomogeneousMultiphaseStateModel> equations);

    /// @brief 分配派生 alpha/T 缓存并检查相与传质组分映射。
    void initialize(const Field& field);

    /// @brief 从当前状态派生 alpha/T，经 TransferLedger 映射守恒源。
    /// 总能量内部源严格为零；潜热只由两相 EOS referenceEnergy 差表达。
    void assemble(Field& field, Residual& residual, double dt);

    const PhaseChange::PhaseChangeDiagnostics& diagnostics() const {
        return diagnostics_;
    }

    /// @brief 最近一个 RK stage 的相间交换账本。
    const InterphaseTransfer::TransferLedger& ledger() const {
        return ledger_;
    }

private:
    Multiphase::MultiPhaseConfig config_;
    std::shared_ptr<const HomogeneousMultiphaseStateModel> equations_;
    ScalarField alpha_;
    ScalarField temperature_;
    int liquidComponent_ = -1;
    int vaporComponent_ = -1;
    int liquidPhase_ = -1;
    int vaporPhase_ = -1;
    int trackedPhase_ = -1;
    PhaseChange::PhaseChangeDiagnostics diagnostics_;
    InterphaseTransfer::TransferLedger ledger_;

    void refreshDerived(const Field& field);
};

} // namespace SF::Physics::FluidStateModel
