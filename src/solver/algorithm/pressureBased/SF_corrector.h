#pragma once

/// @file SF_corrector.h
/// @brief 单流体压力基 HYPRE 压力—速度修正。

#include "SF_config.h"
#include "SF_field.h"
#include "SF_interfaces.h"
#include "solver/linearAlgebra/SF_linearAlgebra.h"

#include <functional>

namespace SF::PressureBased {

/// @brief 一次压力修正的统一 HYPRE 收敛信息。
struct CorrectionSummary {
    int iterations = 0;
    int correctedCells = 0;
    double initialResidual = 0.0;
    double finalResidual = 0.0;
    double maxDivergenceBefore = 0.0;
    double maxDivergenceAfter = 0.0;
    int interfaceFaces = 0;
    double maxTargetPressureJump = 0.0;
    double maxCorrectionPressureJump = 0.0;
    std::int64_t structureRebuilds = 0;
    std::int64_t linearSolves = 0;
};

/// @brief 持久压力校正器；跨时间步复用 HYPRE IJ/ParCSR 结构。
class Corrector {
public:
    using InterfaceJumpProvider = std::function<
        const FDM::IInterfaceJumpCondition*(const Field&)>;

    Corrector(FDM::BoundaryConfig boundaries,
              FDM::PressureCorrectionConfig config,
              double idealGasGamma);

    void setExecutionRuntime(FDM::IExecutionRuntime* runtime) {
        runtime_ = runtime;
    }
    void setInterfaceJumpProvider(InterfaceJumpProvider provider) {
        interfaceJumpProvider_ = std::move(provider);
    }

    /// @brief 求解压力修正并同步更新压力、动量和总能量。
    CorrectionSummary correct(Field& field, double dt);
    CorrectionSummary correct(const std::vector<Field*>& fields, double dt);
    const LinearAlgebra::ReuseStatistics& reuseStatistics() const {
        return linearSolver_.statistics();
    }

private:
    FDM::BoundaryConfig boundaries_;
    FDM::PressureCorrectionConfig config_;
    double idealGasGamma_ = 1.4;
    LinearAlgebra::SolverSession linearSolver_;
    FDM::IExecutionRuntime* runtime_ = nullptr;
    InterfaceJumpProvider interfaceJumpProvider_;
};

} // namespace SF::PressureBased
