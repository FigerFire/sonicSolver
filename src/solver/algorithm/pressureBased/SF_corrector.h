#pragma once

/// @file SF_corrector.h
/// @brief 单流体压力基 HYPRE 压力—速度修正。

#include "SF_config.h"
#include "SF_field.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "core/interfaces/SF_equationCoupling.h"
#include "solver/linearAlgebra/SF_linearAlgebra.h"

#include <functional>
#include <memory>

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
    ~Corrector();

    void setExecutionRuntime(FDM::IExecutionRuntime* runtime) {
        runtime_ = runtime;
    }
    void setInterfaceJumpProvider(InterfaceJumpProvider provider) {
        interfaceJumpProvider_ = std::move(provider);
    }

    /// @brief 求解压力修正并同步更新压力、动量和总能量。
    CorrectionSummary correct(Field& field, double dt);
    CorrectionSummary correct(const std::vector<Field*>& fields, double dt);

    /// @brief Typed operations used by the structured PISO plan.
    void assemble(const std::vector<Field*>& fields, double dt);
    void solve();
    void preparePressureUpdate();
    void correctVelocity();
    void correctFlux();
    CorrectionSummary commitPressureUpdate();
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
    struct Workspace;
    std::unique_ptr<Workspace> workspace_;
};

} // namespace SF::PressureBased
