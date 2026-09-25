/// @file SF_IBM.h
/// @brief IBM 对 application/solver 暴露的薄 facade。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include "SF_ibmConfig.h"
#include "SF_immersedSystem.h"
#include "geoProcessing/SF_STLGeometry.h"
#include "method/ghost/SF_ghost.h"
#include "method/SF_method.h"

#include <string>
#include <vector>

namespace SF::IBM {

/// @brief 统一 IBM 组合根；只负责生命周期和对外端口转发。
///
/// 方法选择只有一个权威来源：`IBMRuntimeConfig::method`。几何加载、方法装配、
/// capability 构造和启动日志分别位于 `common/`，避免 facade 再次膨胀。
class IB final : public FDM::IImmersedSystem {
public:
    /// @brief 从冻结的运行配置装配 IBM。
    bool setup(Field& field,
               const std::string& caseDir,
               const std::vector<std::string>& stlFiles,
               const IBMRuntimeConfig& config);

    bool active() const { return active_; }
    bool usesGhostCells() const;
    bool usesForcing() const;

    const FDM::ImmersedMethodSelection& methodSelection() const override;
    const FDM::ImmersedMethodCapabilities& capabilities() const override;
    const FDM::ImmersedAlgorithmDescriptor& algorithmDescriptor() const override;
    const std::vector<FDM::ImmersedFluidPort>& fluidPorts() const override;

    /// @brief 在指定物理时刻刷新 ghost 状态。
    /// @param field 待闭合的守恒场。
    /// @param time 当前 RK/时间推进 stage 的物理时刻，单位 s。
    /// @param dt 当前物理步长，单位 s；初始化阶段允许为 0。
    void applyGhostCells(Field& field, double time, double dt);

    FDM::ImmersedConstraintResult applyConstraint(
        const std::vector<Field*>& fields, double targetTime, double dt);
    /// @brief 注入 solver 提供的 Runtime；forcing 数学核不直接接触 MPI。
    void setExecutionRuntime(FDM::IExecutionRuntime* runtime);
    const FDM::ImmersedSurfaceSystem& prepareMonolithicSystem(
        Field& field, double targetTime, double dt);
    FDM::ImmersedConstraintResult acceptMonolithicSolution(
        Field& field, double targetTime, double dt,
        const FDM::ImmersedKKTState& state);

    double constraintMask(const Field& field, int i, int j, int k) const;
    Vector3 multiplier(const Field& field, int i, int j, int k) const;

private:
    bool active_ = false;
    GeoProcessing::STLGeometry geometry_;
    GhostIBM::WeightBuilder ghostWeights_;
    GhostIBM::GhostCellIBM ghostCell_;
    Forcing::ImmersedForcingSystem forcing_;
    IBMRuntimeConfig config_;
    FDM::ImmersedMethodSelection selection_;
    FDM::ImmersedMethodCapabilities capabilities_;
    FDM::ImmersedAlgorithmDescriptor descriptor_;
    std::vector<FDM::ImmersedFluidPort> fluidPorts_;
};

} // namespace SF::IBM
