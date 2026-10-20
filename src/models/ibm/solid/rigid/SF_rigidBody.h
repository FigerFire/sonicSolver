#pragma once

/// @file SF_rigidBody.h
/// @brief 表面变分 IBM 的平动/转动刚体状态与构型更新。

#include "solid/SF_bodyModel.h"

#include <array>
#include <vector>

namespace SF::IBM::Forcing {

/// @brief 管理参考 STL 到当前刚体构型的映射及速度状态。
class RigidBodyState final : public IBodyModel {
public:
    void configure(const FDM::IBMForcingConfig& config) override;
    std::vector<GeoProcessing::Triangle> triangles(
        const GeoProcessing::STLGeometry& geometry,
        double targetTime) const override;
    /// @brief 把当前世界坐标点逆映射到 STL 参考构型。
    GeoProcessing::Point referencePoint(
        const Vector3& world, double targetTime) const override;
    /// @brief 返回指定世界坐标点的刚体速度。
    Vector3 velocityAt(
        const Vector3& point, double targetTime) const override;
    Vector3 deformationVelocityAt(
        const Vector3& point, double targetTime) const override;
    /// @brief 写入表面点 G 基向量以及本时间步固体方程 A_s/r_s。
    void prepareEquationView(
        FDM::ImmersedSurfaceSystem& surface,
        double targetTime,
        double dt) const override;
    void advanceCoupled(const Vector3& linearVelocity,
                        const Vector3& angularVelocity,
                        double dt) override;

    Vector3 center(double targetTime) const override;
    std::array<double,9> worldInertia(double targetTime) const;
    const Vector3& linearVelocity() const override { return linearVelocity_; }
    const Vector3& angularVelocity() const override { return angularVelocity_; }

private:
    FDM::IBMForcingConfig config_;
    Vector3 center_;
    Vector3 linearVelocity_;
    Vector3 angularVelocity_;
    std::array<double,9> orientation_{{1,0,0,0,1,0,0,0,1}};

    Vector3 transformed(const Vector3& reference, double targetTime) const;
};

} // namespace SF::IBM::Forcing
