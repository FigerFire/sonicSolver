#pragma once

/// @file SF_bodyModel.h
/// @brief IBM 固体运动与动力学的求解器中立对象接口。

#include "SF_configTypes.h"
#include "SF_immersedConstraint.h"
#include "geoProcessing/SF_STLGeometry.h"

#include <vector>

namespace SF::IBM::Forcing {

/// @brief 给 IBM 方法提供构型、边界速度和固体方程视图。
class IBodyModel {
public:
    virtual ~IBodyModel() = default;
    virtual void configure(const FDM::IBMForcingConfig& config) = 0;
    virtual std::vector<GeoProcessing::Triangle> triangles(
        const GeoProcessing::STLGeometry& geometry,
        double targetTime) const = 0;
    virtual GeoProcessing::Point referencePoint(
        const Vector3& world, double targetTime) const = 0;
    virtual Vector3 velocityAt(
        const Vector3& point, double targetTime) const = 0;
    /// @brief 返回不包含刚体平动/转动的局部形变速度。
    ///
    /// 刚体实现严格返回零；可变形/游动体实现可据此把边界速度写成
    /// `G q + uDef`，供 FTS 与全隐式 KKT 共用同一仿射约束。
    virtual Vector3 deformationVelocityAt(
        const Vector3& point, double targetTime) const = 0;
    virtual void prepareEquationView(
        FDM::ImmersedSurfaceSystem& surface,
        double targetTime, double dt) const = 0;
    virtual void advanceCoupled(
        const Vector3& linearVelocity,
        const Vector3& angularVelocity,
        double dt) = 0;
    virtual Vector3 center(double targetTime) const = 0;
    virtual const Vector3& linearVelocity() const = 0;
    virtual const Vector3& angularVelocity() const = 0;
};

} // namespace SF::IBM::Forcing
