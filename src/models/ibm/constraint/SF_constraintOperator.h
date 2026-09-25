#pragma once

/// @file SF_constraintOperator.h
/// @brief IBM 约束拓扑构造的可替换窄接口。

#include "SF_field.h"
#include "SF_immersedConstraint.h"
#include "geoProcessing/SF_STLGeometry.h"

#include <vector>

namespace SF::IBM::Forcing {

class IBodyModel;

/// @brief body mask 约束拓扑构造接口。
class IBodyConstraintOperator {
public:
    virtual ~IBodyConstraintOperator() = default;
    virtual const FDM::ImmersedBodySystem& build(
        const Field& field,
        const GeoProcessing::STLGeometry& geometry,
        const IBodyModel& bodyModel,
        double targetTime) = 0;
};

/// @brief 表面 Lagrangian 点与 Eulerian raw edge 构造接口。
class ISurfaceConstraintOperator {
public:
    virtual ~ISurfaceConstraintOperator() = default;
    virtual const FDM::ImmersedSurfaceSystem& build(
        const Field& field,
        const std::vector<GeoProcessing::Triangle>& triangles,
        const Vector3& center,
        double supportRadius) = 0;

    /// @brief 返回本地 Eulerian owner edge 的每-marker raw 权重和。
    virtual const std::vector<double>& localNormalizations() const {
        throw std::runtime_error(
            "Surface constraint operator does not expose distributed normalization.");
    }
    /// @brief 使用 Runtime 归并后的权重和刷新 canonical J 行。
    virtual const FDM::ImmersedSurfaceSystem& normalizeDistributed(
        const std::vector<double>&) {
        throw std::runtime_error(
            "Surface constraint operator does not support distributed normalization.");
    }
};

} // namespace SF::IBM::Forcing
