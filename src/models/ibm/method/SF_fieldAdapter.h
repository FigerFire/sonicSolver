#pragma once

/// @file SF_fieldAdapter.h
/// @brief IBM 方法访问 canonical Field 状态的内部适配器。

#include "SF_fluidStateModel.h"
#include "SF_field.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Forcing::FieldAdapter {

inline bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

inline int momentumIndex(const Field& field) {
    return field.hasStateModel()
        ? field.stateModel()->momentumIndex(0) : RU;
}

inline int energyIndex(const Field& field) {
    return field.hasStateModel()
        ? field.stateModel()->energyIndex() : E;
}

inline double density(const Field& field, int i, int j, int k) {
    return field.hasStateModel()
        ? field.thermodynamicState(i,j,k).density
        : field(i,j,k,RHO);
}

inline double volume(const Field& field, int i, int j, int k) {
    const double inverse = field.Jac(i,j,k);
    if (!std::isfinite(inverse) || inverse == 0.0) {
        throw std::runtime_error("IBM method found an invalid cell Jacobian.");
    }
    // 与 FieldOps/CouplingGraph 保持同一个 canonical Eulerian 体积定义。
    // 物理边界点不进入 Peskin 的表面插值边，因此不能以 local patch 边界
    // 决定半对偶体积，避免 MPI split 改变 IBM 力和功。
    return 1.0/std::abs(inverse);
}

inline bool interior(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return i >= ng && i < ng+field.NX()
        && j >= ng && j < ng+field.NY()
        && k >= ng && k < ng+field.NZ();
}

inline void validateState(const Field& field, int momentum, int energy,
                          int i, int j, int k) {
    const double rho = density(field,i,j,k);
    if (!std::isfinite(rho) || rho <= 0.0) {
        throw std::runtime_error("IBM method found a non-positive density.");
    }
    if (!std::isfinite(field(i,j,k,momentum))
        || !std::isfinite(field(i,j,k,momentum+1))
        || !std::isfinite(field(i,j,k,momentum+2))
        || !std::isfinite(field(i,j,k,energy))) {
        throw std::runtime_error("IBM predictor state is non-finite.");
    }
}

inline Vector3 velocity(const Field& field, int momentum,
                        int i, int j, int k) {
    const double rho = density(field,i,j,k);
    return {field(i,j,k,momentum)/rho,
            field(i,j,k,momentum+1)/rho,
            field(i,j,k,momentum+2)/rho};
}

inline void setMomentum(Field& field, int momentum,
                        int i, int j, int k, const Vector3& value) {
    field(i,j,k,momentum) = value.x;
    field(i,j,k,momentum+1) = value.y;
    field(i,j,k,momentum+2) = value.z;
}

} // namespace SF::IBM::Forcing::FieldAdapter
