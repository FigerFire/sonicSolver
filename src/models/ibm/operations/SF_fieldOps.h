#pragma once

/// @file SF_fieldOps.h
/// @brief IBM 对 canonical Field 的无状态读取、校验和写入操作。

#include "SF_fluidStateModel.h"
#include "SF_field.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::FieldOps {

/// @brief 三维结构索引值对象。
struct Index3 {
    int i = 0;
    int j = 0;
    int k = 0;
};

inline Index3 index3(const Field& field, int localCell) {
    if (localCell < 0 || localCell >= field.TotalSize()) {
        throw std::runtime_error("IBM local cell index is outside Field storage.");
    }
    Index3 result;
    field.getIJK(localCell,result.i,result.j,result.k);
    return result;
}

inline int momentumIndex(const Field& field) {
    return field.hasStateModel()
        ? field.stateModel()->momentumIndex(0) : RU;
}

inline int energyIndex(const Field& field) {
    return field.hasStateModel()
        ? field.stateModel()->energyIndex() : E;
}

inline double density(const Field& field, const Index3& index) {
    return field.hasStateModel()
        ? field.thermodynamicState(index.i,index.j,index.k).density
        : field(index.i,index.j,index.k,RHO);
}

inline double volume(const Field& field, const Index3& index) {
    const double inverse = field.Jac(index.i,index.j,index.k);
    if (!std::isfinite(inverse) || inverse == 0.0) {
        throw std::runtime_error("IBM method found an invalid cell Jacobian.");
    }
    // 强迫/约束路径只对已由 IBM 拓扑筛选的 Eulerian unknown 调用本函数。
    // 因此它必须和 CouplingGraph 的 canonical 体积一致，不能把本地 patch
    // 端点当作物理边界再乘半权重；MPI 切分面也是本地端点。
    return 1.0/std::abs(inverse);
}

inline bool interior(const Field& field, const Index3& index) {
    const int ng = field.NG();
    return index.i >= ng && index.i < ng+field.NX()
        && index.j >= ng && index.j < ng+field.NY()
        && index.k >= ng && index.k < ng+field.NZ();
}

inline void validateState(
        const Field& field,
        int momentum,
        int energy,
        const Index3& index) {
    const double rho = density(field,index);
    if (!std::isfinite(rho) || rho <= 0.0) {
        throw std::runtime_error("IBM method found a non-positive density.");
    }
    if (!std::isfinite(field(index.i,index.j,index.k,momentum))
        || !std::isfinite(field(index.i,index.j,index.k,momentum+1))
        || !std::isfinite(field(index.i,index.j,index.k,momentum+2))
        || !std::isfinite(field(index.i,index.j,index.k,energy))) {
        throw std::runtime_error("IBM predictor state is non-finite.");
    }
}

inline Vector3 velocity(
        const Field& field,
        int momentum,
        const Index3& index) {
    const double rho = density(field,index);
    return {field(index.i,index.j,index.k,momentum)/rho,
            field(index.i,index.j,index.k,momentum+1)/rho,
            field(index.i,index.j,index.k,momentum+2)/rho};
}

inline Vector3 momentum(
        const Field& field,
        int momentumIndex,
        const Index3& index) {
    return {field(index.i,index.j,index.k,momentumIndex),
            field(index.i,index.j,index.k,momentumIndex+1),
            field(index.i,index.j,index.k,momentumIndex+2)};
}

inline void setMomentum(
        Field& field,
        int momentumIndex,
        const Index3& index,
        const Vector3& value) {
    field(index.i,index.j,index.k,momentumIndex) = value.x;
    field(index.i,index.j,index.k,momentumIndex+1) = value.y;
    field(index.i,index.j,index.k,momentumIndex+2) = value.z;
}

} // namespace SF::IBM::FieldOps
