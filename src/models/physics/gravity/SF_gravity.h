/// @file SF_gravity.h
/// @brief 重力源项物理模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include "core/residual/SF_residual.h"

namespace SF {
namespace Source {
namespace Gravity {

inline void addTranslation(Field& field, Residual& residual,
                           int i, int j, int k, const Vector3& value) {
    const auto state = field.hasStateModel()
        ? field.thermodynamicState(i,j,k)
        : Physics::FluidStateModel::ThermodynamicState{};
    double rho = field.hasStateModel() ? state.density : field(i, j, k, RHO);
    const int momentum = field.hasStateModel()
        ? field.stateModel()->momentumIndex(0) : RU;
    const int energy = field.hasStateModel()
        ? field.stateModel()->energyIndex() : E;

    double u = field(i, j, k, momentum) / rho;
    double v = field(i, j, k, momentum+1) / rho;
    double w = field(i, j, k, momentum+2) / rho;

    residual.source(i, j, k, momentum) += rho * value.x;
    residual.source(i, j, k, momentum+1) += rho * value.y;
    residual.source(i, j, k, momentum+2) += rho * value.z;
    residual.source(i, j, k, energy) +=
        rho * (u * value.x + v * value.y + w * value.z);
}

} // namespace Gravity
} // namespace Source
} // namespace SF
