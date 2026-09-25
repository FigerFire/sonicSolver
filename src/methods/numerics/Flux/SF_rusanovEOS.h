#pragma once

/// @file SF_rusanovEOS.h
/// @brief FluidStateModel/EOS-aware Rusanov 守恒通量基线。

#include "SF_fluidStateModel.h"
#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"

#include <vector>

namespace SF::Numerics::RusanovEOS {

/// @brief 计算任意 FluidStateModel 的局部 Lax-Friedrichs 面通量。
std::vector<double> faceFlux(
    const Physics::FluidStateModel::Model& equations,
    const double* left,
    const double* right,
    const std::array<double, 3>& unitNormal);

/// @brief 将 FluidStateModel-aware Rusanov 面通量写入调用者提供的连续缓存。
/// @param resultCount 必须等于 FluidStateModel 的变量数；不进行数值降级。
void faceFlux(
    const Physics::FluidStateModel::Model& equations,
    const double* left,
    const double* right,
    const std::array<double, 3>& unitNormal,
    double* result,
    int resultCount);

/// @brief 一阶重构基线：所有压力和声速均来自 FluidStateModel closure。
void div(Field& field, FluxField& fluxField, Residual& residual,
         const Physics::FluidStateModel::Model& equations);

/// @brief 使用 FluidStateModel 返回的声速计算 CFL 时间步。
double deltaT(const Field& field,
              const Physics::FluidStateModel::Model& equations,
              double cfl);

} // namespace SF::Numerics::RusanovEOS
