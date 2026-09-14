#pragma once

/// @file SF_transport.h
/// @brief Eulerian 守恒标量输运方程的统一矩阵装配接口。

#include "solver/linearAlgebra/SF_linearAlgebra.h"
#include "solver/linearAlgebra/SF_distributedRowMap.h"
#include "solver/equation/eulerian/SF_workspace.h"
#include "solver/discretization/eulerian/SF_eulerian.h"

namespace SF::EulerianEulerian {

/// @brief 装配瞬态、迎风对流、扩散、显式源和非负隐式汇。
LinearAlgebra::SparseSystem assembleTransportEquation(
    const Field& field,
    const LinearAlgebra::DistributedRowMap& rowMap,
    const ScalarField& unknown,
    const ScalarField& previousConserved,
    const ScalarField& capacity,
    const ScalarField& diffusivity,
    const ScalarField& explicitSource,
    const std::vector<double>& canonicalMassFlux,
    double timeStep,
    bool axisymmetric,
    int radialCoordinate,
    ScalarField& diagonal,
    const ScalarField* implicitSink = nullptr);

} // namespace SF::EulerianEulerian
