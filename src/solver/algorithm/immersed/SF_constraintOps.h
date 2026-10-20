#pragma once

/// @file SF_constraintOps.h
/// @brief IBM constraint Plan operation 的数值入口。

#include "core/interfaces/SF_boundaryPipeline.h"
#include "core/interfaces/SF_executionRuntime.h"

#include <vector>

namespace SF::ImmersedAlgorithm {

/// @brief 对 predictor state 执行一次既有 variational constraint projection。
FDM::ImmersedConstraintResult projectConstraint(
    const std::vector<Field*>& fields,
    double targetTime,
    double dt,
    FDM::IExecutionRuntime* runtime,
    const FDM::ImmersedCouplingPorts& immersed);

} // namespace SF::ImmersedAlgorithm
