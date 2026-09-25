/// @file SF_constraintOps.cpp
/// @brief `ibm.constraint.project` 的现有数值实现与 provider 校验。

#include "solver/algorithm/immersed/SF_constraintOps.h"
#include "solver/algorithm/immersed/SF_immersedStrategy.h"

#include <cmath>
#include <stdexcept>

namespace SF::ImmersedAlgorithm {

FDM::ImmersedConstraintResult projectConstraint(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt,
        FDM::IExecutionRuntime* runtime,
        const FDM::ImmersedCouplingPorts& immersed) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "ibm.constraint.project requires finite positive dt.");
    }
    if (!immersed.constraint || !immersed.system) {
        throw std::runtime_error(
            "ibm.constraint.project requires matching constraint and "
            "immersed-system providers.");
    }
    if (const auto* provider = immersed.constraint->systemProvider();
        provider && provider != immersed.system) {
        throw std::runtime_error(
            "ibm.constraint.project received providers from different "
            "immersed systems.");
    }
    validateProjectionProvider(*immersed.system);
    immersed.constraint->setExecutionRuntime(runtime);
    if (fields.size() != 1
        && !immersed.system->capabilities().distributed) {
        throw std::runtime_error(
            "ibm.constraint.project received multiple local patches, but "
            "canonical distributed constraint ownership is unavailable.");
    }
    return immersed.constraint->projectPredictedState(
        fields,targetTime,dt);
}

} // namespace SF::ImmersedAlgorithm
