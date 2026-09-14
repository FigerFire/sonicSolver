/// @file SF_correction.cpp
/// @brief 密度基守恒状态校正与物理一致性检查。

#include "solver/algorithm/densityBased/SF_correction.h"
#include "solver/algorithm/immersed/SF_immersedStrategy.h"

#include <cmath>
#include <stdexcept>

namespace SF::DensityBased {

FDM::CapabilitySet Algorithm::capabilities() const {
    return {
        FDM::Capability::PredictedConservativeState,
        FDM::Capability::SingleFieldExecution,
        FDM::Capability::MultiFieldExecution};
}

FDM::FlowAlgorithmResult Algorithm::correct(
        FDM::FlowAlgorithmContext& context,
        double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "DensityBased::Algorithm requires finite positive dt.");
    }
    if (!context.immersed.constraint && !context.immersed.system) return {};
    if (context.immersed.system && !context.immersed.constraint) {
        throw std::runtime_error(
            "densityBased received an immersed system without its "
            "constraint numerical adapter.");
    }
    if (!context.immersed.system) {
        throw std::runtime_error(
            "densityBased received an immersed constraint adapter without "
            "the unified IImmersedSystem selection service.");
    }
    if (const auto* provider = context.immersed.constraint->systemProvider();
        provider && provider != context.immersed.system) {
        throw std::runtime_error(
            "densityBased received constraint and selection services from "
            "different IBM managers.");
    }
    ImmersedAlgorithm::validateForAlgorithm(
        *context.immersed.system,
        FDM::SolverAlgorithm::DensityBased);
    context.immersed.constraint->setExecutionRuntime(context.executionRuntime);
    if (context.fields.size() != 1
        && !context.immersed.system->capabilities().distributed) {
        throw std::runtime_error(
            "densityBased IBM forcing received multiple local patches, but "
            "canonical distributed constraint ownership is not available.");
    }
    const FDM::ImmersedConstraintResult result =
        context.immersed.constraint->projectPredictedState(
            context.fields, context.targetTime, dt);
    return {result.performed, result.performed, result.detail};
}

} // namespace SF::DensityBased
