/// @file SF_method.cpp
/// @brief 统一 IBM forcing facade 的配置校验、算法分派和依赖注入。

#include "method/SF_method.h"
#include "method/SF_fieldAdapter.h"
#include "descriptor/SF_algorithmDescriptor.h"

#include "SF_config.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace SF::IBM::Forcing {

ImmersedForcingSystem::ImmersedForcingSystem()
    : surfaceOperator_(std::make_unique<SurfaceOperator>()),
      bodyOperator_(std::make_unique<BodyOperator>()),
      bodyModel_(std::make_unique<RigidBodyState>()) {}

void ImmersedForcingSystem::setBodyConstraintOperator(
        std::unique_ptr<IBodyConstraintOperator> operatorInstance) {
    if (!operatorInstance) {
        throw std::invalid_argument(
            "ImmersedForcingSystem cannot accept a null body operator.");
    }
    bodyOperator_ = std::move(operatorInstance);
}

void ImmersedForcingSystem::setSurfaceConstraintOperator(
        std::unique_ptr<ISurfaceConstraintOperator> operatorInstance) {
    if (!operatorInstance) {
        throw std::invalid_argument(
            "ImmersedForcingSystem cannot accept a null surface operator.");
    }
    surfaceOperator_ = std::move(operatorInstance);
}

void ImmersedForcingSystem::setBodyModel(
        std::unique_ptr<IBodyModel> bodyModel) {
    if (!bodyModel) {
        throw std::invalid_argument(
            "ImmersedForcingSystem cannot accept a null body model.");
    }
    bodyModel_ = std::move(bodyModel);
}

void ImmersedForcingSystem::configure(
        const GeoProcessing::STLGeometry& geometry,
        IBMRuntimeConfig config) {
    try {
        FDM::validateIBMForcingSelection(config.forcing);
    } catch (const std::exception& error) {
        throw std::runtime_error(
            std::string("IBM method selection is invalid: ") + error.what());
    }
    const bool body = config.forcing.constraintSupport
        == FDM::IBMConstraintSupport::Body;
    if (config.forcing.constraintSupport
            == FDM::IBMConstraintSupport::SurfaceAndBody) {
        throw std::runtime_error(
            "surfaceAndBody requires two explicit constraint assemblies; "
            "this build accepts one support domain per method.");
    }
    if (!config.forcing.fluidPorts.empty()) {
        throw std::runtime_error(
            "Phase-wise IBM fluidPorts are parsed, but multiphase IBM "
            "assembly is not implemented.");
    }
    if (config.forcing.enforcement == FDM::IBMEnforcement::FractionalDLM
        && config.forcing.coupling
            != FDM::IBMConstraintCoupling::IncrementalProjection) {
        throw std::runtime_error(
            "fractionalDLM requires coupling incrementalProjection.");
    }
    if (config.forcing.enforcement == FDM::IBMEnforcement::MonolithicKKT
        && config.forcing.coupling
            != FDM::IBMConstraintCoupling::MonolithicKKT) {
        throw std::runtime_error(
            "fullyImplicitDLM requires coupling monolithicKKT.");
    }
    if (body && !geometry.watertight()) {
        throw std::runtime_error(
            "Body constraints require a watertight STL for inside/outside "
            "classification.");
    }
    if (config.forcing.energyCoupling
        != FDM::IBMEnergyCoupling::MechanicalWork) {
        throw std::runtime_error(
            "IBM forcing requires energyCoupling mechanicalWork.");
    }
    if (!FieldAdapter::finiteVector(config.forcing.centerOfMass)
        || !FieldAdapter::finiteVector(config.forcing.linearVelocity)
        || !FieldAdapter::finiteVector(config.forcing.angularVelocity)
        || !std::isfinite(config.forcing.constraintTolerance)
        || config.forcing.constraintTolerance <= 0.0) {
        throw std::runtime_error(
            "IBM method contains invalid rigid-body data or tolerance.");
    }
    geometry_ = &geometry;
    config_ = std::move(config);
    descriptor_ = Descriptor::variational(config_.forcing);
    laggedMultiplier_.clear();
    bodyModel_->configure(config_.forcing);
}

FDM::ImmersedConstraintResult ImmersedForcingSystem::projectPredictedState(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    switch (config_.forcing.algorithm) {
        case FDM::IBMForcingAlgorithm::PeskinOriginal:
            return applyPeskinOriginal(fields,targetTime,dt);
        case FDM::IBMForcingAlgorithm::DFMExplicitSelfPropelled:
            return applyDFMExplicitSelfPropelled(fields,targetTime,dt);
        case FDM::IBMForcingAlgorithm::DFMFractionalStepSelfPropelled:
            return applyDFMFractionalStepSelfPropelled(
                fields,targetTime,dt);
        case FDM::IBMForcingAlgorithm::DFMFractionalStepPrescribed:
            return applyDFMFractionalStepPrescribed(
                fields,targetTime,dt);
        case FDM::IBMForcingAlgorithm::VelocityForcingFTS:
            return applyVelocityForcingFTS(fields,targetTime,dt);
        case FDM::IBMForcingAlgorithm::VelocityForcingBP:
            return applyVelocityForcingBP(fields,targetTime,dt);
        case FDM::IBMForcingAlgorithm::DFMImplicitPrescribed:
        case FDM::IBMForcingAlgorithm::DFMImplicitSelfPropelled:
        case FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian:
            throw std::runtime_error(
                "implicit DFM must run in the monolithic pressure stage.");
    }
    throw std::runtime_error("Unknown IBM forcing algorithm.");
}

} // namespace SF::IBM::Forcing
