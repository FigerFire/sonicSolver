/// @file SF_algorithmDescriptor.cpp
/// @brief Frozen IBM options; all mathematical inventories live in SystemContribution.
#include "descriptor/SF_algorithmDescriptor.h"
#include "SF_config.h"
#include <stdexcept>

namespace SF::IBM::Descriptor {
using FDM::ImmersedAlgorithmDescriptor;
ImmersedAlgorithmDescriptor ghostCell() {
    ImmersedAlgorithmDescriptor result;
    result.id="ghostCellIBM";result.referenceName="sharp-interface GhostCellIBM";
    result.support=FDM::IBMConstraintSupport::Surface;
    result.representation=FDM::IBMRepresentation::SharpJump;
    result.enforcement=FDM::IBMEnforcement::GhostCell;
    validate(result);
    return result;
}
ImmersedAlgorithmDescriptor variational(const FDM::IBMForcingConfig& config) {
    FDM::validateIBMForcingSelection(config);
    if (config.constraintSupport==FDM::IBMConstraintSupport::SurfaceAndBody)
        throw std::runtime_error("One IBM port cannot merge surface and body multipliers; declare separate constraint blocks.");
    ImmersedAlgorithmDescriptor result;
    result.id=FDM::toString(config.algorithm);result.algorithm=config.algorithm;
    result.referenceName="Bhalla unified IBM formulation";
    result.support=config.constraintSupport;result.representation=config.representation;
    result.enforcement=config.enforcement;result.solid=config.solidModel;
    result.rigidMotionMode=config.rigidMotionMode;
    result.surfaceNormalization=config.surfaceNormalization;
    result.monolithic=config.enforcement==FDM::IBMEnforcement::MonolithicKKT;
    const bool lagged=config.enforcement==FDM::IBMEnforcement::ExplicitIBM;
    const bool penalty=config.enforcement==FDM::IBMEnforcement::BrinkmanPenalty;
    result.variational.fluidKineticIncrement=!lagged && !penalty;
    result.variational.incompressibilityConstraint=result.monolithic;
    result.variational.immersedNoSlipConstraint=!lagged && !penalty;
    if (lagged) result.variational.stationaryFunctional=
        "consume lagged force; diagonal multiplier response for the next physical step";
    else if (penalty) result.variational.stationaryFunctional=
        "post-predictor momentum penalty and midpoint mechanical work; no multiplier constraint";
    else if (!result.monolithic && config.solidModel==FDM::IBMSolidModel::SelfPropelledRigid)
        result.variational.stationaryFunctional=
            "virtual-fluid generalized mass projection followed by pointwise velocity projection";
    else result.variational.stationaryFunctional=result.monolithic
        ?"fluid/solid kinetic increment - solid external work + p^T D u + Lambda^T(Ju-Gq-uDef)"
        :"minimum mass-norm velocity increment subject to Ju=Us";
    if (config.algorithm==FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian)
        result.variational.stationaryFunctional+=" + gamma/2 ||Ju-Gq-uDef||^2_ML";
    validate(result);
    return result;
}
void validate(const ImmersedAlgorithmDescriptor& descriptor) {
    if (descriptor.wallClosure!=FDM::ImmersedWallClosure::EulerSlip
        && descriptor.enforcement!=FDM::IBMEnforcement::GhostCell)
        throw std::runtime_error("Viscous wall closure requires a Ghost boundary implementation, not a forcing selection.");
    if (descriptor.id.empty() || descriptor.referenceName.empty())
        throw std::runtime_error("IBM implementation options require id and referenceName.");
    if (descriptor.monolithic && (!descriptor.variational.fluidKineticIncrement
        || !descriptor.variational.incompressibilityConstraint || !descriptor.variational.immersedNoSlipConstraint))
        throw std::runtime_error("Monolithic IBM must declare kinetic, pressure and no-slip capabilities.");
}
} // namespace SF::IBM::Descriptor
