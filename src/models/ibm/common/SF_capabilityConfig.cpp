/// @file SF_capabilityConfig.cpp
/// @brief IBM capability 快照构造实现。

#include "common/SF_capabilityConfig.h"

namespace SF::IBM::Common {

FDM::ImmersedMethodCapabilities configureCapabilities(
        const FDM::ImmersedMethodSelection& selection) {
    FDM::ImmersedMethodCapabilities capabilities;
    capabilities.surfaceConstraint = selection.hasSurfaceConstraint();
    capabilities.bodyConstraint = selection.hasBodyConstraint();
    capabilities.diffuseTransfer =
        selection.representation == FDM::IBMRepresentation::DiffuseKernel;
    capabilities.eulerianMask =
        selection.representation == FDM::IBMRepresentation::EulerianMask;
    capabilities.sharpJump =
        selection.representation == FDM::IBMRepresentation::SharpJump;
    capabilities.prescribedSolid =
        selection.solid == FDM::IBMSolidModel::Prescribed;
    capabilities.coupledRigid =
        selection.solid == FDM::IBMSolidModel::CoupledRigid;
    capabilities.selfPropelledRigid =
        selection.solid == FDM::IBMSolidModel::SelfPropelledRigid;
    capabilities.densityBased =
        selection.method == FDM::IBMMethod::Ghost
        || selection.enforcement != FDM::IBMEnforcement::MonolithicKKT;
    capabilities.pressureBased = true;
    return capabilities;
}

} // namespace SF::IBM::Common
