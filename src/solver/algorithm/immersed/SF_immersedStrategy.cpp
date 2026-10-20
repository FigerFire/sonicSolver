/// @file SF_immersedStrategy.cpp
/// @brief 把 ghost IBM 或 forcing IBM 接入算法阶段而不暴露具体方法。

#include "solver/algorithm/immersed/SF_immersedStrategy.h"

#include <stdexcept>

namespace SF::ImmersedAlgorithm {

StrategyContract contract(FDM::IBMEnforcement enforcement) {
    StrategyContract result;
    result.kind = enforcement;
    switch (enforcement) {
        case FDM::IBMEnforcement::GhostCell:
            result.name = "GhostCellIBM";
            result.implemented = true;
            break;
        case FDM::IBMEnforcement::ExplicitIBM:
            result.name = "ImmersedExplicitIBM";
            result.requiresPredictedState = true;
            result.producesMultiplier = true;
            result.implemented = true;
            break;
        case FDM::IBMEnforcement::FractionalDLM:
            result.name = "ImmersedFractionalDLM";
            result.requiresPredictedState = true;
            result.producesMultiplier = true;
            result.implemented = true;
            break;
        case FDM::IBMEnforcement::VelocityForcing:
            result.name = "ImmersedVelocityForcing";
            result.requiresPredictedState = true;
            result.producesMultiplier = true;
            result.implemented = true;
            break;
        case FDM::IBMEnforcement::BrinkmanPenalty:
            result.name = "ImmersedBrinkmanPenalty";
            result.requiresPredictedState = true;
            result.implemented = true;
            break;
        case FDM::IBMEnforcement::MonolithicKKT:
            result.name = "ImmersedMonolithicKKT";
            result.requiresPressureCorrection = true;
            result.producesMultiplier = true;
            result.implemented = true;
            break;
    }
    if (result.name.empty()) {
        throw std::runtime_error("Unknown immersed enforcement strategy.");
    }
    return result;
}

namespace {
void validateSystem(const FDM::IImmersedSystem& system) {
    const auto& selection = system.methodSelection();
    const auto& capabilities = system.capabilities();
    const StrategyContract selected = contract(selection.enforcement);
    if (!selected.implemented) {
        throw std::runtime_error(
            "IBM enforcement '" + selected.name
            + "' is declared but its equation assembly is not implemented "
              "in this build.");
    }
    if (selection.method != FDM::IBMMethod::VariationalForcing) {
        throw std::runtime_error(
            "IBM constraint operation received a non-variational family; "
            "ghost-cell IBM must stay in the boundary pipeline.");
    }

    if (selection.hasSurfaceConstraint()
        && !capabilities.surfaceConstraint) {
        throw std::runtime_error(
            "Selected IBM surface constraint is not provided by the "
            "configured immersed system.");
    }
    if (selection.hasBodyConstraint() && !capabilities.bodyConstraint) {
        throw std::runtime_error(
            "Selected IBM body constraint is not provided by the "
            "configured immersed system.");
    }
    if (selection.representation == FDM::IBMRepresentation::DiffuseKernel
        && !capabilities.diffuseTransfer) {
        throw std::runtime_error(
            "Selected IBM diffuseKernel transfer is not available.");
    }
    if (selection.representation == FDM::IBMRepresentation::EulerianMask
        && !capabilities.eulerianMask) {
        throw std::runtime_error(
            "Selected IBM eulerianMask transfer is not available.");
    }
    if (selection.representation == FDM::IBMRepresentation::SharpJump
        && !capabilities.sharpJump) {
        throw std::runtime_error(
            "Selected IBM sharpJump transfer is not available.");
    }
    if (selection.solid == FDM::IBMSolidModel::Prescribed
        && !capabilities.prescribedSolid) {
        throw std::runtime_error(
            "Selected prescribed IBM solid model is not available.");
    }
    if (selection.solid == FDM::IBMSolidModel::CoupledRigid
        && !capabilities.coupledRigid) {
        throw std::runtime_error(
            "Selected coupledRigid IBM solid model is not available.");
    }
    if (selection.solid == FDM::IBMSolidModel::SelfPropelledRigid
        && !capabilities.selfPropelledRigid) {
        throw std::runtime_error(
            "Selected selfPropelledRigid IBM model is not available.");
    }
    if (selection.solid == FDM::IBMSolidModel::Deformable) {
        throw std::runtime_error(
            "Selected deformable IBM solid has no migrated equation block.");
    }
    bool phasePort = false;
    for (const auto& port : system.fluidPorts()) {
        phasePort = phasePort || port.phaseIndex >= 0;
    }
    if (phasePort && !capabilities.multiPhase) {
        throw std::runtime_error(
            "Phase-wise IBM fluidPorts were selected, but multiphase "
            "constraint assembly is not available.");
    }
}
} // namespace

void validateProjectionProvider(const FDM::IImmersedSystem& system) {
    validateSystem(system);
    const StrategyContract selected = contract(
        system.methodSelection().enforcement);
    if (selected.requiresPressureCorrection) {
        throw std::runtime_error(
            "ibm.constraint.project cannot execute enforcement '"
            +selected.name+"'; it requires ibm.kkt.solve.");
    }
    if (!selected.requiresPredictedState) {
        throw std::runtime_error(
            "ibm.constraint.project requires a predicted-state constraint "
            "enforcement provider.");
    }
}

void validateMonolithicProvider(const FDM::IImmersedSystem& system) {
    validateSystem(system);
    const StrategyContract selected = contract(
        system.methodSelection().enforcement);
    if (!selected.requiresPressureCorrection
        || system.methodSelection().enforcement
            != FDM::IBMEnforcement::MonolithicKKT) {
        throw std::runtime_error(
            "ibm.kkt.solve requires a monolithic KKT enforcement provider.");
    }
}

} // namespace SF::ImmersedAlgorithm
