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

void validateForAlgorithm(const FDM::IBMForcingConfig& forcing,
                          FDM::SolverAlgorithm algorithm) {
    FDM::validateIBMForcingSelection(forcing);
    FDM::ImmersedMethodSelection selection;
    selection.method = FDM::IBMMethod::VariationalForcing;
    selection.support = forcing.constraintSupport;
    selection.representation = forcing.representation;
    selection.enforcement = forcing.enforcement;
    selection.solid = forcing.solidModel;
    validateForAlgorithm(selection, algorithm);
}

void validateForAlgorithm(const FDM::ImmersedMethodSelection& selection,
                          FDM::SolverAlgorithm algorithm) {
    const StrategyContract selected = contract(selection.enforcement);
    if (selected.requiresPressureCorrection
        && algorithm != FDM::SolverAlgorithm::PressureBased) {
        throw std::runtime_error(
            "IBM enforcement '" + selected.name
            + "' requires a pressure-based pressure-correction stage.");
    }
    // 顺序型约束只要求一个已经完成流体 predictor/corrector 的状态；该数学
    // 契约与 densityBased/pressureBased 无关。具体算法负责把它放在明确阶段，
    // IBM 模块不再把求解器类型当作方法选择轴。
    if (!selected.implemented) {
        throw std::runtime_error(
            "IBM enforcement '" + selected.name
            + "' is declared but its equation assembly is not implemented "
              "in this build.");
    }
}

void validateForAlgorithm(const FDM::IImmersedSystem& system,
                          FDM::SolverAlgorithm algorithm) {
    const auto& selection = system.methodSelection();
    const auto& capabilities = system.capabilities();
    if (selection.method != FDM::IBMMethod::VariationalForcing) {
        throw std::runtime_error(
            "Flow correction received a non-variational IBM family; "
            "ghost-cell IBM must stay in the boundary pipeline.");
    }
    validateForAlgorithm(selection, algorithm);

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
    if (algorithm == FDM::SolverAlgorithm::DensityBased
        && !capabilities.densityBased) {
        throw std::runtime_error(
            "Configured immersed system does not provide a density-based "
            "enforcement stage.");
    }
    if (algorithm == FDM::SolverAlgorithm::PressureBased
        && !capabilities.pressureBased) {
        throw std::runtime_error(
            "Configured immersed system does not provide a pressure-based "
            "enforcement stage.");
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

} // namespace SF::ImmersedAlgorithm
