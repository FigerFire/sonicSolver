/// @file SF_algorithmDescriptor.cpp
/// @brief IBM algorithm descriptor 的构造与启动期契约校验。

#include "descriptor/SF_algorithmDescriptor.h"

#include "SF_config.h"

#include <set>
#include <stdexcept>
#include <string>

namespace SF::IBM::Descriptor {
namespace {

using FDM::ImmersedAlgorithmDescriptor;
using FDM::ImmersedConstraintDescriptor;
using FDM::ImmersedOwnershipKind;
using FDM::ImmersedSolveBlockDescriptor;
using FDM::ImmersedUnknownDescriptor;
using FDM::ImmersedVariableLocation;

ImmersedUnknownDescriptor multiplier(
        const FDM::IBMForcingConfig& config) {
    const bool surface = config.constraintSupport
        == FDM::IBMConstraintSupport::Surface;
    return {
        surface ? "Lambda_s" : "lambda_b",
        surface ? "surface multiplier" : "body multiplier",
        surface ? ImmersedVariableLocation::SurfaceConstraintDof
                : ImmersedVariableLocation::BodyConstraintDof,
        3,
        ImmersedOwnershipKind::ConstraintOwner};
}

void addSolidUnknowns(ImmersedAlgorithmDescriptor& result) {
    result.unknowns.push_back({
        "U_s", "rigid-body translational velocity",
        ImmersedVariableLocation::SolidGlobalDof, 3,
        ImmersedOwnershipKind::SolidOwner});
    result.unknowns.push_back({
        "omega_s", "rigid-body angular velocity",
        ImmersedVariableLocation::SolidGlobalDof, 3,
        ImmersedOwnershipKind::SolidOwner});
}

std::string strategy(const FDM::IBMForcingConfig& config) {
    if (config.algorithm
        == FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian) {
        return "augmentedLagrangianKKT";
    }
    const FDM::IBMEnforcement enforcement=config.enforcement;
    switch (enforcement) {
        case FDM::IBMEnforcement::ExplicitIBM:
            return "explicitLaggedMultiplier";
        case FDM::IBMEnforcement::FractionalDLM:
            return "fractionalVariationalProjection";
        case FDM::IBMEnforcement::VelocityForcing:
            return "surfaceMassProjection";
        case FDM::IBMEnforcement::BrinkmanPenalty:
            return "dissipativePenalty";
        case FDM::IBMEnforcement::MonolithicKKT:
            return "monolithicKKT";
        case FDM::IBMEnforcement::GhostCell:
            break;
    }
    throw std::runtime_error(
        "Variational IBM descriptor received ghostCell enforcement.");
}

} // namespace

ImmersedAlgorithmDescriptor ghostCell() {
    ImmersedAlgorithmDescriptor result;
    result.id = "ghostCellIBM";
    result.referenceName = "sharp-interface GhostCellIBM";
    result.support = FDM::IBMConstraintSupport::Surface;
    result.representation = FDM::IBMRepresentation::SharpJump;
    result.enforcement = FDM::IBMEnforcement::GhostCell;
    result.solid = FDM::IBMSolidModel::Prescribed;
    result.solveBlocks.push_back({
        "S_IBM_BOUNDARY", "immersed boundary reconstruction",
        "boundaryStencilClosure", {}, {}, {}});
    validate(result);
    return result;
}

ImmersedAlgorithmDescriptor variational(
        const FDM::IBMForcingConfig& config) {
    FDM::validateIBMForcingSelection(config);
    if (config.constraintSupport
        == FDM::IBMConstraintSupport::SurfaceAndBody) {
        throw std::runtime_error(
            "One variational IBM descriptor cannot merge surface and body "
            "multipliers. Declare two explicit constraint blocks.");
    }

    ImmersedAlgorithmDescriptor result;
    result.id = FDM::toString(config.algorithm);
    result.referenceName = "Bhalla unified IBM formulation";
    result.support = config.constraintSupport;
    result.representation = config.representation;
    result.enforcement = config.enforcement;
    result.solid = config.solidModel;
    result.introducesMultiplier =
        config.enforcement != FDM::IBMEnforcement::BrinkmanPenalty;
    result.monolithic =
        config.enforcement == FDM::IBMEnforcement::MonolithicKKT;

    result.variational.fluidKineticIncrement = true;
    // 当前 predictor 已包含黏性项；约束校正块本身使用质量矩阵，不能把
    // predictor 的黏性离散虚报为 KKT Hessian 的组成部分。
    result.variational.viscousDissipation = false;
    result.variational.externalWork = true;
    result.variational.incompressibilityConstraint = result.monolithic;
    result.variational.immersedNoSlipConstraint =
        config.enforcement != FDM::IBMEnforcement::BrinkmanPenalty;
    result.variational.solidKineticIncrement =
        config.solidModel == FDM::IBMSolidModel::CoupledRigid
        || config.solidModel == FDM::IBMSolidModel::SelfPropelledRigid;
    result.variational.solidExternalWork =
        result.variational.solidKineticIncrement;
    result.variational.stationaryFunctional = result.monolithic
        ? "fluid/solid kinetic increment - solid external work "
          "+ p^T D u + Lambda^T(Ju-Gq-uDef)"
        : "minimum mass-norm velocity increment subject to Ju=Us";
    if (config.algorithm
        == FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian) {
        result.variational.stationaryFunctional +=
            " + gamma/2 ||Ju-Gq-uDef||^2_ML";
    }

    if (result.introducesMultiplier) {
        result.unknowns.push_back(multiplier(config));
        result.equations.push_back({
            "E_IBM_STATIONARITY", "immersed momentum stationarity",
            "M_f (u-u*)/dt + G p - J^T Lambda = 0", {}});
        result.constraints.push_back({
            "C_IBM_NO_SLIP", "immersed no-slip",
            config.constraintSupport == FDM::IBMConstraintSupport::Surface
                ? "J u - G q = uDef" : "C u = U_s",
            result.unknowns.front().id});
    }
    if (result.variational.solidKineticIncrement) {
        addSolidUnknowns(result);
        result.equations.push_back({
            "E_RIGID_TRANSLATION", "rigid-body translation",
            "M_s (U_s-U_s^n)/dt = F_ext-C_T Lambda", {"U_s"}});
        result.equations.push_back({
            "E_RIGID_ROTATION", "rigid-body rotation",
            "I_s (omega_s-omega_s^n)/dt + omega_s x I_s omega_s "
            "= T_ext-C_R Lambda", {"omega_s"}});
    }

    ImmersedSolveBlockDescriptor block;
    block.id = result.monolithic ? "S_FLUID_IBM_KKT" : "S_IBM_CONSTRAINT";
    block.name = result.monolithic
        ? "monolithic fluid-pressure-immersed constraint"
        : "immersed variational correction";
    block.strategy = strategy(config);
    for (const auto& equation : result.equations) {
        block.equations.push_back(equation.id);
    }
    for (const auto& unknown : result.unknowns) {
        block.unknowns.push_back(unknown.id);
    }
    for (const auto& constraint : result.constraints) {
        block.constraints.push_back(constraint.id);
    }
    result.solveBlocks.push_back(std::move(block));
    validate(result);
    return result;
}

void validate(const ImmersedAlgorithmDescriptor& descriptor) {
    if (descriptor.id.empty() || descriptor.referenceName.empty()) {
        throw std::runtime_error(
            "IBM algorithm descriptor requires id and referenceName.");
    }
    std::set<std::string> unknownIds;
    for (const auto& unknown : descriptor.unknowns) {
        if (unknown.id.empty() || unknown.name.empty()
            || unknown.components <= 0
            || !unknownIds.insert(unknown.id).second) {
            throw std::runtime_error(
                "IBM algorithm descriptor contains an invalid/duplicate unknown.");
        }
    }
    std::set<std::string> constraintIds;
    std::set<std::string> equationIds;
    for (const auto& equation : descriptor.equations) {
        if (equation.id.empty() || equation.name.empty()
            || equation.form.empty()
            || !equationIds.insert(equation.id).second) {
            throw std::runtime_error(
                "IBM algorithm descriptor contains an invalid/duplicate equation.");
        }
        for (const auto& unknown : equation.solvedUnknowns) {
            if (unknownIds.find(unknown) == unknownIds.end()) {
                throw std::runtime_error(
                    "IBM equation references undeclared unknown '"
                    + unknown + "'.");
            }
        }
    }
    for (const auto& constraint : descriptor.constraints) {
        if (constraint.id.empty() || constraint.name.empty()
            || constraint.equation.empty()
            || !constraintIds.insert(constraint.id).second
            || unknownIds.find(constraint.multiplierUnknown)
                == unknownIds.end()) {
            throw std::runtime_error(
                "IBM constraint descriptor has an invalid multiplier reference.");
        }
    }
    std::set<std::string> blockIds;
    for (const auto& block : descriptor.solveBlocks) {
        if (block.id.empty() || block.name.empty() || block.strategy.empty()
            || !blockIds.insert(block.id).second) {
            throw std::runtime_error(
                "IBM descriptor contains an invalid/duplicate solve block.");
        }
        for (const auto& unknown : block.unknowns) {
            if (unknownIds.find(unknown) == unknownIds.end()) {
                throw std::runtime_error(
                    "IBM solve block references undeclared unknown '"
                    + unknown + "'.");
            }
        }
        for (const auto& equation : block.equations) {
            if (equationIds.find(equation) == equationIds.end()) {
                throw std::runtime_error(
                    "IBM solve block references undeclared equation '"
                    + equation + "'.");
            }
        }
        for (const auto& constraint : block.constraints) {
            if (constraintIds.find(constraint) == constraintIds.end()) {
                throw std::runtime_error(
                    "IBM solve block references undeclared constraint '"
                    + constraint + "'.");
            }
        }
    }
    if (descriptor.monolithic
        && (!descriptor.variational.incompressibilityConstraint
            || !descriptor.variational.immersedNoSlipConstraint)) {
        throw std::runtime_error(
            "Monolithic IBM descriptor must declare pressure and no-slip constraints.");
    }
}

} // namespace SF::IBM::Descriptor
