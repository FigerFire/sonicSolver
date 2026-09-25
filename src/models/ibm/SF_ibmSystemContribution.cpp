/// @file SF_ibmSystemContribution.cpp
/// @brief IBM mathematical metadata contribution implementation.

#include "SF_ibmSystemContribution.h"

#include "SF_immersedSystem.h"
#include "core/system/SF_systemContribution.h"

#include <stdexcept>
#include <utility>

namespace SF::IBM::SystemContribution {
using namespace SF::System;
namespace {

VariableLocation location(FDM::ImmersedVariableLocation value) {
    switch (value) {
        case FDM::ImmersedVariableLocation::EulerianGlobalDof:
            return VariableLocation::EulerianCell;
        case FDM::ImmersedVariableLocation::BodyConstraintDof:
            return VariableLocation::BodyConstraint;
        case FDM::ImmersedVariableLocation::SurfaceConstraintDof:
            return VariableLocation::SurfaceConstraint;
        case FDM::ImmersedVariableLocation::SolidGlobalDof:
            return VariableLocation::SolidGlobal;
    }
    throw std::runtime_error("Unknown immersed variable location.");
}

OwnershipKind ownership(FDM::ImmersedOwnershipKind value) {
    switch (value) {
        case FDM::ImmersedOwnershipKind::EulerianOwner:
            return OwnershipKind::EulerianGlobalDof;
        case FDM::ImmersedOwnershipKind::ConstraintOwner:
            return OwnershipKind::ConstraintGlobalDof;
        case FDM::ImmersedOwnershipKind::SolidOwner:
            return OwnershipKind::SolidGlobalDof;
    }
    throw std::runtime_error("Unknown immersed ownership kind.");
}


} // namespace

void contribute(
        SF::System::SystemContribution& system,
        const FDM::ImmersedAlgorithmDescriptor& immersed) {
    system.recordContribution(
        "model.ibm."+immersed.id,"immersed-boundary contribution");
    for (const auto& source : immersed.unknowns) {
        UnknownDescriptor unknown;
        unknown.id = source.id;
        unknown.name = source.name;
        unknown.location = location(source.location);
        unknown.components = source.components;
        unknown.ownership = ownership(source.ownership);
        unknown.shape = source.components == 1
            ? ValueShape::Scalar : ValueShape::Vector;
        unknown.role = unknown.location == VariableLocation::BodyConstraint
                || unknown.location == VariableLocation::SurfaceConstraint
            ? UnknownRole::Multiplier : UnknownRole::Algebraic;
        unknown.storageBinding = StorageBinding::SpecializedExecutor;
        unknown.nameSpace = "immersed";
        system.addUnknown(std::move(unknown));
    }
    for (const auto& source : immersed.constraints) {
        system.addConstraint({
            source.id,source.name,source.equation,source.multiplierUnknown});
    }
    for (const auto& source : immersed.equations) {
        const std::string unknown = source.solvedUnknowns.empty()
            ? source.id : source.solvedUnknowns.front();
        EquationDescriptor descriptor{
            source.id,source.name,source.form,source.solvedUnknowns};
        descriptor.category = EquationCategory::ConstraintEquation;
        system.addEquation(std::move(descriptor),
            Equation::named(source.id,
                Equation::constraint({source.form.empty()
                    ? source.id : source.form})
                    == Equation::Symbol{unknown}));
    }
    bool hasConstraintTransformation = false;
    for (const auto& source : immersed.solveBlocks) {
        ExecutionPolicy policy;
        policy.id = source.id;
        policy.name = source.name;
        policy.strategyName = source.strategy;
        if (source.strategy == "boundaryStencilClosure") {
            policy.kind = ExecutionPolicyKind::BoundaryClosure;
            policy.strategyKind = FDM::SolveStrategyKind::BoundaryClosure;
        } else if (source.strategy == "monolithicKKT"
                   || source.strategy == "augmentedLagrangianKKT") {
            policy.kind = ExecutionPolicyKind::MonolithicKKT;
            policy.strategyKind = FDM::SolveStrategyKind::MonolithicKKT;
            policy.leafKind = PlanNodeKind::BlockSolve;
            policy.leafOperation = "ibm.kkt.solve";
            hasConstraintTransformation = true;
        } else if (source.strategy == "explicitLaggedMultiplier"
                   || source.strategy == "fractionalVariationalProjection"
                   || source.strategy == "surfaceMassProjection"
                   || source.strategy == "dissipativePenalty") {
            policy.kind = ExecutionPolicyKind::ConstraintProjection;
            policy.strategyKind = FDM::SolveStrategyKind::ConstraintSolve;
            policy.leafKind = PlanNodeKind::Correct;
            policy.leafOperation = "ibm.constraint.project";
            // Brinkman/BP is a direct momentum penalty contribution and has
            // no multiplier constraint for the transformer to rewrite.
            hasConstraintTransformation = hasConstraintTransformation
                || source.strategy != "dissipativePenalty";
        } else {
            throw std::runtime_error(
                "Unknown immersed solve strategy '"+source.strategy+"'.");
        }
        policy.equations = source.equations;
        policy.unknowns = source.unknowns;
        policy.constraints = source.constraints;
        policy.priority = 200;
        if (immersed.monolithic) {
            policy.unknowns.insert(policy.unknowns.begin(),{"pPrime","U"});
            policy.equations.insert(
                policy.equations.begin(),{"E_MOMENTUM","E_PRESSURE"});
        }
        system.addExecutionPolicy(std::move(policy));
    }
    if (hasConstraintTransformation) {
        system.requestTransformation({
            "immersedConstraint","immersed constraint transformation",
            200,true,{}});
    }
}


} // namespace SF::IBM::SystemContribution
