/// @file SF_systemDescription.cpp
/// @brief Explain/IR 枚举的稳定字符串（WHAT/ORDER description formatting）。
///
/// 它们只服务 explain/diagnostic，不参与任何决策。

#include "SF_resolvedSimulationSystem.h"

namespace SF::System {

const ExecutableOperation* findExecutableOperation(
        const ExecutableEquationSystem& system, OperationStage stage) {
    const auto found = std::find_if(
        system.operations.begin(),system.operations.end(),
        [stage](const ExecutableOperation& operation) {
            return operation.stage == stage;
        });
    return found == system.operations.end() ? nullptr : &*found;
}

const char* toString(OperationStage stage) {
    switch (stage) {
        case OperationStage::Prepare: return "prepare";
        case OperationStage::MomentumAssemble: return "momentum-assemble";
        case OperationStage::MomentumSolve: return "momentum-solve";
        case OperationStage::PressureBoundaryPrepare:
            return "pressure-boundary-prepare";
        case OperationStage::PressureAssemble: return "pressure-assemble";
        case OperationStage::PressureSolve: return "pressure-solve";
        case OperationStage::PressureUpdatePrepare:
            return "pressure-update-prepare";
        case OperationStage::VelocityCorrect: return "velocity-correct";
        case OperationStage::FluxCorrect: return "flux-correct";
        case OperationStage::CorrectionCommit: return "correction-commit";
        case OperationStage::StepCommit: return "step-commit";
    }
    return "unknown-operation-stage";
}

const char* toString(OperationCapability capability) {
    switch (capability) {
        case OperationCapability::ConservativeExplicit: return "conservative.explicit";
        case OperationCapability::PressureSchedule: return "pressure.schedule";
        case OperationCapability::MomentumPredictor: return "momentum.predictor";
        case OperationCapability::PressureCorrection: return "pressure.correction";
        case OperationCapability::PressureLinearSolve: return "pressure.linearSystem";
        case OperationCapability::PressureBoundary: return "pressure.boundary";
        case OperationCapability::VelocityCorrection: return "velocity.correction";
        case OperationCapability::FluxCorrection: return "flux.correction";
        case OperationCapability::EulerianPhaseExecution: return "eulerian.phaseExecution";
        case OperationCapability::ImmersedConstraint: return "immersed.constraint";
    }
    return "unknown-operation-capability";
}

const char* toString(PhysicsTemplateKind kind) {
    switch (kind) {
        case PhysicsTemplateKind::SingleFluid: return "singleFluid";
        case PhysicsTemplateKind::HomogeneousMixture:
            return "homogeneousMixture";
        case PhysicsTemplateKind::OneFluidInterface: return "oneFluidInterface";
        case PhysicsTemplateKind::EulerianEulerian: return "eulerianEulerian";
    }
    return "unknown";
}

const char* toString(VariableLocation value) {
    switch (value) {
        case VariableLocation::EulerianCell: return "EulerianCell";
        case VariableLocation::EulerianFace: return "EulerianFace";
        case VariableLocation::BodyConstraint: return "BodyConstraint";
        case VariableLocation::SurfaceConstraint: return "SurfaceConstraint";
        case VariableLocation::SolidGlobal: return "SolidGlobal";
    }
    return "UnknownLocation";
}

const char* toString(OwnershipKind value) {
    switch (value) {
        case OwnershipKind::EulerianGlobalDof: return "EulerianGlobalDof";
        case OwnershipKind::CanonicalFace: return "CanonicalFace";
        case OwnershipKind::ConstraintGlobalDof: return "ConstraintGlobalDof";
        case OwnershipKind::SolidGlobalDof: return "SolidGlobalDof";
    }
    return "UnknownOwnership";
}

const char* toString(OriginKind value) {
    switch (value) {
        case OriginKind::BuiltinDefault: return "builtin-default";
        case OriginKind::BuiltinPreset: return "builtin-preset";
        case OriginKind::Model: return "model";
        case OriginKind::User: return "user";
        case OriginKind::Generated: return "generated";
    }
    return "unknown-origin";
}

const char* toString(ModificationKind value) {
    switch (value) {
        case ModificationKind::Add: return "add";
        case ModificationKind::Extend: return "extend";
        case ModificationKind::Replace: return "replace";
        case ModificationKind::Disable: return "disable";
    }
    return "unknown-modification";
}

const char* toString(EquationCategory value) {
    switch (value) {
        case EquationCategory::PhysicalEquation: return "physical";
        case EquationCategory::ConstraintEquation: return "constraint";
        case EquationCategory::AlgorithmicDerivedEquation:
            return "algorithmic-derived";
        case EquationCategory::AlgebraicRelation: return "algebraic-relation";
    }
    return "unknown-equation-category";
}

const char* toString(TransformationState value) {
    switch (value) {
        case TransformationState::NotRegistered: return "not-registered";
        case TransformationState::RegisteredAndActive:
            return "registered-active";
        case TransformationState::RegisteredButNotApplicable:
            return "registered-not-applicable";
        case TransformationState::RegisteredButInvalid:
            return "registered-invalid";
        case TransformationState::Applied: return "applied";
    }
    return "unknown-transformation-state";
}



const char* toString(ResourceAccessMode value) {
    switch (value) {
        case ResourceAccessMode::Read: return "read";
        case ResourceAccessMode::Write: return "write";
        case ResourceAccessMode::ReadWrite: return "read-write";
    }
    return "unknown-access";
}

const char* toString(SynchronizationRequirement value) {
    switch (value) {
        case SynchronizationRequirement::None: return "none";
        case SynchronizationRequirement::ReadHalo: return "read-halo";
        case SynchronizationRequirement::WriteOwned: return "write-owned";
    }
    return "unknown-synchronization";
}

} // namespace SF::System
