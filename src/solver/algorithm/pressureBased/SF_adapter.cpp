/// @file SF_adapter.cpp
/// @brief 压力基算法与通用 Field/EquationSet 状态之间的适配。

#include "solver/algorithm/pressureBased/SF_adapter.h"
#include "solver/algorithm/immersed/SF_immersedStrategy.h"

#include <sstream>
#include <stdexcept>

namespace SF::PressureBased {

Algorithm::Algorithm(const FDM::SolverConfig& config)
    : corrector_(config.boundaries, config.pressure,
                 config.numerics.idealGasGamma),
      kkt_(config.pressure,config.numerics.idealGasGamma) {}

FDM::CapabilitySet Algorithm::capabilities() const {
    return {
        FDM::Capability::PressureCorrection,
        FDM::Capability::PressureJumpConsumer,
        FDM::Capability::LagrangeMultiplierField,
        FDM::Capability::SingleFieldExecution,
        FDM::Capability::MultiFieldExecution};
}

FDM::FlowAlgorithmResult Algorithm::correct(
        FDM::FlowAlgorithmContext& context,
        double dt) {
    if (context.fields.empty()) {
        throw std::runtime_error(
            "PressureBased::Algorithm requires at least one Field.");
    }
    if (context.prepareBoundaryState) context.prepareBoundaryState();
    if (context.equationSystem) {
        context.equationSystem->preparePressureCorrection(context.fields);
    }

    corrector_.setExecutionRuntime(context.executionRuntime);
    corrector_.setInterfaceJumpProvider(
        [&](const Field& geometry)
            -> const FDM::IInterfaceJumpCondition* {
            return context.equationSystem
                ? context.equationSystem->interfaceJumpCondition(geometry)
                : nullptr;
    });

    if (context.immersed.system && !context.immersed.constraint) {
        throw std::runtime_error(
            "pressureBased received an immersed system without its "
            "constraint numerical adapter.");
    }
    const bool hasImmersedConstraint=context.immersed.constraint!=nullptr;
    if (hasImmersedConstraint) {
        context.immersed.constraint->setExecutionRuntime(context.executionRuntime);
        if (!context.immersed.system) {
            throw std::runtime_error(
                "pressureBased received an immersed constraint adapter "
                "without the unified IImmersedSystem selection service.");
        }
        if (const auto* provider = context.immersed.constraint->systemProvider();
            provider && provider != context.immersed.system) {
            throw std::runtime_error(
                "pressureBased received constraint and selection services "
                "from different IBM managers.");
        }
        ImmersedAlgorithm::validateForAlgorithm(
            *context.immersed.system,
            FDM::SolverAlgorithm::PressureBased);
        const bool monolithic=context.immersed.system->methodSelection().enforcement
            == FDM::IBMEnforcement::MonolithicKKT;
        if (monolithic && context.executionRuntime
            && context.executionRuntime->distributed()) {
            throw std::runtime_error(
                "distributed monolithic IBM KKT is not enabled until the "
                "ConstraintGlobalDof -> owned HYPRE row -> lambda COPY "
                "path is connected to the pressure block.");
        }
        if (monolithic && (context.fields.size()!=1
                           || !context.fields.front())) {
            throw std::runtime_error(
                "surface monolithic KKT currently requires one local Field; "
                "multi-patch/MPI multiplier ownership is not implemented.");
        }
        if (monolithic) {
            const KKTCorrectionSummary result=kkt_.correct(
                *context.fields.front(),context.targetTime,dt,
                *context.immersed.constraint);
            std::ostringstream detail;
            detail << "KKT iterations=" << result.iterations
                   << ", residual=" << result.relativeResidual
                   << ", maxDiv(before/after)="
                   << result.maxDivergenceBefore << "/"
                   << result.maxDivergenceAfter
                   << ", HYPRE(rebuilds/solves)="
                   << result.structureRebuilds << "/" << result.linearSolves
                   << ", " << result.immersed.detail;
            return {true,true,detail.str()};
        }
    }

    const CorrectionSummary pressure = context.fields.size() == 1
        ? corrector_.correct(*context.fields.front(), dt)
        : corrector_.correct(context.fields, dt);
    for (Field* field : context.fields) {
        if (field) field->invalidateThermodynamicCache();
    }

    FDM::ImmersedConstraintResult immersed;
    if (hasImmersedConstraint) {
        // 顺序型 IBM 的公开算法顺序：流体 predictor -> pressure correction
        // -> variational constraint projection。最后一次投影严格满足 IBM 约束；
        // 它不是单体 pressure/IBM KKT，也不会以隐藏迭代伪装成单体算法。
        immersed=context.immersed.constraint->projectPredictedState(
            context.fields,context.targetTime,dt);
        for (Field* field:context.fields) {
            if (field) field->invalidateThermodynamicCache();
        }
    }

    std::ostringstream detail;
    detail << "iterations=" << pressure.iterations
           << ", residual=" << pressure.finalResidual
           << ", maxDiv(before)=" << pressure.maxDivergenceBefore
           << ", maxDiv(after)=" << pressure.maxDivergenceAfter
           << ", interfaceFaces=" << pressure.interfaceFaces
           << ", maxJump(target/correction)="
           << pressure.maxTargetPressureJump << "/"
           << pressure.maxCorrectionPressureJump
               << ", HYPRE(rebuilds/solves)="
               << pressure.structureRebuilds << "/"
               << pressure.linearSolves;
    if (immersed.performed) detail << ", sequentialIBM=" << immersed.detail;
    return {true, true, detail.str()};
}

} // namespace SF::PressureBased
