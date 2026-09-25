/// @file SF_solution.cpp
/// @brief 全隐式 DLM/KKT 解的校验、恢复和事务式提交。

#include "method/SF_method.h"

#include "method/fictitiousDomain/implicit/SF_recovery.h"
#include "operations/SF_validation.h"

namespace SF::IBM::Forcing {

FDM::ImmersedConstraintResult
ImmersedForcingSystem::acceptMonolithicSolution(
        Field& field,
        double targetTime,
        double dt,
        const FDM::ImmersedKKTState& state) {
    Validation::requirePositive(dt,"fullyImplicitDLM dt");
    Monolithic::RecoveryResult recovery = Monolithic::recoverSolution(
        field,surfaceSystem_,state,config_.forcing,*bodyModel_,targetTime);

    if (surfaceSystem_.solidEquation.solveGeneralizedVelocity) {
        bodyModel_->advanceCoupled(
            recovery.solid.linearVelocity,
            recovery.solid.angularVelocity,
            dt);
    }

    lastField_ = &field;
    mask_ = std::move(recovery.mask);
    multiplier_ = std::move(recovery.multiplier);
    lastResult_ = std::move(recovery.constraint);
    return lastResult_;
}

} // namespace SF::IBM::Forcing
