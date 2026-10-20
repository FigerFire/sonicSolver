/// @file SF_recovery.cpp
/// @brief 表面 DLM/KKT 解的状态恢复、伴随传播和作用反作用积分。

#include "method/fictitiousDomain/implicit/SF_recovery.h"

#include "SF_config.h"
#include "immersed/SF_immersed.h"
#include "operations/SF_diagnostics.h"
#include "operations/SF_fieldOps.h"
#include "operations/SF_loads.h"
#include "operations/SF_validation.h"
#include "solid/SF_bodyModel.h"

#include <algorithm>

namespace SF::IBM::Monolithic {

RecoveryResult recoverSolution(
        Field& field,
        const FDM::ImmersedSurfaceSystem& system,
        const FDM::ImmersedKKTState& state,
        const FDM::IBMForcingConfig& config,
        const Forcing::IBodyModel& body,
        double targetTime) {
    Validation::requireMultiplierCount(
        state.surfaceMultiplier.size(),system.points.size());
    Validation::requireFinite(targetTime,"IBM KKT targetTime");

    RecoveryResult result;
    result.mask.assign(static_cast<std::size_t>(field.TotalSize()),0);
    result.multiplier.assign(
        static_cast<std::size_t>(field.TotalSize()),Vector3());
    result.solid = Kinematics::resolve(body,state);
    result.constraint.maximumStationarityResidual =
        state.maximumStationarityResidual;
    result.constraint.kineticIncrementFunctional =
        state.kineticIncrementFunctional;
    const int momentum = FieldOps::momentumIndex(field);

    for (std::size_t marker=0; marker<system.points.size(); ++marker) {
        const auto& point = system.points[marker];
        const Vector3 lambda = state.surfaceMultiplier[marker];
        Validation::requirePositive(point.measure,"IBM surface measure");
        Validation::requireFinite(lambda,"IBM surface multiplier");

        const Vector3 interpolated = FDM::Immersed::interpolate(
            point.interpolation,[&](int cell) {
                Validation::requireIndex(
                    cell,field.TotalSize(),"IBM interpolation cell");
                return FieldOps::velocity(
                    field,momentum,FieldOps::index3(field,cell));
            });
        FDM::Immersed::spread(
            point.interpolation,point.measure,lambda,
            [&](int cell,const Vector3& force) {
                Validation::requireIndex(
                    cell,field.TotalSize(),"IBM spreading cell");
                const FieldOps::Index3 index = FieldOps::index3(field,cell);
                const std::size_t storage = static_cast<std::size_t>(cell);
                result.multiplier[storage] = result.multiplier[storage]
                    +force*(1.0/FieldOps::volume(field,index));
                result.mask[storage] = 1;
            });

        const Vector3 target = Kinematics::targetVelocity(
            result.solid,point,body,targetTime);
        result.constraint.maximumVelocityResidual = std::max(
            result.constraint.maximumVelocityResidual,
            norm(interpolated-target));
        const Vector3 bodyForce = Loads::forceOnBody(lambda,point.measure);
        result.constraint.forceOnBody =
            result.constraint.forceOnBody+bodyForce;
        result.constraint.torqueOnBody = result.constraint.torqueOnBody
            +Loads::torqueOnBody(point.relativePosition,bodyForce);
        result.constraint.fluidMechanicalPower +=
            Loads::fluidMechanicalPower(lambda,target,point.measure);
    }

    Validation::requireConstraintResidual(
        result.constraint.maximumVelocityResidual,
        system.constraintTolerance,"fullyImplicitDLM");
    result.constraint.constrainedCells = system.points.size();
    result.constraint.performed = true;
    result.constraint.detail = Diagnostics::describe(
        FDM::toString(config.algorithm),"surface",
        result.constraint,&result.solid);
    return result;
}

} // namespace SF::IBM::Monolithic
