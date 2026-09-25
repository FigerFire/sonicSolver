/// @file SF_prescribed.cpp
/// @brief 兼容的给定速度 DFM-FTS 体投影；新 case 应显式选择 algorithm。

#include "method/SF_method.h"
#include "constraint/variational/SF_projection.h"
#include "operations/SF_diagnostics.h"
#include "operations/SF_fieldOps.h"
#include "operations/SF_loads.h"
#include "operations/SF_validation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::IBM::Forcing {

FDM::ImmersedConstraintResult
ImmersedForcingSystem::applyDFMFractionalStepPrescribed(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    if (!geometry_) {
        throw std::runtime_error(
            "ImmersedForcingSystem must be configured before constraint projection.");
    }
    if (usesMonolithicKKT()) {
        throw std::runtime_error(
            "surface Lambda_s requires the pressureBased monolithic KKT "
            "stage; densityBased sequential projection is not permitted.");
    }
    if (fields.size() != 1 || !fields.front()) {
        throw std::runtime_error(
            "fractionalDLM currently requires one local structured Field; "
            "multi-patch/MPI multiplier ownership is not implemented.");
    }
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "fractionalDLM requires finite targetTime and positive dt.");
    }

    Field& field = *fields.front();
    lastField_ = &field;
    mask_.assign(static_cast<size_t>(field.TotalSize()), 0);
    multiplier_.assign(static_cast<size_t>(field.TotalSize()), Vector3());
    lastResult_ = {};

    const int momentum = FieldOps::momentumIndex(field);
    const int energy = FieldOps::energyIndex(field);
    if (momentum < 0 || momentum + 2 >= field.NVar()
        || energy < 0 || energy >= field.NVar()) {
        throw std::runtime_error(
            "fractionalDLM requires canonical density, momentum and energy.");
    }

    const auto& body = bodyOperator_->build(
        field,*geometry_,*bodyModel_,targetTime);
    for (const auto& point : body.points) {
        if (point.localCell < 0 || point.localCell >= field.TotalSize()
            || !std::isfinite(point.dualVolume)
            || point.dualVolume <= 0.0
            || !std::isfinite(point.position.x)
            || !std::isfinite(point.position.y)
            || !std::isfinite(point.position.z)
            || !std::isfinite(point.relativePosition.x)
            || !std::isfinite(point.relativePosition.y)
            || !std::isfinite(point.relativePosition.z)) {
            throw std::runtime_error(
                "fractionalDLM received an invalid body constraint point.");
        }
        const FieldOps::Index3 index = FieldOps::index3(
            field,point.localCell);
        FieldOps::validateState(field,momentum,energy,index);
        const double rho = FieldOps::density(field,index);
        const Vector3 oldMomentum = FieldOps::momentum(
            field,momentum,index);
        const double oldEnergy = field(index.i,index.j,index.k,energy);
        const Vector3 oldVelocity = oldMomentum*(1.0/rho);
        const Vector3 targetVelocity =
            bodyModel_->velocityAt(point.position,targetTime);
        const Variational::PointProjectionSolution projection =
            Variational::solvePointProjection({
                oldVelocity,targetVelocity,rho,point.dualVolume,dt});
        const Vector3 newMomentum = projection.correctedVelocity*rho;
        const Vector3 momentumIncrement = newMomentum-oldMomentum;
        const double newEnergy = oldEnergy
            +Variational::mechanicalEnergyIncrement(
                momentumIncrement,oldVelocity,projection.correctedVelocity);
        if (!std::isfinite(newEnergy)) {
            throw std::runtime_error(
                "fractionalDLM mechanical-work energy update is non-finite.");
        }

        FieldOps::setMomentum(field,momentum,index,newMomentum);
        field(index.i,index.j,index.k,energy) = newEnergy;
        const size_t cell = static_cast<size_t>(point.localCell);
        mask_[cell] = 1;
        multiplier_[cell] = projection.multiplier;
        const Vector3 bodyForce = Loads::forceOnBody(
            projection.multiplier,point.dualVolume);
        lastResult_.forceOnBody = lastResult_.forceOnBody+bodyForce;
        lastResult_.torqueOnBody = lastResult_.torqueOnBody
            +Loads::torqueOnBody(point.relativePosition,bodyForce);
        const Vector3 midpointVelocity =
            (oldVelocity+projection.correctedVelocity)*0.5;
        lastResult_.fluidMechanicalPower +=
            Loads::fluidMechanicalPower(
                projection.multiplier,midpointVelocity,point.dualVolume);
        lastResult_.maximumVelocityResidual = std::max(
            lastResult_.maximumVelocityResidual,
            projection.constraintResidual);
        lastResult_.maximumStationarityResidual = std::max(
            lastResult_.maximumStationarityResidual,
            projection.stationarityResidual);
        lastResult_.kineticIncrementFunctional +=
            projection.kineticIncrementFunctional;
        ++lastResult_.constrainedCells;
    }
    finalizeDistributedResult(lastResult_);
    Validation::requireConstraintResidual(
        lastResult_.maximumVelocityResidual,
        config_.forcing.constraintTolerance,"fractionalDLM");
    field.invalidateThermodynamicCache();
    lastResult_.performed = true;
    lastResult_.detail = Diagnostics::describe(
        "dfmFractionalStepPrescribed","body",lastResult_);
    return lastResult_;
}

} // namespace SF::IBM::Forcing
