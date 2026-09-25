/// @file SF_BP.cpp
/// @brief Bhalla Algorithm 5 的 BP 速度强迫：动量方程内施加 Brinkman 罚力。

#include "method/SF_method.h"
#include "method/SF_fieldAdapter.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::IBM::Forcing {

FDM::ImmersedConstraintResult ImmersedForcingSystem::applyVelocityForcingBP(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    if (!geometry_ || fields.size() != 1 || !fields.front()) {
        throw std::runtime_error(
            "brinkmanPenalty currently requires one local structured Field.");
    }
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "brinkmanPenalty requires finite targetTime and positive dt.");
    }
    const double coefficient = config_.forcing.penaltyCoefficient;
    if (!std::isfinite(coefficient) || coefficient <= 0.0) {
        throw std::runtime_error(
            "brinkmanPenalty requires a finite positive penaltyCoefficient.");
    }
    Field& field = *fields.front();
    lastField_ = &field;
    mask_.assign(static_cast<size_t>(field.TotalSize()), 0);
    multiplier_.assign(static_cast<size_t>(field.TotalSize()), Vector3());
    lastResult_ = {};
    const int momentum = FieldAdapter::momentumIndex(field);
    const int energy = FieldAdapter::energyIndex(field);
    if (momentum < 0 || momentum+2 >= field.NVar()
        || energy < 0 || energy >= field.NVar()) {
        throw std::runtime_error(
            "brinkmanPenalty requires canonical momentum and energy fields.");
    }
    const auto& body = bodyOperator_->build(
        field,*geometry_,*bodyModel_,targetTime);
    FDM::ImmersedConstraintResult result;
    for (const auto& point : body.points) {
        if (!std::isfinite(point.dualVolume) || point.dualVolume <= 0.0) {
            throw std::runtime_error(
                "brinkmanPenalty received a non-positive dual volume.");
        }
        int i=0,j=0,k=0;
        field.getIJK(point.localCell,i,j,k);
        FieldAdapter::validateState(field,momentum,energy,i,j,k);
        const double rho = FieldAdapter::density(field,i,j,k);
        const Vector3 oldMomentum(
            field(i,j,k,momentum),field(i,j,k,momentum+1),
            field(i,j,k,momentum+2));
        const Vector3 oldVelocity = oldMomentum*(1.0/rho);
        const Vector3 target =
            bodyModel_->velocityAt(point.position,targetTime);
        const Vector3 forceDensity =
            (target-oldVelocity)*coefficient;
        const Vector3 delta = forceDensity*dt;
        const Vector3 newMomentum = oldMomentum+delta;
        const Vector3 newVelocity = newMomentum*(1.0/rho);
        FieldAdapter::setMomentum(field,momentum,i,j,k,newMomentum);
        const double newEnergy = field(i,j,k,energy)
            + dot(delta,(oldVelocity+newVelocity)*0.5);
        if (!std::isfinite(newEnergy)) {
            throw std::runtime_error(
                "brinkmanPenalty mechanical-work energy update is non-finite.");
        }
        field(i,j,k,energy) = newEnergy;
        const size_t cell = static_cast<size_t>(point.localCell);
        mask_[cell] = 1;
        multiplier_[cell] = forceDensity;
        const Vector3 bodyForce = forceDensity*(-point.dualVolume);
        result.forceOnBody = result.forceOnBody+bodyForce;
        result.torqueOnBody = result.torqueOnBody
            +cross(point.relativePosition,bodyForce);
        result.fluidMechanicalPower +=
            dot(forceDensity,(oldVelocity+newVelocity)*0.5)
            *point.dualVolume;
        result.maximumVelocityResidual = std::max(
            result.maximumVelocityResidual,norm(newVelocity-target));
        ++result.constrainedCells;
    }
    finalizeDistributedResult(result);
    field.invalidateThermodynamicCache();
    result.performed = true;
    std::ostringstream detail;
    detail << "velocityForcingBP body cells=" << result.constrainedCells
           << ", coefficient=" << coefficient << ", bodyForce=("
           << result.forceOnBody.x << "," << result.forceOnBody.y << ","
           << result.forceOnBody.z << ")";
    result.detail = detail.str();
    lastResult_ = result;
    return result;
}

} // namespace SF::IBM::Forcing
