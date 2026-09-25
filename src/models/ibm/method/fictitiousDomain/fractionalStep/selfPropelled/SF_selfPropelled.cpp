/// @file SF_selfPropelled.cpp
/// @brief Bhalla Algorithm 3：一阶 FTS 自推进 DFM 的体约束投影。

#include "method/SF_method.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "method/fictitiousDomain/SF_selfPropelled.h"
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
ImmersedForcingSystem::applyDFMFractionalStepSelfPropelled(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    if (!geometry_ || fields.size() != 1 || !fields.front()) {
        throw std::runtime_error(
            "DFM FTS self-propulsion requires one local structured Field.");
    }
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "DFM FTS self-propulsion requires finite time and positive dt.");
    }
    Field& field = *fields.front();
    lastField_ = &field;
    mask_.assign((size_t)field.TotalSize(),0);
    multiplier_.assign((size_t)field.TotalSize(),Vector3());
    lastResult_ = {};
    const int momentum = FieldOps::momentumIndex(field);
    const int energy = FieldOps::energyIndex(field);
    const auto& body = bodyOperator_->build(
        field,*geometry_,*bodyModel_,targetTime);
    std::vector<SelfPropelled::Sample> samples;
    samples.reserve(body.points.size());
    for (const auto& point : body.points) {
        const FieldOps::Index3 index = FieldOps::index3(
            field,point.localCell);
        FieldOps::validateState(field,momentum,energy,index);
        samples.push_back({
            FieldOps::velocity(field,momentum,index),
            bodyModel_->deformationVelocityAt(point.position,targetTime),
            point.relativePosition,
            FieldOps::density(field,index)*point.dualVolume});
    }
    SelfPropelled::GeneralizedSystem generalizedSystem=
        SelfPropelled::assembleGeneralizedSystem(
            samples,config_.forcing.rigidMotionMode);
    std::vector<double> generalizedEntries(12,0.0);
    for (int n=0;n<9;++n) generalizedEntries[(size_t)n]=
        generalizedSystem.matrix[(size_t)n];
    generalizedEntries[9]=generalizedSystem.rhs.x;
    generalizedEntries[10]=generalizedSystem.rhs.y;
    generalizedEntries[11]=generalizedSystem.rhs.z;
    if (runtime_) runtime_->globalSum(generalizedEntries);
    for (int n=0;n<9;++n) generalizedSystem.matrix[(size_t)n]=
        generalizedEntries[(size_t)n];
    generalizedSystem.rhs={generalizedEntries[9],generalizedEntries[10],
                           generalizedEntries[11]};
    const Vector3 generalized=SelfPropelled::solveGeneralizedSystem(
        generalizedSystem,config_.forcing.rigidMotionMode,
        config_.forcing.externalForce,config_.forcing.externalTorque,dt);
    for (size_t n=0; n<body.points.size(); ++n) {
        const auto& point = body.points[n];
        const FieldOps::Index3 index = FieldOps::index3(
            field,point.localCell);
        const double rho = FieldOps::density(field,index);
        const Vector3 oldMomentum = FieldOps::momentum(
            field,momentum,index);
        const Vector3 oldVelocity = oldMomentum*(1.0/rho);
        const Vector3 target = SelfPropelled::targetVelocity(
            config_.forcing.rigidMotionMode,generalized,
            point.relativePosition,samples[n].deformationVelocity);
        const Variational::PointProjectionSolution projection =
            Variational::solvePointProjection({
                oldVelocity,target,rho,point.dualVolume,dt});
        const Vector3 newMomentum = projection.correctedVelocity*rho;
        const Vector3 momentumIncrement = newMomentum-oldMomentum;
        FieldOps::setMomentum(field,momentum,index,newMomentum);
        const double newEnergy = field(index.i,index.j,index.k,energy)
            +Variational::mechanicalEnergyIncrement(
                momentumIncrement,oldVelocity,projection.correctedVelocity);
        if (!std::isfinite(newEnergy)) {
            throw std::runtime_error(
                "DFM FTS self-propulsion produced non-finite energy.");
        }
        field(index.i,index.j,index.k,energy)=newEnergy;
        const size_t cell=(size_t)point.localCell;
        mask_[cell]=1;
        multiplier_[cell]=projection.multiplier;
        const Vector3 force=Loads::forceOnBody(
            projection.multiplier,point.dualVolume);
        lastResult_.forceOnBody=lastResult_.forceOnBody+force;
        lastResult_.torqueOnBody=lastResult_.torqueOnBody
            +Loads::torqueOnBody(point.relativePosition,force);
        lastResult_.fluidMechanicalPower +=
            Loads::fluidMechanicalPower(
                projection.multiplier,
                (oldVelocity+projection.correctedVelocity)*0.5,
                point.dualVolume);
        lastResult_.maximumVelocityResidual=std::max(
            lastResult_.maximumVelocityResidual,
            projection.constraintResidual);
        lastResult_.maximumStationarityResidual=std::max(
            lastResult_.maximumStationarityResidual,
            projection.stationarityResidual);
        lastResult_.kineticIncrementFunctional +=
            projection.kineticIncrementFunctional;
        ++lastResult_.constrainedCells;
    }
    finalizeDistributedResult(lastResult_);
    Validation::requireConstraintResidual(
        lastResult_.maximumVelocityResidual,
        config_.forcing.constraintTolerance,
        "dfmFractionalStepSelfPropelled");
    Vector3 linear=bodyModel_->linearVelocity();
    Vector3 angular=bodyModel_->angularVelocity();
    if (config_.forcing.rigidMotionMode
        == FDM::IBMRigidMotionMode::Rotate) angular=generalized;
    else linear=generalized;
    bodyModel_->advanceCoupled(linear,angular,dt);
    field.invalidateThermodynamicCache();
    lastResult_.performed=true;
    const Kinematics::SolidKinematics solid{linear,angular};
    lastResult_.detail=Diagnostics::describe(
        "dfmFractionalStepSelfPropelled","body",lastResult_,&solid);
    return lastResult_;
}

} // namespace SF::IBM::Forcing
