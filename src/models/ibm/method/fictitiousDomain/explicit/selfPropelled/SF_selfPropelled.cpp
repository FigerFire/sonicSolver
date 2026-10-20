/// @file SF_selfPropelled.cpp
/// @brief Bhalla Algorithm 2 的显式自推进 DFM：滞后乘子推进并更新下一步约束力。

#include "method/SF_method.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "method/SF_fieldAdapter.h"
#include "method/fictitiousDomain/SF_selfPropelled.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::IBM::Forcing {

FDM::ImmersedConstraintResult
ImmersedForcingSystem::applyDFMExplicitSelfPropelled(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    if (!geometry_ || fields.size()!=1 || !fields.front()) {
        throw std::runtime_error(
            "Explicit self-propelled DFM requires one structured Field.");
    }
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt<=0.0) {
        throw std::runtime_error(
            "Explicit self-propelled DFM requires finite time and positive dt.");
    }
    Field& field=*fields.front();
    lastField_=&field;
    const int total=field.TotalSize();
    if (!laggedMultiplier_.empty()
        && laggedMultiplier_.size()!=(size_t)total) {
        throw std::runtime_error(
            "Explicit self-propelled DFM lagged multiplier topology changed.");
    }
    if (laggedMultiplier_.empty()) {
        laggedMultiplier_.assign((size_t)total,Vector3());
    }
    mask_.assign((size_t)total,0);
    multiplier_.assign((size_t)total,Vector3());
    lastResult_={};
    const int momentum=FieldAdapter::momentumIndex(field);
    const int energy=FieldAdapter::energyIndex(field);
    const auto& body=bodyOperator_->build(
        field,*geometry_,*bodyModel_,targetTime);

    // Algorithm 2 的时间显式性：本步流体只消费上一时间层乘子。
    for (const auto& point : body.points) {
        const size_t cell=(size_t)point.localCell;
        int i=0,j=0,k=0;
        field.getIJK(point.localCell,i,j,k);
        FieldAdapter::validateState(field,momentum,energy,i,j,k);
        const double rho=FieldAdapter::density(field,i,j,k);
        const Vector3 oldMomentum(
            field(i,j,k,momentum),field(i,j,k,momentum+1),
            field(i,j,k,momentum+2));
        const Vector3 oldVelocity=oldMomentum*(1.0/rho);
        const Vector3 delta=laggedMultiplier_[cell]*dt;
        const Vector3 newMomentum=oldMomentum+delta;
        const Vector3 newVelocity=newMomentum*(1.0/rho);
        FieldAdapter::setMomentum(field,momentum,i,j,k,newMomentum);
        const double newEnergy=field(i,j,k,energy)
            +dot(delta,(oldVelocity+newVelocity)*0.5);
        if (!std::isfinite(newEnergy)) {
            throw std::runtime_error(
                "Explicit self-propelled DFM produced non-finite energy.");
        }
        field(i,j,k,energy)=newEnergy;
        const Vector3 appliedForce=laggedMultiplier_[cell]
            *(-point.dualVolume);
        lastResult_.forceOnBody=lastResult_.forceOnBody+appliedForce;
        lastResult_.torqueOnBody=lastResult_.torqueOnBody
            +cross(point.relativePosition,appliedForce);
    }

    std::vector<SelfPropelled::Sample> samples;
    samples.reserve(body.points.size());
    for (const auto& point : body.points) {
        int i=0,j=0,k=0;
        field.getIJK(point.localCell,i,j,k);
        samples.push_back({
            FieldAdapter::velocity(field,momentum,i,j,k),
            bodyModel_->deformationVelocityAt(point.position,targetTime),
            point.relativePosition,
            FieldAdapter::density(field,i,j,k)*point.dualVolume});
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
    std::vector<Vector3> nextMultiplier((size_t)total,Vector3());
    for (size_t n=0; n<body.points.size(); ++n) {
        const auto& point=body.points[n];
        int i=0,j=0,k=0;
        field.getIJK(point.localCell,i,j,k);
        const double rho=FieldAdapter::density(field,i,j,k);
        const Vector3 current=FieldAdapter::velocity(field,momentum,i,j,k);
        const Vector3 target=SelfPropelled::targetVelocity(
            config_.forcing.rigidMotionMode,generalized,
            point.relativePosition,samples[n].deformationVelocity);
        const Vector3 lambda=(target-current)*(rho/dt);
        const size_t cell=(size_t)point.localCell;
        nextMultiplier[cell]=lambda;
        multiplier_[cell]=lambda;
        mask_[cell]=1;
        lastResult_.maximumVelocityResidual=std::max(
            lastResult_.maximumVelocityResidual,norm(current-target));
        ++lastResult_.constrainedCells;
    }
    laggedMultiplier_=std::move(nextMultiplier);
    finalizeDistributedResult(lastResult_);
    Vector3 linear=bodyModel_->linearVelocity();
    Vector3 angular=bodyModel_->angularVelocity();
    if (config_.forcing.rigidMotionMode
        == FDM::IBMRigidMotionMode::Rotate) angular=generalized;
    else linear=generalized;
    bodyModel_->advanceCoupled(linear,angular,dt);
    field.invalidateThermodynamicCache();
    lastResult_.performed=true;
    std::ostringstream detail;
    detail << "dfmExplicitSelfPropelled body cells="
           << lastResult_.constrainedCells
           << ", lagged multiplier advanced, q=("
           << generalized.x << "," << generalized.y << ","
           << generalized.z << "), pre-next-step max|u-Us|="
           << lastResult_.maximumVelocityResidual;
    lastResult_.detail=detail.str();
    return lastResult_;
}

} // namespace SF::IBM::Forcing
