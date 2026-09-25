/// @file SF_methodState.cpp
/// @brief 保存最近一次 IBM 约束 mask、乘子和单体 KKT 选择状态。

#include "method/SF_method.h"
#include "core/interfaces/SF_executionRuntime.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Forcing {

void ImmersedForcingSystem::setExecutionRuntime(
        FDM::IExecutionRuntime* runtime) {
    runtime_ = runtime;
}

void ImmersedForcingSystem::buildSurfaceSystem(
        Field& field, double targetTime, double dt) {
    const Vector3 center=bodyModel_->center(targetTime);
    surfaceSystem_=surfaceOperator_->build(
        field,bodyModel_->triangles(*geometry_,targetTime),center,
        config_.forcing.surfaceSupportRadius);
    surfaceNormalizations_=surfaceOperator_->localNormalizations();
    if (runtime_) runtime_->globalSum(surfaceNormalizations_);
    surfaceSystem_=surfaceOperator_->normalizeDistributed(surfaceNormalizations_);
    bodyModel_->prepareEquationView(surfaceSystem_,targetTime,dt);
}

void ImmersedForcingSystem::reduceSurfaceVectors(
        std::vector<Vector3>& values) const {
    if (values.empty() || !runtime_) return;
    std::vector<double> packed(values.size()*3u,0.0);
    for (std::size_t n=0;n<values.size();++n) {
        packed[3u*n]=values[n].x;
        packed[3u*n+1u]=values[n].y;
        packed[3u*n+2u]=values[n].z;
    }
    runtime_->globalSum(packed);
    for (std::size_t n=0;n<values.size();++n) {
        const Vector3 reduced{packed[3u*n],packed[3u*n+1u],packed[3u*n+2u]};
        if (!std::isfinite(reduced.x) || !std::isfinite(reduced.y)
            || !std::isfinite(reduced.z)) {
            throw std::runtime_error(
                "IBM distributed surface reduction produced a non-finite value.");
        }
        values[n]=reduced;
    }
}

void ImmersedForcingSystem::copySurfaceMultipliers(
        std::vector<Vector3>& values) const {
    if (values.empty() || !runtime_) return;
    if (values.size()!=surfaceSystem_.points.size()) {
        throw std::runtime_error(
            "IBM surface multiplier COPY size differs from marker layout.");
    }
    std::vector<std::int64_t> ids;
    ids.reserve(surfaceSystem_.points.size());
    for (const auto& point : surfaceSystem_.points) {
        if (!point.globalConstraintId.valid()) {
            throw std::runtime_error(
                "IBM surface multiplier COPY received invalid constraint id.");
        }
        ids.push_back(point.globalConstraintId.value());
    }
    for (int component=0; component<3; ++component) {
        std::vector<double> canonical(values.size(),0.0);
        for (std::size_t marker=0; marker<values.size(); ++marker) {
            const auto& point=surfaceSystem_.points[marker];
            const double value=component==0 ? values[marker].x
                : (component==1 ? values[marker].y : values[marker].z);
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "IBM surface multiplier COPY received non-finite lambda.");
            }
            if (ownsSurfacePoint(point)) canonical[marker]=value;
        }
        runtime_->copyCanonicalEntities(ids,canonical);
        for (std::size_t marker=0; marker<values.size(); ++marker) {
            if (component==0) values[marker].x=canonical[marker];
            else if (component==1) values[marker].y=canonical[marker];
            else values[marker].z=canonical[marker];
        }
    }
}

bool ImmersedForcingSystem::ownsSurfacePoint(
        const FDM::ImmersedSurfacePoint& point) const {
    if (!point.globalConstraintId.valid()) {
        throw std::runtime_error(
            "IBM surface point has an invalid canonical constraint identity.");
    }
    return !runtime_ || runtime_->ownsCanonicalEntity(
        point.globalConstraintId.value());
}

void ImmersedForcingSystem::finalizeDistributedResult(
        FDM::ImmersedConstraintResult& result) const {
    if (!runtime_) return;
    std::vector<double> sums{
        static_cast<double>(result.constrainedCells),
        result.kineticIncrementFunctional,
        result.forceOnBody.x,result.forceOnBody.y,result.forceOnBody.z,
        result.torqueOnBody.x,result.torqueOnBody.y,result.torqueOnBody.z,
        result.fluidMechanicalPower};
    runtime_->globalSum(sums);
    result.constrainedCells=static_cast<std::size_t>(std::llround(sums[0]));
    result.kineticIncrementFunctional=sums[1];
    result.forceOnBody={sums[2],sums[3],sums[4]};
    result.torqueOnBody={sums[5],sums[6],sums[7]};
    result.fluidMechanicalPower=sums[8];
    result.maximumVelocityResidual=runtime_->globalMaximum(
        result.maximumVelocityResidual);
    result.maximumStationarityResidual=runtime_->globalMaximum(
        result.maximumStationarityResidual);
}

bool ImmersedForcingSystem::usesMonolithicKKT() const {
    return config_.forcing.algorithm
            == FDM::IBMForcingAlgorithm::DFMImplicitPrescribed
        || config_.forcing.algorithm
            == FDM::IBMForcingAlgorithm::DFMImplicitSelfPropelled
        || config_.forcing.algorithm
            == FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian;
}

double ImmersedForcingSystem::constraintMask(
        const Field& field, int i, int j, int k) const {
    if (&field != lastField_ || mask_.empty()) return 0.0;
    return mask_[static_cast<size_t>(field.getIdx(i, j, k))] ? 1.0 : 0.0;
}

Vector3 ImmersedForcingSystem::multiplier(
        const Field& field, int i, int j, int k) const {
    if (&field != lastField_ || multiplier_.empty()) return {};
    return multiplier_[static_cast<size_t>(field.getIdx(i, j, k))];
}

} // namespace SF::IBM::Forcing
