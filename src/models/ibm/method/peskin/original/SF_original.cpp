/// @file SF_original.cpp
/// @brief Bhalla Algorithm 1：Peskin 原始法的显式滞后表面乘子与规则化传播。

#include "method/SF_method.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "method/SF_fieldAdapter.h"
#include "immersed/SF_immersed.h"
#include "operations/SF_fieldOps.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::IBM::Forcing {

FDM::ImmersedConstraintResult ImmersedForcingSystem::applyPeskinOriginal(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    if (!geometry_ || fields.size()!=1 || !fields.front()) {
        throw std::runtime_error(
            "Peskin original IBM requires one local structured Field.");
    }
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt<=0.0) {
        throw std::runtime_error(
            "Peskin original IBM requires finite time and positive dt.");
    }
    Field& field=*fields.front();
    lastField_=&field;
    const int total=field.TotalSize();
    if (!laggedMultiplier_.empty()
        && laggedMultiplier_.size()!=(size_t)total) {
        throw std::runtime_error(
            "Peskin IBM lagged multiplier topology changed.");
    }
    if (laggedMultiplier_.empty()) {
        laggedMultiplier_.assign((size_t)total,Vector3());
    }
    mask_.assign((size_t)total,0);
    multiplier_.assign((size_t)total,Vector3());
    lastResult_={};
    const int momentum=FieldAdapter::momentumIndex(field);
    const int energy=FieldAdapter::energyIndex(field);
    const Vector3 center=bodyModel_->center(targetTime);

    // 当前流体步只消费上一时间层由表面乘子传播得到的 Eulerian 力密度。
    for (int cell=0; cell<total; ++cell) {
        const Vector3 forceDensity=laggedMultiplier_[(size_t)cell];
        if (norm(forceDensity)==0.0) continue;
        int i=0,j=0,k=0;
        field.getIJK(cell,i,j,k);
        FieldAdapter::validateState(field,momentum,energy,i,j,k);
        const double rho=FieldAdapter::density(field,i,j,k);
        const Vector3 oldMomentum(
            field(i,j,k,momentum),field(i,j,k,momentum+1),
            field(i,j,k,momentum+2));
        const Vector3 oldVelocity=oldMomentum*(1.0/rho);
        const Vector3 delta=forceDensity*dt;
        const Vector3 newMomentum=oldMomentum+delta;
        const Vector3 newVelocity=newMomentum*(1.0/rho);
        FieldAdapter::setMomentum(field,momentum,i,j,k,newMomentum);
        const double newEnergy=field(i,j,k,energy)
            +dot(delta,(oldVelocity+newVelocity)*0.5);
        if (!std::isfinite(newEnergy)) {
            throw std::runtime_error(
                "Peskin IBM mechanical-work update is non-finite.");
        }
        field(i,j,k,energy)=newEnergy;
        const double volume=FieldAdapter::volume(field,i,j,k);
        const Vector3 forceOnBody=forceDensity*(-volume);
        const Vector3 position(
            field.X(i,j,k),field.Y(i,j,k),field.Z(i,j,k));
        lastResult_.forceOnBody=lastResult_.forceOnBody+forceOnBody;
        lastResult_.torqueOnBody=lastResult_.torqueOnBody
            +cross(position-center,forceOnBody);
        lastResult_.fluidMechanicalPower +=
            dot(forceDensity,(oldVelocity+newVelocity)*0.5)*volume;
    }

    buildSurfaceSystem(field,targetTime,dt);
    // 与 FTS/BP 共用同一局部 J 行遍历：各 rank 只读取拥有有效 halo 的
    // Eulerian edge，随后通过 Runtime 显式归并 marker 部分和。这里不能让
    // patch 端点改变 M^-1 的体积权重。
    std::vector<double> inverseMass((size_t)total,0.0);
    std::vector<Vector3> interpolated(
        surfaceSystem_.points.size(),Vector3());
    std::vector<double> diagonals(surfaceSystem_.points.size(),0.0);
    auto velocity = [&](int cell) {
        const FieldOps::Index3 index=FieldOps::index3(field,cell);
        if (!FieldOps::interior(field,index)
            || field.isSolverBoundaryPoint(index.i,index.j,index.k)) {
            throw std::runtime_error(
                "Peskin interpolation references a non-Eulerian cell.");
        }
        FieldOps::validateState(field,momentum,energy,index);
        return FieldOps::velocity(field,momentum,index);
    };
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        for (const auto& edge:point.interpolation) {
            if (edge.cell<0 || edge.cell>=total
                || !std::isfinite(edge.value)) {
                throw std::runtime_error(
                    "Peskin interpolation row contains an invalid edge.");
            }
            const FieldOps::Index3 index=FieldOps::index3(field,edge.cell);
            velocity(edge.cell);
            inverseMass[(size_t)edge.cell]=1.0/(
                FieldOps::density(field,index)*FieldOps::volume(field,index));
        }
        interpolated[marker]=FDM::Immersed::interpolate(
            point.interpolation,velocity);
        for (const auto& edge:point.interpolation) {
            diagonals[marker]+=edge.value*edge.value*point.measure
                *inverseMass[(size_t)edge.cell];
        }
        if (!std::isfinite(diagonals[marker]) || diagonals[marker]<=0.0) {
            throw std::runtime_error(
                "Peskin interpolation produced a singular mass diagonal.");
        }
    }
    reduceSurfaceVectors(interpolated);
    if (runtime_) runtime_->globalSum(diagonals);
    std::vector<Vector3> surfaceMultiplier(
        surfaceSystem_.points.size(),Vector3());
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        const double diagonal=diagonals[marker];
        if (!std::isfinite(diagonal) || diagonal<=0.0) {
            throw std::runtime_error(
                "Peskin IBM found a singular surface mass diagonal.");
        }
        // Explicit Peskin lambda is also a canonical constraint quantity.
        if (ownsSurfacePoint(point)) {
            surfaceMultiplier[marker]=(point.prescribedVelocity
                -interpolated[marker])*(1.0/(dt*diagonal));
        }
    }
    copySurfaceMultipliers(surfaceMultiplier);
    std::vector<Vector3> nextMultiplier((size_t)total,Vector3());
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        const Vector3 target=point.prescribedVelocity;
        const Vector3 lambda=surfaceMultiplier[marker];
        FDM::Immersed::spread(
            point.interpolation,point.measure,lambda,
            [&](int cell,const Vector3& force) {
                int i=0,j=0,k=0;
                field.getIJK(cell,i,j,k);
                const Vector3 densityForce=
                    force*(1.0/FieldAdapter::volume(field,i,j,k));
                nextMultiplier[(size_t)cell]=
                    nextMultiplier[(size_t)cell]+densityForce;
                multiplier_[(size_t)cell]=
                    multiplier_[(size_t)cell]+densityForce;
                mask_[(size_t)cell]=1;
            });
        if (ownsSurfacePoint(point)) {
            lastResult_.maximumVelocityResidual=std::max(
                lastResult_.maximumVelocityResidual,
                norm(interpolated[marker]-target));
            ++lastResult_.constrainedCells;
        }
    }
    laggedMultiplier_=std::move(nextMultiplier);
    finalizeDistributedResult(lastResult_);
    field.invalidateThermodynamicCache();
    lastResult_.performed=true;
    std::ostringstream detail;
    detail << "peskinOriginal surface markers="
           << lastResult_.constrainedCells
           << ", lagged multiplier advanced, pre-next-step max|Ju-Us|="
           << lastResult_.maximumVelocityResidual;
    lastResult_.detail=detail.str();
    return lastResult_;
}

} // namespace SF::IBM::Forcing
