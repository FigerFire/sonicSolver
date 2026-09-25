/// @file SF_FTS.cpp
/// @brief Bhalla Algorithm 5 的 FTS 速度强迫：预测后用表面乘子投影速度。

#include "method/SF_method.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "constraint/variational/SF_projection.h"
#include "immersed/SF_immersed.h"
#include "operations/SF_diagnostics.h"
#include "operations/SF_fieldOps.h"
#include "operations/SF_loads.h"
#include "operations/SF_validation.h"

#include "SF_fluidStateModel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::IBM::Forcing {
namespace {

struct CellIncrement {
    Vector3 momentum;
    double rho = 0.0;
    bool touched = false;
};

} // namespace

FDM::ImmersedConstraintResult ImmersedForcingSystem::applyVelocityForcingFTS(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) {
    if (!geometry_ || fields.size() != 1 || !fields.front()) {
        throw std::runtime_error(
            "velocityForcing currently requires one local structured Field.");
    }
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "velocityForcing requires finite targetTime and positive dt.");
    }
    Field& field = *fields.front();
    lastField_ = &field;
    mask_.assign(static_cast<size_t>(field.TotalSize()), 0);
    multiplier_.assign(static_cast<size_t>(field.TotalSize()), Vector3());
    lastResult_ = {};
    const int momentum = FieldOps::momentumIndex(field);
    const int energy = FieldOps::energyIndex(field);
    if (momentum < 0 || momentum+2 >= field.NVar()
        || energy < 0 || energy >= field.NVar()) {
        throw std::runtime_error(
            "velocityForcing requires canonical momentum and energy fields.");
    }
    buildSurfaceSystem(field,targetTime,dt);

    std::vector<double> inverseMass(
        static_cast<size_t>(field.TotalSize()),0.0);
    std::vector<Vector3> interpolated(
        surfaceSystem_.points.size(),Vector3());
    auto velocity = [&](int cell) {
        const FieldOps::Index3 index = FieldOps::index3(field,cell);
        if (!FieldOps::interior(field,index)
            || field.isSolverBoundaryPoint(index.i,index.j,index.k)) {
            throw std::runtime_error(
                "velocityForcing interpolation references a boundary cell.");
        }
        FieldOps::validateState(field,momentum,energy,index);
        return FieldOps::velocity(field,momentum,index);
    };
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        if (!std::isfinite(point.measure) || point.measure <= 0.0) {
            throw std::runtime_error(
                "velocityForcing received a non-positive surface measure.");
        }
        Validation::requireFinite(
            point.prescribedVelocity,"velocityForcing target velocity");
        for (const auto& weight : point.interpolation) {
            if (weight.cell < 0
                || weight.cell >= field.TotalSize()
                || !std::isfinite(weight.value)) {
                throw std::runtime_error(
                    "velocityForcing received an invalid interpolation row.");
            }
            velocity(weight.cell);
            const FieldOps::Index3 index=FieldOps::index3(
                field,weight.cell);
            inverseMass[static_cast<size_t>(weight.cell)]=
                1.0/(FieldOps::density(field,index)
                     *FieldOps::volume(field,index));
        }
        interpolated[marker] = FDM::Immersed::interpolate(
            point.interpolation, velocity);
    }
    // `diag(S)` 是分区无关的 surface Schur 图不变量；将它与 Ju* 一起
    // 输出，能把 predictor/halo 偏差和 J--M^-1--J^T 组装偏差明确分开。
    std::vector<double> schurDiagonal(surfaceSystem_.points.size(),0.0);
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        double localDiagonal=0.0;
        for (const auto& edge:point.interpolation) {
            localDiagonal += edge.value*edge.value
                *inverseMass[static_cast<std::size_t>(edge.cell)];
        }
        schurDiagonal[marker] = dt*point.measure*localDiagonal;
    }
    if (runtime_) runtime_->globalSum(schurDiagonal);
    double minimumSchurDiagonal=std::numeric_limits<double>::infinity();
    double maximumSchurDiagonal=0.0;
    for (double value:schurDiagonal) {
        if (!std::isfinite(value) || value<=0.0) {
            throw std::runtime_error(
                "velocityForcingFTS constructed a non-positive Schur diagonal.");
        }
        minimumSchurDiagonal=std::min(minimumSchurDiagonal,value);
        maximumSchurDiagonal=std::max(maximumSchurDiagonal,value);
    }
    if (surfaceNormalizations_.size()!=surfaceSystem_.points.size()) {
        throw std::runtime_error(
            "velocityForcingFTS surface normalization diagnostics are stale.");
    }
    double minimumNormalization=std::numeric_limits<double>::infinity();
    double maximumNormalization=0.0;
    for (double value:surfaceNormalizations_) {
        if (!std::isfinite(value) || value<=0.0) {
            throw std::runtime_error(
                "velocityForcingFTS constructed a non-positive J normalization.");
        }
        minimumNormalization=std::min(minimumNormalization,value);
        maximumNormalization=std::max(maximumNormalization,value);
    }
    // 每个 rank 仅计算它拥有的 Eulerian edge 对 Ju* 的部分贡献；Runtime
    // 显式 SUM 后所有 rank 得到同一 marker 约束误差，而 IBM 数学核不识别 MPI。
    reduceSurfaceVectors(interpolated);
    // 记录投影前的 Ju* 误差，作为 serial/MPI 分区一致性诊断。它与
    // 投影后的约束残差分开：若二者前者随 partition 改变，问题在流体
    // predictor/halo；若仅后者改变，才应检查 J M^-1 J^T 或 lambda COPY。
    double maximumPredictedVelocityResidual=0.0;
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        maximumPredictedVelocityResidual=std::max(
            maximumPredictedVelocityResidual,
            norm(interpolated[marker]
                -surfaceSystem_.points[marker].prescribedVelocity));
    }
    std::vector<Vector3> constraintError;
    constraintError.reserve(surfaceSystem_.points.size());
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        constraintError.push_back(
            surfaceSystem_.points[marker].prescribedVelocity-interpolated[marker]);
    }
    Variational::SurfaceSchurProjectionProblem projectionProblem;
    projectionProblem.surface=&surfaceSystem_;
    projectionProblem.inverseEulerianMass=&inverseMass;
    projectionProblem.constraintError=std::move(constraintError);
    projectionProblem.dt=dt;
    projectionProblem.maxIterations=
        config_.forcing.constraintSolverMaxIterations;
    projectionProblem.relativeTolerance=
        config_.forcing.constraintSolverRelativeTolerance;
    projectionProblem.ownedMarkers.reserve(surfaceSystem_.points.size());
    for (const auto& point : surfaceSystem_.points) {
        projectionProblem.ownedMarkers.push_back(
            ownsSurfacePoint(point) ? 1u : 0u);
    }
    projectionProblem.globalSum=[this](std::vector<double>& values) {
        if (runtime_) runtime_->globalSum(values);
    };
    Variational::SurfaceSchurProjectionSolution projection=
        Variational::solveSurfaceSchurProjection(projectionProblem);
    // Schur vectors are replicated after global reductions. Canonical marker
    // owner publishes lambda; Runtime performs sparse COPY to each J^T
    // consumer. This is deliberately neither a SUM nor an average.
    copySurfaceMultipliers(projection.multiplier);

    std::vector<CellIncrement> increments(
        static_cast<size_t>(field.TotalSize()));
    FDM::ImmersedConstraintResult result;
    for (size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        const Vector3 markerVelocity=point.prescribedVelocity
            -projectionProblem.constraintError[marker];
        const Vector3 target=point.prescribedVelocity;
        const Vector3 lambda=projection.multiplier[marker];
        if (ownsSurfacePoint(point)) {
            const Vector3 bodyForce = Loads::forceOnBody(lambda,point.measure);
            result.forceOnBody = result.forceOnBody+bodyForce;
            result.torqueOnBody = result.torqueOnBody
                +Loads::torqueOnBody(point.relativePosition,bodyForce);
            result.fluidMechanicalPower +=
                Loads::fluidMechanicalPower(
                    lambda,(markerVelocity+target)*0.5,point.measure);
            ++result.constrainedCells;
        }
        FDM::Immersed::spread(
            point.interpolation, point.measure, lambda,
            [&](int localCell, const Vector3& force) {
            const size_t cell = static_cast<size_t>(localCell);
            const FieldOps::Index3 index = FieldOps::index3(field,localCell);
            const Vector3 forceDensity =
                force*(1.0/FieldOps::volume(field,index));
            increments[cell].momentum = increments[cell].momentum
                + forceDensity*dt;
            increments[cell].rho = FieldOps::density(field,index);
            increments[cell].touched = true;
            multiplier_[cell] = multiplier_[cell] + forceDensity;
            mask_[cell] = 1;
        });
    }
    if (surfaceSystem_.points.empty()) {
        throw std::runtime_error(
            "velocityForcing found no supported surface markers.");
    }

    for (size_t cell=0; cell<increments.size(); ++cell) {
        if (!increments[cell].touched) continue;
        const FieldOps::Index3 index = FieldOps::index3(
            field,static_cast<int>(cell));
        const Vector3 oldMomentum = FieldOps::momentum(
            field,momentum,index);
        const Vector3 newMomentum = oldMomentum+increments[cell].momentum;
        const Vector3 oldVelocity = oldMomentum*(1.0/increments[cell].rho);
        const Vector3 newVelocity = newMomentum*(1.0/increments[cell].rho);
        FieldOps::setMomentum(field,momentum,index,newMomentum);
        const double newEnergy = field(index.i,index.j,index.k,energy)
            +Variational::mechanicalEnergyIncrement(
                increments[cell].momentum,oldVelocity,newVelocity);
        if (!std::isfinite(newEnergy)) {
            throw std::runtime_error(
                "velocityForcing mechanical-work energy update is non-finite.");
        }
        field(index.i,index.j,index.k,energy) = newEnergy;
        const Vector3 stationarity = (newVelocity-oldVelocity)
            *(increments[cell].rho/dt)
            -increments[cell].momentum*(1.0/dt);
        result.maximumStationarityResidual = std::max(
            result.maximumStationarityResidual,norm(stationarity));
        result.kineticIncrementFunctional += 0.5*increments[cell].rho
            *dot(newVelocity-oldVelocity,newVelocity-oldVelocity)
            *FieldOps::volume(field,index)/dt;
    }
    // 用实际更新后的 Eulerian 状态验证完整 J M^-1 J^T 系统。
    std::vector<Vector3> corrected(surfaceSystem_.points.size(),Vector3());
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        corrected[marker] = FDM::Immersed::interpolate(
            point.interpolation, [&](int cell) {
            return FieldOps::velocity(
                field,momentum,FieldOps::index3(field,cell));
        });
    }
    reduceSurfaceVectors(corrected);
    for (std::size_t marker=0;marker<surfaceSystem_.points.size();++marker) {
        const auto& point=surfaceSystem_.points[marker];
        const Vector3 target = point.prescribedVelocity;
        result.maximumVelocityResidual = std::max(
            result.maximumVelocityResidual,norm(corrected[marker]-target));
    }
    finalizeDistributedResult(result);
    Validation::requireConstraintResidual(
        result.maximumVelocityResidual,
        config_.forcing.constraintTolerance,"velocityForcingFTS");
    field.invalidateThermodynamicCache();
    result.performed = true;
    result.detail = Diagnostics::describe(
        "velocityForcingFTS","surface",result)
        +", preProjection max|Ju*-Us|="
        +std::to_string(maximumPredictedVelocityResidual)
        +", diag(S)=["+std::to_string(minimumSchurDiagonal)
        +","+std::to_string(maximumSchurDiagonal)+"]"
        +", sum(Jraw)=["+std::to_string(minimumNormalization)
        +","+std::to_string(maximumNormalization)+"]"
        +", SchurCG(iterations="+std::to_string(projection.iterations)
        +", relativeResidual="
        +std::to_string(projection.maximumRelativeResidual)+")";
    lastResult_ = result;
    return result;
}

} // namespace SF::IBM::Forcing
