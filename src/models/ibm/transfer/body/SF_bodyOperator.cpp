/// @file SF_bodyOperator.cpp
/// @brief 构造 Eulerian body ownership、对偶体积与相对位置拓扑。

#include "transfer/body/SF_bodyOperator.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Forcing {
namespace {

double nodalWeight(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    double weight = 1.0;
    if (i == ng || i == ng + field.NX() - 1) weight *= 0.5;
    if (j == ng || j == ng + field.NY() - 1) weight *= 0.5;
    if (k == ng || k == ng + field.NZ() - 1) weight *= 0.5;
    return weight;
}

double dualVolume(const Field& field, int i, int j, int k) {
    const double inverseVolume = field.Jac(i,j,k);
    if (!std::isfinite(inverseVolume) || inverseVolume == 0.0) {
        throw std::runtime_error(
            "BodyOperator found an invalid cell Jacobian.");
    }
    return nodalWeight(field,i,j,k)/std::abs(inverseVolume);
}

bool finiteVector(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

} // namespace

const FDM::ImmersedBodySystem& BodyOperator::build(
        const Field& field,
        const GeoProcessing::STLGeometry& geometry,
        const IBodyModel& bodyModel,
        double targetTime) {
    if (!std::isfinite(targetTime)) {
        throw std::runtime_error(
            "BodyOperator requires a finite targetTime.");
    }
    system_.points.clear();
    system_.centerOfMass = bodyModel.center(targetTime);
    const int ng = field.NG();
    for (int k=ng; k<ng+field.NZ(); ++k) {
        for (int j=ng; j<ng+field.NY(); ++j) {
            for (int i=ng; i<ng+field.NX(); ++i) {
                const Vector3 position(
                    field.X(i,j,k),field.Y(i,j,k),field.Z(i,j,k));
                if (!finiteVector(position)) {
                    throw std::runtime_error(
                        "BodyOperator found a non-finite Eulerian point.");
                }
                if (!geometry.contains(
                        bodyModel.referencePoint(position,targetTime))) {
                    continue;
                }
                // 已登记的共享 GlobalDof 仅由唯一 Eulerian owner 生成体约束点。
                // 这使 DFM/Brinkman 的动量更新和刚体载荷无需对副本做平均。
                if (field.globalDofId(i,j,k) >= 0
                    && !field.isGlobalDofOwner(i,j,k)) {
                    continue;
                }
                if (field.isSolverBoundaryPoint(i,j,k)) {
                    throw std::runtime_error(
                        "BodyOperator found a body constraint overlapping "
                        "a fitted physical boundary.");
                }
                FDM::ImmersedBodyPoint point;
                point.localCell = field.getIdx(i,j,k);
                point.position = position;
                point.relativePosition = position-system_.centerOfMass;
                point.dualVolume = dualVolume(field,i,j,k);
                system_.points.push_back(point);
            }
        }
    }
    bool hasDistributedDof=false;
    for (int k=ng;k<ng+field.NZ()&&!hasDistributedDof;++k)
        for (int j=ng;j<ng+field.NY()&&!hasDistributedDof;++j)
            for (int i=ng;i<ng+field.NX();++i)
                if (field.globalDofId(i,j,k)>=0) {
                    hasDistributedDof=true;
                    break;
                }
    // 分布式时一个 rank 不与刚体相交是正常状态；全局 body 是否为空由后续
    // Runtime 归并的总约束数/广义质量矩阵统一判定。串行保留明确几何诊断。
    if (system_.points.empty() && !hasDistributedDof) {
        throw std::runtime_error(
            "BodyOperator found no Eulerian points inside the moving body; "
            "refine the grid or check STL position and SI units.");
    }
    return system_;
}

} // namespace SF::IBM::Forcing
