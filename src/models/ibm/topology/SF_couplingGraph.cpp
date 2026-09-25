/// @file SF_couplingGraph.cpp
/// @brief 串行 IBM 稀疏耦合边、对偶体积和 partition-independent 身份构造。

#include "topology/SF_couplingGraph.h"

#include "immersed/SF_regularizedKernel.h"

#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace SF::IBM::Topology {
namespace {

bool interpolationUnknown(const Field& field, int i, int j, int k) {
    const int ng=field.NG();
    return i>=ng && i<ng+field.NX()
        && j>=ng && j<ng+field.NY()
        && k>=ng && k<ng+field.NZ()
        && field.CellFlag(i,j,k)==FLUID_CELL
        && !field.isSolverBoundaryPoint(i,j,k);
}

double dualVolume(const Field& field, int i, int j, int k) {
    const double inverseVolume=field.Jac(i,j,k);
    if (!std::isfinite(inverseVolume) || inverseVolume==0.0) {
        throw std::runtime_error(
            "IBM coupling graph found an invalid Eulerian Jacobian.");
    }
    // IBM surface support 只允许使用 `interpolationUnknown()` 已筛出的
    // 流体自由度；真实物理边界点不会进入本图。因此这里的对偶体积必须
    // 直接采用 source-canonical metric 的完整 `1/|J|`。不能按当前 patch
    // 的 i/j/k 端点再乘 1/2：MPI 内部分割面同样是本地数组端点，那会把
    // 分区拓扑误当作物理边界，使 J/S 的权重、力和功随 split 改变。
    return 1.0/std::abs(inverseVolume);
}

std::int64_t canonicalDof(const Field& field, int i, int j, int k) {
    const int registered=field.globalDofId(i,j,k);
    return registered>=0 ? static_cast<std::int64_t>(registered)
                         : static_cast<std::int64_t>(field.getIdx(i,j,k));
}

} // namespace

void CouplingGraph::build(
        const Field& field,
        const std::vector<SurfaceMarker>& markers,
        double supportRadius) {
    if (!std::isfinite(supportRadius) || supportRadius<=0.0) {
        throw std::runtime_error(
            "IBM coupling graph requires a positive supportRadius.");
    }
    rawRows_.clear();
    rawRows_.resize(markers.size());
    rows_.clear();
    rawNormalizations_.assign(markers.size(),0.0);
    const int ng=field.NG();
    for (std::size_t markerIndex=0;
         markerIndex<markers.size(); ++markerIndex) {
        const auto& marker=markers[markerIndex];
        auto& row=rawRows_[markerIndex];
        std::unordered_set<std::int64_t> identities;
        double normalization=0.0;
        for (int k=ng;k<ng+field.NZ();++k) {
            for (int j=ng;j<ng+field.NY();++j) {
                for (int i=ng;i<ng+field.NX();++i) {
                    if (!interpolationUnknown(field,i,j,k)) continue;
                    // 未登记 GlobalDof 的单块网格保持串行语义；已登记的共享点只由
                    // 唯一 Eulerian owner 生成 edge，禁止 replica 参与归一化或传播。
                    if (field.globalDofId(i,j,k)>=0
                        && !field.isGlobalDofOwner(i,j,k)) continue;
                    const Vector3 position(
                        field.X(i,j,k),field.Y(i,j,k),field.Z(i,j,k));
                    const double kernel=FDM::Immersed::wendlandC2(
                        norm(position-marker.position),supportRadius);
                    if (kernel<=0.0) continue;
                    const double volume=dualVolume(field,i,j,k);
                    const double raw=FDM::Immersed::volumeWeightedKernel(
                        kernel,volume);
                    const std::int64_t identity=canonicalDof(field,i,j,k);
                    if (!identities.insert(identity).second) {
                        throw std::runtime_error(
                            "IBM coupling graph found a duplicate "
                            "(marker, EulerianDof) edge.");
                    }
                    FDM::ImmersedInterpolationWeight edge;
                    edge.cell=field.getIdx(i,j,k);
                    edge.value=raw;
                    edge.globalEulerianDofId=identity;
                    edge.dualVolume=volume;
                    row.push_back(edge);
                    normalization+=raw;
                }
            }
        }
        if (!std::isfinite(normalization) || normalization < 0.0) {
            throw std::runtime_error(
                "IBM coupling graph generated an invalid local normalization for marker "
                +std::to_string(marker.id.value())
                +".");
        }
        rawNormalizations_[markerIndex]=normalization;
    }
}

void CouplingGraph::normalize(
        const std::vector<double>& globalNormalizations) {
    if (globalNormalizations.size()!=rawRows_.size()) {
        throw std::runtime_error(
            "IBM coupling graph normalization size differs from marker count.");
    }
    rows_=rawRows_;
    for (std::size_t marker=0;marker<rows_.size();++marker) {
        const double normalization=globalNormalizations[marker];
        if (!std::isfinite(normalization)||normalization<=0.0) {
            throw std::runtime_error(
                "IBM coupling graph received a non-positive global marker normalization.");
        }
        double localSum=0.0;
        for(auto& edge:rows_[marker]) {
            edge.value/=normalization;
            localSum+=edge.value;
        }
        if (!std::isfinite(localSum)||localSum<0.0||localSum>1.0+1.0e-12) {
            throw std::runtime_error(
                "IBM coupling graph generated an invalid normalized local J row.");
        }
    }
}

} // namespace SF::IBM::Topology
