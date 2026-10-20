/// @file SF_couplingGraph.cpp
/// @brief 串行 IBM 稀疏耦合边、对偶体积和 partition-independent 身份构造。

#include "topology/SF_couplingGraph.h"

#include "immersed/SF_regularizedKernel.h"
#include "immersed/SF_linearReproducing.h"

#include "core/mesh/SF_nodalQuadrature.h"

#include <cmath>
#include <algorithm>
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

double dualVolume(const Field& field,int i,int j,int k) {
    return StructuredMesh::nodalVolume(field,i,j,k);
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
        double supportRadius, FDM::IBMSurfaceNormalization normalization) {
    if (!std::isfinite(supportRadius) || supportRadius<=0.0) {
        throw std::runtime_error(
            "IBM coupling graph requires a positive supportRadius.");
    }
    (void)FDM::toString(normalization);
    normalization_=normalization;
    const bool linear=normalization_==FDM::IBMSurfaceNormalization::LinearReproducing;
    rawMoments_.assign(linear?markers.size()*16u:0u,0.0);
    rawBasis_.clear();
    if (linear) rawBasis_.resize(markers.size());
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
                    if (linear) {
                        const Vector3 offset=(position-marker.position)*(1.0/supportRadius);
                        const std::array<double,4> basis={1.,offset.x,offset.y,offset.z};
                        rawBasis_[markerIndex].push_back(basis);
                        for (int a=0;a<4;++a) for (int b=0;b<4;++b)
                            rawMoments_[16u*markerIndex+4*a+b]+=raw*basis[a]*basis[b];
                    }
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
        const std::vector<double>& globalNormalizations,
        const std::vector<double>& globalMoments) {
    if (globalNormalizations.size()!=rawRows_.size()) {
        throw std::runtime_error(
            "IBM coupling graph normalization size differs from marker count.");
    }
    const bool linear=normalization_==FDM::IBMSurfaceNormalization::LinearReproducing;
    if (globalMoments.size()!=(linear?rawRows_.size()*16u:0u))
        throw std::runtime_error("IBM transfer moment matrix count differs from selected recipe.");
    rows_=rawRows_;
    for (std::size_t marker=0;marker<rows_.size();++marker) {
        const double normalization=globalNormalizations[marker];
        if (!std::isfinite(normalization)||normalization<=0.0) {
            throw std::runtime_error(
                "IBM coupling graph received a non-positive global marker normalization.");
        }
        if (linear) {
            std::array<double,16> matrix;
            std::copy_n(globalMoments.begin()+16u*marker,16,matrix.begin());
            if (std::abs(matrix[0]-normalization)>64.*std::numeric_limits<double>::epsilon()*normalization)
                throw std::runtime_error("IBM raw moment mass and normalization SUM disagree.");
            try {
                const auto coefficients=FDM::Immersed::linearMomentCoefficients(matrix);
                for (std::size_t edge=0;edge<rows_[marker].size();++edge) {
                    double value=0.;
                    for (int a=0;a<4;++a) value+=rawBasis_[marker][edge][a]*coefficients[a];
                    rows_[marker][edge].value*=value;
                    if (!std::isfinite(rows_[marker][edge].value))
                        throw std::runtime_error("non-finite corrected weight");
                }
            } catch (const std::exception& error) {
                throw std::runtime_error("IBM linearReproducing marker "+std::to_string(marker)+": "+error.what());
            }
            continue;
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
