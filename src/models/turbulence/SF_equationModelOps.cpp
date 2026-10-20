/// @file SF_equationModelOps.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

#include "SF_equationModelOps.h"

#include "core/mesh/SF_meshBoundaryGeometry.h"
#include "core/mesh/SF_dimension.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace SF::Turbulence::ModelOps {

bool isPhysical(const Field& geometry, int i, int j, int k) {
    return StructuredMesh::BoundaryGeometry::isPhysicalPoint(geometry, i, j, k);
}

double strainSquared(
        const Field& geometry,
        const Physics::PhaseSystems::PhaseState& phase,
        int i, int j, int k) {
    const SymmTensor3 strain = Math::symm(CENTRAL2::grad(
        geometry, phase.primitive.velocity, i, j, k));
    return 2.0 * Math::doubleDot(strain, strain);
}

double vorticityMagnitude(
        const Field& geometry,
        const Physics::PhaseSystems::PhaseState& phase,
        int i, int j, int k) {
    return norm(Math::curl(CENTRAL2::grad(
        geometry, phase.primitive.velocity, i, j, k)));
}

double filterWidth(
        const Field& geometry,
        int i, int j, int k,
        double scale) {
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::runtime_error(
            "Smagorinsky filterScale must be finite and positive.");
    }
    double product = 1.0;
    int dimensions = 0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        int di=0,dj=0,dk=0;
        if (axis==0) di=1;
        else if (axis==1) dj=1;
        else dk=1;
        const int lo = geometry.NG();
        const int coordinate = axis==0 ? i : (axis==1 ? j : k);
        const int count = axis==0 ? geometry.NX()
            : (axis==1 ? geometry.NY() : geometry.NZ());
        if (coordinate == lo+count-1) {
            di=-di; dj=-dj; dk=-dk;
        }
        const double dx=geometry.X(i+di,j+dj,k+dk)-geometry.X(i,j,k);
        const double dy=geometry.Y(i+di,j+dj,k+dk)-geometry.Y(i,j,k);
        const double dz=geometry.Z(i+di,j+dj,k+dk)-geometry.Z(i,j,k);
        const double spacing=std::sqrt(dx*dx+dy*dy+dz*dz);
        if (!std::isfinite(spacing) || spacing <= 0.0) {
            throw std::runtime_error(
                "Smagorinsky found invalid active-direction spacing.");
        }
        product *= spacing;
        ++dimensions;
    }
    if (dimensions == 0) {
        throw std::runtime_error(
            "Smagorinsky requires at least one active direction.");
    }
    return scale*std::pow(product,1.0/(double)dimensions);
}

void clearEquationCell(
        TransportEquationState& equation,
        int cell) {
    if (equation.variable.empty()) return;
    equation.diffusivity.values()[(size_t)cell]=0.0;
    equation.explicitSource.values()[(size_t)cell]=0.0;
    equation.implicitSink.values()[(size_t)cell]=0.0;
}

void requirePositiveState(
        const std::string& model,
        const std::string& phase,
        const std::string& secondName,
        double density,
        double phaseMass,
        double kineticEnergy,
        double secondVariable,
        double kineticEnergyFloor,
        double secondFloor,
        int i, int j, int k) {
    if (std::isfinite(density) && density>0.0
        && std::isfinite(phaseMass) && phaseMass>0.0
        && std::isfinite(kineticEnergy)
        && kineticEnergy>=kineticEnergyFloor
        && std::isfinite(secondVariable)
        && secondVariable>=secondFloor) {
        return;
    }
    throw std::runtime_error(
        model+" received invalid state for phase '"+phase+"' at ("
        +std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)
        +"): rho="+std::to_string(density)
        +", phaseMass="+std::to_string(phaseMass)
        +", k="+std::to_string(kineticEnergy)
        +", "+secondName+"="+std::to_string(secondVariable)+".");
}

void buildWallDistance(
        const Physics::PhaseSystems::PhaseSystem& system,
        const FDM::TurbulenceConfig& config,
        PhaseEquationState& state) {
    const Field& geometry=system.geometry();
    if (config.wallPatches.empty()) {
        throw std::runtime_error(
            "Eulerian kOmegaSST requires explicit walls (...) in "
            "constant/turbulenceProperties.");
    }
    state.wallDistance.setupLike(
        geometry,
        "wallDistance."+system.phases()[state.phaseIndex].name,
        -1.0);
    std::unordered_set<int> allWallPoints;
    for (const std::string& patch : config.wallPatches) {
        // 仅在活动计算方向识别壁面，二维 empty 厚度方向不能形成候选壁面。
        const int wallAxis=
            StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(geometry,patch);
        const auto& points=geometry.getSet(patch);
        std::unordered_set<int> physicalPoints;
        int side=-1;
        const int lo=geometry.NG();
        const int hi=lo+StructuredMesh::BoundaryGeometry::axisSize(
            geometry,wallAxis)-1;
        for (int point : points) {
            int i=0,j=0,k=0;
            geometry.getIJK(point,i,j,k);
            if (!isPhysical(geometry,i,j,k)) continue;
            const int coordinate=StructuredMesh::BoundaryGeometry::coordinateIndexForAxis(
                i,j,k,wallAxis);
            const int pointSide=coordinate==lo ? lo
                : (coordinate==hi ? hi : -1);
            if (pointSide<0 || (side>=0 && side!=pointSide)) {
                throw std::runtime_error(
                    "kOmegaSST wall patch '"+patch
                    +"' must occupy exactly one structured boundary side.");
            }
            side=pointSide;
            physicalPoints.insert(point);
            allWallPoints.insert(point);
        }
        if (physicalPoints.empty()) {
            throw std::runtime_error(
                "kOmegaSST wall patch '"+patch+"' has no physical points.");
        }
        for (int k=lo;k<lo+geometry.NZ();++k)
            for (int j=lo;j<lo+geometry.NY();++j)
                for (int i=lo;i<lo+geometry.NX();++i) {
                    std::array<int,3> projected{i,j,k};
                    projected[(size_t)wallAxis]=side;
                    const int wall=geometry.getIdx(
                        projected[0],projected[1],projected[2]);
                    if (!physicalPoints.count(wall)) continue;
                    const double dx=geometry.X(i,j,k)
                        -geometry.X(projected[0],projected[1],projected[2]);
                    const double dy=geometry.Y(i,j,k)
                        -geometry.Y(projected[0],projected[1],projected[2]);
                    const double dz=geometry.Z(i,j,k)
                        -geometry.Z(projected[0],projected[1],projected[2]);
                    const int cell=geometry.getIdx(i,j,k);
                    const double distance=std::sqrt(dx*dx+dy*dy+dz*dz);
                    double& current=
                        state.wallDistance.values()[(size_t)cell];
                    if (current<0.0 || distance<current) current=distance;
                }
    }
    const int ng=geometry.NG();
    for (int k=ng;k<ng+geometry.NZ();++k)
        for (int j=ng;j<ng+geometry.NY();++j)
            for (int i=ng;i<ng+geometry.NX();++i) {
                const int cell=geometry.getIdx(i,j,k);
                const double value=state.wallDistance(i,j,k);
                const bool wallPoint=allWallPoints.count(cell)>0;
                if (!std::isfinite(value)
                    || (wallPoint ? value!=0.0 : value<=0.0)) {
                    throw std::runtime_error(
                        "kOmegaSST could not construct a valid wall distance "
                        "at ("+std::to_string(i)+","+std::to_string(j)+","
                        +std::to_string(k)+").");
                }
            }
}

} // namespace SF::Turbulence::ModelOps
