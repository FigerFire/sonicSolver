/// @file SF_phaseBoundary.cpp
/// @brief PhaseSystem 状态、边界、源项与双欧拉物理实现。

#include "SF_phaseBoundary.h"

#include "core/mesh/SF_meshBoundaryGeometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems {
namespace {

template<class T>
const std::vector<int>& points(const Field& field, const BCSetting<T>& bc) {
    const auto it = field.getAllSets().find(bc.name);
    if (it == field.getAllSets().end()) {
        static const std::vector<int> empty;
        return empty;
    }
    return it->second;
}

std::array<int, 3> inward(const Field& field, int axis, int i, int j, int k) {
    std::array<int, 3> p{i, j, k};
    const int lo = field.NG();
    const int hi = lo + StructuredMesh::BoundaryGeometry::axisSize(field, axis) - 1;
    if (p[(size_t)axis] == lo) ++p[(size_t)axis];
    else if (p[(size_t)axis] == hi) --p[(size_t)axis];
    else throw std::runtime_error("Phase BC point is not on its patch plane.");
    return p;
}

void fillGhost(const Field& f, ScalarField& q, int axis,
               int i, int j, int k, double value, bool fixed) {
    const int coord = axis == 0 ? i : (axis == 1 ? j : k);
    const int lo = f.NG();
    const int hi = lo + StructuredMesh::BoundaryGeometry::axisSize(f, axis) - 1;
    const int sign = coord == lo ? -1 : 1;
    for (int layer = 1; layer <= f.NG(); ++layer) {
        std::array<int, 3> g{i, j, k};
        g[(size_t)axis] += sign * layer;
        const double interior = q(i, j, k);
        q(g[0], g[1], g[2]) = fixed ? 2.0 * value - interior : interior;
    }
}

void applyScalar(const Field& f, ScalarField& q,
                 const std::vector<BCSetting<double>>& settings) {
    for (const auto& bc : settings) {
        const auto set = f.getAllSets().find(bc.name);
        if (set == f.getAllSets().end() || set->second.empty()) continue;
        const int axis = bc.type == EMPTY
            ? StructuredMesh::BoundaryGeometry::boundaryAxisForSet(f, bc.name)
            : StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(f, bc.name);
        for (int id : points(f, bc)) {
            int i = 0, j = 0, k = 0;
            f.getIJK(id, i, j, k);
            if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(f, i, j, k)) continue;
            if (bc.type == FIXED_VALUE) {
                if (!std::isfinite(bc.value))
                    throw std::runtime_error("Phase fixedValue is non-finite.");
                q(i,j,k) = bc.value;
                fillGhost(f, q, axis, i,j,k, bc.value, true);
            } else if (bc.type == ZERO_GRADIENT || bc.type == SYMMETRY) {
                const auto s = inward(f, axis, i,j,k);
                q(i,j,k) = q(s[0],s[1],s[2]);
                fillGhost(f, q, axis, i,j,k, q(i,j,k), false);
            } else if (bc.type == EMPTY) {
                fillGhost(f, q, axis, i,j,k, q(i,j,k), false);
            }
        }
    }
}

void applyVelocity(const Field& f, std::array<ScalarField,3>& U,
                   const std::vector<BCSetting<Vector3>>& settings) {
    for (const auto& bc : settings) {
        const auto set = f.getAllSets().find(bc.name);
        if (set == f.getAllSets().end() || set->second.empty()) continue;
        const int axis = bc.type == EMPTY
            ? StructuredMesh::BoundaryGeometry::boundaryAxisForSet(f, bc.name)
            : StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(f, bc.name);
        const auto normal = StructuredMesh::BoundaryGeometry::analyzeBoundarySet(f, bc.name).normal;
        for (int id : points(f, bc)) {
            int i=0,j=0,k=0; f.getIJK(id,i,j,k);
            if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(f,i,j,k)) continue;
            const auto s = inward(f,axis,i,j,k);
            std::array<double,3> value{
                U[0](s[0],s[1],s[2]), U[1](s[0],s[1],s[2]),
                U[2](s[0],s[1],s[2])};
            if (bc.type == FIXED_VALUE) {
                value = {bc.value.x,bc.value.y,bc.value.z};
            } else if (bc.type == SYMMETRY) {
                const double un = value[0]*normal[0]+value[1]*normal[1]
                                + value[2]*normal[2];
                for (int d=0;d<3;++d) value[(size_t)d]-=un*normal[(size_t)d];
            } else if (bc.type == EMPTY) {
                value[(size_t)axis] = 0.0;
            }
            for (int d=0;d<3;++d) {
                if (!std::isfinite(value[(size_t)d]))
                    throw std::runtime_error("Phase velocity BC is non-finite.");
                U[(size_t)d](i,j,k)=value[(size_t)d];
                fillGhost(f,U[(size_t)d],axis,i,j,k,value[(size_t)d],
                          bc.type==FIXED_VALUE);
            }
        }
    }
}

template<class T, class Setter>
void applyInitial(const Field& f, const std::vector<BCSetting<T>>& settings,
                  Setter&& setter) {
    for (const auto& bc : settings) {
        const auto set = f.getAllSets().find(bc.name);
        if (set == f.getAllSets().end() || set->second.empty()) continue;
        if (bc.type != FIXED_VALUE) {
            throw std::runtime_error("Initial phase fields only accept fixedValue.");
        }
        for (int id : points(f,bc)) {
            int i=0,j=0,k=0; f.getIJK(id,i,j,k);
            setter(i,j,k,bc.value);
        }
    }
}

} // namespace

void PhaseBoundaryApplicator::initialize(
        const Field& f, const Multiphase::PhaseProperties& p,
        double initialPressure,
        PhaseState& s) const {
    auto& primitive = s.primitive;
    std::fill(primitive.alpha.values().begin(),
              primitive.alpha.values().end(), p.volumeFraction);
    std::fill(primitive.temperature.values().begin(),
              primitive.temperature.values().end(),
              p.initialTemperature);
    std::fill(primitive.density.values().begin(),
              primitive.density.values().end(),
              Multiphase::phaseDensity(
                  p, p.initialTemperature, initialPressure));
    const std::array<double,3> initialU{
        p.initialVelocity.x,p.initialVelocity.y,p.initialVelocity.z};
    for(int d=0;d<3;++d) {
        std::fill(primitive.velocity[(size_t)d].values().begin(),
                  primitive.velocity[(size_t)d].values().end(),
                  initialU[(size_t)d]);
    }
    applyInitial(f,p.alphaInitial,[&](int i,int j,int k,double v){
        primitive.alpha(i,j,k)=v;
    });
    applyInitial(f,p.densityInitial,[&](int i,int j,int k,double v){
        primitive.density(i,j,k)=v;
    });
    applyInitial(f,p.temperatureInitial,[&](int i,int j,int k,double v){
        primitive.temperature(i,j,k)=v;
    });
    applyInitial(f,p.velocityInitial,[&](int i,int j,int k,const Vector3& v){
        primitive.velocity[0](i,j,k)=v.x;
        primitive.velocity[1](i,j,k)=v.y;
        primitive.velocity[2](i,j,k)=v.z;
    });
    s.commitPrimitiveToPrimary(p);
}

void PhaseBoundaryApplicator::apply(
        const Field& f, const Multiphase::PhaseProperties& p,
        PhaseState& s) const {
    auto& primitive = s.primitive;
    applyScalar(f,primitive.alpha,p.alphaBoundary);
    applyScalar(f,primitive.density,p.densityBoundary);
    applyScalar(f,primitive.temperature,p.temperatureBoundary);
    applyVelocity(f,primitive.velocity,p.velocityBoundary);
    s.commitPrimitiveToPrimary(p);
}

void PhaseBoundaryApplicator::applySharedPressure(
        const Field& field,
        ScalarField& pressure,
        const std::vector<BCSetting<double>>& settings) const {
    applyScalar(field, pressure, settings);
}

} // namespace SF::Physics::PhaseSystems
