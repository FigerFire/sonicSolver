#include "SF_scalarBoundary.h"
#include "core/mesh/SF_meshBoundaryGeometry.h"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdexcept>
namespace SF::Boundary {
void applyNodeScalarBoundary(const Field& geometry,
                   const std::vector<BCSetting<double>>& settings,
                   ScalarField& field) {
    if (settings.empty()) {
        throw std::runtime_error(
            "Node scalar boundary requires boundaryField entries.");
    }
    if(field.values().size()!=static_cast<size_t>(geometry.TotalSize()))
        throw std::runtime_error("Scalar boundary layout differs from geometry.");
    const int ng = geometry.NG();
    const std::array<int,3> lower{ng, ng, ng};
    const std::array<int,3> upper{
        ng + geometry.NX() - 1,
        ng + geometry.NY() - 1,
        ng + geometry.NZ() - 1};
    for (const auto& setting : settings) {
        if(setting.type!=FIXED_VALUE && setting.type!=ZERO_GRADIENT
            && setting.type!=SYMMETRY && setting.type!=EMPTY)
            throw std::runtime_error("Unsupported node scalar boundary on '"+setting.name+"'.");
        if(setting.type==FIXED_VALUE && !std::isfinite(setting.value))
            throw std::runtime_error("Non-finite fixed node scalar boundary.");
        const int axis = setting.type == EMPTY
            ? StructuredMesh::BoundaryGeometry::boundaryAxisForSet(
                geometry, setting.name)
            : StructuredMesh::BoundaryGeometry::activeBoundaryAxisForSet(
                geometry, setting.name);
        const auto& indices = geometry.getSet(setting.name);
        for (int cell : indices) {
            int i = 0, j = 0, k = 0;
            geometry.getIJK(cell, i, j, k);
            if (!StructuredMesh::BoundaryGeometry::isPhysicalPoint(
                    geometry, i, j, k)) {
                continue;
            }
            std::array<int,3> point{i,j,k};
            const int sign = point[(size_t)axis]
                == lower[(size_t)axis] ? -1 : 1;
            std::array<int,3> adjacent = point;
            adjacent[(size_t)axis] -= sign;
            const int adjacentCell = geometry.getIdx(
                adjacent[0], adjacent[1], adjacent[2]);
            if (setting.type == FIXED_VALUE) {
                field.values()[(size_t)cell] = setting.value;
            } else {
                field.values()[(size_t)cell] =
                    field.values()[(size_t)adjacentCell];
            }
            for (int layer = 1; layer <= ng; ++layer) {
                std::array<int,3> ghost = point;
                ghost[(size_t)axis] += sign * layer;
                const int ghostCell = geometry.getIdx(
                    ghost[0], ghost[1], ghost[2]);
                if (setting.type == FIXED_VALUE) {
                    std::array<int,3> mirror = point;
                    mirror[(size_t)axis] -= sign * layer;
                    const int mirrorCell = geometry.getIdx(
                        mirror[0], mirror[1], mirror[2]);
                    field.values()[(size_t)ghostCell] =
                        2.0 * setting.value
                        - field.values()[(size_t)mirrorCell];
                } else {
                    field.values()[(size_t)ghostCell] =
                        field.values()[(size_t)cell];
                }
            }
        }
    }
}

}

namespace SF::Boundary {
namespace {
int clampToInterior(int v,int lo,int hi){return std::max(lo,std::min(v,hi));}
void setFixedScalar(const Field& flow, std::vector<double>& state,
                    int i, int j, int k,
                    double bcValue,
                    int ng, int nx, int ny, int nz) {
    auto applyFixed = [&](int gI, int gJ, int gK, int rI, int rJ, int rK) {
        double realValue = state[flow.getIdx(rI, rJ, rK)];
        state[flow.getIdx(gI, gJ, gK)] = 2.0 * bcValue - realValue;
    };

    bool isGhost = (i < ng || i >= (nx + ng)
                 || j < ng || j >= (ny + ng)
                 || k < ng || k >= (nz + ng));

    if (isGhost) {
        int ri = clampToInterior(i, ng, nx + ng - 1);
        int rj = clampToInterior(j, ng, ny + ng - 1);
        int rk = clampToInterior(k, ng, nz + ng - 1);
        applyFixed(i, j, k, ri, rj, rk);
        return;
    }

    state[flow.getIdx(i, j, k)] = bcValue;

    if (i == ng) {
        for (int b = 0; b < ng; ++b) applyFixed(b, j, k, i, j, k);
    }
    if (i == nx + ng - 1) {
        for (int b = 1; b <= ng; ++b) applyFixed(i + b, j, k, i, j, k);
    }
    if (j == ng) {
        for (int b = 0; b < ng; ++b) applyFixed(i, b, k, i, j, k);
    }
    if (j == ny + ng - 1) {
        for (int b = 1; b <= ng; ++b) applyFixed(i, j + b, k, i, j, k);
    }
    if (k == ng) {
        for (int b = 0; b < ng; ++b) applyFixed(i, j, b, i, j, k);
    }
    if (k == nz + ng - 1) {
        for (int b = 1; b <= ng; ++b) applyFixed(i, j, k + b, i, j, k);
    }
}

void setZeroGradientScalar(const Field& flow, std::vector<double>& state,
                           int i, int j, int k,
                           int ng, int nx, int ny, int nz) {
    bool isGhost = (i < ng || i >= (nx + ng)
                 || j < ng || j >= (ny + ng)
                 || k < ng || k >= (nz + ng));

    if (isGhost) {
        int ri = clampToInterior(i, ng, nx + ng - 1);
        int rj = clampToInterior(j, ng, ny + ng - 1);
        int rk = clampToInterior(k, ng, nz + ng - 1);
        state[flow.getIdx(i, j, k)] = state[flow.getIdx(ri, rj, rk)];
        return;
    }

    int srcI = i;
    int srcJ = j;
    int srcK = k;
    if (i == ng && nx > 1) srcI = i + 1;
    else if (i == nx + ng - 1 && nx > 1) srcI = i - 1;
    if (j == ng && ny > 1) srcJ = j + 1;
    else if (j == ny + ng - 1 && ny > 1) srcJ = j - 1;
    if (k == ng && nz > 1) srcK = k + 1;
    else if (k == nz + ng - 1 && nz > 1) srcK = k - 1;

    state[flow.getIdx(i, j, k)] = state[flow.getIdx(srcI, srcJ, srcK)];

    if (i == ng) {
        for (int b = 0; b < ng; ++b) state[flow.getIdx(b, j, k)] = state[flow.getIdx(i, j, k)];
    }
    if (i == nx + ng - 1) {
        for (int b = 1; b <= ng; ++b) state[flow.getIdx(i + b, j, k)] = state[flow.getIdx(i, j, k)];
    }
    if (j == ng) {
        for (int b = 0; b < ng; ++b) state[flow.getIdx(i, b, k)] = state[flow.getIdx(i, j, k)];
    }
    if (j == ny + ng - 1) {
        for (int b = 1; b <= ng; ++b) state[flow.getIdx(i, j + b, k)] = state[flow.getIdx(i, j, k)];
    }
    if (k == ng) {
        for (int b = 0; b < ng; ++b) state[flow.getIdx(i, j, b)] = state[flow.getIdx(i, j, k)];
    }
    if (k == nz + ng - 1) {
        for (int b = 1; b <= ng; ++b) state[flow.getIdx(i, j, k + b)] = state[flow.getIdx(i, j, k)];
    }
}

void setEmptyScalar(const Field& flow, std::vector<double>& state,
                    int i, int j, int k,
                    int ng, int nx, int ny, int nz, int axis) {
    int di = 0, dj = 0, dk = 0;
    if (axis == 0) di = 1;
    else if (axis == 1) dj = 1;
    else dk = 1;

    auto fill = [&](int gi, int gj, int gk) {
        state[flow.getIdx(gi, gj, gk)] = state[flow.getIdx(i, j, k)];
    };

    if (axis == 0 && i == ng) {
        for (int b = 1; b <= ng; ++b) fill(i - b * di, j, k);
    }
    if (axis == 0 && i == nx + ng - 1) {
        for (int b = 1; b <= ng; ++b) fill(i + b * di, j, k);
    }
    if (axis == 1 && j == ng) {
        for (int b = 1; b <= ng; ++b) fill(i, j - b * dj, k);
    }
    if (axis == 1 && j == ny + ng - 1) {
        for (int b = 1; b <= ng; ++b) fill(i, j + b * dj, k);
    }
    if (axis == 2 && k == ng) {
        for (int b = 1; b <= ng; ++b) fill(i, j, k - b * dk);
    }
    if (axis == 2 && k == nz + ng - 1) {
        for (int b = 1; b <= ng; ++b) fill(i, j, k + b * dk);
    }
}

}
void applyTransportScalarBoundary(const Field& flow,
                                   const std::vector<BCSetting<double>>& settings,
                                   std::vector<double>& state) {
    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();

    for (const auto& bc : settings) {
        auto applyAt = [&](int i, int j, int k, int axis) {

            switch (bc.type) {
                case FIXED_VALUE:
                    setFixedScalar(flow, state, i, j, k, bc.value, ng, nx, ny, nz);
                    break;
                case ZERO_GRADIENT:
                case SYMMETRY:
                    setZeroGradientScalar(flow, state, i, j, k, ng, nx, ny, nz);
                    break;
                case EMPTY:
                    setEmptyScalar(flow, state, i, j, k, ng, nx, ny, nz, axis);
                    break;
            }
        };

        const auto& allSets = flow.getAllSets();
        if (allSets.find(bc.name) == allSets.end()) {
            std::cerr << "[SF FATAL] Turbulence boundary zone '" << bc.name
                      << "' is not a mesh set." << std::endl;
            std::exit(1);
        }

        const auto& indices = flow.getSet(bc.name);
        const int axis = (bc.type == EMPTY)
            ? StructuredMesh::BoundaryGeometry::boundaryAxisForSet(flow, bc.name)
            : -1;
        for (int flatIndex : indices) {
            int i, j, k;
            flow.getIJK(flatIndex, i, j, k);
            applyAt(i, j, k, axis);
        }
    }
}

}
