// Frozen pre-migration stencil, test oracle only.
#include "core/mesh/SF_meshBoundaryGeometry.h"
#include <array>
#include <stdexcept>
void referenceBoundary(const Field& geometry,
                   const std::vector<BCSetting<double>>& settings,
                   ScalarField& field) {
    if (settings.empty()) {
        throw std::runtime_error(
            "Transported turbulence equation requires boundaryField entries.");
    }
    const int ng = geometry.NG();
    const std::array<int,3> lower{ng, ng, ng};
    const std::array<int,3> upper{
        ng + geometry.NX() - 1,
        ng + geometry.NY() - 1,
        ng + geometry.NZ() - 1};
    for (const auto& setting : settings) {
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

