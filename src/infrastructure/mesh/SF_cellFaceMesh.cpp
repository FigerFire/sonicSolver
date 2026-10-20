/// @file SF_cellFaceMesh.cpp
/// @brief 结构/多块网格拓扑、度量与加载实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.01-----------*/

#include "SF_cellFaceMesh.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>

namespace SF {

namespace {

std::array<int, 4> sortedFaceKey(std::array<int, 4> ids) {
    std::sort(ids.begin(), ids.end());
    return ids;
}

bool faceInPointSet(const MeshFace& face,
                    const std::vector<unsigned char>& mask) {
    for (int pointId : face.pointIds) {
        if (pointId < 0 || pointId >= (int)mask.size()) return false;
        if (!mask[(std::size_t)pointId]) return false;
    }
    return true;
}

void normalizePointSet(std::vector<int>& indices) {
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}

} // namespace

CellFaceMesh CellFaceMesh::fromStructuredPoints(
        int nx, int ny, int nz,
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::map<std::string, std::vector<int>>& pointSets,
        int blockId,
        int ownerRank) {
    CellFaceMesh mesh;
    const int nPoints = nx * ny * nz;
    if (nx <= 0 || ny <= 0 || nz <= 0 ||
        (int)x.size() != nPoints ||
        x.size() != y.size() ||
        x.size() != z.size()) {
        return mesh;
    }

    mesh.points.reserve((std::size_t)nPoints);
    for (int p = 0; p < nPoints; ++p) {
        mesh.points.push_back({x[(std::size_t)p],
                               y[(std::size_t)p],
                               z[(std::size_t)p]});
    }

    const int nxCells = nx - 1;
    const int nyCells = ny - 1;
    const int nzCells = nz - 1;
    if (nxCells <= 0 || nyCells <= 0 || nzCells <= 0) {
        return mesh;
    }

    mesh.cells.reserve((std::size_t)nxCells * nyCells * nzCells);
    std::map<std::array<int, 4>, int> faceMap;

    auto addFace = [&](MeshCell& cell,
                       int localFace,
                       std::array<int, 4> pointIds,
                       int direction) {
        const std::array<int, 4> key = sortedFaceKey(pointIds);
        auto found = faceMap.find(key);
        if (found != faceMap.end()) {
            MeshFace& face = mesh.faces[(std::size_t)found->second];
            face.neighbourCell = (int)mesh.cells.size();
            face.neighbourRank = ownerRank;
            cell.faceIds[(std::size_t)localFace] = found->second;
            return;
        }

        MeshFace face;
        face.pointIds = pointIds;
        face.ownerCell = (int)mesh.cells.size();
        face.neighbourCell = -1;
        face.direction = direction;
        face.ownerRank = ownerRank;
        face.neighbourRank = -1;
        const int faceId = (int)mesh.faces.size();
        mesh.faces.push_back(face);
        faceMap.emplace(key, faceId);
        cell.faceIds[(std::size_t)localFace] = faceId;
    };

    for (int k = 0; k < nzCells; ++k) {
        for (int j = 0; j < nyCells; ++j) {
            for (int i = 0; i < nxCells; ++i) {
                const int p000 = pointIndex(i,     j,     k,     nx, ny);
                const int p100 = pointIndex(i + 1, j,     k,     nx, ny);
                const int p110 = pointIndex(i + 1, j + 1, k,     nx, ny);
                const int p010 = pointIndex(i,     j + 1, k,     nx, ny);
                const int p001 = pointIndex(i,     j,     k + 1, nx, ny);
                const int p101 = pointIndex(i + 1, j,     k + 1, nx, ny);
                const int p111 = pointIndex(i + 1, j + 1, k + 1, nx, ny);
                const int p011 = pointIndex(i,     j + 1, k + 1, nx, ny);

                MeshCell cell;
                cell.pointIds = {p000, p100, p110, p010,
                                 p001, p101, p111, p011};
                cell.blockId = blockId;
                cell.ijk = {i, j, k};

                addFace(cell, 0, {p000, p010, p011, p001}, 0);
                addFace(cell, 1, {p100, p101, p111, p110}, 0);
                addFace(cell, 2, {p000, p001, p101, p100}, 1);
                addFace(cell, 3, {p010, p110, p111, p011}, 1);
                addFace(cell, 4, {p000, p100, p110, p010}, 2);
                addFace(cell, 5, {p001, p011, p111, p101}, 2);

                mesh.cells.push_back(cell);
            }
        }
    }

    for (const auto& [name, ids] : pointSets) {
        if (name == "all") continue;
        std::vector<unsigned char> mask(mesh.points.size(), 0);
        for (int id : ids) {
            if (id >= 0 && id < (int)mask.size()) mask[(std::size_t)id] = 1;
        }

        MeshPatch patch;
        patch.name = name;
        patch.type = "patch";
        for (int faceId = 0; faceId < (int)mesh.faces.size(); ++faceId) {
            MeshFace& face = mesh.faces[(std::size_t)faceId];
            if (face.neighbourCell >= 0) continue;
            if (!faceInPointSet(face, mask)) continue;
            face.patchId = (int)mesh.patches.size();
            patch.faceIds.push_back(faceId);
        }
        if (!patch.faceIds.empty()) {
            mesh.patches.push_back(std::move(patch));
        }
    }

    return mesh;
}

void CellFaceMesh::clear() {
    points.clear();
    cells.clear();
    faces.clear();
    patches.clear();
}

std::map<std::string, std::vector<int>> CellFaceMesh::patchPointSets(
        bool includeAll) const {
    std::map<std::string, std::vector<int>> sets;
    if (includeAll) {
        std::vector<int>& all = sets["all"];
        all.reserve(points.size());
        for (int pointId = 0; pointId < (int)points.size(); ++pointId) {
            all.push_back(pointId);
        }
    }

    for (const MeshPatch& patch : patches) {
        std::vector<int>& ids = sets[patch.name];
        for (int faceId : patch.faceIds) {
            if (faceId < 0 || faceId >= (int)faces.size()) {
                std::cerr << "[SF FATAL] Patch '" << patch.name
                          << "' references invalid faceId " << faceId
                          << ", faceCount=" << faces.size() << std::endl;
                std::exit(1);
            }
            const MeshFace& face = faces[(std::size_t)faceId];
            for (int pointId : face.pointIds) {
                if (pointId < 0 || pointId >= (int)points.size()) {
                    std::cerr << "[SF FATAL] Face " << faceId
                              << " references invalid pointId " << pointId
                              << ", pointCount=" << points.size() << std::endl;
                    std::exit(1);
                }
                ids.push_back(pointId);
            }
        }
        normalizePointSet(ids);
    }

    return sets;
}

int CellFaceMesh::pointIndex(int i, int j, int k, int nx, int ny) {
    return (k * ny + j) * nx + i;
}

} // namespace SF
