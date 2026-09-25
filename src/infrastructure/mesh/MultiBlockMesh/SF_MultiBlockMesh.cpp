/// @file SF_MultiBlockMesh.cpp
/// @brief 结构/多块网格拓扑、度量与加载实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_MultiBlockMesh.h"
#include "SF_meshDecompose.h"
#include "SF_mesh.h"
#include "SF_phaseProperties.h"
#include "core/interfaces/SF_log.h"
#include "methods/numerics/structured/SF_structured.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace SF {

namespace {

static void cleanLine(std::string& line) {
    size_t p = line.find("//");
    if (p != std::string::npos) line = line.substr(0, p);
    p = line.find('#');
    if (p != std::string::npos && (p == 0 || line[p - 1] != '\\')) {
        // Keep SFM section headers intact.
        if (p != 0) line = line.substr(0, p);
    }

    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), notSpace));
    line.erase(std::find_if(line.rbegin(), line.rend(), notSpace).base(), line.end());
}

static void trimOnly(std::string& line) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), notSpace));
    line.erase(std::find_if(line.rbegin(), line.rend(), notSpace).base(), line.end());
}

static std::string upperSection(std::string section) {
    std::transform(section.begin(), section.end(), section.begin(),
                   [](unsigned char c) { return (char)std::toupper(c); });
    return section;
}

static bool isReservedSection(const std::string& section,
                              const std::string& expectedUpper) {
    return upperSection(section) == expectedUpper;
}

static void normalizePointSet(std::vector<int>& indices) {
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}

static void normalizePointSets(RawMeshBlock& raw) {
    for (auto& [name, indices] : raw.pointSets) {
        normalizePointSet(indices);
    }
}

static void appendPointSet(std::vector<int>& dst, const std::vector<int>& src) {
    dst.insert(dst.end(), src.begin(), src.end());
    normalizePointSet(dst);
}

static bool isBuiltinAllSetName(const std::string& name) {
    if (name.size() != 3) return false;
    return std::tolower((unsigned char)name[0]) == 'a'
        && std::tolower((unsigned char)name[1]) == 'l'
        && std::tolower((unsigned char)name[2]) == 'l';
}

static std::vector<int> canonicalAllPointSet(int pointCount) {
    std::vector<int> all((std::size_t)pointCount);
    std::iota(all.begin(), all.end(), 0);
    return all;
}

static bool ensureBuiltinAllPointSet(RawMeshBlock& raw) {
    const int totalPts = raw.nx * raw.ny * raw.nz;
    if (totalPts <= 0) return true;

    const std::vector<int> all = canonicalAllPointSet(totalPts);
    std::vector<std::string> allAliases;
    for (auto& [name, indices] : raw.pointSets) {
        if (!isBuiltinAllSetName(name)) continue;
        normalizePointSet(indices);
        if (!indices.empty() && indices != all) {
            broadcast("Fatal: ",
                      "Set '" + name + "' in " + raw.name
                      + " uses reserved built-in name 'all' but does not "
                      "cover every physical point.");
            return false;
        }
        allAliases.push_back(name);
    }

    raw.pointSets["all"] = all;
    for (const std::string& alias : allAliases) {
        raw.pointSets[alias] = all;
    }
    return true;
}

static bool ensureBuiltinAllPointSets(std::vector<RawMeshBlock>& blocks) {
    for (RawMeshBlock& raw : blocks) {
        if (!ensureBuiltinAllPointSet(raw)) return false;
    }
    return true;
}

struct SourceMetricData {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::vector<std::array<double, 10>> values;
    std::array<std::vector<std::array<double, 4>>, 3> faces;
};

static size_t sourceFaceIndex(
        const SourceMetricData& source,
        int direction, int i, int j, int k) {
    if (direction == 0) {
        if (i < -1 || i >= source.nx ||
            j < 0 || j >= source.ny ||
            k < 0 || k >= source.nz) return (size_t)-1;
        return ((size_t)k * source.ny + (size_t)j)
            * (size_t)(source.nx + 1) + (size_t)(i + 1);
    }
    if (direction == 1) {
        if (i < 0 || i >= source.nx ||
            j < -1 || j >= source.ny ||
            k < 0 || k >= source.nz) return (size_t)-1;
        return ((size_t)k * (size_t)(source.ny + 1)
                + (size_t)(j + 1)) * (size_t)source.nx + (size_t)i;
    }
    if (direction == 2) {
        if (i < 0 || i >= source.nx ||
            j < 0 || j >= source.ny ||
            k < -1 || k >= source.nz) return (size_t)-1;
        return ((size_t)(k + 1) * source.ny + (size_t)j)
            * (size_t)source.nx + (size_t)i;
    }
    return (size_t)-1;
}

static std::array<double, 10> metricAt(
        const Field& field, int i, int j, int k) {
    return {
        field.Jac(i, j, k),
        field.XiX(i, j, k), field.XiY(i, j, k), field.XiZ(i, j, k),
        field.EtX(i, j, k), field.EtY(i, j, k), field.EtZ(i, j, k),
        field.ZeX(i, j, k), field.ZeY(i, j, k), field.ZeZ(i, j, k)};
}

static void setMetricAt(
        Field& field, int i, int j, int k,
        const std::array<double, 10>& metric) {
    field.Jac(i, j, k) = metric[0];
    field.XiX(i, j, k) = metric[1];
    field.XiY(i, j, k) = metric[2];
    field.XiZ(i, j, k) = metric[3];
    field.EtX(i, j, k) = metric[4];
    field.EtY(i, j, k) = metric[5];
    field.EtZ(i, j, k) = metric[6];
    field.ZeX(i, j, k) = metric[7];
    field.ZeY(i, j, k) = metric[8];
    field.ZeZ(i, j, k) = metric[9];
}

static bool buildSourceCanonicalMetrics(
        const std::vector<RawMeshBlock>& zones,
        const MeshRuntimeConfig& config,
        std::vector<SourceMetricData>& metrics) {
    metrics.clear();
    metrics.reserve(zones.size());
    for (const RawMeshBlock& zone : zones) {
        Field source;
        if (!Mesh::setupFieldFromRaw(
                source, zone.nx, zone.ny, zone.nz, zone.ng,
                zone.x, zone.y, zone.z, zone.pointSets,
                config, zone.name + ":canonicalMetric")) {
            return false;
        }
        SourceMetricData data;
        data.nx = zone.nx;
        data.ny = zone.ny;
        data.nz = zone.nz;
        data.values.resize((size_t)zone.nx * zone.ny * zone.nz);
        const int ng = source.NG();
        size_t index = 0;
        for (int k = 0; k < zone.nz; ++k) {
            for (int j = 0; j < zone.ny; ++j) {
                for (int i = 0; i < zone.nx; ++i, ++index) {
                    data.values[index] = metricAt(
                        source, i + ng, j + ng, k + ng);
                }
            }
        }
        data.faces[0].resize(
            (size_t)(zone.nx + 1) * zone.ny * zone.nz);
        data.faces[1].resize(
            (size_t)zone.nx * (zone.ny + 1) * zone.nz);
        data.faces[2].resize(
            (size_t)zone.nx * zone.ny * (zone.nz + 1));
        for (int direction = 0; direction < 3; ++direction) {
            const Math::Dir dir = (Math::Dir)direction;
            if (!Math::isDirectionActive(dir)) continue;
            Math::forFaces(source, dir, [&](int i, int j, int k) {
                const int li = i - ng;
                const int lj = j - ng;
                const int lk = k - ng;
                const size_t face = sourceFaceIndex(
                    data, direction, li, lj, lk);
                if (face == (size_t)-1 ||
                    face >= data.faces[(size_t)direction].size()) {
                    return;
                }
                double metric[4]{0.0, 0.0, 0.0, 0.0};
                Math::faceMetrics(source, i, j, k, dir, metric);
                data.faces[(size_t)direction][face] = {
                    metric[0], metric[1], metric[2], metric[3]};
            });
        }
        metrics.push_back(std::move(data));
    }
    broadcast("Canonical source metrics: ",
              std::to_string(metrics.size())
              + " source zone(s) evaluated before MPI slicing");
    return true;
}

static bool applySourceCanonicalMetrics(
        std::vector<MeshBlockField>& patches,
        const std::vector<MeshDecompose::SourcePatchExtent>& extents,
        const std::vector<SourceMetricData>& metrics) {
    if (patches.size() != extents.size()) return false;
    for (const auto& extent : extents) {
        if (extent.patchId < 0 ||
            extent.patchId >= (int)patches.size() ||
            extent.sourceZoneId < 0 ||
            extent.sourceZoneId >= (int)metrics.size()) {
            return false;
        }
        MeshBlockField& patch = patches[(size_t)extent.patchId];
        const SourceMetricData& source =
            metrics[(size_t)extent.sourceZoneId];
        if (extent.sourceSize !=
            std::array<int, 3>{source.nx, source.ny, source.nz}) {
            return false;
        }
        Field& field = patch.field;
        const int ng = field.NG();
        for (int k = 0; k < field.NZ(); ++k) {
            for (int j = 0; j < field.NY(); ++j) {
                for (int i = 0; i < field.NX(); ++i) {
                    const int si = extent.start[0] + i;
                    const int sj = extent.start[1] + j;
                    const int sk = extent.start[2] + k;
                    const size_t sourceIndex =
                        ((size_t)sk * source.ny + (size_t)sj)
                        * source.nx + (size_t)si;
                    if (sourceIndex >= source.values.size()) return false;
                    setMetricAt(field, i + ng, j + ng, k + ng,
                                source.values[sourceIndex]);
                }
            }
        }
        for (int direction = 0; direction < 3; ++direction) {
            const Math::Dir dir = (Math::Dir)direction;
            if (!Math::isDirectionActive(dir)) continue;
            bool valid = true;
            Math::forFaces(field, dir, [&](int i, int j, int k) {
                const int si = extent.start[0] + i - ng;
                const int sj = extent.start[1] + j - ng;
                const int sk = extent.start[2] + k - ng;
                const size_t sourceFace = sourceFaceIndex(
                    source, direction, si, sj, sk);
                if (sourceFace == (size_t)-1 ||
                    sourceFace >= source.faces[(size_t)direction].size()) {
                    valid = false;
                    return;
                }
                field.setCanonicalFaceMetrics(
                    direction, i, j, k,
                    source.faces[(size_t)direction][sourceFace]);
            });
            if (!valid) return false;
        }
    }
    return true;
}

static bool appendCellFace(MeshCell& cell, int faceId) {
    for (int& slot : cell.faceIds) {
        if (slot < 0) {
            slot = faceId;
            return true;
        }
    }
    return false;
}

static bool finalizeParsedTopology(RawMeshBlock& raw,
                                   int blockId,
                                   std::string* error) {
    if (raw.topology.empty()) return true;

    const int nPoints = raw.nx * raw.ny * raw.nz;
    if (nPoints <= 0) {
        if (error) *error = "topology sections appear before valid #Information.";
        return false;
    }
    if ((int)raw.x.size() != nPoints ||
        raw.x.size() != raw.y.size() ||
        raw.x.size() != raw.z.size()) {
        if (error) *error = "topology sections require matching #Point data.";
        return false;
    }

    raw.topology.points.clear();
    raw.topology.points.reserve((std::size_t)nPoints);
    for (int p = 0; p < nPoints; ++p) {
        raw.topology.points.push_back({raw.x[(std::size_t)p],
                                       raw.y[(std::size_t)p],
                                       raw.z[(std::size_t)p]});
    }

    const int nCells = (int)raw.topology.cells.size();
    const int nxCells = raw.nx - 1;
    const int nyCells = raw.ny - 1;
    const int nzCells = raw.nz - 1;
    if (nxCells <= 0 || nyCells <= 0 || nzCells <= 0) {
        if (error) *error = "cell-face topology requires at least one cell.";
        return false;
    }
    if (nCells != nxCells * nyCells * nzCells) {
        if (error) {
            *error = "cell count is " + std::to_string(nCells)
                   + ", expected "
                   + std::to_string(nxCells * nyCells * nzCells) + ".";
        }
        return false;
    }
    if (raw.topology.faces.empty()) {
        if (error) *error = "cell-face topology contains no #Face entries.";
        return false;
    }

    for (int cellId = 0; cellId < nCells; ++cellId) {
        MeshCell& cell = raw.topology.cells[(std::size_t)cellId];
        for (int pointId : cell.pointIds) {
            if (pointId < 0 || pointId >= nPoints) {
                if (error) {
                    *error = "cell " + std::to_string(cellId)
                           + " references invalid point "
                           + std::to_string(pointId) + ".";
                }
                return false;
            }
        }
        cell.faceIds = {-1, -1, -1, -1, -1, -1};
        cell.blockId = blockId;
        const int i = cellId % nxCells;
        const int j = (cellId / nxCells) % nyCells;
        const int k = cellId / (nxCells * nyCells);
        cell.ijk = {i, j, k};
    }

    for (int faceId = 0; faceId < (int)raw.topology.faces.size(); ++faceId) {
        MeshFace& face = raw.topology.faces[(std::size_t)faceId];
        for (int pointId : face.pointIds) {
            if (pointId < 0 || pointId >= nPoints) {
                if (error) {
                    *error = "face " + std::to_string(faceId)
                           + " references invalid point "
                           + std::to_string(pointId) + ".";
                }
                return false;
            }
        }
        if (face.ownerCell < 0 || face.ownerCell >= nCells) {
            if (error) {
                *error = "face " + std::to_string(faceId)
                       + " has invalid owner cell "
                       + std::to_string(face.ownerCell) + ".";
            }
            return false;
        }
        if (face.neighbourCell < -1 || face.neighbourCell >= nCells) {
            if (error) {
                *error = "face " + std::to_string(faceId)
                       + " has invalid neighbour cell "
                       + std::to_string(face.neighbourCell) + ".";
            }
            return false;
        }
        if (face.patchId >= (int)raw.topology.patches.size()) {
            if (error) {
                *error = "face " + std::to_string(faceId)
                       + " references invalid patch "
                       + std::to_string(face.patchId) + ".";
            }
            return false;
        }
        if (!appendCellFace(raw.topology.cells[(std::size_t)face.ownerCell],
                            faceId)) {
            if (error) {
                *error = "owner cell " + std::to_string(face.ownerCell)
                       + " has more than six faces.";
            }
            return false;
        }
        if (face.neighbourCell >= 0 &&
            !appendCellFace(raw.topology.cells[(std::size_t)face.neighbourCell],
                            faceId)) {
            if (error) {
                *error = "neighbour cell " + std::to_string(face.neighbourCell)
                       + " has more than six faces.";
            }
            return false;
        }
    }

    for (int patchId = 0; patchId < (int)raw.topology.patches.size(); ++patchId) {
        MeshPatch& patch = raw.topology.patches[(std::size_t)patchId];
        for (int faceId : patch.faceIds) {
            if (faceId < 0 || faceId >= (int)raw.topology.faces.size()) {
                if (error) {
                    *error = "patch '" + patch.name
                           + "' references invalid face "
                           + std::to_string(faceId) + ".";
                }
                return false;
            }
            MeshFace& face = raw.topology.faces[(std::size_t)faceId];
            if (face.patchId >= 0 && face.patchId != patchId) {
                if (error) {
                    *error = "face " + std::to_string(faceId)
                           + " belongs to multiple patches.";
                }
                return false;
            }
            face.patchId = patchId;
        }
    }

    for (int cellId = 0; cellId < nCells; ++cellId) {
        const MeshCell& cell = raw.topology.cells[(std::size_t)cellId];
        for (int faceId : cell.faceIds) {
            if (faceId < 0) {
                if (error) {
                    *error = "cell " + std::to_string(cellId)
                           + " is missing one or more faces.";
                }
                return false;
            }
        }
    }

    auto patchSets = raw.topology.patchPointSets(true);
    for (const auto& [name, indices] : patchSets) {
        appendPointSet(raw.pointSets[name], indices);
    }
    return true;
}

bool configuredSetHasPoints(const std::vector<RawMeshBlock>& blocks,
                            const std::string& name) {
    if (isBuiltinAllSetName(name)) return true;
    for (const RawMeshBlock& raw : blocks) {
        auto it = raw.pointSets.find(name);
        if (it != raw.pointSets.end() && !it->second.empty()) return true;
    }
    return false;
}

bool validateConfiguredSetReferences(const std::vector<RawMeshBlock>& blocks,
                                     const std::vector<std::string>& names) {
    for (const std::string& name : names) {
        if (configuredSetHasPoints(blocks, name)) continue;
        broadcast("Fatal: ",
                  "case references set '" + name
                  + "', but no loaded mesh block contains any point in that set.");
        return false;
    }
    return true;
}

bool validateConfiguredPointCoverage(const RawMeshBlock& raw,
                                     const std::vector<std::string>& names) {
    const int totalPts = raw.nx * raw.ny * raw.nz;
    std::vector<unsigned char> covered((size_t)totalPts, 0);

    for (const std::string& name : names) {
        if (isBuiltinAllSetName(name)) {
            std::fill(covered.begin(), covered.end(), 1);
            continue;
        }
        auto it = raw.pointSets.find(name);
        if (it == raw.pointSets.end()) continue;
        for (int idx : it->second) {
            if (idx >= 0 && idx < totalPts) covered[(size_t)idx] = 1;
        }
    }

    std::vector<int> missing;
    for (int idx = 0; idx < totalPts; ++idx) {
        if (!covered[(size_t)idx]) missing.push_back(idx);
    }
    if (missing.empty()) return true;

    std::ostringstream oss;
    oss << missing.size()
        << " point(s) in " << raw.name
        << " are not covered by any configured IC/BC referenced set. First indices:";
    const int show = std::min((int)missing.size(), 12);
    for (int n = 0; n < show; ++n) {
        const int idx = missing[(size_t)n];
        const int i = idx % raw.nx;
        const int j = (idx / raw.nx) % raw.ny;
        const int k = idx / (raw.nx * raw.ny);
        oss << " " << idx << "(" << i << "," << j << "," << k << ")";
    }
    if ((int)missing.size() > show) oss << " ...";
    broadcast("Fatal: ", oss.str());
    return false;
}

bool validateConfiguredSetCoverage(
        const std::vector<RawMeshBlock>& blocks,
        const std::vector<std::string>& names) {
    if (names.empty()) {
        broadcast("Fatal: ",
                  "No IC/BC set is configured. Every physical point must be "
                  "covered by at least one configured IC/BC referenced set.");
        return false;
    }
    if (!validateConfiguredSetReferences(blocks, names)) return false;
    for (const RawMeshBlock& raw : blocks) {
        if (!validateConfiguredPointCoverage(raw, names)) return false;
    }
    return true;
}

struct SourceBoundaryFaceRef {
    int blockId = -1;
    int axis = -1;
    int side = -1;
    std::array<int, 4> globalPointIds{-1, -1, -1, -1};
};

static int sourceLocalCoord(int pointId, int nx, int ny, int axis) {
    if (axis == 0) return pointId % nx;
    if (axis == 1) return (pointId / nx) % ny;
    return pointId / (nx * ny);
}

static int sourceBoundarySide(const RawMeshBlock& raw,
                              const MeshFace& face) {
    if (face.direction < 0 || face.direction > 2) return -1;
    const int maxCoord =
        face.direction == 0 ? raw.nx - 1 :
        face.direction == 1 ? raw.ny - 1 : raw.nz - 1;
    int lowerCount = 0;
    int upperCount = 0;
    for (int pointId : face.pointIds) {
        const int coord =
            sourceLocalCoord(pointId, raw.nx, raw.ny, face.direction);
        if (coord == 0) ++lowerCount;
        if (coord == maxCoord) ++upperCount;
    }
    if (lowerCount == 4) return 0;
    if (upperCount == 4) return 1;
    return -1;
}

static std::string sourceSideLabel(const SourceBoundaryFaceRef& ref) {
    const char axisName = ref.axis == 0 ? 'i' : (ref.axis == 1 ? 'j' : 'k');
    std::ostringstream oss;
    oss << "block" << ref.blockId << "." << axisName
        << (ref.side == 0 ? "Min" : "Max");
    return oss.str();
}

static std::string sourcePairLabel(const SourceBoundaryFaceRef& a,
                                   const SourceBoundaryFaceRef& b) {
    std::string lhs = sourceSideLabel(a);
    std::string rhs = sourceSideLabel(b);
    if (rhs < lhs) std::swap(lhs, rhs);
    return lhs + " <-> " + rhs;
}

static bool reportCompositeSourceInterfaces(
        const std::vector<RawMeshBlock>& rawBlocks) {
    if (rawBlocks.size() <= 1) return true;

    std::map<int, int> pointReferences;
    std::map<int, std::set<int>> pointZones;
    for (size_t blockId = 0; blockId < rawBlocks.size(); ++blockId) {
        const RawMeshBlock& raw = rawBlocks[blockId];
        for (int globalId : raw.globalPointIds) {
            if (globalId < 0) continue;
            ++pointReferences[globalId];
            pointZones[globalId].insert((int)blockId);
        }
    }

    int sharedPointGroups = 0;
    int branchPointGroups = 0;
    int maxPointMultiplicity = 1;
    for (const auto& [globalId, references] : pointReferences) {
        (void)globalId;
        maxPointMultiplicity = std::max(maxPointMultiplicity, references);
        auto zones = pointZones.find(globalId);
        const int zoneCount =
            zones == pointZones.end() ? 0 : (int)zones->second.size();
        if (zoneCount > 1) ++sharedPointGroups;
        if (references >= 3 && zoneCount > 1) ++branchPointGroups;
    }

    std::map<std::array<int, 4>, std::vector<SourceBoundaryFaceRef>> faceGroups;
    for (size_t blockId = 0; blockId < rawBlocks.size(); ++blockId) {
        const RawMeshBlock& raw = rawBlocks[blockId];
        const size_t expected = (size_t)raw.nx * raw.ny * raw.nz;
        if (raw.topology.empty() || raw.globalPointIds.size() != expected) {
            broadcast("Fatal: ",
                      "multi-zone source interface detection requires "
                      "cell-face topology and globalPointIds for block "
                      + std::to_string(blockId));
            return false;
        }

        for (const MeshFace& face : raw.topology.faces) {
            if (face.neighbourCell >= 0) continue;

            SourceBoundaryFaceRef ref;
            ref.blockId = (int)blockId;
            ref.axis = face.direction;
            ref.side = sourceBoundarySide(raw, face);
            if (ref.side < 0) {
                broadcast("Fatal: ",
                          "boundary face in source block "
                          + std::to_string(blockId)
                          + " cannot be mapped to a structured side.");
                return false;
            }

            bool valid = true;
            for (int n = 0; n < 4; ++n) {
                const int pointId = face.pointIds[(size_t)n];
                if (pointId < 0 || pointId >= (int)expected) {
                    valid = false;
                    break;
                }
                const int globalId = raw.globalPointIds[(size_t)pointId];
                if (globalId < 0) {
                    valid = false;
                    break;
                }
                ref.globalPointIds[(size_t)n] = globalId;
            }
            if (!valid) {
                broadcast("Fatal: ",
                          "source interface face references invalid point ids "
                          "in block " + std::to_string(blockId));
                return false;
            }

            auto key = ref.globalPointIds;
            std::sort(key.begin(), key.end());
            faceGroups[key].push_back(ref);
        }
    }

    int sharedFaceGroups = 0;
    int invalidFaceGroups = 0;
    int maxFaceMultiplicity = 1;
    std::map<std::string, int> pairFaceCounts;
    for (const auto& [key, refs] : faceGroups) {
        (void)key;
        if (refs.size() < 2) continue;

        std::set<int> blocks;
        for (const SourceBoundaryFaceRef& ref : refs) {
            blocks.insert(ref.blockId);
        }
        if (blocks.size() < 2) continue;

        maxFaceMultiplicity =
            std::max(maxFaceMultiplicity, (int)refs.size());
        if (refs.size() != 2 || blocks.size() != 2) {
            ++invalidFaceGroups;
            continue;
        }

        ++sharedFaceGroups;
        ++pairFaceCounts[sourcePairLabel(refs[0], refs[1])];
    }

    if (invalidFaceGroups > 0) {
        broadcast("Fatal: ",
                  "source mesh contains boundary faces shared by more than "
                  "two blocks; conservative FDM needs an explicit pairwise "
                  "interface.");
        return false;
    }
    if (sharedPointGroups > 0 && sharedFaceGroups == 0) {
        broadcast("Fatal: ",
                  "source mesh has repeated physical points but no conformal "
                  "shared boundary faces.");
        return false;
    }

    {
        std::ostringstream oss;
        oss << sharedFaceGroups << " shared boundary face groups, "
            << sharedPointGroups << " shared point groups, "
            << "branch point groups=" << branchPointGroups
            << ", max point multiplicity=" << maxPointMultiplicity
            << ", max face multiplicity=" << maxFaceMultiplicity;
        broadcast("Mesh source interfaces: ", oss.str());
    }
    if (!pairFaceCounts.empty()) {
        std::ostringstream oss;
        bool first = true;
        for (const auto& [label, count] : pairFaceCounts) {
            if (!first) oss << "; ";
            first = false;
            oss << label << "=" << count;
        }
        broadcast("Mesh source interface pairs: ", oss.str());
    }
    return true;
}

} // namespace

bool MultiBlockMesh::loadFiles(const std::string& caseDir,
                               const std::vector<std::string>& meshFiles,
                               const MeshRuntimeConfig& config,
                               const CompositePreprocessor& preprocessor) {
    blocks_.clear();
    partitions_.clear();
    haloPlan_.clear();
    std::vector<RawMeshBlock> rawBlocks;
    std::vector<RawMeshBlock> indexedSetBlocks;
    std::vector<MeshDecompose::SourcePatchExtent> sourceExtents;
    std::vector<SourceMetricData> sourceMetrics;

    broadcast("Mesh reading: ", std::to_string(meshFiles.size()) + " file(s)");
    for (const auto& mf : meshFiles) {
        std::string path = resolvePath(caseDir, mf);

        std::vector<RawMeshBlock> segments;
        if (!readSFMFile(path, segments)) {
            broadcast("Warning: ", "Failed to read " + path);
            continue;
        }

        for (size_t s = 0; s < segments.size(); ++s) {
            RawMeshBlock& raw = segments[s];
            bool hasPoints = !raw.x.empty();
            bool hasTopology = !raw.pointSets.empty() || !raw.topology.empty();

            if (hasPoints) {
                raw.name = blockNameFromPath(path, rawBlocks.size());
                if (segments.size() > 1) raw.name += "_" + std::to_string(s);
                raw.sourceFile = path;
                rawBlocks.push_back(std::move(raw));
                continue;
            }

            if (!hasTopology) continue;
            if (raw.nx > 0 && raw.ny > 0 && raw.nz > 0 && raw.ng > 0) {
                indexedSetBlocks.push_back(std::move(raw));
            } else if (!rawBlocks.empty()) {
                appendSetFile(raw, rawBlocks.back());
            } else {
                broadcast("Fatal: ", "Set-only mesh file appears before any #Point block: " + path);
                return false;
            }
        }
    }

    if (!indexedSetBlocks.empty()) {
        if (rawBlocks.empty()) {
            broadcast("Fatal: ", "Indexed topology sections were loaded before any mesh block.");
            return false;
        }
        if (indexedSetBlocks.size() % rawBlocks.size() != 0) {
            broadcast("Fatal: ",
                      "Indexed set sections must contain complete block groups.");
            return false;
        }
        for (size_t i = 0; i < indexedSetBlocks.size(); ++i) {
            appendSetFile(indexedSetBlocks[i], rawBlocks[i % rawBlocks.size()]);
        }
    }

    if (rawBlocks.empty()) {
        broadcast("Fatal: ", "No mesh block with #Point data was loaded.");
        return false;
    }

    if (!ensureBuiltinAllPointSets(rawBlocks)) return false;

    broadcast("Mesh source zones parsed: ", rawBlocks.size());
    for (const RawMeshBlock& raw : rawBlocks) {
        if (!validateRawBlock(raw)) return false;
    }

    MeshDecompose::CompositeMeshSummary topology;
    if (!MeshDecompose::assembleCompositeMesh(rawBlocks,
                                              config.haloTolerance,
                                              topology)) {
        broadcast("Fatal: ", "composite mesh topology assembly failed.");
        return false;
    }
    {
        std::ostringstream oss;
        oss << topology.sourceZoneCount << " source zone(s), "
            << topology.uniquePointCount << " unique physical points, "
            << topology.sharedPointGroupCount << " shared point groups, "
            << "max multiplicity=" << topology.maxPointMultiplicity;
        broadcast("Mesh composite topology: ", oss.str());
    }
    if (!reportCompositeSourceInterfaces(rawBlocks)) return false;
    if (!validateConfiguredSetCoverage(
            rawBlocks, config.configuredSetNames)) return false;
    if (config.parallelEnabled &&
        !buildSourceCanonicalMetrics(rawBlocks, config, sourceMetrics)) {
        broadcast("Fatal: ",
                  "canonical source metric construction failed.");
        return false;
    }
    if (preprocessor && !preprocessor(rawBlocks)) {
        broadcast("Fatal: ", "composite mesh preprocessor failed.");
        return false;
    }

    if (!config.parallelEnabled && rawBlocks.size() > 1) {
        broadcast("Fatal: ",
                  "composite multi-zone mesh currently requires MPI parallel solve; "
                  "serial composite solve is not implemented.");
        return false;
    }

    if (config.parallelEnabled) {
        if (config.automaticPartition) {
            broadcast("Fatal: ",
                      "MPI mesh decomposition requires explicit [parallel].split "
                      "for both single-zone and multi-zone input.");
            return false;
        }
        const int splitProduct =
            config.partitionSplit[0]
            * config.partitionSplit[1]
            * config.partitionSplit[2];
        if (splitProduct != config.parallelProcessCount) {
            broadcast("Fatal: ",
                      "parallel split product="
                      + std::to_string(splitProduct)
                      + " does not match configured nProcs="
                      + std::to_string(config.parallelProcessCount));
            return false;
        }

        std::vector<RawMeshBlock> patches;
        if (!MeshDecompose::decomposeCompositeMesh(rawBlocks,
                                                   config.partitionSplit,
                                                   config.haloTolerance,
                                                   patches,
                                                   sourceExtents,
                                                   partitions_)) {
            return false;
        }
        rawBlocks = std::move(patches);
        if (!ensureBuiltinAllPointSets(rawBlocks)) return false;

        std::ostringstream oss;
        oss << config.partitionSplit[0] << " x "
            << config.partitionSplit[1] << " x "
            << config.partitionSplit[2]
            << " -> " << partitions_.size()
            << " MPI partitions assembled from " << rawBlocks.size()
            << " internal structured storage patches";
        broadcast("Mesh decomposition: ", oss.str());
    } else {
        partitions_.resize(1);
        partitions_[0].rank = 0;
        partitions_[0].patchIds.push_back(0);
    }

    for (size_t blockId = 0; blockId < rawBlocks.size(); ++blockId) {
        RawMeshBlock& raw = rawBlocks[blockId];
        MeshBlockField block;
        block.name = raw.name;
        block.sourceFile = raw.sourceFile;
        block.sourceZoneId = raw.sourceZoneId;
        block.ownerRank = raw.ownerRank;
        block.globalPointIds = raw.globalPointIds;
        if (!raw.topology.empty()) {
            block.topology = raw.topology;
            for (MeshCell& cell : block.topology.cells) {
                cell.blockId = (int)blockId;
            }
            for (MeshFace& face : block.topology.faces) {
                face.ownerRank = raw.ownerRank;
                if (face.neighbourCell >= 0 && face.neighbourRank < 0) {
                    face.neighbourRank = raw.ownerRank;
                }
            }
        } else {
            block.topology = CellFaceMesh::fromStructuredPoints(
                raw.nx, raw.ny, raw.nz,
                raw.x, raw.y, raw.z,
                raw.pointSets,
                (int)blockId,
                raw.ownerRank);
        }
        if (!Mesh::setupFieldFromRaw(block.field,
                                     raw.nx, raw.ny, raw.nz, raw.ng,
                                     raw.x, raw.y, raw.z,
                                     raw.pointSets,
                                     config,
                                     raw.name)) {
            return false;
        }
        if (!raw.ibmPointData.empty()) {
            const size_t expected =
                (size_t)raw.nx * raw.ny * raw.nz;
            if (raw.ibmPointData.size() != expected) {
                broadcast("Fatal: ",
                          "IBM point-data size mismatch in " + raw.name);
                return false;
            }
            block.ibmPointData = std::move(raw.ibmPointData);
            size_t local = 0;
            const int ng = block.field.NG();
            for (int k = 0; k < raw.nz; ++k) {
                for (int j = 0; j < raw.ny; ++j) {
                    for (int i = 0; i < raw.nx; ++i, ++local) {
                        const RawIBMPointData& data =
                            block.ibmPointData[local];
                        const int fi = i + ng;
                        const int fj = j + ng;
                        const int fk = k + ng;
                        block.field.CellFlag(fi, fj, fk) = data.cellType;
                        block.field.setIBMFluidMask(
                            fi, fj, fk, data.cellType == FLUID_CELL);
                        block.field.setWallDistance(
                            fi, fj, fk, std::abs(data.signedDistance));
                    }
                }
            }
        }
        blocks_.push_back(std::move(block));
    }

    if (config.parallelEnabled &&
        !MeshDecompose::copyPartitionGhostCoordinates(blocks_, sourceExtents)) {
        broadcast("Fatal: ", "partition ghost-coordinate synchronization failed.");
        return false;
    }
    if (config.parallelEnabled) {
        // setupFieldFromRaw先按局部外推坐标计算了度规。切分后ghost坐标已被
        // 同源邻patch的真实坐标替换，必须刷新Jacobian/metric，避免新分核面
        // 继续使用切分前的外推几何。
        for (MeshBlockField& block : blocks_) {
            Mesh::refreshMetrics(block.field);
        }
        if (!applySourceCanonicalMetrics(
                blocks_, sourceExtents, sourceMetrics)) {
            broadcast("Fatal: ",
                      "failed to project source canonical metrics onto MPI patches.");
            return false;
        }
    }

    broadcast("Mesh structured patches ready: ", blocks_.size());

    if (config.parallelEnabled) {
        if (blocks_.size() < 2) {
            broadcast("haloExchange preprocess: ", "single patch, skipped.");
        } else {
            MeshCommunication::HaloPreprocessOptions options;
            options.tolerance = config.haloTolerance;
            options.requiredHaloWidth = config.requiredGhostLayers;
            options.allowNonMatchingInterfaces = false;
            options.useDeclaredHaloSets = false;
            options.physicalBoundaryNames =
                config.physicalBoundaryNames;
            if (!MeshCommunication::buildHaloExchangePlan(
                    blocks_, options, haloPlan_)) {
                broadcast("Fatal: ", "haloExchange preprocessing failed.");
                return false;
            }
            for (MeshBlockField& block : blocks_) {
                Mesh::auditMetricIdentity(block.field);
            }
            if (haloPlan_.mappedGhostCells == 0) {
                broadcast("Warning: ",
                          "haloExchange found no overlapping ghost points. "
                          "Check block overlap, declared halo sets, and tolerance.");
            }
        }
        haloPlan_.blockOwnerRanks.clear();
        haloPlan_.blockOwnerRanks.reserve(blocks_.size());
        for (const MeshBlockField& block : blocks_) {
            haloPlan_.blockOwnerRanks.push_back(block.ownerRank);
        }
    }

    broadcast("Mesh ready: ",
              std::to_string(partitions_.size()) + " partition(s), "
              + std::to_string(blocks_.size()) + " structured patch(es)");
    return true;
}

std::string MultiBlockMesh::resolvePath(const std::string& caseDir,
                                        const std::string& meshFile) {
    if (meshFile.empty()) return meshFile;
    if (meshFile[0] == '/' || meshFile[0] == '\\' || meshFile.rfind("..", 0) == 0) {
        return meshFile;
    }
    return caseDir + "/" + meshFile;
}

std::string MultiBlockMesh::blockNameFromPath(const std::string& path,
                                              size_t fallbackId) {
    size_t sep = path.find_last_of("/\\");
    std::string name = (sep == std::string::npos) ? path : path.substr(sep + 1);
    size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) name = name.substr(0, dot);
    if (name.empty()) name = "block" + std::to_string(fallbackId);
    return name;
}

bool MultiBlockMesh::readSFMFile(const std::string& path,
                                 std::vector<RawMeshBlock>& out) {
    std::ifstream file(path);
    if (!file.is_open()) return false;

    out.clear();

    RawMeshBlock current;
    current.sourceFile = path;
    std::string line, section;
    int activePatch = -1;
    int remainingPatchFaces = 0;
    std::string activePointSet;
    int remainingPointSetIds = 0;
    bool parseOk = true;
    std::string parseError;
    auto hasContent = [](const RawMeshBlock& raw) {
        return raw.nx > 0 || raw.ny > 0 || raw.nz > 0 || raw.ng > 0 ||
               !raw.x.empty() || !raw.pointSets.empty() ||
               !raw.topology.empty();
    };
    auto failParse = [&](const std::string& message) {
        parseOk = false;
        parseError = message;
        broadcast("Fatal: ", path + ": " + message);
    };
    auto flushCurrent = [&]() -> bool {
        if (remainingPatchFaces != 0) {
            failParse("patch face list ended early; remaining face count="
                      + std::to_string(remainingPatchFaces) + ".");
            return false;
        }
        if (remainingPointSetIds != 0) {
            failParse("pointSet '" + activePointSet
                      + "' ended early; remaining point count="
                      + std::to_string(remainingPointSetIds) + ".");
            return false;
        }
        if (!hasContent(current)) return true;
        normalizePointSets(current);
        if (!current.topology.empty()) {
            std::string error;
            if (!finalizeParsedTopology(current, (int)out.size(), &error)) {
                failParse(error);
                return false;
            }
        } else if (!current.x.empty()) {
            current.topology = CellFaceMesh::fromStructuredPoints(
                current.nx, current.ny, current.nz,
                current.x, current.y, current.z,
                current.pointSets,
                (int)out.size(),
                current.ownerRank);
        }
        out.push_back(std::move(current));
        current = RawMeshBlock{};
        current.sourceFile = path;
        activePatch = -1;
        remainingPatchFaces = 0;
        activePointSet.clear();
        remainingPointSetIds = 0;
        return true;
    };

    while (parseOk && std::getline(file, line)) {
        trimOnly(line);
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (remainingPatchFaces != 0) {
                failParse("patch '" + current.topology.patches[(std::size_t)activePatch].name
                          + "' declares more faces than provided before section "
                          + line + ".");
                break;
            }
            if (remainingPointSetIds != 0) {
                failParse("pointSet '" + activePointSet
                          + "' declares more points than provided before section "
                          + line + ".");
                break;
            }
            section = line.substr(1);
            trimOnly(section);
            if (isReservedSection(section, "INFORMATION")) {
                if (!flushCurrent()) break;
            } else if (isReservedSection(section, "INTERFACE")) {
                // Legacy interface sections are ignored. Parallel topology is
                // generated from [parallel].split and coordinate matching.
            } else if (isReservedSection(section, "CELL") ||
                       isReservedSection(section, "CELLS") ||
                       isReservedSection(section, "FACE") ||
                       isReservedSection(section, "FACES") ||
                       isReservedSection(section, "PATCH") ||
                       isReservedSection(section, "PATCHES") ||
                       isReservedSection(section, "POINTSET") ||
                       isReservedSection(section, "POINTSETS")) {
                activePatch = -1;
                remainingPatchFaces = 0;
                activePointSet.clear();
                remainingPointSetIds = 0;
            } else if (!isReservedSection(section, "POINT") && !section.empty()) {
                current.pointSets[section];
            }
            continue;
        }

        cleanLine(line);
        if (line.empty()) continue;

        std::stringstream ss(line);
        if (isReservedSection(section, "INFORMATION")) {
            ss >> current.nx >> current.ny >> current.nz >> current.ng;
        } else if (isReservedSection(section, "POINT")) {
            double tx, ty, tz;
            if (ss >> tx >> ty >> tz) {
                current.x.push_back(tx);
                current.y.push_back(ty);
                current.z.push_back(tz);
            }
        } else if (isReservedSection(section, "CELL") ||
                   isReservedSection(section, "CELLS")) {
            MeshCell cell;
            bool ok = true;
            for (int n = 0; n < 8; ++n) {
                if (!(ss >> cell.pointIds[(std::size_t)n])) {
                    ok = false;
                    break;
                }
            }
            if (!ok) {
                failParse("invalid #Cell row; expected eight point indices.");
                break;
            }
            current.topology.cells.push_back(cell);
        } else if (isReservedSection(section, "FACE") ||
                   isReservedSection(section, "FACES")) {
            MeshFace face;
            bool ok = true;
            for (int n = 0; n < 4; ++n) {
                if (!(ss >> face.pointIds[(std::size_t)n])) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                ok = static_cast<bool>(ss >> face.ownerCell
                                           >> face.neighbourCell
                                           >> face.patchId
                                           >> face.direction);
            }
            if (!ok) {
                failParse("invalid #Face row; expected four points plus owner, neighbour, patch, direction.");
                break;
            }
            if (!(ss >> face.ownerRank)) face.ownerRank = current.ownerRank;
            if (!(ss >> face.neighbourRank)) {
                face.neighbourRank =
                    (face.neighbourCell >= 0) ? current.ownerRank : -1;
            }
            current.topology.faces.push_back(face);
        } else if (isReservedSection(section, "PATCH") ||
                   isReservedSection(section, "PATCHES")) {
            if (remainingPatchFaces > 0) {
                int faceId;
                int consumed = 0;
                while (ss >> faceId) {
                    if (remainingPatchFaces <= 0) {
                        failParse("patch face list has more ids than declared.");
                        break;
                    }
                    current.topology.patches[(std::size_t)activePatch]
                        .faceIds.push_back(faceId);
                    --remainingPatchFaces;
                    ++consumed;
                }
                if (!parseOk) break;
                if (consumed == 0) {
                    failParse("patch face list row contains no face ids.");
                    break;
                }
                if (remainingPatchFaces == 0) activePatch = -1;
            } else {
                MeshPatch patch;
                int nFaces = -1;
                if (!static_cast<bool>(ss >> patch.name >> patch.type >> nFaces) ||
                    patch.name.empty() || nFaces < 0) {
                    failParse("invalid #Patch header; expected '<name> <type> <nFaces>'.");
                    break;
                }
                activePatch = (int)current.topology.patches.size();
                remainingPatchFaces = nFaces;
                current.topology.patches.push_back(std::move(patch));
                if (remainingPatchFaces == 0) activePatch = -1;
            }
        } else if (isReservedSection(section, "POINTSET") ||
                   isReservedSection(section, "POINTSETS")) {
            if (remainingPointSetIds > 0) {
                int pointId;
                int consumed = 0;
                while (ss >> pointId) {
                    if (remainingPointSetIds <= 0) {
                        failParse("pointSet list has more ids than declared.");
                        break;
                    }
                    current.pointSets[activePointSet].push_back(pointId);
                    --remainingPointSetIds;
                    ++consumed;
                }
                if (!parseOk) break;
                if (consumed == 0) {
                    failParse("pointSet row contains no point ids.");
                    break;
                }
                if (remainingPointSetIds == 0) activePointSet.clear();
            } else {
                int nPoints = -1;
                if (!(ss >> activePointSet >> nPoints) ||
                    activePointSet.empty() || nPoints < 0) {
                    failParse("invalid #PointSet header; expected '<name> <nPoints>'.");
                    break;
                }
                current.pointSets[activePointSet];
                remainingPointSetIds = nPoints;
                if (remainingPointSetIds == 0) activePointSet.clear();
            }
        } else if (isReservedSection(section, "INTERFACE")) {
            continue;
        } else if (!section.empty()) {
            int idx;
            while (ss >> idx) current.pointSets[section].push_back(idx);
        }
    }
    if (!parseOk) return false;
    if (!flushCurrent()) return false;

    return true;
}

bool MultiBlockMesh::appendSetFile(const RawMeshBlock& setFile,
                                   RawMeshBlock& target) {
    if (setFile.nx > 0 && target.nx > 0 &&
        (setFile.nx != target.nx ||
         setFile.ny != target.ny ||
         setFile.nz != target.nz ||
         setFile.ng != target.ng)) {
        broadcast("Warning: ", "Set block dimensions do not match mesh block " + target.name);
    }
    for (const auto& [name, indices] : setFile.pointSets) {
        appendPointSet(target.pointSets[name], indices);
    }
    return true;
}

bool MultiBlockMesh::validateRawBlock(const RawMeshBlock& raw) {
    if (raw.nx <= 0 || raw.ny <= 0 || raw.nz <= 0 || raw.ng <= 0) {
        broadcast("Fatal: ", "Invalid #Information in " + raw.sourceFile);
        return false;
    }

    int totalPts = raw.nx * raw.ny * raw.nz;
    if ((int)raw.x.size() != totalPts ||
        raw.x.size() != raw.y.size() ||
        raw.x.size() != raw.z.size()) {
        broadcast("Fatal: ", "Point count in " + raw.sourceFile + " is "
                  + std::to_string(raw.x.size()) + ", expected "
                  + std::to_string(totalPts) + ".");
        return false;
    }

    for (const auto& [name, indices] : raw.pointSets) {
        for (int idx : indices) {
            if (idx < 0 || idx >= totalPts) {
                std::ostringstream oss;
                oss << "Set '" << name << "' in " << raw.name
                    << " contains out-of-range point index " << idx
                    << ", valid range is [0, " << (totalPts - 1) << "].";
                broadcast("Fatal: ", oss.str());
                return false;
            }
        }
    }

    return true;
}

} // namespace SF
