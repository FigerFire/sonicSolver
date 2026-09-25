/// @file SF_meshGen.cpp
/// @brief 结构网格 block、edge 与生成流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_meshGen.h"
#include "core/interfaces/SF_log.h"
#include "SF_edges.h"
#include "infrastructure/io/mesh/SF_generatedMeshWriter.h"
#include "SF_MultiBlockMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace SF {
namespace MeshGen {

static inline double N(int i, double t) { return i ? t : (1.0 - t); }

void generateBlockPoints(const double v[8][3],
                         const int    nc[3],
                         const double grade[3],
                         std::vector<double>& x,
                         std::vector<double>& y,
                         std::vector<double>& z) {

    int np[3] = { nc[0] + 1, nc[1] + 1, nc[2] + 1 };
    size_t total = (size_t)np[0] * np[1] * np[2];
    x.resize(total); y.resize(total); z.resize(total);

    auto paramCoord = [](int i, int n, double g) -> double {
        if (n <= 1) return 0.0;
        if (g == 1.0) return (double)i / (n - 1);
        double r = (g == 0.0) ? 1.0 : g;
        return (1.0 - std::pow(r, i)) / (1.0 - std::pow(r, n - 1));
    };

    const int map[2][2][2] = { {{0,1},{3,2}}, {{4,5},{7,6}} };

    size_t idx = 0;
    for (int k = 0; k < np[2]; ++k) {
        double zeta = paramCoord(k, np[2], grade[2]);
        for (int j = 0; j < np[1]; ++j) {
            double eta = paramCoord(j, np[1], grade[1]);
            for (int i = 0; i < np[0]; ++i) {
                double xi = paramCoord(i, np[0], grade[0]);
                double px = 0, py = 0, pz = 0;
                for (int kk = 0; kk < 2; ++kk)
                    for (int jj = 0; jj < 2; ++jj)
                        for (int ii = 0; ii < 2; ++ii) {
                            double w = N(ii, xi) * N(jj, eta) * N(kk, zeta);
                            int vi   = map[kk][jj][ii];
                            px += w * v[vi][0];
                            py += w * v[vi][1];
                            pz += w * v[vi][2];
                        }
                x[idx] = px; y[idx] = py; z[idx] = pz;
                ++idx;
            }
        }
    }
}

void generateStructuredMesh(const MeshParameters& params,
                            std::vector<double>& allX,
                            std::vector<double>& allY,
                            std::vector<double>& allZ,
                            std::vector<int>&    blockSizes) {
    // 有弧边定义 → TFI
    if (!params.edges.empty()) {
        generateStructuredMeshTFI(params, allX, allY, allZ, blockSizes);
        return;
    }
    // 纯直线 → trilinear
    allX.clear(); allY.clear(); allZ.clear();
    blockSizes.clear();

    for (auto& blk : params.blocks) {
        double v[8][3];
        for (int i = 0; i < 8; ++i) {
            auto& mv = params.vertices[blk.verts[i]];
            v[i][0] = mv.x * params.scale;
            v[i][1] = mv.y * params.scale;
            v[i][2] = mv.z * params.scale;
        }
        std::vector<double> bx, by, bz;
        generateBlockPoints(v, blk.cells, blk.grade, bx, by, bz);

        allX.insert(allX.end(), bx.begin(), bx.end());
        allY.insert(allY.end(), by.begin(), by.end());
        allZ.insert(allZ.end(), bz.begin(), bz.end());

        blockSizes.push_back(blk.cells[0] + 1);
        blockSizes.push_back(blk.cells[1] + 1);
        blockSizes.push_back(blk.cells[2] + 1);
    }
}

// ============================================================
//  网格输出
// ============================================================

namespace {

static size_t blockPointCount(const std::vector<int>& blockSizes, int b) {
    return (size_t)blockSizes[b*3] * blockSizes[b*3 + 1] * blockSizes[b*3 + 2];
}

static std::string numberedBlockName(const std::string& baseName, int blockId,
                                     const std::string& ext) {
    std::ostringstream name;
    name << baseName << "_block" << std::setw(3) << std::setfill('0') << blockId << ext;
    return name.str();
}

static std::array<int, 4> sortedFace(std::array<int, 4> f);

static double generatedMeshMatchTolerance(const std::vector<double>& allX,
                                          const std::vector<double>& allY,
                                          const std::vector<double>& allZ) {
    if (allX.empty()) return 1.0e-9;
    double minX = allX[0], maxX = allX[0];
    double minY = allY[0], maxY = allY[0];
    double minZ = allZ[0], maxZ = allZ[0];
    for (size_t i = 1; i < allX.size(); ++i) {
        minX = std::min(minX, allX[i]); maxX = std::max(maxX, allX[i]);
        minY = std::min(minY, allY[i]); maxY = std::max(maxY, allY[i]);
        minZ = std::min(minZ, allZ[i]); maxZ = std::max(maxZ, allZ[i]);
    }
    double dx = maxX - minX;
    double dy = maxY - minY;
    double dz = maxZ - minZ;
    double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
    return std::max(1.0e-9, 1.0e-12 * std::max(1.0, diag));
}

struct GeneratedBucketKey {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;

    bool operator==(const GeneratedBucketKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct GeneratedBucketKeyHash {
    size_t operator()(const GeneratedBucketKey& key) const {
        size_t h = std::hash<std::int64_t>{}(key.x);
        h ^= std::hash<std::int64_t>{}(key.y)
           + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<std::int64_t>{}(key.z)
           + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct GeneratedGlobalPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int references = 0;
};

struct GeneratedBoundaryFaceRef {
    int blockId = -1;
    int axis = -1;
    int side = -1;
    std::array<int, 4> globalPointIds{-1, -1, -1, -1};
};

struct GeneratedInterfaceSummary {
    int uniquePointCount = 0;
    int sharedPointGroups = 0;
    int branchPointGroups = 0;
    int maxPointMultiplicity = 1;
    int sharedFaceGroups = 0;
    int invalidFaceGroups = 0;
    int maxFaceMultiplicity = 1;
    std::map<std::string, int> pairFaceCounts;
};

static GeneratedBucketKey generatedPointBucket(double x,
                                               double y,
                                               double z,
                                               double tolerance) {
    return {
        (std::int64_t)std::floor(x / tolerance),
        (std::int64_t)std::floor(y / tolerance),
        (std::int64_t)std::floor(z / tolerance)
    };
}

static std::vector<std::vector<int>> assignGeneratedGlobalPointIds(
        const std::vector<double>& allX,
        const std::vector<double>& allY,
        const std::vector<double>& allZ,
        const std::vector<int>& blockSizes,
        double tolerance,
        std::vector<GeneratedGlobalPoint>& globalPoints) {
    std::vector<std::vector<int>> blockGlobalIds;
    const int nBlocks = (int)blockSizes.size() / 3;
    blockGlobalIds.resize((size_t)nBlocks);
    globalPoints.clear();
    globalPoints.reserve(allX.size());

    std::unordered_map<GeneratedBucketKey,
                       std::vector<int>,
                       GeneratedBucketKeyHash> buckets;
    const double toleranceSquared = tolerance * tolerance;

    size_t offset = 0;
    for (int b = 0; b < nBlocks; ++b) {
        const size_t nPts = blockPointCount(blockSizes, b);
        blockGlobalIds[(size_t)b].assign(nPts, -1);

        for (size_t localId = 0; localId < nPts; ++localId) {
            const size_t idx = offset + localId;
            const double x = allX[idx];
            const double y = allY[idx];
            const double z = allZ[idx];
            const GeneratedBucketKey base =
                generatedPointBucket(x, y, z, tolerance);

            int globalId = -1;
            for (int dk = -1; dk <= 1 && globalId < 0; ++dk) {
                for (int dj = -1; dj <= 1 && globalId < 0; ++dj) {
                    for (int di = -1; di <= 1 && globalId < 0; ++di) {
                        const GeneratedBucketKey key{
                            base.x + di, base.y + dj, base.z + dk
                        };
                        auto found = buckets.find(key);
                        if (found == buckets.end()) continue;
                        for (int candidateId : found->second) {
                            const GeneratedGlobalPoint& candidate =
                                globalPoints[(size_t)candidateId];
                            const double dx = candidate.x - x;
                            const double dy = candidate.y - y;
                            const double dz = candidate.z - z;
                            if (dx * dx + dy * dy + dz * dz
                                <= toleranceSquared) {
                                globalId = candidateId;
                                break;
                            }
                        }
                    }
                }
            }

            if (globalId < 0) {
                globalId = (int)globalPoints.size();
                GeneratedGlobalPoint point;
                point.x = x;
                point.y = y;
                point.z = z;
                globalPoints.push_back(point);
                buckets[base].push_back(globalId);
            }

            blockGlobalIds[(size_t)b][localId] = globalId;
            ++globalPoints[(size_t)globalId].references;
        }
        offset += nPts;
    }

    return blockGlobalIds;
}

static std::string generatedSideLabel(const GeneratedBoundaryFaceRef& ref) {
    const char axisName = ref.axis == 0 ? 'i' : (ref.axis == 1 ? 'j' : 'k');
    std::ostringstream oss;
    oss << "block" << ref.blockId << "." << axisName
        << (ref.side == 0 ? "Min" : "Max");
    return oss.str();
}

static std::string generatedPairLabel(const GeneratedBoundaryFaceRef& a,
                                      const GeneratedBoundaryFaceRef& b) {
    std::string lhs = generatedSideLabel(a);
    std::string rhs = generatedSideLabel(b);
    if (rhs < lhs) std::swap(lhs, rhs);
    return lhs + " <-> " + rhs;
}

static void appendGeneratedBoundaryFace(
        std::map<std::array<int, 4>,
                 std::vector<GeneratedBoundaryFaceRef>>& faceGroups,
        const std::vector<int>& globalIds,
        int blockId,
        int nx, int ny,
        int axis, int side,
        std::array<int, 4> localPointIds) {
    (void)nx;
    (void)ny;
    GeneratedBoundaryFaceRef ref;
    ref.blockId = blockId;
    ref.axis = axis;
    ref.side = side;
    for (int n = 0; n < 4; ++n) {
        const int pointId = localPointIds[(size_t)n];
        if (pointId < 0 || pointId >= (int)globalIds.size()) return;
        ref.globalPointIds[(size_t)n] = globalIds[(size_t)pointId];
    }
    auto key = ref.globalPointIds;
    std::sort(key.begin(), key.end());
    faceGroups[key].push_back(ref);
}

static void appendGeneratedBlockBoundaryFaces(
        std::map<std::array<int, 4>,
                 std::vector<GeneratedBoundaryFaceRef>>& faceGroups,
        const std::vector<int>& globalIds,
        int blockId,
        int nx, int ny, int nz) {
    auto idx = [nx, ny](int i, int j, int k) {
        return (k * ny + j) * nx + i;
    };
    if (nx < 2 || ny < 2 || nz < 2) return;

    for (int k = 0; k < nz - 1; ++k) {
        for (int j = 0; j < ny - 1; ++j) {
            appendGeneratedBoundaryFace(
                faceGroups, globalIds, blockId, nx, ny, 0, 0,
                {idx(0, j, k), idx(0, j + 1, k),
                 idx(0, j + 1, k + 1), idx(0, j, k + 1)});
            appendGeneratedBoundaryFace(
                faceGroups, globalIds, blockId, nx, ny, 0, 1,
                {idx(nx - 1, j, k), idx(nx - 1, j, k + 1),
                 idx(nx - 1, j + 1, k + 1), idx(nx - 1, j + 1, k)});
        }
    }

    for (int k = 0; k < nz - 1; ++k) {
        for (int i = 0; i < nx - 1; ++i) {
            appendGeneratedBoundaryFace(
                faceGroups, globalIds, blockId, nx, ny, 1, 0,
                {idx(i, 0, k), idx(i + 1, 0, k),
                 idx(i + 1, 0, k + 1), idx(i, 0, k + 1)});
            appendGeneratedBoundaryFace(
                faceGroups, globalIds, blockId, nx, ny, 1, 1,
                {idx(i, ny - 1, k), idx(i, ny - 1, k + 1),
                 idx(i + 1, ny - 1, k + 1), idx(i + 1, ny - 1, k)});
        }
    }

    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx - 1; ++i) {
            appendGeneratedBoundaryFace(
                faceGroups, globalIds, blockId, nx, ny, 2, 0,
                {idx(i, j, 0), idx(i + 1, j, 0),
                 idx(i + 1, j + 1, 0), idx(i, j + 1, 0)});
            appendGeneratedBoundaryFace(
                faceGroups, globalIds, blockId, nx, ny, 2, 1,
                {idx(i, j, nz - 1), idx(i, j + 1, nz - 1),
                 idx(i + 1, j + 1, nz - 1), idx(i + 1, j, nz - 1)});
        }
    }
}

static bool writeInterfaceTopologyReport(
        const std::string& filename,
        const std::vector<double>& allX,
        const std::vector<double>& allY,
        const std::vector<double>& allZ,
        const std::vector<int>& blockSizes,
        double tolerance,
        GeneratedInterfaceSummary& summary) {
    summary = {};
    const int nBlocks = (int)blockSizes.size() / 3;
    if (nBlocks <= 1) return true;
    if (allX.size() != allY.size() || allX.size() != allZ.size()) {
        broadcast("Fatal: ",
                  "createMesh interface report needs matching point arrays.");
        return false;
    }

    std::vector<GeneratedGlobalPoint> globalPoints;
    auto blockGlobalIds = assignGeneratedGlobalPointIds(
        allX, allY, allZ, blockSizes, tolerance, globalPoints);
    summary.uniquePointCount = (int)globalPoints.size();
    for (const GeneratedGlobalPoint& point : globalPoints) {
        summary.maxPointMultiplicity =
            std::max(summary.maxPointMultiplicity, point.references);
        if (point.references > 1) ++summary.sharedPointGroups;
        if (point.references >= 3) ++summary.branchPointGroups;
    }

    std::map<std::array<int, 4>,
             std::vector<GeneratedBoundaryFaceRef>> faceGroups;
    for (int b = 0; b < nBlocks; ++b) {
        const int nx = blockSizes[b * 3];
        const int ny = blockSizes[b * 3 + 1];
        const int nz = blockSizes[b * 3 + 2];
        appendGeneratedBlockBoundaryFaces(faceGroups,
                                          blockGlobalIds[(size_t)b],
                                          b, nx, ny, nz);
    }

    for (const auto& [key, refs] : faceGroups) {
        (void)key;
        if (refs.size() < 2) continue;
        std::set<int> blocks;
        for (const GeneratedBoundaryFaceRef& ref : refs) {
            blocks.insert(ref.blockId);
        }
        if (blocks.size() < 2) continue;

        summary.maxFaceMultiplicity =
            std::max(summary.maxFaceMultiplicity, (int)refs.size());
        if (refs.size() != 2 || blocks.size() != 2) {
            ++summary.invalidFaceGroups;
            continue;
        }

        ++summary.sharedFaceGroups;
        ++summary.pairFaceCounts[generatedPairLabel(refs[0], refs[1])];
    }

    std::ofstream out(filename);
    if (!out) {
        broadcast("Fatal: ", "Cannot open interface report " + filename);
        return false;
    }

    out << "# sonicSolver createMesh interface topology\n";
    out << "tolerance " << std::setprecision(15) << tolerance << "\n";
    out << "sourceBlocks " << nBlocks << "\n";
    out << "uniquePhysicalPoints " << summary.uniquePointCount << "\n";
    out << "sharedPointGroups " << summary.sharedPointGroups << "\n";
    out << "branchPointGroups " << summary.branchPointGroups << "\n";
    out << "maxPointMultiplicity " << summary.maxPointMultiplicity << "\n";
    out << "sharedFaceGroups " << summary.sharedFaceGroups << "\n";
    out << "maxFaceMultiplicity " << summary.maxFaceMultiplicity << "\n";
    out << "invalidFaceGroups " << summary.invalidFaceGroups << "\n\n";
    out << "# InterfacePairs\n";
    for (const auto& [label, count] : summary.pairFaceCounts) {
        out << label << " faces " << count << "\n";
    }

    if (summary.invalidFaceGroups > 0) {
        broadcast("Fatal: ",
                  "generated mesh contains boundary faces shared by more than "
                  "two source blocks; explicit structured interface topology "
                  "cannot be built.");
        return false;
    }
    if (summary.sharedPointGroups > 0 && summary.sharedFaceGroups == 0) {
        broadcast("Fatal: ",
                  "generated mesh has repeated physical points but no "
                  "conformal shared boundary faces.");
        return false;
    }
    return true;
}

static std::map<std::string, std::vector<int>> buildDefaultBlockSets(
    int nx, int ny, int nz) {
    auto idx = [nx, ny](int i, int j, int k) {
        return (k * ny + j) * nx + i;
    };
    auto normalize = [](std::vector<int>& ids) {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    };

    std::vector<int> all;
    all.reserve((size_t)nx * ny * nz);
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int id = idx(i, j, k);
                all.push_back(id);
            }
        }
    }

    std::map<std::string, std::vector<int>> sets;
    sets["all"] = std::move(all);

    auto appendIPlane = [&](const std::string& name, int i) {
        auto& ids = sets[name];
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                ids.push_back(idx(i, j, k));
        normalize(ids);
    };
    auto appendJPlane = [&](const std::string& name, int j) {
        auto& ids = sets[name];
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                ids.push_back(idx(i, j, k));
        normalize(ids);
    };
    auto appendKPlane = [&](const std::string& name, int k) {
        auto& ids = sets[name];
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                ids.push_back(idx(i, j, k));
        normalize(ids);
    };

    appendIPlane("Left", 0);
    appendIPlane("Right", nx - 1);
    appendIPlane("inlet", 0);
    appendIPlane("outlet", nx - 1);

    if (ny > 1) {
        appendJPlane("sideWalls", 0);
        appendJPlane("sideWalls", ny - 1);
        appendJPlane("walls", 0);
        appendJPlane("walls", ny - 1);
    }
    if (nz > 1) {
        appendKPlane("frontBack", 0);
        appendKPlane("frontBack", nz - 1);
        appendKPlane("walls", 0);
        appendKPlane("walls", nz - 1);
    }

    auto& leftHalf = sets["leftHalf"];
    auto& rightHalf = sets["rightHalf"];
    const int split = nx / 2;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                (i < split ? leftHalf : rightHalf).push_back(idx(i, j, k));
            }
        }
    }
    return sets;
}

static std::array<int, 4> sortedFace(std::array<int, 4> f) {
    std::sort(f.begin(), f.end());
    return f;
}

static void appendStructuredFacePoints(std::vector<int>& ids,
                                       int nx, int ny, int nz,
                                       int axis, bool upper) {
    auto idx = [nx, ny](int i, int j, int k) {
        return (k * ny + j) * nx + i;
    };

    if (axis == 0) {
        int i = upper ? nx - 1 : 0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                ids.push_back(idx(i, j, k));
    } else if (axis == 1) {
        int j = upper ? ny - 1 : 0;
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i)
                ids.push_back(idx(i, j, k));
    } else {
        int k = upper ? nz - 1 : 0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                ids.push_back(idx(i, j, k));
    }
}

static std::vector<std::map<std::string, std::vector<int>>> buildFacePatchSets(
    const MeshParameters& params,
    const std::vector<int>& blockSizes) {
    const int nBlocks = (int)blockSizes.size() / 3;
    std::vector<std::map<std::string, std::vector<int>>> sets((size_t)nBlocks);
    if (params.facePatches.empty()) return sets;

    for (auto& blockSets : sets) {
        for (const MeshFacePatch& patch : params.facePatches) {
            blockSets[patch.name];
        }
    }

    struct LocalFace {
        std::array<int, 4> verts;
        int axis = 0;
        bool upper = false;
    };

    for (int b = 0; b < nBlocks && b < (int)params.blocks.size(); ++b) {
        const MeshBlock& block = params.blocks[(size_t)b];
        const int* v = block.verts;
        const std::array<LocalFace, 6> localFaces = {{
            {{{v[0], v[3], v[7], v[4]}}, 0, false},
            {{{v[1], v[2], v[6], v[5]}}, 0, true},
            {{{v[0], v[1], v[5], v[4]}}, 1, false},
            {{{v[3], v[2], v[6], v[7]}}, 1, true},
            {{{v[0], v[1], v[2], v[3]}}, 2, false},
            {{{v[4], v[5], v[6], v[7]}}, 2, true}
        }};

        int nx = blockSizes[(size_t)b * 3];
        int ny = blockSizes[(size_t)b * 3 + 1];
        int nz = blockSizes[(size_t)b * 3 + 2];

        for (const MeshFacePatch& patch : params.facePatches) {
            for (const auto& patchFace : patch.faces) {
                const auto target = sortedFace(patchFace);
                for (const LocalFace& local : localFaces) {
                    if (sortedFace(local.verts) != target) continue;
                    appendStructuredFacePoints(sets[(size_t)b][patch.name],
                                               nx, ny, nz,
                                               local.axis, local.upper);
                }
            }
        }
    }

    for (auto& blockSets : sets) {
        for (auto& [name, ids] : blockSets) {
            std::sort(ids.begin(), ids.end());
            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        }
    }
    return sets;
}

static std::map<std::string, std::vector<int>> buildBlockSets(
    int nx, int ny, int nz,
    const std::map<std::string, std::vector<int>>& extraSets = {}) {
    std::map<std::string, std::vector<int>> sets = buildDefaultBlockSets(nx, ny, nz);
    for (const auto& [name, ids] : extraSets) {
        if (name == "all") continue;
        auto dst = ids;
        std::sort(dst.begin(), dst.end());
        dst.erase(std::unique(dst.begin(), dst.end()), dst.end());
        sets[name] = std::move(dst);
    }
    return sets;
}

static std::vector<std::map<std::string, std::vector<int>>> buildBlockSetSections(
    const std::vector<int>& blockSizes,
    const std::vector<std::map<std::string, std::vector<int>>>& extraSets = {}) {

    const int nBlocks = (int)blockSizes.size() / 3;
    std::vector<std::map<std::string, std::vector<int>>> sections;
    sections.reserve(nBlocks);
    for (int b = 0; b < nBlocks; ++b) {
        const int nx = blockSizes[b * 3];
        const int ny = blockSizes[b * 3 + 1];
        const int nz = blockSizes[b * 3 + 2];
        if ((size_t)b < extraSets.size()) {
            sections.push_back(buildBlockSets(nx, ny, nz, extraSets[(size_t)b]));
        } else {
            sections.push_back(buildBlockSets(nx, ny, nz));
        }
    }
    return sections;
}

} // namespace

bool writeCombinedMesh(const std::string& baseName,
                       const MeshParameters& params,
                       const std::vector<double>& allX,
                       const std::vector<double>& allY,
                       const std::vector<double>& allZ,
                       const std::vector<int>&    blockSizes,
                       int nGhost) {

    int nBlocks = (int)blockSizes.size() / 3;
    if (nBlocks == 0) {
        broadcast("Warning: ", "No mesh blocks to write.");
        return false;
    }

    auto explicitFaceSets = buildFacePatchSets(params, blockSizes);
    auto pointSetSections = buildBlockSetSections(blockSizes, explicitFaceSets);
    GeneratedInterfaceSummary interfaceSummary;
    const std::string interfaceReport = baseName + "_interfaces.txt";
    if (nBlocks > 1 &&
        !writeInterfaceTopologyReport(interfaceReport,
                                      allX, allY, allZ,
                                      blockSizes,
                                      generatedMeshMatchTolerance(allX,
                                                                  allY,
                                                                  allZ),
                                      interfaceSummary)) {
        return false;
    }
    GeneratedMeshWriter::writeMultiZoneSFM(baseName + ".sfm",
                                           allX, allY, allZ,
                                           blockSizes, nGhost,
                                           explicitFaceSets,
                                           pointSetSections);

    // 当前 SFM 已自包含 boundary/point sets。删除旧生成器遗留的 sidecar，
    // 避免同一 case 同时加载两套不同 ghost 容量的集合元数据。
    std::remove((baseName + "_sets.sfm").c_str());

    if (nBlocks == 1) {
        int nx = blockSizes[0];
        int ny = blockSizes[1];
        int nz = blockSizes[2];
        GeneratedMeshWriter::writeBlockVTK(baseName + ".vts", allX, allY, allZ, 0, nx, ny, nz);
        GeneratedMeshWriter::writeVTM(baseName + ".vtm", {baseName + ".vts"});
        broadcast("Generated mesh SFM: ", "wrote one self-contained structured zone.");
        return true;
    }

    std::vector<std::string> vtkFiles;
    vtkFiles.reserve(nBlocks);

    size_t offset = 0;
    for (int b = 0; b < nBlocks; ++b) {
        int nx = blockSizes[b*3];
        int ny = blockSizes[b*3 + 1];
        int nz = blockSizes[b*3 + 2];
        size_t nPts = blockPointCount(blockSizes, b);

        std::string vtkFile = numberedBlockName(baseName, b, ".vts");
        GeneratedMeshWriter::writeBlockVTK(vtkFile, allX, allY, allZ, offset, nx, ny, nz);
        vtkFiles.push_back(vtkFile);
        offset += nPts;
    }

    if (offset != allX.size() || allX.size() != allY.size() || allX.size() != allZ.size()) {
        broadcast("Warning: ", "Generated point arrays do not match block sizes.");
    }

    GeneratedMeshWriter::writeVTM(baseName + ".vtm", vtkFiles);
    std::remove((baseName + ".vts").c_str());
    std::remove((baseName + ".vtu").c_str());
    broadcast("Generated mesh SFM: ",
              "wrote " + std::to_string(nBlocks)
              + " structured zones in one self-contained file.");
    {
        std::ostringstream oss;
        oss << interfaceSummary.sharedFaceGroups
            << " shared boundary face groups, "
            << interfaceSummary.sharedPointGroups
            << " shared point groups, max point multiplicity="
            << interfaceSummary.maxPointMultiplicity
            << ", report=" << interfaceReport;
        broadcast("Generated mesh interfaces: ", oss.str());
    }
    return true;
}

// ============================================================
//  TFI 块生成 (弧边支持)
//  Hex 12 边索引: ξ-dir(0-1,3-2,4-5,7-6), η-dir(0-3,1-2,4-7,5-6), ζ-dir(0-4,1-5,3-7,2-6)
// ============================================================

/// 查找一条边 (v0↔v1) 在 edges 数组中的定义
static const MeshEdge* findEdge(const MeshEdge* edges, int nEdges, int a, int b) {
    for (int i = 0; i < nEdges; ++i)
        if ((edges[i].v0 == a && edges[i].v1 == b) || (edges[i].v0 == b && edges[i].v1 == a))
            return &edges[i];
    return nullptr;
}

/// 构建 block 的 12 条边表
static void buildEdgeTable(const MeshParameters& params, int blkIdx,
                            MeshEdge edgeTable[12]) {
    auto& blk = params.blocks[blkIdx];
    const int* v = blk.verts;
    int nEdges = (int)params.edges.size();

    // 12 条边的 (vi, vj) 对
    const int pairs[12][2] = {
        {v[0],v[1]}, {v[3],v[2]}, {v[4],v[5]}, {v[7],v[6]},   // ξ-dir
        {v[0],v[3]}, {v[1],v[2]}, {v[4],v[7]}, {v[5],v[6]},   // η-dir
        {v[0],v[4]}, {v[1],v[5]}, {v[3],v[7]}, {v[2],v[6]}    // ζ-dir
    };
    for (int i = 0; i < 12; ++i) {
        const MeshEdge* e = findEdge(params.edges.data(), nEdges, pairs[i][0], pairs[i][1]);
        if (e) {
            edgeTable[i] = *e;
            if (edgeTable[i].type == EdgeType::ARC) {
                edgeTable[i].arcP = params.scale * edgeTable[i].arcP;
            }
        } else {
            edgeTable[i] = {};  // default LINE
        }
    }
}

void generateTFIBlock(const double verts[8][3],
                      const int    nc[3],
                      const MeshEdge edges[12],
                      std::vector<double>& x,
                      std::vector<double>& y,
                      std::vector<double>& z) {

    int np[3] = { nc[0]+1, nc[1]+1, nc[2]+1 };
    size_t total = (size_t)np[0] * np[1] * np[2];
    x.resize(total); y.resize(total); z.resize(total);

    // 准备 12 条边两端的顶点坐标
    Vector3 V[8];
    for (int i = 0; i < 8; ++i) V[i] = {verts[i][0], verts[i][1], verts[i][2]};

    // TFI: 先沿边求点, 再混合
    size_t idx = 0;
    for (int k = 0; k < np[2]; ++k) {
        double zeta = (np[2] > 1) ? (double)k / (np[2] - 1) : 0.0;
        for (int j = 0; j < np[1]; ++j) {
            double eta = (np[1] > 1) ? (double)j / (np[1] - 1) : 0.0;
            for (int i = 0; i < np[0]; ++i) {
                double xi = (np[0] > 1) ? (double)i / (np[0] - 1) : 0.0;

                // ── ξ 方向边 (沿 i 方向) ──
                Vector3 e_xi[4] = {
                    edgeEval(V[0], V[1], edges[0], xi),   // j=0, k=0
                    edgeEval(V[3], V[2], edges[1], xi),   // j=1, k=0
                    edgeEval(V[4], V[5], edges[2], xi),   // j=0, k=1
                    edgeEval(V[7], V[6], edges[3], xi)    // j=1, k=1
                };
                // ── η 方向边 ──
                Vector3 e_eta[4] = {
                    edgeEval(V[0], V[3], edges[4], eta),   // i=0, k=0
                    edgeEval(V[1], V[2], edges[5], eta),   // i=1, k=0
                    edgeEval(V[4], V[7], edges[6], eta),   // i=0, k=1
                    edgeEval(V[5], V[6], edges[7], eta)    // i=1, k=1
                };
                // ── ζ 方向边 ──
                Vector3 e_zeta[4] = {
                    edgeEval(V[0], V[4], edges[8],  zeta), // i=0, j=0
                    edgeEval(V[1], V[5], edges[9],  zeta), // i=1, j=0
                    edgeEval(V[3], V[7], edges[10], zeta), // i=0, j=1
                    edgeEval(V[2], V[6], edges[11], zeta)  // i=1, j=1
                };

                // ── 3D TFI: P = P_ξ + P_η + P_ζ - P_ξη - P_ηζ - P_ξζ + T ──
                double px = 0, py = 0, pz = 0;

                // P_ξ: 单向插值沿 ξ 边
                px += (1-eta)*(1-zeta)*e_xi[0].x + eta*(1-zeta)*e_xi[1].x + (1-eta)*zeta*e_xi[2].x + eta*zeta*e_xi[3].x;
                py += (1-eta)*(1-zeta)*e_xi[0].y + eta*(1-zeta)*e_xi[1].y + (1-eta)*zeta*e_xi[2].y + eta*zeta*e_xi[3].y;
                pz += (1-eta)*(1-zeta)*e_xi[0].z + eta*(1-zeta)*e_xi[1].z + (1-eta)*zeta*e_xi[2].z + eta*zeta*e_xi[3].z;

                // P_η: 单向插值沿 η 边
                px += (1-xi)*(1-zeta)*e_eta[0].x + xi*(1-zeta)*e_eta[1].x + (1-xi)*zeta*e_eta[2].x + xi*zeta*e_eta[3].x;
                py += (1-xi)*(1-zeta)*e_eta[0].y + xi*(1-zeta)*e_eta[1].y + (1-xi)*zeta*e_eta[2].y + xi*zeta*e_eta[3].y;
                pz += (1-xi)*(1-zeta)*e_eta[0].z + xi*(1-zeta)*e_eta[1].z + (1-xi)*zeta*e_eta[2].z + xi*zeta*e_eta[3].z;

                // P_ζ: 单向插值沿 ζ 边
                px += (1-xi)*(1-eta)*e_zeta[0].x + xi*(1-eta)*e_zeta[1].x + (1-xi)*eta*e_zeta[2].x + xi*eta*e_zeta[3].x;
                py += (1-xi)*(1-eta)*e_zeta[0].y + xi*(1-eta)*e_zeta[1].y + (1-xi)*eta*e_zeta[2].y + xi*eta*e_zeta[3].y;
                pz += (1-xi)*(1-eta)*e_zeta[0].z + xi*(1-eta)*e_zeta[1].z + (1-xi)*eta*e_zeta[2].z + xi*eta*e_zeta[3].z;

                // - P_ξη: 双线性面 (k=0/k=1 面沿 ξ-η), 在 ζ 中混合
                {
                    double bx_k0 = (1-xi)*(1-eta)*V[0].x + xi*(1-eta)*V[1].x + (1-xi)*eta*V[3].x + xi*eta*V[2].x;
                    double bx_k1 = (1-xi)*(1-eta)*V[4].x + xi*(1-eta)*V[5].x + (1-xi)*eta*V[7].x + xi*eta*V[6].x;
                    double by_k0 = (1-xi)*(1-eta)*V[0].y + xi*(1-eta)*V[1].y + (1-xi)*eta*V[3].y + xi*eta*V[2].y;
                    double by_k1 = (1-xi)*(1-eta)*V[4].y + xi*(1-eta)*V[5].y + (1-xi)*eta*V[7].y + xi*eta*V[6].y;
                    double bz_k0 = (1-xi)*(1-eta)*V[0].z + xi*(1-eta)*V[1].z + (1-xi)*eta*V[3].z + xi*eta*V[2].z;
                    double bz_k1 = (1-xi)*(1-eta)*V[4].z + xi*(1-eta)*V[5].z + (1-xi)*eta*V[7].z + xi*eta*V[6].z;
                    px -= (1-zeta)*bx_k0 + zeta*bx_k1;
                    py -= (1-zeta)*by_k0 + zeta*by_k1;
                    pz -= (1-zeta)*bz_k0 + zeta*bz_k1;
                }

                // - P_ηζ: 双线性面 (i=0/i=1 面沿 η-ζ), 在 ξ 中混合
                {
                    double bx_i0 = (1-eta)*(1-zeta)*V[0].x + eta*(1-zeta)*V[3].x + (1-eta)*zeta*V[4].x + eta*zeta*V[7].x;
                    double bx_i1 = (1-eta)*(1-zeta)*V[1].x + eta*(1-zeta)*V[2].x + (1-eta)*zeta*V[5].x + eta*zeta*V[6].x;
                    double by_i0 = (1-eta)*(1-zeta)*V[0].y + eta*(1-zeta)*V[3].y + (1-eta)*zeta*V[4].y + eta*zeta*V[7].y;
                    double by_i1 = (1-eta)*(1-zeta)*V[1].y + eta*(1-zeta)*V[2].y + (1-eta)*zeta*V[5].y + eta*zeta*V[6].y;
                    double bz_i0 = (1-eta)*(1-zeta)*V[0].z + eta*(1-zeta)*V[3].z + (1-eta)*zeta*V[4].z + eta*zeta*V[7].z;
                    double bz_i1 = (1-eta)*(1-zeta)*V[1].z + eta*(1-zeta)*V[2].z + (1-eta)*zeta*V[5].z + eta*zeta*V[6].z;
                    px -= (1-xi)*bx_i0 + xi*bx_i1;
                    py -= (1-xi)*by_i0 + xi*by_i1;
                    pz -= (1-xi)*bz_i0 + xi*bz_i1;
                }

                // - P_ξζ: 双线性面 (j=0/j=1 面沿 ξ-ζ), 在 η 中混合
                {
                    double bx_j0 = (1-xi)*(1-zeta)*V[0].x + xi*(1-zeta)*V[1].x + (1-xi)*zeta*V[4].x + xi*zeta*V[5].x;
                    double bx_j1 = (1-xi)*(1-zeta)*V[3].x + xi*(1-zeta)*V[2].x + (1-xi)*zeta*V[7].x + xi*zeta*V[6].x;
                    double by_j0 = (1-xi)*(1-zeta)*V[0].y + xi*(1-zeta)*V[1].y + (1-xi)*zeta*V[4].y + xi*zeta*V[5].y;
                    double by_j1 = (1-xi)*(1-zeta)*V[3].y + xi*(1-zeta)*V[2].y + (1-xi)*zeta*V[7].y + xi*zeta*V[6].y;
                    double bz_j0 = (1-xi)*(1-zeta)*V[0].z + xi*(1-zeta)*V[1].z + (1-xi)*zeta*V[4].z + xi*zeta*V[5].z;
                    double bz_j1 = (1-xi)*(1-zeta)*V[3].z + xi*(1-zeta)*V[2].z + (1-xi)*zeta*V[7].z + xi*zeta*V[6].z;
                    px -= (1-eta)*bx_j0 + eta*bx_j1;
                    py -= (1-eta)*by_j0 + eta*by_j1;
                    pz -= (1-eta)*bz_j0 + eta*bz_j1;
                }

                // + T (三线性顶点项)
                px += (1-xi)*(1-eta)*(1-zeta)*V[0].x + xi*(1-eta)*(1-zeta)*V[1].x
                    + (1-xi)*eta*(1-zeta)*V[3].x   + xi*eta*(1-zeta)*V[2].x
                    + (1-xi)*(1-eta)*zeta*V[4].x   + xi*(1-eta)*zeta*V[5].x
                    + (1-xi)*eta*zeta*V[7].x       + xi*eta*zeta*V[6].x;
                py += (1-xi)*(1-eta)*(1-zeta)*V[0].y + xi*(1-eta)*(1-zeta)*V[1].y
                    + (1-xi)*eta*(1-zeta)*V[3].y   + xi*eta*(1-zeta)*V[2].y
                    + (1-xi)*(1-eta)*zeta*V[4].y   + xi*(1-eta)*zeta*V[5].y
                    + (1-xi)*eta*zeta*V[7].y       + xi*eta*zeta*V[6].y;
                pz += (1-xi)*(1-eta)*(1-zeta)*V[0].z + xi*(1-eta)*(1-zeta)*V[1].z
                    + (1-xi)*eta*(1-zeta)*V[3].z   + xi*eta*(1-zeta)*V[2].z
                    + (1-xi)*(1-eta)*zeta*V[4].z   + xi*(1-eta)*zeta*V[5].z
                    + (1-xi)*eta*zeta*V[7].z       + xi*eta*zeta*V[6].z;

                x[idx] = px; y[idx] = py; z[idx] = pz;
                ++idx;
            }
        }
    }
}

void generateStructuredMeshTFI(const MeshParameters& params,
                               std::vector<double>& allX,
                               std::vector<double>& allY,
                               std::vector<double>& allZ,
                               std::vector<int>&    blockSizes) {
    allX.clear(); allY.clear(); allZ.clear();
    blockSizes.clear();

    for (size_t b = 0; b < params.blocks.size(); ++b) {
        auto& blk = params.blocks[b];

        double v[8][3];
        for (int i = 0; i < 8; ++i) {
            auto& mv = params.vertices[blk.verts[i]];
            v[i][0] = mv.x * params.scale;
            v[i][1] = mv.y * params.scale;
            v[i][2] = mv.z * params.scale;
        }

        MeshEdge edgeTable[12];
        buildEdgeTable(params, (int)b, edgeTable);

        std::vector<double> bx, by, bz;
        generateTFIBlock(v, blk.cells, edgeTable, bx, by, bz);

        allX.insert(allX.end(), bx.begin(), bx.end());
        allY.insert(allY.end(), by.begin(), by.end());
        allZ.insert(allZ.end(), bz.begin(), bz.end());

        blockSizes.push_back(blk.cells[0] + 1);
        blockSizes.push_back(blk.cells[1] + 1);
        blockSizes.push_back(blk.cells[2] + 1);
    }
}

} // namespace MeshGen
} // namespace SF
