/// @file SF_meshDecompose.cpp
/// @brief 结构网格 MPI 分解与 ownership 构造实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.13-----------*/

#include "SF_meshDecompose.h"

#include "SF_mesh.h"
#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace SF {
namespace MeshDecompose {
namespace {

struct BucketKey {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;

    bool operator==(const BucketKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct BucketKeyHash {
    size_t operator()(const BucketKey& key) const {
        size_t h = std::hash<std::int64_t>{}(key.x);
        h ^= std::hash<std::int64_t>{}(key.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<std::int64_t>{}(key.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct GlobalPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int references = 0;
    std::unordered_set<std::string> setNames;
};

int rawIndex(int i, int j, int k, int nx, int ny) {
    return (k * ny + j) * nx + i;
}

BucketKey pointBucket(double x, double y, double z, double tolerance) {
    return {
        (std::int64_t)std::floor(x / tolerance),
        (std::int64_t)std::floor(y / tolerance),
        (std::int64_t)std::floor(z / tolerance)
    };
}

void normalizeIndices(std::vector<int>& indices) {
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}

double rawCoordinate(const RawMeshBlock& raw,
                     int i, int j, int k,
                     int axis) {
    const int idx = rawIndex(i, j, k, raw.nx, raw.ny);
    if (axis == 0) return raw.x[(size_t)idx];
    if (axis == 1) return raw.y[(size_t)idx];
    return raw.z[(size_t)idx];
}

bool logicalCuts(int pointCount,int parts,
                 std::vector<int>& cuts,
                 const std::string& zoneName,int axis) {
    cuts.clear();
    if(pointCount<=0||parts<=0) return false;
    if(pointCount==1) {
        if(parts!=1) {
            broadcast("Fatal: ",
                "source-index decomposition cannot split degenerate axis "
                +std::to_string(axis)+" of zone "+zoneName+".");
            return false;
        }
        cuts={0,0};
        return true;
    }
    const int cells=pointCount-1;
    if(parts>cells) {
        broadcast("Fatal: ",
            "source-index split requests more partitions than cells on axis "
            +std::to_string(axis)+" of zone "+zoneName+".");
        return false;
    }
    cuts.reserve((size_t)parts+1);
    for(int part=0;part<=parts;++part) {
        cuts.push_back((cells*part)/parts);
    }
    for(size_t n=1;n<cuts.size();++n) {
        if(cuts[n]<=cuts[n-1]) {
            broadcast("Fatal: ",
                "source-index split produced an empty logical interval.");
            return false;
        }
    }
    return true;
}

const char* coordinateAxisName(int axis) {
    static const char* names[3] = {"x", "y", "z"};
    return axis >= 0 && axis < 3 ? names[axis] : "unknown";
}

bool localAxisIsStraightExtrusion(const RawMeshBlock& zone,
                                  int localAxis,
                                  int coordinateAxis,
                                  double tolerance) {
    const int dims[3] = {zone.nx, zone.ny, zone.nz};
    if (dims[localAxis] <= 1) return false;

    double direction = 0.0;
    for (int k = 0; k < zone.nz; ++k) {
        for (int j = 0; j < zone.ny; ++j) {
            for (int i = 0; i < zone.nx; ++i) {
                const int local[3] = {i, j, k};
                if (local[localAxis] + 1 >= dims[localAxis]) continue;

                int next[3] = {i, j, k};
                ++next[localAxis];
                const double delta[3] = {
                    rawCoordinate(zone, next[0], next[1], next[2], 0)
                        - rawCoordinate(zone, i, j, k, 0),
                    rawCoordinate(zone, next[0], next[1], next[2], 1)
                        - rawCoordinate(zone, i, j, k, 1),
                    rawCoordinate(zone, next[0], next[1], next[2], 2)
                        - rawCoordinate(zone, i, j, k, 2)
                };

                for (int axis = 0; axis < 3; ++axis) {
                    if (axis == coordinateAxis) continue;
                    if (std::abs(delta[axis]) > tolerance) return false;
                }
                if (std::abs(delta[coordinateAxis]) <= tolerance) return false;
                if (direction == 0.0) {
                    direction = delta[coordinateAxis];
                } else if (direction * delta[coordinateAxis] <= 0.0) {
                    return false;
                }
            }
        }
    }
    return direction != 0.0;
}

bool isGlobalExtrusionAxis(const std::vector<RawMeshBlock>& zones,
                           int coordinateAxis,
                           double tolerance) {
    for (const RawMeshBlock& zone : zones) {
        bool zoneMatches = false;
        for (int localAxis = 0; localAxis < 3; ++localAxis) {
            if (localAxisIsStraightExtrusion(zone, localAxis,
                                             coordinateAxis, tolerance)) {
                zoneMatches = true;
                break;
            }
        }
        if (!zoneMatches) return false;
    }
    return true;
}

int extrusionLocalAxis(const RawMeshBlock& zone,
                       int coordinateAxis,double tolerance) {
    for(int localAxis=0;localAxis<3;++localAxis) {
        if(localAxisIsStraightExtrusion(
                zone,localAxis,coordinateAxis,tolerance)) {
            return localAxis;
        }
    }
    return -1;
}

std::int64_t sourceCellWeight(const RawMeshBlock& zone) {
    return static_cast<std::int64_t>(std::max(zone.nx-1,1))
        *static_cast<std::int64_t>(std::max(zone.ny-1,1))
        *static_cast<std::int64_t>(std::max(zone.nz-1,1));
}

int maximumSourceZoneMultiplicity(const std::vector<RawMeshBlock>& zones) {
    int maximumGlobalId = -1;
    for (const RawMeshBlock& zone : zones) {
        for (int globalId : zone.globalPointIds) {
            maximumGlobalId = std::max(maximumGlobalId, globalId);
        }
    }
    if (maximumGlobalId < 0) return 0;

    std::vector<std::unordered_set<int>> sourceZones(
        (size_t)maximumGlobalId + 1);
    for (size_t zoneId = 0; zoneId < zones.size(); ++zoneId) {
        for (int globalId : zones[zoneId].globalPointIds) {
            if (globalId < 0 || globalId > maximumGlobalId) continue;
            sourceZones[(size_t)globalId].insert((int)zoneId);
        }
    }

    int maximumMultiplicity = 0;
    for (const auto& references : sourceZones) {
        maximumMultiplicity =
            std::max(maximumMultiplicity, (int)references.size());
    }
    return maximumMultiplicity;
}

bool validateCompositeSplitTopology(
    const std::vector<RawMeshBlock>& zones,
    const std::array<int, 3>& parts,
    double tolerance) {
    const int maximumMultiplicity = maximumSourceZoneMultiplicity(zones);
    if (maximumMultiplicity < 3) return true;

    std::array<bool, 3> extrusionAxes{false, false, false};
    for (int axis = 0; axis < 3; ++axis) {
        extrusionAxes[(size_t)axis] =
            isGlobalExtrusionAxis(zones, axis, tolerance);
    }

    std::vector<int> transverseSplitAxes;
    for (int axis = 0; axis < 3; ++axis) {
        if (parts[(size_t)axis] > 1 &&
            !extrusionAxes[(size_t)axis]) {
            transverseSplitAxes.push_back(axis);
        }
    }

    const bool multipleTransverseCuts = transverseSplitAxes.size() > 1;
    const bool nonBisectingTransverseCut =
        transverseSplitAxes.size() == 1 &&
        parts[(size_t)transverseSplitAxes.front()] != 2;
    if (!multipleTransverseCuts && !nonBisectingTransverseCut) return true;

    std::ostringstream oss;
    oss << "mesh topology does not support split=("
        << parts[0] << " " << parts[1] << " " << parts[2] << "). "
        << "The composite mesh contains a point shared by "
        << maximumMultiplicity
        << " source zones. Such topology may be bisected along only one "
        << "non-extrusion coordinate axis; globally straight extrusion axes "
        << "may be split independently. Detected extrusion axes:";
    bool wroteAxis = false;
    for (int axis = 0; axis < 3; ++axis) {
        if (!extrusionAxes[(size_t)axis]) continue;
        oss << (wroteAxis ? "," : " ") << coordinateAxisName(axis);
        wroteAxis = true;
    }
    if (!wroteAxis) oss << " none";
    broadcast("Fatal: ", oss.str());
    return false;
}

RawMeshBlock sliceZone(const RawMeshBlock& zone,
                       const SourcePatchExtent& extent,
                       int ownerRank) {
    RawMeshBlock patch;
    patch.name = zone.name + "_patch" + std::to_string(extent.patchId);
    patch.sourceFile = zone.sourceFile;
    patch.sourceZoneId = extent.sourceZoneId;
    patch.nx = extent.size[0];
    patch.ny = extent.size[1];
    patch.nz = extent.size[2];
    patch.ng = zone.ng;
    const size_t pointCount = (size_t)patch.nx * patch.ny * patch.nz;
    patch.x.reserve(pointCount);
    patch.y.reserve(pointCount);
    patch.z.reserve(pointCount);
    patch.globalPointIds.reserve(pointCount);
    if (!zone.ibmPointData.empty()) {
        patch.ibmPointData.reserve(pointCount);
    }

    for (int k = 0; k < patch.nz; ++k) {
        for (int j = 0; j < patch.ny; ++j) {
            for (int i = 0; i < patch.nx; ++i) {
                const int gi = extent.start[0] + i;
                const int gj = extent.start[1] + j;
                const int gk = extent.start[2] + k;
                const int gidx = rawIndex(gi, gj, gk, zone.nx, zone.ny);
                const double x = zone.x[(size_t)gidx];
                const double y = zone.y[(size_t)gidx];
                const double z = zone.z[(size_t)gidx];
                patch.x.push_back(x);
                patch.y.push_back(y);
                patch.z.push_back(z);
                patch.globalPointIds.push_back(zone.globalPointIds[(size_t)gidx]);
                if (!zone.ibmPointData.empty()) {
                    patch.ibmPointData.push_back(
                        zone.ibmPointData[(size_t)gidx]);
                }
            }
        }
    }

    for (const auto& [name, indices] : zone.pointSets) {
        auto& localSet = patch.pointSets[name];
        for (int gidx : indices) {
            if (gidx < 0 || gidx >= zone.nx * zone.ny * zone.nz) continue;
            const int gi = gidx % zone.nx;
            const int gj = (gidx / zone.nx) % zone.ny;
            const int gk = gidx / (zone.nx * zone.ny);
            const bool inside =
                gi >= extent.start[0] && gi < extent.start[0] + extent.size[0] &&
                gj >= extent.start[1] && gj < extent.start[1] + extent.size[1] &&
                gk >= extent.start[2] && gk < extent.start[2] + extent.size[2];
            if (!inside) continue;
            localSet.push_back(rawIndex(gi - extent.start[0],
                                        gj - extent.start[1],
                                        gk - extent.start[2],
                                        patch.nx, patch.ny));
        }
        normalizeIndices(localSet);
    }

    patch.ownerRank=ownerRank;
    return patch;
}

const SourcePatchExtent* findContainingSourceExtent(
    const std::vector<SourcePatchExtent>& extents,
    int sourceZoneId,
    const std::array<int, 3>& sourceIJK,
    int excludedPatch) {
    for (const SourcePatchExtent& extent : extents) {
        if (extent.patchId == excludedPatch) continue;
        if (extent.sourceZoneId != sourceZoneId) continue;
        const bool inside =
            sourceIJK[0] >= extent.start[0] &&
            sourceIJK[0] < extent.start[0] + extent.size[0] &&
            sourceIJK[1] >= extent.start[1] &&
            sourceIJK[1] < extent.start[1] + extent.size[1] &&
            sourceIJK[2] >= extent.start[2] &&
            sourceIJK[2] < extent.start[2] + extent.size[2];
        if (inside) return &extent;
    }
    return nullptr;
}

bool copyCrossSourceGhostCoordinates(std::vector<MeshBlockField>& patches,
                                     double tolerance) {
    (void)patches;
    (void)tolerance;
    // 跨源zone接口不能把邻块坐标线直接拷进本块ghost层；五块圆柱网格
    // 在中心块角点处并不存在唯一光滑结构坐标。保留setupFieldFromRaw给出的
    // 本块外推ghost坐标，使metrics沿本块坐标连续；匹配接口ghost状态由
    // haloExchange按GlobalPointId面拓扑做层向direct-copy。非匹配接口当前
    // fail-fast，不走隐式插值。
    return true;
}

bool validateSourceCellCoverage(
    const std::vector<RawMeshBlock>& zones,
    const std::vector<SourcePatchExtent>& extents) {
    for (size_t zoneId = 0; zoneId < zones.size(); ++zoneId) {
        const RawMeshBlock& zone = zones[zoneId];
        const int cellNx = std::max(zone.nx - 1, 1);
        const int cellNy = std::max(zone.ny - 1, 1);
        const int cellNz = std::max(zone.nz - 1, 1);
        std::vector<int> coverage(
            (size_t)cellNx * cellNy * cellNz, 0);

        auto cellIndex = [&](int i, int j, int k) {
            return ((size_t)k * cellNy + j) * cellNx + i;
        };

        for (const SourcePatchExtent& extent : extents) {
            if (extent.sourceZoneId != (int)zoneId) continue;
            if (extent.sourceSize !=
                std::array<int, 3>{zone.nx, zone.ny, zone.nz}) {
                broadcast("Fatal: ",
                          "source extent size mismatch for zone "
                          + zone.name);
                return false;
            }

            const int cellEnd[3] = {
                extent.cellStart[0] + extent.cellSize[0],
                extent.cellStart[1] + extent.cellSize[1],
                extent.cellStart[2] + extent.cellSize[2]
            };
            for (int axis = 0; axis < 3; ++axis) {
                const int maxCell =
                    axis == 0 ? cellNx : (axis == 1 ? cellNy : cellNz);
                if (extent.cellStart[(size_t)axis] < 0 ||
                    extent.cellSize[(size_t)axis] <= 0 ||
                    cellEnd[axis] > maxCell) {
                    std::ostringstream oss;
                    oss << "source-zone cell extent is invalid in "
                        << zone.name << ": patch=" << extent.patchId
                        << ", axis=" << axis
                        << ", cellStart="
                        << extent.cellStart[(size_t)axis]
                        << ", cellSize="
                        << extent.cellSize[(size_t)axis]
                        << ", maxCell=" << maxCell << ".";
                    broadcast("Fatal: ", oss.str());
                    return false;
                }
            }

            for (int k = extent.cellStart[2];
                 k < cellEnd[2];
                 ++k) {
                for (int j = extent.cellStart[1];
                     j < cellEnd[1];
                     ++j) {
                    for (int i = extent.cellStart[0];
                         i < cellEnd[0];
                         ++i) {
                        ++coverage[cellIndex(i, j, k)];
                    }
                }
            }
        }

        for (int k = 0; k < cellNz; ++k) {
            for (int j = 0; j < cellNy; ++j) {
                for (int i = 0; i < cellNx; ++i) {
                    const int count = coverage[cellIndex(i, j, k)];
                    if (count == 1) continue;
                    std::ostringstream oss;
                    oss << "source-zone cell partition is invalid in "
                        << zone.name << " at source cell ("
                        << i << "," << j << "," << k
                        << "): coverage=" << count
                        << ". Every source cell must belong to exactly one "
                           "MPI partition.";
                    broadcast("Fatal: ", oss.str());
                    return false;
                }
            }
        }
    }
    return true;
}

} // namespace

bool assembleCompositeMesh(std::vector<RawMeshBlock>& zones,
                           double tolerance,
                           CompositeMeshSummary& summary) {
    summary = {};
    summary.sourceZoneCount = (int)zones.size();
    if (zones.empty() || tolerance <= 0.0) return false;

    std::vector<GlobalPoint> globalPoints;
    globalPoints.reserve(1024);
    std::unordered_map<BucketKey, std::vector<int>, BucketKeyHash> buckets;

    const double toleranceSquared = tolerance * tolerance;
    for (size_t zoneId = 0; zoneId < zones.size(); ++zoneId) {
        RawMeshBlock& zone = zones[zoneId];
        zone.sourceZoneId = (int)zoneId;
        if (zone.x.size() != zone.y.size() || zone.x.size() != zone.z.size()) {
            return false;
        }
        zone.globalPointIds.assign(zone.x.size(), -1);

        for (size_t localId = 0; localId < zone.x.size(); ++localId) {
            const double x = zone.x[localId];
            const double y = zone.y[localId];
            const double z = zone.z[localId];
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
                broadcast("Fatal: ",
                          "non-finite mesh coordinate in source zone "
                          + std::to_string(zoneId));
                return false;
            }

            const BucketKey base = pointBucket(x, y, z, tolerance);
            int globalId = -1;
            for (int dk = -1; dk <= 1 && globalId < 0; ++dk) {
                for (int dj = -1; dj <= 1 && globalId < 0; ++dj) {
                    for (int di = -1; di <= 1 && globalId < 0; ++di) {
                        const BucketKey candidateKey{
                            base.x + di, base.y + dj, base.z + dk
                        };
                        auto bucket = buckets.find(candidateKey);
                        if (bucket == buckets.end()) continue;
                        for (int candidateId : bucket->second) {
                            const GlobalPoint& candidate =
                                globalPoints[(size_t)candidateId];
                            const double dx = candidate.x - x;
                            const double dy = candidate.y - y;
                            const double dz = candidate.z - z;
                            if (dx * dx + dy * dy + dz * dz <= toleranceSquared) {
                                globalId = candidateId;
                                break;
                            }
                        }
                    }
                }
            }

            if (globalId < 0) {
                globalId = (int)globalPoints.size();
                GlobalPoint point;
                point.x = x;
                point.y = y;
                point.z = z;
                globalPoints.push_back(std::move(point));
                buckets[base].push_back(globalId);
            }
            zone.globalPointIds[localId] = globalId;
            ++globalPoints[(size_t)globalId].references;
        }
    }

    for (const RawMeshBlock& zone : zones) {
        for (const auto& [name, indices] : zone.pointSets) {
            for (int localId : indices) {
                if (localId < 0 || localId >= (int)zone.globalPointIds.size()) continue;
                const int globalId = zone.globalPointIds[(size_t)localId];
                globalPoints[(size_t)globalId].setNames.insert(name);
            }
        }
    }

    for (RawMeshBlock& zone : zones) {
        for (size_t localId = 0; localId < zone.globalPointIds.size(); ++localId) {
            const int globalId = zone.globalPointIds[localId];
            for (const std::string& name :
                 globalPoints[(size_t)globalId].setNames) {
                zone.pointSets[name].push_back((int)localId);
            }
        }
        for (auto& [name, indices] : zone.pointSets) {
            (void)name;
            normalizeIndices(indices);
        }
    }

    summary.uniquePointCount = (int)globalPoints.size();
    for (const GlobalPoint& point : globalPoints) {
        summary.maxPointMultiplicity =
            std::max(summary.maxPointMultiplicity, point.references);
        if (point.references <= 1) continue;
        ++summary.sharedPointGroupCount;
        summary.sharedPointReferenceCount += point.references;
    }
    return true;
}

bool decomposeCompositeMesh(
    const std::vector<RawMeshBlock>& zones,
    const std::array<int, 3>& parts,
    double tolerance,
    std::vector<RawMeshBlock>& patches,
    std::vector<SourcePatchExtent>& extents,
    std::vector<MeshPartition>& partitions) {
    patches.clear();
    extents.clear();
    partitions.clear();
    if (zones.empty()) return false;

    const int partitionCount = parts[0] * parts[1] * parts[2];
    if (parts[0] <= 0 || parts[1] <= 0 || parts[2] <= 0 ||
        partitionCount <= 0) {
        return false;
    }
    if (!validateCompositeSplitTopology(zones, parts, tolerance)) {
        return false;
    }

    partitions.resize((size_t)partitionCount);
    for (int rank = 0; rank < partitionCount; ++rank) {
        partitions[(size_t)rank].rank = rank;
    }

    std::array<bool,3> extrusionAxes{false,false,false};
    if(zones.size()>1) {
        for(int axis=0;axis<3;++axis) {
            extrusionAxes[(size_t)axis]=
                isGlobalExtrusionAxis(zones,axis,tolerance);
        }
    }
    std::vector<int> transverseAxes;
    int transverseBins=1;
    for(int axis=0;axis<3;++axis) {
        if(parts[(size_t)axis]>1&&!extrusionAxes[(size_t)axis]
            &&zones.size()>1) {
            transverseAxes.push_back(axis);
            transverseBins*=parts[(size_t)axis];
        }
    }
    if(transverseBins>(int)zones.size()) {
        broadcast("Fatal: ",
            "topology-first decomposition has fewer source zones than "
            "requested transverse partitions; splitting through a source "
            "interface branch is not permitted.");
        return false;
    }

    std::vector<int> zoneBin(zones.size(),0);
    if(transverseBins>1) {
        std::vector<size_t> order(zones.size());
        for(size_t zone=0;zone<zones.size();++zone)order[zone]=zone;
        std::sort(order.begin(),order.end(),[&](size_t left,size_t right) {
            return sourceCellWeight(zones[left])
                >sourceCellWeight(zones[right]);
        });
        std::vector<std::int64_t> load((size_t)transverseBins,0);
        for(size_t zone:order) {
            const int bin=(int)std::distance(
                load.begin(),std::min_element(load.begin(),load.end()));
            zoneBin[zone]=bin;
            load[(size_t)bin]+=sourceCellWeight(zones[zone]);
        }
    }

    int nextPatchId = 0;
    for (size_t sourceZoneId = 0; sourceZoneId < zones.size(); ++sourceZoneId) {
        const RawMeshBlock& zone = zones[sourceZoneId];
        if (zone.globalPointIds.size() != zone.x.size()) {
            broadcast("Fatal: ",
                      "source zone has not been assembled into composite topology: "
                      + zone.name);
            return false;
        }

        std::array<std::vector<int>,3> cuts;
        const int dims[3]={zone.nx,zone.ny,zone.nz};
        std::array<int,3> localToRankAxis{-1,-1,-1};
        if(zones.size()==1) {
            for(int axis=0;axis<3;++axis) {
                if(!logicalCuts(
                        dims[axis],parts[(size_t)axis],
                        cuts[(size_t)axis],zone.name,axis)) {
                    return false;
                }
                localToRankAxis[(size_t)axis]=axis;
            }
        } else {
            for(int localAxis=0;localAxis<3;++localAxis) {
                cuts[(size_t)localAxis]=dims[localAxis]>1
                    ?std::vector<int>{0,dims[localAxis]-1}
                    :std::vector<int>{0,0};
            }
            for(int rankAxis=0;rankAxis<3;++rankAxis) {
                if(parts[(size_t)rankAxis]<=1
                    ||!extrusionAxes[(size_t)rankAxis]) continue;
                const int localAxis=extrusionLocalAxis(
                    zone,rankAxis,tolerance);
                if(localAxis<0
                    ||localToRankAxis[(size_t)localAxis]>=0) {
                    broadcast("Fatal: ",
                        "source topology cannot map requested extrusion "
                        "partition to a unique local index axis in zone "
                        +zone.name+".");
                    return false;
                }
                if(!logicalCuts(
                        dims[localAxis],parts[(size_t)rankAxis],
                        cuts[(size_t)localAxis],zone.name,localAxis)) {
                    return false;
                }
                localToRankAxis[(size_t)localAxis]=rankAxis;
            }
        }

        for (size_t kz = 0; kz + 1 < cuts[2].size(); ++kz) {
            for (size_t jy = 0; jy + 1 < cuts[1].size(); ++jy) {
                for (size_t ix = 0; ix + 1 < cuts[0].size(); ++ix) {
                    SourcePatchExtent extent;
                    extent.patchId = nextPatchId++;
                    extent.sourceZoneId = (int)sourceZoneId;
                    const size_t interval[3] = {ix, jy, kz};
                    for (int axis = 0; axis < 3; ++axis) {
                        const std::vector<int>& axisCuts = cuts[(size_t)axis];
                        const size_t p = interval[(size_t)axis];
                        const int lo = axisCuts[p];
                        const int hi = axisCuts[p + 1];

                        extent.cellStart[(size_t)axis] = lo;
                        extent.cellSize[(size_t)axis] =
                            dims[axis] > 1 ? hi - lo : 1;
                        extent.start[(size_t)axis] = lo;
                        extent.size[(size_t)axis] =
                            hi - extent.start[(size_t)axis] + 1;
                    }
                    extent.sourceSize = {zone.nx, zone.ny, zone.nz};
                    if (extent.size[0] <= 0 ||
                        extent.size[1] <= 0 ||
                        extent.size[2] <= 0) {
                        std::ostringstream oss;
                        oss << "Source-index decomposition produced an empty "
                               "owned-point patch in source zone "
                            << zone.name << " at patch "
                            << extent.patchId << ".";
                        broadcast("Fatal: ", oss.str());
                        return false;
                    }
                    if (extent.cellSize[0] <= 0 ||
                        extent.cellSize[1] <= 0 ||
                        extent.cellSize[2] <= 0) {
                        std::ostringstream oss;
                        oss << "Source-index decomposition produced an empty "
                               "owned-cell patch in source zone "
                            << zone.name << " at patch "
                            << extent.patchId << ".";
                        broadcast("Fatal: ", oss.str());
                        return false;
                    }

                    std::array<int,3> rankCoordinate{0,0,0};
                    int encodedBin=zoneBin[sourceZoneId];
                    for(int rankAxis:transverseAxes) {
                        rankCoordinate[(size_t)rankAxis]=
                            encodedBin%parts[(size_t)rankAxis];
                        encodedBin/=parts[(size_t)rankAxis];
                    }
                    const size_t localInterval[3]={ix,jy,kz};
                    for(int localAxis=0;localAxis<3;++localAxis) {
                        const int rankAxis=
                            localToRankAxis[(size_t)localAxis];
                        if(rankAxis>=0) {
                            rankCoordinate[(size_t)rankAxis]=
                                (int)localInterval[(size_t)localAxis];
                        }
                    }
                    const int ownerRank=
                        (rankCoordinate[2]*parts[1]+rankCoordinate[1])
                        *parts[0]+rankCoordinate[0];
                    RawMeshBlock patch =
                        sliceZone(zone,extent,ownerRank);
                    if (ownerRank < 0 || ownerRank >= partitionCount) {
                        return false;
                    }
                    partitions[(size_t)ownerRank].patchIds.push_back(extent.patchId);
                    patches.push_back(std::move(patch));
                    extents.push_back(extent);
                }
            }
        }
    }

    for (const MeshPartition& partition : partitions) {
        if (!partition.patchIds.empty()) continue;
        broadcast("Fatal: ",
                  "Source-index decomposition produced empty MPI partition "
                  + std::to_string(partition.rank));
        return false;
    }
    return validateSourceCellCoverage(zones, extents);
}

bool copyPartitionGhostCoordinates(
    std::vector<MeshBlockField>& patches,
    const std::vector<SourcePatchExtent>& extents) {
    if (patches.size() != extents.size()) return false;

    for (size_t patchId = 0; patchId < patches.size(); ++patchId) {
        Field& field = patches[patchId].field;
        const SourcePatchExtent& owner = extents[patchId];
        const int ng = field.NG();

        for (int k = 0; k < field.MZ(); ++k) {
            for (int j = 0; j < field.MY(); ++j) {
                for (int i = 0; i < field.MX(); ++i) {
                    const bool ghost =
                        i < ng || i >= ng + field.NX() ||
                        j < ng || j >= ng + field.NY() ||
                        k < ng || k >= ng + field.NZ();
                    if (!ghost) continue;

                    const std::array<int, 3> sourceIJK{
                        owner.start[0] + i - ng,
                        owner.start[1] + j - ng,
                        owner.start[2] + k - ng
                    };
                    if (sourceIJK[0] < 0 || sourceIJK[0] >= owner.sourceSize[0] ||
                        sourceIJK[1] < 0 || sourceIJK[1] >= owner.sourceSize[1] ||
                        sourceIJK[2] < 0 || sourceIJK[2] >= owner.sourceSize[2]) {
                        continue;
                    }

                    const SourcePatchExtent* donor =
                        findContainingSourceExtent(extents,
                                                   owner.sourceZoneId,
                                                   sourceIJK,
                                                   owner.patchId);
                    if (!donor) {
                        std::ostringstream oss;
                        oss << "internal partition ghost has no source-index "
                               "donor: patch="
                            << owner.patchId << ", sourceZone="
                            << owner.sourceZoneId << ", sourceIJK=("
                            << sourceIJK[0] << ","
                            << sourceIJK[1] << ","
                            << sourceIJK[2] << ").";
                        broadcast("Fatal: ", oss.str());
                        return false;
                    }

                    Field& donorField = patches[(size_t)donor->patchId].field;
                    const int di = sourceIJK[0] - donor->start[0] + donorField.NG();
                    const int dj = sourceIJK[1] - donor->start[1] + donorField.NG();
                    const int dk = sourceIJK[2] - donor->start[2] + donorField.NG();
                    field.X(i, j, k) = donorField.X(di, dj, dk);
                    field.Y(i, j, k) = donorField.Y(di, dj, dk);
                    field.Z(i, j, k) = donorField.Z(di, dj, dk);
                }
            }
        }
        Mesh::refreshMetrics(field);
    }
    if (!copyCrossSourceGhostCoordinates(patches, 1.0e-10)) {
        return false;
    }
    for (MeshBlockField& patch : patches) {
        Mesh::refreshMetrics(patch.field);
    }
    return true;
}

} // namespace MeshDecompose
} // namespace SF
