/// @file SF_compositeIBM.cpp
/// @brief 多 patch IBM 管理器的配置、几何分发与求解器接口适配。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.14-----------*/

#include "SF_compositeIBM.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace SF {
namespace IBM {
namespace {

using GeoProcessing::Point;
using GeoProcessing::SDFResult;

int rawIndex(int i, int j, int k, int nx, int ny) {
    return (k * ny + j) * nx + i;
}

RawIBMPointData pointDataFromSDF(const SDFResult& sdf) {
    RawIBMPointData data;
    data.signedDistance = sdf.signedDistance;
    if (!sdf.inside) {
        data.cellType = FLUID_CELL;
        return data;
    }

    constexpr double imageEpsilon = 1.0e-10;
    const Point image =
        sdf.closestPoint + sdf.normal * (sdf.distance + imageEpsilon);
    data.cellType = SOLID_CELL;
    data.hasGeometry = true;
    data.wallPoint = {
        sdf.closestPoint.x, sdf.closestPoint.y, sdf.closestPoint.z
    };
    data.imagePoint = {image.x, image.y, image.z};
    data.wallNormal = {sdf.normal.x, sdf.normal.y, sdf.normal.z};
    return data;
}

void writeFieldPointData(Field& field,
                         IBM::IBMGeometry& geometry,
                         int i, int j, int k,
                         const RawIBMPointData& data) {
    field.CellFlag(i, j, k) = data.cellType;
    field.setIBMFluidMask(i, j, k, data.cellType == FLUID_CELL);
    field.setWallDistance(i, j, k, std::abs(data.signedDistance));
    geometry.signedDistance(i, j, k) = data.signedDistance;
    geometry.ghostLayer(i, j, k) = data.ghostLayer;
    if (!data.hasGeometry) return;
    geometry.setGeometry(
        i, j, k,
        Vector3(data.wallPoint[0],
                data.wallPoint[1],
                data.wallPoint[2]),
        Vector3(data.imagePoint[0],
                data.imagePoint[1],
                data.imagePoint[2]),
        Vector3(data.wallNormal[0],
                data.wallNormal[1],
                data.wallNormal[2]),
        Vector3(data.wallVelocity[0],
                data.wallVelocity[1],
                data.wallVelocity[2]));
}

bool availablePoint(const Field& field, int i, int j, int k) {
    if (i < 0 || i >= field.MX() ||
        j < 0 || j >= field.MY() ||
        k < 0 || k >= field.MZ()) {
        return false;
    }
    const int ng = field.NG();
    const bool interior =
        i >= ng && i < ng + field.NX() &&
        j >= ng && j < ng + field.NY() &&
        k >= ng && k < ng + field.NZ();
    return interior || field.isCommunicationHalo(i, j, k);
}

} // namespace

bool CompositeIB::preprocessSourceZones(
    std::vector<RawMeshBlock>& zones,
    const std::string& caseDir,
    const std::vector<std::string>& stlFiles) {
    active_ = false;
    if (config_.method != FDM::IBMMethod::Ghost) {
        broadcast("Fatal: ",
                  "composite MPI IBM currently supports only type=ghost.");
        return false;
    }
    if (zones.empty() || stlFiles.empty()) {
        broadcast("Fatal: ",
                  "composite IBM requires source zones and at least one STL file.");
        return false;
    }
    if (!geometry_.loadSTLFiles(stlFiles, caseDir)) {
        broadcast("Fatal: ", "failed to load composite IBM STL geometry.");
        return false;
    }
    if (!geometry_.watertight()) {
        const auto& report = geometry_.watertightReport();
        broadcast("IBM warning: ",
                  "STL is not watertight. SDF signs may be unreliable.");
        broadcast("IBM boundary edges: ", report.boundaryEdges);
        broadcast("IBM non-manifold edges: ", report.nonManifoldEdges);
    }

    int maximumGlobalId = -1;
    for (const RawMeshBlock& zone : zones) {
        if (zone.globalPointIds.size() != zone.x.size()) {
            broadcast("Fatal: ",
                      "IBM source preprocessing requires assembled globalPointIds.");
            return false;
        }
        for (int globalId : zone.globalPointIds) {
            maximumGlobalId = std::max(maximumGlobalId, globalId);
        }
    }
    if (maximumGlobalId < 0) return false;

    const size_t pointCount = (size_t)maximumGlobalId + 1;
    std::vector<Point> points(pointCount);
    std::vector<unsigned char> seen(pointCount, 0);
    for (const RawMeshBlock& zone : zones) {
        for (size_t local = 0; local < zone.globalPointIds.size(); ++local) {
            const int globalId = zone.globalPointIds[local];
            if (globalId < 0 || (size_t)globalId >= pointCount) return false;
            if (seen[(size_t)globalId]) continue;
            points[(size_t)globalId] = {
                zone.x[local], zone.y[local], zone.z[local]
            };
            seen[(size_t)globalId] = 1;
        }
    }
    if (std::find(seen.begin(), seen.end(), 0) != seen.end()) {
        broadcast("Fatal: ", "composite IBM global-point table has holes.");
        return false;
    }

    broadcast("IBM composite preprocessing: ",
              std::to_string(pointCount) + " unique physical points");
    std::vector<RawIBMPointData> globalData(pointCount);
    std::vector<unsigned char> inside(pointCount, 0);
    for (size_t globalId = 0; globalId < pointCount; ++globalId) {
        const SDFResult sdf = geometry_.signedDistance(points[globalId]);
        if (!std::isfinite(sdf.signedDistance) ||
            !std::isfinite(sdf.distance)) {
            broadcast("Fatal: ",
                      "non-finite IBM signed distance at global point "
                      + std::to_string(globalId));
            return false;
        }
        globalData[globalId] = pointDataFromSDF(sdf);
        inside[globalId] = sdf.inside ? 1 : 0;
    }

    const int ghostDepth = config_.requiredGhostLayers;
    for (int layer = 1; layer <= ghostDepth; ++layer) {
        std::vector<unsigned char> mark(pointCount, 0);
        for (const RawMeshBlock& zone : zones) {
            for (int k = 0; k < zone.nz; ++k) {
                for (int j = 0; j < zone.ny; ++j) {
                    for (int i = 0; i < zone.nx; ++i) {
                        const int local =
                            rawIndex(i, j, k, zone.nx, zone.ny);
                        const int globalId =
                            zone.globalPointIds[(size_t)local];
                        if (!inside[(size_t)globalId] ||
                            globalData[(size_t)globalId].ghostLayer != 0) {
                            continue;
                        }

                        bool touchesTarget = false;
                        for (int dk = -1; dk <= 1 && !touchesTarget; ++dk) {
                            for (int dj = -1; dj <= 1 && !touchesTarget; ++dj) {
                                for (int di = -1; di <= 1; ++di) {
                                    if (di == 0 && dj == 0 && dk == 0) continue;
                                    const int ni = i + di;
                                    const int nj = j + dj;
                                    const int nk = k + dk;
                                    if (ni < 0 || ni >= zone.nx ||
                                        nj < 0 || nj >= zone.ny ||
                                        nk < 0 || nk >= zone.nz) {
                                        continue;
                                    }
                                    const int neighborLocal =
                                        rawIndex(ni, nj, nk,
                                                 zone.nx, zone.ny);
                                    const int neighborGlobal =
                                        zone.globalPointIds[
                                            (size_t)neighborLocal];
                                    if (neighborGlobal == globalId) continue;
                                    if (layer == 1) {
                                        touchesTarget =
                                            !inside[(size_t)neighborGlobal];
                                    } else {
                                        touchesTarget =
                                            globalData[(size_t)neighborGlobal]
                                                .ghostLayer == layer - 1;
                                    }
                                    if (touchesTarget) break;
                                }
                            }
                        }
                        if (touchesTarget) mark[(size_t)globalId] = 1;
                    }
                }
            }
        }
        for (size_t globalId = 0; globalId < pointCount; ++globalId) {
            if (mark[globalId]) {
                globalData[globalId].ghostLayer = layer;
            }
        }
    }

    int fluidCount = 0;
    int ghostCount = 0;
    int solidCount = 0;
    for (RawIBMPointData& data : globalData) {
        if (!data.hasGeometry) {
            data.cellType = FLUID_CELL;
            ++fluidCount;
        } else if (data.ghostLayer > 0) {
            data.cellType = IBM_GHOST_CELL;
            ++ghostCount;
        } else {
            data.cellType = SOLID_CELL;
            ++solidCount;
        }
    }

    for (RawMeshBlock& zone : zones) {
        zone.ibmPointData.resize(zone.globalPointIds.size());
        for (size_t local = 0; local < zone.globalPointIds.size(); ++local) {
            zone.ibmPointData[local] =
                globalData[(size_t)zone.globalPointIds[local]];
        }
    }

    broadcast("IBM composite fluid points: ", fluidCount);
    broadcast("IBM composite ghost points: ", ghostCount);
    broadcast("IBM composite solid points: ", solidCount);
    active_ = true;
    return true;
}

bool CompositeIB::buildBlockGeometries(
    MultiBlockMesh& mesh) {
    for (size_t bid = 0; bid < mesh.size(); ++bid) {
        MeshBlockField& block = mesh.block(bid);
        IBM::IBMGeometry& geometry = blockGeometries_[bid];
        geometry.setup(block.field.MX(),
                       block.field.MY(),
                       block.field.MZ());
        if (block.ibmPointData.empty()) continue;

        const size_t expected =
            (size_t)block.field.NX() * block.field.NY() * block.field.NZ();
        if (block.ibmPointData.size() != expected) {
            broadcast("Fatal: ",
                      "IBM block point-data size mismatch in " + block.name);
            return false;
        }

        const int ng = block.field.NG();
        size_t local = 0;
        for (int k = 0; k < block.field.NZ(); ++k) {
            for (int j = 0; j < block.field.NY(); ++j) {
                for (int i = 0; i < block.field.NX(); ++i, ++local) {
                    const RawIBMPointData& data =
                        block.ibmPointData[local];
                    const int fi = i + ng;
                    const int fj = j + ng;
                    const int fk = k + ng;
                    geometry.signedDistance(fi, fj, fk) =
                        data.signedDistance;
                    geometry.ghostLayer(fi, fj, fk) = data.ghostLayer;
                    if (!data.hasGeometry) continue;
                    geometry.setGeometry(
                        fi, fj, fk,
                        Vector3(data.wallPoint[0],
                                data.wallPoint[1],
                                data.wallPoint[2]),
                        Vector3(data.imagePoint[0],
                                data.imagePoint[1],
                                data.imagePoint[2]),
                        Vector3(data.wallNormal[0],
                                data.wallNormal[1],
                                data.wallNormal[2]),
                        Vector3(data.wallVelocity[0],
                                data.wallVelocity[1],
                                data.wallVelocity[2]));
                }
            }
        }
    }
    return true;
}

bool CompositeIB::classifyPartitionHalos(
    MultiBlockMesh& mesh) {
    if (!active_) return false;

    for (MeshBlockField& block : mesh.blocks()) {
        block.field.clearCommunicationHaloMask();
    }
    for (const MeshCommunication::HaloBlockPlan& blockPlan :
         mesh.haloExchangePlan().blockPlans) {
        if (blockPlan.blockId < 0 ||
            blockPlan.blockId >= (int)mesh.size()) {
            return false;
        }
        Field& field = mesh.block((size_t)blockPlan.blockId).field;
        for (const MeshCommunication::HaloCellMapping& mapping :
             blockPlan.cells) {
            field.setCommunicationHalo(mapping.ownerIJK[0],
                                       mapping.ownerIJK[1],
                                       mapping.ownerIJK[2], true);
        }
    }

    for (size_t bid = 0; bid < mesh.size(); ++bid) {
        Field& field = mesh.block(bid).field;
        for (int k = 0; k < field.MZ(); ++k) {
            for (int j = 0; j < field.MY(); ++j) {
                for (int i = 0; i < field.MX(); ++i) {
                    if (!field.isCommunicationHalo(i, j, k)) continue;
                    const Point point{
                        field.X(i, j, k),
                        field.Y(i, j, k),
                        field.Z(i, j, k)
                    };
                    const SDFResult sdf = geometry_.signedDistance(point);
                    if (!std::isfinite(sdf.signedDistance) ||
                        !std::isfinite(sdf.distance)) {
                        std::cerr
                            << "[SF FATAL] non-finite IBM halo signed distance at block "
                            << mesh.block(bid).name << " (" << i << "," << j << ","
                            << k << ")" << std::endl;
                        return false;
                    }
                    writeFieldPointData(field, blockGeometries_[bid],
                                        i, j, k,
                                        pointDataFromSDF(sdf));
                }
            }
        }
    }

    const int ghostDepth = config_.requiredGhostLayers;
    for (int layer = 1; layer <= ghostDepth; ++layer) {
        for (size_t bid = 0; bid < mesh.size(); ++bid) {
            Field& field = mesh.block(bid).field;
            IBM::IBMGeometry& geometry = blockGeometries_[bid];
            std::vector<int> mark;
            for (int k = 0; k < field.MZ(); ++k) {
                for (int j = 0; j < field.MY(); ++j) {
                    for (int i = 0; i < field.MX(); ++i) {
                        if (!field.isCommunicationHalo(i, j, k) ||
                            field.CellFlag(i, j, k) == FLUID_CELL ||
                            geometry.ghostLayer(i, j, k) != 0) {
                            continue;
                        }

                        bool touchesTarget = false;
                        for (int dk = -1; dk <= 1 && !touchesTarget; ++dk) {
                            for (int dj = -1; dj <= 1 && !touchesTarget; ++dj) {
                                for (int di = -1; di <= 1; ++di) {
                                    if (di == 0 && dj == 0 && dk == 0) continue;
                                    const int ni = i + di;
                                    const int nj = j + dj;
                                    const int nk = k + dk;
                                    if (!availablePoint(field, ni, nj, nk)) continue;
                                    if (layer == 1) {
                                        touchesTarget =
                                            field.CellFlag(ni, nj, nk)
                                            == FLUID_CELL;
                                    } else {
                                        touchesTarget =
                                            geometry.ghostLayer(ni, nj, nk)
                                            == layer - 1;
                                    }
                                    if (touchesTarget) break;
                                }
                            }
                        }
                        if (touchesTarget) {
                            mark.push_back(field.getIdx(i, j, k));
                        }
                    }
                }
            }
            for (int idx : mark) {
                int i = 0, j = 0, k = 0;
                field.getIJK(idx, i, j, k);
                geometry.ghostLayer(i, j, k) = layer;
                field.CellFlag(i, j, k) = IBM_GHOST_CELL;
            }
        }
    }
    return true;
}

bool CompositeIB::setupLocalPatches(
    MultiBlockMesh& mesh,
    const std::vector<int>& localPatchIds,
    int mpiRank) {
    if (!active_) return false;

    blockGeometries_.clear();
    blockGeometries_.resize(mesh.size());
    if (!buildBlockGeometries(mesh)) return false;
    if (!classifyPartitionHalos(mesh)) return false;

    patchWeights_.clear();
    patchWeights_.resize(mesh.size());
    patchReady_.assign(mesh.size(), 0);

    int localFluid = 0;
    int localGhost = 0;
    int localSolid = 0;
    for (int patchId : localPatchIds) {
        if (patchId < 0 || patchId >= (int)mesh.size()) return false;
        Field& field = mesh.block((size_t)patchId).field;
        if (!patchWeights_[(size_t)patchId].buildFromPreclassified(
                field, blockGeometries_[(size_t)patchId], config_)) {
            broadcast("Fatal: ",
                      "IBM local closure preprocessing failed for patch "
                      + std::to_string(patchId));
            return false;
        }
        if (config_.ilwEnabled) {
            const int requestedOrder = config_.requestedTaylorOrder();
            for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
                for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                    for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                        if (field.CellFlag(i, j, k) != IBM_GHOST_CELL) continue;
                        const auto& storage =
                            patchWeights_[(size_t)patchId].ilwStorage();
                        const size_t cell = (size_t)field.getIdx(i, j, k);
                        const bool hasPlan = storage.hasPlan(cell);
                        const IBMILWPointPlan* plan = hasPlan
                            ? &storage.plan(cell)
                            : nullptr;
                        if (plan && plan->maxOrder >= requestedOrder) continue;

                        std::cerr
                            << "[SF FATAL] MPI IBM ILW preprocessing incomplete"
                            << ": rank=" << mpiRank
                            << ", patch=" << patchId
                            << ", ijk=(" << i << "," << j << "," << k << ")"
                            << ", xyz=(" << field.X(i, j, k)
                            << "," << field.Y(i, j, k)
                            << "," << field.Z(i, j, k) << ")"
                            << ", requestedTaylorOrder=" << requestedOrder
                            << ", fluidSamples="
                            << storage.fluidCount(cell)
                            << ", normalSamples="
                            << storage.normalCount(cell)
                            << ", hasPlan=" << hasPlan
                            << ", planMaxOrder="
                            << (plan ? plan->maxOrder : 0)
                            << ", useT2=" << (plan ? plan->useT2 : false)
                            << ", higherFits="
                            << (plan ? plan->higherFits.size() : 0)
                            << std::endl;
                        return false;
                    }
                }
            }
        }
        patchReady_[(size_t)patchId] = 1;
        const auto& counts = patchWeights_[(size_t)patchId].counts();
        localFluid += counts.fluid;
        localGhost += counts.ghost;
        localSolid += counts.solid;
    }
    broadcast("IBM local patch fluid points: ", localFluid);
    broadcast("IBM local patch ghost points: ", localGhost);
    broadcast("IBM local patch solid points: ", localSolid);
    return true;
}

void CompositeIB::applyLocalPatches(
    MultiBlockMesh& mesh,
    const std::vector<int>& localPatchIds,
    double time,
    double dt) {
    if (!active_) return;
    if (!std::isfinite(time) || !std::isfinite(dt) || dt < 0.0) {
        throw std::invalid_argument(
            "Composite IBM apply requires finite time and non-negative dt.");
    }
    for (int patchId : localPatchIds) {
        if (patchId < 0 || patchId >= (int)mesh.size() ||
            patchId >= (int)patchReady_.size() ||
            !patchReady_[(size_t)patchId]) {
            std::cerr
                << "[SF FATAL] IBM apply requested for unprepared patch "
                << patchId << std::endl;
            std::exit(1);
        }
        ghostCell_.apply(mesh.block((size_t)patchId).field,
                         patchWeights_[(size_t)patchId],
                         config_);
    }
}

} // namespace IBM
} // namespace SF
