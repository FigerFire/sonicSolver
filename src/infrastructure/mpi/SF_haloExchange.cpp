/// @file SF_haloExchange.cpp
/// @brief halo、GlobalDof 汇总与并行协调基础设施实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_haloExchange.h"
#include "SF_globalDofResidual.h"
#include "SF_MultiBlockMesh.h"
#include "SF_canonicalFace.h"
#include "SF_lxF.h"
#include "methods/numerics/structured/SF_structured.h"
#include "core/interfaces/SF_log.h"
#include "SF_valueTypes.h"
#include "SF_utility.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>

namespace SF {
namespace Parallel {

using namespace MeshCommunication;

void HaloExchange::assembleGlobalDofResiduals(
        std::vector<MeshBlockField>& blocks,
        const std::vector<Residual*>& residuals) const {
    if (!plan_) {
        throw std::runtime_error(
            "GlobalDof residual assembly requires a communication plan.");
    }
    GlobalDofResidualAssembler::assemble(*plan_, blocks, residuals, backend_);
}

void HaloExchange::assembleGlobalDofResiduals(
        Field& localField, Residual& residual) const {
    if (!plan_ || localBlockId_ < 0) {
        // 单块串行网格没有共享 GlobalDof，也不会构造通信计划。此时
        // AccumulateGlobalDof 的数学作用是恒等映射，空间残差继续读取
        // 已暂存的 local residual。MPI 激活后缺失计划仍必须立即报错。
        if (!mpiEnabled_) return;
        throw std::runtime_error(
            "GlobalDof residual assembly requires a configured local block.");
    }
    GlobalDofResidualAssembler::assemble(
        *plan_, localField, residual, localBlockId_, backend_);
}

void HaloExchange::configure(
        const HaloExchangePlan* plan,
        int localBlockId,
        const Backend::CommunicationBackend* backend) {
    plan_ = plan;
    localBlockId_ = localBlockId;
    backend_ = backend;
    rank_ = backend ? backend->rank() : 0;
    size_ = backend ? std::max(1, backend->size()) : 1;
    mpiRank_ = rank_;
    mpiSize_ = size_;
    mpiEnabled_ = backend && backend->active();
}

int HaloExchange::ownerRank(int blockId) const {
    if (blockId < 0 || size_ <= 0) return -1;
    if (plan_ && blockId < (int)plan_->blockOwnerRanks.size()) {
        return plan_->blockOwnerRanks[(size_t)blockId];
    }
    return blockId % size_;
}

const HaloBlockPlan* HaloExchange::localPlan() const {
    return blockPlan(localBlockId_);
}

const HaloBlockPlan* HaloExchange::blockPlan(int blockId) const {
    if (!plan_) return nullptr;
    for (const auto& blockPlan : plan_->blockPlans) {
        if (blockPlan.blockId == blockId) return &blockPlan;
    }
    return nullptr;
}

std::vector<int> HaloExchange::communicationRanks() const {
    std::set<int> ranks;
    if (!plan_) return {};
    auto connect = [&](int firstBlock, int secondBlock) {
        const int first = ownerRank(firstBlock);
        const int second = ownerRank(secondBlock);
        if (first == rank_ && second >= 0 && second != rank_) {
            ranks.insert(second);
        }
        if (second == rank_ && first >= 0 && first != rank_) {
            ranks.insert(first);
        }
    };
    for (const auto& block : plan_->blockPlans) {
        for (const auto& mapping : block.cells) {
            connect(mapping.ownerBlock, mapping.donorBlock);
        }
    }
    for (const auto& group : plan_->interfaceSyncGroups) {
        for (size_t first = 0; first < group.points.size(); ++first) {
            for (size_t second = first + 1;
                 second < group.points.size(); ++second) {
                connect(group.points[first].blockId,
                        group.points[second].blockId);
            }
        }
    }
    for (const auto& group : plan_->interfaceFluxSyncGroups) {
        for (size_t first = 0; first < group.faces.size(); ++first) {
            for (size_t second = first + 1;
                 second < group.faces.size(); ++second) {
                connect(group.faces[first].blockId,
                        group.faces[second].blockId);
            }
        }
    }
    return {ranks.begin(), ranks.end()};
}

std::vector<int> HaloExchange::rankInteriorCounts(
        int components) const {
    if (!plan_ || components <= 0) {
        throw std::runtime_error(
            "Halo communication requires a positive component count and "
            "a precomputed topology plan.");
    }
    std::vector<int> counts((size_t)size_, 0);
    for (const auto& block : plan_->blockInfos) {
        const int owner = ownerRank(block.blockId);
        if (owner < 0 || owner >= size_
            || block.interiorPointCount < 0) {
            throw std::runtime_error(
                "Halo plan contains an invalid block owner/count.");
        }
        counts[(size_t)owner] +=
            block.interiorPointCount * components;
    }
    return counts;
}

void HaloExchange::exchangeNeighbourPayload(
        const std::vector<double>& send,
        const std::vector<int>& counts,
        const std::vector<int>& displacements,
        int tag,
        std::vector<double>& received) const {
    if (!backend_ || !backend_->active()) return;
    if ((int)counts.size() != size_
        || (int)displacements.size() != size_
        || counts[(size_t)rank_] != (int)send.size()) {
        throw std::runtime_error(
            "Neighbour payload layout differs from the local static plan.");
    }
    std::copy(send.begin(), send.end(),
              received.begin() + displacements[(size_t)rank_]);
    backend_->exchangeNeighbours(
        send, counts, displacements, communicationRanks(), tag, received);
}

namespace {

struct InterfacePoint3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

static InterfacePoint3 operator+(const InterfacePoint3& a,
                                 const InterfacePoint3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

static InterfacePoint3 operator*(const InterfacePoint3& a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}

static double interfaceNorm2(const InterfacePoint3& p) {
    return p.x * p.x + p.y * p.y + p.z * p.z;
}

static InterfacePoint3 interfaceFaceCofactorVector(
        const MeshBlockField& block,
        const HaloInterfaceFace& face,
        bool& ok) {
    ok = false;
    if (face.direction < 0 || face.direction > 2) return {};

    double metrics[4] = {0.0, 0.0, 0.0, 0.0};
    Math::faceMetrics(block.field,
                      face.i,
                      face.j,
                      face.k,
                      (Math::Dir)face.direction,
                      metrics);
    InterfacePoint3 n{metrics[0], metrics[1], metrics[2]};
    const double n2 = interfaceNorm2(n);
    ok = std::isfinite(n.x) && std::isfinite(n.y) &&
         std::isfinite(n.z) && n2 > 1.0e-300;
    return n;
}

/// @brief 仅打包Field的内部物理点（不含ghost层），与OpenCFD做法一致。
/// @return interior-only平铺buffer: [k][j][i][var]，interior索引 = (k*nx + j)*ny + i
static std::vector<double> packInteriorField(const Field& field) {
    int ng = field.NG();
    int nx = field.NX();
    int ny = field.NY();
    int nz = field.NZ();
    int nPoints = nx * ny * nz;
    std::vector<double> buffer((size_t)nPoints * (size_t)field.NVar(), 0.0);
    size_t p = 0;
    for (int k = ng; k < ng + nz; ++k) {
        for (int j = ng; j < ng + ny; ++j) {
            for (int i = ng; i < ng + nx; ++i) {
                for (int v = 0; v < field.NVar(); ++v) {
                    buffer[p++] = field(i, j, k, v);
                }
            }
        }
    }
    return buffer;
}

static void requireScalarCompatible(const Field& field,
                                    const std::vector<double>& values,
                                    const char* context) {
    if ((int)values.size() != field.TotalSize()) {
        broadcast("Fatal: ",
                  std::string(context)
                  + ": scalar value count does not match Field::TotalSize.");
        std::exit(1);
    }
}

static std::vector<double> packInteriorScalar(
        const Field& field,
        const std::vector<double>& values) {
    requireScalarCompatible(field, values, "haloExchange scalar pack");
    int ng = field.NG();
    int nx = field.NX();
    int ny = field.NY();
    int nz = field.NZ();
    int nPoints = nx * ny * nz;
    std::vector<double> buffer((size_t)nPoints, 0.0);
    size_t p = 0;
    for (int k = ng; k < ng + nz; ++k) {
        for (int j = ng; j < ng + ny; ++j) {
            for (int i = ng; i < ng + nx; ++i) {
                const int id = field.getIdx(i, j, k);
                const double value = values[(size_t)id];
                if (!std::isfinite(value)) {
                    broadcast("Fatal: ",
                              "haloExchange scalar pack found non-finite value.");
                    std::exit(1);
                }
                buffer[p++] = value;
            }
        }
    }
    return buffer;
}

static void writeInteriorScalar(const Field& field,
                                std::vector<double>& values,
                                int interiorIndex,
                                double value) {
    requireScalarCompatible(field, values, "haloExchange scalar write");
    if (interiorIndex < 0) return;
    int ng = field.NG();
    int nx = field.NX();
    int ny = field.NY();
    int nz = field.NZ();
    int nPoints = nx * ny * nz;
    if (interiorIndex >= nPoints) return;

    int i = interiorIndex % nx + ng;
    int j = (interiorIndex / nx) % ny + ng;
    int k = interiorIndex / (nx * ny) + ng;
    if (!std::isfinite(value)) {
        broadcast("Fatal: ",
                  "haloExchange scalar write received non-finite value.");
        std::exit(1);
    }
    values[(size_t)field.getIdx(i, j, k)] = value;
}

static void writeInteriorPoint(Field& field,
                               int interiorIndex,
                               const std::vector<double>& values) {
    if (interiorIndex < 0) return;
    int ng = field.NG();
    int nx = field.NX();
    int ny = field.NY();
    int nz = field.NZ();
    int nPoints = nx * ny * nz;
    if (interiorIndex >= nPoints) return;

    int i = interiorIndex % nx + ng;
    int j = (interiorIndex / nx) % ny + ng;
    int k = interiorIndex / (nx * ny) + ng;
    if ((int)values.size() != field.NVar()) return;
    for (int v = 0; v < field.NVar(); ++v) {
        field(i, j, k, v) = values[(size_t)v];
    }
}

static bool scalarInteriorDataOffset(int interior,
                                     int rank,
                                     const std::vector<int>& displacements,
                                     const std::vector<int>& counts,
                                     const std::vector<double>& allData,
                                     size_t& offset) {
    if (rank < 0 || rank >= (int)counts.size()) return false;
    if (interior < 0) return false;
    const size_t localOffset = (size_t)interior;
    if (localOffset + 1 > (size_t)counts[(size_t)rank]) return false;
    offset = (size_t)displacements[(size_t)rank] + localOffset;
    return offset < allData.size();
}

static bool interiorDataOffset(int interior,
                               int rank,
                               int nVar,
                               const std::vector<int>& displacements,
                               const std::vector<int>& counts,
                               const std::vector<double>& allData,
                               size_t& offset) {
    if (rank < 0 || rank >= (int)counts.size()) return false;
    if (interior < 0) return false;
    const size_t localOffset = (size_t)interior * (size_t)nVar;
    if (localOffset + (size_t)nVar > (size_t)counts[(size_t)rank]) return false;
    offset = (size_t)displacements[(size_t)rank] + localOffset;
    return offset + (size_t)nVar <= allData.size();
}

static void ownerCopyScalarInterfaceSyncGroups(
        const HaloExchangePlan& plan,
        const std::function<bool(int, int, size_t&)>& pointDataOffset,
        std::vector<double>& allData,
        const std::function<void(int, int, double)>& localWriter,
        const std::function<bool(const HaloInterfaceSyncGroup&)>&
            participates = {}) {
    for (const HaloInterfaceSyncGroup& group : plan.interfaceSyncGroups) {
        if (participates && !participates(group)) continue;
        if (group.points.empty() || group.canonicalOwner < 0 ||
            group.canonicalOwner >= (int)group.points.size()) {
            broadcast("Fatal: ",
                      "GlobalDof scalar copy has no valid canonical owner.");
            std::exit(1);
        }
        const HaloInterfacePoint& owner =
            group.points[(size_t)group.canonicalOwner];
        size_t ownerOffset = 0;
        if (!pointDataOffset(owner.blockId,
                             owner.interiorIndex,
                             ownerOffset)) {
            broadcast("Fatal: ",
                      "GlobalDof scalar owner is absent from the exchanged payload.");
            std::exit(1);
        }
        const double reference = allData[ownerOffset];
        if (!std::isfinite(reference)) {
            broadcast("Fatal: ",
                      "GlobalDof scalar owner produced a non-finite value.");
            std::exit(1);
        }

        for (const HaloInterfacePoint& point : group.points) {
            size_t offset = 0;
            if (!pointDataOffset(point.blockId,
                                 point.interiorIndex,
                                 offset)) {
                broadcast("Fatal: ",
                          "GlobalDof scalar replica is absent from the exchanged payload.");
                std::exit(1);
            }
            allData[offset] = reference;
            localWriter(point.blockId, point.interiorIndex, reference);
        }
    }
}

static bool validInterfaceFluxFace(const Field& field,
                                   const HaloInterfaceFace& face) {
    const int ng = field.NG();
    const int nx = field.NX();
    const int ny = field.NY();
    const int nz = field.NZ();
    if (face.direction == 0) {
        return face.i >= ng - 1 && face.i < ng + nx &&
               face.j >= ng && face.j < ng + ny &&
               face.k >= ng && face.k < ng + nz;
    }
    if (face.direction == 1) {
        return face.i >= ng && face.i < ng + nx &&
               face.j >= ng - 1 && face.j < ng + ny &&
               face.k >= ng && face.k < ng + nz;
    }
    if (face.direction == 2) {
        return face.i >= ng && face.i < ng + nx &&
               face.j >= ng && face.j < ng + ny &&
               face.k >= ng - 1 && face.k < ng + nz;
    }
    return false;
}

static void ownerCopyInterfaceSyncGroups(
        const HaloExchangePlan& plan,
        int nVar,
        const std::function<bool(int, int, size_t&)>& pointDataOffset,
        std::vector<double>& allData,
        const std::function<void(int, int,
                                 const std::vector<double>&)>&
            localWriter,
        const std::function<std::string(int, int)>& describePoint = {},
        const std::function<bool(const HaloInterfaceSyncGroup&)>&
            participates = {}) {
    for (const HaloInterfaceSyncGroup& group : plan.interfaceSyncGroups) {
        // 单 Field MPI 只接收拓扑邻居的数据。不能让一个不拥有本 group
        // 副本的 rank 读取 allData 中未通信的零填充；它对该 GlobalDof
        // 没有写入职责，也无需执行 owner-to-copy。
        if (participates && !participates(group)) continue;
        if (group.points.empty() || group.canonicalOwner < 0 ||
            group.canonicalOwner >= (int)group.points.size()) {
            broadcast("Fatal: ",
                      "GlobalDof state copy has no valid canonical owner.");
            std::exit(1);
        }
        const HaloInterfacePoint& owner =
            group.points[(size_t)group.canonicalOwner];
        size_t ownerOffset = 0;
        if (!pointDataOffset(owner.blockId,
                             owner.interiorIndex,
                             ownerOffset)) {
            broadcast("Fatal: ",
                      "GlobalDof state owner is absent from the exchanged payload.");
            std::exit(1);
        }
        std::vector<double> reference((size_t)nVar, 0.0);
        for (int v = 0; v < nVar; ++v) {
            reference[(size_t)v] = allData[ownerOffset + (size_t)v];
            if (!std::isfinite(reference[(size_t)v])) {
                broadcast("Fatal: ",
                          "GlobalDof state owner produced a non-finite value.");
                std::exit(1);
            }
        }
        if (nVar == 5) {
            const std::string context = describePoint
                ? std::string("GlobalDof canonical state at ")
                    + describePoint(owner.blockId, owner.interiorIndex)
                : std::string("GlobalDof canonical state");
            Numerics::requirePhysicalState(
                context.c_str(),
                reference[(size_t)RHO], reference[(size_t)RU],
                reference[(size_t)RV], reference[(size_t)RW],
                reference[(size_t)E]);
        }

        for (const HaloInterfacePoint& point : group.points) {
            size_t offset = 0;
            if (!pointDataOffset(point.blockId,
                                 point.interiorIndex,
                                 offset)) {
                broadcast("Fatal: ",
                          "GlobalDof state replica is absent from the exchanged payload.");
                std::exit(1);
            }

            for (int v = 0; v < nVar; ++v) {
                allData[offset + (size_t)v] = reference[(size_t)v];
            }
            localWriter(point.blockId, point.interiorIndex, reference);
        }
    }
}

static const DonorBlockInfo* findDonorBlockInfo(
        const std::vector<DonorBlockInfo>& infos,
        int blockId) {
    for (const DonorBlockInfo& info : infos) {
        if (info.blockId == blockId) return &info;
    }
    return nullptr;
}

template <typename Reader>
static bool tensorInterpolatedValue(const HaloCellMapping& mapping,
                                    const DonorBlockInfo& donorInfo,
                                    int component,
                                    Reader&& readComponent,
                                    double& value) {
    value = 0.0;
    bool hasValue = false;

    for (int axis = 0; axis < 3; ++axis) {
        const int start = mapping.donorStencilStart[(size_t)axis];
        const int size = mapping.donorStencilSize[(size_t)axis];
        const int axisSize =
            axis == 0 ? donorInfo.nx : (axis == 1 ? donorInfo.ny : donorInfo.nz);
        if (size <= 0 || size > kMaxHaloInterpolationStencil ||
            start < 0 || start + size > axisSize) {
            return false;
        }
        for (int n = 0; n < size; ++n) {
            if (!std::isfinite(
                    mapping.donorStencilWeights[(size_t)axis][(size_t)n])) {
                return false;
            }
        }
    }

    for (int kk = 0; kk < mapping.donorStencilSize[2]; ++kk) {
        const int z = mapping.donorStencilStart[2] + kk;
        const double wz =
            mapping.donorStencilWeights[2][(size_t)kk];
        for (int jj = 0; jj < mapping.donorStencilSize[1]; ++jj) {
            const int y = mapping.donorStencilStart[1] + jj;
            const double wy =
                mapping.donorStencilWeights[1][(size_t)jj];
            for (int ii = 0; ii < mapping.donorStencilSize[0]; ++ii) {
                const int x = mapping.donorStencilStart[0] + ii;
                const double wx =
                    mapping.donorStencilWeights[0][(size_t)ii];
                const int interior =
                    (z * donorInfo.ny + y) * donorInfo.nx + x;
                double donorValue = 0.0;
                if (!readComponent(interior, component, donorValue)) {
                    return false;
                }
                const double weight = wx * wy * wz;
                if (!std::isfinite(weight) || !std::isfinite(donorValue)) {
                    return false;
                }
                value += weight * donorValue;
                hasValue = true;
            }
        }
    }
    return hasValue && std::isfinite(value);
}

} // namespace

void HaloExchange::exchange(Field& field) const {
    const HaloBlockPlan* blockPlan = localPlan();
    if (!blockPlan ||
        (blockPlan->cells.empty() && plan_->interfaceSyncGroups.empty())) {
        return;
    }

#if SF_USE_MPI
    if (!mpiEnabled_) return;

    // ── Step 1: 仅打包interior数据（与OpenCFD一致）──
    std::vector<double> sendBuffer = packInteriorField(field);
    int sendCount = (int)sendBuffer.size();

    // 只和拓扑相邻 rank 交换；计数来自静态 halo plan，不在每个 stage
    // 重新 Allgather。
    std::vector<int> counts = rankInteriorCounts(field.NVar());
    if (counts[(size_t)mpiRank_] != sendCount) {
        throw std::runtime_error(
            "Single-field halo pack differs from its static rank count; "
            "use the multi-block exchange overload for multiple local patches.");
    }

    std::vector<int> displacements((size_t)mpiSize_, 0);
    int totalCount = 0;
    for (int r = 0; r < mpiSize_; ++r) {
        displacements[(size_t)r] = totalCount;
        totalCount += counts[(size_t)r];
    }

    std::vector<double> allData((size_t)totalCount, 0.0);
    exchangeNeighbourPayload(
        sendBuffer, counts, displacements, 4100, allData);

    // 重合物理点不是ghost。若plan中显式给出等价类，则先在通信快照
    // 中同步成一个共同状态，再用该快照填充ghost和同步副本。
    ownerCopyInterfaceSyncGroups(
        *plan_,
        field.NVar(),
        [&](int blockId, int interior, size_t& offset) {
            return interiorDataOffset(interior,
                                      ownerRank(blockId),
                                      field.NVar(),
                                      displacements,
                                      counts,
                                      allData,
                                      offset);
        },
        allData,
        [&](int blockId, int interior,
            const std::vector<double>& avg) {
            if (blockId == localBlockId_) {
                writeInteriorPoint(field, interior, avg);
            }
        }, {},
        [&](const HaloInterfaceSyncGroup& group) {
            return std::any_of(
                group.points.begin(), group.points.end(),
                [&](const HaloInterfacePoint& point) {
                    return point.blockId == localBlockId_;
                });
        });

    // ── Step 3: 预提取各donor rank的interior buffer ──
    // 为快速查找，直接存储每个rank的buffer起始指针和大小
    struct RankBuf {
        const double* data = nullptr;
        int count = 0;
    };
    std::vector<RankBuf> rankBufs((size_t)mpiSize_);
    for (int r = 0; r < mpiSize_; ++r) {
        rankBufs[(size_t)r].data = allData.data() + displacements[(size_t)r];
        rankBufs[(size_t)r].count = counts[(size_t)r];
    }

    // ── Step 4: 按mapping填充ghost cell并同步重合物理点 ──
    for (const HaloCellMapping& mapping : blockPlan->cells) {
        int donorRank = ownerRank(mapping.donorBlock);
        if (donorRank < 0 || donorRank >= mpiSize_) continue;

        const double* donorBuf = rankBufs[(size_t)donorRank].data;
        int donorBufCount = rankBufs[(size_t)donorRank].count;
        if (!donorBuf || donorBufCount <= 0) continue;

        int i = mapping.ownerIJK[0];
        int j = mapping.ownerIJK[1];
        int k = mapping.ownerIJK[2];

        for (int v = 0; v < field.NVar(); ++v) {
            if (mapping.kind == HaloMappingKind::DirectCopy) {
                int iIdx = mapping.donorInteriorIndex;
                if (iIdx < 0) continue;
                size_t p = (size_t)iIdx * (size_t)field.NVar() + (size_t)v;
                if (p < (size_t)donorBufCount) {
                    field(i, j, k, v) = donorBuf[p];
                }
                continue;
            }

            double value = 0.0;
            const DonorBlockInfo* donorInfo =
                findDonorBlockInfo(blockPlan->donorBlockInfos,
                                   mapping.donorBlock);
            if (!donorInfo ||
                !tensorInterpolatedValue(
                    mapping,
                    *donorInfo,
                    v,
                    [&](int interior, int component, double& out) {
                        const size_t p =
                            (size_t)interior * (size_t)field.NVar()
                            + (size_t)component;
                        if (p >= (size_t)donorBufCount) return false;
                        out = donorBuf[p];
                        return true;
                    },
                    value)) {
                broadcast("Fatal: ",
                          "haloExchange tensor interpolation failed.");
                std::exit(1);
            }
            field(i, j, k, v) = value;
        }
    }
#else
    (void)field;
#endif
}

void HaloExchange::exchangeScalarValues(
        const Field& field,
        std::vector<double>& values) const {
    const HaloBlockPlan* blockPlan = localPlan();
    if (!blockPlan ||
        (blockPlan->cells.empty() && plan_->interfaceSyncGroups.empty())) {
        return;
    }

    requireScalarCompatible(field, values, "haloExchange scalar");

#if SF_USE_MPI
    if (!mpiEnabled_) return;

    std::vector<double> sendBuffer = packInteriorScalar(field, values);
    int sendCount = (int)sendBuffer.size();

    std::vector<int> counts = rankInteriorCounts(1);
    if (counts[(size_t)mpiRank_] != sendCount) {
        throw std::runtime_error(
            "Single-field scalar halo pack differs from its static rank "
            "count; use the multi-block scalar exchange overload.");
    }

    std::vector<int> displacements((size_t)mpiSize_, 0);
    int totalCount = 0;
    for (int r = 0; r < mpiSize_; ++r) {
        displacements[(size_t)r] = totalCount;
        totalCount += counts[(size_t)r];
    }

    std::vector<double> allData((size_t)totalCount, 0.0);
    exchangeNeighbourPayload(
        sendBuffer, counts, displacements, 4101, allData);

    ownerCopyScalarInterfaceSyncGroups(
        *plan_,
        [&](int blockId, int interior, size_t& offset) {
            return scalarInteriorDataOffset(interior,
                                            ownerRank(blockId),
                                            displacements,
                                            counts,
                                            allData,
                                            offset);
        },
        allData,
        [&](int blockId, int interior, double avg) {
            if (blockId == localBlockId_) {
                writeInteriorScalar(field, values, interior, avg);
            }
        },
        [&](const HaloInterfaceSyncGroup& group) {
            return std::any_of(
                group.points.begin(), group.points.end(),
                [&](const HaloInterfacePoint& point) {
                    return point.blockId == localBlockId_;
                });
        });

    struct RankBuf {
        const double* data = nullptr;
        int count = 0;
    };
    std::vector<RankBuf> rankBufs((size_t)mpiSize_);
    for (int r = 0; r < mpiSize_; ++r) {
        rankBufs[(size_t)r].data = allData.data() + displacements[(size_t)r];
        rankBufs[(size_t)r].count = counts[(size_t)r];
    }

    for (const HaloCellMapping& mapping : blockPlan->cells) {
        int donorRank = ownerRank(mapping.donorBlock);
        if (donorRank < 0 || donorRank >= mpiSize_) {
            broadcast("Fatal: ", "haloExchange scalar found invalid donor rank.");
            std::exit(1);
        }

        const double* donorBuf = rankBufs[(size_t)donorRank].data;
        int donorBufCount = rankBufs[(size_t)donorRank].count;
        if (!donorBuf || donorBufCount <= 0) {
            broadcast("Fatal: ", "haloExchange scalar found empty donor buffer.");
            std::exit(1);
        }

        int i = mapping.ownerIJK[0];
        int j = mapping.ownerIJK[1];
        int k = mapping.ownerIJK[2];

        double value = 0.0;
        if (mapping.kind == HaloMappingKind::DirectCopy) {
            int iIdx = mapping.donorInteriorIndex;
            if (iIdx < 0 || iIdx >= donorBufCount) {
                broadcast("Fatal: ",
                          "haloExchange scalar direct-copy donor index is invalid.");
                std::exit(1);
            }
            value = donorBuf[(size_t)iIdx];
        } else {
            const DonorBlockInfo* donorInfo =
                findDonorBlockInfo(blockPlan->donorBlockInfos,
                                   mapping.donorBlock);
            if (!donorInfo ||
                !tensorInterpolatedValue(
                    mapping,
                    *donorInfo,
                    0,
                    [&](int interior, int, double& out) {
                        if (interior < 0 || interior >= donorBufCount) {
                            return false;
                        }
                        out = donorBuf[(size_t)interior];
                        return true;
                    },
                    value)) {
                broadcast("Fatal: ",
                          "haloExchange scalar tensor interpolation failed.");
                std::exit(1);
            }
        }
        if (!std::isfinite(value)) {
            broadcast("Fatal: ",
                      "haloExchange scalar produced non-finite halo value.");
            std::exit(1);
        }
        values[(size_t)field.getIdx(i, j, k)] = value;
    }
#else
    (void)field;
    (void)values;
#endif
}

void HaloExchange::exchange(std::vector<MeshBlockField>& blocks) const {
    if (!plan_ || blocks.empty()) return;

#if SF_USE_MPI
    if (!mpiEnabled_) return;

    const int nBlocks = (int)blocks.size();
    if (nBlocks == 0) return;
    const int nVar = blocks.front().field.NVar();
    for (const MeshBlockField& block : blocks) {
        if (block.field.NVar() != nVar) {
            broadcast("Fatal: ", "all MPI blocks must use the same FluidStateModel.");
            std::exit(1);
        }
    }

    std::vector<int> blockDataCounts((size_t)nBlocks, 0);
    for (int b = 0; b < nBlocks; ++b) {
        if (b < (int)plan_->blockInfos.size()) {
            blockDataCounts[(size_t)b] =
                plan_->blockInfos[(size_t)b].interiorDataCount;
        } else {
            const Field& f = blocks[(size_t)b].field;
            blockDataCounts[(size_t)b] = f.NX() * f.NY() * f.NZ() * nVar;
        }
    }

    std::vector<int> counts((size_t)mpiSize_, 0);
    for (int b = 0; b < nBlocks; ++b) {
        int owner = ownerRank(b);
        if (owner < 0 || owner >= mpiSize_) continue;
        counts[(size_t)owner] += blockDataCounts[(size_t)b];
    }

    std::vector<int> displacements((size_t)mpiSize_, 0);
    int totalCount = 0;
    for (int r = 0; r < mpiSize_; ++r) {
        displacements[(size_t)r] = totalCount;
        totalCount += counts[(size_t)r];
    }

    std::vector<int> blockOffsets((size_t)nBlocks, -1);
    std::vector<int> rankCursor((size_t)mpiSize_, 0);
    for (int b = 0; b < nBlocks; ++b) {
        int owner = ownerRank(b);
        if (owner < 0 || owner >= mpiSize_) continue;
        blockOffsets[(size_t)b] =
            displacements[(size_t)owner] + rankCursor[(size_t)owner];
        rankCursor[(size_t)owner] += blockDataCounts[(size_t)b];
    }

    std::vector<double> sendBuffer;
    sendBuffer.reserve((size_t)counts[(size_t)mpiRank_]);
    for (int b = 0; b < nBlocks; ++b) {
        if (ownerRank(b) != mpiRank_) continue;
        std::vector<double> packed = packInteriorField(blocks[(size_t)b].field);
        sendBuffer.insert(sendBuffer.end(), packed.begin(), packed.end());
    }
    int sendCount = (int)sendBuffer.size();
    if (sendCount != counts[(size_t)mpiRank_]) {
        broadcast("Fatal: ",
                  "haloExchange pack count mismatch for rank "
                  + std::to_string(mpiRank_));
        return;
    }

    std::vector<double> allData((size_t)totalCount, 0.0);
    exchangeNeighbourPayload(
        sendBuffer, counts, displacements, 4102, allData);

    ownerCopyInterfaceSyncGroups(
        *plan_,
        nVar,
        [&](int blockId, int interior, size_t& offset) {
            if (blockId < 0 || blockId >= nBlocks || interior < 0) {
                return false;
            }
            const int blockOffset = blockOffsets[(size_t)blockId];
            const int blockCount = blockDataCounts[(size_t)blockId];
            if (blockOffset < 0 ||
                (size_t)interior * (size_t)nVar + (size_t)nVar
                    > (size_t)blockCount) {
                return false;
            }
            offset = (size_t)blockOffset + (size_t)interior * (size_t)nVar;
            return offset + (size_t)nVar <= allData.size();
        },
        allData,
        [&](int blockId, int interior,
            const std::vector<double>& avg) {
            if (ownerRank(blockId) == mpiRank_) {
                writeInteriorPoint(blocks[(size_t)blockId].field,
                                   interior,
                                   avg);
            }
        }, {},
        [&](const HaloInterfaceSyncGroup& group) {
            return std::any_of(
                group.points.begin(), group.points.end(),
                [&](const HaloInterfacePoint& point) {
                    return ownerRank(point.blockId) == mpiRank_;
                });
        });

    auto readValue = [&](int blockId, int interiorIndex, int v,
                         double& value) {
        if (blockId < 0 || blockId >= nBlocks) return false;
        if (interiorIndex < 0 || v < 0 || v >= nVar) return false;
        int blockOffset = blockOffsets[(size_t)blockId];
        if (blockOffset < 0) return false;
        size_t p = (size_t)blockOffset
                 + (size_t)interiorIndex * (size_t)nVar
                 + (size_t)v;
        if (p >= allData.size()) return false;
        value = allData[p];
        return true;
    };

    for (int b = 0; b < nBlocks; ++b) {
        if (ownerRank(b) != mpiRank_) continue;
        const HaloBlockPlan* plan = blockPlan(b);
        if (!plan || plan->cells.empty()) continue;

        Field& field = blocks[(size_t)b].field;
        for (const HaloCellMapping& mapping : plan->cells) {
            int i = mapping.ownerIJK[0];
            int j = mapping.ownerIJK[1];
            int k = mapping.ownerIJK[2];

            for (int v = 0; v < nVar; ++v) {
                if (mapping.kind == HaloMappingKind::DirectCopy) {
                    double value = 0.0;
                    if (readValue(mapping.donorBlock,
                                  mapping.donorInteriorIndex,
                                  v,
                                  value)) {
                        field(i, j, k, v) = value;
                    }
                    continue;
                }

                double value = 0.0;
                const DonorBlockInfo* donorInfo =
                    findDonorBlockInfo(plan_->blockInfos,
                                       mapping.donorBlock);
                if (!donorInfo ||
                    !tensorInterpolatedValue(
                        mapping,
                        *donorInfo,
                        v,
                        [&](int interior, int component, double& out) {
                            return readValue(mapping.donorBlock,
                                             interior,
                                             component,
                                             out);
                        },
                        value)) {
                    broadcast("Fatal: ",
                              "haloExchange tensor interpolation lookup failed.");
                    std::exit(1);
                }
                field(i, j, k, v) = value;
            }
        }
    }
#else
    (void)blocks;
#endif
}

void HaloExchange::exchangeScalarValues(
        const std::vector<ScalarBlockValues>& views) const {
    exchangeScalarBlockValues(views, true);
}

void HaloExchange::exchangeScalarIdentifiers(
        const std::vector<ScalarBlockValues>& views) const {
    // IBM flag、ghost layer 等离散 metadata 服从 Eulerian GlobalDof
    // ownership；禁止在重复点上各自产生后再平均或保留相互矛盾的副本。
    exchangeScalarBlockValues(views, true);
}

void HaloExchange::exchangeScalarBlockValues(
        const std::vector<ScalarBlockValues>& views,
        bool synchronizeSharedPoints) const {
    if (!plan_ || views.empty()) return;

#if SF_USE_MPI
    if (!mpiEnabled_) return;

    int nBlocks = 0;
    for (const ScalarBlockValues& view : views) {
        if (view.blockId < 0 || !view.field || !view.values) {
            broadcast("Fatal: ",
                      "haloExchange scalar block view is incomplete.");
            std::exit(1);
        }
        requireScalarCompatible(*view.field, *view.values,
                                "haloExchange scalar block");
        nBlocks = std::max(nBlocks, view.blockId + 1);
    }
    if (nBlocks == 0) return;

    std::vector<const ScalarBlockValues*> byBlock((size_t)nBlocks, nullptr);
    for (const ScalarBlockValues& view : views) {
        if (byBlock[(size_t)view.blockId] != nullptr) {
            broadcast("Fatal: ",
                      "haloExchange scalar received duplicate block view.");
            std::exit(1);
        }
        byBlock[(size_t)view.blockId] = &view;
    }

    std::vector<int> blockDataCounts((size_t)nBlocks, 0);
    for (int b = 0; b < nBlocks; ++b) {
        if (b < (int)plan_->blockInfos.size()) {
            blockDataCounts[(size_t)b] =
                plan_->blockInfos[(size_t)b].interiorPointCount;
        } else if (byBlock[(size_t)b]) {
            const Field& f = *byBlock[(size_t)b]->field;
            blockDataCounts[(size_t)b] = f.NX() * f.NY() * f.NZ();
        }
    }

    std::vector<int> counts((size_t)mpiSize_, 0);
    for (int b = 0; b < nBlocks; ++b) {
        int owner = ownerRank(b);
        if (owner < 0 || owner >= mpiSize_) continue;
        counts[(size_t)owner] += blockDataCounts[(size_t)b];
    }

    std::vector<int> displacements((size_t)mpiSize_, 0);
    int totalCount = 0;
    for (int r = 0; r < mpiSize_; ++r) {
        displacements[(size_t)r] = totalCount;
        totalCount += counts[(size_t)r];
    }

    std::vector<int> blockOffsets((size_t)nBlocks, -1);
    std::vector<int> rankCursor((size_t)mpiSize_, 0);
    for (int b = 0; b < nBlocks; ++b) {
        int owner = ownerRank(b);
        if (owner < 0 || owner >= mpiSize_) continue;
        blockOffsets[(size_t)b] =
            displacements[(size_t)owner] + rankCursor[(size_t)owner];
        rankCursor[(size_t)owner] += blockDataCounts[(size_t)b];
    }

    std::vector<double> sendBuffer;
    sendBuffer.reserve((size_t)counts[(size_t)mpiRank_]);
    for (int b = 0; b < nBlocks; ++b) {
        if (ownerRank(b) != mpiRank_) continue;
        const ScalarBlockValues* view = byBlock[(size_t)b];
        if (!view) {
            broadcast("Fatal: ",
                      "haloExchange scalar missing local owner block view.");
            std::exit(1);
        }
        std::vector<double> packed =
            packInteriorScalar(*view->field, *view->values);
        sendBuffer.insert(sendBuffer.end(), packed.begin(), packed.end());
    }
    int sendCount = (int)sendBuffer.size();
    if (sendCount != counts[(size_t)mpiRank_]) {
        broadcast("Fatal: ",
                  "haloExchange scalar pack count mismatch for rank "
                  + std::to_string(mpiRank_));
        std::exit(1);
    }

    std::vector<double> allData((size_t)totalCount, 0.0);
    exchangeNeighbourPayload(
        sendBuffer, counts, displacements, 4103, allData);

    if (synchronizeSharedPoints) ownerCopyScalarInterfaceSyncGroups(
        *plan_,
        [&](int blockId, int interior, size_t& offset) {
            if (blockId < 0 || blockId >= nBlocks || interior < 0) {
                return false;
            }
            const int blockOffset = blockOffsets[(size_t)blockId];
            const int blockCount = blockDataCounts[(size_t)blockId];
            if (blockOffset < 0 || interior >= blockCount) return false;
            offset = (size_t)blockOffset + (size_t)interior;
            return offset < allData.size();
        },
        allData,
        [&](int blockId, int interior, double avg) {
            if (ownerRank(blockId) != mpiRank_) return;
            const ScalarBlockValues* view = byBlock[(size_t)blockId];
            if (!view) return;
            writeInteriorScalar(*view->field, *view->values, interior, avg);
        },
        [&](const HaloInterfaceSyncGroup& group) {
            return std::any_of(
                group.points.begin(), group.points.end(),
                [&](const HaloInterfacePoint& point) {
                    return ownerRank(point.blockId) == mpiRank_;
                });
        });

    auto readValue = [&](int blockId, int interiorIndex, double& value) {
        if (blockId < 0 || blockId >= nBlocks || interiorIndex < 0) {
            return false;
        }
        int blockOffset = blockOffsets[(size_t)blockId];
        int blockCount = blockDataCounts[(size_t)blockId];
        if (blockOffset < 0 || interiorIndex >= blockCount) return false;
        size_t p = (size_t)blockOffset + (size_t)interiorIndex;
        if (p >= allData.size()) return false;
        value = allData[p];
        return true;
    };

    for (int b = 0; b < nBlocks; ++b) {
        if (ownerRank(b) != mpiRank_) continue;
        const ScalarBlockValues* view = byBlock[(size_t)b];
        if (!view) continue;
        const HaloBlockPlan* plan = blockPlan(b);
        if (!plan || plan->cells.empty()) continue;

        const Field& field = *view->field;
        std::vector<double>& values = *view->values;
        for (const HaloCellMapping& mapping : plan->cells) {
            int i = mapping.ownerIJK[0];
            int j = mapping.ownerIJK[1];
            int k = mapping.ownerIJK[2];

            double value = 0.0;
            if (mapping.kind == HaloMappingKind::DirectCopy) {
                if (!readValue(mapping.donorBlock,
                               mapping.donorInteriorIndex,
                               value)) {
                    broadcast("Fatal: ",
                              "haloExchange scalar direct-copy donor lookup failed.");
                    std::exit(1);
                }
            } else {
                const DonorBlockInfo* donorInfo =
                    findDonorBlockInfo(plan_->blockInfos,
                                       mapping.donorBlock);
                if (!donorInfo ||
                    !tensorInterpolatedValue(
                        mapping,
                        *donorInfo,
                        0,
                        [&](int interior, int, double& out) {
                            return readValue(mapping.donorBlock,
                                             interior,
                                             out);
                        },
                        value)) {
                    broadcast("Fatal: ",
                              "haloExchange scalar tensor interpolation lookup failed.");
                    std::exit(1);
                }
            }
            if (!std::isfinite(value)) {
                broadcast("Fatal: ",
                          "haloExchange scalar produced non-finite multi-block value.");
                std::exit(1);
            }
            values[(size_t)field.getIdx(i, j, k)] = value;
        }
    }
#else
    (void)views;
    (void)synchronizeSharedPoints;
#endif
}

void HaloExchange::synchronizeCanonicalFaceFlux(
        const Field& field, std::vector<double>& values) const {
    if (!plan_ || plan_->interfaceFluxSyncGroups.empty()) return;
    if (values.size() != 3U * static_cast<size_t>(field.TotalSize())) {
        broadcast("Fatal: ",
                  "Eulerian canonical face flux has invalid storage size.");
        std::exit(1);
    }
#if SF_USE_MPI
    if (!mpiEnabled_) return;
    std::vector<int> counts((size_t)mpiSize_, 0);
    int participantCount = 0;
    for (const HaloInterfaceFluxSyncGroup& group :
         plan_->interfaceFluxSyncGroups) {
        for (const HaloInterfaceFace& face : group.faces) {
            const int owner = ownerRank(face.blockId);
            if (owner < 0 || owner >= mpiSize_) {
                broadcast("Fatal: ",
                          "Eulerian face flux references an invalid rank.");
                std::exit(1);
            }
            ++counts[(size_t)owner];
            ++participantCount;
        }
    }
    std::vector<int> displacements((size_t)mpiSize_, 0);
    int totalCount = 0;
    for (int rank = 0; rank < mpiSize_; ++rank) {
        displacements[(size_t)rank] = totalCount;
        totalCount += counts[(size_t)rank];
    }
    std::vector<int> offsets((size_t)participantCount, -1);
    std::vector<int> cursors((size_t)mpiSize_, 0);
    int participant = 0;
    for (const HaloInterfaceFluxSyncGroup& group :
         plan_->interfaceFluxSyncGroups) {
        for (const HaloInterfaceFace& face : group.faces) {
            const int owner = ownerRank(face.blockId);
            offsets[(size_t)participant++] = displacements[(size_t)owner]
                + cursors[(size_t)owner]++;
        }
    }

    std::vector<double> sendBuffer;
    sendBuffer.reserve((size_t)counts[(size_t)mpiRank_]);
    for (const HaloInterfaceFluxSyncGroup& group :
         plan_->interfaceFluxSyncGroups) {
        for (const HaloInterfaceFace& face : group.faces) {
            if (ownerRank(face.blockId) != mpiRank_) continue;
            if (face.blockId != localBlockId_) {
                broadcast("Fatal: ",
                          "Eulerian face-flux overload requires one local "
                          "block per MPI rank.");
                std::exit(1);
            }
            if (!validInterfaceFluxFace(field, face)
                || !std::isfinite(face.orientation)
                || std::abs(std::abs(face.orientation) - 1.0) > 1.0e-12) {
                broadcast("Fatal: ",
                          "Eulerian canonical face descriptor is invalid.");
                std::exit(1);
            }
            const size_t index = static_cast<size_t>(
                face.direction * field.TotalSize()
                + field.getIdx(face.i, face.j, face.k));
            const double value = values[index];
            if (!std::isfinite(value)) {
                broadcast("Fatal: ",
                          "Eulerian canonical face flux is non-finite.");
                std::exit(1);
            }
            sendBuffer.push_back(value);
        }
    }
    if ((int)sendBuffer.size() != counts[(size_t)mpiRank_]) {
        broadcast("Fatal: ",
                  "Eulerian canonical face-flux pack count mismatch.");
        std::exit(1);
    }
    std::vector<double> allData((size_t)totalCount, 0.0);
    exchangeNeighbourPayload(
        sendBuffer, counts, displacements, 4104, allData);

    participant = 0;
    for (const HaloInterfaceFluxSyncGroup& group :
         plan_->interfaceFluxSyncGroups) {
        if (group.faces.empty() || group.canonicalOwner < 0
            || group.canonicalOwner >= (int)group.faces.size()) {
            broadcast("Fatal: ",
                      "Eulerian canonical face owner is invalid.");
            std::exit(1);
        }
        const HaloInterfaceFace& canonical =
            group.faces[(size_t)group.canonicalOwner];
        const double canonicalValue=
            Numerics::CanonicalFace::orient(
                allData[(size_t)offsets[(size_t)(
                    participant+group.canonicalOwner)]],
                canonical.orientation);
        if (!std::isfinite(canonicalValue)) {
            broadcast("Fatal: ",
                      "Eulerian synchronized face flux is non-finite.");
            std::exit(1);
        }
        for (const HaloInterfaceFace& face : group.faces) {
            if (face.blockId == localBlockId_) {
                const size_t index=Numerics::CanonicalFace::index(
                    field,face.direction,face.i,face.j,face.k);
                values[index]=Numerics::CanonicalFace::orient(
                    canonicalValue,face.orientation);
            }
        }
        participant += (int)group.faces.size();
    }
#else
    (void)field;
    (void)values;
#endif
}

void HaloExchange::assembleCanonicalInterfaceFluxes(
        std::vector<MeshBlockField>& blocks,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) const {
    if (!plan_ || plan_->interfaceFluxSyncGroups.empty() || blocks.empty()) {
        return;
    }

#if SF_USE_MPI
    if (!mpiEnabled_) return;

    const int nBlocks = (int)blocks.size();
    if ((int)fluxes.size() != nBlocks || (int)residuals.size() != nBlocks) {
        throw std::runtime_error("canonical flux workspace block count mismatch.");
    }
    const int nVar = blocks.front().field.NVar();
    const int participantPackSize = nVar;
    for (const MeshBlockField& block : blocks) {
        if (block.field.NVar() != nVar) {
            broadcast("Fatal: ", "canonical interface blocks use different FluidStateModels.");
            std::exit(1);
        }
    }
    std::vector<int> counts((size_t)mpiSize_, 0);
    for (const HaloInterfaceFluxSyncGroup& group :
         plan_->interfaceFluxSyncGroups) {
        if (group.faces.empty() || group.canonicalOwner < 0
            || group.canonicalOwner >= (int)group.faces.size()) {
            broadcast("Fatal: ",
                      "canonical interface flux owner is invalid.");
            std::exit(1);
        }
        for (const HaloInterfaceFace& face : group.faces) {
            if (face.blockId < 0 || face.blockId >= nBlocks) {
                broadcast("Fatal: ",
                          "interface flux sync references an invalid block.");
                std::exit(1);
            }
            const int owner = ownerRank(face.blockId);
            if (owner < 0 || owner >= mpiSize_) {
                broadcast("Fatal: ",
                          "interface flux sync references an invalid owner "
                          "rank.");
                std::exit(1);
            }
        }
        const HaloInterfaceFace& canonical =
            group.faces[(size_t)group.canonicalOwner];
        counts[(size_t)ownerRank(canonical.blockId)] += participantPackSize;
    }

    std::vector<int> displacements((size_t)mpiSize_, 0);
    int totalCount = 0;
    for (int r = 0; r < mpiSize_; ++r) {
        displacements[(size_t)r] = totalCount;
        totalCount += counts[(size_t)r];
    }

    std::vector<int> groupOffsets(
        plan_->interfaceFluxSyncGroups.size(), -1);
    std::vector<int> rankCursor((size_t)mpiSize_, 0);
    for (size_t groupId = 0;
         groupId < plan_->interfaceFluxSyncGroups.size(); ++groupId) {
        const auto& group = plan_->interfaceFluxSyncGroups[groupId];
        const auto& canonical = group.faces[(size_t)group.canonicalOwner];
        const int owner = ownerRank(canonical.blockId);
        groupOffsets[groupId] =
            displacements[(size_t)owner] + rankCursor[(size_t)owner];
        rankCursor[(size_t)owner] += participantPackSize;
    }

    std::vector<double> sendBuffer;
    sendBuffer.reserve((size_t)counts[(size_t)mpiRank_]);
    for (const HaloInterfaceFluxSyncGroup& group :
         plan_->interfaceFluxSyncGroups) {
        const HaloInterfaceFace& face =
            group.faces[(size_t)group.canonicalOwner];
        if (ownerRank(face.blockId) != mpiRank_) continue;
        Field& field = blocks[(size_t)face.blockId].field;
        FluxField* fluxField = fluxes[(size_t)face.blockId];
        if (!fluxField) throw std::runtime_error("missing canonical flux workspace.");
        if (!validInterfaceFluxFace(field, face)) {
            broadcast("Fatal: ",
                      "interface flux sync references an invalid "
                      "canonical ConvectiveFlux index.");
            std::exit(1);
        }
        if (!std::isfinite(face.orientation) ||
            std::abs(std::abs(face.orientation) - 1.0) > 1.0e-12) {
            broadcast("Fatal: ",
                      "interface flux sync has invalid canonical orientation.");
            std::exit(1);
        }
        for (int v = 0; v < nVar; ++v) {
            const double flux = (*fluxField)(
                (size_t)face.direction * field.TotalSize()
                + field.getIdx(face.i, face.j, face.k), v);
            if (!std::isfinite(flux)) {
                broadcast("Fatal: ",
                          "canonical owner produced a non-finite WENO "
                          "face flux.");
                std::exit(1);
            }
            sendBuffer.push_back(flux);
        }
    }
    const int sendCount = (int)sendBuffer.size();
    if (sendCount != counts[(size_t)mpiRank_]) {
        broadcast("Fatal: ",
                  "shared-interface flux pack count mismatch for rank "
                  + std::to_string(mpiRank_));
        std::exit(1);
    }

    std::vector<double> allData((size_t)totalCount, 0.0);
    exchangeNeighbourPayload(
        sendBuffer, counts, displacements, 4105, allData);

    for (size_t groupId = 0;
         groupId < plan_->interfaceFluxSyncGroups.size(); ++groupId) {
        const HaloInterfaceFluxSyncGroup& group =
            plan_->interfaceFluxSyncGroups[groupId];
        if (group.faces.empty()) continue;
        if (group.canonicalOwner < 0 ||
            group.canonicalOwner >= (int)group.faces.size()) {
            broadcast("Fatal: ",
                      "canonical interface flux owner is invalid.");
            std::exit(1);
        }
        double canonicalMetricNorm = 0.0;
        for (double component : group.canonicalCofactor) {
            if (!std::isfinite(component)) {
                broadcast("Fatal: ",
                          "canonical interface face geometry is non-finite.");
                std::exit(1);
            }
            canonicalMetricNorm += component * component;
        }
        if (canonicalMetricNorm <= 1.0e-300) {
            broadcast("Fatal: ",
                      "canonical interface face geometry is degenerate.");
            std::exit(1);
        }

        std::vector<double> canonicalFlux((size_t)nVar, 0.0);
        const HaloInterfaceFace& canonicalFace =
            group.faces[(size_t)group.canonicalOwner];
        const size_t canonicalOffset = (size_t)groupOffsets[groupId];
        if (canonicalOffset + participantPackSize > allData.size()) {
            broadcast("Fatal: ",
                      "canonical interface flux offset is invalid.");
            std::exit(1);
        }
        for (int v = 0; v < nVar; ++v) {
            const double ownerFlux =
                allData[canonicalOffset + (size_t)v];
            canonicalFlux[(size_t)v] =
                Numerics::CanonicalFace::orient(
                    ownerFlux,
                    canonicalFace.orientation);
            if (!std::isfinite(canonicalFlux[(size_t)v])) {
                broadcast("Fatal: ",
                          "shared-interface canonical owner flux is "
                          "non-finite.");
                std::exit(1);
            }
        }

        for (int localFace = 0; localFace < (int)group.faces.size();
             ++localFace) {
            const HaloInterfaceFace& face = group.faces[(size_t)localFace];
            if (ownerRank(face.blockId) == mpiRank_) {
                Field& field = blocks[(size_t)face.blockId].field;
                FluxField* fluxField = fluxes[(size_t)face.blockId];
                if (!fluxField) throw std::runtime_error("missing canonical flux workspace.");
                if (!validInterfaceFluxFace(field, face)) {
                    broadcast("Fatal: ",
                              "interface flux sync references an invalid "
                              "ConvectiveFlux index.");
                    std::exit(1);
                }
                for (int v = 0; v < nVar; ++v) {
                    const double localValue =
                        Numerics::CanonicalFace::orient(
                            canonicalFlux[(size_t)v], face.orientation);
                    if (!std::isfinite(localValue)) {
                        broadcast("Fatal: ",
                                  "common interface flux synchronization produced a "
                                  "non-finite local flux.");
                        std::exit(1);
                    }
                    (*fluxField)((size_t)face.direction * field.TotalSize()
                                 + field.getIdx(face.i, face.j, face.k), v) = localValue;
                }
            }
        }
    }

    for (int b = 0; b < nBlocks; ++b) {
        if (ownerRank(b) == mpiRank_) {
            if (!fluxes[(size_t)b] || !residuals[(size_t)b]) {
                throw std::runtime_error("missing local face-flux/residual workspace.");
            }
            Math::assembleConvectiveFluxResidual(
                blocks[(size_t)b].field, *fluxes[(size_t)b], *residuals[(size_t)b]);
        }
    }
#else
    (void)blocks;
#endif
}

} // namespace Parallel
} // namespace SF
