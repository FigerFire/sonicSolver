/// @file SF_communicationPlan.h
/// @brief 网格接口通信计划与 donor/receiver 拓扑描述。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_meshInterface.h"

#include <array>
#include <string>
#include <vector>

namespace SF {

struct MeshBlockField;

namespace MeshCommunication {

constexpr int kMaxHaloInterpolationStencil = 6;

/// @brief halo映射类型。
enum class HaloMappingKind {
    DirectCopy,
    Interpolated
};

struct HaloCellMapping {
    HaloMappingKind kind = HaloMappingKind::Interpolated;

    int ownerBlock = -1;
    int ownerIndex = -1;
    std::array<int, 3> ownerIJK{0, 0, 0};

    int donorBlock = -1;
    int donorIndex = -1;
    std::array<int, 3> donorCellIJK{0, 0, 0};
    std::array<double, 3> localCoord{0.0, 0.0, 0.0};

    /// @brief donor field 的 interior-only 平铺索引（仅物理点，不含ghost）。
    /// 用于 exchange 时从 interior-packed buffer 中取数据。
    int donorInteriorIndex = -1;

    /// @brief donor结构坐标中的高阶张量插值模板起点（interior局部索引）。
    std::array<int, 3> donorStencilStart{0, 0, 0};
    /// @brief donor结构坐标中每个方向采用的一维Lagrange点数。
    std::array<int, 3> donorStencilSize{0, 0, 0};
    /// @brief 每个方向的一维Lagrange权重，张量积在通信填ghost时现场累加。
    std::array<std::array<double, kMaxHaloInterpolationStencil>, 3>
        donorStencilWeights{};
};

struct DonorBlockInfo {
    int blockId = -1;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int ng = 0;
    int interiorPointCount = 0;
    int interiorDataCount = 0;  ///< interiorPointCount * variableCount
    int variableCount = 5;      ///< FluidStateModel 守恒变量数。
};

struct HaloInterfacePoint {
    int blockId = -1;
    int interiorIndex = -1;
    /// @brief 当前结构 patch 对该 GlobalDof 的物理对偶体积片段。
    double dualVolume = 0.0;
};

struct HaloInterfaceSyncGroup {
    /// @brief 同一 GlobalDof 在多个结构 patch 中的冗余表示。
    std::vector<HaloInterfacePoint> points;
    int globalDofId = -1;
    /// @brief points 中唯一产生 canonical state/metadata 的副本。
    int canonicalOwner = -1;
    int ownerRank = -1;
};

struct HaloInterfaceFace {
    int blockId = -1;
    int direction = -1;
    int i = 0;
    int j = 0;
    int k = 0;
    double orientation = 1.0;
    /// @brief 该 face 在本 patch 通量差中的关联符号：高侧 +1，低侧 -1。
    double residualSign = 0.0;
};

struct HaloInterfaceFluxSyncGroup {
    /// @brief 同一物理界面上的结构半面通量DOF。
    std::vector<HaloInterfaceFace> faces;
    /// @brief `faces[canonicalOwner]` 唯一计算并广播数值通量。
    int canonicalOwner = 0;
    /// @brief 预处理阶段去重后分配的稳定 GlobalFace 标识。
    int globalFaceId = -1;
    /// @brief canonical owner 使用的曲面余因子向量。
    std::array<double, 3> canonicalCofactor{0.0, 0.0, 0.0};
    /// @brief canonical owner 面上的逆 Jacobian。
    double canonicalInverseJacobian = 0.0;
};

struct HaloBlockPlan {
    int blockId = -1;
    std::vector<HaloCellMapping> cells;
    /// @brief 本块需要通信的所有donor块的信息（含当前块自身）。
    std::vector<DonorBlockInfo> donorBlockInfos;
};

struct HaloExchangePlan {
    double tolerance = 1.0e-8;
    /// @brief 本次数值格式实际需要交换的 halo 层数。
    int haloWidth = 1;
    int searchedGhostCells = 0;
    int mappedGhostCells = 0;
    int directMappedGhostCells = 0;
    int interpolatedMappedGhostCells = 0;
    int unmappedGhostCells = 0;
    int minInterpolationTensorPoints = 0;
    int maxInterpolationTensorPoints = 0;
    int synchronizedPhysicalPoints = 0;
    int synchronizedInterfaceFluxes = 0;
    /// @brief 每个全局block由哪个MPI rank推进；为空时兼容旧的blockId%rank规则。
    std::vector<int> blockOwnerRanks;
    /// @brief 每个全局block的interior打包尺寸。
    std::vector<DonorBlockInfo> blockInfos;
    std::vector<HaloBlockPlan> blockPlans;
    /// @brief GlobalDof 等价类；交换时执行 owner-to-copy。
    std::vector<HaloInterfaceSyncGroup> interfaceSyncGroups;
    /// @brief 重合物理界面上的面通量等价类；RHS装配时按owner面同步F_hat。
    std::vector<HaloInterfaceFluxSyncGroup> interfaceFluxSyncGroups;
    /// @brief 显式结构接口拓扑；globalPointId 只用于发现候选，不代替该映射。
    std::vector<PatchInterfaceTopology> interfaces;

    void clear() {
        haloWidth = 1;
        searchedGhostCells = 0;
        mappedGhostCells = 0;
        directMappedGhostCells = 0;
        interpolatedMappedGhostCells = 0;
        unmappedGhostCells = 0;
        minInterpolationTensorPoints = 0;
        maxInterpolationTensorPoints = 0;
        synchronizedPhysicalPoints = 0;
        synchronizedInterfaceFluxes = 0;
        blockOwnerRanks.clear();
        blockInfos.clear();
        blockPlans.clear();
        interfaceSyncGroups.clear();
        interfaceFluxSyncGroups.clear();
        interfaces.clear();
    }
};

struct HaloPreprocessOptions {
    double tolerance = 1.0e-8;
    /// @brief 由 convection/ILW numerics policy 决定，不能写死为 Field::NG。
    int requiredHaloWidth = 1;
    /// @brief 当前主线只允许完全对应接口；非匹配接口必须显式拒绝。
    bool allowNonMatchingInterfaces = false;
    bool useDeclaredHaloSets = true;
    std::vector<std::string> physicalBoundaryNames;
};

bool buildHaloExchangePlan(std::vector<MeshBlockField>& blocks,
                           const HaloPreprocessOptions& options,
                           HaloExchangePlan& plan);

} // namespace MeshCommunication
} // namespace SF
