#pragma once

/// @file SF_meshInterface.h
/// @brief 网格 patch 接口的通用拓扑契约；不包含通信后端。

#include <array>

namespace SF::MeshCommunication {

/// @brief patch 边界的拓扑类型；物理边界不进入 coupled communication plan。
enum class PatchBoundaryKind {
    Physical,
    CoupledConformal,
    Periodic,
    NonMatching
};

/// @brief 两个 patch 的接触维数。
enum class InterfaceRelation {
    FaceToFace,
    EdgeToEdge,
    PointToPoint
};

/// @brief 结构 patch 的一个有向边界面。
struct StructuredPatchSide {
    int axis = -1;
    int sign = 0;
    std::array<int, 3> begin{0, 0, 0};
    std::array<int, 3> end{0, 0, 0};
};

/// @brief owner 索引轴到 neighbour 索引轴的置换和方向。
struct StructuredIndexTransform {
    std::array<int, 3> neighbourAxisForOwner{0, 1, 2};
    std::array<int, 3> orientation{1, 1, 1};
};

/// @brief 与通信实现无关的 patch-patch 网格接口。
struct MeshInterface {
    int ownerBlock = -1;
    int neighbourBlock = -1;
    PatchBoundaryKind kind = PatchBoundaryKind::Physical;
    InterfaceRelation relation = InterfaceRelation::PointToPoint;
    StructuredPatchSide ownerSide;
    StructuredPatchSide neighbourSide;
    StructuredIndexTransform indexTransform;
    int haloWidth = 0;
    bool vectorComponentsAreGlobalCartesian = true;
    int canonicalGeometryOwner = -1;
    std::array<double, 3> canonicalCofactor{0.0, 0.0, 0.0};
};

/// @brief 跨执行单元的网格接口；通信由 infrastructure backend 完成。
struct ProcessorInterface : MeshInterface {
    int ownerRank = -1;
    int neighbourRank = -1;
};

/// @brief 同一执行单元内两个结构 block 的接口。
struct BlockInterface : MeshInterface {};

/// @brief 周期映射接口。
struct PeriodicInterface : MeshInterface {
    std::array<double, 3> translation{0.0, 0.0, 0.0};
};

using PatchInterfaceTopology = MeshInterface;

} // namespace SF::MeshCommunication
