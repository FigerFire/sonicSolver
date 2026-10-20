/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.13-----------*/

#pragma once

/// @file SF_meshDecompose.h
/// @brief 多zone结构网格的全局拓扑组装与MPI分区。

#include "SF_MultiBlockMesh.h"

#include <array>
#include <vector>

namespace SF {
namespace MeshDecompose {

/// @brief 输入zone合并后的全局拓扑统计。
struct CompositeMeshSummary {
    int sourceZoneCount = 0;
    int uniquePointCount = 0;
    int sharedPointGroupCount = 0;
    int sharedPointReferenceCount = 0;
    int maxPointMultiplicity = 1;
};

/// @brief 一个结构patch在源zone逻辑索引中的范围。
struct SourcePatchExtent {
    int patchId = -1;
    int sourceZoneId = -1;
    /// @brief 当前patch第一个真实点在源zone中的逻辑索引。
    std::array<int, 3> start{0, 0, 0};
    /// @brief 当前patch拥有的真实点数；相邻MPI patch共享切分面真实点。
    std::array<int, 3> size{0, 0, 0};
    /// @brief 当前patch拥有的第一个源cell逻辑索引。
    std::array<int, 3> cellStart{0, 0, 0};
    /// @brief 当前patch拥有的源cell数量，用于唯一cell覆盖检查。
    std::array<int, 3> cellSize{0, 0, 0};
    /// @brief 源zone中由本zone拥有的真实点起始索引，排除让给相邻源zone的接口点层。
    std::array<int, 3> sourceOwnedStart{0, 0, 0};
    /// @brief 源zone中由本zone拥有的真实点数量，排除让给相邻源zone的接口点层。
    std::array<int, 3> sourceOwnedSize{0, 0, 0};
    /// @brief 源zone的`#Information`点数。
    std::array<int, 3> sourceSize{0, 0, 0};
};

/// @brief 将输入zone组装为统一复合网格拓扑。
///
/// 重合点获得相同globalPointId；同一全局点上的sets取并集并传播到所有
/// 结构副本。zone只保留为有限差分计算坐标patch，不再承担MPI所有权语义。
/// @param zones 输入并原地写回全局点编号和合并后的sets。
/// @param tolerance 物理坐标合并容差。
/// @param summary 输出拓扑统计。
/// @return 成功返回true；无效坐标或容差时返回false。
bool assembleCompositeMesh(std::vector<RawMeshBlock>& zones,
                           double tolerance,
                           CompositeMeshSummary& summary);

/// @brief 按源 zone 结构索引 split 将复合网格分为 MPI partitions。
///
/// 返回的patches仍是结构数组，供有限差分核使用；partitions才是MPI所有权
/// 单位，一个partition可包含多个来自不同源zone的结构patch。切分位置由
/// `#Information` 的 I/J/K cell 范围直接生成，相邻 patch 共享切分点层；
/// 物理坐标不参与切分面反推。
/// @param zones 已完成assembleCompositeMesh的输入zone。
/// @param parts `{px, py, pz}` MPI分区数。
/// @param tolerance 切分面匹配容差。
/// @param patches 输出结构计算patch。
/// @param extents 输出patch在源zone中的逻辑范围。
/// @param partitions 输出MPI partition及其patch列表。
/// @return 拓扑不支持该split、切分面与结构索引不相容或存在空partition时返回false。
bool decomposeCompositeMesh(
    const std::vector<RawMeshBlock>& zones,
    const std::array<int, 3>& parts,
    double tolerance,
    std::vector<RawMeshBlock>& patches,
    std::vector<SourcePatchExtent>& extents,
    std::vector<MeshPartition>& partitions);

/// @brief 用同源zone邻patch的真实坐标刷新分区ghost坐标。
bool copyPartitionGhostCoordinates(
    std::vector<MeshBlockField>& patches,
    const std::vector<SourcePatchExtent>& extents);

} // namespace MeshDecompose
} // namespace SF
