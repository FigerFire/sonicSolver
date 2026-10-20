#pragma once

/// @file SF_fieldComponents.h
/// @brief Field 组合拥有的网格、边界和 IBM 存储对象。

#include "SF_valueTypes.h"

#include <map>
#include <string>
#include <vector>

namespace SF {

/// @brief 只保存坐标、度规、到壁面的有符号距离和 canonical face geometry。
struct GeometryStorage {
    std::vector<double> x, y, z;
    std::vector<double> jac;
    std::vector<double> xiX, xiY, xiZ;
    std::vector<double> etaX, etaY, etaZ;
    std::vector<double> zetaX, zetaY, zetaZ;
    /// @brief 到最近壁面的有符号距离（湍流模型等消费；非 IBM 专用标量）。
    std::vector<double> wallDistance;
    std::vector<double> canonicalFaceMetrics;
    std::vector<unsigned char> canonicalFaceMetricMask;
    /// @brief canonical face 仅允许唯一 owner 计算数值通量。
    std::vector<unsigned char> canonicalFaceOwnerMask;
};

/// @brief 物理边界集合、求解 mask 和通信分类元数据。
struct BoundaryMetadata {
    std::map<std::string, std::vector<int>> sets;
    std::map<std::string, BCType> faceMarkers;
    std::vector<int> cellType;
    std::vector<unsigned char> solverMask;
    std::vector<unsigned char> communicationHaloMask;
    /// @brief interior 点对应的全局 Eulerian 自由度；ghost 默认为 -1。
    std::vector<int> globalDofId;
    /// @brief GlobalDof 的唯一 MPI owner rank。
    std::vector<int> globalDofOwnerRank;
    /// @brief 当前 Field 副本是否为该 GlobalDof 的唯一 canonical owner。
    std::vector<unsigned char> globalDofOwnerMask;
};

/// @brief IBM 分类以及壁面、镜像点、法向和运动速度几何已迁出 Field。
/// 现由 models/ibm/topology/SF_ibmTopology.h 的 IBM::IBMGeometry 权威持有；
/// Field 只保留与 CellFlag 同类的 deferred 分类标记 IBMFluidMask。

} // namespace SF
