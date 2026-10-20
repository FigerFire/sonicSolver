/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.01-----------*/

#pragma once

/// @file SF_cellFaceMesh.h
/// @brief 轻量 points/cells/faces/patches 网格拓扑。

#include <array>
#include <map>
#include <string>
#include <vector>

namespace SF {

/// @brief 三维点坐标。
struct MeshPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/// @brief 六面体单元拓扑。
struct MeshCell {
    std::array<int, 8> pointIds{-1, -1, -1, -1, -1, -1, -1, -1};
    std::array<int, 6> faceIds{-1, -1, -1, -1, -1, -1};
    int blockId = -1;
    std::array<int, 3> ijk{0, 0, 0};
};

/// @brief 四边形面拓扑和owner/neighbour关系。
struct MeshFace {
    std::array<int, 4> pointIds{-1, -1, -1, -1};
    int ownerCell = -1;
    int neighbourCell = -1;
    int patchId = -1;
    int direction = -1;
    int ownerRank = 0;
    int neighbourRank = -1;
};

/// @brief 边界patch,保存一组面。
struct MeshPatch {
    std::string name;
    std::string type;
    std::vector<int> faceIds;
};

/// @brief 轻量cell-face网格。
///
/// 该类是结构FDM求解器和后续face-owned守恒通量之间的拓扑层。
/// 它不替代当前结构stencil视图，而是提供cells/faces/patches和
/// owner-neighbour关系，供输出、MPI界面通量和守恒审计使用。
class CellFaceMesh {
public:
    /// @brief 从结构块节点坐标和点集合生成cell-face拓扑。
    /// @param nx,ny,nz 结构节点数量,不是单元数量。
    /// @param x,y,z 节点坐标数组,布局为[k][j][i]。
    /// @param pointSets 旧SFM点集合,用于生成patch面。
    /// @param blockId 所属结构块编号。
    /// @param ownerRank 所属MPI rank。
    static CellFaceMesh fromStructuredPoints(
        int nx, int ny, int nz,
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::map<std::string, std::vector<int>>& pointSets,
        int blockId = 0,
        int ownerRank = 0);

    /// @brief 清空所有拓扑。
    void clear();

    /// @brief 由patch面反推每个patch覆盖的点索引集合。
    /// @param includeAll 是否额外生成覆盖全块点的`all`集合。
    /// @return key=patch名,value=局部点索引。
    std::map<std::string, std::vector<int>> patchPointSets(
        bool includeAll = true) const;

    /// @brief 网格是否没有拓扑。
    bool empty() const { return points.empty() && cells.empty() && faces.empty(); }

    std::vector<MeshPoint> points;
    std::vector<MeshCell> cells;
    std::vector<MeshFace> faces;
    std::vector<MeshPatch> patches;

private:
    static int pointIndex(int i, int j, int k, int nx, int ny);
};

} // namespace SF
