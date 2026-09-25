/// @file SF_meshGen.h
/// @brief 结构网格 block、edge 与生成流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <array>
#include <vector>
#include <string>
#include "SF_valueTypes.h"

namespace SF {

// ============================================================
//  网格参数数据结构 (解析自 OpenFOAM system/blockMeshDict)
// ============================================================

struct MeshBlock {
    int    verts[8];
    int    cells[3];
    double grade[3];
};

enum class EdgeType { LINE, ARC };

struct MeshEdge {
    EdgeType type   = EdgeType::LINE;
    int      v0 = 0, v1 = 0;     // 端点顶点索引
    Vector3  arcCenter;          // ARC 圆心
    Vector3  arcP;               // ARC 辅助点 p (过三点圆弧)
    double   arcAngle = 0.0;     // ARC 角度 (rad)
};

struct MeshFacePatch {
    std::string name;
    std::vector<std::array<int, 4>> faces;
};

struct MeshParameters {
    double                 scale = 1.0;
    int                    nGhost = 3;
    std::vector<Vector3>   vertices;
    std::vector<MeshBlock> blocks;
    std::vector<MeshEdge>  edges;
    std::vector<MeshFacePatch> facePatches;
};

// ============================================================
//  结构化网格生成
// ============================================================

namespace MeshGen {

/// 从 meshParameters 生成所有 block 的结构化点
/// @param params  解析后的网格参数
/// @param allX,allY,allZ  输出: 拼接后的点坐标
/// @param blockSizes      输出: 每个 block 的 {nx+1, ny+1, nz+1}
void generateStructuredMesh(const MeshParameters& params,
                            std::vector<double>& allX,
                            std::vector<double>& allY,
                            std::vector<double>& allZ,
                            std::vector<int>&    blockSizes);

/// 生成单个 block 的网格点 (trilinear 插值)
/// @param verts  8 个物理顶点坐标
/// @param nc     3 方向的单元数 {nx, ny, nz}
/// @param grade  3 方向的渐变比
/// @param x,y,z  输出: (nc[0]+1)*(nc[1]+1)*(nc[2]+1) 个点
void generateBlockPoints(const double verts[8][3],
                         const int    nc[3],
                         const double grade[3],
                         std::vector<double>& x,
                         std::vector<double>& y,
                         std::vector<double>& z);

/// 用 TFI 生成带弧边的 block (edges[12] 按块12条边排列)
/// v0-v1 映射关系: 0-1(xi), 3-2(xi), 4-5(xi), 7-6(xi)
///                 0-3(eta), 1-2(eta), 4-7(eta), 5-6(eta)
///                 0-4(zeta), 1-5, 3-7, 2-6
void generateTFIBlock(const double verts[8][3],
                      const int    nc[3],
                      const MeshEdge edges[12],
                      std::vector<double>& x,
                      std::vector<double>& y,
                      std::vector<double>& z);

/// 含 edges 的网格生成入口
void generateStructuredMeshTFI(const MeshParameters& params,
                               std::vector<double>& allX,
                               std::vector<double>& allY,
                               std::vector<double>& allZ,
                               std::vector<int>&    blockSizes);

/// @brief 输出生成网格。
///
/// 单zone和多zone都写入一个自包含的base.sfm；多zone可视化仍使用
/// base_blockXXX.vts和base.vtm。
bool writeCombinedMesh(const std::string& baseName,
                       const MeshParameters& params,
                       const std::vector<double>& allX,
                       const std::vector<double>& allY,
                       const std::vector<double>& allZ,
                       const std::vector<int>& blockSizes,
                       int nGhost = 3);

} // namespace MeshGen
} // namespace SF
