/// @file SF_mesh.h
/// @brief 结构/多块网格拓扑、度量与加载实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include "core/field/SF_field.h"
#include "infrastructure/io/SF_IO.h"
#include "infrastructure/io/SF_caseConfig.h"
#include "SF_meshConfig.h"

/// @file SF_mesh.h
/// @brief 结构网格的初始化、度规计算与虚胞生成。
///
/// Mesh类提供静态方法完成: 物理坐标→度规(Jacobian/度量系数)的计算、
/// 虚胞层坐标外推、边界标签同步、以及从原始坐标构建完整Field。

namespace SF {
    class MultiBlockMesh;
    struct RawMeshBlock;

    class Mesh {
    public:
        /// @brief 进入纯网格生成模式(createMesh=true时调用)。
        /// @param io IO管理器；createMesh 固定读取 system/blockMeshDict。
        /// @return 生成成功返回true。
        static bool createMesh(
            const CaseConfig& caseConfig,
            const MeshRuntimeConfig& config);

        /// @brief 核心函数：读取外部网格并立即计算度规(单块模式)。
        /// @param field 输出的Field对象，将被完全初始化。
        /// @param io IO管理器，包含网格文件列表。
        /// @return 成功返回true。
        static bool setupComplexMesh(
            SF::Field& field,
            const CaseConfig& caseConfig,
            const MeshRuntimeConfig& config);

        /// @brief 从IO读取多块网格并填充MultiBlockMesh。
        /// @param mesh 输出的多块网格对象。
        /// @param io IO管理器。
        /// @param preprocessor 复合拓扑装配后、MPI分解前的可选预处理。
        /// @return 成功返回true。
        static bool setupMultiBlockMesh(
            SF::MultiBlockMesh& mesh,
            const CaseConfig& caseConfig,
            const MeshRuntimeConfig& config,
            const std::function<bool(std::vector<RawMeshBlock>&)>&
                preprocessor = {});

        /// @brief 从原始坐标数组构建完整Field(含度规、虚胞、初始条件)。
        /// @param field 输出的Field对象。
        /// @param nx, ny, nz 内部网格维度(不含虚胞)。
        /// @param ng 虚胞层数。
        /// @param px, py, pz 物理坐标数组(长度=nx*ny*nz)。
        /// @param sets 边界集合映射(key=标签名, value=局部索引列表)。
        /// @param label 日志标签(默认"mesh block")。
        /// @return 成功返回true。
        static bool setupFieldFromRaw(SF::Field& field,
                                      int nx, int ny, int nz, int ng,
                                      const std::vector<double>& px,
                                      const std::vector<double>& py,
                                      const std::vector<double>& pz,
                                      const std::map<std::string, std::vector<int>>& sets,
                                      const MeshRuntimeConfig& config,
                                      const std::string& label = "mesh block");

        /// @brief 在更新虚胞坐标后重新计算度规。
        /// @param field 已经完成坐标填充的结构网格场。
        static void refreshMetrics(SF::Field& field);

        /// @brief 强制审计曲线网格保守度规恒等式的离散残差。
        /// @param field 已计算度规的网格场。
        static void auditMetricIdentity(SF::Field& field);

        // ── 虚胞工具函数 ──

        /// @brief 将实胞边界标签同步到新生成的虚胞点。
        /// @param field 网格场。
        /// @param tagMap 实胞索引→标签名的查找表。
        /// @param ri, rj, rk 参考实胞的(i,j,k)全局索引。
        /// @param gi, gj, gk 目标虚胞的(i,j,k)全局索引。
        static void syncBoundaryTags(SF::Field& field,
                                     const std::unordered_map<int, std::vector<std::string>>& tagMap,
                                     int ri, int rj, int rk,
                                     int gi, int gj, int gk);

        /// @brief 线性外推虚胞坐标: P_ghost = 2*P_near - P_far。
        /// @param field 网格场。
        /// @param i, j, k 虚胞的(i,j,k)全局索引。
        /// @param ni, nj, nk 最近内点的(i,j,k)索引。
        /// @param fi, fj, fk 次近内点的(i,j,k)索引。
        static void extrapolateCoord(SF::Field& field,
                                     int i, int j, int k,
                                     int ni, int nj, int nk,
                                     int fi, int fj, int fk);

    private:
        /// @brief 计算全场度规(Jacobian和度量系数)。
        /// @param field 已填充物理坐标的Field。
        static void computeMetrics(SF::Field& field);

        /// @brief 将边界层的度规用最近邻内点值填充。
        /// @param field 网格场(原地修改)。
        static void applyBoundaryMetrics(SF::Field& field);

        /// @brief 拷贝度规从一个网格点到另一个。
        /// @param field 网格场。
        /// @param di, dj, dk 目标点索引。
        /// @param si, sj, sk 源点索引。
        static void copyMetricCell(SF::Field& field,
                                   int di, int dj, int dk,
                                   int si, int sj, int sk);

        /// @brief 生成所有虚胞层坐标(线性外推)。
        /// @param field 网格场(原地修改)。
        static void generateGhostCells(SF::Field& field);

        /// @brief 用最近实胞保守量初始化虚胞状态。
        /// @param field 网格场(原地修改)。
        static void initializeGhostStateFromInterior(SF::Field& field);
    };

} // namespace SF
