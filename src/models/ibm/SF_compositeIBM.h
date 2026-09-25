/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.14-----------*/

#pragma once

/// @file SF_compositeIBM.h
/// @brief 复合结构网格的分解前IBM分类与分块后ghost闭合。

#include "SF_MultiBlockMesh.h"
#include "SF_ibmConfig.h"
#include "SF_ibmTopology.h"
#include "geoProcessing/SF_STLGeometry.h"
#include "method/ghost/SF_ghost.h"

#include <string>
#include <vector>

namespace SF {
namespace IBM {

/// @brief 管理复合网格IBM几何分类和每个结构patch的局部闭合数据。
/*
    原始 multi-block zones
                │
                ↓
    preprocessSourceZones()
    分解前、全局统一处理
                │
                ↓
        MPI partition
                │
                ↓
    setupLocalPatches()
    分解后、本 rank 初始化
                │
                ↓
    time marching
    applyLocalPatches()
*/
class CompositeIB {
public:
    /// @brief 设置由 application 层冻结的 IBM/ILW 运行配置。
    void configure(const IBMRuntimeConfig& config) { config_ = config; }

    /// @brief 在MPI分解前按global point统一完成源zone IBM几何分类。
    /// @param zones 已完成复合拓扑装配的源zone。
    /// @param caseDir 算例目录。
    /// @param stlFiles IBM STL文件列表。
    /// @return 分类成功返回true。
    bool preprocessSourceZones(std::vector<RawMeshBlock>& zones,
                               const std::string& caseDir,
                               const std::vector<std::string>& stlFiles);

    /// @brief 为分解产生的MPI halo补齐IBM几何并建立本rank patch闭合数据。
    /// @param mesh 已完成分解、Field构建和halo计划构建的复合网格。
    /// @param localPatchIds 当前rank实际推进的patch编号。
    /// @param mpiRank 当前MPI rank，用于前处理失败诊断。
    /// @return 所有本地patch完成闭合预处理时返回true。
    bool setupLocalPatches(MultiBlockMesh& mesh,
                           const std::vector<int>& localPatchIds,
                           int mpiRank);

    /// @brief 重建当前rank全部本地patch的IBM ghost状态。
    /// @param mesh 复合网格。
    /// @param localPatchIds 当前rank实际推进的patch编号。
    void applyLocalPatches(MultiBlockMesh& mesh,
                           const std::vector<int>& localPatchIds,
                           double time,
                           double dt);

    /// @brief 查询复合IBM是否已成功加载。
    bool active() const { return active_; }

private:
    bool active_ = false;
    GeoProcessing::STLGeometry geometry_;
    std::vector<GhostIBM::WeightBuilder> patchWeights_;
    std::vector<unsigned char> patchReady_;
    GhostIBM::GhostCellIBM ghostCell_;
    IBMRuntimeConfig config_;
    /// @brief 每个结构化 patch 一份的 IBM 拓扑几何（ghost 层、壁面/镜像点、法向）。
    std::vector<IBM::IBMGeometry> blockGeometries_;

    /// @brief 从分解后 patch 携带的 IBM POD 分类重建 per-patch 拓扑几何与分类。
    bool buildBlockGeometries(MultiBlockMesh& mesh);
    bool classifyPartitionHalos(MultiBlockMesh& mesh);
};

} // namespace IBM
} // namespace SF
