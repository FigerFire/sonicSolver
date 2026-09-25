/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_generatedMeshWriter.h
/// @brief createMesh生成结果的SFM/VTS/VTM写出工具。
///
/// 本文件只负责把已经生成的结构网格坐标和轻量cell-face拓扑写到磁盘，
/// 不推断MPI接口，也不读取case配置。

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace SF {

namespace GeneratedMeshWriter {

/// @brief 写出一个结构网格块的VTS文件。
void writeBlockVTK(const std::string& filename,
                   const std::vector<double>& allX,
                   const std::vector<double>& allY,
                   const std::vector<double>& allZ,
                   size_t offset,
                   int nx, int ny, int nz);

/// @brief 将一个或多个结构zone写入单个自包含SFM文件。
///
/// 每个zone依次写出#Information、#Point、#Cell、#Face、#Patch和#PointSet。
/// zone之间允许共享物理坐标，但始终保留各自的结构索引；加载后由网格拓扑层
/// 建立重合点映射。
void writeMultiZoneSFM(
    const std::string& filename,
    const std::vector<double>& allX,
    const std::vector<double>& allY,
    const std::vector<double>& allZ,
    const std::vector<int>& blockSizes,
    int nGhost,
    const std::vector<std::map<std::string, std::vector<int>>>& facePatchSets,
    const std::vector<std::map<std::string, std::vector<int>>>& pointSets);

/// @brief 写出VTM多块索引文件。
void writeVTM(const std::string& filename,
              const std::vector<std::string>& blockFiles);

} // namespace GeneratedMeshWriter
} // namespace SF
