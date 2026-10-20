/// @file SF_LazyMeshLoader.cpp
/// @brief 结构/多块网格拓扑、度量与加载实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#include "SF_LazyMeshLoader.h"
#include "core/field/SF_field.h"
#include "SF_mesh.h"

#include <iostream>
#include <algorithm>

namespace SF {

// ============================================================
//  文件注册
// ============================================================

void LazyMeshLoader::registerMeshFile(const std::string& filePath) {
    LazyMeshBlock block;
    block.filePath = filePath;
    block.state = MeshBlockState::Unloaded;
    blocks_.push_back(std::move(block));
}

void LazyMeshLoader::registerSetsFile(const std::string& filePath) {
    sharedSetsFile_ = filePath;
    sharedSetsLoaded_ = false;
    sharedSets_.clear();
}

void LazyMeshLoader::registerMeshFiles(const std::vector<std::string>& filePaths) {
    for (const auto& f : filePaths) {
        registerMeshFile(f);
    }
}

// ============================================================
//  按需加载
// ============================================================

bool LazyMeshLoader::loadBlock(int blockIndex) {
    if (blockIndex < 0 || blockIndex >= static_cast<int>(blocks_.size())) {
        std::cerr << "[LazyMesh] 无效块序号: " << blockIndex << std::endl;
        return false;
    }

    auto& block = blocks_[blockIndex];

    // 已加载则直接返回
    if (block.state == MeshBlockState::Loaded) return true;

    // 尝试读取
    std::cout << "[LazyMesh] 加载块 " << blockIndex
              << ": " << block.filePath << std::endl;

    if (meshIO_.readSingleBlock(block.filePath, block.data)) {
        block.state = MeshBlockState::Loaded;
        // 如果是模式A(每块独立sets), 计算全局偏移
        if (!isSharedSetsMode()) {
            computeGlobalOffsets();
        }
        return true;
    }

    block.state = MeshBlockState::Failed;
    std::cerr << "[LazyMesh] 块 " << blockIndex << " 加载失败!" << std::endl;
    return false;
}

bool LazyMeshLoader::loadAllBlocks() {
    bool allOk = true;
    for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
        if (!loadBlock(b)) allOk = false;
    }
    return allOk;
}

bool LazyMeshLoader::loadSharedSets() {
    if (sharedSetsFile_.empty()) return false;
    if (sharedSetsLoaded_) return true;

    std::cout << "[LazyMesh] 加载共享sets: " << sharedSetsFile_ << std::endl;

    if (!meshIO_.readSetsOnly(sharedSetsFile_, sharedSets_)) {
        std::cerr << "[LazyMesh] 共享sets加载失败!" << std::endl;
        return false;
    }

    // 读取共享sets文件的维度信息
    if (!meshIO_.readDimensions(sharedSetsFile_,
                                globalNx_, globalNy_, globalNz_, globalNg_)) {
        std::cerr << "[LazyMesh] 无法读取共享sets的维度信息" << std::endl;
        return false;
    }

    sharedSetsLoaded_ = true;
    return true;
}

void LazyMeshLoader::unloadBlock(int blockIndex) {
    if (blockIndex < 0 || blockIndex >= static_cast<int>(blocks_.size())) return;
    auto& block = blocks_[blockIndex];
    block.data = SFMMeshBlock();
    block.state = MeshBlockState::Unloaded;
}

void LazyMeshLoader::unloadAllBlocks() {
    for (auto& block : blocks_) {
        block.data = SFMMeshBlock();
        block.state = MeshBlockState::Unloaded;
    }
    sharedSets_.clear();
    sharedSetsLoaded_ = false;
}

// ============================================================
//  数据访问
// ============================================================

const SFMMeshBlock& LazyMeshLoader::getBlock(int blockIndex) const {
    static SFMMeshBlock emptyBlock;
    if (blockIndex < 0 || blockIndex >= static_cast<int>(blocks_.size()))
        return emptyBlock;
    return blocks_[blockIndex].data;
}

SFMMeshBlock& LazyMeshLoader::getBlockWritable(int blockIndex) {
    static SFMMeshBlock emptyBlock;
    if (blockIndex < 0 || blockIndex >= static_cast<int>(blocks_.size()))
        return emptyBlock;
    return blocks_[blockIndex].data;
}

bool LazyMeshLoader::getGlobalDimensions(int& nx, int& ny, int& nz, int& ng) const {
    if (isSharedSetsMode()) {
        nx = globalNx_; ny = globalNy_; nz = globalNz_; ng = globalNg_;
        return true;
    }

    // 模式A: 需要先加载所有块才能计算全局维度
    nx = ny = nz = ng = 0;
    return false;
}

// ============================================================
//  全局偏移计算
// ============================================================

void LazyMeshLoader::computeGlobalOffsets() {
    int offset = 0;
    for (size_t b = 0; b < blocks_.size(); ++b) {
        if (blocks_[b].state != MeshBlockState::Loaded) continue;
        const auto& data = blocks_[b].data;
        blocks_[b].globalIndexOffset = offset;
        blocks_[b].globalStartI = (b > 0) ? offset % data.nx : 0;
        blocks_[b].globalStartJ = 0;
        blocks_[b].globalStartK = 0;
        offset += data.cellCount();
    }
}

int LazyMeshLoader::globalOffset(int blockIndex) const {
    if (blockIndex < 0 || blockIndex >= static_cast<int>(blocks_.size())) return 0;
    return blocks_[blockIndex].globalIndexOffset;
}

int LazyMeshLoader::localToGlobalIndex(int blockIndex, int localIdx) const {
    return globalOffset(blockIndex) + localIdx;
}

bool LazyMeshLoader::findBlockForGlobalIndex(int globalIdx,
                                             int& outBlock,
                                             int& outLocalIdx) const {
    for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
        int offset = globalOffset(b);
        int cellCount = blocks_[b].data.cellCount();
        if (globalIdx >= offset && globalIdx < offset + cellCount) {
            outBlock = b;
            outLocalIdx = globalIdx - offset;
            return true;
        }
    }
    return false;
}

// ============================================================
//  局部sets获取
// ============================================================

std::map<std::string, std::vector<int>>
LazyMeshLoader::getLocalSets(int blockIndex) const {
    std::map<std::string, std::vector<int>> result;

    if (isSharedSetsMode()) {
        // 模式B: 从共享sets中筛取属于本块的索引
        int offset = globalOffset(blockIndex);
        int cellCount = blocks_[blockIndex].data.cellCount();

        for (const auto& [name, globalIndices] : sharedSets_) {
            std::vector<int> local;
            for (int gi : globalIndices) {
                if (gi >= offset && gi < offset + cellCount) {
                    local.push_back(gi - offset);
                }
            }
            if (!local.empty()) result[name] = std::move(local);
        }
    } else {
        // 模式A: 直接返回块的局部sets
        result = blocks_[blockIndex].data.pointSets;
    }

    return result;
}

// ============================================================
//  Field构建
// ============================================================

bool LazyMeshLoader::buildFieldFromBlock(
    const SFMMeshBlock& block,
    Field& field,
    const std::map<std::string, std::vector<int>>& sets,
    const std::string& label)
{
    // 借用现有的 Mesh::setupFieldFromRaw
    return Mesh::setupFieldFromRaw(
        field,
        block.nx, block.ny, block.nz, block.ng,
        block.x, block.y, block.z,
        sets,
        config_,
        label);
}

bool LazyMeshLoader::buildField(int blockIndex, Field& field) {
    if (!loadBlock(blockIndex)) return false;

    const auto& block = getBlock(blockIndex);
    if (!block.isValid()) return false;

    auto localSets = getLocalSets(blockIndex);

    return buildFieldFromBlock(block, field, localSets,
                               "block " + std::to_string(blockIndex));
}

bool LazyMeshLoader::buildMergedField(Field& field) {
    if (blocks_.empty()) {
        std::cerr << "[LazyMesh] 无注册的网格块" << std::endl;
        return false;
    }

    // 对于单块,直接构建
    if (blocks_.size() == 1) {
        return buildField(0, field);
    }

    // 对于多块,需要合并坐标和sets
    // 这里暂用第一个块的逻辑,后续扩展
    if (!loadBlock(0)) return false;

    const auto& firstBlock = getBlock(0);
    auto localSets = getLocalSets(0);

    return buildFieldFromBlock(firstBlock, field, localSets,
                               "merged mesh");
}

} // namespace SF
