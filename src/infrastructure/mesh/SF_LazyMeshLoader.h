/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#pragma once

/// @file SF_LazyMeshLoader.h
/// @brief 懒加载网格管理器,类似OpenFOAM的IOobject模式。
///
/// 核心设计理念:
/// - 不一次性将所有网格文件读入内存
/// - 每个SFM文件按需加载: 首次访问时读取并缓存,后续直接返回缓存
/// - 支持分块网格的两种sets索引模式无缝切换
/// - 提供统一的Field构建接口,隐藏SFM解析细节
///
/// 使用示例:
/// @code
///   LazyMeshLoader loader("caseDir");
///   loader.registerMeshFile("mesh.sfm");
///
///   // 按需加载文件中的某个结构zone
///   if (loader.loadBlock(0)) {
///       auto& block = loader.getBlock(0);
///       // 使用block中的坐标和sets...
///   }
///
///   // 构建Field
///   Field field;
///   loader.buildField(0, field);  // 单块
///   loader.buildMergedField(field); // 合并多块
/// @endcode

#include "infrastructure/io/mesh/SF_sfmMesh.h"
#include "SF_meshConfig.h"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <utility>

// 前向声明
namespace SF {
    class Field;
    class MultiBlockMesh;
}

namespace SF {

/// @brief 网格块的懒加载状态。
enum class MeshBlockState {
    /// @brief 未加载,文件路径已注册。
    Unloaded,
    /// @brief 已加载到内存。
    Loaded,
    /// @brief 加载失败。
    Failed
};

/// @brief 单个网格块的懒加载包装。
struct LazyMeshBlock {
    /// @brief SFM文件路径(相对或绝对)。
    std::string filePath;
    /// @brief 加载状态。
    MeshBlockState state = MeshBlockState::Unloaded;
    /// @brief 解析后的网格数据(加载后有效)。
    SFMMeshBlock data;

    /// @brief 块在全局网格中的起始索引偏移(用于模式A)。
    int globalIndexOffset = 0;
    /// @brief 块在全局网格中的维度起始(x方向)。
    int globalStartI = 0;
    /// @brief 块在全局网格中的维度起始(y方向)。
    int globalStartJ = 0;
    /// @brief 块在全局网格中的维度起始(z方向)。
    int globalStartK = 0;
};

/// @brief 懒加载网格管理器。
///
/// 管理多个SFM网格块的注册、按需加载、缓存和Field构建。
/// 兼容两种sets索引模式:
/// - 模式A(每块独立sets): 每块的sets索引相对于该块的局部网格
/// - 模式B(共享sets): 所有块共享一个_sets.sfm文件,索引相对于全局网格
class LazyMeshLoader {
public:
    /// @brief 构造时指定工作目录。
    /// @param caseDir 算例目录,用于解析相对路径。
    explicit LazyMeshLoader(
            const std::string& caseDir = ".",
            MeshRuntimeConfig config = {})
        : caseDir_(caseDir), meshIO_(caseDir), config_(std::move(config)) {}

    /// @brief 设置算例目录。
    /// @param dir 算例目录路径。
    void setCaseDir(const std::string& dir) {
        caseDir_ = dir;
        meshIO_.setCaseDir(dir);
    }

    // ============================================================
    //  文件注册
    // ============================================================

    /// @brief 注册一个网格SFM文件(不会立即读取)。
    /// @param filePath 相对于caseDir的文件路径。
    void registerMeshFile(const std::string& filePath);

    /// @brief 注册共享sets文件(模式B)。
    /// @param filePath 共享sets文件路径。
    void registerSetsFile(const std::string& filePath);

    /// @brief 批量注册网格文件列表。
    /// @param filePaths 文件路径列表。
    void registerMeshFiles(const std::vector<std::string>& filePaths);

    // ============================================================
    //  按需加载 (Lazy Loading)
    // ============================================================

    /// @brief 按需加载指定序号的网格块。
    /// @param blockIndex 块序号(0-based)。
    /// @return 加载成功返回true; 已加载则直接返回true。
    bool loadBlock(int blockIndex);

    /// @brief 加载所有已注册的网格块。
    /// @return 全部加载成功返回true。
    bool loadAllBlocks();

    /// @brief 仅加载共享sets文件(模式B)。
    /// @return 成功返回true。
    bool loadSharedSets();

    /// @brief 卸载指定块(释放内存)。
    /// @param blockIndex 块序号。
    void unloadBlock(int blockIndex);

    /// @brief 卸载所有块。
    void unloadAllBlocks();

    // ============================================================
    //  数据访问
    // ============================================================

    /// @brief 获取已注册的块数量。
    /// @return 块数量。
    int blockCount() const { return static_cast<int>(blocks_.size()); }

    /// @brief 获取指定块的原始数据(必须先加载)。
    /// @param blockIndex 块序号。
    /// @return 块的const引用。
    const SFMMeshBlock& getBlock(int blockIndex) const;

    /// @brief 获取指定块的可写引用。
    /// @param blockIndex 块序号。
    /// @return 块的可写引用。
    SFMMeshBlock& getBlockWritable(int blockIndex);

    /// @brief 获取共享sets(模式B)。
    /// @return sets映射的const引用; 模式A时为空。
    const std::map<std::string, std::vector<int>>& sharedSets() const {
        return sharedSets_;
    }

    /// @brief 判断当前是哪种sets模式。
    /// @return true=模式B(共享sets), false=模式A(每块独立sets)。
    bool isSharedSetsMode() const { return !sharedSetsFile_.empty(); }

    /// @brief 获取全局网格总维度(仅模式B有意义)。
    /// @param nx, ny, nz 输出的全局维度。
    /// @param ng 输出的虚胞层数。
    /// @return 成功获取返回true。
    bool getGlobalDimensions(int& nx, int& ny, int& nz, int& ng) const;

    // ============================================================
    //  Field构建
    // ============================================================

    /// @brief 从单个块构建Field(包含度规计算和虚胞生成)。
    /// @param blockIndex 块序号。
    /// @param field 输出的Field对象。
    /// @return 成功返回true。
    bool buildField(int blockIndex, Field& field);

    /// @brief 构建合并后的单一Field(用于单块求解器)。
    /// @param field 输出的Field对象。
    /// @return 成功返回true。
    bool buildMergedField(Field& field);

    /// @brief 获取块的局部sets(已转换为Field内部索引)。
    /// @param blockIndex 块序号。
    /// @return 局部sets映射; 模式B时自动从共享sets中筛取属于本块的索引。
    std::map<std::string, std::vector<int>>
    getLocalSets(int blockIndex) const;

    // ============================================================
    //  工具函数
    // ============================================================

    /// @brief 获取块在全局网格中的起始索引偏移。
    /// @param blockIndex 块序号。
    /// @return 0-based全局偏移。
    int globalOffset(int blockIndex) const;

    /// @brief 计算块的局部索引在全局网格中的对应索引。
    /// @param blockIndex 块序号。
    /// @param localIdx 局部0-based索引。
    /// @return 全局0-based索引; 无偏移时为localIdx本身。
    int localToGlobalIndex(int blockIndex, int localIdx) const;

    /// @brief 判断某全局索引属于哪个块。
    /// @param globalIdx 全局0-based索引。
    /// @param outBlock 输出的块序号。
    /// @param outLocalIdx 输出的块内局部索引。
    /// @return 找到所属块返回true。
    bool findBlockForGlobalIndex(int globalIdx,
                                 int& outBlock,
                                 int& outLocalIdx) const;

private:
    MeshRuntimeConfig config_;
    std::string caseDir_;
    SFMMeshIO meshIO_;
    std::vector<LazyMeshBlock> blocks_;
    std::string sharedSetsFile_;
    std::map<std::string, std::vector<int>> sharedSets_;
    bool sharedSetsLoaded_ = false;
    int globalNx_ = 0, globalNy_ = 0, globalNz_ = 0, globalNg_ = 0;

    /// @brief 计算块间全局偏移(自动检测并重排)。
    void computeGlobalOffsets();

    /// @brief 从局部坐标构建Field的辅助函数。
    /// @param block 网格块数据。
    /// @param field 输出的Field。
    /// @param sets 边界集合映射(已转为Field内部索引)。
    /// @param label 日志标签。
    /// @return 成功返回true。
    bool buildFieldFromBlock(const SFMMeshBlock& block,
                             Field& field,
                             const std::map<std::string, std::vector<int>>& sets,
                             const std::string& label);
};

} // namespace SF
