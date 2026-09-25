/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#pragma once

/// @file SF_sfmMesh.h
/// @brief SFM格式网格文件的读写接口。
///
/// SFM (Sonic Fluid Mesh) 是sonicSolver的原生网格格式。
/// 格式说明:
/// @code
///   #Information
///   nx ny nz ng              ← 网格尺寸(不含虚胞)和虚胞层数
///
///   #Point
///   x1 y1 z1                 ← 物理坐标，共(nx+1)*(ny+1)*(nz+1)个点
///   ...
///
///   #SetName1
///   idx1 idx2 ...            ← 边界集合点索引(0-based局部索引)
///   ...
/// @endcode
///
/// 分块网格sets兼容两种模式:
/// - 模式A(独立索引): 每个block的SFM文件各自包含#SetName段,
///   索引范围为[0, block.nx*block.ny*block.nz)
/// - 模式B(全局索引): 仅共享的_sets.sfm文件包含#SetName段,
///   索引范围为[0, total_nx*total_ny*total_nz)
///
/// 读取时自动检测模式,写出时通过参数控制模式。

#include <string>
#include <vector>
#include <map>

namespace SF {

/// @brief 解析后的原始网格块数据(不含虚胞)。
struct SFMMeshBlock {
    /// @brief 网格维度(不含虚胞)。
    int nx = 0, ny = 0, nz = 0;
    /// @brief 虚胞层数。
    int ng = 0;
    /// @brief 物理坐标数组, 长度 = (nx+1)*(ny+1)*(nz+1)。
    std::vector<double> x, y, z;
    /// @brief 边界集合点, key=集合名, value=局部索引列表。
    std::map<std::string, std::vector<int>> pointSets;
    /// @brief 文件名(用于日志)。
    std::string sourceFile;
    /// @brief 是否包含有效点数据。
    bool hasPoints = false;
    /// @brief 是否包含sets数据。
    bool hasSets = false;

    /// @brief 检查数据是否有效。
    /// @return nx>0 && ny>0 && nz>0 && 坐标数量正确时返回true。
    bool isValid() const {
        int expectedCount = (nx + 1) * (ny + 1) * (nz + 1);
        return nx > 0 && ny > 0 && nz > 0
            && static_cast<int>(x.size()) == expectedCount
            && x.size() == y.size()
            && x.size() == z.size();
    }

    /// @brief 获取内部网格点数。
    /// @return nx * ny * nz。
    int cellCount() const { return nx * ny * nz; }
};

/// @brief SFM网格IO管理器。
///
/// 提供懒加载能力:不一次性读取所有网格文件,而是在需要时按文件名读取。
/// 类似OpenFOAM的IOobject模式,用到什么读什么。
class SFMMeshIO {
public:
    /// @brief 构造时设置工作目录。
    /// @param caseDir 算例目录路径, 用于解析相对路径。
    explicit SFMMeshIO(const std::string& caseDir = ".")
        : caseDir_(caseDir) {}

    /// @brief 设置工作目录。
    /// @param caseDir 算例目录路径。
    void setCaseDir(const std::string& caseDir) { caseDir_ = caseDir; }

    /// @brief 获取工作目录。
    /// @return 当前算例目录路径。
    const std::string& caseDir() const { return caseDir_; }

    // ============================================================
    //  SFM读取 (懒加载,按需读取)
    // ============================================================

    /// @brief 读取单个SFM文件(需要时调用,不缓存)。
    /// @param filePath SFM文件路径(支持相对路径或绝对路径)。
    /// @param outBlock 输出的网格块数据。
    /// @return 成功读取且数据有效返回true。
    bool readSingleBlock(const std::string& filePath,
                         SFMMeshBlock& outBlock) const;

    /// @brief 读取SFM文件的#Point段(仅坐标,不含sets)。
    /// @param filePath SFM文件路径。
    /// @param outBlock 输出的网格块(仅填充坐标和维度)。
    /// @return 成功返回true。
    bool readPointsOnly(const std::string& filePath,
                        SFMMeshBlock& outBlock) const;

    /// @brief 读取SFM文件的#SetName段(仅sets,不含坐标)。
    /// @param filePath SFM文件路径, 通常是_sets.sfm。
    /// @param outSets 输出的边界集合映射。
    /// @return 成功读取至少一个集合返回true。
    bool readSetsOnly(const std::string& filePath,
                      std::map<std::string, std::vector<int>>& outSets) const;

    /// @brief 检查SFM文件是否包含sets段。
    /// @param filePath SFM文件路径。
    /// @return 包含至少一个#SetName段返回true。
    bool hasSetsSection(const std::string& filePath) const;

    /// @brief 获取SFM文件的维度信息(不读取坐标,极快)。
    /// @param filePath SFM文件路径。
    /// @param nx, ny, nz 输出的网格维度。
    /// @param ng 输出的虚胞层数。
    /// @return 成功读取返回true。
    bool readDimensions(const std::string& filePath,
                        int& nx, int& ny, int& nz, int& ng) const;

    // ============================================================
    //  SFM写出
    // ============================================================

    /// @brief 写出单个网格块为SFM文件。
    /// @param filePath 输出文件路径。
    /// @param block 网格块数据。
    /// @param writeSets 是否写出sets段; 默认true。
    /// @param useGlobalIndex 是否使用全局索引(跨块)写出sets;
    ///        设为true时, 调用者需事先将sets索引转换为全局索引。
    ///        设为false时, writes为局部索引(仅本块内有效)。
    /// @return 成功写出返回true。
    bool writeSingleBlock(const std::string& filePath,
                          const SFMMeshBlock& block,
                          bool writeSets = true,
                          bool useGlobalIndex = false) const;

    /// @brief 写出纯sets文件(不含坐标,用于模式B的共享sets)。
    /// @param filePath 输出文件路径, 通常命名为*_sets.sfm。
    /// @param sets 边界集合映射, key=集合名, value=全局索引列表。
    /// @param nx, ny, nz 全局网格维度(用于格式校验)。
    /// @param ng 虚胞层数。
    /// @return 成功写出返回true。
    bool writeSetsFile(const std::string& filePath,
                       const std::map<std::string, std::vector<int>>& sets,
                       int nx, int ny, int nz, int ng) const;

    /// @brief 写出多块网格(模式A: 每块各自有sets段)。
    /// @param basePath 基础路径, 自动追加_block000.sfm等后缀。
    /// @param blocks 网格块列表。
    /// @param localSets 每块对应的局部sets映射(与blocks一一对应)。
    /// @return 全部写出成功返回true。
    bool writeMultiBlockLocal(
        const std::string& basePath,
        const std::vector<SFMMeshBlock>& blocks,
        const std::vector<std::map<std::string, std::vector<int>>>& localSets)
        const;

    /// @brief 写出多块网格(模式B: 共享sets文件)。
    /// @param basePath 基础路径。
    /// @param blocks 网格块列表(不再写各自的sets段)。
    /// @param globalSets 全局sets映射, 索引范围为所有块的并集。
    /// @param globalNx, globalNy, globalNz 全局网格总维度。
    /// @param ng 虚胞层数。
    /// @return 全部写出成功返回true。
    bool writeMultiBlockGlobal(
        const std::string& basePath,
        const std::vector<SFMMeshBlock>& blocks,
        const std::map<std::string, std::vector<int>>& globalSets,
        int globalNx, int globalNy, int globalNz, int ng) const;

private:
    std::string caseDir_;

    /// @brief 解析相对路径为绝对路径。
    /// @param filePath 原始路径。
    /// @return 绝对路径字符串。
    std::string resolvePath(const std::string& filePath) const;

    /// @brief 清理行尾的注释和空白。
    /// @param line 待清理的行。
    static void cleanLine(std::string& line);

    /// @brief 读取SFM文件内部实现。
    /// @param filePath 文件路径。
    /// @param readPoints 是否读取#Point段。
    /// @param readSets 是否读取#SetName段。
    /// @param outBlock 输出数据。
    /// @return 成功返回true。
    bool readSFMInternal(const std::string& filePath,
                         bool readPoints, bool readSets,
                         SFMMeshBlock& outBlock) const;
};

} // namespace SF
