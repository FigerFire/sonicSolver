#pragma once

/// @file SF_resultWriter.h
/// @brief 求解结果输出会话，独立于 case 解析器。

#include "SF_valueTypes.h"

#include <functional>
#include <string>
#include <vector>

namespace SF {

class CellFaceMesh;
class Field;
namespace FDM { class IParallelCoordinator; }

/// @brief ResultWriter 构造所需的只读输出配置快照。
struct ResultWriterConfig {
    std::string caseDir;
    std::string caseName;
    std::string jobName;
    std::string outputDir = "result";
    std::vector<BCSetting<double>> densityBoundary;
    std::vector<BCSetting<Vector3>> velocityBoundary;
    std::vector<BCSetting<double>> pressureBoundary;
    std::vector<BCSetting<double>> temperatureBoundary;
    bool ibmOutputEnabled = false;
};

/// @brief VTK/VTU/VTM/PVD 求解结果输出器。
///
/// ResultWriter 不解析 case，也不修改 solver 配置。对象拥有一次运行的输出状态，
/// 包括分区布局和 PVD 时间序列。
class ResultWriter {
public:
    struct PieceLayout {
        int blockId = -1;
        int ownerRank = -1;
        int partitionId = -1;
    };

    /// @brief 额外 VTK 点标量的非拥有视图。
    struct ScalarField {
        std::string name;
        std::function<double(int, int, int)> valueAt;
        bool writesTemperatureBoundaryId = false;
        int components = 1;
        std::function<double(int,int,int,int)> componentAt;
        double value(int i,int j,int k,int component) const {
            return componentAt ? componentAt(i,j,k,component) : valueAt(i,j,k);
        }
    };

    /// @brief 一次并行输出中的结构 patch 视图。
    struct FieldPiece {
        const Field* field = nullptr;
        const std::vector<int>* globalPointIds = nullptr;
        const CellFaceMesh* topology = nullptr;
        int blockId = -1;
        std::vector<ScalarField> extraScalars;
    };

    explicit ResultWriter(ResultWriterConfig config);
    /// @brief 模型提供输出数组；writer 只消费名称、组件数及数值视图。
    using FieldProvider = std::function<std::vector<ScalarField>(
        const Field&,const std::vector<ScalarField>&)>;
    void setFieldProvider(FieldProvider provider) { fieldProvider_=std::move(provider); }


    void setParallelCoordinator(FDM::IParallelCoordinator* parallel) {
        parallel_ = parallel;
    }
    void configurePieces(std::vector<PieceLayout> pieces);
    void clearPieces();

    void save(const Field& field, int step);
    void save(const Field& field, double time);
    void save(const Field& field, int step, int blockId);
    void save(const Field& field, double time, int blockId);
    void savePieces(const std::vector<FieldPiece>& pieces, int step);
    void savePieces(const std::vector<FieldPiece>& pieces, double time);
    void save(const Field& field, int step,
              const std::vector<ScalarField>& extraScalars);
    void save(const Field& field, double time,
              const std::vector<ScalarField>& extraScalars);
    void save(const Field& field, int step, int blockId,
              const std::vector<ScalarField>& extraScalars);
    void save(const Field& field, double time, int blockId,
              const std::vector<ScalarField>& extraScalars);

private:
    struct PVDEntry {
        double timestep = 0.0;
        std::string file;
    };

    ResultWriterConfig config_;
    FieldProvider fieldProvider_;
    std::vector<PieceLayout> pieces_;
    FDM::IParallelCoordinator* parallel_ = nullptr;
    std::vector<PVDEntry> pvdEntries_;

    bool useMultiBlockOutput() const;
    void savePartitionPieces(const std::vector<FieldPiece>& pieces,
                             const std::string& base,
                             double timestep);
    bool writePartitionUnstructuredVTK(
        const std::string& fileName,
        const std::vector<FieldPiece>& pieces);
    bool writePartitionPatchVTK(
        const std::string& fileName,
        const std::vector<FieldPiece>& pieces);
    bool writePartitionIndex(
        const std::string& fileName,
        const std::vector<std::string>& partitionFileNames,
        const std::vector<std::string>& patchFileNames);
    void removeLegacyPatchOutputs(const std::string& base) const;
    void writeMultiBlockIndex(
        const std::string& fileName,
        const std::vector<std::string>& blockFileNames);
    void writeToDisk(const Field& field,
                     const std::string& fileName,
                     bool announce = true,
                     const std::vector<ScalarField>& extraScalars = {});
    void recordAndWritePVD(double timestep,
                           const std::string& dataFileName);
    void flushPVD() const;
};

} // namespace SF
