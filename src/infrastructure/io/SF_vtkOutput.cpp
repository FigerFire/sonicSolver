/// @file SF_vtkOutput.cpp
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*----------IO private implementation---------*/

#include "SF_resultWriter.h"
#include "core/interfaces/SF_log.h"
#include "private/SF_vtkUtils.h"
#include "core/interfaces/SF_parallelCoordinatorContract.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

namespace SF {

using namespace IOPrivate;

ResultWriter::ResultWriter(ResultWriterConfig config)
    : config_(std::move(config)) {}

/// @brief 按步数格式化文件名 (不含扩展名)。
static std::string formatStepName(const std::string& prefix, int step) {
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(6) << std::setfill('0') << step;
    return oss.str();
}

/// @brief 按时间格式化文件名 (不含扩展名)。
static std::string formatTimeName(const std::string& prefix, double time) {
    std::ostringstream oss;
    oss << prefix << "_t" << std::fixed << std::setprecision(6) << time;
    std::string s = oss.str();
    // 去掉末尾无意义的零
    if (s.find('.') != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    return s;
}

void ResultWriter::configurePieces(std::vector<PieceLayout> pieces) {
    std::sort(pieces.begin(), pieces.end(),
              [](const PieceLayout& a, const PieceLayout& b) {
                  return a.blockId < b.blockId;
              });
    pieces_ = std::move(pieces);
}

void ResultWriter::clearPieces() {
    pieces_.clear();
}

bool ResultWriter::useMultiBlockOutput() const {
    return parallel_ != nullptr && parallel_->active() && !pieces_.empty();
}

void ResultWriter::save(const SF::Field& field, int step) {
    save(field, step, std::vector<ScalarField>{});
}

void ResultWriter::save(const SF::Field& field,
                 int step,
                 const std::vector<ScalarField>& extraScalars) {
    std::string prefix = outputPrefix(config_.jobName, config_.caseName);
    std::string base   = formatStepName(prefix, step);
    std::string fullPath = config_.caseDir + "/" + config_.outputDir + "/" + base + ".vts";
    writeToDisk(field, fullPath, true, extraScalars);
    recordAndWritePVD((double)step, base + ".vts");
}

void ResultWriter::save(const SF::Field& field, int step, int blockId) {
    save(field, step, blockId, std::vector<ScalarField>{});
}

void ResultWriter::save(const SF::Field& field,
                 int step,
                 int blockId,
                 const std::vector<ScalarField>& extraScalars) {
    // HYPRE 会在单进程压力基计算中初始化 MPI。此时网格仍是单块串行网格，
    // 不应进入遗留的 _block000 输出路径，否则 PVD 时间序列无法注册。
    if (!parallel_ || parallel_->size() == 1) {
        save(field, step, extraScalars);
        return;
    }

    std::string prefix = outputPrefix(config_.jobName, config_.caseName);
    std::string base   = formatStepName(prefix, step);

    if (useMultiBlockOutput()) {
        std::string vtmPath = config_.caseDir + "/" + config_.outputDir + "/" + base + ".vtm";

        std::vector<std::string> blockFiles;
        blockFiles.reserve(pieces_.size());
        for (const auto& piece : pieces_) {
            std::ostringstream blockName;
            blockName << base << "_block"
                      << std::setw(3) << std::setfill('0') << piece.blockId
                      << ".vts";
            std::string blockPath = config_.caseDir + "/" + config_.outputDir + "/" + blockName.str();
            blockFiles.push_back(blockPath);
        }

        if (blockId < 0 || blockId >= (int)blockFiles.size()) {
            broadcast("Error: ", "Invalid VTK block id " + std::to_string(blockId));
            return;
        }
        const auto& piece = pieces_[(size_t)blockId];
        const bool shouldWrite =
            !parallel_ || piece.ownerRank < 0 || parallel_->rank() == piece.ownerRank;
        if (shouldWrite) {
            writeToDisk(field, blockFiles[(size_t)blockId], false, extraScalars);
        }
        if (parallel_) parallel_->barrier();
        if (!parallel_ || parallel_->rank() == 0) {
            writeMultiBlockIndex(vtmPath, blockFiles);
            recordAndWritePVD((double)step, base + ".vtm");
        }
        if (parallel_) parallel_->barrier();
        return;
    }

    // 非多块分块输出: 仍然把 blockId 编码进文件名 (兼容遗留用法)
    std::ostringstream oss;
    oss << config_.caseDir << "/" << config_.outputDir << "/" << base
        << "_block" << std::setw(3) << std::setfill('0') << blockId
        << ".vts";
    writeToDisk(field, oss.str(), true, extraScalars);
}

void ResultWriter::save(const SF::Field& field, double time) {
    save(field, time, std::vector<ScalarField>{});
}

void ResultWriter::save(const SF::Field& field,
                 double time,
                 const std::vector<ScalarField>& extraScalars) {
    std::string prefix   = outputPrefix(config_.jobName, config_.caseName);
    std::string base     = formatTimeName(prefix, time);
    std::string fullPath = config_.caseDir + "/" + config_.outputDir + "/" + base + ".vts";
    writeToDisk(field, fullPath, true, extraScalars);
    recordAndWritePVD(time, base + ".vts");
}

void ResultWriter::save(const SF::Field& field, double time, int blockId) {
    save(field, time, blockId, std::vector<ScalarField>{});
}

void ResultWriter::save(const SF::Field& field,
                 double time,
                 int blockId,
                 const std::vector<ScalarField>& extraScalars) {
    // 单进程 HYPRE 只借用 MPI 线性代数后端，不改变 VTK 的串行拓扑语义。
    if (!parallel_ || parallel_->size() == 1) {
        save(field, time, extraScalars);
        return;
    }

    std::string prefix = outputPrefix(config_.jobName, config_.caseName);
    std::string base   = formatTimeName(prefix, time);

    if (useMultiBlockOutput()) {
        std::string vtmPath = config_.caseDir + "/" + config_.outputDir + "/" + base + ".vtm";

        std::vector<std::string> blockFiles;
        blockFiles.reserve(pieces_.size());
        for (const auto& piece : pieces_) {
            std::ostringstream blockName;
            blockName << base << "_block"
                      << std::setw(3) << std::setfill('0') << piece.blockId
                      << ".vts";
            std::string blockPath = config_.caseDir + "/" + config_.outputDir + "/" + blockName.str();
            blockFiles.push_back(blockPath);
        }

        if (blockId < 0 || blockId >= (int)blockFiles.size()) {
            broadcast("Error: ", "Invalid VTK block id " + std::to_string(blockId));
            return;
        }
        const auto& piece = pieces_[(size_t)blockId];
        const bool shouldWrite =
            !parallel_ || piece.ownerRank < 0 || parallel_->rank() == piece.ownerRank;
        if (shouldWrite) {
            writeToDisk(field, blockFiles[(size_t)blockId], false, extraScalars);
        }
        if (parallel_) parallel_->barrier();
        if (!parallel_ || parallel_->rank() == 0) {
            writeMultiBlockIndex(vtmPath, blockFiles);
            recordAndWritePVD(time, base + ".vtm");
        }
        if (parallel_) parallel_->barrier();
        return;
    }

    // 非多块分块输出
    std::ostringstream oss;
    oss << config_.caseDir << "/" << config_.outputDir << "/" << base
        << "_block" << std::setw(3) << std::setfill('0') << blockId
        << ".vts";
    writeToDisk(field, oss.str(), true, extraScalars);
}

void ResultWriter::savePieces(const std::vector<FieldPiece>& pieces, int step) {
    const std::string prefix = outputPrefix(config_.jobName, config_.caseName);
    const std::string base = formatStepName(prefix, step);
    savePartitionPieces(pieces, base, (double)step);
}

void ResultWriter::savePieces(const std::vector<FieldPiece>& pieces, double time) {
    const std::string prefix = outputPrefix(config_.jobName, config_.caseName);
    const std::string base = formatTimeName(prefix, time);
    savePartitionPieces(pieces, base, time);
}

void ResultWriter::savePartitionPieces(
    const std::vector<FieldPiece>& pieces,
    const std::string& base,
    double timestep) {
    const std::string vtmPath =
        config_.caseDir + "/" + config_.outputDir + "/" + base + ".vtm";

    int partitionCount = 0;
    for (const PieceLayout& layout : pieces_) {
        partitionCount = std::max(partitionCount, layout.partitionId + 1);
    }
    std::vector<std::vector<FieldPiece>> partitionPieces(
        (size_t)partitionCount);

    bool localLayoutOk = partitionCount > 0;
    if (parallel_ && partitionCount != parallel_->size()) {
        localLayoutOk = false;
    }
    for (const FieldPiece& piece : pieces) {
        if (!piece.field || !piece.globalPointIds || piece.blockId < 0 ||
            piece.blockId >= (int)pieces_.size()) {
            localLayoutOk = false;
            continue;
        }
        const PieceLayout& layout = pieces_[(size_t)piece.blockId];
        if (layout.partitionId < 0 ||
            layout.partitionId >= partitionCount) {
            localLayoutOk = false;
            continue;
        }
        partitionPieces[(size_t)layout.partitionId].push_back(piece);
    }
    for (const auto& partition : partitionPieces) {
        if (partition.empty()) localLayoutOk = false;
    }

    const bool layoutOk =
        parallel_ ? parallel_->allRanksAgree(localLayoutOk) : localLayoutOk;
    if (!layoutOk) {
        broadcast("Fatal: ",
                  "Invalid or incomplete VTK partition layout.");
        return;
    }

    std::vector<std::string> partitionFiles;
    std::vector<std::string> patchFiles;
    partitionFiles.reserve((size_t)partitionCount);
    patchFiles.reserve((size_t)partitionCount);
    bool localWriteOk = true;
    for (int partitionId = 0; partitionId < partitionCount; ++partitionId) {
        std::ostringstream fileName;
        fileName << base << "_partition"
                 << std::setw(3) << std::setfill('0') << partitionId
                 << ".vtu";
        partitionFiles.push_back(
            config_.caseDir + "/" + config_.outputDir + "/" + fileName.str());

        std::ostringstream patchFileName;
        patchFileName << base << "_partition"
                      << std::setw(3) << std::setfill('0') << partitionId
                      << "_patches.vtu";
        patchFiles.push_back(
            config_.caseDir + "/" + config_.outputDir + "/" + patchFileName.str());

        if (!parallel_ || parallel_->rank() == partitionId) {
            if (!writePartitionUnstructuredVTK(
                    partitionFiles.back(),
                    partitionPieces[(size_t)partitionId])) {
                localWriteOk = false;
            }
            if (!writePartitionPatchVTK(
                    patchFiles.back(),
                    partitionPieces[(size_t)partitionId])) {
                localWriteOk = false;
            }
        }
    }

    const bool writeOk =
        parallel_ ? parallel_->allRanksAgree(localWriteOk) : localWriteOk;
    if (!writeOk) {
        broadcast("Fatal: ", "Failed to write one or more VTK partitions.");
        return;
    }

    bool localIndexOk = true;
    if (!parallel_ || parallel_->rank() == 0) {
        removeLegacyPatchOutputs(base);
        localIndexOk = writePartitionIndex(vtmPath, partitionFiles, patchFiles);
        if (localIndexOk) {
            recordAndWritePVD(timestep, base + ".vtm");
        }
    }
    const bool indexOk =
        parallel_ ? parallel_->allRanksAgree(localIndexOk) : localIndexOk;
    if (!indexOk) {
        broadcast("Fatal: ", "Failed to write VTK partition index.");
    }
}

} // namespace SF
