/// @file SF_vtkIndexWriter.cpp
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*----------IO private implementation---------*/

#include "SF_resultWriter.h"
#include "core/interfaces/SF_log.h"
#include "private/SF_vtkUtils.h"
#include "SF_cellFaceMesh.h"
#include "core/interfaces/SF_parallelCoordinatorContract.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace SF {

using namespace IOPrivate;

bool ResultWriter::writePartitionPatchVTK(
    const std::string& fileName,
    const std::vector<FieldPiece>& pieces) {
    struct PatchFaceReference {
        const FieldPiece* piece = nullptr;
        const MeshFace* face = nullptr;
        int faceId = -1;
    };

    if (pieces.empty()) return false;

    std::vector<PatchFaceReference> patchFaces;
    for (const FieldPiece& piece : pieces) {
        if (!piece.topology) {
            broadcast("Error: ",
                      "partition VTK patch topology is missing.");
            return false;
        }
        const CellFaceMesh& topology = *piece.topology;
        for (std::size_t faceId = 0; faceId < topology.faces.size();
             ++faceId) {
            const MeshFace& face = topology.faces[faceId];
            if (face.patchId < 0) continue;
            patchFaces.push_back({
                &piece,
                &face,
                (int)faceId
            });
        }
    }

    std::ofstream out(fileName);
    if (!out.is_open()) {
        broadcast("Error: ", "Cannot open file " + fileName);
        return false;
    }
    out << std::setprecision(17);

    const std::size_t pointCount = patchFaces.size() * 4;
    const std::size_t cellCount = patchFaces.size();

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << pointCount
        << "\" NumberOfCells=\"" << cellCount << "\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (const PatchFaceReference& ref : patchFaces) {
        const CellFaceMesh& topology = *ref.piece->topology;
        for (int pointId : ref.face->pointIds) {
            if (pointId < 0 || pointId >= (int)topology.points.size()) {
                broadcast("Error: ",
                          "partition VTK patch face point index mismatch.");
                return false;
            }
            const MeshPoint& point = topology.points[(std::size_t)pointId];
            out << point.x << " "
                << point.y << " "
                << point.z << "\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    out << "      <CellData>\n";
    auto writeIntArray = [&](const char* type,
                             const char* name,
                             auto valueAt) {
        out << "        <DataArray type=\"" << type << "\" Name=\""
            << name << "\" format=\"ascii\">\n";
        for (const PatchFaceReference& ref : patchFaces) {
            out << valueAt(ref) << "\n";
        }
        out << "        </DataArray>\n";
    };
    writeIntArray("Int32", "PatchId",
                  [](const PatchFaceReference& ref) {
                      return ref.face->patchId;
                  });
    writeIntArray("Int32", "OwnerCell",
                  [](const PatchFaceReference& ref) {
                      return ref.face->ownerCell;
                  });
    writeIntArray("Int32", "NeighbourCell",
                  [](const PatchFaceReference& ref) {
                      return ref.face->neighbourCell;
                  });
    writeIntArray("Int32", "FaceId",
                  [](const PatchFaceReference& ref) {
                      return ref.faceId;
                  });
    writeIntArray("Int32", "BlockId",
                  [](const PatchFaceReference& ref) {
                      return ref.piece->blockId;
                  });
    writeIntArray("Int32", "Direction",
                  [](const PatchFaceReference& ref) {
                      return ref.face->direction;
                  });
    writeIntArray("Int32", "OwnerRank",
                  [](const PatchFaceReference& ref) {
                      return ref.face->ownerRank;
                  });
    writeIntArray("Int32", "NeighbourRank",
                  [](const PatchFaceReference& ref) {
                      return ref.face->neighbourRank;
                  });
    out << "      </CellData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
           "format=\"ascii\">\n";
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
        const std::size_t first = cell * 4;
        out << first << " " << first + 1 << " "
            << first + 2 << " " << first + 3 << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int64\" Name=\"offsets\" "
           "format=\"ascii\">\n";
    for (std::size_t cell = 1; cell <= cellCount; ++cell) {
        out << cell * 4 << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" "
           "format=\"ascii\">\n";
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
        out << "9\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}

bool ResultWriter::writePartitionIndex(
    const std::string& fileName,
    const std::vector<std::string>& partitionFileNames,
    const std::vector<std::string>& patchFileNames) {
    (void)patchFileNames;
    std::ofstream out(fileName);
    if (!out.is_open()) {
        broadcast("Error: ", "Cannot open file " + fileName);
        return false;
    }
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
           "byte_order=\"LittleEndian\">\n";
    out << "  <vtkMultiBlockDataSet>\n";
    for (size_t partitionId = 0;
         partitionId < partitionFileNames.size();
         ++partitionId) {
        out << "    <DataSet index=\"" << partitionId
            << "\" name=\"partition" << partitionId
            << "\" file=\""
            << fileNameOnly(partitionFileNames[partitionId])
            << "\"/>\n";
    }
    out << "  </vtkMultiBlockDataSet>\n";
    out << "</VTKFile>\n";
    broadcast("VTK multiblock file generated: ", fileName);
    return true;
}

void ResultWriter::removeLegacyPatchOutputs(const std::string& base) const {
    for (const PieceLayout& layout : pieces_) {
        std::ostringstream fileName;
        fileName << config_.caseDir << "/" << config_.outputDir << "/" << base
                 << "_block" << std::setw(3) << std::setfill('0')
                 << layout.blockId << ".vts";
        std::remove(fileName.str().c_str());
    }
    int partitionCount = 0;
    for (const PieceLayout& layout : pieces_) {
        partitionCount = std::max(partitionCount, layout.partitionId + 1);
    }
    for (int partitionId = 0;
         partitionId < partitionCount;
         ++partitionId) {
        std::ostringstream fileName;
        fileName << config_.caseDir << "/" << config_.outputDir << "/" << base
                 << "_partition" << std::setw(3) << std::setfill('0')
                 << partitionId << ".vts";
        std::remove(fileName.str().c_str());
    }
}

void ResultWriter::writeMultiBlockIndex(const std::string& fileName,
                              const std::vector<std::string>& blockFileNames) {
    std::ofstream out(fileName);
    if (!out.is_open()) {
        broadcast("Error: ", "Cannot open file " + fileName);
        return;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <vtkMultiBlockDataSet>\n";
    const bool groupedByPartition =
        pieces_.size() == blockFileNames.size() &&
        std::all_of(pieces_.begin(), pieces_.end(),
                    [](const PieceLayout& piece) {
                        return piece.partitionId >= 0;
                    });

    if (groupedByPartition) {
        std::map<int, std::vector<size_t>> partitionPieces;
        for (size_t i = 0; i < pieces_.size(); ++i) {
            partitionPieces[pieces_[i].partitionId].push_back(i);
        }

        const size_t dot = fileName.find_last_of('.');
        const std::string base =
            dot == std::string::npos ? fileName : fileName.substr(0, dot);
        for (const auto& [partitionId, pieceIds] : partitionPieces) {
            std::ostringstream legacyPartitionName;
            legacyPartitionName << base << "_partition"
                                << std::setw(3) << std::setfill('0')
                                << partitionId << ".vts";
            std::remove(legacyPartitionName.str().c_str());

            out << "    <Block index=\"" << partitionId
                << "\" name=\"partition" << partitionId << "\">\n";
            for (size_t localPiece = 0;
                 localPiece < pieceIds.size();
                 ++localPiece) {
                const size_t pieceId = pieceIds[localPiece];
                out << "      <DataSet index=\"" << localPiece
                    << "\" name=\"patch"
                    << pieces_[pieceId].blockId
                    << "\" file=\""
                    << fileNameOnly(blockFileNames[pieceId])
                    << "\"/>\n";
            }
            out << "    </Block>\n";
        }
    } else {
        for (size_t i = 0; i < blockFileNames.size(); ++i) {
            out << "    <DataSet index=\"" << i << "\" name=\"block" << i
                << "\" file=\"" << fileNameOnly(blockFileNames[i]) << "\"/>\n";
        }
    }
    out << "  </vtkMultiBlockDataSet>\n";
    out << "</VTKFile>\n";
    out.close();

    broadcast("VTK multiblock file generated: ", fileName);
}

// ============================================================
//  PVD 集合文件
// ============================================================

void ResultWriter::recordAndWritePVD(double timestep, const std::string& dataFileName) {
    // 仅在 root 进程记录 (多块时已由调用方保证)
    if (parallel_ && parallel_->active() && parallel_->rank() != 0) return;

    if (!pvdEntries_.empty() &&
        pvdEntries_.back().file == dataFileName &&
        std::abs(pvdEntries_.back().timestep - timestep) < 1.0e-14) {
        return;
    }

    pvdEntries_.push_back({timestep, dataFileName});
    flushPVD();
}

void ResultWriter::flushPVD() const {
    if (pvdEntries_.empty()) return;

    std::string prefix = outputPrefix(config_.jobName, config_.caseName);
    std::string pvdPath = config_.caseDir + "/" + config_.outputDir + "/" + prefix + ".pvd";

    std::ofstream out(pvdPath);
    if (!out.is_open()) {
        broadcast("Error: ", "Cannot open PVD file " + pvdPath);
        return;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (const auto& entry : pvdEntries_) {
        out << "    <DataSet timestep=\"" << std::scientific << std::setprecision(12)
            << entry.timestep << "\" group=\"\" part=\"0\" file=\""
            << entry.file << "\"/>\n";
    }
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
    out.close();

    broadcast("PVD file updated: ", pvdPath);
}

} // namespace SF
