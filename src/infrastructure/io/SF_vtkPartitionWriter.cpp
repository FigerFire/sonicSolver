/// @file SF_vtkPartitionWriter.cpp
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*----------IO private implementation---------*/

#include "SF_resultWriter.h"
#include "core/field/SF_field.h"
#include "core/interfaces/SF_log.h"

#include "private/SF_vtkUtils.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <unordered_map>

namespace SF {

using namespace IOPrivate;

bool ResultWriter::writePartitionUnstructuredVTK(
    const std::string& fileName,
    const std::vector<FieldPiece>& suppliedPieces) {
    auto pieces=suppliedPieces;
    if(fieldProvider_)for(auto& piece:pieces) {
        if(piece.field)piece.extraScalars=fieldProvider_(*piece.field,piece.extraScalars);
    }
    struct PointReference {
        const FieldPiece* piece = nullptr;
        const Field* field = nullptr;
        int i = 0;
        int j = 0;
        int k = 0;
        int globalPointId = -1;
    };
    struct BoundaryMasks {
        std::vector<int> rho;
        std::vector<int> velocity;
        std::vector<int> pressure;
        std::vector<int> temperature;
    };

    if (pieces.empty()) return false;

    const bool writeTemperatureBoundaryIds =
        !pieces.empty()
        && hasTemperatureBoundaryScalar(pieces.front().extraScalars);
    std::vector<std::string> bcNames;
    appendBoundaryZoneNames(config_.densityBoundary, bcNames);
    appendBoundaryZoneNames(config_.velocityBoundary, bcNames);
    appendBoundaryZoneNames(config_.pressureBoundary, bcNames);
    if (writeTemperatureBoundaryIds) {
        appendBoundaryZoneNames(config_.temperatureBoundary, bcNames);
    }
    std::map<std::string, int> setToId;
    for (size_t sid = 0; sid < bcNames.size(); ++sid) {
        setToId[bcNames[sid]] = (int)sid + 1;
    }

    std::unordered_map<const Field*, BoundaryMasks> masksByField;
    auto buildMasks = [&](const Field& field) -> const BoundaryMasks& {
        auto found = masksByField.find(&field);
        if (found != masksByField.end()) return found->second;

        BoundaryMasks masks;
        masks.rho.assign((size_t)field.TotalSize(), 0);
        masks.velocity.assign((size_t)field.TotalSize(), 0);
        masks.pressure.assign((size_t)field.TotalSize(), 0);
        masks.temperature.assign((size_t)field.TotalSize(), 0);
        const auto& allSets = field.getAllSets();
        auto mark = [&](std::vector<int>& values,
                        const std::string& name) {
            auto id = setToId.find(name);
            auto set = allSets.find(name);
            if (id == setToId.end() || set == allSets.end()) return;
            for (int index : set->second) {
                if (index >= 0 && index < field.TotalSize()) {
                    values[(size_t)index] = id->second;
                }
            }
        };
        for (const auto& bc : config_.densityBoundary) mark(masks.rho, bc.name);
        for (const auto& bc : config_.velocityBoundary) mark(masks.velocity, bc.name);
        for (const auto& bc : config_.pressureBoundary) mark(masks.pressure, bc.name);
        if (writeTemperatureBoundaryIds) {
            for (const auto& bc : config_.temperatureBoundary) {
                mark(masks.temperature, bc.name);
            }
        }
        return masksByField.emplace(&field, std::move(masks)).first->second;
    };

    std::unordered_map<int, int> globalToOutput;
    std::vector<PointReference> points;
    std::size_t cellCount = 0;
    for (const FieldPiece& piece : pieces) {
        const Field& field = *piece.field;
        const std::vector<int>& globalPointIds = *piece.globalPointIds;
        const std::size_t expected =
            (std::size_t)field.NX() * field.NY() * field.NZ();
        if (globalPointIds.size() != expected) {
            broadcast("Error: ",
                      "partition VTK global point topology size mismatch.");
            return false;
        }

        const int ng = field.NG();
        for (int k = 0; k < field.NZ(); ++k) {
            for (int j = 0; j < field.NY(); ++j) {
                for (int i = 0; i < field.NX(); ++i) {
                    const std::size_t local =
                        ((std::size_t)k * field.NY() + j) * field.NX() + i;
                    const int globalId = globalPointIds[local];
                    if (globalId < 0) return false;
                    if (globalToOutput.find(globalId) != globalToOutput.end()) {
                        continue;
                    }
                    const int outputId = (int)points.size();
                    globalToOutput.emplace(globalId, outputId);
                    points.push_back({
                        &piece,
                        &field, i + ng, j + ng, k + ng, globalId
                    });
                }
            }
        }
        if (field.NX() > 1 && field.NY() > 1 && field.NZ() > 1) {
            cellCount += (std::size_t)(field.NX() - 1)
                       * (field.NY() - 1)
                       * (field.NZ() - 1);
        }
    }

    std::ofstream out(fileName);
    if (!out.is_open()) {
        broadcast("Error: ", "Cannot open file " + fileName);
        return false;
    }
    out << std::setprecision(17);

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << points.size()
        << "\" NumberOfCells=\"" << cellCount << "\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (const PointReference& point : points) {
        out << point.field->X(point.i, point.j, point.k) << " "
            << point.field->Y(point.i, point.j, point.k) << " "
            << point.field->Z(point.i, point.j, point.k) << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <PointData>\n";

    auto hidden = [&](const PointReference& point) {
        return config_.ibmOutputEnabled &&
               !point.field->IBMFluidMask(point.i, point.j, point.k);
    };
    auto writeScalar = [&](const char* name, auto valueAt) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name
            << "\" format=\"ascii\">\n";
        for (const PointReference& point : points) {
            const double value = hidden(point) ? 0.0 : valueAt(point);
            out << value << "\n";
        }
        out << "        </DataArray>\n";
    };
    auto writeInteger = [&](const char* type,
                            const char* name,
                            auto valueAt) {
        out << "        <DataArray type=\"" << type << "\" Name=\""
            << name << "\" format=\"ascii\">\n";
        for (const PointReference& point : points) {
            out << valueAt(point) << "\n";
        }
        out << "        </DataArray>\n";
    };

    auto boundaryId = [&](const PointReference& point, int variable) {
        const BoundaryMasks& masks = buildMasks(*point.field);
        const int index = point.field->getIdx(point.i, point.j, point.k);
        if (variable == 0) return masks.rho[(size_t)index];
        if (variable == 1) return masks.velocity[(size_t)index];
        if (variable == 2) return masks.pressure[(size_t)index];
        return masks.temperature[(size_t)index];
    };
    writeInteger("Int32", "BC_ID_rho",
                 [&](const PointReference& point) {
                     return boundaryId(point, 0);
                 });
    writeInteger("Int32", "BC_ID_U",
                 [&](const PointReference& point) {
                     return boundaryId(point, 1);
                 });
    writeInteger("Int32", "BC_ID_p",
                 [&](const PointReference& point) {
                     return boundaryId(point, 2);
                 });
    if (writeTemperatureBoundaryIds) {
        writeInteger("Int32", "BC_ID_T",
                     [&](const PointReference& point) {
                         return boundaryId(point, 3);
                     });
    }
    writeInteger("Int32", "BC_ID",
                 [&](const PointReference& point) {
                     int id = boundaryId(point, 0);
                     const int velocityId = boundaryId(point, 1);
                     const int pressureId = boundaryId(point, 2);
                     if (velocityId != 0) id = velocityId;
                     if (pressureId != 0) id = pressureId;
                     if (writeTemperatureBoundaryIds) {
                         const int temperatureId = boundaryId(point, 3);
                         if (temperatureId != 0) id = temperatureId;
                     }
                     return id;
                 });
    writeInteger("Int64", "GlobalPointId",
                 [](const PointReference& point) {
                     return point.globalPointId;
                 });

    if (!pieces.empty()) {
        for(const auto& array:pieces.front().extraScalars) {
            out << "        <DataArray type=\"Float64\" Name=\"" << xmlEscape(array.name)
                << "\" NumberOfComponents=\"" << array.components << "\" format=\"ascii\">\n";
            for(const auto& point:points) {
                const ScalarField* selected=nullptr;
                for(const auto& candidate:point.piece->extraScalars)
                    if(candidate.name==array.name) { selected=&candidate;break; }
                if(!selected||selected->components!=array.components)
                    throw std::runtime_error("Inconsistent partition output field: "+array.name);
                for(int c=0;c<array.components;++c)
                    out << (hidden(point)?0.0:selected->value(point.i,point.j,point.k,c)) << ' ';
                out << '\n';
            }
            out << "        </DataArray>\n";
        }
    }

    if (config_.ibmOutputEnabled) {
        writeInteger("UInt8", "vtkGhostType",
                     [&](const PointReference& point) {
                         return hidden(point) ? 2 : 0;
                     });
        writeInteger("Int32", "IBMCellType",
                     [](const PointReference& point) {
                         return point.field->CellFlag(
                             point.i, point.j, point.k);
                     });
    }
    out << "      </PointData>\n";

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
           "format=\"ascii\">\n";
    auto outputPointId = [&](const FieldPiece& piece,
                             int i, int j, int k) {
        const Field& field = *piece.field;
        const std::size_t local =
            ((std::size_t)k * field.NY() + j) * field.NX() + i;
        return globalToOutput.at((*piece.globalPointIds)[local]);
    };
    for (const FieldPiece& piece : pieces) {
        const Field& field = *piece.field;
        for (int k = 0; k + 1 < field.NZ(); ++k) {
            for (int j = 0; j + 1 < field.NY(); ++j) {
                for (int i = 0; i + 1 < field.NX(); ++i) {
                    out << outputPointId(piece, i, j, k) << " "
                        << outputPointId(piece, i + 1, j, k) << " "
                        << outputPointId(piece, i + 1, j + 1, k) << " "
                        << outputPointId(piece, i, j + 1, k) << " "
                        << outputPointId(piece, i, j, k + 1) << " "
                        << outputPointId(piece, i + 1, j, k + 1) << " "
                        << outputPointId(piece, i + 1, j + 1, k + 1) << " "
                        << outputPointId(piece, i, j + 1, k + 1) << "\n";
                }
            }
        }
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int64\" Name=\"offsets\" "
           "format=\"ascii\">\n";
    for (std::size_t cell = 1; cell <= cellCount; ++cell) {
        out << cell * 8 << "\n";
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" "
           "format=\"ascii\">\n";
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
        out << "12\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    return true;
}

} // namespace SF
