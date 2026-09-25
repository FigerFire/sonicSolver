/// @file SF_generatedMeshWriter.cpp
/// @brief SFM/生成网格文件的 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_generatedMeshWriter.h"

#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace SF {
namespace GeneratedMeshWriter {
namespace {

std::string fileNameOnly(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

struct WriterFace {
    std::array<int, 4> points{};
    int owner = -1;
    int neighbour = -1;
    int patch = -1;
    int direction = -1;
};

struct WriterPatch {
    std::string name;
    std::vector<int> faces;
};

int pointIndex(int i, int j, int k, int nx, int ny) {
    return (k * ny + j) * nx + i;
}

std::array<int, 4> sortedFaceKey(std::array<int, 4> ids) {
    std::sort(ids.begin(), ids.end());
    return ids;
}

bool faceInSet(const WriterFace& face,
               const std::vector<unsigned char>& mask) {
    for (int pointId : face.points) {
        if (pointId < 0 || pointId >= (int)mask.size()) return false;
        if (!mask[(std::size_t)pointId]) return false;
    }
    return true;
}

void writePointSetSections(
        std::ofstream& sfm,
        const std::vector<std::map<std::string, std::vector<int>>>& pointSets,
        int blockId) {
    if (blockId < 0 || blockId >= (int)pointSets.size()) return;
    const auto& sets = pointSets[(std::size_t)blockId];
    if (sets.empty()) return;

    sfm << "#PointSet\n";
    for (const auto& [name, ids] : sets) {
        sfm << name << " " << ids.size() << "\n";
        for (std::size_t n = 0; n < ids.size(); ++n) {
            sfm << ids[n]
                << ((n + 1) % 10 == 0 || n + 1 == ids.size()
                    ? "\n" : " ");
        }
    }
    sfm << "\n";
}

void writeCellFaceSections(
        std::ofstream& sfm,
        int nx, int ny, int nz,
        const std::vector<std::map<std::string, std::vector<int>>>& facePatchSets,
        int blockId) {
    const int nxCells = nx - 1;
    const int nyCells = ny - 1;
    const int nzCells = nz - 1;
    if (nxCells <= 0 || nyCells <= 0 || nzCells <= 0) return;

    std::vector<std::array<int, 8>> cells;
    std::vector<WriterFace> faces;
    std::map<std::array<int, 4>, int> faceMap;

    auto addFace = [&](int owner,
                       std::array<int, 4> points,
                       int direction) {
        const std::array<int, 4> key = sortedFaceKey(points);
        auto found = faceMap.find(key);
        if (found != faceMap.end()) {
            faces[(std::size_t)found->second].neighbour = owner;
            return found->second;
        }

        WriterFace face;
        face.points = points;
        face.owner = owner;
        face.direction = direction;
        const int faceId = (int)faces.size();
        faces.push_back(face);
        faceMap.emplace(key, faceId);
        return faceId;
    };

    for (int k = 0; k < nzCells; ++k) {
        for (int j = 0; j < nyCells; ++j) {
            for (int i = 0; i < nxCells; ++i) {
                const int p000 = pointIndex(i,     j,     k,     nx, ny);
                const int p100 = pointIndex(i + 1, j,     k,     nx, ny);
                const int p110 = pointIndex(i + 1, j + 1, k,     nx, ny);
                const int p010 = pointIndex(i,     j + 1, k,     nx, ny);
                const int p001 = pointIndex(i,     j,     k + 1, nx, ny);
                const int p101 = pointIndex(i + 1, j,     k + 1, nx, ny);
                const int p111 = pointIndex(i + 1, j + 1, k + 1, nx, ny);
                const int p011 = pointIndex(i,     j + 1, k + 1, nx, ny);

                const int cellId = (int)cells.size();
                cells.push_back({p000, p100, p110, p010,
                                 p001, p101, p111, p011});
                addFace(cellId, {p000, p010, p011, p001}, 0);
                addFace(cellId, {p100, p101, p111, p110}, 0);
                addFace(cellId, {p000, p001, p101, p100}, 1);
                addFace(cellId, {p010, p110, p111, p011}, 1);
                addFace(cellId, {p000, p100, p110, p010}, 2);
                addFace(cellId, {p001, p011, p111, p101}, 2);
            }
        }
    }

    std::vector<WriterPatch> patches;
    if (blockId >= 0 && blockId < (int)facePatchSets.size()) {
        const auto& sets = facePatchSets[(std::size_t)blockId];
        for (const auto& [name, ids] : sets) {
            if (name == "all") continue;
            std::vector<unsigned char> mask((std::size_t)(nx * ny * nz), 0);
            for (int id : ids) {
                if (id >= 0 && id < (int)mask.size()) {
                    mask[(std::size_t)id] = 1;
                }
            }

            WriterPatch patch;
            patch.name = name;
            for (int faceId = 0; faceId < (int)faces.size(); ++faceId) {
                WriterFace& face = faces[(std::size_t)faceId];
                if (face.neighbour >= 0) continue;
                if (!faceInSet(face, mask)) continue;
                face.patch = (int)patches.size();
                patch.faces.push_back(faceId);
            }
            if (!patch.faces.empty()) patches.push_back(std::move(patch));
        }
    }

    sfm << "#Cell\n";
    for (const auto& cell : cells) {
        for (int n = 0; n < 8; ++n) {
            sfm << cell[(std::size_t)n] << (n == 7 ? "\n" : " ");
        }
    }
    sfm << "\n";

    sfm << "#Face\n";
    for (const WriterFace& face : faces) {
        for (int n = 0; n < 4; ++n) {
            sfm << face.points[(std::size_t)n] << " ";
        }
        sfm << face.owner << " " << face.neighbour << " "
            << face.patch << " " << face.direction << "\n";
    }
    sfm << "\n";

    sfm << "#Patch\n";
    for (const WriterPatch& patch : patches) {
        sfm << patch.name << " patch " << patch.faces.size() << "\n";
        for (std::size_t n = 0; n < patch.faces.size(); ++n) {
            sfm << patch.faces[n]
                << ((n + 1) % 10 == 0 || n + 1 == patch.faces.size()
                    ? "\n" : " ");
        }
    }
    sfm << "\n";
}

} // namespace

void writeBlockVTK(const std::string& filename,
                   const std::vector<double>& allX,
                   const std::vector<double>& allY,
                   const std::vector<double>& allZ,
                   size_t offset,
                   int nx, int ny, int nz) {
    std::ofstream out(filename);
    if (!out) {
        broadcast("Error: ", "Cannot open " + filename);
        return;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <StructuredGrid WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    out << "    <Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                size_t idx = offset + ((size_t)k * ny + j) * nx + i;
                out << (float)allX[idx] << " "
                    << (float)allY[idx] << " "
                    << (float)allZ[idx] << " ";
            }
            out << "\n";
        }
    }

    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "    </Piece>\n";
    out << "  </StructuredGrid>\n";
    out << "</VTKFile>\n";
}

void writeMultiZoneSFM(
        const std::string& filename,
        const std::vector<double>& allX,
        const std::vector<double>& allY,
        const std::vector<double>& allZ,
        const std::vector<int>& blockSizes,
        int nGhost,
        const std::vector<std::map<std::string, std::vector<int>>>& facePatchSets,
        const std::vector<std::map<std::string, std::vector<int>>>& pointSets) {
    std::ofstream sfm(filename);
    if (!sfm) {
        broadcast("Error: ", "Cannot open " + filename);
        return;
    }

    const int nBlocks = (int)blockSizes.size() / 3;
    size_t offset = 0;
    sfm << std::setprecision(15);
    for (int b = 0; b < nBlocks; ++b) {
        const int nx = blockSizes[b * 3];
        const int ny = blockSizes[b * 3 + 1];
        const int nz = blockSizes[b * 3 + 2];

        sfm << "#Information\n" << nx << " " << ny << " " << nz << " " << nGhost << "\n\n";
        sfm << "#Point\n";
        const size_t nPts = (size_t)nx * ny * nz;
        for (size_t n = 0; n < nPts; ++n) {
            const size_t idx = offset + n;
            sfm << allX[idx] << " " << allY[idx] << " " << allZ[idx] << "\n";
        }
        sfm << "\n";

        writeCellFaceSections(sfm, nx, ny, nz, facePatchSets, b);
        writePointSetSections(sfm, pointSets, b);
        offset += nPts;
    }
}

void writeVTM(const std::string& filename,
              const std::vector<std::string>& blockFiles) {
    std::ofstream vtm(filename);
    if (!vtm) {
        broadcast("Error: ", "Cannot open " + filename);
        return;
    }

    vtm << "<?xml version=\"1.0\"?>\n";
    vtm << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    vtm << "  <vtkMultiBlockDataSet>\n";
    for (size_t i = 0; i < blockFiles.size(); ++i) {
        vtm << "    <DataSet index=\"" << i << "\" name=\"block" << i
            << "\" file=\"" << fileNameOnly(blockFiles[i]) << "\"/>\n";
    }
    vtm << "  </vtkMultiBlockDataSet>\n";
    vtm << "</VTKFile>\n";
}

} // namespace GeneratedMeshWriter
} // namespace SF
