/// @file SF_vtkStructuredWriter.cpp
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*----------IO private implementation---------*/

#include "SF_resultWriter.h"
#include "core/field/SF_field.h"
#include "core/interfaces/SF_log.h"

#include "private/SF_vtkUtils.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <map>

namespace SF {

using namespace IOPrivate;

void ResultWriter::writeToDisk(const SF::Field& field,
                     const std::string& fileName,
                     bool announce,
                     const std::vector<ScalarField>& suppliedFields) {
    const auto extraScalars=fieldProvider_ ? fieldProvider_(field,suppliedFields) : suppliedFields;
    std::ofstream out(fileName);
    if (!out.is_open()) {
        broadcast("Error: ", "Cannot open file " + fileName);
        return;
    }
    out << std::setprecision(17);
    int N_GHOST = field.NG();
    int NX = field.NX(), NY = field.NY(), NZ = field.NZ();

    int istart = N_GHOST, iend = NX + N_GHOST;
    int jstart = N_GHOST, jend = NY + N_GHOST;
    int kstart = N_GHOST, kend = NZ + N_GHOST;

    int nx_out = iend - istart;
    int ny_out = jend - jstart;
    int nz_out = kend - kstart;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <StructuredGrid WholeExtent=\"0 " << nx_out-1 << " 0 " << ny_out-1 << " 0 " << nz_out-1 << "\">\n";
    out << "    <Piece Extent=\"0 " << nx_out-1 << " 0 " << ny_out-1 << " 0 " << nz_out-1 << "\">\n";

    // VTK/ParaView 的 ghost 标志在 PointData 和 CellData 中使用不同位值:
    // - PointData: 2 表示 HIDDENPOINT, 可直接隐藏 IBM ghost/solid 点；
    // - CellData : 32 表示 HIDDENCELL, 不能把 2 当作隐藏单元使用。
    // 这里保持旧版输出语义: 对 IBM 非流体点写 PointData vtkGhostType=2,
    // 让 ParaView 自动 blank 点；IBMCellType 仅作为诊断字段保留。
    constexpr int VTK_VISIBLE_POINT = 0;
    constexpr int VTK_HIDDEN_POINT = 2;

    auto isHiddenIBMPoint = [&](int i, int j, int k) {
        return config_.ibmOutputEnabled && !field.IBMFluidMask(i, j, k);
    };

    // Points
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int k = kstart; k < kend; ++k) {
        for (int j = jstart; j < jend; ++j) {
            for (int i = istart; i < iend; ++i) {
                out << field.X(i,j,k) << " " << field.Y(i,j,k) << " " << field.Z(i,j,k) << " ";
            }
            out << "\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // PointData
    out << "      <PointData>\n";

    // ── BC_ID: 输出参与计算的BC配置集合ID ──
    // 0=内点/内部共享点/未被 OpenFOAM boundaryField 引用的集合, 正整数=引用的 set。
    // 单个point可同时属于多个物理set，因此额外输出变量级ID数组；
    // 变量级数组按该变量在 boundaryField 中的配置顺序覆盖，贴近实际 BC 应用过程。
    {
        const bool writeTemperatureBoundaryIds =
            hasTemperatureBoundaryScalar(extraScalars);
        const auto& allSets = field.getAllSets();

        // 只收集 boundaryField 实际引用的集合名。set 名称只是标签，不携带物理语义。
        std::vector<std::string> bcNames;
        appendBoundaryZoneNames(config_.densityBoundary, bcNames);
        appendBoundaryZoneNames(config_.velocityBoundary, bcNames);
        appendBoundaryZoneNames(config_.pressureBoundary, bcNames);
        if (writeTemperatureBoundaryIds) {
            appendBoundaryZoneNames(config_.temperatureBoundary, bcNames);
        }

        // 建立name→ID映射(1-based)
        std::map<std::string, int> setToId;
        for (size_t sid = 0; sid < bcNames.size(); ++sid)
            setToId[bcNames[sid]] = (int)(sid + 1);

        auto markSetId = [&](std::vector<int>& mask, const std::string& name) {
            auto sidIt = setToId.find(name);
            if (sidIt == setToId.end()) return;
            auto setIt = allSets.find(name);
            if (setIt == allSets.end()) return;
            for (int idx : setIt->second) {
                if (idx < 0 || idx >= field.TotalSize()) continue;
                mask[(size_t)idx] = sidIt->second;
            }
        };

        auto writeIdArray = [&](const char* arrayName,
                                const std::vector<int>& mask) {
            out << "        <DataArray type=\"Int32\" Name=\"" << arrayName
                << "\" format=\"ascii\">\n";
            for (int k = kstart; k < kend; ++k) {
                for (int j = jstart; j < jend; ++j) {
                    for (int i = istart; i < iend; ++i) {
                        int fidx = field.getIdx(i, j, k);
                        out << mask[(size_t)fidx] << " ";
                    }
                    out << "\n";
                }
            }
            out << "        </DataArray>\n";
        };

        std::vector<int> rhoMask((size_t)field.TotalSize(), 0);
        std::vector<int> uMask((size_t)field.TotalSize(), 0);
        std::vector<int> pMask((size_t)field.TotalSize(), 0);
        std::vector<int> tMask((size_t)field.TotalSize(), 0);
        for (const auto& bc : config_.densityBoundary) markSetId(rhoMask, bc.name);
        for (const auto& bc : config_.velocityBoundary) markSetId(uMask, bc.name);
        for (const auto& bc : config_.pressureBoundary) markSetId(pMask, bc.name);
        if (writeTemperatureBoundaryIds) {
            for (const auto& bc : config_.temperatureBoundary) {
                markSetId(tMask, bc.name);
            }
        }

        std::vector<int> bcIdMask = rhoMask;
        for (size_t idx = 0; idx < bcIdMask.size(); ++idx) {
            if (uMask[idx] != 0) bcIdMask[idx] = uMask[idx];
            if (pMask[idx] != 0) bcIdMask[idx] = pMask[idx];
            if (writeTemperatureBoundaryIds && tMask[idx] != 0) {
                bcIdMask[idx] = tMask[idx];
            }
        }

        writeIdArray("BC_ID", bcIdMask);
        writeIdArray("BC_ID_rho", rhoMask);
        writeIdArray("BC_ID_U", uMask);
        writeIdArray("BC_ID_p", pMask);
        if (writeTemperatureBoundaryIds) writeIdArray("BC_ID_T", tMask);

        // 输出ID→名称图例
        out << "        <!-- BC_ID legend: 0=interior/internal-shared/unconfigured";
        for (size_t sid = 0; sid < bcNames.size(); ++sid)
            out << ", " << (sid + 1) << "=" << bcNames[sid];
        out << "; variable BC_ID arrays use the same legend and follow per-variable boundaryField order on overlaps -->\n";
    }

    for (const auto& array : extraScalars) {
        if ((!array.valueAt&&!array.componentAt) || array.name.empty())
            throw std::runtime_error("Invalid registered output array");
        out << "        <DataArray type=\"Float64\" Name=\"" << xmlEscape(array.name)
            << "\" NumberOfComponents=\"" << array.components << "\" format=\"ascii\">\n";
        for(int k=kstart;k<kend;++k)for(int j=jstart;j<jend;++j)for(int i=istart;i<iend;++i) {
            for(int c=0;c<array.components;++c)
                out << (isHiddenIBMPoint(i,j,k)?0.0:array.value(i,j,k,c)) << ' ';
            out << '\n';
        }
        out << "        </DataArray>\n";
    }

    if (config_.ibmOutputEnabled) {
        out << "        <DataArray type=\"UInt8\" Name=\"vtkGhostType\" format=\"ascii\">\n";
        for (int k = kstart; k < kend; ++k) {
            for (int j = jstart; j < jend; ++j) {
                for (int i = istart; i < iend; ++i) {
                    out << (isHiddenIBMPoint(i, j, k) ? VTK_HIDDEN_POINT : VTK_VISIBLE_POINT) << " ";
                }
                out << "\n";
            }
        }
        out << "        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"IBMCellType\" format=\"ascii\">\n";
        for (int k = kstart; k < kend; ++k) {
            for (int j = jstart; j < jend; ++j) {
                for (int i = istart; i < iend; ++i) {
                    out << field.CellFlag(i, j, k) << " ";
                }
                out << "\n";
            }
        }
        out << "        </DataArray>\n";
    }

    out << "      </PointData>\n";

    out << "    </Piece>\n";
    out << "  </StructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close();
    if (announce) broadcast("VTK file generated: ", fileName);
}

} // namespace SF
