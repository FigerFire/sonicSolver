/// @file SF_mesh.cpp
/// @brief 结构/多块网格拓扑、度量与加载实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/



#include "SF_mesh.h"
#include "infrastructure/io/mesh/SF_meshParameters.h"
#include "SF_MultiBlockMesh.h"
#include "SF_meshGen.h"
#include "SF_init.h"
#include "core/interfaces/SF_log.h"
#include "methods/numerics/structured/SF_structured.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <utility>

namespace   SF{
    namespace {
        bool isDirectPath(const std::string& path)
        {
            return !path.empty() &&
                   (path[0] == '/' || path[0] == '\\' || path.rfind("..", 0) == 0);
        }

        std::string meshOutputBaseFromTarget(std::string target)
        {
            size_t dot = target.find_last_of('.');
            if (dot != std::string::npos) target = target.substr(0, dot);

            const std::string setSuffix = "_sets";
            if (target.size() > setSuffix.size() &&
                target.compare(target.size() - setSuffix.size(), setSuffix.size(), setSuffix) == 0) {
                target.erase(target.size() - setSuffix.size());
            }

            size_t blockPos = target.rfind("_block");
            if (blockPos != std::string::npos && blockPos + 6 < target.size()) {
                bool numericSuffix = true;
                for (size_t i = blockPos + 6; i < target.size(); ++i) {
                    if (!std::isdigit(static_cast<unsigned char>(target[i]))) {
                        numericSuffix = false;
                        break;
                    }
                }
                if (numericSuffix) target.erase(blockPos);
            }

            return target;
        }

        bool validateGhostLayerCount(
            int nGhost,
            const std::string& label,
            const MeshRuntimeConfig& config)
        {
            const int required = config.requiredGhostLayers;
            if (nGhost >= required) return true;

            std::ostringstream oss;
            oss << label << " declares nGhost=" << nGhost
                << ", but convection=" << config.convectionScheme
                << " and ILW=" << config.ilwOrder
                << " require at least " << required
                << " ghost layers. Please regenerate or edit the SFM "
                << "#Information ghost-count field explicitly.";
            SF::broadcast("Fatal: ", oss.str());
            return false;
        }

        std::map<std::string, std::vector<int>> toFieldIndexSets(
            const SF::Field& field,
            const std::map<std::string, std::vector<int>>& rawSets)
        {
            std::map<std::string, std::vector<int>> converted;
            int nx = field.NX();
            int ny = field.NY();
            int nz = field.NZ();
            int ng = field.NG();

            for (const auto& [name, indices] : rawSets) {
                auto& out = converted[name];
                out.reserve(indices.size());
                for (int rawIdx : indices) {
                    if (rawIdx < 0 || rawIdx >= nx * ny * nz) continue;
                    int i = rawIdx % nx;
                    int j = (rawIdx / nx) % ny;
                    int k = rawIdx / (nx * ny);
                    out.push_back(field.getIdx(i + ng, j + ng, k + ng));
                }
            }
            return converted;
        }

    }

    bool Mesh::createMesh(
        const CaseConfig& caseConfig,
        const MeshRuntimeConfig& config) {
        const std::string paramPath = caseConfig.caseDir
            + "/" + caseConfig.meshParameterFile;
        SF::broadcast("createMesh mode: ", paramPath);

        SF::MeshParameters params;
        if (!SF::MeshIO::readParameters(paramPath, params)) {
            SF::broadcast("Fatal: ", "Failed to parse mesh parameters.");
            return false;
        }
        if (!validateGhostLayerCount(
                params.nGhost, "Mesh parameters", config)) {
            return false;
        }

        std::vector<double> allX, allY, allZ;
        std::vector<int> blockSizes;
        SF::MeshGen::generateStructuredMesh(params, allX, allY, allZ, blockSizes);

        int nBlocks = (int)blockSizes.size() / 3;
        std::string targetMesh = caseConfig.meshFiles.empty()
                               ? "constant/mesh.sfm"
                               : caseConfig.meshFiles.front();
        targetMesh = meshOutputBaseFromTarget(targetMesh);
        std::string base = isDirectPath(targetMesh) ? targetMesh
            : caseConfig.caseDir + "/" + targetMesh;

        if (!SF::MeshGen::writeCombinedMesh(base, params, allX, allY, allZ,
                                            blockSizes, params.nGhost)) {
            return false;
        }
        SF::broadcast("Done. ", std::to_string(nBlocks) + " blocks -> " + base + ".*");
        return true;
    }

    bool Mesh::setupComplexMesh(
        SF::Field& field,
        const CaseConfig& caseConfig,
        const MeshRuntimeConfig& config) {
        MultiBlockMesh mesh;
        if (!setupMultiBlockMesh(mesh, caseConfig, config)) {
            broadcast("Fatal! " , " Mesh loading failed. Check your .SFM file." );
            return false;
        }

        if (!mesh.isSingleBlock()) {
            broadcast("Fatal: ",
                      "main solver still advances one Field; use MultiBlockMesh API for "
                      + std::to_string(mesh.size()) + " imported blocks.");
            return false;
        }

        field = std::move(mesh.block(0).field);
        broadcast("Mesh ready: ", "single block");
        return true;
    }

    bool Mesh::setupMultiBlockMesh(
        SF::MultiBlockMesh& mesh,
        const CaseConfig& caseConfig,
        const MeshRuntimeConfig& config,
        const std::function<bool(std::vector<RawMeshBlock>&)>& preprocessor) {
        return mesh.loadFiles(caseConfig.caseDir, caseConfig.meshFiles,
                              config, preprocessor);
    }

    bool Mesh::setupFieldFromRaw(SF::Field& field,
                                 int nx, int ny, int nz, int ng,
                                 const std::vector<double>& px,
                                 const std::vector<double>& py,
                                 const std::vector<double>& pz,
                                 const std::map<std::string, std::vector<int>>& sets,
                                 const MeshRuntimeConfig& config,
                                 const std::string& label) {
        if ((int)px.size() != nx * ny * nz ||
            px.size() != py.size() ||
            px.size() != pz.size()) {
            broadcast("Fatal: ", "Invalid point count while setting " + label);
            return false;
        }
        if (!validateGhostLayerCount(
                ng, label.empty() ? "SFM mesh" : label, config)) {
            return false;
        }

        field.setup(nx, ny, nz, ng);
        field.setBoundarySets(toFieldIndexSets(field, sets));

        int count = 0;
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int mi = i + ng;
                    int mj = j + ng;
                    int mk = k + ng;
                    field.X(mi, mj, mk) = px[count];
                    field.Y(mi, mj, mk) = py[count];
                    field.Z(mi, mj, mk) = pz[count];
                    ++count;
                }
            }
        }

        (void)label;
        Init::setupAll(field, config.initialConditions,
                       config.idealGasGamma,
                       config.idealGasConstant);
        generateGhostCells(field);
        initializeGhostStateFromInterior(field);
        computeMetrics(field);
        return true;
    }

    void Mesh::refreshMetrics(SF::Field& field) {
        computeMetrics(field);
    }

    //  计算度规
    void Mesh::computeMetrics(SF::Field& field) {
        int endX = field.MX() - 1;
        int endY = field.MY() - 1;
        int endZ = field.MZ() - 1;

        auto coordinate = [&](int component, int i, int j, int k) {
            if (component == 0) return field.X(i, j, k);
            if (component == 1) return field.Y(i, j, k);
            return field.Z(i, j, k);
        };
        auto shifted = [](int i, int j, int k, int axis, int offset) {
            std::array<int, 3> index{i, j, k};
            index[(size_t)axis] += offset;
            return index;
        };
        auto coordinateDerivative = [&](int component,
                                        int axis,
                                        int i, int j, int k) {
            const auto plus = shifted(i, j, k, axis, 1);
            const auto minus = shifted(i, j, k, axis, -1);
            return 0.5 * (
                coordinate(component, plus[0], plus[1], plus[2]) -
                coordinate(component, minus[0], minus[1], minus[2]));
        };
        auto productDerivative = [&](int firstComponent,
                                     int secondComponent,
                                     int innerAxis,
                                     int outerAxis,
                                     int i, int j, int k) {
            const auto plus = shifted(i, j, k, outerAxis, 1);
            const auto minus = shifted(i, j, k, outerAxis, -1);
            const double plusProduct =
                coordinate(firstComponent, plus[0], plus[1], plus[2])
                * coordinateDerivative(secondComponent, innerAxis,
                                       plus[0], plus[1], plus[2]);
            const double minusProduct =
                coordinate(firstComponent, minus[0], minus[1], minus[2])
                * coordinateDerivative(secondComponent, innerAxis,
                                       minus[0], minus[1], minus[2]);
            return 0.5 * (plusProduct - minusProduct);
        };
        auto conservativeCofactor = [&](int coordinateAxis,
                                        int component,
                                        int i, int j, int k) {
            // Thomas-Lombard 守恒 metric 的 symmetric conservative form
            // (Vinokur-Yee/SCMM)。两种等价 curl 写法对称平均，并始终使用
            // 与当前中心通量一致的同一离散导数算子。
            static constexpr int firstAxis[3] = {1, 2, 0};
            static constexpr int secondAxis[3] = {2, 0, 1};
            static constexpr int firstComponent[3] = {1, 2, 0};
            static constexpr int secondComponent[3] = {2, 0, 1};
            const int a = firstAxis[coordinateAxis];
            const int b = secondAxis[coordinateAxis];
            const int p = firstComponent[component];
            const int q = secondComponent[component];
            const double firstCurl =
                productDerivative(p, q, b, a, i, j, k)
              - productDerivative(p, q, a, b, i, j, k);
            const double symmetricCurl =
                productDerivative(q, p, a, b, i, j, k)
              - productDerivative(q, p, b, a, i, j, k);
            return 0.5 * (firstCurl + symmetricCurl);
        };

        for (int k = 1; k < endZ; ++k) {
            for (int j = 1; j < endY; ++j) {
                for (int i = 1; i < endX; ++i) {
                    const std::array<double, 3> rXi{
                        coordinateDerivative(0, 0, i, j, k),
                        coordinateDerivative(1, 0, i, j, k),
                        coordinateDerivative(2, 0, i, j, k)
                    };
                    const std::array<double, 3> rEta{
                        coordinateDerivative(0, 1, i, j, k),
                        coordinateDerivative(1, 1, i, j, k),
                        coordinateDerivative(2, 1, i, j, k)
                    };
                    const std::array<double, 3> rZeta{
                        coordinateDerivative(0, 2, i, j, k),
                        coordinateDerivative(1, 2, i, j, k),
                        coordinateDerivative(2, 2, i, j, k)
                    };
                    const std::array<double, 3> cXi{
                        conservativeCofactor(0, 0, i, j, k),
                        conservativeCofactor(0, 1, i, j, k),
                        conservativeCofactor(0, 2, i, j, k)
                    };
                    const std::array<double, 3> cEta{
                        conservativeCofactor(1, 0, i, j, k),
                        conservativeCofactor(1, 1, i, j, k),
                        conservativeCofactor(1, 2, i, j, k)
                    };
                    const std::array<double, 3> cZeta{
                        conservativeCofactor(2, 0, i, j, k),
                        conservativeCofactor(2, 1, i, j, k),
                        conservativeCofactor(2, 2, i, j, k)
                    };
                    auto dot = [](const std::array<double, 3>& a,
                                  const std::array<double, 3>& b) {
                        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                    };
                    const double J =
                        (dot(rXi, cXi) + dot(rEta, cEta)
                         + dot(rZeta, cZeta)) / 3.0;
                    if (!std::isfinite(J) || std::abs(J) < 1e-18) {
                        std::cerr << "[SF FATAL] Degenerate mesh Jacobian at ("
                                  << i << "," << j << "," << k
                                  << "), J=" << J << std::endl;
                        std::exit(1);
                    }
                    field.Jac(i,j,k) = 1.0 / J;
                    field.XiX(i,j,k) = cXi[0] * field.Jac(i,j,k);
                    field.XiY(i,j,k) = cXi[1] * field.Jac(i,j,k);
                    field.XiZ(i,j,k) = cXi[2] * field.Jac(i,j,k);
                    field.EtX(i,j,k) = cEta[0] * field.Jac(i,j,k);
                    field.EtY(i,j,k) = cEta[1] * field.Jac(i,j,k);
                    field.EtZ(i,j,k) = cEta[2] * field.Jac(i,j,k);
                    field.ZeX(i,j,k) = cZeta[0] * field.Jac(i,j,k);
                    field.ZeY(i,j,k) = cZeta[1] * field.Jac(i,j,k);
                    field.ZeZ(i,j,k) = cZeta[2] * field.Jac(i,j,k);
                }
            }
        }
        // --- 关键修复：边界层的一阶偏置差分 (补全最外面一圈 0 和 M-1) ---
        // 为了防止 updateField 在边界点溢出，最外层那一个点可以用一阶差分
        applyBoundaryMetrics(field);
        auditMetricIdentity(field);
    }

    void Mesh::auditMetricIdentity(SF::Field& field) {
        double maxMetricDivergence = 0.0;
        double maxCofactorMagnitude = 0.0;
        int maxI = 0;
        int maxJ = 0;
        int maxK = 0;
        int maxComponent = 0;
        const int ng = field.NG();
        for (int k = ng; k < field.NZ() + ng; ++k) {
            for (int j = ng; j < field.NY() + ng; ++j) {
                for (int i = ng; i < field.NX() + ng; ++i) {
                    for (int component = 0; component < 3; ++component) {
                        double divergence = 0.0;
                        if (SF::Math::isDirectionActive(SF::Math::XI)) {
                            double fp[4], fm[4];
                            SF::Math::faceMetrics(field, i, j, k,
                                                  SF::Math::XI, fp);
                            SF::Math::faceMetrics(field, i - 1, j, k,
                                                  SF::Math::XI, fm);
                            maxCofactorMagnitude = std::max({
                                maxCofactorMagnitude,
                                std::abs(fp[component]),
                                std::abs(fm[component])});
                            divergence += fp[component] - fm[component];
                        }
                        if (SF::Math::isDirectionActive(SF::Math::ETA)) {
                            double fp[4], fm[4];
                            SF::Math::faceMetrics(field, i, j, k,
                                                  SF::Math::ETA, fp);
                            SF::Math::faceMetrics(field, i, j - 1, k,
                                                  SF::Math::ETA, fm);
                            maxCofactorMagnitude = std::max({
                                maxCofactorMagnitude,
                                std::abs(fp[component]),
                                std::abs(fm[component])});
                            divergence += fp[component] - fm[component];
                        }
                        if (SF::Math::isDirectionActive(SF::Math::ZETA)) {
                            double fp[4], fm[4];
                            SF::Math::faceMetrics(field, i, j, k,
                                                  SF::Math::ZETA, fp);
                            SF::Math::faceMetrics(field, i, j, k - 1,
                                                  SF::Math::ZETA, fm);
                            maxCofactorMagnitude = std::max({
                                maxCofactorMagnitude,
                                std::abs(fp[component]),
                                std::abs(fm[component])});
                            divergence += fp[component] - fm[component];
                        }
                        if (!std::isfinite(divergence)) {
                            std::cerr << "[SF FATAL] non-finite conservative "
                                      << "metric identity residual at ("
                                      << i << "," << j << "," << k
                                      << "), component=" << component
                                      << std::endl;
                            std::exit(1);
                        }
                        const double absDivergence = std::abs(divergence);
                        if (absDivergence > maxMetricDivergence) {
                            maxMetricDivergence = absDivergence;
                            maxI = i;
                            maxJ = j;
                            maxK = k;
                            maxComponent = component;
                        }
                    }
                }
            }
        }
        std::ostringstream oss;
        const double tolerance =
            65536.0 * std::numeric_limits<double>::epsilon()
            * std::max(maxCofactorMagnitude,
                       std::numeric_limits<double>::min());
        oss << "max |D_xi C_xi + D_eta C_eta + D_zeta C_zeta|="
            << maxMetricDivergence << " at (" << maxI << "," << maxJ
            << "," << maxK << "), component=" << maxComponent
            << ", tolerance=" << tolerance;
        SF::broadcast("Metric identity audit: ", oss.str());
        if (!std::isfinite(tolerance) || tolerance <= 0.0 ||
            maxMetricDivergence > tolerance) {
            std::cerr
                << "[SF FATAL] discrete GCL/metric identity is violated. "
                << "SCMM geometry must satisfy the discrete divergence "
                << "theorem before flux evaluation; metric repair or "
                << "averaging is forbidden. " << oss.str() << std::endl;
            std::exit(1);
        }
    }
    
    // 边界条件度规
    void Mesh::applyBoundaryMetrics(SF::Field& field) {
        // 强制让最外层一圈的度规等于紧邻的内层点
        // 这样 updateField 无论如何都不会读到 0
        // --- 1. 补齐 I 方向 (X) 的最外壳 (左右两端) ---
        int MX = field.MX() ;
        int MY = field.MY() ;
        int MZ = field.MZ() ;
        for (int k = 0; k < MZ; ++k) {
            for (int j = 0; j < MY; ++j) {
                // 左边界点 (0) 拷贝自其内侧邻居 (1)
                copyMetricCell(field, 0, j, k, 1, j, k);
                // 右边界点 (MX-1) 拷贝自其内侧邻居 (MX-2)
                copyMetricCell(field, MX - 1, j, k, MX - 2, j, k);
            }
        }

        // --- 2. 补齐 J 方向 (Y) 的最外壳 (上下两端) ---
        for (int k = 0; k < MZ; ++k) {
            for (int i = 0; i < MX; ++i) {
                // 下边界点 (j=0) 拷贝自其内侧邻居 (j=1)
                copyMetricCell(field, i, 0, k, i, 1, k);
                // 上边界点 (j=MY-1) 拷贝自其内侧邻居 (j=MY-2)
                copyMetricCell(field, i, MY - 1, k, i, MY - 2, k);
            }
        }

        // --- 3. 补齐 K 方向 (Z) 的最外壳 (前后两端) ---
        for (int j = 0; j < MY; ++j) {
            for (int i = 0; i < MX; ++i) {
                // 前边界点 (k=0) 拷贝自其内侧邻居 (k=1)
                copyMetricCell(field, i, j, 0, i, j, 1);
                // 后边界点 (k=MZ-1) 拷贝自其内侧邻居 (k=MZ-2)
                copyMetricCell(field, i, j, MZ - 1, i, j, MZ - 2);
            }
        }
    }

    // 度规单元赋予
    void Mesh::copyMetricCell(SF::Field& field, int di, int dj, int dk, int si, int sj, int sk) { 
        // 拷贝 Jacobian (此时你存的是 1/J 还是 J 都要拷) 
        field.Jac(di, dj, dk) = field.Jac(si, sj, sk); 
        field.XiX(di, dj, dk) = field.XiX(si, sj, sk); 
        field.XiY(di, dj, dk) = field.XiY(si, sj, sk); 
        field.XiZ(di, dj, dk) = field.XiZ(si, sj, sk); 
        field.EtX(di, dj, dk) = field.EtX(si, sj, sk); 
        field.EtY(di, dj, dk) = field.EtY(si, sj, sk); 
        field.EtZ(di, dj, dk) = field.EtZ(si, sj, sk); 
        field.ZeX(di, dj, dk) = field.ZeX(si, sj, sk); 
        field.ZeY(di, dj, dk) = field.ZeY(si, sj, sk); 
        field.ZeZ(di, dj, dk) = field.ZeZ(si, sj, sk); 
        }

    void Mesh::generateGhostCells(SF::Field& field) {
        int MX = field.MX(); int MY = field.MY(); int MZ = field.MZ();
        int NX = field.NX(); int NY = field.NY(); int NZ = field.NZ();
        int ng = field.NG();

        // --- 1. 预处理：构建临时查找表 ---
        std::unordered_map<int, std::vector<std::string>> tagMap;
        for (auto const& [name, indices] : field.getAllSets()) {
            for (int idx : indices) tagMap[idx].push_back(name);
        }
        // --- 1. 处理 I 方向 ---
        for (int k = 0; k < MZ; ++k) {
            for (int j = 0; j < MY; ++j) {
                // 向左外推 (从 ng-1 到 0)
                for (int i = ng - 1; i >= 0; --i) {
                    extrapolateCoord(field, i, j, k, i + 1, j, k, i + 2, j, k);
                    syncBoundaryTags(field, tagMap, ng, j, k, i, j, k);
                }
                // 向右外推 (从 MX-ng 到 MX-1)
                for (int i = MX - ng; i < MX; ++i) {
                    extrapolateCoord(field, i, j, k, i - 1, j, k, i - 2, j, k);
                    syncBoundaryTags(field, tagMap, NX + ng - 1, j, k, i, j, k);
                }
            }
        }
        // --- 2. 处理 J 方向 (逻辑同上) ---
        for (int k = 0; k < MZ; ++k) {
            for (int i = 0; i < MX; ++i) {
                for (int j = ng - 1; j >= 0; --j) { // 往下推
                    extrapolateCoord(field, i, j, k, i, j + 1, k, i, j + 2, k);
                    syncBoundaryTags(field, tagMap, i, ng, k, i, j, k);
                }
                for (int j = MY - ng; j < MY; ++j) { // 往上推
                    extrapolateCoord(field, i, j, k, i, j - 1, k, i, j - 2, k);
                    syncBoundaryTags(field, tagMap, i, NY + ng - 1, k, i, j, k);
                }
            }
        }
        // --- 3. 处理 K 方向 (逻辑同上) ---
        for (int j = 0; j < MY; ++j) {
            for (int i = 0; i < MX; ++i) {
                for (int k = ng - 1; k >= 0; --k) { // 往侧推
                    extrapolateCoord(field, i, j, k, i, j, k + 1, i, j, k + 2);
                    syncBoundaryTags(field, tagMap, i, j, ng, i, j, k);
                }
                for (int k = MZ - ng; k < MZ; ++k) { // 往另一侧推
                    extrapolateCoord(field, i, j, k, i, j, k - 1, i, j, k - 2);
                    syncBoundaryTags(field, tagMap, i, j, NZ + ng - 1, i, j, k);
                }
            }
        }
    }

    void Mesh::initializeGhostStateFromInterior(SF::Field& field) {
        int ng = field.NG();
        int iMin = ng;
        int jMin = ng;
        int kMin = ng;
        int iMax = field.NX() + ng - 1;
        int jMax = field.NY() + ng - 1;
        int kMax = field.NZ() + ng - 1;

        for (int k = 0; k < field.MZ(); ++k) {
            for (int j = 0; j < field.MY(); ++j) {
                for (int i = 0; i < field.MX(); ++i) {
                    bool isGhost = i < iMin || i > iMax ||
                                   j < jMin || j > jMax ||
                                   k < kMin || k > kMax;
                    if (!isGhost) continue;

                    int ri = std::max(iMin, std::min(i, iMax));
                    int rj = std::max(jMin, std::min(j, jMax));
                    int rk = std::max(kMin, std::min(k, kMax));
                    for (int v = 0; v < 5; ++v) {
                        field(i, j, k, v) = field(ri, rj, rk, v);
                    }
                }
            }
        }
    }

    
    /// @brief 线性外推坐标
    /// @param targetIdx (i,j,k) 待生成的虚胞坐标
    /// @param nearIdx (ni,nj,nk) 靠近边界的第一个点（供体1）
    /// @param farIdx (fi,fj,fk) 远离边界的第二个点（供体2）
    void Mesh::extrapolateCoord(SF::Field& field, 
                                int i,  int j,  int k, 
                                int ni, int nj, int nk, 
                                int fi, int fj, int fk) 
    {
        // 线性外推公式: P_target = P_near + (P_near - P_far) = 2 * P_near - P_far
        field.X(i, j, k) = 2.0 * field.X(ni, nj, nk) - field.X(fi, fj, fk);
        field.Y(i, j, k) = 2.0 * field.Y(ni, nj, nk) - field.Y(fi, fj, fk);
        field.Z(i, j, k) = 2.0 * field.Z(ni, nj, nk) - field.Z(fi, fj, fk);
    }

    void Mesh::syncBoundaryTags(
        SF::Field& field, const std::unordered_map<int, std::vector<std::string>>& tagMap,
        int ri, int rj, int rk, int gi, int gj, int gk)
    {
        int ng = field.NG();
        ri = std::max(ng, std::min(ri, field.NX() + ng - 1));
        rj = std::max(ng, std::min(rj, field.NY() + ng - 1));
        rk = std::max(ng, std::min(rk, field.NZ() + ng - 1));

        int rIdx = field.getIdx(ri, rj, rk);
        int gIdx = field.getIdx(gi, gj, gk);

        // 检查实胞是否有标签
        auto it = tagMap.find(rIdx);
        if (it != tagMap.end()) {
            auto& allSets = field.getAllSets();
            for (const auto& tagName : it->second) {
                allSets[tagName].push_back(gIdx);
            }
        }
    }    
}
