/// @file SF_IO.cpp
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_compatibility.h"
#include "app/application/model/SF_configParser.h"
#include "core/interfaces/SF_log.h"
#include "SF_cellFaceMesh.h"
#include "core/mesh/SF_dimension.h"
#include "SF_meshGen.h"
#include "SF_phaseChange.h"
#include "private/SF_casePath.h"
#include <fstream>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <regex>
#include <initializer_list>
#include <unordered_map>
#include <sys/stat.h>

namespace SF {

using namespace IOPrivate;

 
// ============================================================
//  构造 & 顶层加载
// ============================================================

CaseAdapter::CaseAdapter(const std::string& caseFilePath) : caseFilePath_(caseFilePath) {
    if (logFromThisProcess()) {
        std::cout << std::string(60, '=') << "\n"
                  << "Welcome To SonicSolver\n"
                  << "*copyright by Li Pengfei*\n"
                  << std::string(60, '=') << "\n";
    }

    const std::string cleanPath = stripTrailingSlash(caseFilePath);
    if (isDirectoryPath(cleanPath)) {
        caseDir_ = cleanPath.empty() ? "." : cleanPath;
        caseName_ = basePathNameNoExtension(caseDir_);
        caseFilePath_ = caseDir_;
    } else {
        // 分离目录和文件名
        size_t sep = cleanPath.find_last_of("/\\");
        if (sep != std::string::npos) {
            caseDir_  = cleanPath.substr(0, sep);
            caseName_ = cleanPath.substr(sep + 1);
        } else {
            caseDir_  = ".";
            caseName_ = cleanPath;
        }
        // 去掉扩展名
        size_t dot = caseName_.find_last_of('.');
        if (dot != std::string::npos) caseName_ = caseName_.substr(0, dot);
    }
}

void CaseAdapter::buildCaseConfig() {
    caseConfig_.caseDir = caseDir_;
    caseConfig_.caseName = caseName_;
    caseConfig_.jobName = jobName_;
    caseConfig_.meshParameterFile = meshGenerator_;
    caseConfig_.meshFiles = meshFiles_;
    caseConfig_.ibmGeometryFiles = ibmGeometryFiles_;
    caseConfig_.multiPhase = multiPhaseConfig_;
    caseConfig_.multiPhaseEnabled = multiPhaseLoaded_;

    // typed solver 配置由 decode*() 逐语义块写入（runtime/numerics/
    // algorithm/turbulence/ILW/IBM/sources/fields）。这里只补齐跨语义块的
    // 派生量，不构造任何"typed → legacy globals → typed"往返。
    FDM::SolverConfig& solver = caseConfig_.solver;
    solver.numerics.startTime = caseConfig_.time.startTime;
    solver.ibm.method = ibmMethod_;
    solver.ibm.forcing = ibmForcingConfig_;
    // transport 常量解码完成后再投影到消费它的 typed 位置。
    solver.boundaries.thermalDynamicViscosity = solver.numerics.dynamicViscosity;
    solver.boundaries.thermalPrandtl = solver.numerics.prandtl;
    solver.turbulence.laminarDynamicViscosity =
        solver.numerics.dynamicViscosity;
    if (!solver.numerics.termRecipesDeclared) {
        // scheme 选择是 typed authority，recipe 由它派生。
        solver.numerics.recipes.convection = FDM::builtInConvectionRecipe(
            solver.numerics.convection, solver.numerics.flux);
        solver.numerics.recipes.diffusion = solver.numerics.viscousEnabled
            ? std::optional<FDM::TermRecipe>(
                  FDM::builtInDiffusionRecipe(solver.numerics.viscous))
            : std::nullopt;
    }
    solver.numerics.recipes.time = solver.numerics.timeRecipe;
    solver.sources.enabled = FDM::parseSourceKinds(sourceScheme_);
    if (multiPhaseLoaded_
        && multiPhaseConfig_.phaseChange.enabled
        && !multiPhaseConfig_.phaseChange.wallBoiling.empty()) {
        solver.sources.wallHeat =
            multiPhaseConfig_.phaseChange.wallBoiling;
        if (std::find(solver.sources.enabled.begin(),
                      solver.sources.enabled.end(),
                      FDM::SourceKind::WallHeat)
            == solver.sources.enabled.end()) {
            solver.sources.enabled.push_back(
                FDM::SourceKind::WallHeat);
        }
    }
    if (solverPropertiesLoaded_) {
        solver.pressure = solverProperties_;
        solver.pressure.relaxation =
            solverProperties_.coupling.pressureRelaxation;
        solver.pressure.velocityRelaxation =
            solverProperties_.coupling.momentumRelaxation;
    }
    FDM::validateNumericsConfig(solver.numerics);

    caseConfig_.output = resultWriterConfig();
}

ResultWriterConfig CaseAdapter::resultWriterConfig() const {
    ResultWriterConfig config;
    config.caseDir = caseDir_;
    config.caseName = caseName_;
    config.jobName = jobName_;
    config.outputDir = outputDir_;
    config.densityBoundary = densityOutputBC_;
    config.velocityBoundary = velocityOutputBC_;
    config.pressureBoundary = pressureOutputBC_;
    config.temperatureBoundary = temperatureOutputBC_;
    config.ibmOutputEnabled = ibmOutputEnabled_;
    return config;
}

} // namespace SF
