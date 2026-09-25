#pragma once

/// @file SF_caseConfig.h
/// @brief 完整 case 的强类型只读配置。

#include "core/config/SF_config.h"
#include "core/config/SF_equationComposition.h"
#include "core/config/SF_runControl.h"
#include "SF_multiphase.h"
#include "SF_resultWriter.h"

#include <array>
#include <memory>
#include "core/model/SF_model.h"
#include <string>
#include <vector>

namespace SF {

/// @brief MPI 与网格分区配置。
struct CaseParallelConfig {
    bool enabled = false;
    int processCount = 1;
    bool automaticPartition = true;
    std::array<int, 3> partitionSplit{1, 1, 1};
    double haloTolerance = 1.0e-8;
};

/// @brief reader 完成解析后交给 application 的统一 case 值对象。
///
/// 原生 ModelDescription 或旧格式兼容适配器只负责构造本对象；solver、mesh 和
/// writer 均消费本值对象，不查询 parser 全局变量或 IO 对象内部状态。
struct CaseConfig {
    std::shared_ptr<const Model::Description> modelDescription;
    std::string caseDir;
    std::string caseName;
    std::string jobName;
    /// @brief Legacy 兼容输入标签（`algorithm.yaml: type:
    ///        densityBase|pressureBase`）。
    ///
    /// 它只在 composition root 被翻译成**显式 equation request**
    /// （conservative preset 或 pressure-constraint 方程族），既不进入
    /// ResolvedSimulationSystem，也不参与任何 runtime dispatch。
    /// 空字符串表示 case 没有声明该标签。
    std::string compatFlowLabel;
    /// @brief case 是否显式声明了 pressure-coupling preset（algorithm.yaml）。
    bool pressureCouplingDeclared = false;
    std::string meshParameterFile;
    std::vector<std::string> meshFiles;
    std::vector<std::string> ibmGeometryFiles;

    FDM::SolverConfig solver;
    /// User-facing equations/models/algorithms declaration.  The legacy
    /// solver config remains a compatibility input until all runtimes consume
    /// the compiled system directly.
    EquationCompositionConfig composition;
    Physics::Multiphase::MultiPhaseConfig multiPhase;
    bool multiPhaseEnabled = false;

    Time::RunControl time;
    CaseParallelConfig parallel;
    bool writeInitial = true;
    bool createMesh = false;

    ResultWriterConfig output;
};

} // namespace SF
