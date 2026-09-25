/// @file SF_IO.h
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <string>
#include <vector>
#include <optional>
#include <map>
#include "core/model/SF_model.h"
#include "SF_valueTypes.h"
#include "SF_multiphase.h"
#include "SF_config.h"
#include "SF_caseConfig.h"

// 前向声明 (避免循环依赖 SF_IO ↔ SF_meshGen)
namespace SF { struct MeshParameters; }

namespace SF {

/// @brief 原生 case 的 typed section 集合。
///
/// 每个 slot 对应一个已声明的语义对象（solvers/*.yaml、models/*.yaml）。
/// 这里不构造 controlDict / outputDict / solverProperties / fvSchemes /
/// fvSolution / constant/* 形状的中间字典：那些名字只描述 OpenFOAM 的磁盘
/// 布局而不是数值语义，把它们当作解码入口会重新引入一层翻译。
struct NativeCaseSections {
    Model::Parameters runtime{Model::Parameters::object()};
    Model::Parameters output{Model::Parameters::object()};
    Model::Parameters parallel{Model::Parameters::object()};
    Model::Parameters algorithm{Model::Parameters::object()};
    Model::Parameters numerics{Model::Parameters::object()};
    Model::Parameters thermoDynamics{Model::Parameters::object()};
    Model::Parameters phaseSystem{Model::Parameters::object()};
    Model::Parameters phaseChange{Model::Parameters::object()};
    Model::Parameters turbulence{Model::Parameters::object()};
    Model::Parameters ibm{Model::Parameters::object()};
    Model::Parameters ilw{Model::Parameters::object()};
    Model::Parameters gravity{Model::Parameters::object()};
    Model::Parameters mrf{Model::Parameters::object()};
    Model::Parameters wallHeat{Model::Parameters::object()};
    /// field 的语义名字就是它的解码入口。`boundaries` 已按 mesh patch set
    /// 展开 `default`，因此解码阶段不再做目录名或文件路径推断。
    std::vector<Model::FieldDescriptor> fields;

    bool hasAlgorithm = false;
    bool hasOutput = false;
    bool hasParallel = false;
    bool hasThermoDynamics = false;
    bool hasPhaseSystem = false;
    bool hasPhaseChange = false;
    bool hasTurbulence = false;
    bool hasIbm = false;
    bool hasIlw = false;
    bool hasGravity = false;
    bool hasMrf = false;
    bool hasWallHeat = false;

    /// @brief 按语义名字查 field；未声明返回 nullptr。
    const Model::FieldDescriptor* field(const std::string& name) const {
        for (const auto& f : fields) {
            if (f.name == name) return &f;
        }
        return nullptr;
    }
};

/// @brief 原生 case reader facade。
///
/// 只负责把已解析的原生语义对象绑定成 CaseConfig；网格加载和结果输出由
/// 独立模块负责。绑定过程不产生第二套文档形状的 authority。
class CaseAdapter {
public:
    explicit CaseAdapter(const std::string& caseFilePath);

    /// @brief 将原生语义对象直接绑定成 CaseConfig。
    CaseConfig build(const Model::Description& description);
    /// @brief 已绑定的 typed section 集合。
    ///
    /// OWNERSHIP: 只读观察入口，供 explain 与回归测试确认"语义对象 -> typed
    /// section"是直接映射；它不构成第二份配置 authority。
    const NativeCaseSections& sections() const { return sections_; }
    bool fileExists(const std::string& path) const;
    bool isDirectoryPath(const std::string& path) const;

    // ── 初始化 (解析 case 目录, 但不立即加载网格/大文件) ──
    /// @brief 从 build() 存入的 typed section 解码完整强类型配置。
    /// @return 成功时返回 CaseConfig；输入缺失时返回 std::nullopt。
    std::optional<CaseConfig> decodeNativeCase();

private:
    NativeCaseSections sections_;
    bool native_ = false;
    std::string caseFilePath_;
    std::string caseDir_;
    std::string caseName_;
    std::string jobName_;       ///< 来自 #job 词条, 用于PVD/VTM/VTS文件命名
    /// @brief Legacy 兼容标签（densityBase/pressureBase）；只用于翻译成显式
    ///        equation request，不是求解器身份。
    std::string legacyFlowLabel_;
    std::string outputDir_;
    /// mesh generator 的相对路径；由 mesh/mesh.yaml 的 generator 提供。
    std::string meshGenerator_;

    std::vector<std::string> meshFiles_;
    std::vector<std::string> ibmGeometryFiles_;
    std::vector<WallHeatSetting> wallHeatBoundarySettings_;
    std::vector<BCSetting<double>> temperatureOutputBC_;
    std::vector<BCSetting<double>> densityOutputBC_;
    std::vector<BCSetting<Vector3>> velocityOutputBC_;
    std::vector<BCSetting<double>> pressureOutputBC_;
    bool ibmOutputEnabled_ = false;
    FDM::IBMMethod ibmMethod_ = FDM::IBMMethod::Ghost;
    FDM::IBMForcingConfig ibmForcingConfig_;
    Physics::Multiphase::MultiPhaseConfig multiPhaseConfig_;
    FDM::PressureCorrectionConfig solverProperties_;
    CaseConfig caseConfig_;
    /// @brief case-local source 选择串（数值块 + 边界 wallHeat 追加）。
    ///
    /// 它只服务本 case 的解码过程，因此是 adapter 成员而不是 parser 全局
    /// 暂存值；最终选择由 FDM::parseSourceKinds(sourceScheme_) 转成 typed
    /// SourceConfig。
    std::string sourceScheme_;

    // ── 懒加载状态标记 ──
    bool icLoaded_ = false;
    bool bcLoaded_ = false;
    bool scLoaded_ = false;
    bool multiPhaseLoaded_ = false;
    bool solverPropertiesLoaded_ = false;

    void buildCaseConfig();
    ResultWriterConfig resultWriterConfig() const;
    void decodeNativeSections();
    void decodeRuntimeSection();
    void decodeOutputSection();
    void decodeNumericsSection();
    void decodeParallelSection();
    void decodeAlgorithmSection();
    void decodeThermoDynamicsSection();
    void decodePhaseSystemSection();
    void decodePhaseChangeSection();
    void decodeTurbulenceSection();
    void decodeIBMSection();
    void decodeILWSection();
    void decodeFields();
    void decodePhiField(const Model::FieldDescriptor& field);
    void decodeVectorField(const Model::FieldDescriptor& field,
                           std::vector<BCSetting<Vector3>>& internal,
                           std::vector<BCSetting<Vector3>>& boundary,
                           Vector3* uniformValue);
    void decodeScalarField(const Model::FieldDescriptor& field,
                           std::vector<BCSetting<double>>& internal,
                           std::vector<BCSetting<double>>& boundary,
                           double* uniformValue);
    void decodeTemperatureField(const Model::FieldDescriptor& field);
    void decodeSourceSettings();
    void decodeGravitySection();
    void decodeMRFSection();
    void decodeWallHeatSection();
    void resolvePhaseChangeWallBoilingHeatFlux();
    void addIBMGeometryFile(const std::string& fileName,
                            bool explicitPath);

};

} // namespace SF
