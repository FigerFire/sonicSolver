#pragma once

/// @file SF_interfaceModel.h
/// @brief 锐界面表示模型的统一物理接口与工厂。

#include "core/interfaces/SF_equationCoupling.h"
#include "SF_phaseProperties.h"

#include <memory>
#include <vector>

namespace SF::Physics::Multiphase {
class LevelSetField;
}

namespace SF::Physics::InterfaceModels {

/// @brief 界面表示类型；运动/平衡模型由 PhaseSystem 另行决定。
enum class Representation {
    LevelSet,
    VolumeOfFluid,
    CoupledLevelSetVolumeOfFluid
};

/// @brief 求解器无关的锐界面状态与物理贡献接口。
///
/// InterfaceModel 拥有 alpha/phi 和几何缓存，但不拥有主流守恒状态、不选择
/// 密度基/压力基算法，也不提交主流时间层。
class Model : public FDM::ITransportModel,
              public FDM::IInterfaceJumpCondition {
public:
    virtual ~Model() = default;

    /// @brief 返回界面表示类型。
    virtual Representation representation() const noexcept = 0;
    /// @brief 返回模型是否完成初始化。
    virtual bool initialized() const noexcept = 0;
    /// @brief 返回创建该界面模型的只读物理配置。
    virtual const Multiphase::MultiPhaseConfig& config() const noexcept = 0;
    /// @brief 根据网格和 case 初值创建 canonical 界面状态。
    virtual void initialize(const Field& field) = 0;
    /// @brief 从 canonical 界面状态刷新几何、物性和诊断缓存。
    virtual void refresh(const Field& field) = 0;
    /// @brief 分阶段几何计算，调用方负责两阶段之间的 normal stencil 发布。
    virtual void computeNormals(const Field& field) = 0;
    virtual void refreshFromNormals(const Field& field) = 0;
    /// @brief 向主方程装配表面张力等界面贡献。
    virtual void addSourceTerms(Field& field, Residual& residual) const = 0;
    /// @brief 返回 canonical 界面标量，供状态注册和输出使用。
    virtual ScalarField& primaryScalar() = 0;
    virtual const ScalarField& primaryScalar() const = 0;
    /// @brief Level Set 专用状态访问；非 Level Set 模型返回 nullptr。
    virtual Multiphase::LevelSetField* levelSetState() noexcept {
        return nullptr;
    }
    virtual const Multiphase::LevelSetField* levelSetState() const noexcept {
        return nullptr;
    }
    /// @brief 返回相指标和流动质量诊断。
    virtual const Multiphase::ConservationDiagnostics& conservation() const = 0;
    /// @brief 接收 runtime 对 owner 物理积分的全局 SUM；只更新诊断快照。
    virtual void publishConservation(double phaseIndicator,double flowMass) = 0;
};

/// @brief 按 phaseProperties 创建已配置但尚未初始化的界面模型。
/// @throws std::runtime_error 请求尚未实现的 VOF/CLSVOF 或未知表示时失败。
std::unique_ptr<Model> makeInterfaceModel(
    const Multiphase::MultiPhaseConfig& config);

} // namespace SF::Physics::InterfaceModels
