/// @file SF_multiphase.h
/// @brief 多相配置、物性与模型一致性校验实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "core/interfaces/SF_transportModel.h"
#include "SF_phaseProperties.h"
#include "SF_scalarField.h"
#include "SF_phaseChange.h"

#include <functional>
#include <vector>

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief legacy mixture/thermal 物理模型调度入口。
///
/// 锐界面状态由 InterfaceModels::Model 独立拥有，本类不接受 Level Set。
class MultiPhaseModel : public FDM::ITransportModel {
public:
    MultiPhaseModel();
    ~MultiPhaseModel();
    MultiPhaseModel(const MultiPhaseModel&) = delete;
    MultiPhaseModel& operator=(const MultiPhaseModel&) = delete;
    MultiPhaseModel(MultiPhaseModel&&) noexcept;
    MultiPhaseModel& operator=(MultiPhaseModel&&) noexcept;

    /// @brief 写入并校验配置。
    /// @param config IO 层解析出的多相配置。
    void configure(const MultiPhaseConfig& config);

    /// @brief 注入求解器拥有的辅助标量边界离散设置。
    void configureBoundaryNumerics(bool ilwEnabled, int ilwAccuracyOrder) {
        scalarILWEnabled_ = ilwEnabled;
        scalarILWAccuracyOrder_ = ilwAccuracyOrder;
    }

    /// @brief 当前是否启用多相模型。
    bool enabled() const { return config_.enabled; }

    /// @brief 当前是否使用 mixture-like 单速度 alpha 模型。
    bool isMixture() const { return enabled() && isMixtureType(config_.type); }

    /// @brief 当前是否使用单温度热输运模型。
    bool isThermal() const { return enabled() && isThermalType(config_.type); }

    /// @brief legacy 状态是否已完成初始化。
    bool initialized() const { return initialized_; }

    /// @brief 访问只读配置。
    const MultiPhaseConfig& config() const { return config_; }

    /// @brief 访问 alpha 标量场。
    ScalarField& alpha() { return alpha_; }
    const ScalarField& alpha() const { return alpha_; }

    /// @brief 时间积分器推进的守恒相质量（单位体积内被跟踪相质量）。
    ScalarField& phaseMass() { return phaseMass_; }
    const ScalarField& phaseMass() const { return phaseMass_; }

    /// @brief 访问单温度标量场。
    ScalarField& temperature() { return temperature_; }
    const ScalarField& temperature() const { return temperature_; }

    /// @brief 单温度标量是否已初始化。
    bool hasTemperature() const { return hasTemperature_; }

    /// @brief 查询相变能量源缓存是否已更新。
    bool hasPhaseChangeSource() const { return hasPhaseChangeSource_; }

    /// @brief 访问相变能量源缓存。
    double phaseChangeEnergySource(int i, int j, int k) const {
        return phaseChangeEnergySource_[(size_t)alpha_.getIdx(i, j, k)];
    }

    /// @brief 查询 mixture 物性缓存是否已更新。
    bool hasMixtureProperties() const { return hasMixtureProperties_; }

    /// @brief 访问 mixture 密度缓存。
    double mixtureDensity(int i, int j, int k) const {
        return mixtureDensity_[(size_t)alpha_.getIdx(i, j, k)];
    }

    /// @brief 访问 mixture 动力黏度缓存。
    double mixtureViscosity(int i, int j, int k) const {
        return mixtureViscosity_[(size_t)alpha_.getIdx(i, j, k)];
    }

    /// @brief 访问守恒诊断缓存。
    const ConservationDiagnostics& conservation() const;

    /// @brief 根据当前配置初始化具体多相模型。
    /// @param field 已加载网格和 set 的流场。
    void initialize(const Field& field);

    /// @brief 初始化 mixture-like alpha 标量场。
    /// @param field 已加载网格和 set 的流场。
    void initializeAlpha(const Field& field);

    /// @brief 推进 legacy thermal 物理时间步；mixture 由主积分器推进。
    /// @param field 提供速度和网格度规的主流场。
    /// @param dt 主求解器本步时间步长。
    void advance(Field& field, double dt);

    /// @brief 应用多相辅助场边界条件。
    /// @param field 主流场, 提供网格 set 和边界定位。
    void applyAuxiliaryBoundaryConditions(const Field& field);

    /// @brief 根据当前辅助场更新几何、物性和守恒诊断。
    /// @param field 主流场。
    void updateFlowCoupling(const Field& field);

    ///  从守恒量 rhoE 和 alpha 反推混合物温度。
    ///
    /// 使用混合物 EOS：T = (rhoE - 0.5.rho.|u|.2) / (.alpha.*.rho_l.*Cp_l + (1-.alpha).*.rho_v.*Cp_v)。
    ///  field 守恒量场（只读），从中读取 rho/rhoU/rhoV/rhoW/rhoE。
    void refreshDerivedFromConserved(const Field& field);

    /// @brief 冷启动闭合：以 alpha、T 和当前速度构造 rho、rhoE。
    /// @param field 待写入的密度基守恒状态。
    void initializeConservedFromPrimitive(Field& field);

    /// @brief 从当前统一状态装配相质量与能量 RHS，不提交时间层。
    ///
    /// 此方法供 assembleRHS 回调在 RK4 各 stage 内调用：
    /// 1. 从 rhoE 和相质量推导 alpha、温度
    /// 2. 调用 computeRates 获取候选 mdot（legacy 非相变主线）
    /// 3. 总能量不添加重复潜热源
    /// 4. 将相质量的对流、扩散和相变源写入 phaseMassRHS
    ///
    /// @param field 守恒量场（只读状态，累加 Source(E)）。
    /// @param phaseMassRHS 相质量时间导数，与场同尺寸。
    /// @param dt        时间步长
    void assembleIntegratedRHS(Field& field, Residual& residual,
                               std::vector<double>& phaseMassRHS,
                               double dt);

    /// @brief 时间层提交后刷新派生状态和守恒诊断。
    void commitTimeLevel(const Field& field);


    /// @brief 将多相物理源项装配到主守恒方程源项容器。
    /// @param field 主流场, 其 Field::Source 将被累加。
    void addSourceTerms(Field& field, Residual& residual) const;

    /// @brief ITransportModel: 应用多相辅助场边界。
    void applyBoundary(const Field& field) override;

    /// @brief ITransportModel: 更新多相物性缓存；不推进辅助场。
    void correct(const Field& field, double dt) override;

    /// @brief ITransportModel: 返回多相局部动力黏度。
    double dynamicViscosity(const Field& field,
                            int i,
                            int j,
                            int k,
                            double laminarMu) const override;

private:
    MultiPhaseConfig config_;
    ScalarField alpha_;
    ScalarField phaseMass_;
    ScalarField temperature_;
    std::vector<double> mixtureDensity_;
    std::vector<double> mixtureViscosity_;
    std::vector<double> phaseChangeEnergySource_;
    PhaseChange::PhaseChangeDiagnostics phaseChangeDiagnostics_;
    ConservationDiagnostics diagnostics_;
    bool initialized_ = false;
    bool hasMixtureProperties_ = false;
    bool hasTemperature_ = false;
    bool hasPhaseChangeSource_ = false;
    bool scalarILWEnabled_ = false;
    int scalarILWAccuracyOrder_ = 0;

    void updateMixtureProperties(const Field& field);
    void initializeTemperature(const Field& field);
    void advanceTemperature(Field& field, double dt);
    void updateAlphaFromPhaseMass(const Field& field);
    void updatePhaseMassFromAlpha(const Field& field);
    void applyPhaseChangeClosure(Field& field, double dt);
    void updateConservationDiagnostics(const Field& field,
                                       bool resetReference);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
