#pragma once

/// @file SF_equationSystem.h
/// @brief 可注册到求解工作流的湍流输运方程系统。

#include "SF_config.h"
#include "SF_phaseSystem.h"

#include <string>
#include <vector>

namespace SF::Turbulence {

/// @brief 一个守恒标量输运方程的状态和装配系数。
struct TransportEquationState {
    ScalarField variable;
    ScalarField previousConserved;
    ScalarField diffusivity;
    ScalarField explicitSource;
    ScalarField implicitSink;
};

/// @brief 一个 Eulerian 相拥有的湍流状态。
struct PhaseEquationState {
    size_t phaseIndex = 0;
    TransportEquationState kineticEnergy;
    TransportEquationState dissipation;
    TransportEquationState specificDissipation;
    ScalarField eddyViscosity;
    ScalarField wallDistance;
};

/// @brief 湍流物理只提供状态与系数，线性方程由 workflow 统一求解。
class EquationSystem {
public:
    explicit EquationSystem(FDM::TurbulenceConfig config);

    /// @brief 按显式 phase 列表分配并初始化湍流主状态。
    void initialize(
        const Physics::PhaseSystems::PhaseSystem& phaseSystem);

    bool active() const { return active_; }
    bool hasTransportEquations() const;
    std::string description() const;

    /// @brief 应用湍流标量物理边界；halo 由 Solver Algorithm 统一请求。
    void applyBoundary(
        const Physics::PhaseSystems::PhaseSystem& phaseSystem);
    /// @brief 从当前相状态计算 mu_t、生产耗散和扩散系数。
    void prepare(
        const Physics::PhaseSystems::PhaseSystem& phaseSystem);
    /// @brief 成功完成物理时间步后提交旧时间层。
    void commit(
        const Physics::PhaseSystems::PhaseSystem& phaseSystem);
    /// @brief 校验所有被求解湍流状态为有限正值。
    void validate(
        const Physics::PhaseSystems::PhaseSystem& phaseSystem,
        const std::string& stage) const;

    double momentumEddyViscosity(size_t phase, int cell) const;
    double energyEddyDiffusivity(size_t phase, int cell) const;
    double sourceTimeStep(
        const Physics::PhaseSystems::PhaseSystem& phaseSystem,
        double sourceCfl) const;

    std::vector<PhaseEquationState>& states() { return states_; }
    const std::vector<PhaseEquationState>& states() const { return states_; }
    const FDM::TurbulenceConfig& config() const { return config_; }
    std::vector<TransportEquationState*> transportEquations(
        PhaseEquationState& state);
    std::vector<const TransportEquationState*> transportEquations(
        const PhaseEquationState& state) const;

private:
    FDM::TurbulenceConfig config_;
    std::vector<PhaseEquationState> states_;
    bool active_ = false;

    PhaseEquationState* find(size_t phase);
    const PhaseEquationState* find(size_t phase) const;
};

} // namespace SF::Turbulence
