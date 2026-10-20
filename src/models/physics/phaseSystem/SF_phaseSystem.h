#pragma once

/// @file SF_phaseSystem.h
/// @brief 单速度与多速度相系统共享的抽象边界。

#include "SF_phaseProperties.h"
#include "SF_phaseState.h"
#include "SF_transferLedger.h"

#include <memory>
#include <vector>

namespace SF::Physics::PhaseSystems {

/// @brief 一个相对的半隐式动量耦合系数。
struct PhasePairMomentumCoupling {
    size_t first = 0;
    size_t second = 0;
    ScalarField coefficient;     ///< drag + virtual-mass 的线性系数。
    ScalarField dragCoefficient; ///< 湍流扩散速度闭式使用的纯 drag 系数。
};

/// @brief 相方程贡献容器。
///
/// 相间模型先写入并校验成对守恒项；随后外部 source provider 可继续累加重力、
/// MRF 和 wall heat 等非零总源。求解器只消费这个统一容器。
struct PhaseEquationSources {
    std::vector<ScalarField> mass;
    std::vector<std::array<ScalarField, 3>> momentum;
    std::vector<ScalarField> energy;
    ScalarField wallBoilingMass;
    ScalarField wallBoilingConvectiveHeat;
    ScalarField wallBoilingQuenchingHeat;
    ScalarField wallBoilingEvaporativeHeat;
    ScalarField wallBoilingDepartureDiameter;
    ScalarField interphaseMechanicalHeating;
    InterphaseTransfer::TransferLedger transferLedger;
    std::vector<PhasePairMomentumCoupling> momentumCouplings;
    void setupLike(
        const Field& field,
        const std::vector<PhaseState>& phases,
        const std::vector<std::array<size_t, 2>>& phasePairs);
    void clear();
    /// @brief 外部 provider 每次装配前只清空其诊断 ledger。
    void clearExternalLedger();
};

/// @brief PhaseSystem 只拥有相状态和物理闭式，不依赖压力基求解器。
class PhaseSystem {
public:
    virtual ~PhaseSystem() = default;
    virtual void initialize(const Field& geometry) = 0;
    /// @brief 由唯一主状态和共享压力刷新所有 primitive cache。
    virtual void recoverPrimitiveState() = 0;
    /// @brief 用 N-1 个独立 phaseMass 闭合参考相主状态。
    virtual void reconstructReferencePhaseMass() = 0;
    /// @brief 计算物理相间闭式；previousVelocity 由 solver workspace 提供。
    virtual void computeInterphase(
        double dt,
        const std::vector<PhaseVectorField>& previousVelocity) = 0;
    virtual void validateState(const std::string& stage) const = 0;

    virtual const Field& geometry() const = 0;
    virtual ScalarField& sharedPressure() = 0;
    virtual const ScalarField& sharedPressure() const = 0;
    virtual std::vector<PhaseState>& phases() = 0;
    virtual const std::vector<PhaseState>& phases() const = 0;
    virtual PhaseEquationSources& sources() = 0;
    virtual const PhaseEquationSources& sources() const = 0;
    virtual size_t referencePhaseIndex() const = 0;
    virtual const Multiphase::PhaseProperties& phaseProperties(
        size_t phase) const = 0;
    virtual const Multiphase::MultiPhaseConfig& config() const = 0;
};

/// @brief 构造由 phaseProperties 声明的独立相系统。
std::unique_ptr<PhaseSystem> makePhaseSystem(
    const Multiphase::MultiPhaseConfig& config);

} // namespace SF::Physics::PhaseSystems
