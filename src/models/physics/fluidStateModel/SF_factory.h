/// @file SF_factory.h
/// @brief 守恒变量 FluidStateModel 与工厂实现。

#pragma once

#include "SF_homogeneousMultiphaseStateModel.h"
#include "SF_singleFluidStateModel.h"
#include "SF_perfectGasEOS.h"
#include "SF_phaseProperties.h"
#include "SF_field.h"

#include <memory>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::FluidStateModel {

/// @brief 构造默认单流体 perfect-gas FluidStateModel。
inline std::shared_ptr<const SingleFluidStateModel>
makeSingleFluidPerfectGas(
        double gamma,double gasConstant,
        double dynamicViscosity,double prandtl) {
    if (!std::isfinite(prandtl)||prandtl<=0.0) {
        throw std::runtime_error(
            "Single-fluid FluidStateModel requires positive finite Prandtl.");
    }
    const double cp=gamma*gasConstant/(gamma-1.0);
    const double conductivity=dynamicViscosity*cp/prandtl;
    return std::make_shared<const SingleFluidStateModel>(
        std::make_shared<EOS::PerfectGasEOS>(gamma,gasConstant),
        dynamicViscosity,conductivity);
}

/// @brief 用 rho/U 与显式 p 或 T 初值冷启动单流体守恒能量。
inline void initializeSingleFluidField(
        Field& field,
        const std::shared_ptr<const SingleFluidStateModel>& equations,
        const std::vector<BCSetting<double>>& pressureInitial,
        const std::vector<BCSetting<double>>& temperatureInitial) {
    if (!equations||field.NVar()!=equations->variableCount()) {
        throw std::runtime_error(
            "Single-fluid cold start requires a five-variable Field.");
    }
    std::vector<double> pressure((size_t)field.TotalSize(),0.0);
    std::vector<double> temperature((size_t)field.TotalSize(),0.0);
    std::vector<unsigned char> hasPressure(
        (size_t)field.TotalSize(),0);
    std::vector<unsigned char> hasTemperature(
        (size_t)field.TotalSize(),0);
    auto collect=[&](
            const std::vector<BCSetting<double>>& settings,
            std::vector<double>& values,
            std::vector<unsigned char>& assigned) {
        for (const auto& setting : settings) {
            const auto& cells=field.getSet(setting.name);
            for (int cell : cells) {
                if (cell<0||cell>=field.TotalSize()) {
                    throw std::runtime_error(
                        "Single-fluid initial-condition set '"+setting.name
                        +"' contains an out-of-range Field index "
                        +std::to_string(cell)+".");
                }
                values[(size_t)cell]=setting.value;
                assigned[(size_t)cell]=1;
            }
        }
    };
    collect(pressureInitial,pressure,hasPressure);
    collect(temperatureInitial,temperature,hasTemperature);
    const int ng=field.NG();
    for (int k=ng;k<ng+field.NZ();++k)
        for (int j=ng;j<ng+field.NY();++j)
            for (int i=ng;i<ng+field.NX();++i) {
                const int cell=field.getIdx(i,j,k);
                double q[5];
                for (int v=0;v<5;++v) q[v]=field(i,j,k,v);
                const double T=temperature[(size_t)cell];
                const double p=pressure[(size_t)cell];
                if (hasTemperature[(size_t)cell]) {
                    if (!std::isfinite(T)||T<=0.0) {
                        throw std::runtime_error(
                            "Single-fluid cold start received a non-positive "
                            "temperature in a physical cell.");
                    }
                    field(i,j,k,E)=
                        equations->totalEnergyFromTemperature(q,5,T);
                } else if (hasPressure[(size_t)cell]) {
                    if (!std::isfinite(p)||p<=0.0) {
                        throw std::runtime_error(
                            "Single-fluid cold start received a non-positive "
                            "pressure in a physical cell.");
                    }
                    field(i,j,k,E)=
                        equations->totalEnergyFromPressure(q,5,p);
                } else {
                    throw std::runtime_error(
                        "Single-fluid cold start leaves a physical cell "
                        "without positive explicit p or T.");
                }
            }
    field.setStateModel(equations);
    for (int k=ng;k<ng+field.NZ();++k)
        for (int j=ng;j<ng+field.NY();++j)
            for (int i=ng;i<ng+field.NX();++i) {
                const int cell=field.getIdx(i,j,k);
                const auto state=field.thermodynamicState(i,j,k);
                const double T=temperature[(size_t)cell];
                const double p=pressure[(size_t)cell];
                const bool temperatureSpecified=
                    hasTemperature[(size_t)cell]!=0;
                const double expected=temperatureSpecified ? T : p;
                const double actual=temperatureSpecified
                    ? state.temperature : state.pressure;
                if (!std::isfinite(actual)
                    ||std::abs(actual-expected)
                        >1.0e-10*std::max(std::abs(expected),1.0)) {
                    throw std::runtime_error(
                        "Single-fluid EOS cold-start inversion mismatch at ("
                        +std::to_string(i)+","+std::to_string(j)+","
                        +std::to_string(k)+").");
                }
            }
}

/// @brief 绑定 restart 守恒 Q，并逐点完成 EOS 派生状态闭合检查。
inline void initializeSingleFluidRestart(
        Field& field,
        const std::shared_ptr<const SingleFluidStateModel>& equations) {
    if (!equations||field.NVar()!=equations->variableCount()) {
        throw std::runtime_error(
            "Single-fluid restart Q layout does not match FluidStateModel.");
    }
    field.setStateModel(equations);
    const int ng=field.NG();
    for (int k=ng;k<ng+field.NZ();++k)
        for (int j=ng;j<ng+field.NY();++j)
            for (int i=ng;i<ng+field.NX();++i)
                (void)field.thermodynamicState(i,j,k);
}

inline EOS::ModelPtr makeEOS(const Multiphase::PhaseProperties& phase) {
    const std::string model = Multiphase::normalizeModelType(phase.thermoModel);
    if (model == "stiffenedgas") {
        return std::make_shared<EOS::StiffenedGasEOS>(
            phase.gamma, phase.Cv, phase.pInfinity, phase.e0);
    }
    if (model == "perfectgas" || model == "idealgas") {
        return std::make_shared<EOS::PerfectGasEOS>(
            phase.gamma, phase.gasConstant, phase.e0);
    }
    throw std::runtime_error(
        "FluidStateModel EOS factory does not map legacy thermo model '"
        + phase.thermoModel + "'. Configure an explicit EOS model.");
}

inline std::shared_ptr<const HomogeneousMultiphaseStateModel>
makeHomogeneous(const Multiphase::MultiPhaseConfig& config) {
    if (config.phaseChange.enabled) {
        const Multiphase::PhaseProperties* liquid = nullptr;
        const Multiphase::PhaseProperties* vapor = nullptr;
        for (const auto& phase : config.phases) {
            if (phase.role == Multiphase::PhaseRole::Liquid) liquid = &phase;
            if (phase.role == Multiphase::PhaseRole::Gas) vapor = &phase;
        }
        if (!liquid || !vapor)
            throw std::runtime_error(
                "Phase-change FluidStateModel requires explicit liquid/gas roles.");
        const double encodedLatentHeat = vapor->e0 - liquid->e0;
        const double scale = std::max(std::abs(config.phaseChange.latentHeat), 1.0);
        if (!std::isfinite(encodedLatentHeat)
            || std::abs(encodedLatentHeat-config.phaseChange.latentHeat)
                > 1.0e-10*scale)
            throw std::runtime_error(
                "EOS referenceEnergy difference must equal phaseChange latentHeat; "
                "do not add a separate -mdot*L source.");
    }
    std::vector<Component> components;
    components.reserve(config.phases.size());
    for (const auto& phase : config.phases) {
        if (phase.species.empty()) {
            components.push_back({phase.name, "bulk", makeEOS(phase), 1.0,
                                  phase.viscosity,
                                  phase.thermalConductivity});
        } else {
            for (const auto& species : phase.species) {
                components.push_back({phase.name, species.name, makeEOS(phase),
                                      species.massFraction,
                                      phase.viscosity,
                                      phase.thermalConductivity});
            }
        }
    }
    return std::make_shared<const HomogeneousMultiphaseStateModel>(
        std::move(components));
}

/// @brief 用 case 的共压、共温和相体积分数初始化一个已完成网格装配的 Field。
/// @details 先保存旧五变量中的速度，再显式重分配为 FluidStateModel 布局；不复制或
/// 投影旧 rho/rhoE，避免把固定-gamma 状态伪装成多相热力学状态。
inline void initializeHomogeneousField(
        Field& field,
        const Multiphase::MultiPhaseConfig& config,
        const std::shared_ptr<const HomogeneousMultiphaseStateModel>& equations) {
    if (!equations || field.NVar() != 5) {
        throw std::runtime_error("Homogeneous initialization requires an unconverted five-variable Field.");
    }
    const int cells = field.TotalSize();
    std::vector<std::array<double, 3>> velocity((size_t)cells);
    for (int k = 0; k < field.MZ(); ++k) for (int j = 0; j < field.MY(); ++j)
        for (int i = 0; i < field.MX(); ++i) {
            const int id = field.getIdx(i, j, k);
            const double rho = field(i, j, k, RHO);
            if (!std::isfinite(rho) || rho <= 0.0) {
                throw std::runtime_error("Homogeneous initialization found invalid legacy density.");
            }
            velocity[(size_t)id] = {field(i,j,k,RU)/rho, field(i,j,k,RV)/rho, field(i,j,k,RW)/rho};
        }
    std::vector<double> alpha;
    alpha.reserve(config.phases.size());
    for (const auto& phase : config.phases) alpha.push_back(phase.volumeFraction);
    field.resizeConservedVariables(equations->variableCount());
    field.setStateModel(equations);
    for (int k = 0; k < field.MZ(); ++k) for (int j = 0; j < field.MY(); ++j)
        for (int i = 0; i < field.MX(); ++i) {
            const std::vector<double> q = equations->conservativeFromPressureTemperature(
                alpha, config.initialPressure, config.temperature.defaultValue,
                velocity[(size_t)field.getIdx(i,j,k)]);
            for (int v = 0; v < field.NVar(); ++v) field(i,j,k,v) = q[(size_t)v];
        }
}

} // namespace SF::Physics::FluidStateModel
