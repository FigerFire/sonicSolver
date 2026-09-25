/// @file SF_homogeneousMultiphaseStateModel.h
/// @brief 守恒变量 FluidStateModel 与工厂实现。

#pragma once

#include "SF_fluidStateModel.h"
#include "SF_stiffenedGasEOS.h"
#include "SF_transferLedger.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::Physics::FluidStateModel {

struct Component {
    std::string phase;
    std::string species;
    EOS::ModelPtr eos;
    double initialMassFraction = 1.0;
    double dynamicViscosity = 0.0;
    double thermalConductivity = 0.0;
};

/// @brief 单压、单温、单速度的均质多相方程组。
///
/// 每个相—组分只保存 `partialDensity=alpha_phase*rho_phase*Y_species`。
/// 相质量和总密度分别由索引求和得到，不另存冗余 rho/phaseMass。
class HomogeneousMultiphaseStateModel final : public Model {
public:
    explicit HomogeneousMultiphaseStateModel(std::vector<Component> components)
        : components_(std::move(components)) {
        if (components_.empty())
            throw std::runtime_error("HomogeneousMultiphaseStateModel needs components.");
        for (size_t n = 0; n < components_.size(); ++n) {
            const Component& c = components_[n];
            if (!c.eos) throw std::runtime_error("Homogeneous component has no EOS.");
            if (!dynamic_cast<const EOS::StiffenedGasEOS*>(c.eos.get()))
                throw std::runtime_error(
                    "Homogeneous baseline currently requires explicit stiffenedGas EOS for every component.");
            if (!std::isfinite(c.initialMassFraction) || c.initialMassFraction < 0.0) {
                throw std::runtime_error("Homogeneous component initial mass fraction is invalid.");
            }
            variables_.push_back({"partialDensity." + c.phase + "." + c.species,
                                  phaseIndex(c.phase), (int)n});
            layoutVariables_.push_back({
                variables_.back().name,
                State::FieldLocation::Cell,
                State::ConservationKind::Conservative,
                State::UpdatePolicy::Explicit,
                {}, "kg/m3"});
        }
        variables_.push_back({"rhoU"});
        variables_.push_back({"rhoV"});
        variables_.push_back({"rhoW"});
        variables_.push_back({"rhoE"});
        layoutVariables_.push_back({"rhoU", State::FieldLocation::Cell,
            State::ConservationKind::Conservative,
            State::UpdatePolicy::Explicit, {}, "kg/(m2 s)"});
        layoutVariables_.push_back({"rhoV", State::FieldLocation::Cell,
            State::ConservationKind::Conservative,
            State::UpdatePolicy::Explicit, {}, "kg/(m2 s)"});
        layoutVariables_.push_back({"rhoW", State::FieldLocation::Cell,
            State::ConservationKind::Conservative,
            State::UpdatePolicy::Explicit, {}, "kg/(m2 s)"});
        layoutVariables_.push_back({"rhoE", State::FieldLocation::Cell,
            State::ConservationKind::Conservative,
            State::UpdatePolicy::Explicit, {}, "J/m3"});
        layout_ = State::StateLayout(std::move(layoutVariables_));
    }

    Family family() const override { return Family::HomogeneousMultiphase; }
    const State::StateLayout& layout() const override { return layout_; }
    const std::vector<Variable>& primaryVariables() const override { return variables_; }
    int densityVariableCount() const override { return (int)components_.size(); }
    int momentumIndex(int component) const override {
        return densityVariableCount() + component;
    }
    int energyIndex() const override { return densityVariableCount() + 3; }

    /// @brief 将相级 TransferLedger 映射为本方程组的守恒源项。
    /// @param ledger 与 Field 无关的相间交换账本。
    /// @param transferComponentByPhase 每个相接收传质的 partialDensity 索引。
    /// @param sources 输出的 cell-major FluidStateModel 源项。
    void mapTransfersToSources(
            const InterphaseTransfer::TransferLedger& ledger,
            const std::vector<int>& transferComponentByPhase,
            InterphaseTransfer::SourceBundle& sources) const {
        if (ledger.phaseCount() != phaseCount()
            || sources.cellCount() != ledger.cellCount()
            || sources.variableCount() != variableCount()
            || static_cast<int>(transferComponentByPhase.size()) != phaseCount()) {
            throw std::runtime_error(
                "Homogeneous transfer mapping dimensions do not match FluidStateModel.");
        }
        ledger.validate();
        for (int phase = 0; phase < phaseCount(); ++phase) {
            const int variable = transferComponentByPhase[(size_t)phase];
            if (variable < 0 || variable >= densityVariableCount()
                || variables_[(size_t)variable].phase != phase) {
                throw std::runtime_error(
                    "Homogeneous transfer component route does not belong to its phase.");
            }
        }
        for (int cell = 0; cell < ledger.cellCount(); ++cell) {
            for (int phase = 0; phase < phaseCount(); ++phase) {
                sources(cell, transferComponentByPhase[(size_t)phase]) +=
                    ledger.phaseMassSource(cell, phase);
            }
            for (int d = 0; d < 3; ++d) {
                double mixtureMomentum = 0.0;
                for (int phase = 0; phase < phaseCount(); ++phase) {
                    mixtureMomentum += ledger.phaseMomentumSource(cell, phase, d);
                }
                sources(cell, momentumIndex(d)) += mixtureMomentum;
            }
            double mixtureEnergy = ledger.mixtureEnergySource(cell);
            for (int phase = 0; phase < phaseCount(); ++phase) {
                mixtureEnergy += ledger.phaseEnergySource(cell, phase);
            }
            sources(cell, energyIndex()) += mixtureEnergy;
        }
    }

    /// @brief 相数量。
    int phaseCount() const { return static_cast<int>(phases_.size()); }

    /// @brief 按 FluidStateModel 相索引返回相名。
    const std::string& phaseName(int phase) const {
        if (phase < 0 || phase >= phaseCount()) {
            throw std::runtime_error("Homogeneous FluidStateModel phase index is out of range.");
        }
        return phases_[(size_t)phase];
    }

    /// @brief 按相名、组分名查找 partialDensity 主变量索引。
    int componentIndex(const std::string& phase,
                       const std::string& species) const {
        for (int k = 0; k < densityVariableCount(); ++k) {
            if (components_[(size_t)k].phase == phase
                && components_[(size_t)k].species == species) return k;
        }
        return -1;
    }

    /// @brief 返回相中的第一个 partialDensity 变量索引。
    int firstComponentIndex(int phase) const { return firstComponent(phase); }

    /// @brief 返回相的唯一组分索引；多组分相必须由调用者显式指定 species。
    int uniqueComponentIndex(const std::string& phase) const {
        int result = -1;
        for (int k = 0; k < densityVariableCount(); ++k) {
            if (components_[(size_t)k].phase != phase) continue;
            if (result >= 0) {
                throw std::runtime_error(
                    "Phase '" + phase
                    + "' has multiple species; phase change must name the transfer species explicitly.");
            }
            result = k;
        }
        if (result < 0) throw std::runtime_error("Unknown phase '" + phase + "'.");
        return result;
    }

    /// @brief 返回相名对应索引。
    int findPhaseIndex(const std::string& phase) const {
        const auto it = std::find(phases_.begin(), phases_.end(), phase);
        return it == phases_.end() ? -1 : static_cast<int>(std::distance(phases_.begin(), it));
    }

    ThermodynamicState close(const double* q, int count) const override {
        if (!q || count != variableCount())
            throw std::runtime_error("Homogeneous equation-set state size mismatch.");
        ThermodynamicState result;
        const int n = densityVariableCount();
        for (int k = 0; k < n; ++k) {
            if (!std::isfinite(q[k]) || q[k] < 0.0)
                throw std::runtime_error("Homogeneous state has invalid partial density.");
            result.density += q[k];
        }
        if (!std::isfinite(result.density) || result.density <= 0.0)
            throw std::runtime_error("Homogeneous state has zero total density.");
        for (int d = 0; d < 3; ++d)
            result.velocity[(size_t)d] = q[n + d] / result.density;
        const double kinetic = 0.5 * result.density
            * (result.velocity[0]*result.velocity[0]
               + result.velocity[1]*result.velocity[1]
               + result.velocity[2]*result.velocity[2]);
        result.internalEnergyDensity = q[energyIndex()] - kinetic;
        if (!std::isfinite(result.internalEnergyDensity))
            throw std::runtime_error("Homogeneous state has invalid internal energy.");
        closeStiffenedMixture(q, result);
        return result;
    }

    double totalEnergyFromPressure(const double* q, int count,
                                   double pressure) const override {
        requireConservedState(q, count);
        double temperature = 0.0;
        const double internal = internalEnergyAtPressure(
            q, pressure, temperature);
        return internal + kineticEnergy(q);
    }

    double totalEnergyFromTemperature(const double* q, int count,
                                      double temperature) const override {
        requireConservedState(q, count);
        if (!std::isfinite(temperature) || temperature <= 0.0) {
            throw std::runtime_error(
                "Homogeneous temperature boundary requires finite T>0.");
        }
        const double pressure = pressureAtTemperature(q, temperature);
        double verifiedTemperature = 0.0;
        const double internal = internalEnergyAtPressure(
            q, pressure, verifiedTemperature);
        if (std::abs(verifiedTemperature-temperature)
                > 1.0e-10*std::max(std::abs(temperature), 1.0)) {
            throw std::runtime_error(
                "Homogeneous temperature inversion violated the volume constraint.");
        }
        return internal + kineticEnergy(q);
    }

    /// @brief 由共压、共温和相体积分数构造无冗余守恒状态。
    /// @details 输入不满足体积分数和为一或 EOS 体积约束时直接报错，绝不归一化。
    std::vector<double> conservativeFromPressureTemperature(
            const std::vector<double>& volumeFraction,
            double pressure, double temperature,
            const std::array<double, 3>& velocity) const {
        if ((int)volumeFraction.size() != (int)phases_.size()) {
            throw std::runtime_error("Homogeneous initial phase volumeFraction size mismatch.");
        }
        std::vector<double> q((size_t)variableCount(), 0.0);
        double alphaSum = 0.0;
        double totalDensity = 0.0;
        double internalEnergyDensity = 0.0;
        for (int phase = 0; phase < (int)phases_.size(); ++phase) {
            const double alpha = volumeFraction[(size_t)phase];
            if (!std::isfinite(alpha) || alpha < 0.0) {
                throw std::runtime_error("Homogeneous initial volume fraction is invalid.");
            }
            int first = -1;
            double ySum = 0.0;
            for (int k = 0; k < densityVariableCount(); ++k) {
                if (variables_[(size_t)k].phase == phase) {
                    if (first < 0) first = k;
                    ySum += components_[(size_t)k].initialMassFraction;
                }
            }
            if (first < 0 || std::abs(ySum - 1.0) > 1.0e-12) {
                throw std::runtime_error("Homogeneous initial phase species fractions must sum exactly to one.");
            }
            const EOS::State state = components_[(size_t)first].eos
                ->fromPressureTemperature(pressure, temperature);
            const double phaseMass = alpha * state.density;
            for (int k = 0; k < densityVariableCount(); ++k) {
                if (variables_[(size_t)k].phase == phase) {
                    q[(size_t)k] = phaseMass * components_[(size_t)k].initialMassFraction;
                }
            }
            alphaSum += alpha;
            totalDensity += phaseMass;
            internalEnergyDensity += phaseMass * state.internalEnergy;
        }
        if (std::abs(alphaSum - 1.0) > 1.0e-12 || totalDensity <= 0.0) {
            throw std::runtime_error("Homogeneous initial phase volume fractions must sum exactly to one.");
        }
        for (int d = 0; d < 3; ++d) {
            q[(size_t)momentumIndex(d)] = totalDensity * velocity[(size_t)d];
        }
        q[(size_t)energyIndex()] = internalEnergyDensity
            + 0.5 * totalDensity * (velocity[0]*velocity[0]
                                  + velocity[1]*velocity[1]
                                  + velocity[2]*velocity[2]);
        const ThermodynamicState verified = close(q.data(), variableCount());
        // 刚化压力使总内能达到 O(1e9)，压力反解的舍入误差会被放大；
        // 这里只放宽一致性判据，不修改或投影守恒状态。
        if (std::abs(verified.pressure - pressure) > 1.0e-9 * std::max(std::abs(pressure), 1.0)
            || std::abs(verified.temperature - temperature) > 1.0e-9 * std::max(std::abs(temperature), 1.0)) {
            throw std::runtime_error(
                "Homogeneous initial state is not EOS volume-consistent: requested p="
                + std::to_string(pressure) + ", T=" + std::to_string(temperature)
                + ", closed p=" + std::to_string(verified.pressure)
                + ", T=" + std::to_string(verified.temperature) + ".");
        }
        return q;
    }

private:
    std::vector<Component> components_;
    std::vector<Variable> variables_;
    std::vector<State::VariableDescriptor> layoutVariables_;
    State::StateLayout layout_;
    std::vector<std::string> phases_;

    int phaseIndex(const std::string& phase) {
        auto it = std::find(phases_.begin(), phases_.end(), phase);
        if (it != phases_.end()) return (int)std::distance(phases_.begin(), it);
        phases_.push_back(phase);
        return (int)phases_.size() - 1;
    }

    void requireConservedState(const double* q, int count) const {
        if (!q || count != variableCount()) {
            throw std::runtime_error(
                "Homogeneous boundary state size differs from FluidStateModel.");
        }
        double rho = 0.0;
        for (int k = 0; k < densityVariableCount(); ++k) {
            if (!std::isfinite(q[k]) || q[k] < 0.0) {
                throw std::runtime_error(
                    "Homogeneous boundary has invalid partial density.");
            }
            rho += q[k];
        }
        if (!std::isfinite(rho) || rho <= 0.0) {
            throw std::runtime_error(
                "Homogeneous boundary has zero total density.");
        }
    }

    double totalDensity(const double* q) const {
        double rho = 0.0;
        for (int k = 0; k < densityVariableCount(); ++k) rho += q[k];
        return rho;
    }

    double kineticEnergy(const double* q) const {
        const double rho = totalDensity(q);
        double momentum2 = 0.0;
        for (int d = 0; d < 3; ++d) {
            const double momentum = q[momentumIndex(d)];
            if (!std::isfinite(momentum)) {
                throw std::runtime_error(
                    "Homogeneous boundary has non-finite momentum.");
            }
            momentum2 += momentum*momentum;
        }
        return 0.5*momentum2/rho;
    }

    double temperatureAtPressure(const double* partialDensity,
                                 double pressure) const {
        double inverseT = 0.0;
        for (int phase = 0; phase < static_cast<int>(phases_.size()); ++phase) {
            const int k = firstComponent(phase);
            const double mass = phaseMass(partialDensity, phase);
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(
                *components_[(size_t)k].eos);
            const double denominator = pressure + eos.pInfinity();
            if (!std::isfinite(denominator) || denominator <= 0.0) {
                throw std::runtime_error(
                    "Homogeneous pressure boundary crossed -pInfinity.");
            }
            inverseT += mass*(eos.gamma()-1.0)*eos.cv()/denominator;
        }
        if (!std::isfinite(inverseT) || inverseT <= 0.0) {
            throw std::runtime_error(
                "Homogeneous pressure boundary has singular volume constraint.");
        }
        return 1.0/inverseT;
    }

    double internalEnergyAtPressure(const double* partialDensity,
                                    double pressure,
                                    double& temperature) const {
        temperature = temperatureAtPressure(partialDensity, pressure);
        double energy = 0.0;
        for (int phase = 0; phase < static_cast<int>(phases_.size()); ++phase) {
            const int k = firstComponent(phase);
            const double mass = phaseMass(partialDensity, phase);
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(
                *components_[(size_t)k].eos);
            const double rho = (pressure+eos.pInfinity())
                / ((eos.gamma()-1.0)*eos.cv()*temperature);
            const double alpha = mass/rho;
            energy += mass*(eos.referenceEnergy()+eos.cv()*temperature)
                    + alpha*eos.pInfinity();
        }
        if (!std::isfinite(energy)) {
            throw std::runtime_error(
                "Homogeneous pressure boundary produced non-finite energy.");
        }
        return energy;
    }

    double pressureAtTemperature(const double* partialDensity,
                                 double temperature) const {
        double lower = -std::numeric_limits<double>::max();
        for (const Component& component : components_) {
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(
                *component.eos);
            lower = std::max(lower, -eos.pInfinity());
        }
        lower += std::max(1.0, std::abs(lower)*1.0e-12);
        auto volumeResidual = [&](double pressure) {
            double alpha = 0.0;
            for (int phase = 0; phase < static_cast<int>(phases_.size()); ++phase) {
                const int k = firstComponent(phase);
                const double mass = phaseMass(partialDensity, phase);
                const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(
                    *components_[(size_t)k].eos);
                const double rho = (pressure+eos.pInfinity())
                    / ((eos.gamma()-1.0)*eos.cv()*temperature);
                alpha += mass/rho;
            }
            return alpha-1.0;
        };
        double upper = std::max(1.0e5, lower+1.0e5);
        double fLower = volumeResidual(lower);
        double fUpper = volumeResidual(upper);
        for (int expand = 0; fLower*fUpper > 0.0 && expand < 80; ++expand) {
            upper = 2.0*upper+1.0;
            fUpper = volumeResidual(upper);
        }
        if (!std::isfinite(fLower) || !std::isfinite(fUpper)
            || fLower*fUpper > 0.0) {
            throw std::runtime_error(
                "Homogeneous temperature boundary could not bracket pressure.");
        }
        double pressure = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            pressure = 0.5*(lower+upper);
            const double f = volumeResidual(pressure);
            if (!std::isfinite(f)) {
                throw std::runtime_error(
                    "Homogeneous temperature boundary produced non-finite residual.");
            }
            if (std::abs(f) <= 1.0e-12) return pressure;
            if (fLower*f <= 0.0) {
                upper = pressure;
                fUpper = f;
            } else {
                lower = pressure;
                fLower = f;
            }
        }
        throw std::runtime_error(
            "Homogeneous temperature boundary pressure solve did not converge.");
    }

    double residualAtPressure(const double* partialDensity,
                              double internalEnergyDensity,
                              double pressure,
                              double& temperature) const {
        double inverseT = 0.0;
        for (int phase = 0; phase < (int)phases_.size(); ++phase) {
            const int k = firstComponent(phase);
            const double mass = phaseMass(partialDensity, phase);
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(*components_[(size_t)k].eos);
            const double denominator = pressure + eos.pInfinity();
            if (!std::isfinite(denominator) || denominator <= 0.0)
                throw std::runtime_error("Mixture EOS pressure crossed -pInfinity.");
            inverseT += mass * (eos.gamma() - 1.0) * eos.cv()
                      / denominator;
        }
        if (!std::isfinite(inverseT) || inverseT <= 0.0)
            throw std::runtime_error("Mixture EOS volume constraint is singular.");
        temperature = 1.0 / inverseT;
        double energy = 0.0;
        for (int phase = 0; phase < (int)phases_.size(); ++phase) {
            const int k = firstComponent(phase);
            const double mass = phaseMass(partialDensity, phase);
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(*components_[(size_t)k].eos);
            const double rho = (pressure + eos.pInfinity())
                / ((eos.gamma() - 1.0) * eos.cv() * temperature);
            const double alpha = mass / rho;
            energy += mass
                    * (eos.referenceEnergy() + eos.cv() * temperature)
                    + alpha * eos.pInfinity();
        }
        return energy - internalEnergyDensity;
    }

    void closeStiffenedMixture(const double* partialDensity,
                               ThermodynamicState& result) const {
        double lower = -std::numeric_limits<double>::max();
        for (const Component& c : components_) {
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(*c.eos);
            lower = std::max(lower, -eos.pInfinity());
        }
        lower += std::max(1.0, std::abs(lower) * 1.0e-12);
        double upper = std::max(1.0e5, lower + 1.0e5);
        double tLower = 0.0, tUpper = 0.0;
        double fLower = residualAtPressure(partialDensity,
                                           result.internalEnergyDensity,
                                           lower, tLower);
        double fUpper = residualAtPressure(partialDensity,
                                           result.internalEnergyDensity,
                                           upper, tUpper);
        for (int expand = 0; fLower * fUpper > 0.0 && expand < 80; ++expand) {
            upper = 2.0 * upper + 1.0;
            fUpper = residualAtPressure(partialDensity,
                                        result.internalEnergyDensity,
                                        upper, tUpper);
        }
        if (!std::isfinite(fLower) || !std::isfinite(fUpper)
            || fLower * fUpper > 0.0)
            throw std::runtime_error("Mixture stiffened-gas closure could not bracket pressure.");

        double pressure = 0.0, temperature = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            pressure = 0.5 * (lower + upper);
            const double f = residualAtPressure(partialDensity,
                                                 result.internalEnergyDensity,
                                                 pressure, temperature);
            if (!std::isfinite(f))
                throw std::runtime_error("Mixture stiffened-gas closure produced non-finite residual.");
            if (std::abs(f) <= 1.0e-12
                    * std::max(std::abs(result.internalEnergyDensity), 1.0))
                break;
            if (fLower * f <= 0.0) {
                upper = pressure;
                fUpper = f;
            } else {
                lower = pressure;
                fLower = f;
            }
            if (iteration == 99)
                throw std::runtime_error("Mixture stiffened-gas pressure solve did not converge.");
        }

        result.pressure = pressure;
        result.temperature = temperature;
        result.phaseMass.assign(phases_.size(), 0.0);
        result.volumeFraction.assign(phases_.size(), 0.0);
        double woodCompressibility = 0.0;
        result.dynamicViscosity = 0.0;
        result.thermalConductivity = 0.0;
        for (int phase = 0; phase < (int)phases_.size(); ++phase) {
            const int k = firstComponent(phase);
            const double mass = phaseMass(partialDensity, phase);
            const auto& eos = static_cast<const EOS::StiffenedGasEOS&>(*components_[(size_t)k].eos);
            const double rho = (pressure + eos.pInfinity())
                / ((eos.gamma() - 1.0) * eos.cv() * temperature);
            const double alpha = mass / rho;
            result.phaseMass[(size_t)phase] = mass;
            result.volumeFraction[(size_t)phase] += alpha;
            const Component& component = components_[(size_t)k];
            if (!std::isfinite(component.dynamicViscosity)
                || component.dynamicViscosity < 0.0
                || !std::isfinite(component.thermalConductivity)
                || component.thermalConductivity < 0.0) {
                throw std::runtime_error(
                    "Homogeneous phase transport properties are invalid.");
            }
            result.dynamicViscosity += alpha*component.dynamicViscosity;
            result.thermalConductivity += alpha*component.thermalConductivity;
            const double c2 = eos.gamma() * (pressure + eos.pInfinity()) / rho;
            woodCompressibility += alpha / (rho * c2);
        }
        if (!std::isfinite(woodCompressibility) || woodCompressibility <= 0.0)
            throw std::runtime_error("Mixture EOS produced invalid Wood compressibility.");
        result.soundSpeed = std::sqrt(1.0 / (result.density * woodCompressibility));
        double alphaSum = 0.0;
        for (double a : result.volumeFraction) alphaSum += a;
        if (!std::isfinite(alphaSum) || std::abs(alphaSum - 1.0) > 1.0e-10)
            throw std::runtime_error("Mixture EOS violates the volume constraint.");
    }

    int firstComponent(int phase) const {
        for (int k = 0; k < densityVariableCount(); ++k) {
            if (variables_[(size_t)k].phase == phase) return k;
        }
        throw std::runtime_error("Homogeneous FluidStateModel has an empty phase.");
    }

    double phaseMass(const double* partialDensity, int phase) const {
        double mass = 0.0;
        for (int k = 0; k < densityVariableCount(); ++k) {
            if (variables_[(size_t)k].phase == phase) mass += partialDensity[k];
        }
        if (!std::isfinite(mass) || mass < 0.0) {
            throw std::runtime_error("Homogeneous phase mass is invalid.");
        }
        return mass;
    }
};

} // namespace SF::Physics::FluidStateModel
