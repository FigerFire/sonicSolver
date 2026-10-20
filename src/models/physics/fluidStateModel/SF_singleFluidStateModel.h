/// @file SF_singleFluidStateModel.h
/// @brief 守恒变量 FluidStateModel 与工厂实现。

#pragma once

#include "SF_fluidStateModel.h"
#include "SF_EOS.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace SF::Physics::FluidStateModel {

class SingleFluidStateModel final : public Model {
public:
    explicit SingleFluidStateModel(
        EOS::ModelPtr eos,
        double dynamicViscosity=0.0,
        double thermalConductivity=0.0);

    Family family() const override { return Family::SingleFluid; }
    const State::StateLayout& layout() const override { return layout_; }
    const std::vector<Variable>& primaryVariables() const override { return variables_; }
    int densityVariableCount() const override { return 1; }
    int momentumIndex(int component) const override { return 1 + component; }
    int energyIndex() const override { return 4; }

    ThermodynamicState close(
        const double* q,
        int count) const override;

    double totalEnergyFromPressure(
        const double* q,
        int count,
        double pressure) const override;

    double totalEnergyFromTemperature(
        const double* q,
        int count,
        double temperature) const override;
    std::optional<double> perfectGasGamma() const override;

private:
    EOS::ModelPtr eos_;
    double dynamicViscosity_=0.0;
    double thermalConductivity_=0.0;
    std::vector<Variable> variables_;
    State::StateLayout layout_;

    static double totalEnergyFromThermo(
        const double* q,
        int count,
        const EOS::State& thermo);
    static void requireEnergyInput(const double* q,int count);
};

} // namespace SF::Physics::FluidStateModel
