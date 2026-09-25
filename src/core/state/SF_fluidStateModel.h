#pragma once

/// @file SF_fluidStateModel.h
/// @brief 可变守恒状态布局与热力学闭合的核心契约。

#include "SF_stateLayout.h"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::Physics::FluidStateModel {

/// @brief 已实现与预留的方程组族；枚举存在不代表已有可运行实现。
enum class Family {
    SingleFluid,
    HomogeneousMultiphase,
    LowMachSingleFluid,
    OneFluidCLSVOF,
    CompressiblePhaseMass,
    FiveEquation,
    EulerianEulerian
};

struct Variable {
    std::string name;
    int phase = -1;
    int species = -1;
};

struct ThermodynamicState {
    double density = 0.0;
    std::array<double, 3> velocity{0.0, 0.0, 0.0};
    double internalEnergyDensity = 0.0;
    double pressure = 0.0;
    double temperature = 0.0;
    double soundSpeed = 0.0;
    double dynamicViscosity = 0.0;
    double thermalConductivity = 0.0;
    std::vector<double> phaseMass;
    std::vector<double> volumeFraction;
};

/// @brief Equation 层向 Field 和算法层提供的方程组接口。
class Model {
public:
    virtual ~Model() = default;
    virtual Family family() const = 0;
    virtual const State::StateLayout& layout() const = 0;
    virtual const std::vector<Variable>& primaryVariables() const = 0;
    virtual int densityVariableCount() const = 0;
    virtual int momentumIndex(int component) const = 0;
    virtual int energyIndex() const = 0;
    virtual ThermodynamicState close(const double* conserved,
                                     int variableCount) const = 0;
    virtual double totalEnergyFromPressure(const double*, int, double) const {
        throw std::runtime_error(
            "This FluidStateModel/EOS does not implement fixed-density pressure inversion.");
    }
    virtual double totalEnergyFromTemperature(const double*, int, double) const {
        throw std::runtime_error(
            "This FluidStateModel/EOS does not implement fixed-density temperature inversion.");
    }
    /// @brief 返回现有特征高阶通量可消费的 PerfectGas gamma；其他 EOS 显式拒绝。
    virtual std::optional<double> perfectGasGamma() const { return std::nullopt; }
    int variableCount() const {
        return static_cast<int>(primaryVariables().size());
    }
};

} // namespace SF::Physics::FluidStateModel
