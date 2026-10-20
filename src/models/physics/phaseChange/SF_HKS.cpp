/// @file SF_HKS.cpp
/// @brief 相变模型、饱和性质与守恒传递实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_phaseChange.h"
#include "SF_utility.h"

#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {

void addHKSRates(const ModelContext& ctx, std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double a = coefficient(pc, "accommodationCoefficient", 1.0);
    const double molarMass = coefficient(pc, "molarMass", 0.01801528);
    const double gasConstant = coefficient(pc, "gasConstant", 8.314462618);
    const double areaDensity = coefficient(pc, "interfaceAreaDensity", 1.0);
    if (a <= 0.0 || a >= 2.0 || molarMass <= 0.0
        || gasConstant <= 0.0 || areaDensity <= 0.0) {
        throw std::runtime_error(
            "HKS phaseChange requires 0<a<2, molarMass>0, gasConstant>0, "
            "and interfaceAreaDensity>0.");
    }
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double prefactor =
        (2.0 * a / (2.0 - a)) * std::sqrt(molarMass / (2.0 * pi * gasConstant));
    Math::forInterior(ctx.field, [&](int i, int j, int k) {
        if (ctx.field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double T = ctx.temperature(i, j, k);
        if (!std::isfinite(T) || T <= 0.0) {
            throw std::runtime_error("HKS phaseChange: invalid temperature.");
        }
        const double p =
            Numerics::requirePhysicalState("HKS phaseChange", ctx.field, i, j, k);
        const double psat = pc.saturationPressure > 0.0
            ? pc.saturationPressure
            : coefficient(pc, "pSatCoefficient", 1.0) * pc.saturationTemperature;
        const double surfaceMdot =
            prefactor * (psat / std::sqrt(T) - p / std::sqrt(T));
        mdot[(size_t)ctx.alpha.getIdx(i, j, k)] += areaDensity * surfaceMdot;
    });
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
