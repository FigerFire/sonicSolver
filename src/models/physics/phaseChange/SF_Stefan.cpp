/// @file SF_Stefan.cpp
/// @brief 相变模型、饱和性质与守恒传递实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_phaseChange.h"

#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {

void addStefanRates(const ModelContext& ctx, std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double areaDensity = coefficient(pc, "interfaceAreaDensity", 0.0);
    const double kl = coefficient(pc, "liquidThermalConductivity", 0.0);
    const double kv = coefficient(pc, "vaporThermalConductivity", 0.0);
    const double gradL = coefficient(pc, "liquidTemperatureGradient", 0.0);
    const double gradV = coefficient(pc, "vaporTemperatureGradient", 0.0);
    if (areaDensity <= 0.0 || pc.latentHeat <= 0.0) {
        throw std::runtime_error(
            "Stefan phaseChange requires interfaceAreaDensity > 0 and latentHeat > 0.");
    }
    const double surfaceMdot = (kl * gradL - kv * gradV) / pc.latentHeat;
    if (!std::isfinite(surfaceMdot)) {
        throw std::runtime_error("Stefan phaseChange produced non-finite mdot.");
    }
    Math::forInterior(ctx.field, [&](int i, int j, int k) {
        if (ctx.field.CellFlag(i, j, k) != FLUID_CELL) return;
        mdot[(size_t)ctx.alpha.getIdx(i, j, k)] += areaDensity * surfaceMdot;
    });
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
