/// @file SF_enthalpyPorosity.cpp
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

void addEnthalpyPorosityRates(const ModelContext& ctx,
                              std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double Tm = coefficient(pc, "meltingTemperature",
                                  pc.saturationTemperature);
    const double width = coefficient(pc, "mushyTemperatureWidth", 1.0);
    const double rate = coefficient(pc, "relaxationCoefficient",
                                    pc.evaporationCoefficient);
    if (Tm <= 0.0 || width <= 0.0 || rate <= 0.0) {
        throw std::runtime_error(
            "EnthalpyPorosity phaseChange requires meltingTemperature>0, "
            "mushyTemperatureWidth>0 and relaxationCoefficient>0.");
    }
    Math::forInterior(ctx.field, [&](int i, int j, int k) {
        if (ctx.field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double T = ctx.temperature(i, j, k);
        if (!std::isfinite(T) || T <= 0.0) {
            throw std::runtime_error(
                "EnthalpyPorosity phaseChange: invalid temperature.");
        }
        const double superheat = (T - Tm) / width;
        mdot[(size_t)ctx.alpha.getIdx(i, j, k)] +=
            rate * ctx.liquidPhase.density * superheat;
    });
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
