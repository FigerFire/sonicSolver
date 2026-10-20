/// @file SF_phaseField.cpp
/// @brief 相变模型、饱和性质与守恒传递实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_phaseChange.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {

void addPhaseFieldRates(const ModelContext& ctx, std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double mobility = coefficient(pc, "mobility", 0.0);
    const double areaDensity = coefficient(pc, "interfaceAreaDensity", 1.0);
    if (mobility <= 0.0 || areaDensity <= 0.0) {
        throw std::runtime_error(
            "PhaseField phaseChange requires mobility>0 and interfaceAreaDensity>0.");
    }

    Math::forInterior(ctx.field, [&](int i, int j, int k) {
        if (ctx.field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double a = ctx.alpha(i, j, k);
        const double T = ctx.temperature(i, j, k);
        if (!std::isfinite(a) || !std::isfinite(T) || T <= 0.0) {
            throw std::runtime_error("PhaseField phaseChange: invalid alpha/T.");
        }
        const double alphaL = ctx.alphaIsLiquid ? a : 1.0 - a;
        const double alphaV = ctx.alphaIsLiquid ? 1.0 - a : a;
        const double interfaceWindow =
            std::max(0.0, alphaL) * std::max(0.0, alphaV);
        const double drive =
            (T - pc.saturationTemperature) / pc.saturationTemperature;
        mdot[(size_t)ctx.alpha.getIdx(i, j, k)] +=
            mobility * areaDensity * ctx.liquidPhase.density
            * interfaceWindow * drive;
    });
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
