/// @file SF_ZGB.cpp
/// @brief 相变模型、饱和性质与守恒传递实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_phaseChange.h"
#include "SF_utility.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {

void addZGBRates(const ModelContext& ctx, std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double Fvap = coefficient(pc, "vaporizationCoefficient", 50.0);
    const double Fcond = coefficient(pc, "condensationCoefficient", 0.01);
    const double alphaNuc = coefficient(pc, "nucleationVolumeFraction", 5.0e-4);
    const double Rb = coefficient(pc, "bubbleRadius", 1.0e-6);
    const auto& vaporPhase = ctx.alphaIsLiquid ? ctx.otherPhase
                                               : ctx.alphaPhase;
    if (Fvap < 0.0 || Fcond < 0.0 || alphaNuc < 0.0 || Rb <= 0.0) {
        throw std::runtime_error(
            "ZGB phaseChange requires non-negative coefficients and bubbleRadius>0.");
    }

    Math::forInterior(ctx.field, [&](int i, int j, int k) {
        if (ctx.field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double p = Numerics::requirePhysicalState(
            "ZGB phaseChange", ctx.field, i, j, k);
        const double a = ctx.alpha(i, j, k);
        const double alphaV = ctx.alphaIsLiquid ? 1.0 - a : a;
        const double evap = p < pc.saturationPressure
            ? Fvap * 3.0 * alphaNuc * (1.0 - alphaV)
              * vaporPhase.density / Rb
              * std::sqrt(2.0 * (pc.saturationPressure - p)
                          / (3.0 * ctx.liquidPhase.density))
            : 0.0;
        const double cond = p > pc.saturationPressure
            ? Fcond * 3.0 * alphaV * vaporPhase.density / Rb
              * std::sqrt(2.0 * (p - pc.saturationPressure)
                          / (3.0 * ctx.liquidPhase.density))
            : 0.0;
        mdot[(size_t)ctx.alpha.getIdx(i, j, k)] += evap - cond;
    });
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
