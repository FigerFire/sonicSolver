/// @file SF_SchnerrSauer.cpp
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

void addSchnerrSauerRates(const ModelContext& ctx,
                          std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double n0 = coefficient(pc, "bubbleNumberDensity", 1.0e12);
    const double coeff = coefficient(pc, "coefficient", 1.0);
    const auto& vaporPhase = ctx.alphaIsLiquid ? ctx.otherPhase
                                               : ctx.alphaPhase;
    if (n0 <= 0.0 || coeff < 0.0) {
        throw std::runtime_error(
            "SchnerrSauer phaseChange requires bubbleNumberDensity>0 and coefficient>=0.");
    }

    Math::forInterior(ctx.field, [&](int i, int j, int k) {
        if (ctx.field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double p = Numerics::requirePhysicalState(
            "SchnerrSauer phaseChange", ctx.field, i, j, k);
        const double a = ctx.alpha(i, j, k);
        const double alphaL = ctx.alphaIsLiquid ? a : 1.0 - a;
        const double alphaV = ctx.alphaIsLiquid ? 1.0 - a : a;
        const double limitedV =
            std::max(1.0e-12, std::min(1.0 - 1.0e-12, alphaV));
        constexpr double pi = 3.141592653589793238462643383279502884;
        const double bubbleRadius =
            std::pow(3.0 * limitedV / (4.0 * pi * n0), 1.0 / 3.0);
        const double dp = pc.saturationPressure - p;
        const double speed =
            std::sqrt(2.0 * std::abs(dp) / (3.0 * ctx.liquidPhase.density));
        const double rhoMix =
            alphaL * ctx.liquidPhase.density
            + alphaV * vaporPhase.density;
        if (rhoMix <= 0.0 || bubbleRadius <= 0.0) {
            throw std::runtime_error(
                "SchnerrSauer phaseChange produced invalid mixture state.");
        }
        const double magnitude =
            coeff * (ctx.liquidPhase.density * vaporPhase.density / rhoMix)
            * (alphaL * alphaV / bubbleRadius) * speed;
        mdot[(size_t)ctx.alpha.getIdx(i, j, k)] += dp > 0.0
            ? magnitude : -magnitude;
    });
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
