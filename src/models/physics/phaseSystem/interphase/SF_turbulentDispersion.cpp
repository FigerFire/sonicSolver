/// @file SF_turbulentDispersion.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems::Interphase {

void addTurbulentDispersion(
        const PairContext& pair,
        PhaseEquationSources& sources) {
    const std::string model = Multiphase::normalizeModelType(
        pair.options.turbulentDispersionModel);
    if (model == "none") return;
    if (model != "gdb" && model != "favreaveraged") {
        throw std::runtime_error(
            "Eulerian turbulentDispersion supports only none or GDB.");
    }

    const Context& context = pair.system;
    const auto& dispersed = context.phases[pair.dispersed];
    const auto& drag = sources.momentumCouplings[pair.couplingIndex]
                           .dragCoefficient;
    const double diffusivity =
        pair.options.turbulentDispersionCoefficient
        * pair.options.turbulentKinematicViscosity
        / pair.options.turbulentSchmidtNumber;
    const int ghost = context.geometry.NG();
    for (int k = ghost; k < ghost + context.geometry.NZ(); ++k) {
        for (int j = ghost; j < ghost + context.geometry.NY(); ++j) {
            for (int i = ghost; i < ghost + context.geometry.NX(); ++i) {
                if (context.geometry.CellFlag(i,j,k) != FLUID_CELL
                    || context.geometry.isSolverBoundaryPoint(i,j,k)) {
                    continue;
                }
                const int cell = context.geometry.getIdx(i,j,k);
                const double alpha =
                    dispersed.primitive.alpha.values()[(size_t)cell];
                if (!std::isfinite(alpha) || alpha <= 0.0) {
                    throw std::runtime_error(
                        "GDB turbulentDispersion requires alpha_d > 0.");
                }
                const auto alphaGradient =
                    CENTRAL2::grad(
                        context.geometry, dispersed.primitive.alpha, i, j, k);
                std::array<double, 3> forceOnContinuous{};
                for (int component = 0; component < 3; ++component) {
                    const double dispersionVelocity =
                        -diffusivity * alphaGradient[(size_t)component] / alpha;
                    // 分散相受力为 beta*Vtd；统一接口接收连续相受力。
                    forceOnContinuous[(size_t)component] =
                        -drag.values()[(size_t)cell] * dispersionVelocity;
                }
                addConservativePairForce(
                    pair, sources, cell, forceOnContinuous);
            }
        }
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
