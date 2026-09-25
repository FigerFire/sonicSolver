/// @file SF_lift.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems::Interphase {

void addLift(
        const PairContext& pair,
        PhaseEquationSources& sources) {
    const auto& options = pair.options;
    if (Multiphase::normalizeModelType(options.liftModel)
        != "constantcoefficient") {
        throw std::runtime_error(
            "Eulerian lift supports only constantCoefficient.");
    }
    const Context& context = pair.system;
    const auto& continuous = context.phases[pair.continuous];
    const auto& dispersed = context.phases[pair.dispersed];
    const int ghost = context.geometry.NG();
    for (int k = ghost; k < ghost + context.geometry.NZ(); ++k) {
        for (int j = ghost; j < ghost + context.geometry.NY(); ++j) {
            for (int i = ghost; i < ghost + context.geometry.NX(); ++i) {
                if (context.geometry.CellFlag(i,j,k) != FLUID_CELL) continue;
                const auto vorticity = CENTRAL2::curl(
                    context.geometry, continuous.primitive.velocity, i, j, k);
                const int cell = continuous.primitive.alpha.getIdx(i,j,k);
                const double factor = options.liftCoefficient
                    * continuous.primitive.density.values()[(size_t)cell]
                    * dispersed.primitive.alpha.values()[(size_t)cell];
                std::array<double, 3> relative{};
                for (int component = 0; component < 3; ++component) {
                    relative[(size_t)component] =
                        dispersed.primitive.velocity[(size_t)component]
                            .values()[(size_t)cell]
                        - continuous.primitive.velocity[(size_t)component]
                            .values()[(size_t)cell];
                }
                // 原式给出分散相受力 C_L rho alpha (Ur x curl(Uc))；
                // 统一接口接收连续相受力，因此取反。
                const std::array<double, 3> forceOnContinuous{
                    -factor * (relative[1] * vorticity[2]
                               - relative[2] * vorticity[1]),
                    -factor * (relative[2] * vorticity[0]
                               - relative[0] * vorticity[2]),
                    -factor * (relative[0] * vorticity[1]
                               - relative[1] * vorticity[0])};
                addConservativePairForce(
                    pair, sources, cell, forceOnContinuous);
            }
        }
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
