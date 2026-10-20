/// @file SF_RPI.cpp
/// @brief RPI 壁面沸腾闭式模型与热流分配实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_RPI.h"

#include "SF_closure.h"
#include "SF_heatFluxPartition.h"
#include "SF_wallMapping.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

void validateConfig(const Multiphase::PhaseChangeOptions& config,
                    const std::string& context) {
    validateClosureConfig(config, context);
    validateHeatFluxConfig(config, context);
}

void computeRates(const ModelContext& ctx, std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    for (const auto& setting : pc.wallBoiling) {
        if (!std::isfinite(setting.heatFlux) || setting.heatFlux <= 0.0) {
            throw std::runtime_error(
                "RPI phaseChange requires finite positive wall heatFlux; "
                "set it on the matching 0/T wallHeatFlux patch or in phaseChange.");
        }

        const auto samples = mapWallPatch(ctx, setting);
        for (const WallCellSample& sample : samples) {
            const double wallTemperature =
                std::isfinite(setting.wallTemperature)
                ? setting.wallTemperature
                : ctx.temperature(
                      sample.wallI, sample.wallJ, sample.wallK);
            const double liquidTemperature =
                ctx.temperature(sample.fluidI, sample.fluidJ, sample.fluidK);
            if (!std::isfinite(wallTemperature) || wallTemperature <= 0.0
                || !std::isfinite(liquidTemperature)
                || liquidTemperature <= 0.0) {
                throw std::runtime_error(
                    "RPI phaseChange invalid wall/near-wall temperature.");
            }
            if (wallTemperature <= pc.saturationTemperature) {
                continue;
            }

            const WallHeatBalance balance = evaluateWallHeatBalance(
                ctx, liquidTemperature, wallTemperature, setting.heatFlux);
            const BubbleClosure& closure = balance.closure;
            const HeatFluxPartition& q = balance.heat;
            if (budgetCheckEnabled(pc)
                && q.total > setting.heatFlux * (1.0 + 1.0e-10)) {
                std::ostringstream oss;
                oss << "RPI phaseChange heat partition q_total=" << q.total
                    << " exceeds wall heatFlux=" << setting.heatFlux
                    << " on patch '" << setting.patch << "' (q_c="
                    << q.convective << ", q_q=" << q.quenching
                    << ", q_e=" << q.evaporative << ").";
                throw std::runtime_error(oss.str());
            }

            const double volumetricMdot =
                closure.wallMassFlux * sample.areaOverVolume;
            if (!std::isfinite(volumetricMdot) || volumetricMdot <= 0.0) {
                throw std::runtime_error(
                    "RPI phaseChange volumetric source is invalid.");
            }
            mdot[(size_t)ctx.alpha.getIdx(sample.fluidI,
                                          sample.fluidJ,
                                          sample.fluidK)] += volumetricMdot;
        }
    }
}

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
