/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

/// @file SF_massTransfer.h
/// @brief Lee 模型质量传递速率计算模块。
///
/// 使用 Lee 蒸发/冷凝模型计算质量传递源项：
///   蒸发 (T > Tsat): mdot = Ce * alpha_l * rho_l * (T - Tsat) / Tsat
///   冷凝 (T < Tsat): mdot = Cc * alpha_v * rho_v * (Tsat - T) / Tsat
///
/// 调用 numerics 和 math 模块进行状态检查和迭代。

#pragma once

#include "SF_phaseChange.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_utility.h"
#include "SF_wallMapping.h"

#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace Lee {

namespace Detail {

inline std::string leeSourcePatch(
        const Multiphase::PhaseChangeOptions& pc) {
    std::string patch = selection(pc, "sourcePatch", "");
    if (patch.empty()) patch = selection(pc, "wallPatch", "");
    if (patch.empty()) patch = selection(pc, "patch", "");
    if (patch.empty()) {
        const std::string multiPatch = selection(pc, "sourcePatches", "");
        const std::string wallPatches = selection(pc, "wallPatches", "");
        const std::string patches = selection(pc, "patches", "");
        const std::string region = selection(pc, "sourceRegion", "");
        if (!multiPatch.empty() || !wallPatches.empty()
            || !patches.empty() || !region.empty()) {
            throw std::runtime_error(
                "Lee::addMassSource: only single sourcePatch/wallPatch/patch "
                "selection is implemented.");
        }
    }
    return patch;
}

inline int leeSourceLayers(
        const Multiphase::PhaseChangeOptions& pc) {
    const double raw = coefficient(pc, "sourceLayers", 1.0);
    if (!std::isfinite(raw) || raw < 1.0) {
        throw std::runtime_error(
            "Lee::addMassSource: sourceLayers must be finite and >= 1.");
    }
    const int layers = static_cast<int>(std::llround(raw));
    if (std::abs(raw - (double)layers) > 1.0e-12) {
        throw std::runtime_error(
            "Lee::addMassSource: sourceLayers must be an integer.");
    }
    return layers;
}

inline bool isAllPatchName(std::string patch) {
    for (char& c : patch) c = (char)std::tolower((unsigned char)c);
    return patch == "all";
}

inline std::vector<unsigned char> buildSourceMask(const ModelContext& ctx) {
    const auto& pc = ctx.config.phaseChange;
    const std::string patch = leeSourcePatch(pc);
    if (patch.empty()) return {};
    if (isAllPatchName(patch)) {
        throw std::runtime_error(
            "Lee::addMassSource: sourcePatch=all is not a wall-localized "
            "selection; omit sourcePatch for bulk Lee source.");
    }

    WallHeatSetting setting;
    setting.patch = patch;
    const int layers = leeSourceLayers(pc);
    const auto samples = RPI::mapWallPatch(ctx, setting);
    std::vector<unsigned char> mask((size_t)ctx.field.TotalSize(), 0);
    int applied = 0;

    for (const auto& sample : samples) {
        const int di = sample.fluidI - sample.wallI;
        const int dj = sample.fluidJ - sample.wallJ;
        const int dk = sample.fluidK - sample.wallK;
        if (std::abs(di) + std::abs(dj) + std::abs(dk) != 1) {
            throw std::runtime_error(
                "Lee::addMassSource: invalid wall-normal source direction.");
        }

        for (int layer = 1; layer <= layers; ++layer) {
            const int i = sample.wallI + layer * di;
            const int j = sample.wallJ + layer * dj;
            const int k = sample.wallK + layer * dk;
            if (ctx.field.isSolverBoundaryPoint(i, j, k)
                || ctx.field.CellFlag(i, j, k) != FLUID_CELL) {
                throw std::runtime_error(
                    "Lee::addMassSource: sourceLayers extends outside solved "
                    "fluid cells next to patch '" + patch + "'.");
            }
            const int id = ctx.alpha.getIdx(i, j, k);
            if (mask[(size_t)id] == 0) {
                mask[(size_t)id] = 1;
                ++applied;
            }
        }
    }

    if (applied == 0) {
        throw std::runtime_error(
            "Lee::addMassSource: sourcePatch '" + patch
            + "' did not select any solved fluid cells.");
    }
    return mask;
}

} // namespace Detail

/// @brief 计算一个网格单元上的 Lee 模型质量传递速率。
///
/// @param alphaL 液相体积分数。
/// @param alphaV 气相体积分数。
/// @param T 局部温度 (K)。
/// @param Tsat 饱和温度 (K)。
/// @param rhoL 液相密度 (kg/m^3)。
/// @param rhoV 气相密度 (kg/m^3)。
/// @param Ce 蒸发系数 (1/s)。
/// @param Cc 冷凝系数 (1/s)。
/// @return mdot (kg/(m^3*s))，正值为 liquid->vapor。
inline double cellMassSource(double alphaL, double alphaV,
                             double T, double Tsat,
                             double rhoL, double rhoV,
                             double Ce, double Cc) {
    if (!std::isfinite(alphaL) || !std::isfinite(alphaV) ||
        !std::isfinite(T) || !std::isfinite(Tsat) ||
        !std::isfinite(rhoL) || !std::isfinite(rhoV) ||
        !std::isfinite(Ce) || !std::isfinite(Cc)) {
        throw std::runtime_error(
            "Lee::cellMassSource: non-finite input.");
    }
    if (rhoL <= 0.0 || rhoV <= 0.0 || Tsat <= 0.0) {
        throw std::runtime_error(
            "Lee::cellMassSource: density and Tsat must be > 0.");
    }

    double mdot = 0.0;
    if (T > Tsat && Ce > 0.0) {
        // 蒸发: liquid -> vapor
        mdot = Ce * alphaL * rhoL * (T - Tsat) / Tsat;
    } else if (T < Tsat && Cc > 0.0) {
        // 冷凝: vapor -> liquid
        mdot = -Cc * alphaV * rhoV * (Tsat - T) / Tsat;
    }

    if (!std::isfinite(mdot)) {
        throw std::runtime_error(
            "Lee::cellMassSource: non-finite mdot.");
    }
    return mdot;
}

/// @brief 遍历流体域填充 mdot 数组。
///
/// 使用 Math::forFluidInterior 遍历，用 Numerics::requirePhysicalState
/// 获取局部压力，再用 saturationTemperature 计算 Tsat，
/// 最后用 cellMassSource 计算质量传递速率。
///
/// @param ctx 相变模型上下文。
/// @param mdot 输出的质量源项数组 (kg/(m^3*s))。
inline void addMassSource(const ModelContext& ctx,
                          std::vector<double>& mdot) {
    const auto& pc = ctx.config.phaseChange;
    const double TsatConfig = pc.saturationTemperature;
    const double Ce = pc.evaporationCoefficient;
    const double Cc = pc.condensationCoefficient;

    if (Ce <= 0.0 && Cc <= 0.0) {
        throw std::runtime_error(
            "Lee::addMassSource: both Ce and Cc are zero.");
    }

    const std::vector<unsigned char> sourceMask =
        Detail::buildSourceMask(ctx);

    Math::forFluidInterior(ctx.field, [&](int i, int j, int k) {
        const int id = ctx.alpha.getIdx(i, j, k);
        if (!sourceMask.empty()
            && sourceMask[(size_t)id] == 0) {
            return;
        }

        const double alpha = ctx.alpha(i, j, k);
        const double T = ctx.temperature(i, j, k);
        if (!std::isfinite(alpha) || !std::isfinite(T) || T <= 0.0) {
            throw std::runtime_error(
                "Lee::addMassSource: invalid alpha/T at ("
                + std::to_string(i) + "," + std::to_string(j) + ","
                + std::to_string(k) + ").");
        }

        // 局部压力（用守恒量计算）
        const double p = Numerics::requirePhysicalState(
            "Lee::addMassSource", ctx.field, i, j, k);

        // 饱和温度（压力依赖或常值）
        const double Tsat = saturationTemperature(p, pc);
        if (!std::isfinite(Tsat) || Tsat <= 0.0) {
            throw std::runtime_error(
                "Lee::addMassSource: invalid Tsat at ("
                + std::to_string(i) + "," + std::to_string(j) + ","
                + std::to_string(k) + ").");
        }

        const double alphaL = ctx.alphaIsLiquid ? alpha : 1.0 - alpha;
        const double alphaV = ctx.alphaIsLiquid ? 1.0 - alpha : alpha;
        const double rhoL = ctx.liquidPhase.density;
        const double rhoV = ctx.otherPhase.density;

        const double rate = cellMassSource(alphaL, alphaV, T, Tsat,
                                           rhoL, rhoV, Ce, Cc);
        mdot[(size_t)id] += rate;
    });
}

} // namespace Lee
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
