/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_ghostCell.h
/// @brief IBM ghost-cell 守恒量刷新。

#include "SF_weightBuilder.h"
#include "SF_slipWall.h"

#include <cstdlib>
#include <iostream>

namespace SF {
namespace IBM {
namespace GhostIBM {

/// @brief 根据 ghost-cell 权重刷新 IBM ghost 单元守恒状态。
class GhostCellIBM {
public:
    /// @brief 使用ILW边界闭合刷新所有IBM ghost-cell值。
    /// @param field 结构网格场；IBM_GHOST_CELL中的守恒量会被覆盖。
    /// @param weights WeightBuilder生成的ghost点插值权重。
    /// @param config application 层冻结的 IBM/ILW 配置。
    void apply(Field& field,
               const WeightBuilder& weights,
               const IBMRuntimeConfig& config) const {
        const int requestedTaylorOrder = config.requestedTaylorOrder();
        for (const auto& w : weights.weights()) {
            double qFluid[5] = {
                interpolate(field, w, RHO),
                interpolate(field, w, RU),
                interpolate(field, w, RV),
                interpolate(field, w, RW),
                interpolate(field, w, E)
            };
            double normal[3] = {w.wallNormal.x, w.wallNormal.y, w.wallNormal.z};
            double qGhost[5];

            if (config.ilwEnabled) {
                if (!SF::IBM::GhostILW::buildGhostStateForCell(
                        field, weights.ilwStorage(), weights.geometry(),
                        w.ghostI, w.ghostJ, w.ghostK,
                        Math::XI, requestedTaylorOrder,
                        config, qGhost)) {
                    const size_t cell =
                        (size_t)field.getIdx(w.ghostI, w.ghostJ, w.ghostK);
                    const auto& storage = weights.ilwStorage();
                    const bool hasPlan =
                        storage.hasPlan(cell);
                    const auto* plan = hasPlan
                        ? &storage.plan(cell)
                        : nullptr;
                    std::cerr
                        << "[SF FATAL] ILW ghost-cell closure failed at ("
                        << w.ghostI << "," << w.ghostJ << ","
                        << w.ghostK << "), flag="
                        << field.CellFlag(w.ghostI, w.ghostJ, w.ghostK)
                        << ", requestedTaylorOrder=" << requestedTaylorOrder
                        << ", fluidSamples="
                        << storage.fluidCount(cell)
                        << ", normalSamples="
                        << storage.normalCount(cell)
                        << ", hasPlan=" << hasPlan
                        << ", useT2=" << (plan ? plan->useT2 : false)
                        << ", planMaxOrder=" << (plan ? plan->maxOrder : 0)
                        << ", higherFits=" << (plan ? plan->higherFits.size() : 0)
                        << ", targetDistance=" << (plan ? plan->targetDistance : 0.0)
                        << ". LowOrder fallback is disabled; "
                        << "set [numerics].ILW=0 and select a low-order IBM "
                        << "scheme explicitly if that behavior is desired."
                        << std::endl;
                    std::exit(1);
                }
            } else {
                SF::IBM::GhostILW::firstOrderEulerWallGhostState(
                    qFluid, normal, config.gamma, qGhost);
            }

            field(w.ghostI, w.ghostJ, w.ghostK, RHO) = qGhost[0];
            field(w.ghostI, w.ghostJ, w.ghostK, RU) = qGhost[1];
            field(w.ghostI, w.ghostJ, w.ghostK, RV) = qGhost[2];
            field(w.ghostI, w.ghostJ, w.ghostK, RW) = qGhost[3];
            field(w.ghostI, w.ghostJ, w.ghostK, E) = qGhost[4];
        }
    }

private:
    static double interpolate(const Field& field, const GhostCellWeight& weight, int var) {
        double value = 0.0;
        for (const auto& donor : weight.donors) {
            value += donor.weight * field(donor.i, donor.j, donor.k, var);
        }
        return value;
    }
};

} // namespace GhostIBM
} // namespace IBM
} // namespace SF
