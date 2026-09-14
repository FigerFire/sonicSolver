/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_applicator.h
/// @brief 将边界配置值应用到 Field 的 boundary 模块适配器。
///
/// 低层边界核仍由 `SF_boundary.*`、`algebraic/` 和 `ILW/` 提供。
/// 本类只保存一份 `FDM::BoundaryConfig`，并按求解步骤统一调用这些边界核，
/// 避免 Solver Algorithm 直接读取 parser 全局边界数组。

#include "SF_boundary.h"
#include "SF_config.h"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace SF::Boundary {

/// @brief Applies configured density, velocity, and pressure-derived energy boundaries.
///
/// `Applicator` owns a copy of `BoundaryConfig`. This is deliberate:
/// once a Solver Algorithm is constructed, its boundary behavior is deterministic and
/// independent of later mutations to global parser state.
class Applicator final {
public:
    /// @brief Construct an empty applicator.
    ///
    /// Applying an empty applicator is a no-op, which is useful for tests or
    /// cases with periodic/externally synchronized boundaries.
    Applicator() = default;

    /// @brief Construct from parsed boundary settings.
    /// @param config Boundary settings copied into the applicator.
    explicit Applicator(FDM::BoundaryConfig config)
        : config_(std::move(config)) {}

    /// @brief Apply all configured physical boundary conditions.
    ///
    /// Density and velocity are applied first. Pressure input is then converted
    /// into conservative energy so it can use the updated density/momentum state.
    ///
    /// @param field Field whose real boundary points and ghost cells are updated.
    void apply(Field& field) {
        markSolverBoundaryMask(field);
        if (field.hasEquationSet() && config_.ilwEnabled
            && !field.equationSet()->perfectGasGamma()) {
            throw std::runtime_error(
                "EquationSet physical-boundary ILW currently requires a "
                "PerfectGas thermodynamic capability.");
        }
        if (field.hasEquationSet()) {
            applyEquationSetDensity(field);
        } else {
            SF::Boundary::update<double>(
                field, config_.density, RHO,
                config_.ilwEnabled, config_.ilwOrder);
        }
        const int momentum = field.hasEquationSet()
            ? field.equationSet()->momentumIndex(0) : RU;
        SF::Boundary::update<Vector3>(
            field, config_.velocity, momentum,
            config_.ilwEnabled, config_.ilwOrder);
        SF::Boundary::updateEnergyFromPressure(
            field, config_.energyFromPressure,
            config_.ilwEnabled, config_.ilwOrder);
        SF::Boundary::updateEnergyFromThermalBoundary(
            field, config_.thermal, config_.thermalDynamicViscosity,
            config_.thermalPrandtl, config_.ilwEnabled);
    }

private:
    static double totalDensity(const Field& field,
                               int i, int j, int k) {
        double rho = 0.0;
        const int count = field.equationSet()->densityVariableCount();
        for (int v = 0; v < count; ++v) rho += field(i,j,k,v);
        if (!std::isfinite(rho) || rho <= 0.0) {
            throw std::runtime_error(
                "EquationSet density boundary found non-positive total density.");
        }
        return rho;
    }

    static void setTotalDensityFromComposition(
            Field& field,
            int i, int j, int k,
            int compositionI, int compositionJ, int compositionK,
            double targetDensity) {
        if (!std::isfinite(targetDensity) || targetDensity <= 0.0) {
            throw std::runtime_error(
                "EquationSet fixed-density ghost state is non-positive; "
                "the prescribed boundary state is incompatible with its interior reflection.");
        }
        const double referenceDensity = totalDensity(
            field, compositionI, compositionJ, compositionK);
        const int count = field.equationSet()->densityVariableCount();
        for (int v = 0; v < count; ++v) {
            const double fraction =
                field(compositionI,compositionJ,compositionK,v)
                / referenceDensity;
            if (!std::isfinite(fraction) || fraction < 0.0) {
                throw std::runtime_error(
                    "EquationSet boundary composition contains an invalid partial density.");
            }
            field(i,j,k,v) = targetDensity*fraction;
        }
    }

    void applyEquationSetDensity(Field& field) const {
        const int count = field.equationSet()->densityVariableCount();
        for (const auto& bc : config_.density) {
            // A one-density PerfectGas EquationSet has the same conservative
            // density variable as the established ILW kernel.  Reuse that
            // kernel so binding thermodynamics does not alter the boundary
            // reconstruction or its stage/order.
            if (count == 1) {
                SF::Boundary::update<double>(
                    field, std::vector<BCSetting<double>>{bc}, RHO,
                    config_.ilwEnabled, config_.ilwOrder);
                continue;
            }
            if (bc.type != FIXED_VALUE) {
                const std::vector<BCSetting<double>> one{bc};
                for (int v = 0; v < count; ++v)
                    SF::Boundary::update<double>(field, one, v, false);
                continue;
            }
            forBoundarySettingPoints(
                field, bc, [&](int i, int j, int k, int axis) {
                    setTotalDensityFromComposition(
                        field, i,j,k, i,j,k, bc.value);
                    auto applyFixed = [&](int gi, int gj, int gk,
                                          int ri, int rj, int rk) {
                        const double reflectedDensity =
                            2.0*bc.value-totalDensity(field,ri,rj,rk);
                        setTotalDensityFromComposition(
                            field, gi,gj,gk, ri,rj,rk,
                            reflectedDensity);
                    };
                    forBoundaryGhostsAlongAxis(
                        field, i,j,k, axis, applyFixed);
                });
        }
    }

    /// @brief 根据BC配置声明真实物理边界求解mask。
    /// @param field 将被写入mask的Field。
    void markSolverBoundaryMask(Field& field) const {
        field.clearSolverBoundaryMask();
        auto mark = [&](const auto& bc) {
            if (bc.type == EMPTY) return;
            const auto& sets = field.getAllSets();
            if (sets.find(bc.name) != sets.end()) {
                field.markSolverBoundarySet(bc.name);
            }
        };
        for (const auto& bc : config_.density) {
            mark(bc);
        }
        for (const auto& bc : config_.velocity) {
            mark(bc);
        }
        for (const auto& bc : config_.energyFromPressure) {
            mark(bc);
        }
        for (const auto& bc : config_.thermal) {
            if (bc.type == ThermalBCType::Empty) continue;
            const auto& sets = field.getAllSets();
            if (sets.find(bc.name) != sets.end()) {
                field.markSolverBoundarySet(bc.name);
            }
        }
    }

    FDM::BoundaryConfig config_;
};

} // namespace SF::Boundary
