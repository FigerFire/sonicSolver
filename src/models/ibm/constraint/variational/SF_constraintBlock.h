#pragma once

/// @file SF_constraintBlock.h
/// @brief 与压力/密度算法无关的 IBM `J`/`J^T` 约束块。

#include "SF_immersedConstraint.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace SF::IBM::Variational {

/// @brief 一个表面约束分量到 Eulerian velocity DOF 的耦合项。
struct SurfaceConstraintTerm {
    int marker = -1;
    int component = -1;
    int localCell = -1;
    std::int64_t globalEulerianDof = -1;
    double interpolation = 0.0;
    double adjoint = 0.0;
};

/// @brief 独立于 flow algorithm 的表面 `J` 与伴随 `J^T M_L` 数据块。
///
/// `localCellResolver` 只属于当前离散适配器；该类不读取 Field，不知道压力、
/// EOS、HYPRE row 或 MPI。pressure/density adapter 只消费 terms 并自行组装方程。
class SurfaceConstraintBlock final {
public:
    /// @brief 从统一表面系统生成每个 marker/component 的 J 与 adjoint 项。
    static SurfaceConstraintBlock build(
        const FDM::ImmersedSurfaceSystem& surface,
        const std::vector<int>& velocityComponents,
        const std::function<int(int)>& localCellResolver) {
        if (velocityComponents.empty() || !localCellResolver) {
            throw std::runtime_error(
                "Surface constraint block requires velocity components and resolver.");
        }
        SurfaceConstraintBlock result;
        result.componentCount_=static_cast<int>(velocityComponents.size());
        result.markerCount_=static_cast<int>(surface.points.size());
        result.rows_.resize(static_cast<std::size_t>(
            result.markerCount_*result.componentCount_));
        for (int marker=0;marker<result.markerCount_;++marker) {
            const auto& point=surface.points[static_cast<std::size_t>(marker)];
            if (!point.globalConstraintId.valid()
                || !std::isfinite(point.measure) || point.measure<=0.0) {
                throw std::runtime_error(
                    "Surface constraint block received invalid marker metadata.");
            }
            for (int component=0;component<result.componentCount_;++component) {
                auto& row=result.rows_[static_cast<std::size_t>(
                    marker*result.componentCount_+component)];
                for (const auto& weight:point.interpolation) {
                    if (!std::isfinite(weight.value)
                        || weight.globalEulerianDofId<0) {
                        throw std::runtime_error(
                            "Surface constraint block received invalid J weight.");
                    }
                    const int localCell=localCellResolver(weight.cell);
                    row.push_back({marker,component,localCell,
                        weight.globalEulerianDofId,weight.value,
                        -point.measure*weight.value});
                }
            }
        }
        return result;
    }

    const std::vector<SurfaceConstraintTerm>& row(
            int marker,int component) const {
        if (marker<0 || marker>=markerCount_
            || component<0 || component>=componentCount_) {
            throw std::runtime_error(
                "Surface constraint block row index is outside its layout.");
        }
        return rows_[static_cast<std::size_t>(
            marker*componentCount_+component)];
    }

    const std::vector<std::vector<SurfaceConstraintTerm>>& rows() const {
        return rows_;
    }

private:
    int markerCount_ = 0;
    int componentCount_ = 0;
    std::vector<std::vector<SurfaceConstraintTerm>> rows_;
};

} // namespace SF::IBM::Variational
