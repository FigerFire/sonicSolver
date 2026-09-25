/// @file SF_tabulatedEOS.h
/// @brief 热力学状态方程模型实现。

#pragma once

#include "SF_EOS.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::Physics::EOS {

struct TabulatedPoint { double rho, e, p, T, c; };

/// @brief 一维能量曲线表；仅允许同一 density 的显式线性插值。
class TabulatedEOS final : public Model {
public:
    explicit TabulatedEOS(std::vector<TabulatedPoint> points)
        : points_(std::move(points)) {
        if (points_.size() < 2) throw std::runtime_error("TabulatedEOS needs at least two points.");
        std::sort(points_.begin(), points_.end(),
                  [](const auto& a, const auto& b) { return a.e < b.e; });
        rho_ = points_.front().rho;
        for (const auto& p : points_) {
            if (!std::isfinite(p.rho) || std::abs(p.rho-rho_) > 1.0e-12*std::max(rho_,1.0)
                || !std::isfinite(p.e) || !std::isfinite(p.p) || p.p <= 0.0
                || !std::isfinite(p.T) || p.T <= 0.0 || !std::isfinite(p.c) || p.c <= 0.0)
                throw std::runtime_error("TabulatedEOS table is invalid or not a single-density curve.");
        }
    }
    std::string type() const override { return "tabulated"; }
    State fromDensityEnergy(double rho, double e) const override {
        if (std::abs(rho-rho_) > 1.0e-12*std::max(rho_,1.0))
            throw std::runtime_error("TabulatedEOS does not extrapolate in density.");
        auto hi = std::lower_bound(points_.begin(), points_.end(), e,
            [](const TabulatedPoint& p, double x) { return p.e < x; });
        if (hi == points_.begin() || hi == points_.end())
            throw std::runtime_error("TabulatedEOS does not extrapolate in energy.");
        const auto& b = *hi; const auto& a = *(hi-1);
        const double w = (e-a.e)/(b.e-a.e);
        return {rho, e, a.p+w*(b.p-a.p), a.T+w*(b.T-a.T), a.c+w*(b.c-a.c)};
    }
    State fromPressureTemperature(double, double) const override {
        throw std::runtime_error("TabulatedEOS p,T inversion requires a 2D table.");
    }
    double referenceEnergy() const override { return points_.front().e; }
private:
    std::vector<TabulatedPoint> points_;
    double rho_ = 0.0;
};

} // namespace SF::Physics::EOS
