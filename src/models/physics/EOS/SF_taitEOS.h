/// @file SF_taitEOS.h
/// @brief 热力学状态方程模型实现。

#pragma once

#include "SF_EOS.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::Physics::EOS {

/// @brief barotropic Tait EOS: p=B[(rho/rho0)^n-1]+p0。
class TaitEOS final : public Model {
public:
    TaitEOS(double rho0, double p0, double bulk, double exponent,
            double referenceEnergy)
        : rho0_(rho0), p0_(p0), bulk_(bulk), exponent_(exponent),
          referenceEnergy_(referenceEnergy) {
        if (rho0_ <= 0.0 || bulk_ <= 0.0 || exponent_ <= 1.0
            || !std::isfinite(p0_) || !std::isfinite(referenceEnergy_))
            throw std::runtime_error("TaitEOS parameters are invalid.");
    }
    std::string type() const override { return "Tait"; }
    State fromDensityEnergy(double rho, double e) const override {
        if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(e))
            throw std::runtime_error("TaitEOS requires rho>0 and finite e.");
        const double ratio = rho / rho0_;
        const double p = bulk_ * (std::pow(ratio, exponent_) - 1.0) + p0_;
        const double c = std::sqrt(bulk_ * exponent_
                                   * std::pow(ratio, exponent_ - 1.0) / rho0_);
        return {rho, e, p, std::numeric_limits<double>::quiet_NaN(), c};
    }
    State fromPressureTemperature(double p, double) const override {
        const double base = 1.0 + (p - p0_) / bulk_;
        if (!std::isfinite(base) || base <= 0.0)
            throw std::runtime_error("TaitEOS pressure is outside its domain.");
        return fromDensityEnergy(rho0_ * std::pow(base, 1.0 / exponent_),
                                 referenceEnergy_);
    }
    double referenceEnergy() const override { return referenceEnergy_; }
private:
    double rho0_, p0_, bulk_, exponent_, referenceEnergy_;
};

} // namespace SF::Physics::EOS
