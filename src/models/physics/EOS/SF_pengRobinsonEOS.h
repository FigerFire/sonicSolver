/// @file SF_pengRobinsonEOS.h
/// @brief 热力学状态方程模型实现。

#pragma once

#include "SF_EOS.h"

#include <cmath>
#include <stdexcept>

namespace SF::Physics::EOS {

/// @brief Peng-Robinson 的 p(rho,T) 实现；caloric closure 需显式 Cv 和 q。
class PengRobinsonEOS final : public Model {
public:
    PengRobinsonEOS(double gasConstant, double criticalTemperature,
                    double criticalPressure, double acentricFactor,
                    double cv, double referenceEnergy)
        : R_(gasConstant), Tc_(criticalTemperature), Pc_(criticalPressure),
          omega_(acentricFactor), cv_(cv), q_(referenceEnergy) {
        if (R_ <= 0.0 || Tc_ <= 0.0 || Pc_ <= 0.0 || cv_ <= 0.0)
            throw std::runtime_error("PengRobinsonEOS parameters are invalid.");
    }
    std::string type() const override { return "PengRobinson"; }
    State fromDensityEnergy(double rho, double e) const override {
        if (rho <= 0.0 || !std::isfinite(rho) || !std::isfinite(e))
            throw std::runtime_error("PengRobinsonEOS requires rho>0 and finite e.");
        const double T = (e - q_) / cv_;
        if (!std::isfinite(T) || T <= 0.0)
            throw std::runtime_error("PengRobinsonEOS produced T<=0.");
        const double p = pressure(rho, T);
        const double dr = std::max(rho * 1.0e-6, 1.0e-9);
        const double dpdr = (pressure(rho + dr, T) - pressure(rho - dr, T)) / (2.0 * dr);
        if (!std::isfinite(dpdr) || dpdr <= 0.0)
            throw std::runtime_error("PengRobinsonEOS is outside a mechanically stable branch.");
        return {rho, e, p, T, std::sqrt(dpdr)};
    }
    State fromPressureTemperature(double, double) const override {
        throw std::runtime_error(
            "PengRobinsonEOS p,T inversion requires an explicit phase/root selection and is not implicit.");
    }
    double referenceEnergy() const override { return q_; }
private:
    double R_, Tc_, Pc_, omega_, cv_, q_;
    double pressure(double rho, double T) const {
        const double kappa = 0.37464 + 1.54226*omega_ - 0.26992*omega_*omega_;
        const double alpha = std::pow(1.0 + kappa*(1.0-std::sqrt(T/Tc_)), 2.0);
        const double a = 0.45724 * R_*R_*Tc_*Tc_ / Pc_ * alpha;
        const double b = 0.07780 * R_*Tc_ / Pc_;
        const double v = 1.0 / rho;
        const double d1 = v - b;
        const double d2 = v*v + 2.0*b*v - b*b;
        if (d1 <= 0.0 || d2 <= 0.0)
            throw std::runtime_error("PengRobinsonEOS density is outside its volume domain.");
        return R_*T/d1 - a/d2;
    }
};

} // namespace SF::Physics::EOS
