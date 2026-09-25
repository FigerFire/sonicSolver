/// @file SF_perfectGasEOS.h
/// @brief 热力学状态方程模型实现。

#pragma once

#include "SF_EOS.h"

#include <cmath>
#include <stdexcept>

namespace SF::Physics::EOS {

class PerfectGasEOS final : public Model {
public:
    PerfectGasEOS(double gamma, double gasConstant, double referenceEnergy = 0.0)
        : gamma_(gamma), gasConstant_(gasConstant), referenceEnergy_(referenceEnergy) {
        if (!std::isfinite(gamma_) || gamma_ <= 1.0
            || !std::isfinite(gasConstant_) || gasConstant_ <= 0.0
            || !std::isfinite(referenceEnergy_)) {
            throw std::runtime_error("PerfectGasEOS requires gamma>1, R>0 and finite q.");
        }
    }

    std::string type() const override { return "perfectGas"; }

    State fromDensityEnergy(double rho, double e) const override {
        requireDensity(rho);
        const double thermalEnergy = e - referenceEnergy_;
        const double p = (gamma_ - 1.0) * rho * thermalEnergy;
        const double T = thermalEnergy * (gamma_ - 1.0) / gasConstant_;
        return checked(rho, e, p, T,
                       std::sqrt(gamma_ * p / rho));
    }

    State fromPressureTemperature(double p, double T) const override {
        if (!std::isfinite(p) || p <= 0.0 || !std::isfinite(T) || T <= 0.0)
            throw std::runtime_error("PerfectGasEOS requires p>0 and T>0.");
        const double rho = p / (gasConstant_ * T);
        const double e = referenceEnergy_ + gasConstant_ * T / (gamma_ - 1.0);
        return checked(rho, e, p, T, std::sqrt(gamma_ * p / rho));
    }

    State fromDensityPressure(double rho, double p) const override {
        requireDensity(rho);
        if (!std::isfinite(p)||p<=0.0) {
            throw std::runtime_error(
                "PerfectGasEOS requires positive finite pressure.");
        }
        const double thermalEnergy=p/((gamma_-1.0)*rho);
        const double e=referenceEnergy_+thermalEnergy;
        const double T=p/(rho*gasConstant_);
        return checked(rho,e,p,T,std::sqrt(gamma_*p/rho));
    }

    State fromDensityTemperature(double rho, double T) const override {
        requireDensity(rho);
        if (!std::isfinite(T)||T<=0.0) {
            throw std::runtime_error(
                "PerfectGasEOS requires positive finite temperature.");
        }
        const double p=rho*gasConstant_*T;
        const double e=referenceEnergy_
            +gasConstant_*T/(gamma_-1.0);
        return checked(rho,e,p,T,std::sqrt(gamma_*p/rho));
    }

    double referenceEnergy() const override { return referenceEnergy_; }
    double gamma() const { return gamma_; }
    double gasConstant() const { return gasConstant_; }

private:
    double gamma_;
    double gasConstant_;
    double referenceEnergy_;

    static void requireDensity(double rho) {
        if (!std::isfinite(rho) || rho <= 0.0)
            throw std::runtime_error("PerfectGasEOS requires rho>0.");
    }
    static State checked(double rho, double e, double p, double T, double c) {
        if (!std::isfinite(e) || !std::isfinite(p) || p <= 0.0
            || !std::isfinite(T) || T <= 0.0
            || !std::isfinite(c) || c <= 0.0)
            throw std::runtime_error("PerfectGasEOS produced an invalid state.");
        return {rho, e, p, T, c};
    }
};

} // namespace SF::Physics::EOS
