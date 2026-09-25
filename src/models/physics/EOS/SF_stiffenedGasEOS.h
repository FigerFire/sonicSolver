/// @file SF_stiffenedGasEOS.h
/// @brief 热力学状态方程模型实现。

#pragma once

#include "SF_EOS.h"

#include <cmath>
#include <stdexcept>

namespace SF::Physics::EOS {

/// @brief stiffened-gas: p=(gamma-1)rho(e-q)-gamma*pInf。
class StiffenedGasEOS final : public Model {
public:
    StiffenedGasEOS(double gamma, double cv, double pInfinity,
                    double referenceEnergy)
        : gamma_(gamma), cv_(cv), pInfinity_(pInfinity),
          referenceEnergy_(referenceEnergy) {
        if (!std::isfinite(gamma_) || gamma_ <= 1.0
            || !std::isfinite(cv_) || cv_ <= 0.0
            || !std::isfinite(pInfinity_) || pInfinity_ < 0.0
            || !std::isfinite(referenceEnergy_))
            throw std::runtime_error("StiffenedGasEOS parameters are invalid.");
    }

    std::string type() const override { return "stiffenedGas"; }

    State fromDensityEnergy(double rho, double e) const override {
        if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(e))
            throw std::runtime_error("StiffenedGasEOS requires rho>0 and finite e.");
        const double p = (gamma_ - 1.0) * rho * (e - referenceEnergy_)
                       - gamma_ * pInfinity_;
        const double T = (e - referenceEnergy_ - pInfinity_ / rho) / cv_;
        const double c2 = gamma_ * (p + pInfinity_) / rho;
        if (!std::isfinite(p) || p + pInfinity_ <= 0.0
            || !std::isfinite(T) || T <= 0.0
            || !std::isfinite(c2) || c2 <= 0.0)
            throw std::runtime_error("StiffenedGasEOS produced an invalid state.");
        return {rho, e, p, T, std::sqrt(c2)};
    }

    State fromPressureTemperature(double p, double T) const override {
        if (!std::isfinite(p) || p + pInfinity_ <= 0.0
            || !std::isfinite(T) || T <= 0.0)
            throw std::runtime_error("StiffenedGasEOS requires p+pInf>0 and T>0.");
        const double rho = (p + pInfinity_) / ((gamma_ - 1.0) * cv_ * T);
        const double e = referenceEnergy_ + cv_ * T + pInfinity_ / rho;
        const double c = std::sqrt(gamma_ * (p + pInfinity_) / rho);
        return {rho, e, p, T, c};
    }

    State fromDensityPressure(double rho, double p) const override {
        if (!std::isfinite(rho)||rho<=0.0
            ||!std::isfinite(p)||p+pInfinity_<=0.0) {
            throw std::runtime_error(
                "StiffenedGasEOS fixed-density pressure state is invalid.");
        }
        const double e=referenceEnergy_
            +(p+gamma_*pInfinity_)/((gamma_-1.0)*rho);
        return fromDensityEnergy(rho,e);
    }

    State fromDensityTemperature(double rho, double T) const override {
        if (!std::isfinite(rho)||rho<=0.0
            ||!std::isfinite(T)||T<=0.0) {
            throw std::runtime_error(
                "StiffenedGasEOS fixed-density temperature state is invalid.");
        }
        const double e=referenceEnergy_+cv_*T+pInfinity_/rho;
        return fromDensityEnergy(rho,e);
    }

    double referenceEnergy() const override { return referenceEnergy_; }
    double gamma() const { return gamma_; }
    double cv() const { return cv_; }
    double pInfinity() const { return pInfinity_; }

private:
    double gamma_;
    double cv_;
    double pInfinity_;
    double referenceEnergy_;
};

} // namespace SF::Physics::EOS
