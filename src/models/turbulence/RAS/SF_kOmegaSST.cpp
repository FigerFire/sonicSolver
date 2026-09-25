/// @file SF_kOmegaSST.cpp
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_kOmegaSST.h"

#include "methods/numerics/structured/SF_structured.h"
#include "SF_viscous.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <algorithm>
#include <cmath>

namespace SF {
namespace Turbulence {
namespace RAS {

namespace {

struct SSTConstants {
    double betaStar = 0.09;
    double sigmaK1 = 0.85;
    double sigmaW1 = 0.5;
    double beta1 = 0.075;
    double gamma1 = 5.0 / 9.0;
    double sigmaK2 = 1.0;
    double sigmaW2 = 0.856;
    double beta2 = 0.0828;
    double gamma2 = 0.44;
    double a1 = 0.31;
};

Vector3 scalarPhysicalGradient(const Field& flow,
                               const ScalarFields& state,
                               ScalarSlot slot,
                               int i, int j, int k) {
    return CENTRAL2::gradSampled(
        flow,
        [&](int ii, int jj, int kk) {
            return state.clamped(slot, ii, jj, kk);
        },
        i, j, k);
}

double vorticityMagnitude(const Field& flow, int i, int j, int k) {
    return norm(Math::curl(Viscous::velocityGradientAt(
        flow, i, j, k, Viscous::CentralOrder::SECOND)));
}

double laminarKinematicViscosity(double rho, double laminarMu) {
    return std::max(laminarMu, 0.0) / std::max(rho, 1.0e-12);
}

double wallDistance(const Field& flow, int i, int j, int k) {
    const double d = flow.wallDistance(i, j, k);
    if (std::isfinite(d) && d > 1.0e-14) return d;
    return 1.0e20;
}

double crossDiffusionCD(const Field& flow,
                        const ScalarFields& state,
                        int i, int j, int k,
                        double rho,
                        double omega,
                        const SSTConstants& c) {
    const Vector3 gradK = scalarPhysicalGradient(
        flow, state, ScalarSlot::K, i, j, k);
    const Vector3 gradW = scalarPhysicalGradient(
        flow, state, ScalarSlot::Omega, i, j, k);
    return std::max(2.0 * rho * c.sigmaW2 * SF::dot(gradK, gradW)
                    / std::max(omega, 1.0e-30),
                    1.0e-20);
}

double blendingF1(const Field& flow,
                  const ScalarFields& state,
                  int i, int j, int k,
                  double rho,
                  double kVal,
                  double omega,
                  double laminarMu,
                  const SSTConstants& c) {
    const double d = wallDistance(flow, i, j, k);
    const double nu = laminarKinematicViscosity(rho, laminarMu);
    const double sqrtK = std::sqrt(std::max(kVal, 0.0));
    const double cd = crossDiffusionCD(flow, state, i, j, k, rho, omega, c);

    const double term1 = sqrtK / std::max(c.betaStar * omega * d, 1.0e-30);
    const double term2 = 500.0 * nu / std::max(d * d * omega, 1.0e-30);
    const double term3 = 4.0 * rho * c.sigmaW2 * kVal / std::max(cd * d * d, 1.0e-30);
    const double arg1 = std::min(std::max(term1, term2), term3);
    return std::tanh(std::pow(std::max(arg1, 0.0), 4.0));
}

double blendingF2(const Field& flow,
                  int i, int j, int k,
                  double rho,
                  double kVal,
                  double omega,
                  double laminarMu,
                  const SSTConstants& c) {
    const double d = wallDistance(flow, i, j, k);
    const double nu = laminarKinematicViscosity(rho, laminarMu);
    const double sqrtK = std::sqrt(std::max(kVal, 0.0));
    const double term1 = 2.0 * sqrtK / std::max(c.betaStar * omega * d, 1.0e-30);
    const double term2 = 500.0 * nu / std::max(d * d * omega, 1.0e-30);
    const double arg2 = std::max(term1, term2);
    return std::tanh(std::pow(std::max(arg2, 0.0), 2.0));
}

double sstScalarDiffusion(const Field& flow,
                          const ScalarFields& state,
                          ScalarSlot slot,
                          double sigma,
                          double laminarMu,
                          int i, int j, int k) {
    double result = 0.0;
    double Jcell = 1.0 / std::max(std::abs(flow.Jac(i, j, k)), 1.0e-30);

    auto gammaAt = [&](int ii, int jj, int kk) -> double {
        return std::max(laminarMu, 0.0)
             + sigma * std::max(state.EddyMu(ii, jj, kk), 0.0);
    };

    auto phiAt = [&](int ii, int jj, int kk) -> double {
        return state.slot(slot, ii, jj, kk);
    };

    auto addDirection = [&](Math::Dir d) {
        if (!Math::isDirectionActive(d)) return;

        int di, dj, dk;
        Math::dirOffset(d, di, dj, dk);

        double cfL[4], cfR[4];
        Math::faceMetrics(flow, i - di, j - dj, k - dk, d, cfL);
        Math::faceMetrics(flow, i,      j,      k,      d, cfR);

        double JgL = (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
                   / std::max(cfL[3], 1.0e-30);
        double JgR = (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
                   / std::max(cfR[3], 1.0e-30);

        double gamL = 0.5 * (gammaAt(i - di, j - dj, k - dk) + gammaAt(i, j, k));
        double gamR = 0.5 * (gammaAt(i, j, k) + gammaAt(i + di, j + dj, k + dk));

        double fluxL = gamL * JgL * (phiAt(i, j, k) - phiAt(i - di, j - dj, k - dk));
        double fluxR = gamR * JgR * (phiAt(i + di, j + dj, k + dk) - phiAt(i, j, k));

        result += (fluxR - fluxL) * Jcell;
    };

    addDirection(Math::XI);
    addDirection(Math::ETA);
    addDirection(Math::ZETA);

    return result;
}

} // namespace

void KOmegaSSTModel::initialize(const Field& flow,
                                ScalarFields& state,
                                const FDM::TurbulenceConfig& config) {
    initializeScalar(flow, state, config.scalars.kInitial, ScalarSlot::K);
    initializeScalar(flow, state, config.scalars.omegaInitial, ScalarSlot::Omega);
    applyBoundary(flow, state, config);
    refreshEddyMu(flow, state, config);
}

void KOmegaSSTModel::applyBoundary(const Field& flow,
                                   ScalarFields& state,
                                   const FDM::TurbulenceConfig& config) {
    applyScalarBC(flow, state, config.scalars.kBoundary, ScalarSlot::K);
    applyScalarBC(flow, state, config.scalars.omegaBoundary, ScalarSlot::Omega);

    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();
    for (int k = ng; k < nz + ng; ++k) {
        for (int j = ng; j < ny + ng; ++j) {
            for (int i = ng; i < nx + ng; ++i) {
                state.K(i, j, k) = clampPositive(state.K(i, j, k), config.coefficients.kFloor);
                state.Omega(i, j, k) = clampPositive(state.Omega(i, j, k), config.coefficients.omegaFloor);
            }
        }
    }
}

void KOmegaSSTModel::correct(const Field& flow,
                             ScalarFields& state,
                             const FDM::TurbulenceConfig& config,
                             double dt) {
    const auto& coeff = config.coefficients;
    const SSTConstants c;

    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();

    for (int k = ng; k < nz + ng; ++k) {
        for (int j = ng; j < ny + ng; ++j) {
            for (int i = ng; i < nx + ng; ++i) {
                if (flow.CellFlag(i, j, k) != FLUID_CELL) continue;

                double rho   = std::max(flow(i, j, k, RHO), 1.0e-12);
                double kVal  = clampPositive(state.K(i, j, k),     coeff.kFloor);
                double wVal  = clampPositive(state.Omega(i, j, k), coeff.omegaFloor);

                double sMag = strainRateMagnitude(flow, i, j, k);
                double muT  = eddyDynamicViscosity(flow, state, config, i, j, k);
                double f1 = blendingF1(
                    flow, state, i, j, k, rho, kVal, wVal,
                    config.laminarDynamicViscosity, c);
                double sigmaK = f1 * c.sigmaK1 + (1.0 - f1) * c.sigmaK2;
                double sigmaW = f1 * c.sigmaW1 + (1.0 - f1) * c.sigmaW2;
                double beta = f1 * c.beta1 + (1.0 - f1) * c.beta2;
                double gamma = f1 * c.gamma1 + (1.0 - f1) * c.gamma2;

                double production = std::min(muT * sMag * sMag,
                                             20.0 * c.betaStar * rho * kVal * wVal);

                double dk = production / rho
                          - c.betaStar * kVal * wVal;
                double domega = gamma * production / std::max(muT, 1.0e-30)
                              - beta * wVal * wVal;

                const Vector3 gradK = scalarPhysicalGradient(
                    flow, state, ScalarSlot::K, i, j, k);
                const Vector3 gradW = scalarPhysicalGradient(
                    flow, state, ScalarSlot::Omega, i, j, k);
                double crossDiff = 2.0 * (1.0 - f1) * rho * c.sigmaW2
                                 * SF::dot(gradK, gradW)
                                 / std::max(wVal, coeff.omegaFloor);

                // ---- SST diffusion  ∇·((μ + σ μ_t) ∇φ) / ρ ----
                double diffK = sstScalarDiffusion(
                    flow, state, ScalarSlot::K, sigmaK,
                    config.laminarDynamicViscosity, i, j, k);
                double diffW = sstScalarDiffusion(
                    flow, state, ScalarSlot::Omega, sigmaW,
                    config.laminarDynamicViscosity, i, j, k);

                dk     += diffK / rho;
                domega += (diffW + crossDiff) / rho;

                state.K(i, j, k)     = clampPositive(kVal + dt * dk,     coeff.kFloor);
                state.Omega(i, j, k) = clampPositive(wVal + dt * domega, coeff.omegaFloor);
            }
        }
    }

    applyBoundary(flow, state, config);
    refreshEddyMu(flow, state, config);
}

double KOmegaSSTModel::eddyDynamicViscosity(const Field& flow,
                                            const ScalarFields& state,
                                            const FDM::TurbulenceConfig& config,
                                            int i, int j, int k) const {
    const SSTConstants c;
    double rho = std::max(flow(i, j, k, RHO), 1.0e-12);
    double kVal = clampPositive(state.clamped(ScalarSlot::K, i, j, k), config.coefficients.kFloor);
    double omegaVal = clampPositive(state.clamped(ScalarSlot::Omega, i, j, k), config.coefficients.omegaFloor);
    double f2 = blendingF2(
        flow, i, j, k, rho, kVal, omegaVal,
        config.laminarDynamicViscosity, c);
    double omegaMag = vorticityMagnitude(flow, i, j, k);
    return rho * c.a1 * kVal
         / std::max(c.a1 * omegaVal, omegaMag * f2);
}

} // namespace RAS
} // namespace Turbulence
} // namespace SF
