/// @file SF_kEpsilon.cpp
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_kEpsilon.h"

#include <algorithm>

namespace SF {
namespace Turbulence {
namespace RAS {

void KEpsilonModel::initialize(const Field& flow,
                               ScalarFields& state,
                               const FDM::TurbulenceConfig& config) {
    initializeScalar(flow, state, config.scalars.kInitial, ScalarSlot::K);
    initializeScalar(flow, state, config.scalars.epsilonInitial, ScalarSlot::Epsilon);
    applyBoundary(flow, state, config);
    refreshEddyMu(flow, state, config);
}

void KEpsilonModel::applyBoundary(const Field& flow,
                                  ScalarFields& state,
                                  const FDM::TurbulenceConfig& config) {
    applyScalarBC(flow, state, config.scalars.kBoundary, ScalarSlot::K);
    applyScalarBC(flow, state, config.scalars.epsilonBoundary, ScalarSlot::Epsilon);

    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();
    for (int k = ng; k < nz + ng; ++k) {
        for (int j = ng; j < ny + ng; ++j) {
            for (int i = ng; i < nx + ng; ++i) {
                state.K(i, j, k) = clampPositive(state.K(i, j, k), config.coefficients.kFloor);
                state.Epsilon(i, j, k) = clampPositive(state.Epsilon(i, j, k), config.coefficients.epsilonFloor);
            }
        }
    }
}

void KEpsilonModel::correct(const Field& flow,
                            ScalarFields& state,
                            const FDM::TurbulenceConfig& config,
                            double dt) {
    const auto& coeff = config.coefficients;
    const double sigmaK    = std::max(coeff.sigmaK, 1.0e-12);
    const double sigmaEps  = std::max(coeff.sigmaEpsilon, 1.0e-12);
    const double invSigmaK = 1.0 / sigmaK;
    const double invSigmaEps = 1.0 / sigmaEps;

    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();

    for (int k = ng; k < nz + ng; ++k) {
        for (int j = ng; j < ny + ng; ++j) {
            for (int i = ng; i < nx + ng; ++i) {
                if (flow.CellFlag(i, j, k) != FLUID_CELL) continue;

                double rho = std::max(flow(i, j, k, RHO), 1.0e-12);
                double kVal   = clampPositive(state.K(i, j, k),       coeff.kFloor);
                double epsVal = clampPositive(state.Epsilon(i, j, k), coeff.epsilonFloor);

                // Eddy viscosity: standard k-epsilon  (Cmu * rho * k^2 / epsilon)
                double muT = coeff.cMu * rho * kVal * kVal / epsVal;

                double sMag = strainRateMagnitude(flow, i, j, k);
                double production = std::min(muT * sMag * sMag,
                                             50.0 * rho * epsVal);

                // ---- source terms (σ_k / σ_ε are diffusion coefficients,
                //      they must NOT divide the dissipation terms) ----
                double dk   = production / rho
                            - epsVal;
                double deps = coeff.c1 * epsVal * production
                              / std::max(rho * kVal, coeff.kFloor)
                            - coeff.c2 * epsVal * epsVal
                              / std::max(kVal, coeff.kFloor);

                // ---- turbulent diffusion  ∇·((μ_t / σ) ∇φ) / ρ ----
                double diffK   = scalarTurbulentDiffusion(
                    flow, state, ScalarSlot::K, invSigmaK, i, j, k);
                double diffEps = scalarTurbulentDiffusion(
                    flow, state, ScalarSlot::Epsilon, invSigmaEps, i, j, k);

                dk   += diffK   / rho;
                deps += diffEps / rho;

                state.K(i, j, k)       = clampPositive(kVal   + dt * dk,   coeff.kFloor);
                state.Epsilon(i, j, k) = clampPositive(epsVal + dt * deps, coeff.epsilonFloor);
            }
        }
    }

    applyBoundary(flow, state, config);
    refreshEddyMu(flow, state, config);
}

double KEpsilonModel::eddyDynamicViscosity(const Field& flow,
                                           const ScalarFields& state,
                                           const FDM::TurbulenceConfig& config,
                                           int i, int j, int k) const {
    double rho = std::max(flow(i, j, k, RHO), 1.0e-12);
    double kVal = clampPositive(state.clamped(ScalarSlot::K, i, j, k), config.coefficients.kFloor);
    double epsVal = clampPositive(state.clamped(ScalarSlot::Epsilon, i, j, k), config.coefficients.epsilonFloor);
    return config.coefficients.cMu * rho * kVal * kVal / epsVal;
}

} // namespace RAS
} // namespace Turbulence
} // namespace SF
