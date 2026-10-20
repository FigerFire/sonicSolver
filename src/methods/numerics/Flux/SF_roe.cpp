/// @file SF_roe.cpp
/// @brief 可压缩通量或 Riemann 数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_roe.h"

#include "SF_riemann.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace Flux {
namespace RoeFlux {

void flux(const double qL[5], const double qR[5],
          const double normal[3], double fluxOut[5], double gamma) {
    const double gm1 = gamma - 1.0;

    const double pL = Numerics::requirePhysicalState("RoeFlux(left)",
                                                     qL[0], qL[1], qL[2],
                                                     qL[3], qL[4], gamma);
    const double pR = Numerics::requirePhysicalState("RoeFlux(right)",
                                                     qR[0], qR[1], qR[2],
                                                     qR[3], qR[4], gamma);

    const double invRL = 1.0 / qL[0];
    const double uL = qL[1] * invRL;
    const double vL = qL[2] * invRL;
    const double wL = qL[3] * invRL;
    const double hL = (qL[4] + pL) * invRL;

    const double invRR = 1.0 / qR[0];
    const double uR = qR[1] * invRR;
    const double vR = qR[2] * invRR;
    const double wR = qR[3] * invRR;
    const double hR = (qR[4] + pR) * invRR;

    const double sqrtRL = std::sqrt(qL[0]);
    const double sqrtRR = std::sqrt(qR[0]);
    const double sumSqrt = sqrtRL + sqrtRR;

    const double uRoe = (sqrtRL * uL + sqrtRR * uR) / sumSqrt;
    const double vRoe = (sqrtRL * vL + sqrtRR * vR) / sumSqrt;
    const double wRoe = (sqrtRL * wL + sqrtRR * wR) / sumSqrt;
    const double hRoe = (sqrtRL * hL + sqrtRR * hR) / sumSqrt;
    const double q2Roe = uRoe * uRoe + vRoe * vRoe + wRoe * wRoe;
    const double a2Roe = gm1 * (hRoe - 0.5 * q2Roe);

    if (!std::isfinite(a2Roe) || a2Roe <= 0.0) {
        std::cerr << "[SF FATAL] RoeFlux: invalid Roe averaged sound speed"
                  << " a2=" << a2Roe
                  << " h=" << hRoe
                  << " q2=" << q2Roe
                  << std::endl;
        std::exit(1);
    }

    const double aRoe = std::sqrt(a2Roe);
    const double unRoe = uRoe * normal[0] + vRoe * normal[1] + wRoe * normal[2];
    const double lambda[5] = {unRoe, unRoe, unRoe, unRoe + aRoe, unRoe - aRoe};

    double absLambda[5];
    const double entropyEps = 0.1 * aRoe;
    for (int m = 0; m < 5; ++m) {
        const double value = std::abs(lambda[m]);
        if (value < 2.0 * entropyEps) {
            absLambda[m] = 0.5 * (value * value / (2.0 * entropyEps) + entropyEps);
        } else {
            absLambda[m] = value;
        }
    }

    double dQ[5];
    for (int v = 0; v < 5; ++v) dQ[v] = qR[v] - qL[v];

    const double rhoRoe = sqrtRL * sqrtRR;
    double qRoe[5];
    qRoe[0] = rhoRoe;
    qRoe[1] = rhoRoe * uRoe;
    qRoe[2] = rhoRoe * vRoe;
    qRoe[3] = rhoRoe * wRoe;
    qRoe[4] = rhoRoe * (0.5 * q2Roe + a2Roe / gm1);

    double leftEigen[25], rightEigen[25];
    Riemann::buildCharacteristicMatrix(
        qRoe, normal, gamma, leftEigen, rightEigen);

    double alpha[5];
    Riemann::gemv5(leftEigen, dQ, alpha);

    double dissipation[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (int wave = 0; wave < 5; ++wave) {
        const double strength = absLambda[wave] * alpha[wave];
        for (int v = 0; v < 5; ++v) {
            dissipation[v] += rightEigen[v * 5 + wave] * strength;
        }
    }

    double fL[5], fR[5];
    physicalEulerFlux(qL, normal, gamma, fL);
    physicalEulerFlux(qR, normal, gamma, fR);

    for (int v = 0; v < 5; ++v) {
        fluxOut[v] = 0.5 * (fL[v] + fR[v] - dissipation[v]);
    }
}

void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k, Math::Dir d,
                   double gamma) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);

    const FaceGeometry geom = makeFaceGeometry(field, i, j, k, d);
    double qL[5], qR[5], faceFlux[5];
    loadConservative(field, i, j, k, qL);
    loadConservative(field, i + di, j + dj, k + dk, qR);
    flux(qL, qR, geom.normal, faceFlux, gamma);
    storePhysicalFlux(fluxField, field, i, j, k, d, geom, faceFlux);
}

void computeAllFluxes(Field& field, FluxField& fluxField, Residual& residual,
                      double gamma) {
    auto computeDirection = [&](Math::Dir d) {
        Math::forFaces(field, d, [&](int i, int j, int k) {
            if (!Math::shouldCalculateFaceFlux(field, d, i, j, k)) return;
            if (SF::IBM::isIbmNonFluidCell(field, i, j, k)) {
                int di, dj, dk;
                Math::dirOffset(d, di, dj, dk);
                if (SF::IBM::isIbmNonFluidCell(field, i + di, j + dj, k + dk)) {
                    double zero[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
                    Math::storeFlux(fluxField, field, i, j, k, d, zero);
                    return;
                }
            }
            storeFaceFlux(fluxField, field, i, j, k, d, gamma);
        });
    };

    computeDirection(Math::XI);
    computeDirection(Math::ETA);
    computeDirection(Math::ZETA);
    Math::assembleConvectiveFluxResidual(field, fluxField, residual);
}

} // namespace RoeFlux
} // namespace Flux
} // namespace SF
