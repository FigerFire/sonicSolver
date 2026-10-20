/// @file SF_lxF.cpp
/// @brief 可压缩通量或 Riemann 数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_lxF.h"

#include "SF_riemann.h"

namespace SF {
namespace Flux {
namespace LaxFriedrichsFlux {

void split(const double q[5], const double normal[3],
           double fPos[5], double fNeg[5], double gamma) {
    Riemann::splitLaxFriedrichs(q, normal, gamma, fPos, fNeg);
}

void flux(const double qL[5], const double qR[5],
          const double normal[3], double fluxOut[5], double gamma) {
    double fPos[5], fNeg[5], unused[5];
    split(qL, normal, fPos, unused, gamma);
    split(qR, normal, unused, fNeg, gamma);
    for (int v = 0; v < 5; ++v) fluxOut[v] = fPos[v] + fNeg[v];
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
            storeFaceFlux(fluxField, field, i, j, k, d, gamma);
        });
    };

    computeDirection(Math::XI);
    computeDirection(Math::ETA);
    computeDirection(Math::ZETA);
    Math::assembleConvectiveFluxResidual(field, fluxField, residual);
}

} // namespace LaxFriedrichsFlux
} // namespace Flux
} // namespace SF
