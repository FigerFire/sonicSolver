/// @file SF_lw.cpp
/// @brief 可压缩通量或 Riemann 数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_lw.h"

#include <cstdlib>
#include <iostream>

namespace SF {
namespace Flux {
namespace LaxWendroffFlux {

void flux(const double qL[5], const double qR[5],
          const double normal[3],
          double dt, double spacing,
          double fluxOut[5], double gamma) {
    if (!(dt > 0.0) || !(spacing > 0.0)) {
        std::cerr << "[SF FATAL] LaxWendroffFlux: dt and spacing must be positive"
                  << " dt=" << dt
                  << " spacing=" << spacing
                  << std::endl;
        std::exit(1);
    }

    double fL[5], fR[5];
    physicalEulerFlux(qL, normal, gamma, fL);
    physicalEulerFlux(qR, normal, gamma, fR);

    const double lambda = dt / spacing;
    double qHalf[5];
    for (int v = 0; v < 5; ++v) {
        qHalf[v] = 0.5 * (qL[v] + qR[v]) - 0.5 * lambda * (fR[v] - fL[v]);
    }

    physicalEulerFlux(qHalf, normal, gamma, fluxOut);
}

void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k, Math::Dir d,
                   double dt, double gamma) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);

    const FaceGeometry geom = makeFaceGeometry(field, i, j, k, d);
    double qL[5], qR[5], faceFlux[5];
    loadConservative(field, i, j, k, qL);
    loadConservative(field, i + di, j + dj, k + dk, qR);
    flux(qL, qR, geom.normal, dt,
         faceSpacing(field, i, j, k, d, geom.normal), faceFlux, gamma);
    storePhysicalFlux(fluxField, field, i, j, k, d, geom, faceFlux);
}

void computeAllFluxes(Field& field, FluxField& fluxField, Residual& residual,
                      double dt, double gamma) {
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
            storeFaceFlux(fluxField, field, i, j, k, d, dt, gamma);
        });
    };

    computeDirection(Math::XI);
    computeDirection(Math::ETA);
    computeDirection(Math::ZETA);
    Math::assembleConvectiveFluxResidual(field, fluxField, residual);
}

} // namespace LaxWendroffFlux
} // namespace Flux
} // namespace SF
