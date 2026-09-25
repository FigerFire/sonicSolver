/// @file SF_lowOrder.cpp
/// @brief 可压缩通量或 Riemann 数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_lowOrder.h"

#include "SF_lxF.h"
#include "SF_roe.h"
#include "SF_sw.h"

#include <stdexcept>

namespace SF {
namespace Flux {
namespace LowOrderFlux {

Method fallbackFor(FDM::FluxSplitter method) {
    switch (method) {
        case FDM::FluxSplitter::LaxFriedrichs:
            return Method::FirstOrderLaxFriedrichs;
        case FDM::FluxSplitter::StegerWarming:
            return Method::FirstOrderStegerWarming;
        case FDM::FluxSplitter::Rusanov:
            return Method::FirstOrderRusanov;
        case FDM::FluxSplitter::Roe:
        case FDM::FluxSplitter::LaxWendroff:
            return Method::FirstOrderRoe;
    }
    throw std::invalid_argument("Unknown low-order numerical flux method.");
}

void flux(const double qL[5], const double qR[5],
          const double normal[3],
          Method method,
          double fluxOut[5],
          double gamma) {
    if (method == Method::FirstOrderRusanov) {
        throw std::invalid_argument(
            "FirstOrderRusanov requires the active FluidStateModel; use the "
            "FluidStateModel-aware storeFaceFlux overload.");
    }
    if (method == Method::FirstOrderRoe) {
        RoeFlux::flux(qL, qR, normal, fluxOut, gamma);
        return;
    }

    double fPos[5], fNeg[5], unused[5];
    if (method == Method::FirstOrderLaxFriedrichs) {
        LaxFriedrichsFlux::split(qL, normal, fPos, unused, gamma);
        LaxFriedrichsFlux::split(qR, normal, unused, fNeg, gamma);
    } else {
        StegerWarmingFlux::split(qL, normal, fPos, unused, gamma);
        StegerWarmingFlux::split(qR, normal, unused, fNeg, gamma);
    }

    for (int v = 0; v < 5; ++v) fluxOut[v] = fPos[v] + fNeg[v];
}

void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k,
                   Math::Dir d, Method method,
                   double gamma) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    if (SF::IBM::isIbmNonFluidCell(field, i, j, k)
        && SF::IBM::isIbmNonFluidCell(field, i + di, j + dj, k + dk)) {
        double zero[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        Math::storeFlux(fluxField, field, i, j, k, d, zero);
        return;
    }

    const FaceGeometry geom = makeFaceGeometry(field, i, j, k, d);
    double qL[5], qR[5], faceFlux[5];
    loadConservative(field, i, j, k, qL);
    loadConservative(field, i + di, j + dj, k + dk, qR);
    flux(qL, qR, geom.normal, method, faceFlux, gamma);
    storePhysicalFlux(fluxField, field, i, j, k, d, geom, faceFlux);
}

} // namespace LowOrderFlux
} // namespace Flux
} // namespace SF
