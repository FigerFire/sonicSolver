/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_WENO.h
/// @brief WENO界面重构和Flux模块通量方法的统一装配入口。
///
/// 本文件属于对流离散层: WENO负责把单元中心状态重构到界面左右状态。
/// 真正的面通量算法(Roe/Lax-Wendroff/Lax-Friedrichs/Steger-Warming/ILW)
/// 全部通过数值通量核调用。

#include "core/config/SF_config.h"
#include "core/field/SF_field.h"
#include "methods/numerics/Flux/SF_face.h"
#include "methods/numerics/Flux/SF_lowOrder.h"
#include "methods/numerics/Flux/SF_lw.h"
#include "methods/numerics/Flux/SF_lxF.h"
#include "methods/numerics/structured/SF_structured.h"
#include "methods/numerics/Flux/SF_roe.h"
#include "methods/numerics/Flux/SF_rusanovEOS.h"
#include "methods/numerics/Flux/SF_sw.h"
#include "methods/numerics/convection/SF_riemann.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace Flux {

inline void storeLowOrderFaceFlux(
        FluxField& fluxField,
        Field& field,
        int i,
        int j,
        int k,
        Math::Dir direction,
        LowOrderFlux::Method method,
        double gamma,
        const Physics::EquationSet::Model& thermodynamics) {
    if (method != LowOrderFlux::Method::FirstOrderRusanov) {
        LowOrderFlux::storeFaceFlux(
            fluxField, field, i, j, k, direction, method, gamma);
        return;
    }

    int di = 0, dj = 0, dk = 0;
    Math::dirOffset(direction, di, dj, dk);
    const FaceGeometry geometry = makeFaceGeometry(
        field, i, j, k, direction);
    double qLeft[5], qRight[5], faceFlux[5];
    loadConservative(field, i, j, k, qLeft);
    loadConservative(field, i + di, j + dj, k + dk, qRight);
    Numerics::RusanovEOS::faceFlux(
        thermodynamics, qLeft, qRight,
        {geometry.normal[0], geometry.normal[1], geometry.normal[2]},
        faceFlux, 5);
    storePhysicalFlux(
        fluxField, field, i, j, k, direction, geometry, faceFlux);
}

template <int NStencil>
inline constexpr int scalarWENOStencilWidth() {
    if constexpr (NStencil == 4) {
        return 4;
    } else if constexpr (NStencil == 6) {
        return 5;
    } else if constexpr (NStencil == 8) {
        return 7;
    } else {
        static_assert(NStencil == 4 || NStencil == 6 || NStencil == 8,
                      "Unsupported WENO stencil size");
        return 0;
    }
}

template <int NStencil>
inline void fillScalarWENOStencils(const double characteristicStencil[NStencil][5],
                                   int component,
                                   double left[8],
                                   double right[8]) {
    constexpr int width = scalarWENOStencilWidth<NStencil>();
    for (int s = 0; s < width; ++s) {
        left[s] = characteristicStencil[s][component];
        right[s] = characteristicStencil[NStencil - 1 - s][component];
        if (!std::isfinite(left[s]) || !std::isfinite(right[s])) {
            std::cerr << "[SF FATAL] non-finite characteristic WENO stencil "
                      << "value, NStencil=" << NStencil
                      << ", component=" << component
                      << ", scalarIndex=" << s
                      << ", left=" << left[s]
                      << ", right=" << right[s] << std::endl;
            std::exit(1);
        }
    }
}

/// @brief 特征空间WENO重构并使用通量分裂计算一个面通量。
/// @tparam NStencil 模板点数，WENO3=4、WENO5=6、WENO7=8。
/// @tparam WenoKernel 标量WENO核函数类型。
/// @tparam RiemannKernel 通量分裂函数类型。
/// @param qStencil 守恒变量模板数组。
/// @param metrics 半网格面度规。
/// @param flux 输出已乘面积的离散面通量。
/// @param wenoCore 标量WENO核函数。
/// @param riemann 正/负通量分裂函数。
template <int NStencil, typename WenoKernel, typename RiemannKernel>
inline void characteristicWENO(const double qStencil[NStencil][5],
                               const double metrics[4],
                               double* flux,
                               WenoKernel&& wenoCore,
                               RiemannKernel&& riemann,
                               double gamma = 1.4) {
    const double sx = metrics[0];
    const double sy = metrics[1];
    const double sz = metrics[2];
    const double smag = std::sqrt(sx * sx + sy * sy + sz * sz);
    if (!std::isfinite(smag) || smag <= 1.0e-300) {
        std::cerr << "[SF FATAL] characteristic WENO received a degenerate "
                  << "face cofactor, area=" << smag << std::endl;
        std::exit(1);
    }
    const double normal[3] = {
        sx / smag,
        sy / smag,
        sz / smag
    };
    const double area = smag;

    constexpr int centerLeft = NStencil / 2 - 1;
    double qAverage[5];
    for (int v = 0; v < 5; ++v) {
        qAverage[v] = 0.5 * (qStencil[centerLeft][v] + qStencil[centerLeft + 1][v]);
    }

    double leftEigen[25], rightEigen[25];
    Riemann::buildCharacteristicMatrix(
        qAverage, normal, gamma, leftEigen, rightEigen);

    double characteristicStencil[NStencil][5];
    for (int s = 0; s < NStencil; ++s) {
        Riemann::gemv5(leftEigen, qStencil[s], characteristicStencil[s]);
    }

    double wLeft[5], wRight[5];
    for (int v = 0; v < 5; ++v) {
        double stencilLeft[8] = {};
        double stencilRight[8] = {};
        fillScalarWENOStencils<NStencil>(characteristicStencil, v,
                                         stencilLeft, stencilRight);
        wLeft[v] = wenoCore(stencilLeft);
        wRight[v] = wenoCore(stencilRight);
        if (!std::isfinite(wLeft[v]) || !std::isfinite(wRight[v])) {
            std::cerr << "[SF FATAL] scalar WENO reconstruction produced "
                      << "non-finite characteristic value, NStencil="
                      << NStencil << ", component=" << v
                      << ", left=" << wLeft[v]
                      << ", right=" << wRight[v] << std::endl;
            std::exit(1);
        }
    }

    double qLeft[5], qRight[5];
    Riemann::gemv5(rightEigen, wLeft, qLeft);
    Riemann::gemv5(rightEigen, wRight, qRight);

    double fPos[5], unused[5], fNeg[5];
    riemann(qLeft, normal, fPos, unused);
    riemann(qRight, normal, unused, fNeg);
    for (int v = 0; v < 5; ++v) flux[v] = (fPos[v] + fNeg[v]) * area;
}

/// @brief 特征空间WENO重构左右界面状态，再交给直接面通量函数。
/// @tparam NStencil 模板点数，WENO3=4、WENO5=6、WENO7=8。
/// @tparam WenoKernel 标量WENO核函数类型。
/// @tparam FaceFlux 直接面通量函数类型。
/// @param qStencil 守恒变量模板数组。
/// @param metrics 半网格面度规。
/// @param flux 输出已乘面积的离散面通量。
/// @param wenoCore 标量WENO核函数。
/// @param faceFlux 直接面通量函数: `(qL, qR, normal, fluxNoArea)`。
template <int NStencil, typename WenoKernel, typename FaceFlux>
inline void characteristicWENOFaceFlux(const double qStencil[NStencil][5],
                                       const double metrics[4],
                                       double* flux,
                                       WenoKernel&& wenoCore,
                                       FaceFlux&& faceFlux,
                                       double gamma = 1.4) {
    const double sx = metrics[0];
    const double sy = metrics[1];
    const double sz = metrics[2];
    const double smag = std::sqrt(sx * sx + sy * sy + sz * sz);
    if (!std::isfinite(smag) || smag <= 1.0e-300) {
        std::cerr << "[SF FATAL] characteristic WENO face flux received a "
                  << "degenerate face cofactor, area=" << smag << std::endl;
        std::exit(1);
    }
    const double normal[3] = {
        sx / smag,
        sy / smag,
        sz / smag
    };
    const double area = smag;

    constexpr int centerLeft = NStencil / 2 - 1;
    double qAverage[5];
    for (int v = 0; v < 5; ++v) {
        qAverage[v] = 0.5 * (qStencil[centerLeft][v] + qStencil[centerLeft + 1][v]);
    }

    double leftEigen[25], rightEigen[25];
    Riemann::buildCharacteristicMatrix(
        qAverage, normal, gamma, leftEigen, rightEigen);

    double characteristicStencil[NStencil][5];
    for (int s = 0; s < NStencil; ++s) {
        Riemann::gemv5(leftEigen, qStencil[s], characteristicStencil[s]);
    }

    double wLeft[5], wRight[5];
    for (int v = 0; v < 5; ++v) {
        double stencilLeft[8] = {};
        double stencilRight[8] = {};
        fillScalarWENOStencils<NStencil>(characteristicStencil, v,
                                         stencilLeft, stencilRight);
        wLeft[v] = wenoCore(stencilLeft);
        wRight[v] = wenoCore(stencilRight);
        if (!std::isfinite(wLeft[v]) || !std::isfinite(wRight[v])) {
            std::cerr << "[SF FATAL] scalar WENO reconstruction produced "
                      << "non-finite characteristic value, NStencil="
                      << NStencil << ", component=" << v
                      << ", left=" << wLeft[v]
                      << ", right=" << wRight[v] << std::endl;
            std::exit(1);
        }
    }

    double qLeft[5], qRight[5], faceFluxNoArea[5];
    Riemann::gemv5(rightEigen, wLeft, qLeft);
    Riemann::gemv5(rightEigen, wRight, qRight);
    faceFlux(qLeft, qRight, normal, faceFluxNoArea);
    for (int v = 0; v < 5; ++v) flux[v] = faceFluxNoArea[v] * area;
}

/// @brief 计算一个方向的WENO面通量。
/// @tparam NStencil 模板点数，WENO3=4、WENO5=6、WENO7=8。
/// @tparam Kernel 面通量装配核函数类型。
/// @param field 结构网格场，ConvectiveFlux会被写入。
/// @param d 计算坐标方向。
/// @param kernel 负责从模板计算面通量的函数。
/// @param fallbackMethod WENO模板碰到IBM非流体点时采用的低阶回退方法。
/// @param ibmBoundary IBM近壁模板处理方式；默认降阶，ILW时先闭合模板再WENO。
/// @param requestedILWOrder 用户请求的ILW Taylor阶数；0表示按模块默认能力处理。
template <int NStencil, typename Kernel>
inline void computeDirectionalFlux(Field& field, FluxField& fluxField,
                                   Math::Dir d,
                                   Kernel&& kernel,
                                   LowOrderFlux::Method fallbackMethod =
                                       LowOrderFlux::Method::FirstOrderRoe,
                                   FDM::IBMBoundaryScheme ibmBoundary =
                                       FDM::IBMBoundaryScheme::LowOrder,
                                   int requestedILWOrder = 0,
                                   double gamma = 1.4,
                                   const Physics::EquationSet::Model* thermodynamics = nullptr) {
    const int* offsets = nullptr;
    if constexpr (NStencil == 4) offsets = Math::WENO3_OFFSETS;
    else if constexpr (NStencil == 6) offsets = Math::WENO5_OFFSETS;
    else if constexpr (NStencil == 8) offsets = Math::WENO7_OFFSETS;
    else static_assert(NStencil == 0, "Unsupported WENO stencil size");

    Math::forFaces(field, d, [&](int i, int j, int k) {
        if (!Math::shouldCalculateFaceFlux(field, d, i, j, k)) return;
        int di = 0, dj = 0, dk = 0;
        Math::dirOffset(d, di, dj, dk);

        if (!faceTouchesPhysicalFluid(field, i, j, k, d)) {
            double zero[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
            Math::storeFlux(fluxField, field, i, j, k, d, zero);
            return;
        }

        if (faceIsIBMInterface(field, i, j, k, d)
            && ibmBoundary != FDM::IBMBoundaryScheme::ILW) {
            if (!thermodynamics) {
                throw std::runtime_error(
                    "Reconstructed convection requires an active EquationSet.");
            }
            storeLowOrderFaceFlux(
                fluxField, field, i, j, k, d, fallbackMethod, gamma,
                *thermodynamics);
            return;
        }

        if (SF::IBM::isIbmNonFluidCell(field, i, j, k)
            && SF::IBM::isIbmNonFluidCell(field, i + di, j + dj, k + dk)) {
            double zero[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
            Math::storeFlux(fluxField, field, i, j, k, d, zero);
            return;
        }

        if (stencilTouchesIBM(field, i, j, k, d, offsets, NStencil)) {
            if (ibmBoundary == FDM::IBMBoundaryScheme::ILW) {
                double stencil[NStencil][5];
                double metrics[4];
                double flux[5];
                const int ilwTaylorOrder =
                    (requestedILWOrder > 0)
                        ? FDM::ilwTaylorOrderFromAccuracy(requestedILWOrder)
                        : ((NStencil == 4) ? 2 : 4);
                if (stencilClosedByIBMGhosts(field, i, j, k, d,
                                             offsets, NStencil)) {
                    Math::extractStencil(field, i, j, k, d, offsets,
                                         NStencil, &stencil[0][0]);
                    Math::faceMetrics(field, i, j, k, d, metrics);
                    kernel(stencil, metrics, flux, i, j, k, d);
                    bool fluxValid = true;
                    for (double v : flux) {
                        if (!std::isfinite(v)) { fluxValid = false; break; }
                    }
                    if (fluxValid) {
                        Math::storeFlux(fluxField, field, i, j, k, d, flux);
                        return;
                    }
                    std::cerr
                        << "[SF FATAL] ILW closed WENO stencil produced "
                        << "invalid flux at face (" << i << "," << j
                        << "," << k << "). LowOrder fallback is disabled; "
                        << "set [numerics].ILW=0 and select a low-order IBM "
                        << "scheme explicitly if that behavior is desired."
                        << std::endl;
                    std::exit(1);
                }
                std::cerr
                    << "[SF FATAL] ILW closed WENO stencil failed at face ("
                    << i << "," << j << "," << k
                    << "), requestedTaylorOrder=" << ilwTaylorOrder
                    << ". Stencil reached SOLID_CELL. LowOrder fallback is disabled; "
                    << "set [numerics].ILW=0 "
                    << "and select a low-order IBM scheme explicitly if that "
                    << "behavior is desired."
                    << std::endl;
                std::exit(1);
            }
            if (!thermodynamics) {
                throw std::runtime_error(
                    "Reconstructed convection requires an active EquationSet.");
            }
            storeLowOrderFaceFlux(
                fluxField, field, i, j, k, d, fallbackMethod, gamma,
                *thermodynamics);
            return;
        }

        double stencil[NStencil][5];
        double metrics[4];
        double flux[5];
        Math::extractStencil(field, i, j, k, d, offsets, NStencil, &stencil[0][0]);
        Math::faceMetrics(field, i, j, k, d, metrics);
        kernel(stencil, metrics, flux, i, j, k, d);
        Math::storeFlux(fluxField, field, i, j, k, d, flux);
    });
}

/// @brief 计算三个方向的WENO面通量。
/// @tparam NStencil 模板点数，WENO3=4、WENO5=6、WENO7=8。
/// @tparam Kernel 面通量装配核函数类型。
/// @param field 结构网格场，ConvectiveFlux会被写入。
/// @param kernel 负责从模板计算面通量的函数。
/// @param fallbackMethod WENO模板碰到IBM非流体点时采用的低阶回退方法。
/// @param ibmBoundary IBM近壁模板处理方式。
/// @param requestedILWOrder 用户请求的ILW Taylor阶数。
template <int NStencil = 6, typename Kernel>
inline void computeAllFluxes(Field& field, FluxField& fluxField,
                             Kernel&& kernel,
                             LowOrderFlux::Method fallbackMethod =
                                 LowOrderFlux::Method::FirstOrderRoe,
                             FDM::IBMBoundaryScheme ibmBoundary =
                                 FDM::IBMBoundaryScheme::LowOrder,
                             int requestedILWOrder = 0,
                             double gamma = 1.4,
                             const Physics::EquationSet::Model* thermodynamics = nullptr) {
    computeDirectionalFlux<NStencil>(
        field, fluxField, Math::XI, kernel, fallbackMethod, ibmBoundary,
        requestedILWOrder, gamma, thermodynamics);
    computeDirectionalFlux<NStencil>(
        field, fluxField, Math::ETA, kernel, fallbackMethod, ibmBoundary,
        requestedILWOrder, gamma, thermodynamics);
    computeDirectionalFlux<NStencil>(
        field, fluxField, Math::ZETA, kernel, fallbackMethod, ibmBoundary,
        requestedILWOrder, gamma, thermodynamics);
}

} // namespace Flux
} // namespace SF
