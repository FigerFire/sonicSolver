/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_div.h
/// @brief 对流通量散度离散算子 (convective flux divergence).
///
/// 每个对流格式是一个独立的命名空间，内部提供统一的 div() 接口。
/// 调用方式:
/// @code
///   WENO5::div(field, FDM::FluxSplitter::Roe, dt); // WENO5 + Roe面通量
/// @endcode
///
/// div 计算 ∇·F_c(Q) 时先将界面数值通量写入FluxField，
/// 再统一装配到结构方向残差缓存供 Math::spatialResidual 取散度。

#include "SF_field.h"
#include "SF_config.h"
#include "methods/numerics/structured/SF_structured.h"
#include "solver/discretization/convection/SF_WENO.h"

// 对流底层模块
#include "methods/numerics/convection/SF_WENO3.h"
#include "methods/numerics/convection/SF_WENO5.h"
#include "methods/numerics/convection/SF_TENO5.h"
#include "methods/numerics/convection/SF_WENO7.h"
#include "methods/numerics/convection/SF_riemann.h"

#include <cstdlib>
#include <stdexcept>

namespace SF {

namespace DivDetail {

/// @brief 通用WENO对流装配：WENO负责界面重构，flux负责面通量。
/// @param field 结构网格场。
/// @param fluxMethod 用户选择的面通量方法。
/// @param dt 当前时间步；仅Lax-Wendroff使用。
/// @param ibmBoundary IBM近壁WENO模板处理方式。
/// @param requestedILWOrder 用户请求的ILW Taylor阶数；0表示关闭或默认。
/// @param schemeName 日志/报错用格式名。
/// @param requiredGhost 当前WENO阶数需要的ghost层数。
/// @param wenoCore 标量WENO核函数。
template <int NStencil, typename WenoCore>
inline void wenoDiv(Field& field, FluxField& fluxField, Residual& residual,
                    FDM::FluxSplitter fluxMethod,
                    double dt,
                    FDM::IBMBoundaryScheme ibmBoundary,
                    int requestedILWOrder,
                    const char* schemeName,
                    int requiredGhost,
                    double gamma,
                    const Physics::EquationSet::Model& thermodynamics,
                    WenoCore&& wenoCore) {
    if (!Math::checkGhostDepth(field, requiredGhost, schemeName)) std::exit(1);

    const auto fallback = Flux::LowOrderFlux::fallbackFor(fluxMethod);

    Flux::computeAllFluxes<NStencil>(field, fluxField,
        [&](const auto& s, const auto& m, auto* f,
            int i, int j, int k, Math::Dir d) {
            Flux::characteristicWENOFaceFlux<NStencil>(
                s, m, f, wenoCore,
                [&](const double qLeft[5], const double qRight[5],
                    const double normal[3], double out[5]) {
                    switch (fluxMethod) {
                    case FDM::FluxSplitter::StegerWarming:
                        Flux::StegerWarmingFlux::flux(
                            qLeft, qRight, normal, out, gamma);
                        return;
                    case FDM::FluxSplitter::Rusanov:
                        Numerics::RusanovEOS::faceFlux(
                            thermodynamics, qLeft, qRight,
                            {normal[0], normal[1], normal[2]}, out, 5);
                        return;
                    case FDM::FluxSplitter::LaxFriedrichs:
                        Flux::LaxFriedrichsFlux::flux(
                            qLeft, qRight, normal, out, gamma);
                        return;
                    case FDM::FluxSplitter::Roe:
                        Flux::RoeFlux::flux(
                            qLeft, qRight, normal, out, gamma);
                        return;
                    case FDM::FluxSplitter::LaxWendroff:
                        Flux::LaxWendroffFlux::flux(
                            qLeft, qRight, normal, dt,
                            Flux::faceSpacing(field, i, j, k, d, normal),
                            out, gamma);
                        return;
                    }
                    throw std::runtime_error(
                        "Unknown reconstructed-face numerical flux.");
                }, gamma);
        },
        fallback,
        ibmBoundary,
        requestedILWOrder,
        gamma,
        &thermodynamics);
    Math::assembleConvectiveFluxResidual(field, fluxField, residual);
}

} // namespace DivDetail

// ─────────────────────────────────────────────
//  WENO3 — 三阶 WENO 对流
// ─────────────────────────────────────────────
namespace WENO3 {

/// WENO3 对流项 (需要 ghost ≥ 2)
/// @param field 物理场。
/// @param fluxMethod 面通量方法：StegerWarming/LaxFriedrichs/Roe/LaxWendroff。
/// @param dt 当前时间步；仅LaxWendroff使用。
/// @param ibmBoundary IBM近壁WENO模板处理方式。
inline void div(Field& field, FluxField& fluxField, Residual& residual,
                FDM::FluxSplitter fluxMethod,
                double dt,
                FDM::IBMBoundaryScheme ibmBoundary,
                int requestedILWOrder,
                double gamma,
                const Physics::EquationSet::Model& thermodynamics) {
    DivDetail::wenoDiv<4>(field, fluxField, residual, fluxMethod, dt, ibmBoundary,
                          requestedILWOrder, "WENO3",
                          Math::minGhostWENO3, gamma, thermodynamics,
                          ::SF::WENO3::weno3_core);
}

} // namespace WENO3

// ─────────────────────────────────────────────
//  WENO5 — 五阶 WENO 对流
// ─────────────────────────────────────────────
namespace WENO5 {

/// WENO5 对流项 (需要 ghost ≥ 3)
/// @param field 物理场。
/// @param fluxMethod 面通量方法：StegerWarming/LaxFriedrichs/Roe/LaxWendroff。
/// @param dt 当前时间步；仅LaxWendroff使用。
/// @param ibmBoundary IBM近壁WENO模板处理方式。
inline void div(Field& field, FluxField& fluxField, Residual& residual,
                FDM::FluxSplitter fluxMethod,
                double dt,
                FDM::IBMBoundaryScheme ibmBoundary,
                int requestedILWOrder,
                double gamma,
                const Physics::EquationSet::Model& thermodynamics) {
    DivDetail::wenoDiv<6>(field, fluxField, residual, fluxMethod, dt, ibmBoundary,
                          requestedILWOrder, "WENO5",
                          Math::minGhostWENO5, gamma, thermodynamics,
                          ::SF::WENO5::weno5_core);
}

} // namespace WENO5

// ─────────────────────────────────────────────
//  TENO5 — 五阶 TENO 对流
// ─────────────────────────────────────────────
namespace TENO5 {

/// TENO5 对流项 (需要 ghost ≥ 3)
/// @param field 物理场。
/// @param fluxMethod 面通量方法：StegerWarming/LaxFriedrichs/Roe/LaxWendroff。
/// @param dt 当前时间步；仅LaxWendroff使用。
/// @param ibmBoundary IBM近壁WENO/TENO模板处理方式。
inline void div(Field& field, FluxField& fluxField, Residual& residual,
                FDM::FluxSplitter fluxMethod,
                double dt,
                FDM::IBMBoundaryScheme ibmBoundary,
                int requestedILWOrder,
                double gamma,
                const Physics::EquationSet::Model& thermodynamics) {
    DivDetail::wenoDiv<6>(field, fluxField, residual, fluxMethod, dt, ibmBoundary,
                          requestedILWOrder, "TENO5",
                          Math::minGhostWENO5, gamma, thermodynamics,
                          ::SF::TENO5::teno5_core);
}

} // namespace TENO5

// ─────────────────────────────────────────────
//  WENO7 — 七阶 WENO 对流
// ─────────────────────────────────────────────
namespace WENO7 {

/// WENO7 对流项 (需要 ghost ≥ 4)
/// @param field 物理场。
/// @param fluxMethod 面通量方法：StegerWarming/LaxFriedrichs/Roe/LaxWendroff。
/// @param dt 当前时间步；仅LaxWendroff使用。
/// @param ibmBoundary IBM近壁WENO模板处理方式。
inline void div(Field& field, FluxField& fluxField, Residual& residual,
                FDM::FluxSplitter fluxMethod,
                double dt,
                FDM::IBMBoundaryScheme ibmBoundary,
                int requestedILWOrder,
                double gamma,
                const Physics::EquationSet::Model& thermodynamics) {
    DivDetail::wenoDiv<8>(field, fluxField, residual, fluxMethod, dt, ibmBoundary,
                          requestedILWOrder, "WENO7",
                          Math::minGhostWENO7, gamma, thermodynamics,
                          ::SF::WENO7::weno7_core);
}

} // namespace WENO7

/// @brief 按配置把抽象 `div(flux)` 绑定到具体对流格式。
inline void divDispatch(
        Field& field, FluxField& fluxField, Residual& residual,
        FDM::ConvectionScheme scheme,
        FDM::FluxSplitter flux,
        double timeStep,
        FDM::IBMBoundaryScheme ibmBoundary,
        int ilwOrder,
        double idealGasGamma,
        const Physics::EquationSet::Model& thermodynamics) {
    switch (scheme) {
        case FDM::ConvectionScheme::WENO3:
            WENO3::div(field, fluxField, residual, flux, timeStep, ibmBoundary,
                       ilwOrder, idealGasGamma, thermodynamics);
            return;
        case FDM::ConvectionScheme::WENO5:
            WENO5::div(field, fluxField, residual, flux, timeStep, ibmBoundary,
                       ilwOrder, idealGasGamma, thermodynamics);
            return;
        case FDM::ConvectionScheme::TENO5:
            TENO5::div(field, fluxField, residual, flux, timeStep, ibmBoundary,
                       ilwOrder, idealGasGamma, thermodynamics);
            return;
        case FDM::ConvectionScheme::WENO7:
            WENO7::div(field, fluxField, residual, flux, timeStep, ibmBoundary,
                       ilwOrder, idealGasGamma, thermodynamics);
            return;
    }
    throw std::runtime_error("Unknown convection discretization.");
}

} // namespace SF
