/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_laplacian.h
/// @brief 粘性扩散项离散算子 (viscous flux divergence).
///
/// 每个中心差分阶数是一个独立的命名空间，内部提供统一的 laplacian() 接口。
/// 调用方式:
/// @code
///   CENTRAL2::laplacian(field, mu, Pr);
///   CENTRAL4::laplacian(field, mu, Pr);
/// @endcode
///
/// laplacian 计算 ∇·F_v(Q, ∇Q) 并将粘性通量减入 ResX/ResY/ResZ，
/// 后续由 Math::spatialResidual 取散度，符号为 −∇·F_v。

#include "SF_field.h"
#include "core/residual/SF_residual.h"
#include "SF_transportModel.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_viscous.h"

#include <cstdlib>
#include <string>

namespace SF {

// ─────────────────────────────────────────────
//  CENTRAL2 — 二阶中心差分粘性项
// ─────────────────────────────────────────────
namespace CENTRAL2 {

/// @brief Assemble second-order central viscous flux divergence.
/// @param field Field whose residual flux arrays are updated in-place.
/// @param mu Dynamic viscosity.
/// @param Pr Prandtl number.
/// @param enabled Enables/disables viscous assembly without reading globals.
/// @param transportModel 可选的有效粘度/湍流闭合接口。
inline void laplacian(Field& field, Residual& residual, double mu, double Pr, bool enabled,
                      double idealGasGamma,
                      double idealGasConstant,
                      const FDM::ITransportModel* transportModel = nullptr);

/// 二阶中心差分粘性项 (需要 ghost ≥ 1)
/// @param field  物理场
/// @param mu     动力粘度
/// @param Pr     普朗特数 (默认 0.72)
inline void laplacian(Field& field, Residual& residual, double mu, double Pr = 0.72) {
    laplacian(field, residual, mu, Pr, true, 1.4, 287.05);
}

/// @brief Assemble second-order central viscous flux divergence.
/// @param field Field whose residual flux arrays are updated in-place.
/// @param mu Dynamic viscosity.
/// @param Pr Prandtl number.
/// @param enabled Enables/disables viscous assembly without reading globals.
/// @param transportModel 可选的有效粘度/湍流闭合接口。
inline void laplacian(Field& field, Residual& residual, double mu, double Pr, bool enabled,
                      double idealGasGamma,
                      double idealGasConstant,
                      const FDM::ITransportModel* transportModel) {
    if (!enabled) return;

    if (!Math::checkGhostDepth(field, 1, "CENTRAL2 viscous")) {
        std::exit(1);
    }

    Viscous::TransportProperties transport;
    transport.mu = mu;
    transport.pr = Pr;
    transport.idealGasGamma = idealGasGamma;
    transport.idealGasConstant = idealGasConstant;
    transport.transportModel = transportModel;

    Viscous::computeCentralViscousRHS(field, residual, Viscous::CentralOrder::SECOND, transport);
}

} // namespace CENTRAL2

// ─────────────────────────────────────────────
//  CENTRAL4 — 四阶中心差分粘性项
// ─────────────────────────────────────────────
namespace CENTRAL4 {

/// @brief Assemble fourth-order central viscous flux divergence.
/// @param field Field whose residual flux arrays are updated in-place.
/// @param mu Dynamic viscosity.
/// @param Pr Prandtl number.
/// @param enabled Enables/disables viscous assembly without reading globals.
/// @param transportModel 可选的有效粘度/湍流闭合接口。
inline void laplacian(Field& field, Residual& residual, double mu, double Pr, bool enabled,
                      double idealGasGamma,
                      double idealGasConstant,
                      const FDM::ITransportModel* transportModel = nullptr);

/// 四阶中心差分粘性项 (需要 ghost ≥ 2)
/// @param field  物理场
/// @param mu     动力粘度
/// @param Pr     普朗特数 (默认 0.72)
inline void laplacian(Field& field, Residual& residual, double mu, double Pr = 0.72) {
    laplacian(field, residual, mu, Pr, true, 1.4, 287.05);
}

/// @brief Assemble fourth-order central viscous flux divergence.
/// @param field Field whose residual flux arrays are updated in-place.
/// @param mu Dynamic viscosity.
/// @param Pr Prandtl number.
/// @param enabled Enables/disables viscous assembly without reading globals.
/// @param transportModel 可选的有效粘度/湍流闭合接口。
inline void laplacian(Field& field, Residual& residual, double mu, double Pr, bool enabled,
                      double idealGasGamma,
                      double idealGasConstant,
                      const FDM::ITransportModel* transportModel) {
    if (!enabled) return;

    if (!Math::checkGhostDepth(field, 2, "CENTRAL4 viscous")) {
        std::exit(1);
    }

    Viscous::TransportProperties transport;
    transport.mu = mu;
    transport.pr = Pr;
    transport.idealGasGamma = idealGasGamma;
    transport.idealGasConstant = idealGasConstant;
    transport.transportModel = transportModel;

    Viscous::computeCentralViscousRHS(field, residual, Viscous::CentralOrder::FOURTH, transport);
}

} // namespace CENTRAL4

/// @brief 按配置把抽象 `diffusion(coefficient, unknown)` 绑定到中心格式。
inline void laplacianDispatch(
        Field& field, Residual& residual,
        FDM::ViscousScheme scheme,
        bool enabled,
        double dynamicViscosity,
        double prandtl = 0.72,
        double idealGasGamma = 1.4,
        double idealGasConstant = 287.05,
        const FDM::ITransportModel* transport = nullptr) {
    if (scheme == FDM::ViscousScheme::Central4) {
        CENTRAL4::laplacian(
            field, residual, dynamicViscosity, prandtl, enabled,
            idealGasGamma, idealGasConstant, transport);
        return;
    }
    CENTRAL2::laplacian(
        field, residual, dynamicViscosity, prandtl, enabled,
        idealGasGamma, idealGasConstant, transport);
}

} // namespace SF
