/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_viscous.h
/// @brief 中心差分黏性通量离散核。
///
/// 本文件实现 CENTRAL2/CENTRAL4 共享的黏性通量组装工具，由
/// `discretization/diffusion/SF_laplacian.h` 调度调用。

#include "core/interfaces/SF_transportModel.h"
#include "SF_field.h"
#include "core/residual/SF_residual.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_utility.h"
#include "methods/math/SF_newtonian.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF {
namespace Viscous {

enum class CentralOrder { SECOND = 2, FOURTH = 4 };

struct TransportProperties {
    double mu = 0.0;
    double pr = 0.72;
    double idealGasGamma = 1.4;
    double idealGasConstant = 287.05;
    const FDM::ITransportModel* transportModel = nullptr;
};

struct Primitive {
    double rho = 1.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    double p = 0.0;
    double T = 0.0;
    double thermalConductivity = 0.0;
};

inline Primitive primitiveAt(const SF::Field& field, int i, int j, int k,
                             const TransportProperties& transport) {
    Primitive q;
    if (field.hasStateModel()) {
        const auto state = field.thermodynamicState(i,j,k);
        q.rho = state.density;
        q.p = state.pressure;
        q.u = state.velocity[0];
        q.v = state.velocity[1];
        q.w = state.velocity[2];
        q.T = state.temperature;
        q.thermalConductivity = state.thermalConductivity;
        return q;
    }
    q.rho = field(i, j, k, RHO);
    q.p = Numerics::requirePhysicalState(
        "viscous primitiveAt", q.rho,
        field(i, j, k, RU), field(i, j, k, RV), field(i, j, k, RW),
        field(i, j, k, E), transport.idealGasGamma);
    q.u   = field(i, j, k, RU) / q.rho;
    q.v   = field(i, j, k, RV) / q.rho;
    q.w   = field(i, j, k, RW) / q.rho;
    q.T = q.p / (q.rho * transport.idealGasConstant);
    return q;
}

inline double primitiveComponent(const SF::Field& field, int i, int j, int k,
                                 int c,
                                 const TransportProperties& transport) {
    Primitive q = primitiveAt(field, i, j, k, transport);
    if (c == 0) return q.u;
    if (c == 1) return q.v;
    if (c == 2) return q.w;
    return q.T;
}

inline double dynamicViscosityAt(const SF::Field& field,
                                 int i, int j, int k,
                                 const TransportProperties& transport) {
    double value = transport.mu;
    if (field.hasStateModel()) {
        value = field.thermodynamicState(i,j,k).dynamicViscosity;
    }
    if (transport.transportModel) {
        value = transport.transportModel->dynamicViscosity(field, i, j, k,
                                                           value);
    }
    if (!std::isfinite(value) || value < 0.0) {
        std::cerr << "[SF FATAL] invalid dynamic viscosity at ("
                  << i << "," << j << "," << k << "): mu=" << value
                  << std::endl;
        std::exit(1);
    }
    return value;
}

inline bool hasStencil(const SF::Field& field, int i, int j, int k,
                       Math::Dir d, int radius) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    int im = i - radius * di, ip = i + radius * di;
    int jm = j - radius * dj, jp = j + radius * dj;
    int km = k - radius * dk, kp = k + radius * dk;
    return im >= 0 && ip < field.MX()
        && jm >= 0 && jp < field.MY()
        && km >= 0 && kp < field.MZ();
}

inline double centralDerivative2(const SF::Field& field, int i, int j, int k,
                                 Math::Dir d, int component,
                                 const TransportProperties& transport) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);

    double fp = primitiveComponent(
        field, i + di, j + dj, k + dk, component, transport);
    double fm = primitiveComponent(
        field, i - di, j - dj, k - dk, component, transport);
    return 0.5 * (fp - fm);
}

inline double centralDerivative4(const SF::Field& field, int i, int j, int k,
                                 Math::Dir d, int component,
                                 const TransportProperties& transport) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);

    double fpp = primitiveComponent(
        field, i + 2 * di, j + 2 * dj, k + 2 * dk, component, transport);
    double fp  = primitiveComponent(
        field, i + di, j + dj, k + dk, component, transport);
    double fm  = primitiveComponent(
        field, i - di, j - dj, k - dk, component, transport);
    double fmm = primitiveComponent(
        field, i - 2 * di, j - 2 * dj, k - 2 * dk, component, transport);
    return (-fpp + 8.0 * fp - 8.0 * fm + fmm) / 12.0;
}

inline double derivative(const SF::Field& field, int i, int j, int k,
                         Math::Dir d, int component, CentralOrder order,
                         const TransportProperties& transport) {
    if (!Math::isDirectionActive(d)) return 0.0;

    if (order == CentralOrder::FOURTH) {
        if (!hasStencil(field, i, j, k, d, 2)) {
            throw std::runtime_error(
                "CENTRAL4 stencil unavailable at ("
                + std::to_string(i) + "," + std::to_string(j) + ","
                + std::to_string(k) + "), direction="
                + std::to_string(static_cast<int>(d))
                + ". Select viscous type CENTRAL2 explicitly for a "
                  "second-order boundary closure.");
        }
        return centralDerivative4(
            field, i, j, k, d, component, transport);
    }
    return centralDerivative2(field, i, j, k, d, component, transport);
}

inline Vector3 physicalGradientAt(const SF::Field& field,
                                  int i, int j, int k,
                                  int component,
                                  CentralOrder order,
                                  const TransportProperties& transport) {
    const double dXi = derivative(
        field, i, j, k, Math::XI, component, order, transport);
    const double dEt = derivative(
        field, i, j, k, Math::ETA, component, order, transport);
    const double dZe = derivative(
        field, i, j, k, Math::ZETA, component, order, transport);

    return Numerics::VectorCalculus::physicalGradient(
        field, dXi, dEt, dZe, i, j, k);
}

/// @brief 仅速度梯度调用的无状态 ideal-gas 兼容入口。
inline Vector3 physicalGradientAt(const SF::Field& field,
                                  int i, int j, int k,
                                  int component,
                                  CentralOrder order) {
    return physicalGradientAt(
        field, i, j, k, component, order, TransportProperties{});
}

inline Vector3 faceGradient(const SF::Field& field,
                            int i, int j, int k,
                            Math::Dir faceDir,
                            int component,
                            CentralOrder order,
                            const TransportProperties& transport) {
    int di, dj, dk;
    Math::dirOffset(faceDir, di, dj, dk);
    const Vector3 left =
        physicalGradientAt(field, i, j, k, component, order, transport);
    const Vector3 right =
        physicalGradientAt(field, i + di, j + dj, k + dk,
                           component, order, transport);
    return {
        0.5 * (left.x + right.x),
        0.5 * (left.y + right.y),
        0.5 * (left.z + right.z)
    };
}

/// @brief 单元中心速度梯度，第一下标为速度分量。
inline Tensor3 velocityGradientAt(
        const SF::Field& field,
        int i, int j, int k,
        CentralOrder order,
        const TransportProperties& transport = TransportProperties{}) {
    Tensor3 result;
    for (int component = 0; component < 3; ++component) {
        result[(size_t)component] = physicalGradientAt(
            field, i, j, k, component, order, transport);
    }
    return result;
}

/// @brief 面中心速度梯度，由相邻单元中心梯度算术平均获得。
inline Tensor3 faceVelocityGradient(
        const SF::Field& field,
        int i, int j, int k,
        Math::Dir faceDir,
        CentralOrder order,
        const TransportProperties& transport) {
    Tensor3 result;
    for (int component = 0; component < 3; ++component) {
        result[(size_t)component] = faceGradient(
            field, i, j, k, faceDir, component, order, transport);
    }
    return result;
}

inline Primitive facePrimitive(const SF::Field& field, int i, int j, int k,
                               Math::Dir faceDir,
                               const TransportProperties& transport) {
    int di, dj, dk;
    Math::dirOffset(faceDir, di, dj, dk);
    Primitive a = primitiveAt(field, i, j, k, transport);
    Primitive b = primitiveAt(
        field, i + di, j + dj, k + dk, transport);
    Primitive f;
    f.rho = 0.5 * (a.rho + b.rho);
    f.u   = 0.5 * (a.u + b.u);
    f.v   = 0.5 * (a.v + b.v);
    f.w   = 0.5 * (a.w + b.w);
    f.p   = 0.5 * (a.p + b.p);
    f.T   = 0.5 * (a.T + b.T);
    f.thermalConductivity =
        0.5 * (a.thermalConductivity + b.thermalConductivity);
    return f;
}

inline void viscousFluxAtFace(const SF::Field& field, int i, int j, int k,
                              Math::Dir faceDir, CentralOrder order,
                              const TransportProperties& transport,
                              std::vector<double>& flux) {
    flux.assign((size_t)field.NVar(), 0.0);
    Primitive q = facePrimitive(field, i, j, k, faceDir, transport);

    const Tensor3 gradU = faceVelocityGradient(
        field, i, j, k, faceDir, order, transport);
    const Vector3 gradT =
        faceGradient(field, i, j, k, faceDir, 3, order, transport);

    int di, dj, dk;
    Math::dirOffset(faceDir, di, dj, dk);
    double localMu = 0.5 * (dynamicViscosityAt(field, i, j, k, transport)
                          + dynamicViscosityAt(field, i + di, j + dj, k + dk, transport));
    if (!std::isfinite(transport.pr) || transport.pr <= 1.0e-12) {
        std::cerr << "[SF FATAL] invalid Prandtl number for heat conduction: Pr="
                  << transport.pr << std::endl;
        std::exit(1);
    }
    double kappa = 0.0;
    if (field.hasStateModel()) {
        kappa = q.thermalConductivity;
        if (!std::isfinite(kappa) || kappa <= 0.0) {
            std::cerr << "[SF FATAL] FluidStateModel viscous heat flux requires "
                      << "finite positive phase thermal conductivity, got k="
                      << kappa << std::endl;
            std::exit(1);
        }
    } else {
        double cp = transport.idealGasGamma * transport.idealGasConstant
            / (transport.idealGasGamma - 1.0);
        kappa = localMu * cp / transport.pr;
    }

    const SymmTensor3 stress =
        Constitutive::Newtonian::stress(localMu, gradU);

    double cofactor[4];
    Math::faceMetrics(field, i, j, k, faceDir, cofactor);
    const double cx = cofactor[0];
    const double cy = cofactor[1];
    const double cz = cofactor[2];
    const Vector3 area(cx, cy, cz);
    const Vector3 viscousTraction =
        Constitutive::Newtonian::traction(stress, area);
    const int momentum = field.hasStateModel()
        ? field.stateModel()->momentumIndex(0) : RU;
    const int energy = field.hasStateModel()
        ? field.stateModel()->energyIndex() : E;
    flux[(size_t)momentum] = viscousTraction.x;
    flux[(size_t)momentum+1] = viscousTraction.y;
    flux[(size_t)momentum+2] = viscousTraction.z;

    const Vector3 velocity(q.u, q.v, q.w);
    flux[(size_t)energy] = dot(velocity, viscousTraction)
                         + kappa * dot(area, gradT);
    for (int n = 0; n < field.NVar(); ++n) {
        if (!std::isfinite(flux[(size_t)n])) {
            std::cerr << "[SF FATAL] non-finite viscous conservative flux at face ("
                      << i << "," << j << "," << k << "), component="
                      << n << ", value=" << flux[(size_t)n] << std::endl;
            std::exit(1);
        }
    }
}

inline void subtractFlux(SF::Field& field, SF::Residual& residual,
                         int i, int j, int k,
                         Math::Dir d, const std::vector<double>& flux) {
    if ((int)flux.size() != field.NVar()) {
        throw std::runtime_error("Viscous flux size differs from Field NVar.");
    }
    for (int v = 0; v < field.NVar(); ++v) {
        if (d == Math::XI)       residual.x(i, j, k, v) -= flux[(size_t)v];
        else if (d == Math::ETA) residual.y(i, j, k, v) -= flux[(size_t)v];
        else                     residual.z(i, j, k, v) -= flux[(size_t)v];
    }
}

inline void computeCentralViscousRHS(SF::Field& field, SF::Residual& residual,
                                     CentralOrder order,
                                     const TransportProperties& transport) {
    auto computeDir = [&](Math::Dir d) {
        if (!Math::isDirectionActive(d)) return;
        Math::forFaces(field, d, [&](int i, int j, int k) {
            std::vector<double> flux;
            viscousFluxAtFace(field, i, j, k, d, order, transport, flux);
            subtractFlux(field, residual, i, j, k, d, flux);
        });
    };

    computeDir(Math::XI);
    computeDir(Math::ETA);
    computeDir(Math::ZETA);
}

} // namespace Viscous
} // namespace SF
