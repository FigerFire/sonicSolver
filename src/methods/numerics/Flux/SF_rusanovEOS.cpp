/// @file SF_rusanovEOS.cpp
/// @brief 可压缩通量或 Riemann 数值格式实现。

#include "SF_rusanovEOS.h"
#include "SF_canonicalFace.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::Numerics::RusanovEOS {
namespace {

std::vector<double> physicalFlux(
        const Physics::FluidStateModel::Model& equations,
        const double* q,
        const Physics::FluidStateModel::ThermodynamicState& state,
        const std::array<double, 3>& normal) {
    std::vector<double> flux((size_t)equations.variableCount(), 0.0);
    const double un = state.velocity[0]*normal[0]
                    + state.velocity[1]*normal[1]
                    + state.velocity[2]*normal[2];
    for (int k = 0; k < equations.densityVariableCount(); ++k)
        flux[(size_t)k] = q[k] * un;
    for (int d = 0; d < 3; ++d) {
        const int index = equations.momentumIndex(d);
        flux[(size_t)index] = q[index] * un + state.pressure * normal[(size_t)d];
    }
    const int energy = equations.energyIndex();
    flux[(size_t)energy] = (q[energy] + state.pressure) * un;
    return flux;
}

} // namespace

std::vector<double> faceFlux(
        const Physics::FluidStateModel::Model& equations,
        const double* left,
        const double* right,
        const std::array<double, 3>& unitNormal) {
    const int nVar = equations.variableCount();
    std::vector<double> result((size_t)nVar, 0.0);
    faceFlux(equations, left, right, unitNormal, result.data(), nVar);
    return result;
}

void faceFlux(
        const Physics::FluidStateModel::Model& equations,
        const double* left,
        const double* right,
        const std::array<double, 3>& unitNormal,
        double* result,
        int resultCount) {
    const int nVar = equations.variableCount();
    if (!left || !right || !result || resultCount != nVar) {
        throw std::invalid_argument(
            "EOS-aware Rusanov face flux received an invalid state/output layout.");
    }
    const auto l = equations.close(left, nVar);
    const auto r = equations.close(right, nVar);
    const auto fL = physicalFlux(equations, left, l, unitNormal);
    const auto fR = physicalFlux(equations, right, r, unitNormal);
    const double unL = l.velocity[0]*unitNormal[0]
                     + l.velocity[1]*unitNormal[1]
                     + l.velocity[2]*unitNormal[2];
    const double unR = r.velocity[0]*unitNormal[0]
                     + r.velocity[1]*unitNormal[1]
                     + r.velocity[2]*unitNormal[2];
    const double lambda = std::max(std::abs(unL) + l.soundSpeed,
                                   std::abs(unR) + r.soundSpeed);
    if (!std::isfinite(lambda) || lambda <= 0.0)
        throw std::runtime_error("EOS-aware Rusanov received invalid wave speed.");
    for (int v = 0; v < nVar; ++v) {
        result[v] = 0.5 * (fL[(size_t)v] + fR[(size_t)v])
                  - 0.5 * lambda * (right[v] - left[v]);
    }
}

void div(Field& field, FluxField& fluxField, Residual& residual,
         const Physics::FluidStateModel::Model& equations) {
    if (field.NVar() != equations.variableCount())
        throw std::runtime_error("Field and FluidStateModel variable counts differ.");
    fluxField.clear();
    auto direction = [&](Math::Dir dir) {
        Math::forFaces(field, dir, [&](int i, int j, int k) {
            if (!Math::shouldCalculateFaceFlux(field, dir, i, j, k)) return;
            int di=0, dj=0, dk=0;
            Math::dirOffset(dir, di, dj, dk);
            std::vector<double> left((size_t)field.NVar());
            std::vector<double> right((size_t)field.NVar());
            for (int v = 0; v < field.NVar(); ++v) {
                left[(size_t)v] = field(i, j, k, v);
                right[(size_t)v] = field(i+di, j+dj, k+dk, v);
            }
            const auto geometry=CanonicalFace::geometry(
                field,static_cast<int>(dir),i,j,k);
            std::vector<double> flux((size_t)field.NVar(), 0.0);
            faceFlux(equations, left.data(), right.data(),
                     geometry.unitNormal, flux.data(), field.NVar());
            for(double& value:flux)value*=geometry.area;
            Math::storeFlux(fluxField, field, i, j, k, dir, flux);
        });
    };
    direction(Math::XI);
    direction(Math::ETA);
    direction(Math::ZETA);
    Math::assembleConvectiveFluxResidual(field, fluxField, residual);
}

double deltaT(const Field& field,
              const Physics::FluidStateModel::Model& equations,
              double cfl) {
    if (field.NVar() != equations.variableCount()
        || !std::isfinite(cfl) || cfl <= 0.0)
        throw std::runtime_error("EOS-aware deltaT received invalid input.");
    double result = std::numeric_limits<double>::max();
    std::vector<double> q((size_t)field.NVar(), 0.0);
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) q[(size_t)v] = field(i,j,k,v);
        const auto state = equations.close(q.data(), field.NVar());
        const double ux = std::abs(state.velocity[0]*field.XiX(i,j,k)
                                 + state.velocity[1]*field.XiY(i,j,k)
                                 + state.velocity[2]*field.XiZ(i,j,k));
        const double ue = std::abs(state.velocity[0]*field.EtX(i,j,k)
                                 + state.velocity[1]*field.EtY(i,j,k)
                                 + state.velocity[2]*field.EtZ(i,j,k));
        const double uz = std::abs(state.velocity[0]*field.ZeX(i,j,k)
                                 + state.velocity[1]*field.ZeY(i,j,k)
                                 + state.velocity[2]*field.ZeZ(i,j,k));
        double spectralRadius = 0.0;
        if (Math::isDirectionActive(Math::XI))
            spectralRadius += ux + state.soundSpeed*Math::vecMag(
                field.XiX(i,j,k), field.XiY(i,j,k), field.XiZ(i,j,k));
        if (Math::isDirectionActive(Math::ETA))
            spectralRadius += ue + state.soundSpeed*Math::vecMag(
                field.EtX(i,j,k), field.EtY(i,j,k), field.EtZ(i,j,k));
        if (Math::isDirectionActive(Math::ZETA))
            spectralRadius += uz + state.soundSpeed*Math::vecMag(
                field.ZeX(i,j,k), field.ZeY(i,j,k), field.ZeZ(i,j,k));
        if (!std::isfinite(spectralRadius) || spectralRadius <= 0.0)
            throw std::runtime_error("EOS-aware deltaT found zero/invalid spectral radius.");
        result = std::min(result, cfl/spectralRadius);
    });
    if (!std::isfinite(result) || result == std::numeric_limits<double>::max())
        throw std::runtime_error("EOS-aware deltaT found no valid fluid cell.");
    return result;
}

} // namespace SF::Numerics::RusanovEOS
