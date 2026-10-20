#include "core/config/SF_config.h"
#include "app/application/model/SF_configParser.h"
#include "methods/numerics/Flux/SF_lxF.h"
#include "methods/numerics/Flux/SF_roe.h"
#include "methods/numerics/Flux/SF_rusanovEOS.h"
#include "methods/numerics/Flux/SF_sw.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "methods/numerics/convection/SF_WENO3.h"
#include "methods/numerics/convection/SF_WENO5.h"
#include "methods/numerics/convection/SF_TENO5.h"
#include "methods/numerics/convection/SF_WENO7.h"
#include "methods/numerics/convection/SF_riemann.h"
#include "solver/discretization/convection/SF_WENO.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>

namespace {

bool near(double a, double b, double tolerance = 1.0e-12) {
    return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

bool sameFlux(const double a[5], const double b[5], double tolerance = 1.0e-12) {
    for (int v = 0; v < 5; ++v) {
        if (!near(a[v], b[v], tolerance)) return false;
    }
    return true;
}

void perfectGasState(double rho, double u, double v, double w, double p,
                     double gamma, double q[5]) {
    q[0] = rho;
    q[1] = rho * u;
    q[2] = rho * v;
    q[3] = rho * w;
    q[4] = p / (gamma - 1.0)
         + 0.5 * rho * (u * u + v * v + w * w);
}

template <int NStencil, typename Kernel>
bool reconstructedConstantState(const double q[5], Kernel&& kernel,
                                double qLeft[5], double qRight[5],
                                const char* scheme) {
    double stencil[NStencil][5];
    for (auto& sample : stencil) {
        std::copy(q, q + 5, sample);
    }
    const double normal[3] = {1.0, 0.0, 0.0};
    double leftEigen[25], rightEigen[25];
    SF::Riemann::buildCharacteristicMatrix(
        q, normal, 1.4, leftEigen, rightEigen);

    double characteristic[NStencil][5];
    for (int s = 0; s < NStencil; ++s) {
        SF::Riemann::gemv5(leftEigen, stencil[s], characteristic[s]);
    }
    double wLeft[5], wRight[5];
    for (int v = 0; v < 5; ++v) {
        double scalarLeft[8] = {};
        double scalarRight[8] = {};
        SF::Flux::fillScalarWENOStencils<NStencil>(
            characteristic, v, scalarLeft, scalarRight);
        wLeft[v] = kernel(scalarLeft);
        wRight[v] = kernel(scalarRight);
    }
    SF::Riemann::gemv5(rightEigen, wLeft, qLeft);
    SF::Riemann::gemv5(rightEigen, wRight, qRight);

    bool preserved = true;
    for (int v = 0; v < 5; ++v) {
        preserved = preserved && near(qLeft[v], q[v]) && near(qRight[v], q[v]);
    }
    if (!preserved) {
        std::cerr << scheme << " does not preserve a constant conservative stencil.\n";
        std::cerr << "Q0:";
        for (double value : std::array<double, 5>{q[0], q[1], q[2], q[3], q[4]})
            std::cerr << ' ' << value;
        std::cerr << "\nL*Q0:";
        for (double value : wLeft) std::cerr << ' ' << value;
        std::cerr << "\nR*Wleft / R*Wright / difference:\n";
        for (int v = 0; v < 5; ++v) {
            std::cerr << "  " << v << ": " << qLeft[v] << " / "
                      << qRight[v] << " / " << qLeft[v] - q[v] << '\n';
        }
    }
    return preserved;
}

template <int NStencil, typename Kernel, typename Flux>
void reconstructedConstantFlux(const double q[5], Kernel&& kernel,
                               Flux&& flux, double out[5]) {
    double stencil[NStencil][5];
    for (auto& sample : stencil) std::copy(q, q + 5, sample);
    const double metrics[4] = {1.0, 0.0, 0.0, 1.0};
    SF::Flux::characteristicWENOFaceFlux<NStencil>(
        stencil, metrics, out, std::forward<Kernel>(kernel),
        std::forward<Flux>(flux), 1.4);
}

} // namespace

int main() {
    constexpr double gamma = 1.4;
    constexpr double gasConstant = 287.05;
    auto equations = SF::Physics::FluidStateModel::makeSingleFluidPerfectGas(
        gamma, gasConstant, 0.0, 0.72);

    double uniform[5];
    perfectGasState(1.2, 30.0, -2.0, 1.0, 101325.0, gamma, uniform);
    const double normal[3] = {1.0, 0.0, 0.0};
    const std::array<double, 3> normalArray = {1.0, 0.0, 0.0};

    double rusanovUniform[5];
    SF::Numerics::RusanovEOS::faceFlux(
        *equations, uniform, uniform, normalArray, rusanovUniform, 5);
    double physicalUniform[5] = {
        uniform[1],
        uniform[1] * 30.0 + 101325.0,
        uniform[2] * 30.0,
        uniform[3] * 30.0,
        (uniform[4] + 101325.0) * 30.0};

    double reconstructedLeft[5], reconstructedRight[5];
    if (!reconstructedConstantState<4>(
            uniform, SF::WENO3::weno3_core,
            reconstructedLeft, reconstructedRight, "WENO3")
        || !reconstructedConstantState<6>(
            uniform, SF::WENO5::weno5_core,
            reconstructedLeft, reconstructedRight, "WENO5")
        || !reconstructedConstantState<6>(
            uniform, SF::TENO5::teno5_core,
            reconstructedLeft, reconstructedRight, "TENO5")
        || !reconstructedConstantState<8>(
            uniform, SF::WENO7::weno7_core,
            reconstructedLeft, reconstructedRight, "WENO7")) {
        return 1;
    }

    if (!sameFlux(rusanovUniform, physicalUniform)) {
        std::cerr << "Rusanov consistency failed for a uniform state.\n";
        return 1;
    }

    double reconstructedRusanov[5];
    reconstructedConstantFlux<6>(
        uniform,
        SF::WENO5::weno5_core,
        [&](const double left[5], const double right[5],
            const double n[3], double out[5]) {
            SF::Numerics::RusanovEOS::faceFlux(
                *equations, left, right, {n[0], n[1], n[2]}, out, 5);
        },
        reconstructedRusanov);
    if (!sameFlux(reconstructedRusanov, physicalUniform)) {
        std::cerr << "Reconstructed Rusanov is inconsistent on a constant stencil.\n";
        for (int component = 0; component < 5; ++component) {
            std::cerr << "  component " << component
                      << ": reconstructed=" << reconstructedRusanov[component]
                      << ", physical=" << physicalUniform[component] << "\n";
        }
        return 1;
    }

    double left[5], right[5];
    perfectGasState(1.0, 2.0, 0.1, 0.0, 1.0, gamma, left);
    perfectGasState(0.125, -0.5, 0.0, 0.0, 0.1, gamma, right);

    double rusanov[5], steger[5], roe[5], laxFriedrichs[5];
    SF::Numerics::RusanovEOS::faceFlux(
        *equations, left, right, normalArray, rusanov, 5);
    SF::Flux::StegerWarmingFlux::flux(
        left, right, normal, steger, gamma);
    SF::Flux::RoeFlux::flux(left, right, normal, roe, gamma);
    SF::Flux::LaxFriedrichsFlux::flux(
        left, right, normal, laxFriedrichs, gamma);
    if (sameFlux(rusanov, steger, 1.0e-10)
        || sameFlux(rusanov, laxFriedrichs, 1.0e-10)
        || sameFlux(steger, laxFriedrichs, 1.0e-10)) {
        std::cerr << "Configured numerical fluxes are not observably distinct.\n";
        return 1;
    }

    double stegerUniform[5], roeUniform[5], laxFriedrichsUniform[5];
    SF::Flux::StegerWarmingFlux::flux(
        uniform, uniform, normal, stegerUniform, gamma);
    SF::Flux::RoeFlux::flux(
        uniform, uniform, normal, roeUniform, gamma);
    SF::Flux::LaxFriedrichsFlux::flux(
        uniform, uniform, normal, laxFriedrichsUniform, gamma);
    if (!sameFlux(stegerUniform, physicalUniform)
        || !sameFlux(roeUniform, physicalUniform)
        || !sameFlux(laxFriedrichsUniform, physicalUniform)) {
        std::cerr << "A direct equal-state numerical flux is inconsistent.\n";
        return 1;
    }

    double reconstructedSW[5], reconstructedLF[5];
    reconstructedConstantFlux<6>(
        uniform,
        SF::WENO5::weno5_core,
        [&](const double qLeft[5], const double qRight[5],
            const double n[3], double out[5]) {
            SF::Flux::StegerWarmingFlux::flux(qLeft, qRight, n, out, gamma);
        },
        reconstructedSW);
    reconstructedConstantFlux<6>(
        uniform,
        SF::WENO5::weno5_core,
        [&](const double qLeft[5], const double qRight[5],
            const double n[3], double out[5]) {
            SF::Flux::LaxFriedrichsFlux::flux(qLeft, qRight, n, out, gamma);
        },
        reconstructedLF);
    if (!sameFlux(reconstructedSW, physicalUniform)
        || !sameFlux(reconstructedLF, physicalUniform)) {
        std::cerr << "Reconstructed high-order flux is inconsistent on a constant stencil.\n";
        return 1;
    }

    if (SF::FDM::parseFluxSplitter("Rusanov")
                != SF::FDM::FluxSplitter::Rusanov
            || std::string(SF::FDM::toString(SF::FDM::FluxSplitter::Rusanov))
                != "Rusanov") {
        std::cerr << "Rusanov configuration dispatch failed.\n";
        return 1;
    }

    return 0;
}
