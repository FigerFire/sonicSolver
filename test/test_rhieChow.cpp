#include "solver/discretization/pressure/SF_rhieChow.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace {
void requireClose(double actual,double expected) {
    if (std::abs(actual-expected)>1e-14)
        throw std::runtime_error("Rhie-Chow face contract failed.");
}
}

int main() {
    using SF::Pressure::RhieChow::predict;
    using SF::Pressure::RhieChow::correct;
    const std::array<double,3> area{2.0,0.0,0.0};
    const std::array<double,3> zero{0.0,0.0,0.0};
    const std::array<double,3> uniformU{3.0,0.0,0.0};
    const auto uniform=predict(uniformU,uniformU,zero,zero,
                               area,7.0,7.0,0.25,0.5);
    requireClose(uniform.flux,6.0);
    requireClose(uniform.correctionCoefficient,1.0);
    requireClose(correct(uniform.flux,uniform.correctionCoefficient,0,0),
                 uniform.flux);

    // A linear pressure field is already represented by the cell gradient.
    const std::array<double,3> linearGradient{2.0,0.0,0.0};
    const auto linear=predict(zero,zero,linearGradient,linearGradient,
                              area,1.0,2.0,0.25,0.5);
    requireClose(linear.flux,0.0);

    // A collocated checkerboard gradient vanishes at cell centers, but the
    // direct face difference remains in the pressure response.
    const auto checkerboard=predict(zero,zero,zero,zero,
                                    area,-1.0,1.0,0.25,0.5);
    requireClose(checkerboard.flux,-2.0);
    requireClose(correct(0.0,checkerboard.correctionCoefficient,0.0,1.0),
                 -1.0);
}
