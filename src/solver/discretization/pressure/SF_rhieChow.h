#pragma once

/// @file SF_rhieChow.h
/// @brief Collocated pressure/velocity face coupling, independent of PISO.

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace SF::PressureBased::RhieChow {

struct FaceResult {
    double flux = 0.0;
    double correctionCoefficient = 0.0;
};

inline FaceResult predict(const std::array<double,3>& lowerVelocity,
                          const std::array<double,3>& upperVelocity,
                          const std::array<double,3>& lowerGradient,
                          const std::array<double,3>& upperGradient,
                          const std::array<double,3>& area,
                          double lowerPressure, double upperPressure,
                          double inverseMomentumDiagonal,
                          double spacing) {
    if (!(spacing>0.0) || !(inverseMomentumDiagonal>0.0))
        throw std::runtime_error("Rhie-Chow face geometry/coefficient is invalid.");
    double areaSquared=0.0, velocityFlux=0.0, interpolatedGradient=0.0;
    for (int c=0;c<3;++c) {
        areaSquared+=area[(size_t)c]*area[(size_t)c];
        velocityFlux+=0.5*(lowerVelocity[(size_t)c]
            +upperVelocity[(size_t)c])*area[(size_t)c];
        interpolatedGradient+=0.5*(lowerGradient[(size_t)c]
            +upperGradient[(size_t)c])*area[(size_t)c];
    }
    const double directGradient=(upperPressure-lowerPressure)
        *std::sqrt(areaSquared)/spacing;
    return {velocityFlux-inverseMomentumDiagonal
        *(directGradient-interpolatedGradient),
        inverseMomentumDiagonal*std::sqrt(areaSquared)/spacing};
}

inline double correct(double predictedFlux,double response,
                      double lowerCorrection,double upperCorrection) {
    return predictedFlux-response*(upperCorrection-lowerCorrection);
}

} // namespace SF::PressureBased::RhieChow
