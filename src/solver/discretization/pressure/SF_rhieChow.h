#pragma once

/// @file SF_rhieChow.h
/// @brief Collocated pressure/velocity face coupling, independent of PISO.

#include "core/state/SF_distributedField.h"
#include <vector>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace SF::Pressure::RhieChow {

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

/// @brief 绑定后的面几何/边界信息；不含压力算法、时钟或模型对象。
struct Face {
    int lower,upper,boundary;std::size_t index;std::array<double,3> area;
    double spacing;bool twoSided,pressureOpen;
};
/// @brief 独立通量数值操作；调用方绑定实际输入并负责 canonical COPY/发布。
template<class Gradient>
void construct(const State::DistributedFieldView& velocity,const State::DistributedFieldView& pressure,
        const std::vector<double>& response,const std::vector<Face>& faces,Gradient gradient,
        std::vector<double>& flux,std::vector<double>& faceResponse) {
    std::fill(flux.begin(),flux.end(),0.0);std::fill(faceResponse.begin(),faceResponse.end(),0.0);
    for(const auto& f:faces) {
        double q=0;
        if(f.twoSided || f.pressureOpen) {
            std::array<double,3> a{},b{};for(int c=0;c<3;++c){a[c]=velocity.read(f.lower,c);b[c]=velocity.read(f.upper,c);}
            const auto r=f.twoSided?0.5*(response.at(f.lower)+response.at(f.upper)):response.at(f.boundary==f.lower?f.upper:f.lower);
            const auto result=predict(a,b,gradient(f.lower),gradient(f.upper),f.area,pressure.read(f.lower,0),pressure.read(f.upper,0),r,f.spacing);
            q=result.flux;faceResponse.at(f.index)=result.correctionCoefficient;
        } else for(int c=0;c<3;++c)q+=velocity.read(f.boundary,c)*f.area[c];
        flux.at(f.index)=q;
    }
}

} // namespace SF::Pressure::RhieChow
