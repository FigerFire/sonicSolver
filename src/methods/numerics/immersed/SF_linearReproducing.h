#pragma once
/// @file SF_linearReproducing.h
/// @brief 正权 raw kernel 的 affine moment 解；修正后权重可有符号，不做 clipping。
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::FDM::Immersed {
/// @brief 解 sum(r_i b_i b_i^T) a=e0，使 w_i=r_i b_i^T a 重现常量与三维线性场。
/// @param matrix 完整归并的 4x4 row-major Gram matrix，b=[1,(x-X)/R]。
/// @return affine 权重系数。支撑不满秩或近退化时直接报错，不正则化/降阶。
inline std::array<double,4> linearMomentCoefficients(const std::array<double,16>& matrix) {
    std::array<double,4> scale{},rhs{},solution{};
    std::array<double,16> lower{};
    for (double value:matrix) if (!std::isfinite(value))
        throw std::runtime_error("linearReproducing requires finite global moments.");
    for (int i=0;i<4;++i) {
        if (matrix[4*i+i]<=0.)
            throw std::runtime_error("linearReproducing requires full-rank 3D affine support; zero moment diagonal.");
        scale[i]=std::sqrt(matrix[4*i+i]);
    }
    // Equilibration removes units/thin EMPTY-axis scaling, not a singular mode.
    for (int i=0;i<4;++i) for (int j=0;j<=i;++j) {
        double value=matrix[4*i+j]/scale[i]/scale[j];
        if (std::abs(value-matrix[4*j+i]/scale[j]/scale[i])>1e-12)
            throw std::runtime_error("linearReproducing Gram matrix is not symmetric.");
        for (int k=0;k<j;++k) value-=lower[4*i+k]*lower[4*j+k];
        if (i==j) {
            if (value<=128.*std::numeric_limits<double>::epsilon())
                throw std::runtime_error("linearReproducing affine support is rank deficient or ill-conditioned.");
            lower[4*i+j]=std::sqrt(value);
        } else lower[4*i+j]=value/lower[4*j+j];
    }
    rhs[0]=1./scale[0];
    for (int i=0;i<4;++i) {
        for (int j=0;j<i;++j) rhs[i]-=lower[4*i+j]*rhs[j];
        rhs[i]/=lower[4*i+i];
    }
    for (int i=3;i>=0;--i) {
        double value=rhs[i];
        for (int j=i+1;j<4;++j) value-=lower[4*j+i]*solution[j];
        solution[i]=value/lower[4*i+i];
    }
    for (int i=0;i<4;++i) solution[i]/=scale[i];
    for (int i=0;i<4;++i) {
        double residual=-(i==0?1.:0.);
        for (int j=0;j<4;++j) residual+=matrix[4*i+j]*solution[j];
        if (!std::isfinite(solution[i]) || !std::isfinite(residual) || std::abs(residual)>1e-10)
            throw std::runtime_error("linearReproducing failed its global affine moment equations.");
    }
    return solution;
}
}
