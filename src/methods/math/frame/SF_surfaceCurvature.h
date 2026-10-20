/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#pragma once

/// @file SF_surfaceCurvature.h
/// @brief 局部曲面曲率张量拟合工具。

#include <algorithm>
#include <cmath>
#include <vector>

namespace SF {
namespace Math {
namespace SurfaceCurvature {

/// @brief 局部切平面曲率张量。
struct Tensor {
    double k11 = 0.0;
    double k12 = 0.0;
    double k22 = 0.0;
};

/// @brief 参与曲率拟合的局部法向样本。
struct NormalSample {
    double a = 0.0;
    double b = 0.0;
    double dn1 = 0.0;
    double dn2 = 0.0;
    double r2 = 0.0;
};

/// @brief 由邻近法向变化拟合曲率张量。
/// @param normals 已过滤后的局部法向样本。
/// @param useT2 true表示使用二维切平面张量，false表示二维计算平面。
/// @param curvature 输出曲率张量。
/// @return 拟合矩阵非奇异时返回true。
inline bool fitTensor(const std::vector<NormalSample>& normals,
                      bool useT2,
                      Tensor& curvature) {
    if (normals.empty()) return false;

    if (!useT2) {
        double ata = 0.0;
        double atb = 0.0;
        for (const auto& sample : normals) {
            ata += sample.a * sample.a;
            atb += sample.a * sample.dn1;
        }
        if (ata <= 1e-20) return false;
        curvature.k11 = atb / ata;
        curvature.k12 = 0.0;
        curvature.k22 = 0.0;
        return true;
    }

    double ata00 = 0.0, ata01 = 0.0, ata11 = 0.0;
    double atb10 = 0.0, atb11 = 0.0;
    double atb20 = 0.0, atb21 = 0.0;
    for (const auto& sample : normals) {
        ata00 += sample.a * sample.a;
        ata01 += sample.a * sample.b;
        ata11 += sample.b * sample.b;
        atb10 += sample.a * sample.dn1;
        atb11 += sample.b * sample.dn1;
        atb20 += sample.a * sample.dn2;
        atb21 += sample.b * sample.dn2;
    }

    double fit00 = ata00;
    double fit11 = ata11;
    double det = fit00 * fit11 - ata01 * ata01;
    if (std::abs(det) <= 1e-20) {
        const double ridge =
            std::max(1.0e-20, 1.0e-8 * std::max(ata00 + ata11, 1.0));
        fit00 += ridge;
        fit11 += ridge;
        det = fit00 * fit11 - ata01 * ata01;
    }
    if (std::abs(det) <= 1e-30) return false;

    const double b11 = (fit11 * atb10 - ata01 * atb11) / det;
    const double b12 = (fit00 * atb11 - ata01 * atb10) / det;
    const double b21 = (fit11 * atb20 - ata01 * atb21) / det;
    const double b22 = (fit00 * atb21 - ata01 * atb20) / det;

    curvature.k11 = b11;
    curvature.k12 = 0.5 * (b12 + b21);
    curvature.k22 = b22;
    return true;
}

} // namespace SurfaceCurvature
} // namespace Math
} // namespace SF
