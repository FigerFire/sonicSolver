#pragma once

/// @file SF_newtonian.h
/// @brief 牛顿流体应力的无状态本构算子。

#include "methods/math/SF_tensor.h"

namespace SF::Constitutive::Newtonian {

/// @brief Stokes 假设下的可压缩牛顿流体黏性应力。
/// @param dynamicViscosity 动力黏度。
/// @param velocityGradient 速度梯度，第一下标是速度分量。
/// @return `mu*(gradU+gradU^T)-2/3*mu*div(U)*I`。
inline SymmTensor3 stress(double dynamicViscosity,
                         const Tensor3& velocityGradient) {
    SymmTensor3 result = Math::twoSymm(velocityGradient);
    const double dilation = (2.0 / 3.0) * Math::trace(velocityGradient);
    result.xx = dynamicViscosity * (result.xx - dilation);
    result.yy = dynamicViscosity * (result.yy - dilation);
    result.zz = dynamicViscosity * (result.zz - dilation);
    result.xy *= dynamicViscosity;
    result.xz *= dynamicViscosity;
    result.yz *= dynamicViscosity;
    return result;
}

/// @brief 对称应力张量作用于面积向量后的牵引力。
inline Vector3 traction(const SymmTensor3& stressTensor,
                        const Vector3& area) {
    return {
        area.x * stressTensor.xx
            + area.y * stressTensor.xy
            + area.z * stressTensor.xz,
        area.x * stressTensor.xy
            + area.y * stressTensor.yy
            + area.z * stressTensor.yz,
        area.x * stressTensor.xz
            + area.y * stressTensor.yz
            + area.z * stressTensor.zz};
}

} // namespace SF::Constitutive::Newtonian
