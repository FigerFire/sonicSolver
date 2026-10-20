/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_hjWeno.h
/// @brief Level Set Hamilton-Jacobi WENO一侧导数核。

#include "SF_field.h"
#include "SF_state.h"

#include <limits>
#include <vector>

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief 一个计算坐标方向上的左右一侧导数。
struct OneSidedDerivative {
    double minus = 0.0; ///< D^- phi, 从负向模板重构。
    double plus = 0.0;  ///< D^+ phi, 从正向模板重构。
};

/// @brief 三个计算坐标方向上的一侧导数。
struct OneSidedGradient {
    OneSidedDerivative xi;
    OneSidedDerivative eta;
    OneSidedDerivative zeta;
};

/// @brief 用户显式给出的 HJ-WENO 非线性权重参数。
struct HJWenoWeightOptions {
    double epsilon = std::numeric_limits<double>::quiet_NaN(); ///< 正则量，必须 > 0。
    double power = std::numeric_limits<double>::quiet_NaN();   ///< 权重指数，必须 > 0。
};

/// @brief HJ-WENO公共工具。
class HJWeno {
public:
    /// @brief 校验WENO阶数。
    /// @param order 允许3/5/7。
    /// @param context 错误消息上下文。
    static void validateOrder(int order, const char* context);

    /// @brief 给定阶数需要的ghost层数。
    /// @param order WENO阶数。
    /// @return WENO3/5/7分别需要2/3/4层。
    static int requiredGhostLayers(int order);

    /// @brief 检查Field是否满足某阶HJ-WENO模板。
    /// @param field 参考网格。
    /// @param order WENO阶数。
    /// @param context 错误消息上下文。
    static void requireStencil(const Field& field,
                               int order,
                               const char* context);

    /// @brief 判断某计算坐标方向是否参与Level Set离散。
    /// @param field 参考网格。
    /// @param axis 0=xi, 1=eta, 2=zeta。
    static bool activeAxis(const Field& field, int axis);

    /// @brief 计算一个方向的一侧导数。
    /// @param levelSet Level Set场。
    /// @param values phi快照。
    /// @param i,j,k 存储索引。
    /// @param axis 0=xi, 1=eta, 2=zeta。
    /// @param order WENO阶数。
    /// @param weights 用户显式输入的非线性权重参数。
    static OneSidedDerivative derivative(const LevelSetField& levelSet,
                                         const std::vector<double>& values,
                                         int i, int j, int k,
                                         int axis,
                                         int order,
                                         const HJWenoWeightOptions& weights);

    /// @brief 按速度符号选取上风导数。
    /// @param d 一侧导数。
    /// @param speed 对应计算坐标方向的contravariant速度。
    static double upwind(const OneSidedDerivative& d, double speed);

    /// @brief 速度在某计算坐标方向上的contravariant系数。
    static double contravariantSpeed(const Field& field,
                                     int i, int j, int k,
                                     int axis,
                                     const Vector3& velocity);

    /// @brief 将计算坐标梯度转换成物理空间梯度。
    static Vector3 physicalGradient(const Field& field,
                                    int i, int j, int k,
                                    double dXi,
                                    double dEta,
                                    double dZeta);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
