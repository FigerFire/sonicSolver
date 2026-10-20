/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_taylor.h
/// @brief 一维Taylor展开算子。
///
/// 该文件只负责数学展开，不包含任何IBM、边界条件或物理模型。
/// ILW模块把边界状态及其法向导数传入这里，再得到ghost侧外推值。

#include <algorithm>
#include <array>
#include <vector>

namespace SF {
namespace Math {
namespace Taylor {

/// @brief 计算阶乘。
/// @param n 非负整数。
/// @return n!。
inline double factorial(int n) {
    double value = 1.0;
    for (int i = 2; i <= n; ++i) value *= (double)i;
    return value;
}

/// @brief 将用户请求阶数截断到已有导数表支持的最高阶。
/// @param requestedOrder 用户请求的最高Taylor阶。
/// @param derivativeCount 已有导数数量。
/// @return 实际可使用的最高阶；导数表为空时返回-1。
inline int usableOrder(int requestedOrder, int derivativeCount) {
    if (derivativeCount <= 0) return -1;
    const int requested = std::max(0, requestedOrder);
    return std::min(requested, derivativeCount - 1);
}

/// @brief 顺序构造一维Taylor导数表。
///
/// `builder(order, built, derivative)` 会按 `order=0..requestedOrder`
/// 依次调用。`built` 保存已经构造好的低阶导数，因此高阶构造可以自然依赖
/// 0阶壁面状态或1阶边界方程结果。
///
/// @tparam NVAR 变量个数。
/// @tparam Builder 导数构造器类型。
/// @param requestedOrder 需要构造的最高Taylor阶。
/// @param builder 顺序导数构造器。
/// @param derivatives 输出导数表，`derivatives[m][v] = q_v^(m)(0)`。
/// @return 全部目标阶构造成功时返回true。
template <int NVAR, typename Builder>
inline bool buildDerivativeTable(
        int requestedOrder,
        Builder&& builder,
        std::vector<std::array<double, NVAR>>& derivatives) {
    derivatives.clear();
    const int targetOrder = std::max(0, requestedOrder);
    derivatives.reserve((size_t)targetOrder + 1);

    for (int order = 0; order <= targetOrder; ++order) {
        std::array<double, NVAR> derivative{};
        if (!builder(order, derivatives, derivative)) return false;
        derivatives.push_back(derivative);
    }
    return true;
}

/// @brief 对N个变量做一维Taylor展开。
///
/// 展开形式:
/// `q(s) = q(0) + s q'(0) + s^2/2 q''(0) + ...`
///
/// @tparam NVAR 变量个数。
/// @param derivatives `derivatives[m][v]`表示第v个变量的m阶导数。
/// @param displacement 从展开点到目标点的有符号距离。
/// @param maxOrder 使用的最高导数阶数；会自动截断到derivatives可用范围。
/// @return 目标点处的外推值。
template <int NVAR>
inline std::array<double, NVAR>
evaluate1D(const std::vector<std::array<double, NVAR>>& derivatives,
           double displacement,
           int maxOrder) {
    std::array<double, NVAR> value{};
    const int order = usableOrder(maxOrder, (int)derivatives.size());
    if (order < 0) return value;

    double power = 1.0;
    for (int m = 0; m <= order; ++m) {
        const double coeff = power / factorial(m);
        for (int v = 0; v < NVAR; ++v) {
            value[v] += coeff * derivatives[m][v];
        }
        power *= displacement;
    }
    return value;
}

/// @brief 对已含阶乘分母的Taylor系数做一维展开。
///
/// `coefficients[m] = q^(m)(0)/m!`，因此展开只需累乘位移。
/// @tparam NCOEFF 系数数组长度。
/// @param coefficients Taylor系数。
/// @param displacement 从展开点到目标点的有符号距离。
/// @param maxOrder 使用的最高阶数。
/// @return 目标点外推值。
template <size_t NCOEFF>
inline double evaluateScaled1D(const std::array<double, NCOEFF>& coefficients,
                               double displacement,
                               int maxOrder) {
    const int order = usableOrder(maxOrder, (int)coefficients.size());
    if (order < 0) return 0.0;

    double value = coefficients[0];
    double power = 1.0;
    for (int m = 1; m <= order; ++m) {
        power *= displacement;
        value += power * coefficients[(size_t)m];
    }
    return value;
}

} // namespace Taylor
} // namespace Math
} // namespace SF
