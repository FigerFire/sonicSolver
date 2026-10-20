/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_derivative.h
/// @brief 一维非均匀点列求导算子。
///
/// 该文件实现通用有限差分权重求解: 给定展开点`s=0`附近的样本
/// `(s_i, q_i)`，求`q^(m)(0)`。它不读取Field、不处理IBM拓扑；
/// ILW模块负责提供沿壁面法向/近似法向的一侧样本。

#include "SF_linearSystem.h"
#include "SF_taylor.h"

#include <array>
#include <cmath>
#include <vector>

namespace SF {
namespace Math {
namespace Derivative {

/// @brief 一维样本点。
/// @tparam NVAR 变量个数。
template <int NVAR>
struct Sample1D {
    /// @brief 样本到展开点的有符号距离。
    double s = 0.0;
    /// @brief 样本变量值。
    std::array<double, NVAR> value{};
};

/// @brief 解线性方程组 `A x = b`。
/// @param a 系数矩阵，函数内部会修改。
/// @param b 右端项，函数内部会修改。
/// @param x 解向量。
/// @return 成功求解返回true；矩阵奇异返回false。
inline bool solveLinearSystem(std::vector<std::vector<double>> a,
                              std::vector<double> b,
                              std::vector<double>& x) {
    return LinearSystem::solve(a, b, x);
}

/// @brief 计算展开点`s=0`处某阶导数的有限差分权重。
///
/// 权重满足:
/// `sum_i w_i s_i^m = delta_{m,d} d!`
///
/// @param nodes 样本点坐标`s_i`。
/// @param derivativeOrder 目标导数阶数`d`。
/// @param weights 输出权重`w_i`。
/// @return 成功返回true；样本不足或矩阵奇异返回false。
inline bool finiteDifferenceWeights(const std::vector<double>& nodes,
                                    int derivativeOrder,
                                    std::vector<double>& weights) {
    const int n = (int)nodes.size();
    if (n <= derivativeOrder || derivativeOrder < 0) return false;

    std::vector<std::vector<double>> a(n, std::vector<double>(n, 0.0));
    std::vector<double> b(n, 0.0);
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            double power = 1.0;
            for (int p = 0; p < row; ++p) power *= nodes[col];
            a[row][col] = power;
        }
        b[row] = (row == derivativeOrder) ? Taylor::factorial(derivativeOrder) : 0.0;
    }

    return solveLinearSystem(a, b, weights);
}

/// @brief 根据一维样本计算从0阶到maxDerivative阶的导数。
/// @tparam NVAR 变量个数。
/// @param samples 一维样本；通常第一个样本为边界点`s=0`。
/// @param maxDerivative 希望计算的最高导数阶数。
/// @param derivatives 输出导数数组，`derivatives[m][v]=q_v^(m)(0)`。
/// @return 至少成功得到0阶导数时返回true。
template <int NVAR>
inline bool derivativesAtOrigin(const std::vector<Sample1D<NVAR>>& samples,
                                int maxDerivative,
                                std::vector<std::array<double, NVAR>>& derivatives) {
    derivatives.clear();
    if (samples.empty()) return false;

    const int usableOrder = std::max(0, std::min(maxDerivative, (int)samples.size() - 1));
    derivatives.assign((size_t)usableOrder + 1, {});

    std::vector<double> nodes(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) nodes[i] = samples[i].s;

    for (int order = 0; order <= usableOrder; ++order) {
        std::vector<double> weights;
        if (!finiteDifferenceWeights(nodes, order, weights)) return order > 0;
        for (int v = 0; v < NVAR; ++v) {
            double value = 0.0;
            for (size_t s = 0; s < samples.size(); ++s) {
                value += weights[s] * samples[s].value[v];
            }
            derivatives[order][v] = value;
        }
    }
    return true;
}

} // namespace Derivative
} // namespace Math
} // namespace SF
