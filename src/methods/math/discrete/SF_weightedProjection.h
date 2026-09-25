/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_weightedProjection.h
/// @brief 向量值样本的线性权重代入与投影残差指标。
///
/// 本文件提供预计算投影矩阵对向量样本的纯代数代入工具。不依赖Field、
/// IO、MPI或求解器状态。

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace SF {
namespace Math {
namespace WeightedProjection {

/// @brief 对向量值样本应用一组预计算线性权重。
/// @param weights 样本线性权重。
/// @param values 样本向量值。
/// @return 权重和对应的向量值。
template <std::size_t N>
inline std::array<double, N> applyWeights(
        const std::vector<double>& weights,
        const std::vector<std::array<double, N>>& values) {
    std::array<double, N> result{};
    const int n = std::min<int>((int)weights.size(), (int)values.size());
    for (int i = 0; i < n; ++i) {
        for (std::size_t v = 0; v < N; ++v) {
            result[v] += weights[(std::size_t)i] * values[(std::size_t)i][v];
        }
    }
    return result;
}

/// @brief 用预计算投影矩阵计算多项式残差平滑指标。
/// @param projectionRows `coeff = M * values` 的投影矩阵行。
/// @param basisRows 样本点处的多项式基函数值。
/// @param values 样本向量值。
/// @return 残差与高阶系数粗糙度构成的平滑指标。
template <std::size_t N>
inline double projectedPolynomialSmoothness(
        const std::vector<std::vector<double>>& projectionRows,
        const std::vector<std::vector<double>>& basisRows,
        const std::vector<std::array<double, N>>& values) {
    if (values.empty() || projectionRows.empty() ||
        basisRows.size() != values.size()) {
        return 1.0e30;
    }

    std::vector<std::array<double, N>> coeff(projectionRows.size());
    for (std::size_t c = 0; c < projectionRows.size(); ++c) {
        coeff[c] = applyWeights<N>(projectionRows[c], values);
    }

    double residual = 0.0;
    for (std::size_t sample = 0; sample < values.size(); ++sample) {
        if (basisRows[sample].size() != coeff.size()) return 1.0e30;
        for (std::size_t v = 0; v < N; ++v) {
            double predicted = 0.0;
            for (std::size_t c = 0; c < coeff.size(); ++c) {
                predicted += coeff[c][v] * basisRows[sample][c];
            }
            const double diff = predicted - values[sample][v];
            residual += diff * diff;
        }
    }

    double roughness = 0.0;
    for (std::size_t c = 1; c < coeff.size(); ++c) {
        for (std::size_t v = 0; v < N; ++v) {
            roughness += coeff[c][v] * coeff[c][v];
        }
    }
    return residual / std::max<std::size_t>(1, values.size())
         + 1.0e-8 * roughness;
}

} // namespace WeightedProjection
} // namespace Math
} // namespace SF
