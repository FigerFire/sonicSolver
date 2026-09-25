/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_polynomial.h
/// @brief 局部多项式基与最小二乘投影矩阵。
///
/// 本文件只处理纯数学对象：样本局部坐标、多项式基、正规方程和
/// `coeff = M * values` 投影矩阵。它不依赖Field、不读取IBM拓扑，
/// 供ILW/WENO-type外推在前处理阶段预计算几何矩阵。

#include "SF_derivative.h"
#include "SF_taylor.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace SF {
namespace Math {
namespace Polynomial {

/// @brief 局部张量多项式样本点。
struct TensorPoint {
    /// @brief 法向局部坐标。
    double s = 0.0;
    /// @brief 第一切向局部坐标。
    double a = 0.0;
    /// @brief 第二切向局部坐标；二维退化时可忽略。
    double b = 0.0;
};

/// @brief 最小二乘投影预计算计划。
struct LeastSquaresPlan {
    /// @brief 是否使用第二切向坐标。
    bool useT2 = false;
    /// @brief 多项式总阶数。
    int order = 0;
    /// @brief 每个样本点对应的多项式基函数值。
    std::vector<std::vector<double>> basisRows;
    /// @brief `projectionRows[c][i]` 表示第c个系数对第i个样本值的线性权重。
    std::vector<std::vector<double>> projectionRows;
};

/// @brief 整数幂。
inline double powInt(double x, int p) {
    double value = 1.0;
    for (int i = 0; i < p; ++i) value *= x;
    return value;
}

/// @brief 返回总阶不超过order的局部张量多项式基函数数量。
inline int tensorBasisSize(bool useT2, int order) {
    int count = 0;
    for (int total = 0; total <= order; ++total) {
        for (int ps = total; ps >= 0; --ps) {
            count += useT2 ? (total - ps + 1) : 1;
        }
    }
    return count;
}

/// @brief 计算局部张量多项式基函数值。
///
/// 基函数包含阶乘分母，因此纯法向项 `s^m/m!` 的系数直接对应
/// m阶法向导数。
inline std::vector<double> tensorBasis(double s, double a, double b,
                                       bool useT2, int order) {
    std::vector<double> phi;
    phi.reserve((size_t)tensorBasisSize(useT2, order));

    for (int total = 0; total <= order; ++total) {
        for (int ps = total; ps >= 0; --ps) {
            if (useT2) {
                for (int pa = total - ps; pa >= 0; --pa) {
                    const int pb = total - ps - pa;
                    const double denom = Taylor::factorial(ps)
                                       * Taylor::factorial(pa)
                                       * Taylor::factorial(pb);
                    phi.push_back(powInt(s, ps) * powInt(a, pa) * powInt(b, pb)
                                  / denom);
                }
            } else {
                const int pa = total - ps;
                const double denom = Taylor::factorial(ps) * Taylor::factorial(pa);
                phi.push_back(powInt(s, ps) * powInt(a, pa) / denom);
            }
        }
    }
    return phi;
}

/// @brief 返回纯法向导数对应的多项式系数位置。
inline int normalDerivativeIndex(bool useT2, int order, int derivativeOrder) {
    int idx = 0;
    for (int total = 0; total <= order; ++total) {
        for (int ps = total; ps >= 0; --ps) {
            if (useT2) {
                for (int pa = total - ps; pa >= 0; --pa) {
                    const int pb = total - ps - pa;
                    if (ps == derivativeOrder && pa == 0 && pb == 0) return idx;
                    ++idx;
                }
            } else {
                const int pa = total - ps;
                if (ps == derivativeOrder && pa == 0) return idx;
                ++idx;
            }
        }
    }
    return -1;
}

/// @brief 对正规方程矩阵一次性求逆。
///
/// `buildLeastSquaresPlan` 需要 `(A^T A)^{-1}` 的每一行。逐个右端项求解
/// 会重复消元；这里用 Gauss-Jordan 对增广矩阵 `[A^T A | I]` 一次完成，
/// 若主元退化则返回 false，让调用方按 ILW fail-fast 诊断处理。
inline bool invertNormalMatrix(const std::vector<std::vector<double>>& matrix,
                               std::vector<std::vector<double>>& inverse) {
    const int n = (int)matrix.size();
    if (n <= 0) return false;
    for (const auto& row : matrix) {
        if ((int)row.size() != n) return false;
    }

    double normMax = 0.0;
    std::vector<std::vector<double>> aug(
        (size_t)n, std::vector<double>((size_t)(2 * n), 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double value = matrix[(size_t)i][(size_t)j];
            aug[(size_t)i][(size_t)j] = value;
            normMax = std::max(normMax, std::abs(value));
        }
        aug[(size_t)i][(size_t)(n + i)] = 1.0;
    }
    (void)normMax;
    constexpr double pivotTol = 1.0e-14;

    for (int col = 0; col < n; ++col) {
        int pivotRow = col;
        double pivotAbs = std::abs(aug[(size_t)col][(size_t)col]);
        for (int row = col + 1; row < n; ++row) {
            const double candidate = std::abs(aug[(size_t)row][(size_t)col]);
            if (candidate > pivotAbs) {
                pivotAbs = candidate;
                pivotRow = row;
            }
        }
        if (pivotAbs <= pivotTol) return false;
        if (pivotRow != col) {
            std::swap(aug[(size_t)pivotRow], aug[(size_t)col]);
        }

        const double pivot = aug[(size_t)col][(size_t)col];
        for (int j = 0; j < 2 * n; ++j) {
            aug[(size_t)col][(size_t)j] /= pivot;
        }
        for (int row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = aug[(size_t)row][(size_t)col];
            if (std::abs(factor) < 1.0e-30) continue;
            for (int j = 0; j < 2 * n; ++j) {
                aug[(size_t)row][(size_t)j] -=
                    factor * aug[(size_t)col][(size_t)j];
            }
        }
    }

    inverse.assign((size_t)n, std::vector<double>((size_t)n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse[(size_t)i][(size_t)j] =
                aug[(size_t)i][(size_t)(n + j)];
        }
    }
    return true;
}

/// @brief 预计算最小二乘投影矩阵。
///
/// 给定样本点局部坐标，构造 `A` 和 `M=(A^T A)^{-1}A^T`。
/// 后续时间推进中只需执行 `coeff = M * value`。
///
/// @param points 样本局部坐标。
/// @param useT2 是否启用第二切向坐标。
/// @param order 多项式总阶数。
/// @param plan 输出投影计划。
/// @return 样本不足或正规方程奇异时返回false。
inline bool buildLeastSquaresPlan(const std::vector<TensorPoint>& points,
                                  bool useT2,
                                  int order,
                                  LeastSquaresPlan& plan) {
    const int nBasis = tensorBasisSize(useT2, order);
    const int nSample = (int)points.size();
    if (nSample < nBasis) return false;

    plan.useT2 = useT2;
    plan.order = order;
    plan.basisRows.assign((size_t)nSample, {});
    for (int i = 0; i < nSample; ++i) {
        plan.basisRows[(size_t)i] = tensorBasis(points[(size_t)i].s,
                                                points[(size_t)i].a,
                                                points[(size_t)i].b,
                                                useT2, order);
    }

    std::vector<std::vector<double>> ata(
        (size_t)nBasis, std::vector<double>((size_t)nBasis, 0.0));
    for (const auto& phi : plan.basisRows) {
        for (int r = 0; r < nBasis; ++r) {
            for (int c = 0; c < nBasis; ++c) {
                ata[(size_t)r][(size_t)c] += phi[(size_t)r] * phi[(size_t)c];
            }
        }
    }

    std::vector<std::vector<double>> invAta;
    if (!invertNormalMatrix(ata, invAta)) return false;

    plan.projectionRows.assign(
        (size_t)nBasis, std::vector<double>((size_t)nSample, 0.0));
    for (int coeff = 0; coeff < nBasis; ++coeff) {
        for (int sample = 0; sample < nSample; ++sample) {
            double weight = 0.0;
            for (int b = 0; b < nBasis; ++b) {
                weight += invAta[(size_t)coeff][(size_t)b]
                        * plan.basisRows[(size_t)sample][(size_t)b];
            }
            plan.projectionRows[(size_t)coeff][(size_t)sample] = weight;
        }
    }

    return true;
}

} // namespace Polynomial
} // namespace Math
} // namespace SF
