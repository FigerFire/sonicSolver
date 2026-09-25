/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_linearSystem.h
/// @brief 小型稠密线性方程组求解器。
///
/// 本文件只包含纯数学求解器，不依赖Field、不读取全局配置。
/// 主要供局部差分权重、最小二乘投影和ILW边界闭合的小矩阵使用。

#include <cmath>
#include <vector>

namespace SF {
namespace Math {
namespace LinearSystem {

/// @brief 解线性方程组 `A x = b`。
/// @param a 系数矩阵，函数内部会修改副本。
/// @param b 右端项，函数内部会修改副本。
/// @param x 输出解向量。
/// @return 成功求解返回true；矩阵奇异或维度不匹配返回false。
inline bool solve(std::vector<std::vector<double>> a,
                  std::vector<double> b,
                  std::vector<double>& x) {
    const int n = (int)b.size();
    if (n <= 0 || (int)a.size() != n) return false;
    for (const auto& row : a) {
        if ((int)row.size() != n) return false;
    }

    x.assign((size_t)n, 0.0);
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(a[(size_t)row][(size_t)col])
                > std::abs(a[(size_t)pivot][(size_t)col])) {
                pivot = row;
            }
        }
        if (std::abs(a[(size_t)pivot][(size_t)col]) < 1.0e-14) {
            return false;
        }
        if (pivot != col) {
            std::swap(a[(size_t)pivot], a[(size_t)col]);
            std::swap(b[(size_t)pivot], b[(size_t)col]);
        }

        const double invPivot = 1.0 / a[(size_t)col][(size_t)col];
        for (int j = col; j < n; ++j) {
            a[(size_t)col][(size_t)j] *= invPivot;
        }
        b[(size_t)col] *= invPivot;

        for (int row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = a[(size_t)row][(size_t)col];
            if (std::abs(factor) < 1.0e-30) continue;
            for (int j = col; j < n; ++j) {
                a[(size_t)row][(size_t)j] -= factor * a[(size_t)col][(size_t)j];
            }
            b[(size_t)row] -= factor * b[(size_t)col];
        }
    }

    x = b;
    return true;
}

} // namespace LinearSystem
} // namespace Math
} // namespace SF
