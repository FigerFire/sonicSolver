/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_thermodynamicClosure.h"

/// @file SF_boundaryClosure.h
/// @brief ILW 物理边界高阶闭合工具。
///
/// 该文件只服务普通物理边界 ghost 层。IBM 曲面 ghost-cell 闭合由
/// `src/models/ibm/method/ghost/SF_ilwClosure.*` 中的局部几何预处理和 Taylor 构造负责。
///
/// 提供多维张量积 WENO 系数构造：
/// 2D 时收集 N×N 模板，3D 时收集 N×N×N 模板（N = ILW 阶数），
/// 对张量积多项式做最小二乘拟合，再以 WENO 型非线性权重
/// 合成法向 Taylor 系数。样本不足或拟合失败时 fail-fast。

#include "SF_ILW.h"
#include "core/mesh/SF_dimension.h"
#include "methods/math/discrete/SF_polynomial.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace SF {
namespace Boundary {
namespace ILW {
namespace BoundaryClosure {

/// @brief 支持的最高 ILW 精度阶。
constexpr int kMaxAccuracyOrder = 9;

/// @brief 支持的最高 Taylor 法向导数阶。
constexpr int kMaxTaylorOrder = 8;

/// @brief 物理边界法向信息。
struct BoundaryNormal {
    /// @brief 法向轴，0/1/2 对应 x/y/z。
    int axis = -1;
    /// @brief 外法向索引符号，-1 为低端边界，+1 为高端边界。
    int sign = 0;
};

/// @brief 普通物理边界的活跃切向方向。
struct TangentialAxes {
    /// @brief 第一切向轴。
    int t1 = -1;
    /// @brief 第二切向轴；2D时为-1。
    int t2 = -1;
    /// @brief true表示存在第二切向轴，需要使用三变量张量拟合。
    bool useT2 = false;
};

/// @brief 根据当前求解维度和边界法向选择活跃切向轴。
inline TangentialAxes activeTangentialAxes(const BoundaryNormal& normal) {
    TangentialAxes axes;
    for (int axis = 0; axis < 3; ++axis) {
        if (axis == normal.axis) continue;
        if (!Math::isDirectionActiveIndex(axis)) continue;
        if (axes.t1 < 0) {
            axes.t1 = axis;
        } else if (axes.t2 < 0) {
            axes.t2 = axis;
        }
    }
    axes.useT2 = axes.t2 >= 0;
    if (axes.t1 < 0) {
        std::cerr << "[SF FATAL] ILW boundary normal axis "
                  << normal.axis
                  << " has no active tangential direction. At least two "
                  << "active directions are required for physical-boundary ILW."
                  << std::endl;
        std::exit(1);
    }
    return axes;
}

/// @brief 校验并返回当前边界策略指定的 ILW 精度阶。
/// @return 3/5/7/9。
int requireAccuracyOrder(int accuracyOrder);

/// @brief 返回指定 ILW 精度对应的 Taylor 最高导数阶。
/// @return 2/4/6/8。
int taylorOrder(int accuracyOrder);

/// @brief 判断边界点或 ghost 点所在的活跃法向。
/// @param field 结构网格场。
/// @param i i 索引。
/// @param j j 索引。
/// @param k k 索引。
/// @return 法向信息；无法定位时 sign=0。
BoundaryNormal boundaryNormal(const Field& field, int i, int j, int k);

/// @brief 取法向局部网格间距。
/// @param field 结构网格场。
/// @param i 边界侧内点 i。
/// @param j 边界侧内点 j。
/// @param k 边界侧内点 k。
/// @param normal 边界法向信息。
/// @return 正法向网格距。
double normalSpacing(const Field& field,
                     int i, int j, int k,
                     const BoundaryNormal& normal);

/// @brief 返回 ghost 层号。
/// @param gi ghost 点 i。
/// @param gj ghost 点 j。
/// @param gk ghost 点 k。
/// @param ri 最近边界内点 i。
/// @param rj 最近边界内点 j。
/// @param rk 最近边界内点 k。
/// @param normal 边界法向信息。
/// @return 1 为紧邻边界的 ghost 层；无法判断时返回 0。
int ghostLayer(int gi, int gj, int gk,
               int ri, int rj, int rk,
               const BoundaryNormal& normal);

/// @brief 从边界面到 ghost 中心的正距离。
/// @param layer ghost 层号。
/// @param h 法向网格距。
/// @return ghost 中心的局部法向坐标。
double ghostDistance(int layer, double h);

/// @brief 多维张量积 WENO 外推，从样本值和局部坐标计算法向 Taylor 系数。
///
/// 样本按内点优先（s 越小越靠内）排序，每个样本包含 `(value, s, a, b)`。
/// `useT2 == true` 时使用 `(s,a,b)` 三变量拟合；否则只用 `(s,a)` 平面。
///
/// @param samples 样本值及局部 `(s,a,b)` 坐标。
/// @param accuracyOrder ILW 精度阶，3/5/7/9。
/// @param useT2 是否启用第二切向坐标。
/// @param coefficients 输出 Taylor 系数 `U^(m)(0)/m!`。
/// @return 拟合成功时返回 true。
bool extrapolateMultiDim(
    const std::vector<std::pair<double, Math::Polynomial::TensorPoint>>& samples,
    int accuracyOrder,
    bool useT2,
    std::array<double, kMaxAccuracyOrder>& coefficients);

/// @brief 计算 Taylor 系数在 ghost 点的值。
/// @param coefficients `U^(m)(0)/m!` 系数。
/// @param taylorOrder 使用的最高 Taylor 阶。
/// @param distance ghost 点局部法向坐标。
/// @return 外推值。
double evaluateTaylor(const std::array<double, kMaxAccuracyOrder>& coefficients,
                      int taylorOrder,
                      double distance);

/// @brief ILW 物理边界闭合失败时终止。
/// @param bcName 边界类型名。
/// @param i i 索引。
/// @param j j 索引。
/// @param k k 索引。
/// @param reason 失败原因。
void fatalClosure(const char* bcName,
                  int i, int j, int k,
                  const std::string& reason);

/// @brief 沿物理边界收集多维内点样本。
///
/// 模板在法向、第一切向和第二切向（若 3D）上分别收集 `accuracyOrder` 个点，
/// 构造大小为 N×N（2D）或 N×N×N（3D）的规则模板。
/// 每个样本携带其在边界局部坐标系 `(s,a,b)` 中的坐标和标量值。
///
/// @param field 结构网格场。
/// @param i 边界侧内点 i。
/// @param j 边界侧内点 j。
/// @param k 边界侧内点 k。
/// @param normal 边界法向信息。
/// @param accuracyOrder ILW 精度阶，3/5/7/9。
/// @param getter 标量取值函数，签名为 `(field,i,j,k)`。
/// @param h 法向网格间距（用于构造局部坐标）。
/// @param samples 输出样本 `(value, {s,a,b})`。
/// @return 实际收集到的有效样本数。
template <typename Getter>
int collectMultiDimSamples(
    const Field& field,
    int i, int j, int k,
    const BoundaryNormal& normal,
    int accuracyOrder,
    Getter&& getter,
    double h,
    std::vector<std::pair<double, Math::Polynomial::TensorPoint>>& samples) {
    samples.clear();
    const int N = accuracyOrder;
    const int half = N / 2;

    /// 确定活跃切向轴。2D empty 问题只保留一个切向方向。
    const TangentialAxes tangential = activeTangentialAxes(normal);
    const int t1 = tangential.t1;
    const int t2 = tangential.t2;
    const bool useT2 = tangential.useT2;

    auto addSample = [&](int pi, int pj, int pk,
                         double s, double a, double b) {
        const double value = getter(field, pi, pj, pk);
        if (!std::isfinite(value) ||
            !std::isfinite(s) ||
            !std::isfinite(a) ||
            !std::isfinite(b)) {
            fatalClosure("BoundaryClosure", pi, pj, pk,
                         "non-finite sample in multi-dimensional ILW stencil.");
        }
        samples.push_back({value, {s, a, b}});
    };

    /// 计算边界处各方向的步长。
    auto step = [&](int axisIndex, int baseI, int baseJ, int baseK,
                    int& ni, int& nj, int& nk) -> double {
        ni = baseI;
        nj = baseJ;
        nk = baseK;
        double dp[3] = {field.X(baseI, baseJ, baseK),
                         field.Y(baseI, baseJ, baseK),
                         field.Z(baseI, baseJ, baseK)};
        if (axisIndex == 0) {
            ni = baseI + 1;
        } else if (axisIndex == 1) {
            nj = baseJ + 1;
        } else {
            nk = baseK + 1;
        }
        if (!Boundary::isPhysicalCell(field, ni, nj, nk)) return 0.0;
        double dn[3] = {field.X(ni, nj, nk),
                         field.Y(ni, nj, nk),
                         field.Z(ni, nj, nk)};
        double dx = dn[0] - dp[0];
        double dy = dn[1] - dp[1];
        double dz = dn[2] - dp[2];
        double mag = std::sqrt(dx * dx + dy * dy + dz * dz);
        return (mag > 1e-14) ? mag : 0.0;
    };

    /// 计算三个方向的步长。
    double hT1 = h;
    double hT2 = h;
    {
        int dummyNi, dummyNj, dummyNk;
        double hTemp = step(t1, i, j, k, dummyNi, dummyNj, dummyNk);
        if (hTemp > 1e-14) hT1 = hTemp;
        if (useT2) {
            hTemp = step(t2, i, j, k, dummyNi, dummyNj, dummyNk);
            if (hTemp > 1e-14) hT2 = hTemp;
        }
    }

    samples.reserve(useT2 ? (size_t)N * N * N : (size_t)N * N);

    /// 沿法向向内步进。
    for (int ns = 0; ns < N; ++ns) {
        int si = i;
        int sj = j;
        int sk = k;
        if (normal.axis == 0) si -= normal.sign * ns;
        else if (normal.axis == 1) sj -= normal.sign * ns;
        else sk -= normal.sign * ns;

        if (!Boundary::isPhysicalCell(field, si, sj, sk)) continue;

        const double s = -(ns + 0.5) * h;  ///< 负值表示域内。

        /// 沿第一切向偏移。
        for (int na = -half; na <= half; ++na) {
            int ai = si;
            int aj = sj;
            int ak = sk;
            if (t1 == 0) ai += na;
            else if (t1 == 1) aj += na;
            else ak += na;

            if (!Boundary::isPhysicalCell(field, ai, aj, ak)) continue;

            if (!useT2) {
                const double a = na * hT1;
                addSample(ai, aj, ak, s, a, 0.0);
                continue;
            }

            /// 沿第二切向偏移。
            for (int nb = -half; nb <= half; ++nb) {
                int bi = ai;
                int bj = aj;
                int bk = ak;
                if (t2 == 0) bi += nb;
                else if (t2 == 1) bj += nb;
                else bk += nb;

                if (!Boundary::isPhysicalCell(field, bi, bj, bk)) continue;

                const double a = na * hT1;
                const double b = nb * hT2;
                addSample(bi, bj, bk, s, a, b);
            }
        }
    }

    return (int)samples.size();
}

/// @brief 用多维张量积 WENO 构造一个标量变量的法向 Taylor 系数。
///
/// 对 2D 问题收集 N×N 模板，3D 问题收集 N×N×N 模板，
/// 用张量积多项式最小二乘拟合和多维 WENO 型非线性权重
/// 合成法向 Taylor 系数。样本不足或拟合失败时直接报错，不做隐式降维兜底。
///
/// @param field 结构网格场。
/// @param i 边界侧内点 i。
/// @param j 边界侧内点 j。
/// @param k 边界侧内点 k。
/// @param normal 边界法向信息。
/// @param getter 标量取值函数。
/// @param bcName 边界类型名。
/// @param coefficients 输出 Taylor 系数。
template <typename Getter>
void buildMultiDimScalarTaylorCoefficients(
        const Field& field,
        int i, int j, int k,
        const BoundaryNormal& normal,
        int accuracyOrder,
        Getter&& getter,
        const char* bcName,
        std::array<double, kMaxAccuracyOrder>& coefficients) {
    const int accuracy = requireAccuracyOrder(accuracyOrder);
    const double h = normalSpacing(field, i, j, k, normal);

    std::vector<std::pair<double, Math::Polynomial::TensorPoint>> samples;
    const int collected = collectMultiDimSamples(
        field, i, j, k, normal, accuracy,
        std::forward<Getter>(getter), h, samples);

    const bool useT2 = activeTangentialAxes(normal).useT2;
    const int minRequired = useT2
        ? Math::Polynomial::tensorBasisSize(true, accuracy - 1)
        : Math::Polynomial::tensorBasisSize(false, accuracy - 1);

    if (collected < minRequired) {
        fatalClosure(bcName, i, j, k,
                     "insufficient multi-dimensional samples for ILW"
                     + std::to_string(accuracy)
                     + " boundary extrapolation: collected "
                     + std::to_string(collected)
                     + ", required at least "
                     + std::to_string(minRequired)
                     + ". Reduce [numerics].ILW explicitly or enlarge the "
                       "valid multi-dimensional stencil.");
    }

    if (!extrapolateMultiDim(samples, accuracy, useT2, coefficients)) {
        fatalClosure(bcName, i, j, k,
                     "multi-dimensional ILW WENO-type extrapolation "
                     "failed for ILW" + std::to_string(accuracy)
                     + ". No implicit dimension fallback is applied.");
    }
}

/// @brief 默认 Field 变量取值器。
struct FieldVariableGetter {
    int vIdx = 0;

    /// @brief 读取 Field 标量变量。
    double operator()(const Field& field, int i, int j, int k) const {
        return field(i, j, k, vIdx);
    }
};

/// @brief 压力标量取值器。
struct PressureGetter {
    /// @brief 从守恒量计算压力。
    double operator()(const Field& field, int i, int j, int k) const {
        return Boundary::pressureAt(field, i, j, k);
    }
};

} // namespace BoundaryClosure
} // namespace ILW
} // namespace Boundary
} // namespace SF
