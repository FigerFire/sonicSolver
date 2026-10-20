/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_face.h
/// @brief 对流通量面几何、状态和IBM邻接工具。
///
/// 该文件只提供无状态工具函数，不拥有 Field，也不读取全局配置。
/// Roe、LW、Lax-Friedrichs、Steger-Warming和WENO装配通过这些函数
/// 共享面法向、面面积、左右状态和IBM邻接判定。

#include "SF_cellType.h"
#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_utility.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace Flux {

/// @brief 半网格面上的几何量。
struct FaceGeometry {
    /// @brief 面余因子向量及面中心逆Jacobian [J*xi_x,J*xi_y,J*xi_z,1/J]。
    double metrics[4] = {0.0, 0.0, 0.0, 1.0};
    /// @brief 物理空间单位法向，方向从左单元指向右单元。
    double normal[3] = {1.0, 0.0, 0.0};
    /// @brief 面面积权重，用于把法向物理通量转成离散面通量。
    double area = 1.0;
};

/// @brief 从Field拷贝一个守恒状态。
/// @param field 结构网格场。
/// @param i i方向全局存储索引。
/// @param j j方向全局存储索引。
/// @param k k方向全局存储索引。
/// @param q 输出守恒量 [rho,rho*u,rho*v,rho*w,E]。
inline void loadConservative(const Field& field, int i, int j, int k, double q[5]) {
    q[0] = field(i, j, k, RHO);
    q[1] = field(i, j, k, RU);
    q[2] = field(i, j, k, RV);
    q[3] = field(i, j, k, RW);
    q[4] = field(i, j, k, E);
}

/// @brief 计算半网格面的单位法向和面积权重。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @return FaceGeometry，包含单位法向和面积权重。
inline FaceGeometry makeFaceGeometry(const Field& field, int i, int j, int k, Math::Dir d) {
    FaceGeometry geom;
    Math::faceMetrics(field, i, j, k, d, geom.metrics);

    const double sx = geom.metrics[0];
    const double sy = geom.metrics[1];
    const double sz = geom.metrics[2];
    const double smag = std::sqrt(sx * sx + sy * sy + sz * sz);

    if (!std::isfinite(smag) || smag <= 1.0e-300) {
        std::cerr << "[SF FATAL] Degenerate face cofactor at face ("
                  << i << "," << j << "," << k << "), dir=" << d
                  << ", area=" << smag << std::endl;
        std::exit(1);
    }

    // 面积不能加到单位法向的分母。对小尺度曲线网格，这会把
    // p*n*area 改成 p*C*area/(area+eps)，从而破坏离散度规恒等式和
    // 均匀自由流。退化面应显式报错，正常面必须使用精确归一化。
    geom.normal[0] = sx / smag;
    geom.normal[1] = sy / smag;
    geom.normal[2] = sz / smag;
    geom.area = smag;
    return geom;
}

/// @brief 计算两个相邻点沿面法向的距离。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @param normal 面单位法向。
/// @return 沿法向投影距离；退化时返回两点欧氏距离。
inline double faceSpacing(const Field& field, int i, int j, int k,
                          Math::Dir d, const double normal[3]) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    const int ir = i + di;
    const int jr = j + dj;
    const int kr = k + dk;

    const double dx = field.X(ir, jr, kr) - field.X(i, j, k);
    const double dy = field.Y(ir, jr, kr) - field.Y(i, j, k);
    const double dz = field.Z(ir, jr, kr) - field.Z(i, j, k);
    const double projected = std::abs(dx * normal[0] + dy * normal[1] + dz * normal[2]);
    const double euclidean = std::sqrt(dx * dx + dy * dy + dz * dz);
    return std::max(projected > 1e-14 ? projected : euclidean, 1e-14);
}

/// @brief 计算Euler方程法向物理通量。
/// @param q 守恒状态 [rho,rho*u,rho*v,rho*w,E]。
/// @param normal 面单位法向。
/// @param flux 输出未乘面积的法向物理通量。
inline void physicalEulerFlux(const double q[5], const double normal[3],
                              double gamma, double flux[5]) {
    const double rho = q[0];
    const double u = q[1] / rho;
    const double v = q[2] / rho;
    const double w = q[3] / rho;
    const double p = Numerics::requirePhysicalState(
        "physicalEulerFlux", q[0], q[1], q[2], q[3], q[4], gamma);
    const double un = u * normal[0] + v * normal[1] + w * normal[2];

    flux[0] = rho * un;
    flux[1] = rho * u * un + p * normal[0];
    flux[2] = rho * v * un + p * normal[1];
    flux[3] = rho * w * un + p * normal[2];
    flux[4] = (q[4] + p) * un;
}

/// @brief 判断一个WENO模板是否跨入IBM非流体区域。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @param offsets WENO模板偏移数组。
/// @param nStencil 模板点数。
/// @return 模板包含IBM ghost/solid点时返回true。
inline bool stencilTouchesIBM(const Field& field,
                              int i, int j, int k, Math::Dir d,
                              const int* offsets, int nStencil) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    for (int s = 0; s < nStencil; ++s) {
        const int ii = i + offsets[s] * di;
        const int jj = j + offsets[s] * dj;
        const int kk = k + offsets[s] * dk;
        if (SF::IBM::isIbmNonFluidCell(field, ii, jj, kk)) return true;
    }
    return false;
}

/// @brief 判断IBM相交模板是否已全部由可读ghost状态闭合。
/// @return 非流体模板点全部为IBM_GHOST_CELL时返回true；进入SOLID_CELL返回false。
inline bool stencilClosedByIBMGhosts(const Field& field,
                                     int i, int j, int k, Math::Dir d,
                                     const int* offsets, int nStencil) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    for (int s = 0; s < nStencil; ++s) {
        const int ii = i + offsets[s] * di;
        const int jj = j + offsets[s] * dj;
        const int kk = k + offsets[s] * dk;
        if (SF::IBM::isIbmNonFluidCell(field, ii, jj, kk)
            && field.CellFlag(ii, jj, kk) != IBM_GHOST_CELL) {
            return false;
        }
    }
    return true;
}

/// @brief 判断一个面是否直接连接流体和IBM非流体单元。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @return 一侧为FLUID_CELL且另一侧为IBM_GHOST_CELL/SOLID_CELL时返回true。
inline bool faceIsIBMInterface(const Field& field, int i, int j, int k, Math::Dir d) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    const bool leftFluid = SF::IBM::isFluidCell(field, i, j, k);
    const bool rightFluid = SF::IBM::isFluidCell(field, i + di, j + dj, k + dk);
    const bool leftIBM = SF::IBM::isIbmNonFluidCell(field, i, j, k);
    const bool rightIBM = SF::IBM::isIbmNonFluidCell(field, i + di, j + dj, k + dk);
    return (leftFluid && rightIBM) || (leftIBM && rightFluid);
}

/// @brief 判断一个面是否连接当前patch的物理流体点。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @return 任一侧是当前patch物理区域内的FLUID_CELL时返回true。
///
/// 分块边界外侧的通信halo只提供模板数据，不在当前patch推进。若面两侧
/// 都没有本地物理流体点，该面不会进入任何流体单元残差，无需装配WENO。
inline bool faceTouchesPhysicalFluid(const Field& field,
                                     int i, int j, int k,
                                     Math::Dir d) {
    int di, dj, dk;
    Math::dirOffset(d, di, dj, dk);
    return SF::IBM::isFluidCell(field, i, j, k)
        || SF::IBM::isFluidCell(field, i + di, j + dj, k + dk);
}

/// @brief 把未乘面积的面通量写入Field面通量数组。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @param geom 面几何量。
/// @param flux 未乘面积的法向通量。
inline void storePhysicalFlux(FluxField& fluxField, const Field& field,
                              int i, int j, int k, Math::Dir d,
                              const FaceGeometry& geom, const double flux[5]) {
    double areaFlux[5];
    for (int v = 0; v < 5; ++v) areaFlux[v] = flux[v] * geom.area;
    Math::storeFlux(fluxField, field, i, j, k, d, areaFlux);
}

} // namespace Flux
} // namespace SF
