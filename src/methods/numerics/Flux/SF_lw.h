/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_lw.h
/// @brief 二阶Lax-Wendroff/Richtmyer预测-校正对流通量。
///
/// LW通量需要当前时间步长，因此它通过divDispatch的dt参数进入，
/// 不通过全局变量读取时间步。IBM边界面仍由ILWBoundary默认接管。

#include "SF_field.h"
#include "SF_face.h"

namespace SF {
namespace Flux {
namespace LaxWendroffFlux {

/// @brief 计算一个面的Lax-Wendroff通量。
/// @param qL 左侧守恒状态。
/// @param qR 右侧守恒状态。
/// @param normal 面单位法向。
/// @param dt 当前显式时间步长。
/// @param spacing 左右相邻点沿法向的距离。
/// @param flux 输出未乘面积的法向通量。
void flux(const double qL[5], const double qR[5],
          const double normal[3],
          double dt, double spacing,
          double flux[5], double gamma = 1.4);

/// @brief 计算指定半网格面的LW通量并写入FluxField。
/// @param field 结构网格场，内部面通量数组会被写入。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @param dt 当前显式时间步长。
void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k, Math::Dir d,
                   double dt, double gamma = 1.4);

/// @brief 计算全场三个方向的LW通量。
/// @param field 结构网格场，内部面通量数组会被写入并装配到残差缓存。
/// @param dt 当前显式时间步长。
void computeAllFluxes(Field& field, FluxField& fluxField, Residual& residual,
                      double dt, double gamma = 1.4);

} // namespace LaxWendroffFlux
} // namespace Flux
} // namespace SF
