/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_roe.h
/// @brief Roe近似Riemann通量模块。
///
/// 该模块只负责Roe面通量，不负责求解器时间推进、边界条件或配置解析。
/// IBM边界附近的特殊处理由ILWBoundary模块提供；Roe模块可以作为主体格式，
/// 也可以作为WENO模板跨入IBM区域时的局部低阶回退格式。

#include "SF_field.h"
#include "SF_face.h"

namespace SF {
namespace Flux {
namespace RoeFlux {

/// @brief 计算单个面上的Roe数值通量。
/// @param qL 左侧守恒状态 [rho,rho*u,rho*v,rho*w,E]。
/// @param qR 右侧守恒状态 [rho,rho*u,rho*v,rho*w,E]。
/// @param normal 面单位法向，方向从左状态指向右状态。
/// @param flux 输出未乘面积的法向数值通量。
void flux(const double qL[5], const double qR[5],
          const double normal[3], double flux[5], double gamma = 1.4);

/// @brief 计算指定半网格面的Roe通量并写入FluxField。
/// @param field 结构网格场，内部面通量数组会被写入。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k, Math::Dir d,
                   double gamma = 1.4);

/// @brief 计算全场三个方向的Roe通量。
/// @param field 结构网格场，内部面通量数组会被写入并装配到残差缓存。
void computeAllFluxes(Field& field, FluxField& fluxField, Residual& residual,
                      double gamma = 1.4);

} // namespace RoeFlux
} // namespace Flux
} // namespace SF
