/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_lxF.h
/// @brief Lax-Friedrichs通量分裂与一阶面通量接口。
///
/// 该模块是用户配置`flux = "LaxFriedrichs"`对应的Flux层实现。
/// WENO层可调用`split()`完成特征WENO后的正/负通量分裂；近壁低阶回退
/// 或测试代码可调用`flux()`/`storeFaceFlux()`完成一阶左右常值面通量。

#include "SF_field.h"
#include "SF_face.h"

namespace SF {
namespace Flux {
namespace LaxFriedrichsFlux {

/// @brief Lax-Friedrichs正/负通量分裂。
/// @param q 守恒状态 [rho,rho*u,rho*v,rho*w,E]。
/// @param normal 面单位法向。
/// @param fPos 输出正通量。
/// @param fNeg 输出负通量。
void split(const double q[5], const double normal[3],
           double fPos[5], double fNeg[5], double gamma = 1.4);

/// @brief 计算一阶Lax-Friedrichs面通量。
/// @param qL 左侧守恒状态。
/// @param qR 右侧守恒状态。
/// @param normal 面单位法向。
/// @param flux 输出未乘面积的法向数值通量。
void flux(const double qL[5], const double qR[5],
          const double normal[3], double flux[5], double gamma = 1.4);

/// @brief 计算指定面的Lax-Friedrichs通量并写入FluxField。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k, Math::Dir d,
                   double gamma = 1.4);

/// @brief 计算全场三个方向的一阶Lax-Friedrichs通量。
/// @param field 结构网格场，内部面通量数组会被写入并装配到残差缓存。
void computeAllFluxes(Field& field, FluxField& fluxField, Residual& residual,
                      double gamma = 1.4);

} // namespace LaxFriedrichsFlux
} // namespace Flux
} // namespace SF
