/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_lowOrder.h
/// @brief IBM近壁WENO模板不可用时的一阶低阶面通量。
///
/// 低阶通量不是用户层的主体格式，而是一个显式命名的局部回退模块。
/// 当WENO模板碰到IBM_GHOST_CELL/SOLID_CELL时，求解器不能继续跨固体区域
/// 做高阶重构，此时使用本模块的一阶左右状态面通量闭合该面。

#include "SF_field.h"
#include "SF_face.h"
#include "SF_config.h"

namespace SF {
namespace Flux {
namespace LowOrderFlux {

/// @brief 近IBM回退通量类型。
enum class Method {
    /// @brief 一阶左右常值 + Roe Riemann通量。
    FirstOrderRoe,
    /// @brief 一阶左右常值 + Lax-Friedrichs通量分裂。
    FirstOrderLaxFriedrichs,
    /// @brief 一阶左右常值 + Steger-Warming通量分裂。
    FirstOrderStegerWarming,
    /// @brief 一阶左右常值 + FluidStateModel-aware Rusanov 通量。
    FirstOrderRusanov
};

/// @brief 从用户通量方法映射到近IBM低阶回退方法。
/// @param method 用户配置的主体flux方法。
/// @return 对应的低阶回退方法；Lax-Wendroff近壁默认回退到Roe。
Method fallbackFor(FDM::FluxSplitter method);

/// @brief 计算一个一阶低阶面通量。
/// @param qL 左侧守恒状态。
/// @param qR 右侧守恒状态。
/// @param normal 面单位法向，方向从左侧指向右侧。
/// @param method 低阶回退通量类型。
/// @param flux 输出未乘面积的法向数值通量。
void flux(const double qL[5], const double qR[5],
          const double normal[3],
          Method method,
          double flux[5],
          double gamma = 1.4);

/// @brief 计算指定面的一阶低阶通量并写入Field残差数组。
/// @param field 结构网格场。
/// @param i 面左侧单元的i索引。
/// @param j 面左侧单元的j索引。
/// @param k 面左侧单元的k索引。
/// @param d 面所属的计算坐标方向。
/// @param method 低阶回退通量类型。
void storeFaceFlux(FluxField& fluxField, Field& field, int i, int j, int k,
                   Math::Dir d, Method method,
                   double gamma = 1.4);

} // namespace LowOrderFlux
} // namespace Flux
} // namespace SF
