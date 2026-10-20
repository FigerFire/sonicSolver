/// @file SF_WENO3.h
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <cmath>

namespace SF {
namespace WENO3 {

/// WENO3 重构核心 (4 点模板, 2 子模板各 3 点)
double weno3_core(const double* v);

} // namespace WENO3
} // namespace SF
