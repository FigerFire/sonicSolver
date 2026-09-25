/// @file SF_WENO7.h
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <cmath>

namespace SF {
namespace WENO7 {

/// WENO7 重构核心
double weno7_core(const double* v);   // v[7]

} // namespace WENO7
} // namespace SF
