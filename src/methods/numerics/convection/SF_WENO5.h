/// @file SF_WENO5.h
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <cmath>

namespace SF {
namespace WENO5 {

/// WENO5 重构核心 (Jiang & Shu 1996)
double weno5_core(const double* v);   // v[5]

} // namespace WENO5
} // namespace SF
