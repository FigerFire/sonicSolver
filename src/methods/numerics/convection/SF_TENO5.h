/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_TENO5.h
/// @brief 五阶 TENO 标量界面重构核。

namespace SF {
namespace TENO5 {

/// @brief TENO5 重构核心。
/// @param v 五点上风模板值。
/// @return 面左侧重构值。
double teno5_core(const double* v);

} // namespace TENO5
} // namespace SF
