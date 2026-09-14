/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_discretization.h
/// @brief FDM 离散算子总头文件。
///
/// 使用方案: 每个数值格式是一个独立命名空间，通过统一的算子名调用。
///
/// 算子速查:
/// ┌──────────────┬────────────────────────────────────────────┐
/// │ 算子         │ 可用命名空间                               │
/// ├──────────────┼────────────────────────────────────────────┤
/// │ ddt          │ Euler, SSPRK3, RK4                         │
/// │ div          │ WENO3, WENO5, TENO5, WENO7 + flux method   │
/// │ grad/curl    │ CENTRAL2                                   │
/// │ laplacian    │ CENTRAL2, CENTRAL4                         │
/// │ Sp           │ SourceTerm                                 │
/// │ deltaT       │ SF::deltaT (自由函数)                      │
/// └──────────────┴────────────────────────────────────────────┘
///
/// 典型求解器写法:
/// @code
///   // ∂Q/∂t = -∇·F_c + ∇·F_v + S
///   RK4::ddt(field, dt, [&](Field& f) {
///       f.clearResidual();
///       WENO5::div(f, FDM::FluxSplitter::StegerWarming);
///       CENTRAL2::laplacian(f, mu, Pr);
///       SourceTerm::Sp(f, sourceConfig);
///   });
/// @endcode
///
/// NS 方程符号约定 (可压缩守恒形式):
///   ∂Q/∂t + ∇·F_c = ∇·F_v + S
///   → ddt = -div + laplacian + Sp

#include "SF_derivative.h"
#include "solver/discretization/time/SF_time.h"
#include "solver/discretization/convection/SF_convection.h"
#include "solver/discretization/diffusion/SF_diffusion.h"
#include "solver/discretization/source/SF_source.h"
#include "solver/discretization/structured/SF_structured.h"
#include "SF_deltaT.h"
#include "SF_taylor.h"
