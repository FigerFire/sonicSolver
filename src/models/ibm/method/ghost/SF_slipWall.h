/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#pragma once

/// @file SF_slipWall.h
/// @brief IBM Euler滑移壁面ghost状态构造。
///
/// 本文件只保存滑移壁面的局部状态构造，不读取Field、不做WENO权重、
/// 不负责IBM几何搜索。

namespace SF {
namespace IBM {
namespace GhostILW {

/// @brief 构造一阶Euler滑移壁面ghost状态。
/// @param qFluid 流体侧守恒状态 `[rho,rho*u,rho*v,rho*w,E]`。
/// @param wallNormal 从固体指向流体的壁面单位法向。
/// @param wallVelocity 全局壁面速度。
/// @param gamma 比热比。
/// @param qGhost 输出ghost守恒状态。
void firstOrderEulerWallGhostStateWithVelocity(const double qFluid[5],
                                               const double wallNormal[3],
                                               const double wallVelocity[3],
                                               double gamma,
                                               double qGhost[5]);

/// @brief 构造静止Euler滑移壁面ghost状态。
/// @param qFluid 流体侧守恒状态 `[rho,rho*u,rho*v,rho*w,E]`。
/// @param wallNormal 从固体指向流体的壁面单位法向。
/// @param gamma 比热比。
/// @param qGhost 输出ghost守恒状态。
void firstOrderEulerWallGhostState(const double qFluid[5],
                                   const double wallNormal[3],
                                   double gamma,
                                   double qGhost[5]);

} // namespace GhostILW
} // namespace IBM
} // namespace SF
