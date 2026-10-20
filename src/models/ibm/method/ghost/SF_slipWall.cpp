/// @file SF_slipWall.cpp
/// @brief ghost IBM 的 Euler 滑移壁面反射与守恒状态重建。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "method/ghost/SF_slipWall.h"

#include "SF_eulerState.h"
#include "SF_localFrame.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace IBM {
namespace GhostILW {

void firstOrderEulerWallGhostStateWithVelocity(const double qFluid[5],
                                               const double wallNormal[3],
                                               const double wallVelocity[3],
                                               double gamma,
                                               double qGhost[5]) {
    const Math::LocalFrame frame = Math::makeLocalFrame(wallNormal);
    const Math::Euler::Primitive fluidLocal =
        Math::Euler::localPrimitiveFromConservative(qFluid, frame, gamma);
    const Math::Euler::Primitive wallVelLocal =
        Math::Euler::localWallVelocity(wallVelocity, frame);

    Math::Euler::Primitive ghostLocal = fluidLocal;
    ghostLocal[Math::Euler::LUN] =
        2.0 * wallVelLocal[Math::Euler::LUN] - fluidLocal[Math::Euler::LUN];

    if (!std::isfinite(ghostLocal[Math::Euler::LRHO]) ||
        ghostLocal[Math::Euler::LRHO] <= 0.0 ||
        !std::isfinite(ghostLocal[Math::Euler::LP]) ||
        ghostLocal[Math::Euler::LP] <= 0.0) {
        std::cerr << "[SF FATAL] IBM ILW slip-wall primitive is invalid"
                  << std::endl;
        std::exit(1);
    }

    const Math::Euler::Primitive ghostGlobal =
        Math::Euler::toGlobalPrimitive(ghostLocal, frame);
    Math::Euler::conservativeFromPrimitive(ghostGlobal, gamma, qGhost);

    if (!std::isfinite(qGhost[0]) || !std::isfinite(qGhost[1]) ||
        !std::isfinite(qGhost[2]) || !std::isfinite(qGhost[3]) ||
        !std::isfinite(qGhost[4]) || qGhost[0] <= 0.0) {
        std::cerr << "[SF FATAL] IBM ILW slip-wall ghost state is invalid"
                  << std::endl;
        std::exit(1);
    }
}

void firstOrderEulerWallGhostState(const double qFluid[5],
                                   const double wallNormal[3],
                                   double gamma,
                                   double qGhost[5]) {
    const double wallVelocity[3] = {0.0, 0.0, 0.0};
    firstOrderEulerWallGhostStateWithVelocity(
        qFluid, wallNormal, wallVelocity, gamma, qGhost);
}

} // namespace GhostILW
} // namespace IBM
} // namespace SF
