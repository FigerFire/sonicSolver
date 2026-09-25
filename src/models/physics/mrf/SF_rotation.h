/// @file SF_rotation.h
/// @brief MRF 平动/旋转参考系源项模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_sourceAssembly.h"

namespace SF {
namespace Source {
namespace Rotating {

inline Vector3 relativeVelocity(const Field& field, int i, int j, int k,
                                const Vector3& frameVelocity) {
    return MRF::cellVelocity(field, i, j, k) - frameVelocity;
}

inline Vector3 centrifugalAcceleration(const Vector3& radius,
                                       const Vector3& omegaVector) {
    return -1.0 * cross(omegaVector, cross(omegaVector, radius));
}

inline Vector3 coriolisAcceleration(const Vector3& velocityRelative,
                                    const Vector3& omegaVector) {
    return -2.0 * cross(omegaVector, velocityRelative);
}

inline Vector3 rotatingFrameAcceleration(const Field& field, int i, int j, int k,
                                         const Vector3& origin,
                                         const Vector3& omegaVector,
                                         const Vector3& frameVelocity) {
    Vector3 radius = MRF::cellPosition(field, i, j, k) - origin;
    Vector3 velocityRelative = relativeVelocity(field, i, j, k, frameVelocity);
    return coriolisAcceleration(velocityRelative, omegaVector)
         + centrifugalAcceleration(radius, omegaVector);
}

inline void addRotatingFrame(Field& field, Residual& residual,
                             int i, int j, int k,
                             const Vector3& origin,
                             const Vector3& omegaVector,
                             const Vector3& frameVelocity = Vector3()) {
    Vector3 accel = rotatingFrameAcceleration(field, i, j, k,
                                              origin, omegaVector, frameVelocity);
    MRF::addBodyForceSource(field, residual, i, j, k, accel,
                            MRF::cellVelocity(field, i, j, k));
}

inline void addSource(Field& field, Residual& residual, int i, int j, int k,
                      const std::vector<RotatingSetting>& settings) {
    for (const auto& setting : settings) {
        if (MRF::appliesToZone(field, setting.zone, i, j, k)) {
            Vector3 frameVelocity = setting.hasVelocity ? setting.velocity : Vector3();
            addRotatingFrame(field, residual, i, j, k,
                             setting.center,
                             MRF::angularVelocity(setting),
                             frameVelocity);
        }
    }
}

} // namespace Rotating
} // namespace Source
} // namespace SF
