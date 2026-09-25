/// @file SF_translation.h
/// @brief MRF 平动/旋转参考系源项模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_sourceAssembly.h"

namespace SF {
namespace Source {
namespace Translation {

inline Vector3 velocity{0.0, 0.0, 0.0};
inline Vector3 acceleration{0.0, 0.0, 0.0};

inline void addTranslatingFrame(Field& field, Residual& residual,
                                int i, int j, int k,
                                const Vector3& frameAcceleration) {
    MRF::addBodyForceSource(field, residual, i, j, k,
                            -1.0 * frameAcceleration,
                            MRF::cellVelocity(field, i, j, k));
}

inline void addSource(Field& field, Residual& residual, int i, int j, int k) {
    addTranslatingFrame(field, residual, i, j, k, acceleration);
}

} // namespace Translation
} // namespace Source
} // namespace SF
