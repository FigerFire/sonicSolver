/// @file SF_mrf.h
/// @brief MRF 平动/旋转参考系源项模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_sourceAssembly.h"
#include "SF_rotation.h"
#include "SF_translation.h"

namespace SF {
namespace Source {
namespace MRF {

inline void addSource(Field& field, Residual& residual, int i, int j, int k,
                      const std::vector<RotatingSetting>& settings) {
    Rotating::addSource(field, residual, i, j, k, settings);
    Translation::addSource(field, residual, i, j, k);
}

} // namespace MRF
} // namespace Source
} // namespace SF
