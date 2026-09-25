/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_properties.h
/// @brief 由Level Set更新两相混合物性缓存。

#include "SF_field.h"
#include "SF_phaseProperties.h"
#include "SF_state.h"

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief Level Set物性更新工具。
class Properties {
public:
    /// @brief 根据phi更新LevelSetField内的密度和动力黏度缓存。
    /// @param field 参考网格。
    /// @param levelSet Level Set场。
    /// @param config 多相配置。
    static void update(const Field& field,
                       LevelSetField& levelSet,
                       const MultiPhaseConfig& config);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
