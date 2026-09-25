/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_band.h
/// @brief Level Set窄带标记接口。

#include "SF_field.h"
#include "SF_state.h"

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief 窄带更新工具。
class NarrowBand {
public:
    /// @brief 根据 |phi|<=width 更新窄带掩码。
    /// @param field 参考网格。
    /// @param levelSet Level Set场。
    /// @param width 窄带半宽, 必须为正。
    static void update(const Field& field,
                       LevelSetField& levelSet,
                       double width);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
