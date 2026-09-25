/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_mass.h
/// @brief Level Set相质量统计和质量修正接口。

#include "SF_field.h"
#include "SF_state.h"

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief 质量修正接口。
class MassCorrection {
public:
    /// @brief 统计液相Heaviside指示函数之和。
    /// @param field 参考网格。
    /// @param levelSet Level Set场。
    /// @param epsilon Heaviside正则化半厚度；0表示锐界面统计。
    static double liquidIndicatorSum(const Field& field,
                                     const LevelSetField& levelSet,
                                     double epsilon);

    /// @brief 对phi施加全局质量修正。
    ///
    /// 该接口预留给后续Sussman/Fatemi质量约束步骤。当前不实现自动改phi，
    /// 显式调用会fail-fast，避免隐藏地平移界面。
    static void applyGlobalShift(Field& field,
                                 LevelSetField& levelSet,
                                 double targetLiquidIndicatorSum,
                                 double epsilon);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
