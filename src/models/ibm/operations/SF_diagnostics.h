#pragma once

/// @file SF_diagnostics.h
/// @brief IBM 结果的纯格式化接口。

#include "SF_immersedConstraint.h"
#include "operations/SF_kinematics.h"

#include <string>
#include <string_view>

namespace SF::IBM::Diagnostics {

std::string describe(
    std::string_view algorithm,
    std::string_view support,
    const FDM::ImmersedConstraintResult& result,
    const Kinematics::SolidKinematics* solid = nullptr);

} // namespace SF::IBM::Diagnostics
