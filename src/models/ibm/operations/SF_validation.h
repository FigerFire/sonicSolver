#pragma once

/// @file SF_validation.h
/// @brief IBM 动态数值状态的 fail-fast 校验操作。

#include "SF_valueTypes.h"

#include <cstddef>
#include <string_view>

namespace SF::IBM::Validation {

void requireFinite(double value, std::string_view name);
void requireFinite(const Vector3& value, std::string_view name);
void requirePositive(double value, std::string_view name);
void requireIndex(int index, int size, std::string_view name);
void requireMultiplierCount(
    std::size_t actual,
    std::size_t expected);
void requireConstraintResidual(
    double residual,
    double tolerance,
    std::string_view algorithm);

} // namespace SF::IBM::Validation
