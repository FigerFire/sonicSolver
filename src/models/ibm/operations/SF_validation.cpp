/// @file SF_validation.cpp
/// @brief IBM 数值状态校验；无效状态直接拒绝且不做隐式兜底。

#include "operations/SF_validation.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace SF::IBM::Validation {

void requireFinite(double value, std::string_view name) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string(name)+" must be finite.");
    }
}

void requireFinite(const Vector3& value, std::string_view name) {
    if (!std::isfinite(value.x) || !std::isfinite(value.y)
        || !std::isfinite(value.z)) {
        throw std::runtime_error(std::string(name)+" must be finite.");
    }
}

void requirePositive(double value, std::string_view name) {
    requireFinite(value,name);
    if (value <= 0.0) {
        throw std::runtime_error(std::string(name)+" must be positive.");
    }
}

void requireIndex(int index, int size, std::string_view name) {
    if (index < 0 || index >= size) {
        throw std::runtime_error(std::string(name)+" is outside its storage.");
    }
}

void requireMultiplierCount(std::size_t actual, std::size_t expected) {
    if (actual != expected) {
        throw std::runtime_error(
            "IBM KKT multiplier count does not match constraint points.");
    }
}

void requireConstraintResidual(
        double residual,
        double tolerance,
        std::string_view algorithm) {
    requireFinite(residual,"IBM constraint residual");
    requirePositive(tolerance,"IBM constraint tolerance");
    if (residual > tolerance) {
        throw std::runtime_error(
            std::string(algorithm)+" failed its velocity constraint: residual="
            +std::to_string(residual)+", tolerance="
            +std::to_string(tolerance)+".");
    }
}

} // namespace SF::IBM::Validation
