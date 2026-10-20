/// @file SF_loads.cpp
/// @brief IBM 作用反作用、力矩与离散机械功率。

#include "operations/SF_loads.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Loads {
namespace {

bool finite(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

void requireLoadInput(const Vector3& value, double measure) {
    if (!finite(value) || !std::isfinite(measure) || measure <= 0.0) {
        throw std::runtime_error(
            "IBM load integration requires finite multiplier and positive measure.");
    }
}

} // namespace

Vector3 forceOnBody(const Vector3& fluidMultiplier, double measure) {
    requireLoadInput(fluidMultiplier,measure);
    return fluidMultiplier*(-measure);
}

Vector3 torqueOnBody(
        const Vector3& relativePosition,
        const Vector3& bodyForce) {
    if (!finite(relativePosition) || !finite(bodyForce)) {
        throw std::runtime_error("IBM torque received non-finite input.");
    }
    return cross(relativePosition,bodyForce);
}

double fluidMechanicalPower(
        const Vector3& fluidMultiplier,
        const Vector3& velocity,
        double measure) {
    requireLoadInput(fluidMultiplier,measure);
    if (!finite(velocity)) {
        throw std::runtime_error("IBM power received non-finite velocity.");
    }
    const double result = dot(fluidMultiplier,velocity)*measure;
    if (!std::isfinite(result)) {
        throw std::runtime_error("IBM mechanical power is non-finite.");
    }
    return result;
}

} // namespace SF::IBM::Loads
