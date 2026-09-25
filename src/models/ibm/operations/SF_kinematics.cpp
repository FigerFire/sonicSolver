/// @file SF_kinematics.cpp
/// @brief IBM 刚体和变形表面的目标速度恢复。

#include "operations/SF_kinematics.h"

#include "solid/SF_bodyModel.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Kinematics {
namespace {

bool finite(const Vector3& value) {
    return std::isfinite(value.x)
        && std::isfinite(value.y)
        && std::isfinite(value.z);
}

} // namespace

SolidKinematics resolve(
        const Forcing::IBodyModel& body,
        const FDM::ImmersedKKTState& state) {
    SolidKinematics result{body.linearVelocity(),body.angularVelocity()};
    if (state.hasGeneralizedVelocity) {
        for (double value : state.generalizedVelocity) {
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "IBM KKT returned non-finite generalized velocity.");
            }
        }
        result.linearVelocity = {state.generalizedVelocity[0],
                                 state.generalizedVelocity[1],
                                 state.generalizedVelocity[2]};
        result.angularVelocity = {state.generalizedVelocity[3],
                                  state.generalizedVelocity[4],
                                  state.generalizedVelocity[5]};
    }
    if (!finite(result.linearVelocity) || !finite(result.angularVelocity)) {
        throw std::runtime_error(
            "IBM solid kinematics contains a non-finite velocity.");
    }
    return result;
}

Vector3 targetVelocity(
        const SolidKinematics& solid,
        const FDM::ImmersedSurfacePoint& point,
        const Forcing::IBodyModel& body,
        double time) {
    if (!std::isfinite(time)) {
        throw std::runtime_error("IBM target velocity requires finite time.");
    }
    const Vector3 rigid = solid.linearVelocity
        +cross(solid.angularVelocity,point.relativePosition);
    const Vector3 result = rigid
        +body.deformationVelocityAt(point.position,time);
    if (!finite(result)) {
        throw std::runtime_error("IBM target velocity is non-finite.");
    }
    return result;
}

} // namespace SF::IBM::Kinematics
