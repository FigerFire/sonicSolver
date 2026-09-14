/// @file SF_pressureJump.cpp
/// @brief 锐利界面压力跳跃对压力方程的离散修正。

#include "solver/algorithm/pressureBased/SF_pressureJump.h"

#include "SF_thermodynamicClosure.h"

#include <cmath>
#include <stdexcept>

namespace SF::PressureBased {

FaceCorrectionJump pressureCorrectionJump(
        const FDM::IInterfaceJumpCondition* condition,
        const Field& field,
        int leftI, int leftJ, int leftK,
        int axis,
        double relaxation) {
    FaceCorrectionJump result;
    if (!condition) return result;
    if (!std::isfinite(relaxation) || relaxation <= 0.0
        || relaxation > 1.0) {
        throw std::runtime_error(
            "pressureCorrectionJump requires relaxation in (0,1].");
    }

    const FDM::InterfacePressureJump physical =
        condition->pressureJump(field, leftI, leftJ, leftK, axis);
    if (!physical.crossesInterface) return result;
    if (!std::isfinite(physical.fractionFromLeft)
        || physical.fractionFromLeft < 0.0
        || physical.fractionFromLeft > 1.0
        || !std::isfinite(physical.targetRightMinusLeft)) {
        throw std::runtime_error(
            "pressureCorrectionJump received an invalid interface jump.");
    }

    int rightI = leftI, rightJ = leftJ, rightK = leftK;
    if (axis == 0) ++rightI;
    else if (axis == 1) ++rightJ;
    else if (axis == 2) ++rightK;
    else {
        throw std::runtime_error(
            "pressureCorrectionJump axis must be 0, 1, or 2.");
    }
    const double current =
        Boundary::pressureAt(field,rightI,rightJ,rightK)
        - Boundary::pressureAt(field,leftI,leftJ,leftK);
    const double correction =
        (physical.targetRightMinusLeft - current) / relaxation;
    if (!std::isfinite(current) || !std::isfinite(correction)) {
        throw std::runtime_error(
            "pressureCorrectionJump produced a non-finite correction jump.");
    }

    result.active = true;
    result.fractionFromLeft = physical.fractionFromLeft;
    result.targetRightMinusLeft = physical.targetRightMinusLeft;
    result.correctionRightMinusLeft = correction;
    return result;
}

} // namespace SF::PressureBased
