/// @file SF_diagnostics.cpp
/// @brief IBM 约束、载荷和固体运动日志格式化。

#include "operations/SF_diagnostics.h"

#include <sstream>

namespace SF::IBM::Diagnostics {

std::string describe(
        std::string_view algorithm,
        std::string_view support,
        const FDM::ImmersedConstraintResult& result,
        const Kinematics::SolidKinematics* solid) {
    std::ostringstream output;
    output << algorithm << " " << support << " dofs="
           << result.constrainedCells
           << ", max|Ju-Us|=" << result.maximumVelocityResidual
           << ", max|dJ/du|=" << result.maximumStationarityResidual
           << ", kineticFunctional="
           << result.kineticIncrementFunctional
           << ", bodyForce=(" << result.forceOnBody.x << ","
           << result.forceOnBody.y << "," << result.forceOnBody.z << ")"
           << ", bodyTorque=(" << result.torqueOnBody.x << ","
           << result.torqueOnBody.y << "," << result.torqueOnBody.z << ")"
           << ", fluidPower=" << result.fluidMechanicalPower;
    if (solid) {
        output << ", solidU=(" << solid->linearVelocity.x << ","
               << solid->linearVelocity.y << ","
               << solid->linearVelocity.z << ")"
               << ", solidOmega=(" << solid->angularVelocity.x << ","
               << solid->angularVelocity.y << ","
               << solid->angularVelocity.z << ")";
    }
    return output.str();
}

} // namespace SF::IBM::Diagnostics
