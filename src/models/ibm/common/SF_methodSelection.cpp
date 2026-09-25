/// @file SF_methodSelection.cpp
/// @brief IBM 正交方法选择实现。

#include "common/SF_methodSelection.h"

namespace SF::IBM::Common {

FDM::ImmersedMethodSelection selectMethod(
        const IBMRuntimeConfig& config) {
    FDM::ImmersedMethodSelection selection;
    selection.method = config.method;
    if (config.method == FDM::IBMMethod::Ghost) {
        selection.support = FDM::IBMConstraintSupport::Surface;
        selection.representation = FDM::IBMRepresentation::SharpJump;
        selection.enforcement = FDM::IBMEnforcement::GhostCell;
        selection.solid = FDM::IBMSolidModel::Prescribed;
        selection.fluidPorts.push_back(
            {"momentum", "energy", -1, true, true});
        return selection;
    }

    selection.support = config.forcing.constraintSupport;
    selection.representation = config.forcing.representation;
    selection.enforcement = config.forcing.enforcement;
    selection.solid = config.forcing.solidModel;
    if (config.forcing.fluidPorts.empty()) {
        selection.fluidPorts.push_back(
            {"momentum", "energy", -1, true, true});
    } else {
        for (std::size_t phase = 0;
             phase < config.forcing.fluidPorts.size(); ++phase) {
            const std::string& name = config.forcing.fluidPorts[phase];
            selection.fluidPorts.push_back(
                {"momentum." + name, "energy." + name,
                 static_cast<int>(phase), true, true});
        }
    }
    return selection;
}

} // namespace SF::IBM::Common
