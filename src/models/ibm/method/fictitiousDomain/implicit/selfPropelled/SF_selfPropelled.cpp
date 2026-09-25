/// @file SF_selfPropelled.cpp
/// @brief Bhalla Algorithm 7：全隐式自推进 DFM 的固体未知量选择与契约校验。

#include "method/SF_method.h"

#include <stdexcept>

namespace SF::IBM::Forcing {

void ImmersedForcingSystem::validateImplicitAlgorithmContract() const {
    if (config_.forcing.algorithm
        == FDM::IBMForcingAlgorithm::DFMImplicitSelfPropelled) {
        if (config_.forcing.solidModel
                != FDM::IBMSolidModel::SelfPropelledRigid
            || config_.forcing.motion
                != FDM::IBMSolidMotion::CoupledRigid) {
            throw std::runtime_error(
                "dfmImplicitSelfPropelled requires a coupled "
                "selfPropelledRigid generalized velocity in the KKT system.");
        }
        return;
    }
    if (config_.forcing.algorithm
            != FDM::IBMForcingAlgorithm::DFMImplicitPrescribed
        && config_.forcing.algorithm
            != FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian) {
        throw std::runtime_error(
            "Monolithic IBM stage received a non-implicit DFM algorithm.");
    }
}

} // namespace SF::IBM::Forcing
