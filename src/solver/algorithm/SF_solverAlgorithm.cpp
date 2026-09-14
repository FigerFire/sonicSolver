/// @file SF_solverAlgorithm.cpp
/// @brief 根据 SolverConfig 创建密度基或压力基流动算法。

#include "solver/algorithm/SF_solverAlgorithm.h"

#include "solver/algorithm/densityBased/SF_correction.h"

#include <memory>
#include <stdexcept>

#include "solver/algorithm/pressureBased/SF_adapter.h"

namespace SF::FDM {

std::unique_ptr<IFlowAlgorithm> makeFlowAlgorithm(
        const SolverConfig& config) {
    if (config.numerics.solver == SolverAlgorithm::DensityBased) {
        return std::make_unique<DensityBased::Algorithm>();
    }
    if (config.numerics.solver == SolverAlgorithm::PressureBased) {
        return std::make_unique<PressureBased::Algorithm>(config);
    }
    throw std::runtime_error("Unknown flow algorithm family.");
}

} // namespace SF::FDM
