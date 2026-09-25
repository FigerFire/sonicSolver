/// @file SF_interfaceModel.cpp
/// @brief 界面表示模型的创建与统一接口实现。

#include "SF_interfaceModel.h"

#include "levelSet/SF_levelSet.h"

#include <stdexcept>

namespace SF::Physics::InterfaceModels {

std::unique_ptr<Model> makeInterfaceModel(
        const Multiphase::MultiPhaseConfig& config) {
    if (Multiphase::isLevelSetType(config.type)) {
        return std::make_unique<LevelSetModel>(config);
    }
    throw std::runtime_error(
        "InterfaceModel factory: representation '" + config.type
        + "' is not implemented. Available representation: levelSet.");
}

} // namespace SF::Physics::InterfaceModels
