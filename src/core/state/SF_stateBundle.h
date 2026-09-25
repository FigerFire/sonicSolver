#pragma once

/// @file SF_stateBundle.h
/// @brief workflow 可提交的非拥有物理状态包。

#include "SF_variableRegistry.h"
#include "SF_distributedField.h"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SF::Physics::FluidStateModel { class Model; }

namespace SF::State {

/// @brief workflow 可提交的状态包；不接管 Field 所有权。
///
/// `patches` 是参与求解的唯一 patch 枚举。这里的 state 不是 MPI canonical
/// entity 的同义词；唯一 owner/COPY/SUM 语义仍由 distributed runtime 管理。
struct StateBundle {
    std::vector<Field*> patches;
    std::shared_ptr<const Physics::FluidStateModel::Model> stateModel;
    VariableRegistry transported;
    /// @brief 主守恒量、模型辅助量和算法临时量的统一分布式注册表。
    DistributedFieldRegistry distributed;
    double time = 0.0;
    double dt = 0.0;
    int step = 0;

    void validatePatches() const {
        if (patches.empty()) {
            throw std::runtime_error(
                "StateBundle requires at least one participating patch.");
        }
        for (const Field* patch : patches) {
            if (!patch) {
                throw std::runtime_error(
                    "StateBundle contains a null participating patch.");
            }
            if (stateModel && patch->stateModel()
                && patch->stateModel() != stateModel) {
                throw std::runtime_error(
                    "StateBundle equation binding differs from a participating Field.");
            }
        }
    }

    /// @brief 验证 density-based timestep 的唯一热力学 binding。
    void validateDensityEquationBinding() const {
        validatePatches();
        if (!stateModel) {
            throw std::runtime_error(
                "Density-based StateBundle requires a FluidStateModel binding.");
        }
        for (size_t index = 0; index < patches.size(); ++index) {
            const auto& actual = patches[index]->stateModel();
            if (actual == stateModel) continue;
            std::ostringstream message;
            message << "Density-based StateBundle FluidStateModel mismatch at patch "
                    << index << "; expected=" << stateModel.get()
                    << ", actual=" << actual.get() << ".";
            throw std::runtime_error(message.str());
        }
    }

    Field& singlePatch() {
        validatePatches();
        if (patches.size() != 1) {
            throw std::runtime_error(
                "StateBundle::singlePatch requires exactly one patch.");
        }
        return *patches.front();
    }

    const Field& singlePatch() const {
        return const_cast<StateBundle*>(this)->singlePatch();
    }

    /// @brief 注册全部 patch 的主守恒状态；状态包不接管 Field 所有权。
    void registerConservativeState() {
        validatePatches();
        for (size_t patch = 0; patch < patches.size(); ++patch) {
            distributed.add(conservativeView(
                "conservative", (int)patch, *patches[patch]));
        }
    }

    /// @brief 注册一个模型辅助标量。
    void registerAuxiliary(
            const VariableDescriptor& descriptor,
            int blockId, Field& geometry, ScalarField& value,
            ExchangeKind exchange = ExchangeKind::State) {
        distributed.add(scalarView(
            descriptor.name, blockId, geometry, value,
            descriptor.halo.depth, descriptor.halo.stages, exchange));
    }
};

} // namespace SF::State
