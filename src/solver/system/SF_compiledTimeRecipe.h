#pragma once

/// @file SF_compiledTimeRecipe.h
/// @brief Temporal-method output consumed by the existing stage backend.

#include "core/config/types/SF_timeRecipe.h"

#include <array>
#include <cstddef>
#include <stdexcept>

namespace SF::System {

struct TimeStage {
    double abscissa = 0.0;
    double baseWeight = 0.0;
    double incrementWeight = 0.0;
};

/// A frozen backend selected by the TemporalMethod object. The existing
/// numerical stage kernels consume this tag during the migration.
enum class ExplicitStageBackend { ForwardEuler, SSPRK3, ClassicalRK4 };

/// Frozen method mathematics. Construction is owned by ITemporalMethod;
/// execution retains the validated numerical kernels during migration.
class CompiledTimeRecipe {
public:
    static CompiledTimeRecipe fromMethod(
        FDM::TimeRecipe selected, std::array<TimeStage,4> stages,
        int count, ExplicitStageBackend backend,
        std::array<double,4> finalWeights = {},
        double finalDivisor = 1.0) {
        if (count < 1 || count > 4 || count != selected.stageCount()
            || finalDivisor == 0.0)
            throw std::runtime_error("Temporal method stage contract is invalid.");
        CompiledTimeRecipe result;
        result.selected_ = selected;
        result.stages_ = stages;
        result.count_ = count;
        result.backend_ = backend;
        result.finalWeights_ = finalWeights;
        result.finalDivisor_ = finalDivisor;
        return result;
    }

    FDM::TimeRecipeId id() const { return selected_.id(); }
    FDM::TimeRecipe selected() const { return selected_; }
    ExplicitStageBackend backend() const { return backend_; }
    FDM::TimeFamily family() const { return selected_.family(); }
    FDM::TimeTopology topology() const { return selected_.topology(); }
    int order() const { return selected_.order(); }
    int stageCount() const { return count_; }
    void requireProviderStages(int providerStages) const {
        if (providerStages!=stageCount())
            throw std::runtime_error("Explicit time provider stage count disagrees with compiled recipe.");
    }
    const TimeStage& stage(int index) const {
        if (index < 0 || index >= count_)
            throw std::out_of_range("Compiled time stage index is invalid.");
        return stages_[(size_t)index];
    }
    std::array<double,4> finalWeights() const { return finalWeights_; }
    double finalDivisor() const { return finalDivisor_; }

private:
    FDM::TimeRecipe selected_ = FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);
    std::array<TimeStage,4> stages_{};
    int count_ = 1;
    ExplicitStageBackend backend_ = ExplicitStageBackend::ForwardEuler;
    std::array<double,4> finalWeights_{};
    double finalDivisor_ = 1.0;
};

} // namespace SF::System
