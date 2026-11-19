#pragma once
/// @file SF_timeRecipe.h
/// @brief 不可变的 built-in 时间离散配方 contract。
///
/// 时间方法只以完整 recipe 选择；调用方不单独设置 stage 数或 implicit 标志。

#include <stdexcept>

namespace SF::FDM {
/// @brief 已实现的 built-in 时间离散配方标识；验收范围见各方法报告。
enum class TimeRecipeId { ForwardEuler, SSPRK3, ClassicalRK4, PartitionedFixedEuler };
enum class TimeFamily { RungeKutta, FixedPoint };
enum class TimeTopology { ExplicitStages, FixedInterval };

/// @brief 不可变的时间数值 contract；调用方只能整体选择 built-in 配方。
class TimeRecipe {
public:
    constexpr TimeRecipe() = default;

    static constexpr TimeRecipe builtIn(TimeRecipeId id) {
        switch (id) {
            case TimeRecipeId::PartitionedFixedEuler:
                return {id,TimeFamily::FixedPoint,TimeTopology::FixedInterval,1,1};
            case TimeRecipeId::ForwardEuler:
                return {id,TimeFamily::RungeKutta,
                        TimeTopology::ExplicitStages,1,1};
            case TimeRecipeId::SSPRK3:
                return {id,TimeFamily::RungeKutta,
                        TimeTopology::ExplicitStages,3,3};
            case TimeRecipeId::ClassicalRK4:
                return {id,TimeFamily::RungeKutta,
                        TimeTopology::ExplicitStages,4,4};
        }
        throw std::invalid_argument("Unknown built-in time recipe id.");
    }

    constexpr TimeRecipeId id() const { return id_; }
    constexpr TimeFamily family() const { return family_; }
    constexpr TimeTopology topology() const { return topology_; }
    constexpr int order() const { return order_; }
    constexpr int stageCount() const { return stageCount_; }

private:
    constexpr TimeRecipe(TimeRecipeId id, TimeFamily family,
                         TimeTopology topology, int order, int stageCount)
        : id_(id), family_(family), topology_(topology),
          order_(order), stageCount_(stageCount) {}

    TimeRecipeId id_ = TimeRecipeId::ClassicalRK4;
    TimeFamily family_ = TimeFamily::RungeKutta;
    TimeTopology topology_ = TimeTopology::ExplicitStages;
    int order_ = 4;
    int stageCount_ = 4;
};

constexpr TimeRecipe builtInTimeRecipe(TimeRecipeId id) {
    return TimeRecipe::builtIn(id);
}
/// @brief Convert a time scheme to the canonical config/log token.
/// @param scheme Strongly typed time scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(TimeRecipeId id) {
    switch (id) {
        case TimeRecipeId::PartitionedFixedEuler: return "partitionedFixedEuler";
        case TimeRecipeId::ForwardEuler: return "forwardEuler";
        case TimeRecipeId::SSPRK3: return "SSPRK3";
        case TimeRecipeId::ClassicalRK4: return "classicalRK4";
    }
    return "unknown";
}
inline const char* toString(TimeFamily family) {
    switch (family) {
        case TimeFamily::FixedPoint: return "FixedPoint";
        case TimeFamily::RungeKutta: return "RungeKutta";
    }
    return "unknown";
}
inline const char* toString(TimeTopology topology) {
    switch (topology) {
        case TimeTopology::FixedInterval: return "FixedInterval";
        case TimeTopology::ExplicitStages: return "ExplicitStages";
    }
    return "unknown";
}
} // namespace SF::FDM

