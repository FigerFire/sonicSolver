#pragma once

/// @file SF_numericalSystem.h
/// @brief HOW — 数学 term 到 built-in numerical recipe 的冻结绑定。
///        NumericalCompiler 的冻结输出；不持有 Field、MPI 或执行循环。

#include "core/config/SF_configTypes.h"
#include "core/system/SF_expression.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace SF::System {

/// @brief 编译后的时间推进策略：HOW 的时间方向部分。
///
/// runtime 只从这里读取 time recipe，不再回到 raw SolverConfig 重读选择。
struct CompiledTimeIntegration {
    FDM::TimeRecipe recipe;
};

/// @brief 编译后的时间步长限制：HOW 的时间步长部分。
struct TimeStepPolicy {
    double cfl = 0.0;
    double maxDeltaT = std::numeric_limits<double>::max();
};

/// @brief Eulerian phase transport 的冻结数值输入；不是 pressure coupling preset。
struct CompiledPhaseTransport {
    FDM::PhaseConvectionScheme convection = FDM::PhaseConvectionScheme::Upwind;
    double sourceCfl = 0.5;
};

/// @brief 一个 mathematical term 到 built-in numerical recipe 的冻结绑定。
struct BoundTerm {
    std::string equationId;
    Equation::TermKind kind;
    std::string primary;
    std::string secondary;
    std::size_t ordinal = 0;
    FDM::TermRecipe recipe;

    BoundTerm(std::string equation, const Equation::Term& term,
              std::size_t termOrdinal, FDM::TermRecipe boundRecipe)
        : equationId(std::move(equation)), kind(term.kind),
          primary(term.primary.name), secondary(term.secondary.name),
          ordinal(termOrdinal), recipe(std::move(boundRecipe)) {}
};

enum class RecipeConsumerKind { EquationTerm, ExecutableOperation };

/// @brief Selected recipe 的实际编译消费者，供 validation/explain 使用。
struct CompiledRecipeBinding {
    FDM::TermRecipe recipe;
    RecipeConsumerKind kind;
    std::string consumerId;
};

struct CompiledNumericalSystem {
    FDM::NumericalRecipeSet recipes;
    CompiledTimeIntegration time;
    TimeStepPolicy dt;
    CompiledPhaseTransport phaseTransport;
    std::vector<BoundTerm> terms;
    std::vector<CompiledRecipeBinding> recipeBindings;
    int requiredHaloWidth = 0;
    std::vector<std::string> workspaceRequirements;
    std::vector<std::string> providerRequirements;
};

} // namespace SF::System
