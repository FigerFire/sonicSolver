#pragma once

/// @file SF_methodObjects.h
/// @brief Numerical method contracts; HOW supplies the output, methods supply
///        mathematical realization and transitional backend lowering.

#include "SF_compiledTimeRecipe.h"
#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

class TermProviderCatalog;
struct FormulaOperatorBinding;
struct CompiledNumericalSystem;
/// @brief 已迁移的局部数值方法入口；仅按 WHICH 方法选择，不访问模型或运行资源。
std::optional<CompiledEquationCall> compileDirectEquation(const ExecutableEquationSystem&,
    const EquationCall&,const NumericalBinding&);

/// @brief 已迁移方法的有限过程/布局约束，不重新匹配动量 AST。
std::string validateCompiledConservativeStage(const ExecutableEquationSystem&,const CompiledSolvePlan&);
void validateDirectProgram(const ExecutableEquationSystem&,const CompiledNumericalSystem&,const CompiledSolvePlan&);

class ITemporalMethod {
public:
    virtual ~ITemporalMethod() = default;
    virtual FDM::TimeRecipeId id() const = 0;
    virtual CompiledTimeRecipe compile(FDM::TimeRecipe selected) const = 0;
    virtual SolvePlanNode compileFragment(
        const CompiledTimeRecipe& compiled,
        const CompiledEquationCall& residual,
        const std::vector<CompiledMathRef>& calls,
        const std::vector<TemporalParticipant>& participants,
        const std::vector<SolvePlanNode>& physicalStepPrefix) const = 0;
};

class TemporalMethodRegistry {
public:
    void add(const ITemporalMethod& method);
    const ITemporalMethod& at(FDM::TimeRecipeId id) const;
private:
    std::vector<const ITemporalMethod*> methods_;
};

class IProvider {
public:
    virtual ~IProvider() = default;
    virtual std::string_view id() const = 0;
    /// Implementation ownership is part of the selected numerical contract.
    virtual std::string_view runtimeProvider() const { return {}; }
    virtual std::string_view operationProvider(const OpId&) const {
        return runtimeProvider();
    }
    /// @brief 所选实现的空间编译入口。空值表示此方法不提供空间编译。
    virtual SpatialLowering spatialLowering() const { return nullptr; }
    /// @brief 整体残差/未知量集合由所选实现lower；默认拒绝，不展开为逐项执行。
    virtual CompiledEquationCall compileBlock(const ExecutableEquationSystem&,
        const ExecutionScope&,const NumericalBinding&) const {
        throw std::runtime_error("Unsupported implicit block solve =>: selected "+std::string(id())+" has no steady block capability");
    }
    virtual CompiledEquationCall compile(
        const ExecutableEquationSystem& system,
        const EquationCall& step,
        const NumericalBinding& binding) const = 0;
};

class ProviderRegistry {
public:
    void add(const IProvider& method);
    const IProvider& at(std::string_view id) const;
private:
    std::vector<const IProvider*> methods_;
};

TemporalMethodRegistry builtinTemporalMethods();
ProviderRegistry builtinProviders();
} // namespace SF::System
