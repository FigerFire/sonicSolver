#pragma once

/// @file SF_termProviderCatalog.h
/// @brief Deterministic numerical-provider selection for compiled equation terms.

#include "SF_runtimeRequirements.h"
#include "core/system/SF_sourceProvider.h"

#include <functional>
#include <optional>

namespace SF::System {

struct FormulaOperatorBinding {
    FormulaOperatorBinding()=default;
    FormulaOperatorBinding(std::string equation,std::string address,std::string implementation)
        :formula(std::move(equation)),occurrence(std::move(address)),provider(std::move(implementation)) {
        if(formula.empty() && !occurrence.empty())operation=Math::resolveOperator(occurrence);
    }
    std::optional<Math::Op> operation;
    /// Empty formula with an operator name selects the global default;
    /// empty occurrence selects a formula-level default. Both set selects
    /// one stable operator occurrence. More specific bindings win.
    std::string formula;
    std::string occurrence;
    std::string provider;
};

/// @brief Resolve occurrence > formula default > global operator default.
const FormulaOperatorBinding* selectFormulaBinding(
    const Equation& formula, const FormulaExpr& expression,
    std::string_view occurrence,
    const std::vector<FormulaOperatorBinding>& bindings);

struct TermProviderDescriptor {
    std::string id;
    std::function<bool(const TermMatchContext&)> match;
    std::optional<FDM::TermRecipe> sourceRecipe;
    std::function<PrimitiveMomentumSource()> compilePrimitiveSource;
    PrimitiveMomentumSpatialTerm primitiveSpatial = nullptr;
    std::string owner = "solver.discretization";
    std::function<ConservativeSourceKernel()> compileConservativeSource;
    std::function<StageConservativeSource(const TermMatchContext&)> compileStageSource;
    /// @brief 数值实现支持的数学身份；matcher继续校验布局/系数/具体结构。
    std::vector<Math::Op> operations;
};

struct ResolvedTermProvider {
    std::string id;
    BindingStatus status = BindingStatus::Unsupported;
    std::string reason;
    std::optional<FDM::TermRecipe> recipe;
    PrimitiveMomentumSource primitiveSource;
    ConservativeSourceKernel conservativeSource;
    StageConservativeSource stageSource;
    PrimitiveMomentumSpatialTerm primitiveSpatial = nullptr;
    std::string owner;
    bool compiledDataAvailable = false;
};

/// @brief Open registrations; zero/one/multiple matches have explicit outcomes.
class TermProviderCatalog {
public:
    void add(TermProviderDescriptor descriptor);
    void addSource(SourceTermProviderDescriptor descriptor);
    ResolvedTermProvider resolve(const TermMatchContext& context,
        const std::vector<FormulaOperatorBinding>& bindings = {}) const;
    static TermProviderCatalog builtIn();
private:
    std::vector<TermProviderDescriptor> descriptors_;
};

} // namespace SF::System
