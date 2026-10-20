#pragma once

/// @file SF_sourceProvider.h
/// @brief Model-owned, solver-independent source term registration value.

#include "core/interfaces/SF_termKernel.h"
#include "core/system/SF_formula.h"
#include "core/config/SF_configTypes.h"

#include <functional>
#include <optional>
#include <string>

namespace SF::System {

struct TermMatchContext {
    const Equation& formula;
    const FormulaExpr& expression;
    std::string_view occurrence;
    std::string_view output;
    std::string_view equationMethod;
    const FDM::TermRecipe* selectedRecipe = nullptr;
};

/// @brief A model supplies a matcher, recipe, and value-owning kernel factories.
struct SourceTermProviderDescriptor {
    std::string id;
    std::function<bool(const TermMatchContext&)> match;
    std::optional<FDM::TermRecipe> sourceRecipe;
    std::function<PrimitiveMomentumSource()> compilePrimitiveSource;
    std::string owner;
    std::function<ConservativeSourceKernel()> compileConservativeSource;
};

} // namespace SF::System
