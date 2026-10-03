#pragma once

/// @file SF_methodObjects.h
/// @brief Numerical method contracts; HOW supplies the output, methods supply
///        mathematical realization and transitional backend lowering.

#include "SF_compiledTimeRecipe.h"
#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"

#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

class ITemporalMethod {
public:
    virtual ~ITemporalMethod() = default;
    virtual FDM::TimeRecipeId id() const = 0;
    virtual CompiledTimeRecipe compile(FDM::TimeRecipe selected) const = 0;
    virtual SolvePlanNode compileFragment(
        const CompiledTimeRecipe& compiled,
        const CompiledEquationCall& residual,
        const std::vector<SolvePlanNode>& physicalStepPrefix = {}) const = 0;
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
    /// Global term recipes are consumed only by providers exposing spatial assembly.
    virtual bool usesSpatialRecipes() const { return false; }
    virtual std::string_view operationProvider(const OpId&) const {
        return runtimeProvider();
    }
    /// Provider-local lifecycle only decorates compiled scopes; source stays immutable.
    virtual void lowerLifecycle(const ExecutionScope&, bool,
        const std::vector<CompiledEquationCall>&, SolvePlanNode&) const {}
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
/// @brief Explicit WHAT/STATE/HOW/WHICH inputs; STATE never infers call order.
CompiledExecutionProgram compileExecutionProgram(
    const ExecutableEquationSystem& system,const StateRegistry& state,
    const ExecutionProgram& program,
    const std::vector<NumericalBinding>& bindings,
    const ProviderRegistry& registry,
    const CompiledTimeRecipe* time = nullptr,
    const ITemporalMethod* temporalMethod = nullptr);

/// @brief Composition-snapshot convenience; delegates to the independent channels.
inline CompiledExecutionProgram compileExecutionProgram(
        const ExecutableEquationSystem& system,const ExecutionProgram& program,
        const std::vector<NumericalBinding>& bindings,const ProviderRegistry& registry,
        const CompiledTimeRecipe* time=nullptr,const ITemporalMethod* temporalMethod=nullptr) {
    return compileExecutionProgram(system,system.state,program,bindings,registry,time,temporalMethod);
}

/// Lower equation-level HOW to executable control flow for method providers
/// that have already supplied their backend OpIds. Temporal fragments are
/// compiled by ITemporalMethod and take precedence for transient bodies.
SolvePlanNode compileMethodProgram(const CompiledExecutionProgram& program, const ProviderRegistry* providers = nullptr);

} // namespace SF::System
