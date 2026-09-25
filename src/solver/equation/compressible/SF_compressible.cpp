/// @file SF_compressible.cpp
/// @brief 单流体可压缩守恒方程项及源项表达式的组装。

#include "solver/equation/compressible/SF_compressible.h"

#include "SF_config.h"
#include "solver/discretization/SF_discretization.h"
#include "SF_field.h"

#include <algorithm>

namespace SF::Equation::Compressible {

System::System(const Equation::System& definition)
    : definition_(&definition),
      assemblyPlan_(Equation::makeAssemblyPlan(
          definition,{"E_MASS","E_MOMENTUM","E_ENERGY"})) {
}

bool System::contains(TermKind kind) const {
    return std::any_of(
        assemblyPlan_.begin(),assemblyPlan_.end(),
        [kind](const Equation::AssemblyPlan& plan) {
            return plan.contains(kind);
        });
}

void System::begin(FluxField& fluxField, Residual& residual) const {
    fluxField.clear();
    residual.clear();
}

void System::convection(
        Field& field, FluxField& fluxField, Residual& residual,
        const AssemblyContext& context) const {
    if (!contains(TermKind::Divergence)) return;
    if (!context.convection
        || context.convection->role() != FDM::TermRole::Convection) {
        throw std::runtime_error(
            "Compressible convection has no compiled convection TermRecipe.");
    }
    if (!context.thermodynamics) {
        throw std::runtime_error(
            "Reconstructed convection requires an authoritative FluidStateModel binding.");
    }
    const auto equationGamma = context.thermodynamics->perfectGasGamma();
    if (!equationGamma) {
        throw std::runtime_error(
            "Configured reconstructed convection requires a PerfectGas FluidStateModel; "
            "generic FluidStateModel thermodynamics are not supported by the "
            "current characteristic reconstruction implementation.");
    }
    divDispatch(field, fluxField, residual,
                context.convection->convection(),
                context.convection->flux(),
                context.timeStep,
                context.ibmBoundary,
                context.ilwOrder,
                *equationGamma,
                *context.thermodynamics);
}

void System::diffusionAndSources(
        Field& field, Residual& residual, const AssemblyContext& context) const {
    if (contains(TermKind::Diffusion)) {
        if (!context.diffusion
            || context.diffusion->role() != FDM::TermRole::Diffusion) {
            throw std::runtime_error(
                "Compressible diffusion has no compiled diffusion TermRecipe.");
        }
        laplacianDispatch(
            field, residual,
            context.diffusion->diffusion(),
            context.dynamicViscosity,
            context.prandtl,
            context.idealGasGamma,
            context.idealGasConstant,
            context.transport);
    }
    if (contains(TermKind::Source)) {
        if (!context.sources || !context.sourceParameters) {
            throw std::runtime_error(
                "Compressible sources have no compiled source bindings.");
        }
        SourceTerm::Sp(
            field,residual,*context.sourceParameters,*context.sources);
    }
}

void System::assemble(
        Field& field, FluxField& fluxField, Residual& residual,
        const AssemblyContext& context) const {
    begin(fluxField, residual);
    convection(field, fluxField, residual, context);
    diffusionAndSources(field, residual, context);
}

} // namespace SF::Equation::Compressible
