/// @file SF_compressible.cpp
/// @brief 单流体可压缩守恒方程项及源项表达式的组装。

#include "solver/equation/compressible/SF_compressible.h"

#include "SF_config.h"
#include "solver/discretization/SF_discretization.h"
#include "SF_field.h"

#include <algorithm>

namespace SF::Equation::Compressible {

void System::begin(FluxField& fluxField, Residual& residual) const {
    fluxField.clear();
    residual.clear();
}

void System::convection(
        Field& field, FluxField& fluxField, Residual& residual,
        const AssemblyContext& context) const {
    if (!context.convection) return;
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
    if (context.diffusion) {
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
    if (context.sources && !context.sources->empty()) {
        SourceTerm::Sp(field,residual,*context.sources);
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
