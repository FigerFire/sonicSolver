#pragma once

#include "SF_parserCommon.h"
#include "core/config/SF_configTypes.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::FDM {
inline ConvectionScheme parseConvectionScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty()) return ConvectionScheme::WENO5;
    if (t == "weno3") return ConvectionScheme::WENO3;
    if (t == "weno5") return ConvectionScheme::WENO5;
    if (t == "teno5") return ConvectionScheme::TENO5;
    if (t == "weno7") return ConvectionScheme::WENO7;
    throw std::invalid_argument(
        "Unsupported convection scheme '" + value
        + "'. Supported schemes: WENO3, WENO5, TENO5, WENO7. "
          "Recommended pairings: WENO5/TENO5 with StegerWarming or "
          "LaxFriedrichs for smooth shock-capturing cases; WENO3 with "
          "LaxFriedrichs for IBM/ILW debug runs; WENO7 only when ghost depth "
          "and boundary closure are known to be sufficient.");
}
/// @brief Parse the governing-equation formulation.
/// @param value User-facing string such as `"conservativeFluxDifference"`.
/// @return Strongly typed formulation; empty input keeps the conservative default.
inline EquationFormulation parseEquationFormulation(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "conservative" || t == "conservativefd"
        || t == "conservativefinitedifference"
        || t == "conservativefluxdifference"
        || t == "fluxdifference" || t == "fluxdifferencing") {
        return EquationFormulation::ConservativeFluxDifference;
    }
    if (t == "primitive" || t == "primitivedifferential"
        || t == "nonconservative" || t == "differential") {
        return EquationFormulation::PrimitiveDifferential;
    }
    throw std::invalid_argument(
        "Unsupported numerics formulation '" + value
        + "'. Use conservativeFluxDifference for density-based compressible flow.");
}
/// @brief Parse the variable family reconstructed at faces.
/// @param value User-facing string such as `"characteristic"`.
/// @return Strongly typed reconstruction variable; empty input keeps characteristic.
inline ReconstructionVariable parseReconstructionVariable(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "characteristic" || t == "char") {
        return ReconstructionVariable::Characteristic;
    }
    if (t == "conservative" || t == "u") return ReconstructionVariable::Conservative;
    if (t == "primitive" || t == "prim") return ReconstructionVariable::Primitive;
    throw std::invalid_argument(
        "Unsupported reconstruction variable '" + value
        + "'. Current WENO flux assembly supports characteristic reconstruction.");
}
/// @brief Parse the common-face interface flux policy.
/// @param value User-facing string such as `"shared"` or `"sharedInterfaceFlux"`.
/// @return Strongly typed interface policy.
inline InterfaceFluxPolicy parseInterfaceFluxPolicy(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "shared" || t == "sharedflux"
        || t == "sharedinterfaceflux" || t == "commonfaceflux") {
        return InterfaceFluxPolicy::SharedInterfaceFlux;
    }
    throw std::invalid_argument(
        "Unsupported interface flux policy '" + value
        + "'. Current conservative solver requires sharedInterfaceFlux.");
}
/// @brief Parse an inviscid flux method name.
/// @param value User-facing string such as `"StegerWarming"`, `"Rusanov"` or `"Roe"`.
/// @return Strongly typed flux method; empty input keeps `StegerWarming`.
inline FluxSplitter parseFluxSplitter(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty()) return FluxSplitter::StegerWarming;
    if (t == "rusanov" || t == "localaxfriedrichs" || t == "llf") return FluxSplitter::Rusanov;
    if (t == "laxfriedrichs" || t == "lf" || t == "lxf" || t == "laxf") return FluxSplitter::LaxFriedrichs;
    if (t == "roe") return FluxSplitter::Roe;
    if (t == "laxwendroff" || t == "lw") return FluxSplitter::LaxWendroff;
    if (t == "stegerwarming" || t == "sw") return FluxSplitter::StegerWarming;
    throw std::invalid_argument(
        "Unsupported inviscid flux method '" + value
        + "'. Supported methods: StegerWarming, Rusanov, LaxFriedrichs, Roe, "
          "LaxWendroff. Recommended pairings: WENO5/TENO5+StegerWarming for "
          "regular cases, WENO3+LaxFriedrichs for IBM/ILW debug, Roe for "
          "direct face-flux checks.");
}
/// @brief Parse the IBM near-wall WENO closure mode.
/// @param value User-facing string such as `"Downgrade"`, `"LowOrder"` or `"ILW"`.
/// @return Strongly typed IBM boundary scheme; empty input keeps low-order.
inline IBMBoundaryScheme parseIBMBoundaryScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "downgrade" || t == "loworder"
        || t == "firstorder") {
        return IBMBoundaryScheme::LowOrder;
    }
    if (t == "ilw" || t == "highorder" || t == "highorderilw") {
        return IBMBoundaryScheme::ILW;
    }
    throw std::invalid_argument(
        "Unsupported IBM WENO closure '" + value
        + "'. Supported closures: Downgrade/LowOrder or ILW. "
          "Use ILW only when [numerics].ILW is configured with the requested "
          "WENO order; use Downgrade for explicit low-order IBM debugging.");
}
/// @brief 解析 IBM 主方法；未实现的 direct/continuous forcing 不会被替换。
/// @brief Resolve one user token to a complete built-in time recipe.
inline TimeRecipe resolveTimeRecipe(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "classicalrk4")
        return builtInTimeRecipe(TimeRecipeId::ClassicalRK4);
    if (t == "forwardeuler")
        return builtInTimeRecipe(TimeRecipeId::ForwardEuler);
    if (t == "ssprk3") return builtInTimeRecipe(TimeRecipeId::SSPRK3);
    throw std::invalid_argument(
        "Unsupported time recipe '" + value
        + "'. Supported recipes: forwardEuler, SSPRK3, classicalRK4.");
}
/// @brief Parse a viscous central-difference scheme name.
/// @param value User-facing string such as `"CENTRAL2"` or `"CENTRAL4"`.
/// @return Strongly typed viscous scheme; empty input keeps `CENTRAL2`.
inline ViscousScheme parseViscousScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "central2" || t == "cd2") return ViscousScheme::Central2;
    if (t == "central4" || t == "cd4") return ViscousScheme::Central4;
    throw std::invalid_argument(
        "Unsupported viscous scheme '" + value
        + "'. Supported schemes: CENTRAL2, CENTRAL4. Recommended pairing: "
          "CENTRAL2 for current viscous/turbulence validation; CENTRAL4 only "
          "after boundary stencils are verified.");
}
/// @brief Resolve only source-registered built-in convection recipe names.
inline TermRecipe resolveConvectionTermRecipe(const std::string& value) {
    const std::string token = normalizeToken(value);
    for (ConvectionScheme scheme : {ConvectionScheme::WENO3,
                                    ConvectionScheme::WENO5,
                                    ConvectionScheme::TENO5,
                                    ConvectionScheme::WENO7}) {
        for (FluxSplitter flux : {FluxSplitter::StegerWarming,
                                  FluxSplitter::Rusanov,
                                  FluxSplitter::LaxFriedrichs,
                                  FluxSplitter::Roe,
                                  FluxSplitter::LaxWendroff}) {
            const TermRecipe recipe = builtInConvectionRecipe(scheme,flux);
            if (token == normalizeToken(toString(recipe.id()))) return recipe;
        }
    }
    throw std::invalid_argument(
        "Unsupported convection term recipe '"+value
        +"'. Select a registered built-in recipe such as weno7Rusanov or "
         "weno7Steger.");
}
inline TermRecipe resolveDiffusionTermRecipe(const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "central2explicit") {
        return builtInDiffusionRecipe(ViscousScheme::Central2);
    }
    if (token == "central4explicit") {
        return builtInDiffusionRecipe(ViscousScheme::Central4);
    }
    throw std::invalid_argument(
        "Unsupported diffusion term recipe '"+value
        +"'. Supported recipe: central2Explicit or central4Explicit.");
}
// 校验已移入 core/config/SF_numericsPolicy.h（typed 值对象自校验）。
} // namespace SF::FDM
