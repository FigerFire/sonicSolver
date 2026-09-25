#pragma once
/// @file SF_termRecipe.h
/// @brief 数学 term 的空间离散职责与 built-in recipe contract。

#include "core/config/types/SF_timeRecipe.h"

#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace SF::FDM {
/// @brief 无黏对流通量散度格式。
enum class ConvectionScheme { WENO3, WENO5, TENO5, WENO7 };

/// @brief 控制方程离散形式。
enum class EquationFormulation {
    ConservativeFluxDifference,
    PrimitiveDifferential
};

/// @brief 界面左右状态的重构变量族。
enum class ReconstructionVariable {
    Characteristic,
    Conservative,
    Primitive
};

/// @brief 公共面的数值通量所有权策略。
enum class InterfaceFluxPolicy { SharedInterfaceFlux };

/// @brief 无黏面通量或通量分裂方法。
enum class FluxSplitter {
    StegerWarming,
    Rusanov,
    LaxFriedrichs,
    Roe,
    LaxWendroff
};

/// @brief 黏性通量中心差分阶数。
enum class ViscousScheme { Central2, Central4 };
/// @brief 显式源项族。
enum class SourceKind { Gravity, MRF, WallHeat };

/// @brief 数学 term 的空间离散职责；不描述全局执行顺序。
enum class TermRole { Convection, Diffusion, Source };
enum class TemporalRole { ExplicitResidual };

/// @brief 已注册且具有 production provider 的 built-in term recipe。
enum class TermRecipeId {
    Weno3Steger, Weno3Rusanov, Weno3LaxFriedrichs, Weno3Roe,
    Weno3LaxWendroff,
    Weno5Steger, Weno5Rusanov, Weno5LaxFriedrichs, Weno5Roe,
    Weno5LaxWendroff,
    Teno5Steger, Teno5Rusanov, Teno5LaxFriedrichs, Teno5Roe,
    Teno5LaxWendroff,
    Weno7Steger, Weno7Rusanov, Weno7LaxFriedrichs, Weno7Roe,
    Weno7LaxWendroff,
    Central2Explicit, Central4Explicit,
    GravityExplicit, MRFExplicit, WallHeatExplicit
};

/// @brief 不可逐项修改的空间 term 数值 contract。
class TermRecipe {
public:
    constexpr TermRecipeId id() const { return id_; }
    constexpr TermRole role() const { return role_; }
    constexpr TemporalRole temporalRole() const { return temporalRole_; }
    constexpr ConvectionScheme convection() const { return convection_; }
    constexpr ReconstructionVariable reconstruction() const {
        return reconstruction_;
    }
    constexpr FluxSplitter flux() const { return flux_; }
    constexpr ViscousScheme diffusion() const { return diffusion_; }
    constexpr SourceKind source() const { return source_; }
    constexpr int haloWidth() const { return haloWidth_; }

    static constexpr TermRecipe convectionRecipe(
            TermRecipeId id, ConvectionScheme convection,
            FluxSplitter flux, int haloWidth) {
        return {id,TermRole::Convection,TemporalRole::ExplicitResidual,
                convection,ReconstructionVariable::Characteristic,flux,
                ViscousScheme::Central2,SourceKind::Gravity,haloWidth};
    }
    static constexpr TermRecipe diffusionRecipe(
            TermRecipeId id, ViscousScheme diffusion) {
        return {id,TermRole::Diffusion,TemporalRole::ExplicitResidual,
                ConvectionScheme::WENO5,
                ReconstructionVariable::Characteristic,
                FluxSplitter::StegerWarming,diffusion,
                SourceKind::Gravity,0};
    }
    static constexpr TermRecipe sourceRecipe(
            TermRecipeId id, SourceKind source) {
        return {id,TermRole::Source,TemporalRole::ExplicitResidual,
                ConvectionScheme::WENO5,
                ReconstructionVariable::Characteristic,
                FluxSplitter::StegerWarming,ViscousScheme::Central2,
                source,0};
    }

private:
    constexpr TermRecipe(
            TermRecipeId id, TermRole role, TemporalRole temporalRole,
            ConvectionScheme convection,
            ReconstructionVariable reconstruction, FluxSplitter flux,
            ViscousScheme diffusion, SourceKind source, int haloWidth)
        : id_(id), role_(role), temporalRole_(temporalRole),
          convection_(convection), reconstruction_(reconstruction), flux_(flux),
          diffusion_(diffusion), source_(source), haloWidth_(haloWidth) {}

    TermRecipeId id_;
    TermRole role_;
    TemporalRole temporalRole_;
    ConvectionScheme convection_;
    ReconstructionVariable reconstruction_;
    FluxSplitter flux_;
    ViscousScheme diffusion_;
    SourceKind source_;
    int haloWidth_;
};

/// @brief 用户解析完成后的 numerical recipe selection。
struct NumericalRecipeSet {
    TimeRecipe time;
    std::optional<TermRecipe> convection;
    std::optional<TermRecipe> diffusion;
};
/// @brief Convert a convective scheme to the canonical config/log token.
/// @param scheme Strongly typed scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(ConvectionScheme scheme) {
    switch (scheme) {
        case ConvectionScheme::WENO3: return "WENO3";
        case ConvectionScheme::WENO5: return "WENO5";
        case ConvectionScheme::TENO5: return "TENO5";
        case ConvectionScheme::WENO7: return "WENO7";
    }
    return "WENO5";
}
/// @brief Convert formulation to the canonical config/log token.
inline const char* toString(EquationFormulation formulation) {
    switch (formulation) {
        case EquationFormulation::ConservativeFluxDifference:
            return "conservativeFluxDifference";
        case EquationFormulation::PrimitiveDifferential:
            return "primitiveDifferential";
    }
    return "conservativeFluxDifference";
}
/// @brief Convert reconstruction variable to the canonical config/log token.
inline const char* toString(ReconstructionVariable variable) {
    switch (variable) {
        case ReconstructionVariable::Characteristic: return "characteristic";
        case ReconstructionVariable::Conservative: return "conservative";
        case ReconstructionVariable::Primitive: return "primitive";
    }
    return "characteristic";
}
/// @brief Convert interface flux policy to the canonical config/log token.
inline const char* toString(InterfaceFluxPolicy policy) {
    switch (policy) {
        case InterfaceFluxPolicy::SharedInterfaceFlux:
            return "sharedInterfaceFlux";
    }
    return "sharedInterfaceFlux";
}
/// @brief Convert a flux method to the canonical config/log token.
/// @param splitter Strongly typed flux method.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(FluxSplitter splitter) {
    switch (splitter) {
        case FluxSplitter::StegerWarming: return "StegerWarming";
        case FluxSplitter::Rusanov: return "Rusanov";
        case FluxSplitter::LaxFriedrichs: return "LaxFriedrichs";
        case FluxSplitter::Roe: return "Roe";
        case FluxSplitter::LaxWendroff: return "LaxWendroff";
    }
    return "unknown";
}
inline TermRecipeId convectionRecipeId(
        ConvectionScheme scheme, FluxSplitter flux) {
    const int schemeOffset = scheme == ConvectionScheme::WENO3 ? 0
        : scheme == ConvectionScheme::WENO5 ? 5
        : scheme == ConvectionScheme::TENO5 ? 10 : 15;
    const int fluxOffset = flux == FluxSplitter::StegerWarming ? 0
        : flux == FluxSplitter::Rusanov ? 1
        : flux == FluxSplitter::LaxFriedrichs ? 2
        : flux == FluxSplitter::Roe ? 3 : 4;
    return static_cast<TermRecipeId>(schemeOffset+fluxOffset);
}
/// @brief Resolve an already parsed legacy convection pair to one built-in recipe.
inline TermRecipe builtInConvectionRecipe(
        ConvectionScheme scheme, FluxSplitter flux) {
    const int halo = scheme == ConvectionScheme::WENO3 ? 2
        : (scheme == ConvectionScheme::WENO7 ? 4 : 3);
    return TermRecipe::convectionRecipe(
        convectionRecipeId(scheme,flux),scheme,flux,halo);
}
inline TermRecipe builtInDiffusionRecipe(ViscousScheme scheme) {
    return TermRecipe::diffusionRecipe(
        scheme == ViscousScheme::Central2
            ? TermRecipeId::Central2Explicit
            : TermRecipeId::Central4Explicit,
        scheme);
}
inline TermRecipe builtInSourceRecipe(SourceKind source) {
    switch (source) {
        case SourceKind::Gravity:
            return TermRecipe::sourceRecipe(
                TermRecipeId::GravityExplicit,source);
        case SourceKind::MRF:
            return TermRecipe::sourceRecipe(TermRecipeId::MRFExplicit,source);
        case SourceKind::WallHeat:
            return TermRecipe::sourceRecipe(
                TermRecipeId::WallHeatExplicit,source);
    }
    throw std::invalid_argument("Unknown built-in source term recipe.");
}
inline const char* toString(TermRecipeId id) {
    static constexpr const char* names[] = {
        "weno3Steger", "weno3Rusanov", "weno3LaxFriedrichs", "weno3Roe",
        "weno3LaxWendroff",
        "weno5Steger", "weno5Rusanov", "weno5LaxFriedrichs", "weno5Roe",
        "weno5LaxWendroff",
        "teno5Steger", "teno5Rusanov", "teno5LaxFriedrichs", "teno5Roe",
        "teno5LaxWendroff",
        "weno7Steger", "weno7Rusanov", "weno7LaxFriedrichs", "weno7Roe",
        "weno7LaxWendroff",
        "central2Explicit", "central4Explicit",
        "gravityExplicit", "mrfExplicit", "wallHeatExplicit"
    };
    const auto index = static_cast<std::size_t>(id);
    if (index >= sizeof(names)/sizeof(names[0])) return "unknown";
    return names[index];
}
inline const char* toString(TermRole role) {
    switch (role) {
        case TermRole::Convection: return "convection";
        case TermRole::Diffusion: return "diffusion";
        case TermRole::Source: return "source";
    }
    return "unknown";
}
inline const char* toString(TemporalRole role) {
    if (role == TemporalRole::ExplicitResidual) return "ExplicitResidual";
    return "unknown";
}
/// @brief Convert a viscous scheme to the canonical config/log token.
/// @param scheme Strongly typed viscous scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(ViscousScheme scheme) {
    switch (scheme) {
        case ViscousScheme::Central2: return "CENTRAL2";
        case ViscousScheme::Central4: return "CENTRAL4";
    }
    return "CENTRAL2";
}
/// @brief Convert source kinds to a compact source-list token.
/// @param kinds Ordered source kinds.
/// @return `"None"` when empty, otherwise names joined with `+`.
inline std::string toString(const std::vector<SourceKind>& kinds) {
    if (kinds.empty()) return "None";

    std::ostringstream out;
    for (size_t i = 0; i < kinds.size(); ++i) {
        if (i > 0) out << "+";
        switch (kinds[i]) {
            case SourceKind::Gravity: out << "Gravity"; break;
            case SourceKind::MRF: out << "MRF"; break;
            case SourceKind::WallHeat: out << "WallHeat"; break;
        }
    }
    return out.str();
}
} // namespace SF::FDM

