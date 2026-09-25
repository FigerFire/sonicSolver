#pragma once

/// @file SF_buildRequest.h
/// @brief Composition-only input for constructing a resolved system.

#include "SF_physicsTemplate.h"
#include "SF_couplingStatus.h"
#include "core/system/SF_equationIR.h"
#include "core/config/SF_equationComposition.h"
#include "core/interfaces/SF_buildCapabilities.h"
#include "models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"

#include <optional>
#include <string>
#include <vector>

namespace SF::FDM { struct ImmersedAlgorithmDescriptor; }

namespace SF::System {

/// Resolved physical content for the built-in density single-fluid preset.
struct SingleFluidPresetSpec {
    bool diffusion = false;
};

/// @brief 单流体压力约束方程族请求。
///
/// 这是 equation-source 请求，不是 solver family：它声明 WHAT
/// （质量/动量/能量 + 一个 div 类约束，压力作为约束乘子），
/// 不选择 time recipe、coupling preset 或 backend。
struct PressureConstraintSpec {
    bool diffusion = false;
};

/// Startup composition request. Runtime/execution headers do not include this
/// type, which keeps model contribution details out of the resolved contract.
struct BuildRequest {
    PhysicsTemplateKind templateOrigin = PhysicsTemplateKind::SingleFluid;
    std::optional<SingleFluidPresetSpec> singleFluidPreset;
    std::optional<PressureConstraintSpec> pressureConstraint;
    std::vector<std::string> phaseNames;
    std::optional<Physics::InterfaceModels::LevelSetContribution::Spec> levelSet;
    bool homogeneousThermodynamics = false;
    bool legacyMixture = false;
    bool phaseChange = false;
    bool transportedLegacyAlpha = false;
    std::optional<Turbulence::SystemContributionSpec> turbulence;
    /// @brief 显式注册的压力耦合 preset；生效与否由 resolved equation system
    ///        决定（active/inactive/invalid/unsupported），不做静默忽略。
    std::optional<CouplingPresetRequest> coupling;
    /// @brief 用户对已有 equation/未知量的显式修改请求。
    ///
    /// Precedence（§15，typed system layer）：
    ///   builtin defaults → model contributions → user add/extend
    ///   → user replace/disable → validation
    /// 当前只有 add/extend 有 typed lowering；replace/disable 会显式 fail-fast，
    /// 不允许同名后写覆盖前写，也不允许静默忽略。
    std::vector<SystemModification> userModifications;
    bool parallel = false;
    /// @brief 本 binary 真实提供的 backend 能力；由 composition root 填入。
    ///
    /// 默认值是“没有任何可选 backend”，因此单元测试或任何忘记填写的调用点
    /// 在需要 MPI/HYPRE/KKT 的 case 上会得到明确的 unavailable capability，
    /// 而不是静默假设 backend 存在。
    BuildCapabilities capabilities;
    EquationCompositionConfig composition;
    const FDM::ImmersedAlgorithmDescriptor* immersed = nullptr;
};

} // namespace SF::System
