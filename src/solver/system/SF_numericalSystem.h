#pragma once

/// @file SF_numericalSystem.h
/// @brief WHICH — 数学 term 到 built-in numerical recipe 的冻结绑定。
///        NumericalCompiler 的冻结输出；不持有 Field、MPI 或执行循环。

#include "core/config/SF_configTypes.h"
#include "core/system/SF_equationIR.h"
#include "core/interfaces/SF_termKernel.h"
#include "SF_compiledTimeRecipe.h"
#include "solver/discretization/pressure/SF_momentumTermKernels.h"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace SF { class Field; }

namespace SF::System {

using PrimitiveMomentumSpatialTerm = Discretization::PressureMomentum::Kernel;

/// @brief 编译后的时间推进策略：WHICH 的 temporal lowering。
///
/// runtime 只从这里读取 time recipe，不再回到 raw SolverConfig 重读选择。
struct CompiledTimeIntegration {
    CompiledTimeRecipe recipe;
};

/// @brief 编译后的时间步长限制：WHICH 的时间步长限制。
struct TimeStepPolicy {
    double cfl = 0.0;
    double maxDeltaT = std::numeric_limits<double>::max();
};

/// @brief Eulerian phase transport 的冻结数值输入；不是 pressure coupling preset。
struct CompiledPhaseTransport {
    FDM::PhaseConvectionScheme convection = FDM::PhaseConvectionScheme::Upwind;
    double sourceCfl = 0.5;
};

enum class PressureFaceCoupling { None, RhieChow };

/// @brief Frozen inputs for constant-density pressure operations.
/// rhoConst remains in CompiledStateRealization; this value contains only numerical inputs.
struct PressureNumericalConfig {
    double dynamicViscosity = 0.0;
    double velocityRelaxation = 1.0;
    double pressureRelaxation = 1.0;
    double relativeTolerance = 1.0e-2;
    double absoluteTolerance = 1.0e-10;
    FDM::PressureReference reference;
    FDM::LinearSolverConfig pressureLinear;
    std::vector<BCSetting<Vector3>> velocityBoundary;
    std::vector<BCSetting<double>> pressureBoundary;
};

/// @brief A frozen Equation AST occurrence bound to its spatial provider.
/// The source path is independent of legacy SF::Equation::Term ordinals/roles.
struct CompiledSpatialBinding {
    std::string formulaId;
    std::string executionEquationId;
    std::string output;
    std::string equationMethod;
    std::string occurrence;
    std::string operatorName;
    std::string primary;
    std::string secondary;
    enum class Side { Left, Right } side = Side::Left;
    std::optional<FDM::TermRecipe> recipe;
    std::string provider;
    std::string providerOwner;
    bool compiledDataAvailable = false;
    PrimitiveMomentumSource primitiveSource;
    ConservativeSourceKernel conservativeSource;
    StageConservativeSource stageSource;
    PrimitiveMomentumSpatialTerm primitiveSpatial = nullptr;

    CompiledSpatialBinding(std::string formula, std::string executionEquation,
              std::string target, std::string method, std::string path,
              std::string mathOperator, std::string symbol,
              std::string coefficient, Side termSide,
              std::optional<FDM::TermRecipe> boundRecipe,
              std::string providerId, std::string owner,
              bool dataAvailable,
              PrimitiveMomentumSource sourceEvaluator = {},
              ConservativeSourceKernel conservativeEvaluator = {},
              PrimitiveMomentumSpatialTerm spatialEvaluator = nullptr)
        : formulaId(std::move(formula)),
          executionEquationId(std::move(executionEquation)),
          output(std::move(target)),equationMethod(std::move(method)),
          occurrence(std::move(path)),operatorName(std::move(mathOperator)),
          primary(std::move(symbol)),secondary(std::move(coefficient)),
          side(termSide),
          recipe(std::move(boundRecipe)),provider(std::move(providerId)),
          providerOwner(std::move(owner)),
          compiledDataAvailable(dataAvailable),
          primitiveSource(std::move(sourceEvaluator)),
          conservativeSource(std::move(conservativeEvaluator)),
          primitiveSpatial(spatialEvaluator) {}
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
    PressureFaceCoupling pressureFaceCoupling = PressureFaceCoupling::None;
    std::optional<PressureNumericalConfig> pressureOperator;
    std::vector<CompiledSpatialBinding> operators;
    std::vector<CompiledRecipeBinding> recipeBindings;
    int requiredHaloWidth = 0;
    std::vector<std::string> workspaceRequirements;
    std::vector<std::string> providerRequirements;
};

} // namespace SF::System
