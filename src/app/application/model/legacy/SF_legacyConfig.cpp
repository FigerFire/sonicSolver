/// @file SF_legacyConfig.cpp
/// @brief 初值设置与 legacy 配置迁移实现。
///
/// STATUS: Legacy —— 本文件不属于任何 CMake target，也不被任何生产路径
/// include。它引用的 `SolverAlgorithm` / parser 全局暂存值已从 runtime
/// authority 中删除，仅作为历史映射留在磁盘上供人工比对。

#include "SF_legacyConfig.h"

#include "SF_config.h"
#include "SF_parameter.h"
#include "core/state/SF_idealGasDefaults.h"

namespace SF::Legacy {

FDM::SolverConfig makeSolverConfig() {
    FDM::SolverConfig config;
    config.numerics.solver = solverApplication == "pressureBase"
        ? FDM::SolverAlgorithm::PressureBased
        : FDM::SolverAlgorithm::DensityBased;
    config.numerics.startTime = startTime;
    config.numerics.cfl = CFL;
    config.numerics.maxDeltaT = maxDeltaT;
    config.numerics.formulation =
        FDM::parseEquationFormulation(solverFormulation);
    config.numerics.convection =
        FDM::parseConvectionScheme(convectionScheme);
    config.numerics.reconstruction =
        FDM::parseReconstructionVariable(reconstructionVariable);
    config.numerics.flux = FDM::parseFluxSplitter(fluxSplitter);
    config.numerics.interfaceFlux =
        FDM::parseInterfaceFluxPolicy(interfaceFluxPolicy);
    config.numerics.ibmBoundary = enableILW
        ? FDM::IBMBoundaryScheme::ILW
        : FDM::parseIBMBoundaryScheme(ibmBoundaryScheme);
    config.numerics.ilwOrder = ilwOrder;
    config.numerics.timeRecipe = FDM::resolveTimeRecipe(timeRecipeName);
    config.numerics.timeRecipeDeclared =
        !FDM::normalizeToken(timeRecipeName).empty();
    config.numerics.recipes.time = config.numerics.timeRecipe;
    config.numerics.viscous = FDM::parseViscousScheme(viscousScheme);
    config.numerics.termRecipesDeclared = termRecipesDeclared;
    if (termRecipesDeclared) {
        // Native term input is a complete selection boundary.  Missing roles
        // must remain absent so NumericalCompiler can reject a mathematical
        // term without a recipe instead of inheriting a hidden default.
        config.numerics.recipes.convection.reset();
        config.numerics.recipes.diffusion.reset();
        if (!FDM::normalizeToken(convectionTermRecipe).empty()) {
            config.numerics.recipes.convection =
                FDM::resolveConvectionTermRecipe(convectionTermRecipe);
            config.numerics.convection =
                config.numerics.recipes.convection->convection();
            config.numerics.reconstruction =
                config.numerics.recipes.convection->reconstruction();
            config.numerics.flux = config.numerics.recipes.convection->flux();
        }
        if (!FDM::normalizeToken(diffusionTermRecipe).empty()) {
            config.numerics.recipes.diffusion =
                FDM::resolveDiffusionTermRecipe(diffusionTermRecipe);
            config.numerics.viscous =
                config.numerics.recipes.diffusion->diffusion();
        }
    } else {
        config.numerics.recipes.convection = FDM::builtInConvectionRecipe(
            config.numerics.convection,config.numerics.flux);
        if (enableViscous) {
            config.numerics.recipes.diffusion =
                FDM::builtInDiffusionRecipe(config.numerics.viscous);
        }
    }
    config.numerics.viscousEnabled = enableViscous;
    config.numerics.dynamicViscosity = mu;
    config.numerics.prandtl = Prandtl;
    config.numerics.idealGasGamma = DefaultIdealGasGamma;
    config.numerics.idealGasConstant = DefaultIdealGasConstant;

    config.pressure.maxIterations = pressureMaxIterations;
    config.pressure.relativeTolerance = pressureRelativeTolerance;
    config.pressure.absoluteTolerance = pressureAbsoluteTolerance;
    config.pressure.relaxation = pressureRelaxation;
    config.pressure.velocityRelaxation = pressureVelocityRelaxation;

    config.initial.velocity = U_IC;
    config.initial.pressure = p_IC;
    config.initial.density = rho_IC;
    config.initial.energy = rhoE_IC;
    config.initial.temperature = T_IC;

    config.boundaries.ilwEnabled = enableILW;
    config.boundaries.ilwOrder = ilwOrder;
    config.boundaries.velocity = U_BC;
    config.boundaries.energyFromPressure = p_BC;
    config.boundaries.density = rho_BC;
    config.boundaries.thermal = T_BC;
    config.boundaries.thermalDynamicViscosity = mu;
    config.boundaries.thermalPrandtl = Prandtl;

    config.ibm.enabled = enableIBM;
    config.ibm.normalAngleDegrees = ilwNormalAngleDegrees;
    config.ibm.normalSearchMinLayers = ilwNormalSearchMinLayers;
    config.ibm.normalSearchMaxLayers = ilwNormalSearchMaxLayers;
    config.ibm.normalSearchTargetCandidates =
        ilwNormalSearchTargetCandidates;
    config.ibm.normalSearchKeepSamples = ilwNormalSearchKeepSamples;

    config.sources.enabled = FDM::parseSourceKinds(sourceScheme);
    config.sources.gravity = gravitySettings;
    config.sources.rotating = rotatingSettings;
    config.sources.wallHeat = wallHeatSettings;

    config.turbulence.enabled = enableTurbulence;
    config.turbulence.laminarDynamicViscosity = mu;
    config.turbulence.family =
        FDM::parseTurbulenceFamily(turbulenceFamily);
    config.turbulence.model =
        FDM::parseTurbulenceModelKind(turbulenceModel);
    config.turbulence.phaseNames = turbulencePhases;
    config.turbulence.wallPatches = turbulenceWallPatches;
    config.turbulence.turbulentPrandtl = turbulentPrandtl;
    config.turbulence.coupleMomentum = turbulenceCoupleMomentum;
    config.turbulence.coupleEnergy = turbulenceCoupleEnergy;
    config.turbulence.scalars.kInitial = k_IC;
    config.turbulence.scalars.epsilonInitial = epsilon_IC;
    config.turbulence.scalars.omegaInitial = omega_IC;
    config.turbulence.scalars.kBoundary = k_BC;
    config.turbulence.scalars.epsilonBoundary = epsilon_BC;
    config.turbulence.scalars.omegaBoundary = omega_BC;
    config.turbulence.coefficients.cMu = turbulenceCMu;
    config.turbulence.coefficients.c1 = turbulenceC1;
    config.turbulence.coefficients.c2 = turbulenceC2;
    config.turbulence.coefficients.sigmaK = turbulenceSigmaK;
    config.turbulence.coefficients.sigmaEpsilon = turbulenceSigmaEpsilon;
    config.turbulence.coefficients.betaStar = turbulenceBetaStar;
    config.turbulence.coefficients.beta1 = turbulenceBeta1;
    config.turbulence.coefficients.gamma1 = turbulenceGamma1;
    config.turbulence.coefficients.a1 = turbulenceA1;
    config.turbulence.coefficients.sigmaOmega1 = turbulenceSigmaOmega1;
    config.turbulence.coefficients.cSmagorinsky = turbulenceCSmagorinsky;
    config.turbulence.coefficients.filterScale = turbulenceFilterScale;
    config.turbulence.coefficients.kFloor = turbulenceKFloor;
    config.turbulence.coefficients.epsilonFloor = turbulenceEpsilonFloor;
    config.turbulence.coefficients.omegaFloor = turbulenceOmegaFloor;

    FDM::validateNumericsConfig(config.numerics);
    return config;
}

} // namespace SF::Legacy
