#pragma once

/// @file SF_parameter.h
/// @brief legacy case parser 的暂存状态；新求解代码应消费强类型配置值对象。

#include "core/state/SF_valueTypes.h"

#include <array>
#include <string>
#include <vector>

namespace SF {

// controlDict legacy 暂存值。
extern double CFL;
extern double maxDeltaT;
extern double startTime;
extern double endTime;
extern bool output;
extern int saveByStep;
extern double saveByTime;
extern bool writeInitial;
extern int endStep;

// fvSchemes 和 solver 配置 legacy 暂存值。
extern std::string convectionScheme;
extern std::string timeRecipeName;
extern std::string fluxSplitter;
extern std::string solverFormulation;
extern std::string reconstructionVariable;
extern std::string interfaceFluxPolicy;
extern std::string viscousScheme;
extern std::string convectionTermRecipe;
extern std::string diffusionTermRecipe;
extern bool termRecipesDeclared;
extern std::string ibmBoundaryScheme;
extern int ilwOrder;
extern double ilwNormalAngleDegrees;
extern int ilwNormalSearchMinLayers;
extern int ilwNormalSearchMaxLayers;
extern int ilwNormalSearchTargetCandidates;
extern int ilwNormalSearchKeepSamples;
extern std::string sourceScheme;
extern bool enableViscous;
extern double mu;
extern double dynamicViscosity;
extern double Prandtl;
extern bool createMesh;
extern std::string solverApplication;
extern int pressureMaxIterations;
extern double pressureRelativeTolerance;
extern double pressureAbsoluteTolerance;
extern double pressureRelaxation;
extern double pressureVelocityRelaxation;

// turbulence legacy 暂存值。
extern bool enableTurbulence;
extern std::string turbulenceFamily;
extern std::string turbulenceModel;
extern double turbulenceCMu;
extern double turbulenceC1;
extern double turbulenceC2;
extern double turbulenceSigmaK;
extern double turbulenceSigmaEpsilon;
extern double turbulenceBetaStar;
extern double turbulenceBeta1;
extern double turbulenceGamma1;
extern double turbulenceA1;
extern double turbulenceSigmaOmega1;
extern double turbulenceCSmagorinsky;
extern double turbulenceFilterScale;
extern double turbulenceKFloor;
extern double turbulenceEpsilonFloor;
extern double turbulenceOmegaFloor;
extern std::vector<std::string> turbulencePhases;
extern std::vector<std::string> turbulenceWallPatches;
extern double turbulentPrandtl;
extern bool turbulenceCoupleMomentum;
extern bool turbulenceCoupleEnergy;

// MPI、IBM 和物理模块 legacy 暂存值。
extern bool parallelEnabled;
extern int parallelNProcs;
extern bool parallelPartitionsAuto;
extern std::array<int, 3> parallelPartitions;
extern double parallelHaloTolerance;
extern bool enableIBM;
extern bool enableILW;
extern bool enableGravity;
extern bool enableMRF;
extern bool enableWallHeatSource;
extern std::vector<ZoneVectorSetting> gravitySettings;
extern std::vector<RotatingSetting> rotatingSettings;
extern std::vector<WallHeatSetting> wallHeatSettings;

// 边界条件 legacy 暂存值。
extern std::vector<BCSetting<Vector3>> U_BC;
extern std::vector<BCSetting<double>> p_BC;
extern std::vector<BCSetting<double>> rho_BC;
extern std::vector<BCSetting<double>> k_BC;
extern std::vector<BCSetting<double>> epsilon_BC;
extern std::vector<BCSetting<double>> omega_BC;
extern std::vector<ThermalBCSetting> T_BC;

// 初始条件 legacy 暂存值。
extern std::vector<BCSetting<Vector3>> U_IC;
extern std::vector<BCSetting<double>> p_IC;
extern std::vector<BCSetting<double>> rho_IC;
extern std::vector<BCSetting<double>> rhoE_IC;
extern std::vector<BCSetting<double>> k_IC;
extern std::vector<BCSetting<double>> epsilon_IC;
extern std::vector<BCSetting<double>> omega_IC;
extern std::vector<BCSetting<double>> T_IC;
extern Vector3 U_0;
extern double p_0;
extern double rho_0;
extern double T_0;

} // namespace SF
