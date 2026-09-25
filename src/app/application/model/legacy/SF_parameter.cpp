/// @file SF_parameter.cpp
/// @brief 初值设置与 legacy 配置迁移实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_parameter.h"
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace SF {

    // ── 求解器控制 (全由 system/controlDict 覆盖) ──
    double CFL          = 0.0;
    double maxDeltaT    = std::numeric_limits<double>::max();
    double startTime    = 0.0;
    double endTime      = 0.0;
    bool   output       = false;
    int    saveByStep   = 0;
    double saveByTime   = 0.0;
    bool   writeInitial = true;
    int    endStep      = 0;

    // ── 数值方法 (全由 system/fvSchemes 覆盖) ──
    std::string convectionScheme;
    std::string timeRecipeName;
    std::string fluxSplitter;
    std::string solverFormulation = "conservativeFluxDifference";
    std::string reconstructionVariable = "characteristic";
    std::string interfaceFluxPolicy = "shared";
    std::string viscousScheme = "CENTRAL2";
    std::string convectionTermRecipe;
    std::string diffusionTermRecipe;
    bool        termRecipesDeclared = false;
    std::string ibmBoundaryScheme = "Downgrade";
    int         ilwOrder = 0;
    double      ilwNormalAngleDegrees = 69.51268488527785;
    int         ilwNormalSearchMinLayers = 1;
    int         ilwNormalSearchMaxLayers = 0;
    int         ilwNormalSearchTargetCandidates = 512;
    int         ilwNormalSearchKeepSamples = 256;
    std::string sourceScheme = "None";
    bool        enableViscous = false;
    double      mu = 0.0;
    double      dynamicViscosity = 0.0;
    double      Prandtl = 0.72;
    bool        createMesh = false;
    std::string solverApplication = "densityBase";
    int         pressureMaxIterations = 3000;
    double      pressureRelativeTolerance = 1.0e-2;
    double      pressureAbsoluteTolerance = 1.0e-10;
    double      pressureRelaxation = 1.0e-3;
    double      pressureVelocityRelaxation = 0.1;
    bool        enableTurbulence = false;
    std::string turbulenceFamily = "DNS";
    std::string turbulenceModel = "DNS";
    double      turbulenceCMu = 0.09;
    double      turbulenceC1 = 1.44;
    double      turbulenceC2 = 1.92;
    double      turbulenceSigmaK = 1.0;
    double      turbulenceSigmaEpsilon = 1.3;
    double      turbulenceBetaStar = 0.09;
    double      turbulenceBeta1 = 0.075;
    double      turbulenceGamma1 = 5.0 / 9.0;
    double      turbulenceA1 = 0.31;
    double      turbulenceSigmaOmega1 = 0.5;
    double      turbulenceCSmagorinsky = 0.17;
    double      turbulenceFilterScale = 1.0;
    double      turbulenceKFloor = 1.0e-10;
    double      turbulenceEpsilonFloor = 1.0e-10;
    double      turbulenceOmegaFloor = 1.0e-10;
    std::vector<std::string> turbulencePhases;
    std::vector<std::string> turbulenceWallPatches;
    double      turbulentPrandtl = 0.9;
    bool        turbulenceCoupleMomentum = true;
    bool        turbulenceCoupleEnergy = true;
    bool        parallelEnabled = false;
    int         parallelNProcs = 1;
    bool        parallelPartitionsAuto = true;
    std::array<int, 3> parallelPartitions = {1, 1, 1};
    double      parallelHaloTolerance = 1.0e-8;
    bool        enableIBM = false;
    bool        enableILW = false;
    bool        enableGravity = false;
    bool        enableMRF = false;
    bool        enableWallHeatSource = false;
    std::vector<ZoneVectorSetting> gravitySettings;
    std::vector<RotatingSetting> rotatingSettings;
    std::vector<WallHeatSetting> wallHeatSettings;

    // ── 边界条件 (全由 0/ 场文件 boundaryField 覆盖) ──
    std::vector<BCSetting<Vector3>> U_BC;
    std::vector<BCSetting<double>>  p_BC;
    std::vector<BCSetting<double>>  rho_BC;
    std::vector<BCSetting<double>>  k_BC;
    std::vector<BCSetting<double>>  epsilon_BC;
    std::vector<BCSetting<double>>  omega_BC;
    std::vector<ThermalBCSetting>   T_BC;

    // ── 初始条件 (全由 0/ 场文件 internalField/internalSets 覆盖) ──
    std::vector<BCSetting<Vector3>> U_IC;
    std::vector<BCSetting<double>>  p_IC;
    std::vector<BCSetting<double>>  rho_IC;
    std::vector<BCSetting<double>>  rhoE_IC;
    std::vector<BCSetting<double>>  k_IC;
    std::vector<BCSetting<double>>  epsilon_IC;
    std::vector<BCSetting<double>>  omega_IC;
    std::vector<BCSetting<double>>  T_IC;

    // ── 内部场 ──
    Vector3 U_0{0, 0, 0};
    double  p_0   = 0.0;
    double  rho_0 = 0.0;
    double  T_0   = 0.0;

}
