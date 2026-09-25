#pragma once
/// @file SF_configTypes.h
/// @brief 求解器配置的强类型值对象入口。
///
/// 各模块的配置类型已按 ownership 拆分，本文件只保留 Solver 直接消费的
/// 聚合配置（NumericsConfig / BoundaryConfig / InitialConditionConfig /
/// SourceConfig / SolverConfig），并重新导出各模块类型：
///   IBM        -> core/config/types/SF_ibmConfigTypes.h
///   湍流       -> core/config/types/SF_turbulenceConfigTypes.h
///   线性求解   -> core/config/types/SF_linearSolverConfigTypes.h
///   压力基算法 -> core/config/types/SF_pressureConfigTypes.h
///   时间配方   -> core/config/types/SF_timeRecipe.h
///   term 配方  -> core/config/types/SF_termRecipe.h

#include "core/state/SF_valueTypes.h"

#include "core/config/types/SF_ibmConfigTypes.h"
#include "core/config/types/SF_turbulenceConfigTypes.h"
#include "core/config/types/SF_pressureConfigTypes.h"
#include "core/config/types/SF_timeRecipe.h"
#include "core/config/types/SF_termRecipe.h"
#include "core/config/types/SF_linearSolverConfigTypes.h"

#include <limits>
#include <optional>
#include <vector>

namespace SF::FDM {

/// @brief 单个求解器实例的空间、时间和输运离散配置。
struct NumericsConfig {
    double startTime = 0.0;
    double cfl = 1.0;
    double maxDeltaT = std::numeric_limits<double>::max();
    EquationFormulation formulation =
        EquationFormulation::ConservativeFluxDifference;
    ConvectionScheme convection = ConvectionScheme::WENO5;
    ReconstructionVariable reconstruction =
        ReconstructionVariable::Characteristic;
    FluxSplitter flux = FluxSplitter::StegerWarming;
    InterfaceFluxPolicy interfaceFlux =
        InterfaceFluxPolicy::SharedInterfaceFlux;
    IBMBoundaryScheme ibmBoundary = IBMBoundaryScheme::LowOrder;
    int ilwOrder = 0;
    TimeRecipe timeRecipe = builtInTimeRecipe(TimeRecipeId::ClassicalRK4);
    bool timeRecipeDeclared = false;
    NumericalRecipeSet recipes{
        timeRecipe,
        TermRecipe::convectionRecipe(
            TermRecipeId::Weno5Steger,ConvectionScheme::WENO5,
            FluxSplitter::StegerWarming,3),
        std::nullopt};
    bool termRecipesDeclared = false;
    ViscousScheme viscous = ViscousScheme::Central2;
    bool viscousEnabled = false;
    double dynamicViscosity = 0.0;
    double prandtl = 0.72;
    double idealGasGamma = 1.4;
    double idealGasConstant = 287.05;
};

/// @brief 求解器边界值配置。
struct BoundaryConfig {
    bool ilwEnabled = false;
    int ilwOrder = 0;
    std::vector<BCSetting<Vector3>> velocity;
    std::vector<BCSetting<double>> energyFromPressure;
    std::vector<BCSetting<double>> density;
    std::vector<ThermalBCSetting> thermal;
    double thermalDynamicViscosity = 0.0;
    double thermalPrandtl = 0.72;
};

/// @brief 单流体守恒状态初始化输入。
struct InitialConditionConfig {
    std::vector<BCSetting<Vector3>> velocity;
    std::vector<BCSetting<double>> pressure;
    std::vector<BCSetting<double>> density;
    std::vector<BCSetting<double>> energy;
    std::vector<BCSetting<double>> temperature;
};

/// @brief 显式源项配置。
struct SourceConfig {
    std::vector<SourceKind> enabled;
    std::vector<ZoneVectorSetting> gravity;
    std::vector<RotatingSetting> rotating;
    std::vector<WallHeatSetting> wallHeat;
};

/// @brief Solver Algorithm 消费的完整强类型配置。
struct SolverConfig {
    NumericsConfig numerics;
    PressureCorrectionConfig pressure;
    InitialConditionConfig initial;
    BoundaryConfig boundaries;
    IBMConfig ibm;
    SourceConfig sources;
    TurbulenceConfig turbulence;
};

/// @brief 所有配置块都缺失时的 typed 基线。
///
/// native case 解码从这个值开始，再逐个语义块覆盖；它不读取任何 parser
/// 全局变量。这里显式写出的三个字段是历史上来自 legacy 暂存对象初值的
/// 默认量（CFL 暂存初值为 0，湍流 family/model 暂存初值为 DNS），其余字段
/// 与 typed 默认值一致。
/// @return 未配置任何语义块时的求解器配置。
inline SolverConfig defaultSolverConfig() {
    SolverConfig config;
    config.numerics.cfl = 0.0;
    config.turbulence.family = TurbulenceFamily::DNS;
    config.turbulence.model = TurbulenceModelKind::DNS;
    return config;
}

} // namespace SF::FDM
