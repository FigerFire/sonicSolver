#pragma once
/// @file SF_turbulenceConfigTypes.h
/// @brief 湍流模型选择与常数配置值对象。

#include "core/state/SF_valueTypes.h"

#include <string>
#include <vector>

namespace SF::FDM {
/// @brief 湍流模型大类。
enum class TurbulenceFamily { None, RAS, LES, DNS };

/// @brief 具体湍流模型。
enum class TurbulenceModelKind {
    None,
    kEpsilon,
    kOmegaSST,
    Smagorinsky,
    DNS
};

/// @brief 湍流输运标量的初值和边界输入。
struct TurbulenceScalarConfig {
    std::vector<BCSetting<double>> kInitial;
    std::vector<BCSetting<double>> epsilonInitial;
    std::vector<BCSetting<double>> omegaInitial;
    std::vector<BCSetting<double>> kBoundary;
    std::vector<BCSetting<double>> epsilonBoundary;
    std::vector<BCSetting<double>> omegaBoundary;
};

/// @brief 湍流模型常数。
struct TurbulenceCoefficients {
    double cMu = 0.09;
    double c1 = 1.44;
    double c2 = 1.92;
    double sigmaK = 1.0;
    double sigmaEpsilon = 1.3;
    double betaStar = 0.09;
    double beta1 = 0.075;
    double gamma1 = 5.0 / 9.0;
    double a1 = 0.31;
    double sigmaOmega1 = 0.5;
    double cSmagorinsky = 0.17;
    double filterScale = 1.0;
    double kFloor = 1.0e-10;
    double epsilonFloor = 1.0e-10;
    double omegaFloor = 1.0e-10;
};

/// @brief 完整湍流配置值对象。
struct TurbulenceConfig {
    bool enabled = false;
    double laminarDynamicViscosity = 0.0;
    TurbulenceFamily family = TurbulenceFamily::None;
    TurbulenceModelKind model = TurbulenceModelKind::None;
    std::vector<std::string> phaseNames;
    std::vector<std::string> wallPatches;
    double turbulentPrandtl = 0.9;
    bool coupleMomentum = true;
    bool coupleEnergy = true;
    TurbulenceScalarConfig scalars;
    TurbulenceCoefficients coefficients;
};
/// @brief 转换湍流大类到稳定日志字符串。
/// @param family 强类型湍流大类。
/// @return 稳定字符串。
inline const char* toString(TurbulenceFamily family) {
    switch (family) {
        case TurbulenceFamily::None: return "None";
        case TurbulenceFamily::RAS: return "RAS";
        case TurbulenceFamily::LES: return "LES";
        case TurbulenceFamily::DNS: return "DNS";
    }
    return "None";
}
/// @brief 转换具体湍流模型到稳定日志字符串。
/// @param model 强类型模型枚举。
/// @return 稳定字符串。
inline const char* toString(TurbulenceModelKind model) {
    switch (model) {
        case TurbulenceModelKind::None: return "None";
        case TurbulenceModelKind::kEpsilon: return "kEpsilon";
        case TurbulenceModelKind::kOmegaSST: return "kOmegaSST";
        case TurbulenceModelKind::Smagorinsky: return "Smagorinsky";
        case TurbulenceModelKind::DNS: return "DNS";
    }
    return "None";
}
} // namespace SF::FDM

