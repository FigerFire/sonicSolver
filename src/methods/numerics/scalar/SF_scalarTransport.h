/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.02-----------*/

#pragma once

/// @file SF_scalarTransport.h
/// @brief 通用标量输运、边界和有界性诊断工具。

#include "SF_boundaryGeometry.h"
#include "reconstruction/ILW/SF_boundaryClosure.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_scalarField.h"
#include "SF_utility.h"
#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace SF {
namespace FDM {
namespace Scalar {

/// @brief 标量越界处理模式。
enum class BoundMode {
    None,      ///< 不检查也不修改。
    Diagnose,  ///< 只输出越界诊断，不修改数值。
    Project    ///< 显式投影到上下界。
};

/// @brief 标量上下界设置。
struct Bounds {
    bool enabled = false;              ///< true 时启用上下界逻辑。
    double lower = 0.0;                ///< 下界。
    double upper = 1.0;                ///< 上界。
    BoundMode mode = BoundMode::Diagnose; ///< 越界处理模式。
};

/// @brief 一个显式标量输运方程的基础设置。
struct TransportConfig {
    std::string name = "scalar";       ///< 标量名，用于诊断输出。
    bool convectionEnabled = true;      ///< 是否组装一阶迎风守恒对流。
    bool diffusionEnabled = false;      ///< 是否组装中心扩散。
    double diffusivity = 0.0;           ///< 常标量扩散系数。
    bool sourceEnabled = false;         ///< 是否启用常源项。
    double source = 0.0;                ///< 常源项。
    bool ilwBoundaryEnabled = false;     ///< 辅助标量是否使用 ILW 边界闭合。
    int ilwAccuracyOrder = 0;            ///< ILW 精度阶；0 表示关闭。
    Bounds bounds;                      ///< 可选上下界诊断/投影。
};

/// @brief 标量越界诊断结果。
struct BoundReport {
    bool sampled = false;
    int below = 0;
    int above = 0;
    double minValue = std::numeric_limits<double>::max();
    double maxValue = -std::numeric_limits<double>::max();

    /// @brief 是否发现越界。
    bool violated() const { return below > 0 || above > 0; }
};

/// @brief 解析有界处理模式。
/// @param value 用户配置文本。
/// @return 对应模式；空字符串返回 Diagnose。
inline BoundMode parseBoundMode(std::string value) {
    value.erase(std::remove_if(value.begin(), value.end(),
                               [](unsigned char c) {
                                   return c == '_' || c == '-' || std::isspace(c);
                               }),
                value.end());
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    if (value.empty() || value == "diagnose" || value == "warn") {
        return BoundMode::Diagnose;
    }
    if (value == "none" || value == "off") return BoundMode::None;
    if (value == "project" || value == "clip" || value == "bounded") {
        return BoundMode::Project;
    }
    std::cerr << "[SF FATAL] unknown scalar boundMode '" << value
              << "'. Use none, diagnose, or project/clip." << std::endl;
    std::exit(1);
}

/// @brief 有界处理模式名称。
inline const char* toString(BoundMode mode) {
    switch (mode) {
    case BoundMode::None: return "none";
    case BoundMode::Diagnose: return "diagnose";
    case BoundMode::Project: return "project";
    }
    return "unknown";
}

inline void requireCompatible(const Field& flow,
                              const ScalarField& scalar,
                              const char* context) {
    if (!scalar.isCompatibleWith(flow)) {
        std::cerr << "[SF FATAL] " << context
                  << ": scalar field '" << scalar.name()
                  << "' dimensions do not match Field." << std::endl;
        std::exit(1);
    }
}

inline void requireFiniteScalar(const std::string& context,
                                double value,
                                int i,
                                int j,
                                int k) {
    if (!std::isfinite(value)) {
        std::cerr << "[SF FATAL] " << context
                  << ": non-finite scalar at (" << i << "," << j << ","
                  << k << ") value=" << value << std::endl;
        std::exit(1);
    }
}

inline bool isGhostCell(const Field& flow, int i, int j, int k) {
    const int ng = flow.NG();
    return i < ng || i >= flow.NX() + ng
        || j < ng || j >= flow.NY() + ng
        || k < ng || k >= flow.NZ() + ng;
}

inline void nearestInteriorCell(const Field& flow,
                                int i,
                                int j,
                                int k,
                                int& ri,
                                int& rj,
                                int& rk) {
    const int ng = flow.NG();
    ri = std::max(ng, std::min(i, flow.NX() + ng - 1));
    rj = std::max(ng, std::min(j, flow.NY() + ng - 1));
    rk = std::max(ng, std::min(k, flow.NZ() + ng - 1));
}

inline void nearestInsideSample(const Field& flow,
                                int i,
                                int j,
                                int k,
                                int axis,
                                int& si,
                                int& sj,
                                int& sk) {
    const int ng = flow.NG();
    si = i;
    sj = j;
    sk = k;
    if (axis == 0) {
        if (i == ng && flow.NX() > 1) si = i + 1;
        else if (i == flow.NX() + ng - 1 && flow.NX() > 1) si = i - 1;
    } else if (axis == 1) {
        if (j == ng && flow.NY() > 1) sj = j + 1;
        else if (j == flow.NY() + ng - 1 && flow.NY() > 1) sj = j - 1;
    } else if (axis == 2) {
        if (k == ng && flow.NZ() > 1) sk = k + 1;
        else if (k == flow.NZ() + ng - 1 && flow.NZ() > 1) sk = k - 1;
    }
}

template <typename Fn>
inline void forBoundaryGhostsAlongAxis(const Field& flow,
                                       int i,
                                       int j,
                                       int k,
                                       int axis,
                                       Fn&& fn) {
    const int ng = flow.NG();
    if (!Math::isDirectionActiveIndex(axis)) return;

    if (axis == 0) {
        if (i == ng) {
            for (int b = 0; b < ng; ++b) fn(b, j, k, i, j, k);
        }
        if (i == flow.NX() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i + b, j, k, i, j, k);
        }
    } else if (axis == 1) {
        if (j == ng) {
            for (int b = 0; b < ng; ++b) fn(i, b, k, i, j, k);
        }
        if (j == flow.NY() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j + b, k, i, j, k);
        }
    } else if (axis == 2) {
        if (k == ng) {
            for (int b = 0; b < ng; ++b) fn(i, j, b, i, j, k);
        }
        if (k == flow.NZ() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j, k + b, i, j, k);
        }
    }
}

/// @brief 应用标量初始条件。
inline void applyInitialConditions(const Field& flow,
                                   ScalarField& scalar,
                                   const std::vector<BCSetting<double>>& settings) {
    requireCompatible(flow, scalar, "Scalar initial conditions");
    for (const auto& ic : settings) {
        if (ic.type != FIXED_VALUE) {
            std::cerr << "[SF FATAL] scalar initial condition for '"
                      << scalar.name() << "' on set '" << ic.name
                      << "' supports only fixedValue." << std::endl;
            std::exit(1);
        }
        if (!std::isfinite(ic.value)) {
            std::cerr << "[SF FATAL] scalar initial condition for '"
                      << scalar.name() << "' on set '" << ic.name
                      << "' is non-finite: " << ic.value << std::endl;
            std::exit(1);
        }
        const auto& indices = flow.getSet(ic.name);
        for (int idx : indices) {
            int i = 0, j = 0, k = 0;
            flow.getIJK(idx, i, j, k);
            scalar(i, j, k) = ic.value;
        }
    }
}

inline void applyFixedBoundary(const Field& flow,
                               ScalarField& scalar,
                               int i,
                               int j,
                               int k,
                               int axis,
                               double value) {
    auto applyFixed = [&](int gi, int gj, int gk,
                          int ri, int rj, int rk) {
        scalar(gi, gj, gk) = 2.0 * value - scalar(ri, rj, rk);
    };

    if (isGhostCell(flow, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(flow, i, j, k, ri, rj, rk);
        applyFixed(i, j, k, ri, rj, rk);
        return;
    }

    scalar(i, j, k) = value;
    forBoundaryGhostsAlongAxis(flow, i, j, k, axis, applyFixed);
}

inline void applyZeroGradientBoundary(const Field& flow,
                                      ScalarField& scalar,
                                      int i,
                                      int j,
                                      int k,
                                      int axis) {
    if (isGhostCell(flow, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(flow, i, j, k, ri, rj, rk);
        scalar(i, j, k) = scalar(ri, rj, rk);
        return;
    }

    int si = i, sj = j, sk = k;
    nearestInsideSample(flow, i, j, k, axis, si, sj, sk);
    scalar(i, j, k) = scalar(si, sj, sk);

    auto broadcast = [&](int gi, int gj, int gk,
                         int ri, int rj, int rk) {
        scalar(gi, gj, gk) = scalar(ri, rj, rk);
    };
    forBoundaryGhostsAlongAxis(flow, i, j, k, axis, broadcast);
}

inline void applyEmptyBoundary(const Field& flow,
                               ScalarField& scalar,
                               int i,
                               int j,
                               int k,
                               int axis) {
    int di = 0, dj = 0, dk = 0;
    if (axis == 0) di = 1;
    else if (axis == 1) dj = 1;
    else dk = 1;

    auto fill = [&](int gi, int gj, int gk) {
        scalar(gi, gj, gk) = scalar(i, j, k);
    };

    if (axis == 0 && i == flow.NG()) {
        for (int b = 1; b <= flow.NG(); ++b) fill(i - b * di, j, k);
    }
    if (axis == 0 && i == flow.NX() + flow.NG() - 1) {
        for (int b = 1; b <= flow.NG(); ++b) fill(i + b * di, j, k);
    }
    if (axis == 1 && j == flow.NG()) {
        for (int b = 1; b <= flow.NG(); ++b) fill(i, j - b * dj, k);
    }
    if (axis == 1 && j == flow.NY() + flow.NG() - 1) {
        for (int b = 1; b <= flow.NG(); ++b) fill(i, j + b * dj, k);
    }
    if (axis == 2 && k == flow.NG()) {
        for (int b = 1; b <= flow.NG(); ++b) fill(i, j, k - b * dk);
    }
    if (axis == 2 && k == flow.NZ() + flow.NG() - 1) {
        for (int b = 1; b <= flow.NG(); ++b) fill(i, j, k + b * dk);
    }
}

inline SF::Boundary::ILW::BoundaryClosure::BoundaryNormal
scalarILWNormalForAxis(const Field& flow,
                       int i,
                       int j,
                       int k,
                       int axis,
                       const std::string& scalarName,
                       const char* bcName) {
    namespace Closure = SF::Boundary::ILW::BoundaryClosure;
    Closure::BoundaryNormal normal;
    normal.axis = axis;
    const int ng = flow.NG();
    if (axis == 0) {
        if (i <= ng) normal.sign = -1;
        else if (i >= flow.NX() + ng - 1) normal.sign = 1;
    } else if (axis == 1) {
        if (j <= ng) normal.sign = -1;
        else if (j >= flow.NY() + ng - 1) normal.sign = 1;
    } else if (axis == 2) {
        if (k <= ng) normal.sign = -1;
        else if (k >= flow.NZ() + ng - 1) normal.sign = 1;
    }
    if (normal.axis < 0 || normal.axis > 2 || normal.sign == 0) {
        Closure::fatalClosure(
            bcName, i, j, k,
            "scalar '" + scalarName
            + "' boundary point is not on the requested active boundary axis.");
    }
    return normal;
}

inline void nearestScalarILWInteriorCell(const Field& flow,
                                         int i,
                                         int j,
                                         int k,
                                         int& ri,
                                         int& rj,
                                         int& rk) {
    nearestInteriorCell(flow, i, j, k, ri, rj, rk);
}

inline void applyILWFixedBoundary(const Field& flow,
                                  ScalarField& scalar,
                                  int i,
                                  int j,
                                  int k,
                                  int axis,
                                  double value,
                                  int accuracyOrder) {
    namespace Closure = SF::Boundary::ILW::BoundaryClosure;
    const auto normal =
        scalarILWNormalForAxis(flow, i, j, k, axis,
                               scalar.name(), "scalarFixedValue");
    auto getter = [&](const Field&, int gi, int gj, int gk) {
        return scalar(gi, gj, gk);
    };
    auto writeGhost = [&](int gi, int gj, int gk,
                          int ri, int rj, int rk) {
        int refI = ri, refJ = rj, refK = rk;
        nearestScalarILWInteriorCell(flow, ri, rj, rk, refI, refJ, refK);
        std::array<double, Closure::kMaxAccuracyOrder> coeff{};
        Closure::buildMultiDimScalarTaylorCoefficients(
            flow, refI, refJ, refK, normal, accuracyOrder, getter,
            "scalarFixedValue", coeff);
        coeff[0] = value;
        const int layer =
            Closure::ghostLayer(gi, gj, gk, refI, refJ, refK, normal);
        if (layer < 1 || layer > flow.NG()) {
            Closure::fatalClosure("scalarFixedValue", gi, gj, gk,
                                  "invalid ghost layer for scalar ILW fixedValue.");
        }
        const double h = Closure::normalSpacing(flow, refI, refJ, refK, normal);
        scalar(gi, gj, gk) =
            Closure::evaluateTaylor(
                coeff, Closure::taylorOrder(accuracyOrder),
                Closure::ghostDistance(layer, h));
    };

    if (isGhostCell(flow, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestScalarILWInteriorCell(flow, i, j, k, ri, rj, rk);
        writeGhost(i, j, k, ri, rj, rk);
        return;
    }

    scalar(i, j, k) = value;
    forBoundaryGhostsAlongAxis(flow, i, j, k, axis, writeGhost);
}

inline void applyILWZeroGradientBoundary(const Field& flow,
                                         ScalarField& scalar,
                                         int i,
                                         int j,
                                         int k,
                                         int axis,
                                         int accuracyOrder) {
    namespace Closure = SF::Boundary::ILW::BoundaryClosure;
    const auto normal =
        scalarILWNormalForAxis(flow, i, j, k, axis,
                               scalar.name(), "scalarZeroGradient");
    auto getter = [&](const Field&, int gi, int gj, int gk) {
        return scalar(gi, gj, gk);
    };
    auto writeGhost = [&](int gi, int gj, int gk,
                          int ri, int rj, int rk) {
        int refI = ri, refJ = rj, refK = rk;
        nearestScalarILWInteriorCell(flow, ri, rj, rk, refI, refJ, refK);
        std::array<double, Closure::kMaxAccuracyOrder> coeff{};
        Closure::buildMultiDimScalarTaylorCoefficients(
            flow, refI, refJ, refK, normal, accuracyOrder, getter,
            "scalarZeroGradient", coeff);
        coeff[1] = 0.0;
        const int layer =
            Closure::ghostLayer(gi, gj, gk, refI, refJ, refK, normal);
        if (layer < 1 || layer > flow.NG()) {
            Closure::fatalClosure("scalarZeroGradient", gi, gj, gk,
                                  "invalid ghost layer for scalar ILW zeroGradient.");
        }
        const double h = Closure::normalSpacing(flow, refI, refJ, refK, normal);
        scalar(gi, gj, gk) =
            Closure::evaluateTaylor(
                coeff, Closure::taylorOrder(accuracyOrder),
                Closure::ghostDistance(layer, h));
    };

    if (isGhostCell(flow, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestScalarILWInteriorCell(flow, i, j, k, ri, rj, rk);
        writeGhost(i, j, k, ri, rj, rk);
        return;
    }

    forBoundaryGhostsAlongAxis(flow, i, j, k, axis, writeGhost);
}

/// @brief 应用通用标量边界条件。
inline void applyBoundaryConditions(const Field& flow,
                                    ScalarField& scalar,
                                    const std::vector<BCSetting<double>>& settings,
                                    bool ilwEnabled = false,
                                    int ilwAccuracyOrder = 0) {
    requireCompatible(flow, scalar, "Scalar boundary conditions");
    for (const auto& bc : settings) {
        if (bc.type == FIXED_VALUE && !std::isfinite(bc.value)) {
            std::cerr << "[SF FATAL] scalar fixedValue boundary for '"
                      << scalar.name() << "' on set '" << bc.name
                      << "' is non-finite: " << bc.value << std::endl;
            std::exit(1);
        }

        const auto& allSets = flow.getAllSets();
        auto it = allSets.find(bc.name);
        if (it == allSets.end()) {
            std::cerr << "[SF FATAL] scalar boundary set '" << bc.name
                      << "' for '" << scalar.name()
                      << "' is not a mesh set." << std::endl;
            std::exit(1);
        }

        const int axis = (bc.type == EMPTY)
            ? Boundary::Geometry::boundaryAxisForSet(flow, bc.name)
            : Boundary::Geometry::activeBoundaryAxisForSet(flow, bc.name);

        for (int idx : it->second) {
            int i = 0, j = 0, k = 0;
            flow.getIJK(idx, i, j, k);
            if (!Boundary::Geometry::isPhysicalPoint(flow, i, j, k)) continue;

            switch (bc.type) {
            case FIXED_VALUE:
                if (ilwEnabled && ilwAccuracyOrder > 0) {
                    applyILWFixedBoundary(flow, scalar, i, j, k,
                                          axis, bc.value,
                                          ilwAccuracyOrder);
                } else {
                    applyFixedBoundary(flow, scalar, i, j, k, axis, bc.value);
                }
                break;
            case ZERO_GRADIENT:
            case SYMMETRY:
                if (ilwEnabled && ilwAccuracyOrder > 0) {
                    applyILWZeroGradientBoundary(
                        flow, scalar, i, j, k, axis,
                        ilwAccuracyOrder);
                } else {
                    applyZeroGradientBoundary(flow, scalar, i, j, k, axis);
                }
                break;
            case EMPTY:
                applyEmptyBoundary(flow, scalar, i, j, k, axis);
                break;
            }
        }
    }
}

inline double scalarAt(const ScalarField& scalar,
                       const std::vector<double>& values,
                       int i,
                       int j,
                       int k) {
    return values[(size_t)scalar.getIdx(i, j, k)];
}

inline double normalVelocityAtFace(const Field& flow,
                                   int i,
                                   int j,
                                   int k,
                                   Math::Dir dir) {
    int di = 0, dj = 0, dk = 0;
    Math::dirOffset(dir, di, dj, dk);

    const double pL =
        Numerics::requirePhysicalState("scalar transport face left",
                                       flow, i, j, k);
    const double pR =
        Numerics::requirePhysicalState("scalar transport face right",
                                       flow, i + di, j + dj, k + dk);
    (void)pL;
    (void)pR;

    const double rhoL = flow(i, j, k, RHO);
    const double rhoR = flow(i + di, j + dj, k + dk, RHO);
    const double u = 0.5 * (flow(i, j, k, RU) / rhoL
                          + flow(i + di, j + dj, k + dk, RU) / rhoR);
    const double v = 0.5 * (flow(i, j, k, RV) / rhoL
                          + flow(i + di, j + dj, k + dk, RV) / rhoR);
    const double w = 0.5 * (flow(i, j, k, RW) / rhoL
                          + flow(i + di, j + dj, k + dk, RW) / rhoR);

    double metric[4];
    Math::faceMetrics(flow, i, j, k, dir, metric);
    return u * metric[0] + v * metric[1] + w * metric[2];
}

inline double convectiveFaceFlux(const Field& flow,
                                 const ScalarField& scalar,
                                 const std::vector<double>& values,
                                 int i,
                                 int j,
                                 int k,
                                 Math::Dir dir) {
    int di = 0, dj = 0, dk = 0;
    Math::dirOffset(dir, di, dj, dk);
    const double un = normalVelocityAtFace(flow, i, j, k, dir);
    const double left = scalarAt(scalar, values, i, j, k);
    const double right = scalarAt(scalar, values, i + di, j + dj, k + dk);
    requireFiniteScalar("scalar transport face left", left, i, j, k);
    requireFiniteScalar("scalar transport face right",
                        right, i + di, j + dj, k + dk);
    return un >= 0.0 ? un * left : un * right;
}

inline double convectiveDivergence(const Field& flow,
                                   const ScalarField& scalar,
                                   const std::vector<double>& values,
                                   int i,
                                   int j,
                                   int k) {
    double div = 0.0;
    if (Math::isDirectionActive(Math::XI)) {
        div += convectiveFaceFlux(flow, scalar, values, i, j, k, Math::XI)
             - convectiveFaceFlux(flow, scalar, values, i - 1, j, k, Math::XI);
    }
    if (Math::isDirectionActive(Math::ETA)) {
        div += convectiveFaceFlux(flow, scalar, values, i, j, k, Math::ETA)
             - convectiveFaceFlux(flow, scalar, values, i, j - 1, k, Math::ETA);
    }
    if (Math::isDirectionActive(Math::ZETA)) {
        div += convectiveFaceFlux(flow, scalar, values, i, j, k, Math::ZETA)
             - convectiveFaceFlux(flow, scalar, values, i, j, k - 1, Math::ZETA);
    }
    return flow.Jac(i, j, k) * div;
}

inline double centralDiffusion(const Field& flow,
                               const ScalarField& scalar,
                               const std::vector<double>& values,
                               int i,
                               int j,
                               int k,
                               double diffusivity) {
    if (!std::isfinite(diffusivity) || diffusivity < 0.0) {
        std::cerr << "[SF FATAL] scalar diffusivity for '" << scalar.name()
                  << "' must be finite and >= 0, got " << diffusivity
                  << std::endl;
        std::exit(1);
    }
    if (diffusivity == 0.0) return 0.0;

    double result = 0.0;
    const double Jcell = flow.Jac(i, j, k);

    auto addDir = [&](Math::Dir dir, int im, int jm, int km,
                      int ip, int jp, int kp) {
        double cfL[4], cfR[4];
        if (dir == Math::XI) {
            Math::faceMetrics(flow, i - 1, j, k, dir, cfL);
            Math::faceMetrics(flow, i, j, k, dir, cfR);
        } else if (dir == Math::ETA) {
            Math::faceMetrics(flow, i, j - 1, k, dir, cfL);
            Math::faceMetrics(flow, i, j, k, dir, cfR);
        } else {
            Math::faceMetrics(flow, i, j, k - 1, dir, cfL);
            Math::faceMetrics(flow, i, j, k, dir, cfR);
        }

        if (!std::isfinite(cfL[3]) || !std::isfinite(cfR[3])
            || std::abs(cfL[3]) <= 1.0e-30
            || std::abs(cfR[3]) <= 1.0e-30) {
            std::cerr << "[SF FATAL] scalar diffusion found invalid face "
                      << "Jacobian for '" << scalar.name() << "' at ("
                      << i << "," << j << "," << k << ")" << std::endl;
            std::exit(1);
        }
        const double JgL =
            (cfL[0] * cfL[0] + cfL[1] * cfL[1] + cfL[2] * cfL[2])
            / cfL[3];
        const double JgR =
            (cfR[0] * cfR[0] + cfR[1] * cfR[1] + cfR[2] * cfR[2])
            / cfR[3];
        if (!std::isfinite(JgL) || !std::isfinite(JgR)) {
            std::cerr << "[SF FATAL] scalar diffusion metric is non-finite for '"
                      << scalar.name() << "' at (" << i << "," << j
                      << "," << k << ")" << std::endl;
            std::exit(1);
        }

        const double center = scalarAt(scalar, values, i, j, k);
        const double left = scalarAt(scalar, values, im, jm, km);
        const double right = scalarAt(scalar, values, ip, jp, kp);
        result += diffusivity * Jcell
                * (JgR * (right - center) - JgL * (center - left));
    };

    if (Math::isDirectionActive(Math::XI)) {
        addDir(Math::XI, i - 1, j, k, i + 1, j, k);
    }
    if (Math::isDirectionActive(Math::ETA)) {
        addDir(Math::ETA, i, j - 1, k, i, j + 1, k);
    }
    if (Math::isDirectionActive(Math::ZETA)) {
        addDir(Math::ZETA, i, j, k - 1, i, j, k + 1);
    }
    return result;
}

/// @brief 只装配标量输运 RHS，不修改标量时间层。
/// @param flow 当前 stage 的流动守恒状态。
/// @param scalar 当前 stage 的守恒标量。
/// @param config 对流、扩散与常源配置。
/// @param rhs 输出 d(scalar)/dt，尺寸与 Field 一致。
inline void assembleRHS(const Field& flow,
                        const ScalarField& scalar,
                        const TransportConfig& config,
                        std::vector<double>& rhs) {
    requireCompatible(flow, scalar, "Scalar RHS");
    rhs.assign((size_t)flow.TotalSize(), 0.0);
    const std::vector<double>& values = scalar.values();
    Math::forInterior(flow, [&](int i, int j, int k) {
        if (flow.CellFlag(i, j, k) != FLUID_CELL) return;
        double value = 0.0;
        if (config.convectionEnabled) {
            value -= convectiveDivergence(flow, scalar, values, i, j, k);
        }
        if (config.diffusionEnabled) {
            value += centralDiffusion(flow, scalar, values, i, j, k,
                                      config.diffusivity);
        }
        if (config.sourceEnabled) {
            if (!std::isfinite(config.source)) {
                std::cerr << "[SF FATAL] scalar source for '"
                          << scalar.name() << "' is non-finite: "
                          << config.source << std::endl;
                std::exit(1);
            }
            value += config.source;
        }
        requireFiniteScalar("scalar RHS", value, i, j, k);
        rhs[(size_t)scalar.getIdx(i, j, k)] = value;
    });
}

/// @brief 检查或投影标量上下界。
inline BoundReport applyBounds(const Field& flow,
                               ScalarField& scalar,
                               const Bounds& bounds) {
    BoundReport report;
    if (!bounds.enabled || bounds.mode == BoundMode::None) return report;
    if (!std::isfinite(bounds.lower) || !std::isfinite(bounds.upper)
        || bounds.lower > bounds.upper) {
        std::cerr << "[SF FATAL] invalid scalar bounds for '"
                  << scalar.name() << "': [" << bounds.lower
                  << ", " << bounds.upper << "]" << std::endl;
        std::exit(1);
    }

    Math::forInterior(flow, [&](int i, int j, int k) {
        if (flow.CellFlag(i, j, k) != FLUID_CELL) return;
        double& value = scalar(i, j, k);
        requireFiniteScalar("scalar bounds", value, i, j, k);
        report.sampled = true;
        report.minValue = std::min(report.minValue, value);
        report.maxValue = std::max(report.maxValue, value);
        if (value < bounds.lower) {
            ++report.below;
            if (bounds.mode == BoundMode::Project) value = bounds.lower;
        } else if (value > bounds.upper) {
            ++report.above;
            if (bounds.mode == BoundMode::Project) value = bounds.upper;
        }
    });

    if (report.violated()) {
        SF::broadcast("Scalar bounds     : ",
                      scalar.name() + " mode=" + toString(bounds.mode)
                      + ", below=" + std::to_string(report.below)
                      + ", above=" + std::to_string(report.above)
                      + ", min=" + std::to_string(report.minValue)
                      + ", max=" + std::to_string(report.maxValue));
    }
    return report;
}

/// @brief 一阶显式推进通用标量输运方程。
inline void advanceForwardEuler(const Field& flow,
                                ScalarField& scalar,
                                const std::vector<BCSetting<double>>& boundaries,
                                const TransportConfig& config,
                                double dt) {
    requireCompatible(flow, scalar, "Scalar transport");
    if (!std::isfinite(dt) || dt <= 0.0) {
        std::cerr << "[SF FATAL] scalar transport for '" << scalar.name()
                  << "' requires finite positive dt, got " << dt
                  << std::endl;
        std::exit(1);
    }

    applyBoundaryConditions(
        flow, scalar, boundaries,
        config.ilwBoundaryEnabled, config.ilwAccuracyOrder);
    applyBounds(flow, scalar, config.bounds);

    const std::vector<double> oldValues = scalar.values();
    std::vector<double> nextValues = oldValues;

    Math::forInterior(flow, [&](int i, int j, int k) {
        if (flow.CellFlag(i, j, k) != FLUID_CELL) return;
        const int id = scalar.getIdx(i, j, k);
        const double oldValue = oldValues[(size_t)id];
        requireFiniteScalar("scalar transport old value", oldValue, i, j, k);

        double rhs = 0.0;
        if (config.convectionEnabled) {
            rhs -= convectiveDivergence(flow, scalar, oldValues, i, j, k);
        }
        if (config.diffusionEnabled) {
            rhs += centralDiffusion(flow, scalar, oldValues, i, j, k,
                                    config.diffusivity);
        }
        if (config.sourceEnabled) {
            if (!std::isfinite(config.source)) {
                std::cerr << "[SF FATAL] scalar source for '" << scalar.name()
                          << "' is non-finite: " << config.source << std::endl;
                std::exit(1);
            }
            rhs += config.source;
        }

        const double updated = oldValue + dt * rhs;
        requireFiniteScalar("scalar transport updated value",
                            updated, i, j, k);
        nextValues[(size_t)id] = updated;
    });

    scalar.values().swap(nextValues);
    applyBoundaryConditions(
        flow, scalar, boundaries,
        config.ilwBoundaryEnabled, config.ilwAccuracyOrder);
    const BoundReport report = applyBounds(flow, scalar, config.bounds);
    if (report.violated() && config.bounds.mode == BoundMode::Project) {
        applyBoundaryConditions(
            flow, scalar, boundaries,
            config.ilwBoundaryEnabled, config.ilwAccuracyOrder);
    }
}

} // namespace Scalar
} // namespace FDM
} // namespace SF
