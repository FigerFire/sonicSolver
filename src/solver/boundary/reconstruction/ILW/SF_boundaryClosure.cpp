/// @file SF_boundaryClosure.cpp
/// @brief ILW 特征边界方程与 Taylor ghost 状态闭合。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_boundaryClosure.h"

#include "core/mesh/SF_dimension.h"
#include "methods/math/discrete/SF_polynomial.h"
#include "SF_taylor.h"
#include "SF_numericsPolicy.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace SF {
namespace Boundary {
namespace ILW {
namespace BoundaryClosure {
namespace {

constexpr double kWenoEpsilon = 1.0e-6;
constexpr double kWenoPower = 3.0;

bool supportedAccuracy(int accuracyOrder) {
    return accuracyOrder == 3 || accuracyOrder == 5 ||
           accuracyOrder == 7 || accuracyOrder == 9;
}

double wenoWeightScale(double h) {
    if (!std::isfinite(h) || h <= 0.0) return 1.0e-8;
    return std::min(0.5, std::max(1.0e-8, h));
}

} // namespace

int requireAccuracyOrder(int accuracyOrder) {
    if (!FDM::isSupportedILWOrder(accuracyOrder) || accuracyOrder == 0) {
        std::cerr << "[SF FATAL] ILW physical boundary closure requires "
                  << "[numerics].ILW = 3, 5, 7, or 9; got "
                  << accuracyOrder << "." << std::endl;
        std::exit(1);
    }
    return accuracyOrder;
}

int taylorOrder(int accuracyOrder) {
    return FDM::ilwTaylorOrderFromAccuracy(
        requireAccuracyOrder(accuracyOrder));
}

BoundaryNormal boundaryNormal(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    BoundaryNormal normal;

    if (Math::isDirectionActiveIndex(0) &&
        (i < ng || i >= field.NX() + ng)) {
        normal.axis = 0;
        normal.sign = (i < ng) ? -1 : 1;
        return normal;
    }
    if (Math::isDirectionActiveIndex(1) &&
        (j < ng || j >= field.NY() + ng)) {
        normal.axis = 1;
        normal.sign = (j < ng) ? -1 : 1;
        return normal;
    }
    if (Math::isDirectionActiveIndex(2) &&
        (k < ng || k >= field.NZ() + ng)) {
        normal.axis = 2;
        normal.sign = (k < ng) ? -1 : 1;
        return normal;
    }

    if (Math::isDirectionActiveIndex(0) &&
        (i == ng || i == field.NX() + ng - 1)) {
        normal.axis = 0;
        normal.sign = (i == ng) ? -1 : 1;
        return normal;
    }
    if (Math::isDirectionActiveIndex(1) &&
        (j == ng || j == field.NY() + ng - 1)) {
        normal.axis = 1;
        normal.sign = (j == ng) ? -1 : 1;
        return normal;
    }
    if (Math::isDirectionActiveIndex(2) &&
        (k == ng || k == field.NZ() + ng - 1)) {
        normal.axis = 2;
        normal.sign = (k == ng) ? -1 : 1;
        return normal;
    }

    return normal;
}

double normalSpacing(const Field& field,
                     int i, int j, int k,
                     const BoundaryNormal& normal) {
    int ni = i;
    int nj = j;
    int nk = k;
    const int ng = field.NG();
    if (normal.axis == 0) {
        if (normal.sign > 0 && i > ng) {
            ni = i - 1;
        } else {
            ni = i + 1;
        }
    } else if (normal.axis == 1) {
        if (normal.sign > 0 && j > ng) {
            nj = j - 1;
        } else {
            nj = j + 1;
        }
    } else if (normal.axis == 2) {
        if (normal.sign > 0 && k > ng) {
            nk = k - 1;
        } else {
            nk = k + 1;
        }
    }

    if (!Boundary::isPhysicalCell(field, ni, nj, nk)) {
        fatalClosure("BoundaryClosure", i, j, k,
                     "normal grid-spacing neighbor is outside the physical domain.");
    }

    const double dx = field.X(ni, nj, nk) - field.X(i, j, k);
    const double dy = field.Y(ni, nj, nk) - field.Y(i, j, k);
    const double dz = field.Z(ni, nj, nk) - field.Z(i, j, k);
    const double h = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (!std::isfinite(h) || h <= 1.0e-14) {
        fatalClosure("BoundaryClosure", i, j, k,
                     "invalid normal grid spacing on physical boundary.");
    }
    return h;
}

int ghostLayer(int gi, int gj, int gk,
               int ri, int rj, int rk,
               const BoundaryNormal& normal) {
    if (normal.axis == 0) return std::abs(gi - ri);
    if (normal.axis == 1) return std::abs(gj - rj);
    if (normal.axis == 2) return std::abs(gk - rk);
    return 0;
}

double ghostDistance(int layer, double h) {
    return ((double)layer - 0.5) * h;
}

namespace {

/// @brief 从样本中估算局部坐标缩放因子。
double estimateCoordScale(
    const std::vector<std::pair<double, Math::Polynomial::TensorPoint>>& samples,
    int coordIndex) {
    std::vector<double> coords;
    coords.reserve(samples.size());
    for (const auto& sample : samples) {
        double v = (coordIndex == 0) ? sample.second.s
                 : (coordIndex == 1) ? sample.second.a
                 : sample.second.b;
        if (std::isfinite(v)) coords.push_back(v);
    }
    if (coords.size() < 2) return 1.0;

    std::sort(coords.begin(), coords.end());
    std::vector<double> gaps;
    gaps.reserve(coords.size() - 1);
    double previous = coords.front();
    for (size_t n = 1; n < coords.size(); ++n) {
        const double current = coords[n];
        const double gap = current - previous;
        if (gap > 1.0e-14) {
            gaps.push_back(gap);
            previous = current;
        }
    }

    if (gaps.empty()) return 1.0;
    std::sort(gaps.begin(), gaps.end());
    return std::max(1.0e-14, gaps[gaps.size() / 2]);
}

/// @brief 为指定阶数的子模板收集样本。
std::vector<Math::Polynomial::TensorPoint> gatherSubStencilPoints(
    const std::vector<std::pair<double, Math::Polynomial::TensorPoint>>& samples,
    int polyOrder,
    double sScale,
    double aScale,
    double bScale,
    bool useT2) {
    std::vector<Math::Polynomial::TensorPoint> pts;
    pts.reserve(samples.size());

    /// 子模板的法向范围为前 (polyOrder+1) 层，切向为中心 ±polyOrder/2。
    const double sMax = (polyOrder + 0.5 + 1e-6) * sScale;
    const double aMax = (polyOrder * 0.5 + 1e-6) * aScale;
    const double bMax = (polyOrder * 0.5 + 1e-6) * bScale;

    for (const auto& sample : samples) {
        const double absS = std::abs(sample.second.s);
        const double absA = std::abs(sample.second.a);
        const double absB = std::abs(sample.second.b);
        if (absS > sMax) continue;
        if (absA > aMax) continue;
        if (useT2 && absB > bMax) continue;
        pts.push_back(sample.second);
    }
    return pts;
}

/// @brief 在指定局部坐标处求多项式导数值。
double evaluatePolynomialDerivative(
    const std::vector<double>& coeff,
    int polyOrder,
    bool useT2,
    int ds, int da, int db,
    double s, double a, double b,
    double sScale, double aScale, double bScale) {
    double value = 0.0;
    int idx = 0;
    for (int total = 0; total <= polyOrder; ++total) {
        for (int ps = total; ps >= 0; --ps) {
            if (useT2) {
                for (int pa = total - ps; pa >= 0; --pa) {
                    const int pb = total - ps - pa;
                    double term = coeff[(size_t)idx];
                    if (ps >= ds && pa >= da && pb >= db) {
                        const double denom =
                            Math::Taylor::factorial(ps - ds)
                          * Math::Taylor::factorial(pa - da)
                          * Math::Taylor::factorial(pb - db);
                        term *= Math::Polynomial::powInt(s / sScale, ps - ds)
                              * Math::Polynomial::powInt(a / aScale, pa - da)
                              * Math::Polynomial::powInt(b / bScale, pb - db)
                              / denom;
                    } else {
                        term = 0.0;
                    }
                    value += term;
                    ++idx;
                }
            } else {
                const int pa = total - ps;
                double term = coeff[(size_t)idx];
                if (ps >= ds && pa >= da) {
                    const double denom =
                        Math::Taylor::factorial(ps - ds)
                      * Math::Taylor::factorial(pa - da);
                    term *= Math::Polynomial::powInt(s / sScale, ps - ds)
                          * Math::Polynomial::powInt(a / aScale, pa - da)
                          / denom;
                } else {
                    term = 0.0;
                }
                value += term;
                ++idx;
            }
        }
    }
    return value;
}

/// @brief 从拟合多项式系数中提取法向导数。
double extractNormalDerivative(
    const std::vector<double>& coeff,
    int polyOrder,
    int derivativeOrder,
    bool useT2,
    double sScale) {
    const int idx = Math::Polynomial::normalDerivativeIndex(
        useT2, polyOrder, derivativeOrder);
    if (idx < 0 || idx >= (int)coeff.size()) return 0.0;

    const double invScale = 1.0 / std::pow(std::max(1e-14, sScale), derivativeOrder);
    return coeff[(size_t)idx] * invScale;
}

/// @brief 对多项式系数计算多维光滑指示子。
///
/// 使用 4×4（2D）或 4×4×4（3D）Gauss 积分计算各阶导数平方的加权和。
double multiDimSmoothnessIndicator(
    const std::vector<double>& coeff,
    int polyOrder,
    bool useT2,
    double sScale,
    double aScale,
    double bScale) {
    if (polyOrder <= 0) {
        return 2.0 * sScale * sScale;
    }

    constexpr double gp[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526
    };
    constexpr double gw[4] = {
        0.3478548451374538, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374538
    };

    const double halfS = 0.5 * sScale;
    const double halfA = 0.5 * aScale;
    const double halfB = 0.5 * bScale;
    const double jac3D = halfS * halfA * halfB;
    const double jac2D = halfS * halfA;

    double beta = 0.0;

    for (int total = 1; total <= polyOrder; ++total) {
        const double scale = std::pow(sScale * sScale, total - 1);
        for (int ds = total; ds >= 0; --ds) {
            const int taRemaining = total - ds;
            if (useT2) {
                for (int da = taRemaining; da >= 0; --da) {
                    const int db = taRemaining - da;
                    double integral = 0.0;
                    for (int qs = 0; qs < 4; ++qs) {
                        const double sv = halfS * gp[qs];
                        for (int qa = 0; qa < 4; ++qa) {
                            const double av = halfA * gp[qa];
                            for (int qb = 0; qb < 4; ++qb) {
                                const double bv = halfB * gp[qb];
                                double dv = evaluatePolynomialDerivative(
                                    coeff, polyOrder, useT2,
                                    ds, da, db, sv, av, bv, sScale, aScale, bScale);
                                integral += gw[qs] * gw[qa] * gw[qb] * dv * dv;
                            }
                        }
                    }
                    beta += scale * integral;
                }
            } else {
                const int da = taRemaining;
                double integral = 0.0;
                for (int qs = 0; qs < 4; ++qs) {
                    const double sv = halfS * gp[qs];
                    for (int qa = 0; qa < 4; ++qa) {
                        const double av = halfA * gp[qa];
                        double dv = evaluatePolynomialDerivative(
                            coeff, polyOrder, useT2,
                            ds, da, 0, sv, av, 0.0, sScale, aScale, 1.0);
                        integral += gw[qs] * gw[qa] * dv * dv;
                    }
                }
                beta += scale * integral;
            }
        }
    }

    if (useT2) beta *= jac3D;
    else beta *= jac2D;

    return std::isfinite(beta) ? beta : 1.0e30;
}

/// @brief 从子模板样本拟合多项式。
bool fitPolynomialToSubStencil(
    const std::vector<std::pair<double, Math::Polynomial::TensorPoint>>& samples,
    const std::vector<Math::Polynomial::TensorPoint>& subPoints,
    int polyOrder,
    bool useT2,
    double sScale,
    double aScale,
    double bScale,
    std::vector<double>& coeff) {
    /// 缩放坐标以改善条件数。
    std::vector<Math::Polynomial::TensorPoint> scaledPoints;
    scaledPoints.reserve(subPoints.size());
    for (const auto& pt : subPoints) {
        scaledPoints.push_back({
            pt.s / sScale,
            pt.a / aScale,
            pt.b / (useT2 ? bScale : 1.0)
        });
    }

    Math::Polynomial::LeastSquaresPlan plan;
    if (!Math::Polynomial::buildLeastSquaresPlan(
            scaledPoints, useT2, polyOrder, plan)) {
        return false;
    }

    /// 收集对应样本值。
    std::vector<double> values;
    values.reserve(subPoints.size());
    for (const auto& pt : subPoints) {
        /// 在原始样本中查找匹配的坐标。
        bool found = false;
        for (const auto& sample : samples) {
            if (std::abs(sample.second.s - pt.s) < 1e-12
             && std::abs(sample.second.a - pt.a) < 1e-12
             && std::abs(sample.second.b - pt.b) < 1e-12) {
                values.push_back(sample.first);
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    const int nBasis = (int)plan.projectionRows.size();
    coeff.assign((size_t)nBasis, 0.0);
    for (int c = 0; c < nBasis; ++c) {
        double sum = 0.0;
        for (size_t i = 0; i < values.size(); ++i) {
            sum += plan.projectionRows[(size_t)c][i] * values[i];
        }
        coeff[(size_t)c] = sum;
    }
    return true;
}

} // namespace

bool extrapolateMultiDim(
    const std::vector<std::pair<double, Math::Polynomial::TensorPoint>>& samples,
    int accuracyOrder,
    bool useT2,
    std::array<double, kMaxAccuracyOrder>& coefficients) {
    coefficients.fill(0.0);
    if (!supportedAccuracy(accuracyOrder)) return false;
    if (samples.empty()) return false;

    const int N = accuracyOrder;
    const double sScale = estimateCoordScale(samples, 0);
    const double aScale = estimateCoordScale(samples, 1);
    const double bScale = useT2
        ? estimateCoordScale(samples, 2)
        : 1.0;
    const double hEst = std::max(1e-14, sScale);

    /// 线性权重：低阶模板权重按 h^(N-1-r) 缩放。
    std::array<double, kMaxAccuracyOrder> linearWeights{};
    const double hw = wenoWeightScale(hEst);
    double lowerWeightSum = 0.0;
    for (int r = 0; r < N - 1; ++r) {
        linearWeights[(size_t)r] = std::pow(hw, N - 1 - r);
        lowerWeightSum += linearWeights[(size_t)r];
    }
    linearWeights[(size_t)(N - 1)] =
        std::max(1.0e-12, 1.0 - lowerWeightSum);

    /// candidateDerivs[d][r] = d 阶导数的第 r 个候选（未加权）。
    std::array<std::array<double, kMaxAccuracyOrder>, kMaxAccuracyOrder>
        candidateDerivs{};
    for (auto& row : candidateDerivs) row.fill(0.0);

    std::array<double, kMaxAccuracyOrder> alphaVals{};
    double alphaSum = 0.0;

    for (int r = 0; r < N; ++r) {
        std::vector<Math::Polynomial::TensorPoint> subPoints =
            gatherSubStencilPoints(samples, r, sScale, aScale, bScale, useT2);
        const int minPts = Math::Polynomial::tensorBasisSize(useT2, r);
        if ((int)subPoints.size() < minPts) {
            alphaVals[(size_t)r] = 0.0;
            continue;
        }

        std::vector<double> coeff;
        if (!fitPolynomialToSubStencil(samples, subPoints, r, useT2,
                                       sScale, aScale, bScale, coeff)) {
            alphaVals[(size_t)r] = 0.0;
            continue;
        }

        const double beta = multiDimSmoothnessIndicator(
            coeff, r, useT2, sScale, aScale, bScale);
        const double denom = std::pow(kWenoEpsilon + beta, kWenoPower);
        alphaVals[(size_t)r] = linearWeights[(size_t)r] / denom;
        if (!std::isfinite(alphaVals[(size_t)r]) || alphaVals[(size_t)r] < 0.0) {
            alphaVals[(size_t)r] = 0.0;
        }
        alphaSum += alphaVals[(size_t)r];

        /// 暂存各阶导数的原始（未加权）候选。
        for (int d = 0; d <= r; ++d) {
            candidateDerivs[(size_t)d][(size_t)r] =
                extractNormalDerivative(coeff, r, d, useT2, sScale);
        }
    }

    if (!(alphaSum > 0.0) || !std::isfinite(alphaSum)) return false;

    std::array<double, kMaxAccuracyOrder> omega{};
    for (int r = 0; r < N; ++r) {
        omega[(size_t)r] = alphaVals[(size_t)r] / alphaSum;
    }

    for (int d = 0; d < N; ++d) {
        double value = 0.0;
        for (int r = d; r < N; ++r) {
            value += omega[(size_t)r] * candidateDerivs[(size_t)d][(size_t)r];
        }
        /// 除以 d! 存储在 coefficients[d] 中，便于直接做 Taylor 展开。
        coefficients[(size_t)d] = value / Math::Taylor::factorial(d);
    }

    return true;
}

double evaluateTaylor(const std::array<double, kMaxAccuracyOrder>& coefficients,
                      int taylorOrder,
                      double distance) {
    return Math::Taylor::evaluateScaled1D(
        coefficients, distance, std::min(taylorOrder, kMaxTaylorOrder));
}

void fatalClosure(const char* bcName,
                  int i, int j, int k,
                  const std::string& reason) {
    std::cerr << "[SF FATAL] ILW " << bcName
              << " boundary closure failed at ("
              << i << "," << j << "," << k << "): "
              << reason << std::endl;
    std::exit(1);
}

} // namespace BoundaryClosure
} // namespace ILW
} // namespace Boundary
} // namespace SF
