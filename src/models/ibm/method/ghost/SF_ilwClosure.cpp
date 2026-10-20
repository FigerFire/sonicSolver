/// @file SF_ilwClosure.cpp
/// @brief ghost IBM 的 ILW 样本预处理与高阶 ghost 状态闭合。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#include "method/ghost/SF_ilwClosure.h"

#include "SF_cellType.h"
#include "core/mesh/SF_dimension.h"
#include "SF_eulerState.h"
#include "SF_localFrame.h"
#include "methods/math/discrete/SF_polynomial.h"
#include "SF_slipWall.h"
#include "SF_surfaceCurvature.h"
#include "SF_taylor.h"
#include "SF_utility.h"
#include "SF_weightedProjection.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace SF {
namespace IBM {
namespace GhostILW {

// ═══════════════════════════════════════════════════════════════
//  IBM ILW — 匿名空间：内部辅助函数
// ═══════════════════════════════════════════════════════════════

namespace {

constexpr int kMaxSupportedTaylorOrder = 8;
constexpr double kWenoExtrapolationPower = 3.0;
constexpr double kWenoExtrapolationEpsilon = 1.0e-6;
constexpr int kMaxStructured2DVariantsPerOrder = 6;
constexpr int kMaxStructured3DVariantsPerOrder = 8;
constexpr int kMaxStructured3DWindowStartsPerAxis = 10;
constexpr int kMaxStructured3DTemplatesExamined =
    3 * kMaxStructured3DVariantsPerOrder;
constexpr int kMaxRanked3DVariantsPerOrder = 5;

using Primitive = Math::Euler::Primitive;
using LocalFrame = Math::LocalFrame;
using Clock = std::chrono::steady_clock;
using namespace Math::Euler;
using Math::localCoordinates;
using Math::makeLocalFrame;
using Math::makeLocalFrame2D;
using Math::signedDistanceFromWall;

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

bool isIBMReadableCell(const Field& field, int i, int j, int k) {
    return SF::IBM::isIBMReadableCell(field, i, j, k);
}

bool isIbmNonFluidCell(const Field& field, int i, int j, int k) {
    return SF::IBM::isIbmNonFluidCell(field, i, j, k);
}

void loadConservative(const Field& field, int i, int j, int k, double q[5]) {
    q[0] = field(i, j, k, RHO);
    q[1] = field(i, j, k, RU);
    q[2] = field(i, j, k, RV);
    q[3] = field(i, j, k, RW);
    q[4] = field(i, j, k, E);
}

using CurvatureTensor = Math::SurfaceCurvature::Tensor;

struct CurvatureFilterStats {
    int cached = 0;
    int skippedNonPhysical = 0;
    int skippedMissingGeometry = 0;
    int skippedAlignment = 0;
    int skippedCoincident = 0;
    int accepted = 0;
    double minAbsAlignment = std::numeric_limits<double>::max();
    double maxAbsAlignment = 0.0;
    double minTangentialR2 = std::numeric_limits<double>::max();
};

bool ilwFail(const char* reason) {
    static int emitted = 0;
    if (emitted < 20) {
        std::cerr << "[SF ILW] " << reason << std::endl;
    } else if (emitted == 20) {
        std::cerr << "[SF ILW] further ILW diagnostics suppressed" << std::endl;
    }
    ++emitted;
    return false;
}

void printCurvatureFilterDiagnostic(int ibmI, int ibmJ, int ibmK,
                                    const CurvatureFilterStats& stats,
                                    const double wallPoint[3],
                                    const LocalFrame& frame,
                                    const IBMRuntimeConfig& config,
                                    bool useT2) {
    static int emitted = 0;
    if (emitted >= 12) return;
    ++emitted;

    std::cerr << "[SF ILW] curvature normal filter at ("
              << ibmI << "," << ibmJ << "," << ibmK
              << "): cached=" << stats.cached
              << ", accepted=" << stats.accepted
              << ", skippedNonPhysical=" << stats.skippedNonPhysical
              << ", skippedMissingGeometry=" << stats.skippedMissingGeometry
              << ", skippedAlignment=" << stats.skippedAlignment
              << ", skippedCoincident=" << stats.skippedCoincident
              << ", useT2=" << useT2
              << ", normalAngleLimit=" << config.normalAngleDegrees
              << "deg"
              << ", alignmentThreshold="
              << config.normalAlignmentThreshold();
    if (stats.minAbsAlignment != std::numeric_limits<double>::max()) {
        std::cerr << ", absAlignmentRange=["
                  << stats.minAbsAlignment << ","
                  << stats.maxAbsAlignment << "]";
    }
    if (stats.minTangentialR2 != std::numeric_limits<double>::max()) {
        std::cerr << ", minTangentialR2=" << stats.minTangentialR2;
    }
    std::cerr << ", wallPoint=("
              << wallPoint[0] << "," << wallPoint[1] << ","
              << wallPoint[2] << ")"
              << ", frame.n=(" << frame.n[0] << "," << frame.n[1]
              << "," << frame.n[2] << ")"
              << std::endl;
}

bool buildFluidVisibleBoundaryGeometry(const Field& field,
                                       const IBMGeometry& geometry,
                                       int ibmI, int ibmJ, int ibmK,
                                       double wallPoint[3],
                                       double exteriorNormal[3],
                                       double wallVelocity[3],
                                       double& targetDistance);

void cellPoint(const Field& field, int i, int j, int k, double p[3]) {
    p[0] = field.X(i, j, k);
    p[1] = field.Y(i, j, k);
    p[2] = field.Z(i, j, k);
}

bool estimateCurvatureTensor(const Field& field,
                             const Storage& storage,
                             const IBMGeometry& geometry,
                             int ibmI, int ibmJ, int ibmK,
                             const double wallPoint[3],
                             const LocalFrame& frame,
                             const IBMRuntimeConfig& config,
                             bool useT2,
                             CurvatureTensor& curvature) {
    std::vector<Math::SurfaceCurvature::NormalSample> normals;
    normals.reserve(64);

    const size_t cell = (size_t)field.getIdx(ibmI, ibmJ, ibmK);
    const int nCached = storage.normalCount(cell);
    const double alignmentThreshold = config.normalAlignmentThreshold();
    CurvatureFilterStats filterStats;
    filterStats.cached = nCached;
    if (nCached <= 0) {
        printCurvatureFilterDiagnostic(ibmI, ibmJ, ibmK, filterStats,
                                       wallPoint, frame, config, useT2);
        return ilwFail("curvature tensor has no cached IBM normals");
    }

    for (int c = 0; c < nCached; ++c) {
        int ii = 0, jj = 0, kk = 0;
        field.getIJK(storage.normalSample(cell, c),
                     ii, jj, kk);
        if (!isIBMReadableCell(field, ii, jj, kk)) {
            ++filterStats.skippedNonPhysical;
            continue;
        }
        if (!geometry.hasGeometry(ii, jj, kk)) {
            ++filterStats.skippedMissingGeometry;
            continue;
        }

        const Vector3 wp = geometry.wallPoint(ii, jj, kk);
        const Vector3 wn = geometry.wallNormal(ii, jj, kk);
        double p[3] = {wp.x, wp.y, wp.z};
        double n[3] = {-wn.x, -wn.y, -wn.z};
        Math::normalizeArray(n);
        const Math::Vector3 normalVector = Math::vectorFromArray(n);
        double normalAlignment = Math::dot(normalVector, frame.n);
        const double absAlignment = std::abs(normalAlignment);
        filterStats.minAbsAlignment =
            std::min(filterStats.minAbsAlignment, absAlignment);
        filterStats.maxAbsAlignment =
            std::max(filterStats.maxAbsAlignment, absAlignment);
        if (absAlignment < alignmentThreshold) {
            ++filterStats.skippedAlignment;
            continue;
        }
        if (normalAlignment < 0.0) {
            n[0] = -n[0];
            n[1] = -n[1];
            n[2] = -n[2];
            normalAlignment = -normalAlignment;
        }

        double s = 0.0, a = 0.0, b = 0.0;
        localCoordinates(p, wallPoint, frame, s, a, b);
        const Math::Vector3 dn = Math::vectorFromArray(n) - frame.n;
        const double dn1 = Math::dot(dn, frame.t1);
        const double dn2 = Math::dot(dn, frame.t2);
        const double r2 = a * a + b * b;
        filterStats.minTangentialR2 =
            std::min(filterStats.minTangentialR2, r2);
        if (r2 <= 1e-20) {
            ++filterStats.skippedCoincident;
            continue;
        }
        ++filterStats.accepted;
        normals.push_back({a, b, dn1, dn2, r2});
    }

    if (normals.empty()) {
        printCurvatureFilterDiagnostic(ibmI, ibmJ, ibmK, filterStats,
                                       wallPoint, frame, config, useT2);
        return ilwFail("curvature tensor has no neighboring IBM normals");
    }
    std::sort(normals.begin(), normals.end(),
              [](const Math::SurfaceCurvature::NormalSample& lhs,
                 const Math::SurfaceCurvature::NormalSample& rhs) {
                  return lhs.r2 < rhs.r2;
              });
    if ((int)normals.size() > 64) normals.resize(64);

    if (!Math::SurfaceCurvature::fitTensor(normals, useT2, curvature)) {
        return ilwFail(useT2
            ? "curvature tensor 3D fit is singular"
            : "curvature tensor 2D fit is singular");
    }
    return true;
}

struct PlanSampleGeometry {
    int idx = 0;
    double s = 0.0;
    double a = 0.0;
    double b = 0.0;
    double r2 = 0.0;
};

enum class PlanFailure {
    None,
    FluidSamples,
    Curvature,
    Fit
};

const char* planFailureName(PlanFailure failure) {
    switch (failure) {
        case PlanFailure::None: return "None";
        case PlanFailure::FluidSamples: return "FluidSamples";
        case PlanFailure::Curvature: return "Curvature";
        case PlanFailure::Fit: return "Fit";
    }
    return "Unknown";
}

void printRuntimePlanFailure(const Field& field,
                             const Storage& storage,
                             const IBMGeometry& geometry,
                             int ibmI, int ibmJ, int ibmK,
                             int requestedOrder,
                             PlanFailure failure) {
    static int emitted = 0;
    if (emitted >= 12) return;
    ++emitted;

    std::cerr << "[SF ILW] runtime plan build failed at ("
              << ibmI << "," << ibmJ << "," << ibmK
              << "): reason=" << planFailureName(failure)
              << ", flag=" << field.CellFlag(ibmI, ibmJ, ibmK)
              << ", requestedTaylorOrder=" << requestedOrder
              << ", fluidSamples="
              << storage.fluidCount((size_t)field.getIdx(ibmI, ibmJ, ibmK))
              << ", normalSamples="
              << storage.normalCount((size_t)field.getIdx(ibmI, ibmJ, ibmK))
              << ", hasGeometry="
              << (geometry.hasGeometry(ibmI, ibmJ, ibmK) ? 1 : 0)
              << ", hasPrecomputedPlan="
              << (storage.hasPlan((size_t)field.getIdx(ibmI, ibmJ, ibmK)) ? 1 : 0)
              << std::endl;
}

bool hasTangentialSpread(const std::vector<PlanSampleGeometry>& samples,
                         bool t2Direction) {
    if (samples.size() < 2) return false;
    double lo = std::numeric_limits<double>::max();
    double hi = -std::numeric_limits<double>::max();
    for (const auto& sample : samples) {
        const double v = t2Direction ? sample.b : sample.a;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    return (hi - lo) > 1e-10;
}

bool gatherTensorFluidSampleGeometry(const Field& field,
                                     const Storage& storage,
                                     int ibmI, int ibmJ, int ibmK,
                                     const double wallPoint[3],
                                     const LocalFrame& frame,
                                     std::vector<PlanSampleGeometry>& samples) {
    samples.clear();
    const size_t cell = (size_t)field.getIdx(ibmI, ibmJ, ibmK);
    const int nCached = storage.fluidCount(cell);
    if (nCached <= 0) return false;

    for (int n = 0; n < nCached; ++n) {
        int ii = 0, jj = 0, kk = 0;
        field.getIJK(storage.fluidSample(cell, n),
                     ii, jj, kk);
        if (!isIBMReadableCell(field, ii, jj, kk)) continue;
        if (field.CellFlag(ii, jj, kk) != FLUID_CELL) continue;

        double pt[3];
        cellPoint(field, ii, jj, kk, pt);
        double s = 0.0, a = 0.0, b = 0.0;
        localCoordinates(pt, wallPoint, frame, s, a, b);
        if (s >= -1e-12) continue;

        samples.push_back({field.getIdx(ii, jj, kk),
                           s, a, b, s * s + a * a + b * b});
    }

    std::sort(samples.begin(), samples.end(),
              [](const PlanSampleGeometry& lhs, const PlanSampleGeometry& rhs) {
                  return lhs.r2 < rhs.r2;
              });

    constexpr int maxSamples = 512;
    if ((int)samples.size() > maxSamples) samples.resize(maxSamples);
    return !samples.empty();
}

std::vector<Math::Polynomial::TensorPoint>
pointsFromGeometry(const std::vector<PlanSampleGeometry>& samples,
                   double scale) {
    const double invScale = 1.0 / std::max(1.0e-12, scale);
    std::vector<Math::Polynomial::TensorPoint> points;
    points.reserve(samples.size());
    for (const auto& sample : samples) {
        points.push_back({
            sample.s * invScale,
            sample.a * invScale,
            sample.b * invScale
        });
    }
    return points;
}

bool buildFitCandidate(const std::vector<PlanSampleGeometry>& samples,
                       bool useT2,
                       int polynomialOrder,
                       int derivativeOrder,
                       double linearWeight,
                       double localSpacing,
                       IBMILWFitCandidate& candidate) {
    Math::Polynomial::LeastSquaresPlan lsPlan;
    if (!Math::Polynomial::buildLeastSquaresPlan(
            pointsFromGeometry(samples, localSpacing),
            useT2, polynomialOrder, lsPlan)) {
        return false;
    }

    candidate.derivativeOrder = derivativeOrder;
    candidate.polynomialOrder = polynomialOrder;
    candidate.linearWeight = linearWeight;
    candidate.localSpacing = localSpacing;
    candidate.sampleCells.clear();
    candidate.sampleCells.reserve(samples.size());
    for (const auto& sample : samples) candidate.sampleCells.push_back(sample.idx);
    candidate.basisRows = lsPlan.basisRows;
    candidate.projectionRows = lsPlan.projectionRows;

    if (derivativeOrder > 0) {
        if (derivativeOrder > polynomialOrder) {
            candidate.derivativeWeights.assign(samples.size(), 0.0);
            return true;
        }
        const int coeffIndex = Math::Polynomial::normalDerivativeIndex(
            useT2, polynomialOrder, derivativeOrder);
        if (coeffIndex < 0 ||
            coeffIndex >= (int)candidate.projectionRows.size()) {
            return false;
        }
        candidate.derivativeWeights = candidate.projectionRows[(size_t)coeffIndex];
        const double invScaleDerivative =
            1.0 / std::pow(std::max(1.0e-12, localSpacing), derivativeOrder);
        for (double& weight : candidate.derivativeWeights) {
            weight *= invScaleDerivative;
        }
        double weightSum = 0.0;
        for (double weight : candidate.derivativeWeights) weightSum += weight;
        const double correction =
            weightSum / std::max<size_t>(1, candidate.derivativeWeights.size());
        for (double& weight : candidate.derivativeWeights) {
            weight -= correction;
        }
    }
    return true;
}

bool buildFitCandidateWithGrowingStencil(
        const std::vector<PlanSampleGeometry>& samples,
        bool useT2,
        int polynomialOrder,
        int derivativeOrder,
        int preferredSamples,
        double linearWeight,
        double localSpacing,
        bool minimizeAmplification,
        IBMILWFitCandidate& candidate) {
    const int nBasis = Math::Polynomial::tensorBasisSize(useT2, polynomialOrder);
    const int nSamples = (int)samples.size();
    if (nSamples < nBasis) return false;

    auto rowAmplification = [](std::vector<double> weights,
                               double requiredSum) {
        if (weights.empty()) return std::numeric_limits<double>::max();
        double sum = 0.0;
        for (double weight : weights) sum += weight;
        const double correction =
            (requiredSum - sum) / (double)weights.size();
        double amplification = 0.0;
        for (double weight : weights) {
            amplification += std::abs(weight + correction);
        }
        return amplification;
    };
    auto fitAmplification = [&](const IBMILWFitCandidate& fit) {
        if (fit.projectionRows.empty()) {
            return std::numeric_limits<double>::max();
        }
        double amplification =
            rowAmplification(fit.projectionRows[0], 1.0);
        const int d1Index = Math::Polynomial::normalDerivativeIndex(
            useT2, polynomialOrder, 1);
        if (d1Index >= 0 &&
            d1Index < (int)fit.projectionRows.size()) {
            amplification = std::max(
                amplification,
                rowAmplification(
                    fit.projectionRows[(size_t)d1Index], 0.0));
        }
        return amplification;
    };

    bool built = false;
    double bestAmplification = std::numeric_limits<double>::max();
    int nFit = std::min(nSamples, std::max(nBasis, preferredSamples));
    while (true) {
        std::vector<PlanSampleGeometry> substencil(samples.begin(),
                                                   samples.begin() + nFit);
        IBMILWFitCandidate trial;
        if (buildFitCandidate(substencil, useT2, polynomialOrder,
                              derivativeOrder, linearWeight, localSpacing,
                              trial)) {
            if (!minimizeAmplification) {
                candidate = std::move(trial);
                return true;
            }
            const double amplification = fitAmplification(trial);
            if (std::isfinite(amplification) &&
                amplification < bestAmplification) {
                bestAmplification = amplification;
                candidate = std::move(trial);
                built = true;
            }
        }
        if (nFit == nSamples) break;
        nFit = std::min(nSamples, std::max(nFit + nBasis, 2 * nFit));
    }

    return built;
}

bool appendFitCandidatesWithGrowingStencils(
        const std::vector<PlanSampleGeometry>& samples,
        bool useT2,
        int polynomialOrder,
        int derivativeOrder,
        int preferredSamples,
        double linearWeight,
        double localSpacing,
        std::vector<IBMILWFitCandidate>& candidates) {
    const int nBasis = Math::Polynomial::tensorBasisSize(useT2, polynomialOrder);
    const int nSamples = (int)samples.size();
    if (nSamples < nBasis) return false;

    std::vector<int> stencilSizes;
    int nFit = std::max(nBasis, preferredSamples);
    const int grow = std::max(4, nBasis);
    while (nFit < nSamples && (int)stencilSizes.size() < 6) {
        stencilSizes.push_back(nFit);
        nFit = std::min(nSamples, nFit + grow);
    }
    if (stencilSizes.empty() || stencilSizes.back() != nSamples) {
        stencilSizes.push_back(nSamples);
    }

    bool appended = false;
    bool hasDerivativeCapableCandidate = false;
    for (int size : stencilSizes) {
        std::vector<PlanSampleGeometry> substencil(samples.begin(),
                                                   samples.begin() + size);
        IBMILWFitCandidate candidate;
        if (!buildFitCandidate(substencil, useT2, polynomialOrder,
                               derivativeOrder, linearWeight, localSpacing,
                               candidate)) {
            continue;
        }
        candidates.push_back(candidate);
        appended = true;
        if (polynomialOrder >= derivativeOrder) {
            hasDerivativeCapableCandidate = true;
        }
    }
    return appended && hasDerivativeCapableCandidate;
}

/// @brief 用空间分箱法估计2D局部网格间距（中位数最近邻距离）。
///
/// 将样本按 (s, a) 坐标分箱，每个样本只在其所在箱及相邻箱内搜索最近邻，
/// 将 O(K²) 降为近似 O(K)。箱尺寸由包围盒和样本数自动估计。
/// @param samples 流体侧样本几何数组。
/// @return 样本间最近邻距离的中位数（≥ 1e-8）。
double estimateLocalSpacing2D(const std::vector<PlanSampleGeometry>& samples) {
    const int n = (int)samples.size();
    if (n < 2) return 1.0;

    // 计算 (s, a) 包围盒
    double sMin = std::numeric_limits<double>::max();
    double sMax = -std::numeric_limits<double>::max();
    double aMin = std::numeric_limits<double>::max();
    double aMax = -std::numeric_limits<double>::max();
    for (const auto& sample : samples) {
        sMin = std::min(sMin, sample.s);
        sMax = std::max(sMax, sample.s);
        aMin = std::min(aMin, sample.a);
        aMax = std::max(aMax, sample.a);
    }

    const double sRange = std::max(1.0e-12, sMax - sMin);
    const double aRange = std::max(1.0e-12, aMax - aMin);

    // 箱尺寸估计：假设样本近似均匀分布
    const double areaPerSample = sRange * aRange / (double)n;
    const double binSize = std::max(1.0e-12, std::sqrt(areaPerSample));

    const int nBinsS = std::max(1, (int)std::ceil(sRange / binSize));
    const int nBinsA = std::max(1, (int)std::ceil(aRange / binSize));
    const double invBinSize = 1.0 / binSize;

    // 将样本分配到箱
    // bins[idx] = 该箱内的样本索引列表
    std::vector<std::vector<int>> bins((size_t)(nBinsS * nBinsA));
    for (int i = 0; i < n; ++i) {
        const int bs = std::min(nBinsS - 1,
            std::max(0, (int)((samples[i].s - sMin) * invBinSize)));
        const int ba = std::min(nBinsA - 1,
            std::max(0, (int)((samples[i].a - aMin) * invBinSize)));
        bins[(size_t)(bs * nBinsA + ba)].push_back(i);
    }

    // 在每个样本的邻域箱内搜索最近邻
    std::vector<double> nearest;
    nearest.reserve((size_t)n);

    for (int i = 0; i < n; ++i) {
        const int bs = std::min(nBinsS - 1,
            std::max(0, (int)((samples[i].s - sMin) * invBinSize)));
        const int ba = std::min(nBinsA - 1,
            std::max(0, (int)((samples[i].a - aMin) * invBinSize)));

        double best = std::numeric_limits<double>::max();
        // 检查自身箱及相邻箱（3×3邻域）
        for (int dbs = -1; dbs <= 1; ++dbs) {
            const int nbs = bs + dbs;
            if (nbs < 0 || nbs >= nBinsS) continue;
            for (int dba = -1; dba <= 1; ++dba) {
                const int nba = ba + dba;
                if (nba < 0 || nba >= nBinsA) continue;
                for (int j : bins[(size_t)(nbs * nBinsA + nba)]) {
                    if (i == j) continue;
                    const double ds = samples[i].s - samples[j].s;
                    const double da = samples[i].a - samples[j].a;
                    const double d = ds * ds + da * da;
                    if (d > 1.0e-20) best = std::min(best, d);
                }
            }
        }

        if (std::isfinite(best) && best < std::numeric_limits<double>::max()) {
            nearest.push_back(std::sqrt(best));
        }
    }

    if (nearest.empty()) return 1.0;
    std::sort(nearest.begin(), nearest.end());
    return std::max(1.0e-8, nearest[nearest.size() / 2]);
}

std::vector<PlanSampleGeometry> uniqueSamples2D(
        const std::vector<PlanSampleGeometry>& samples,
        double h) {
    std::vector<PlanSampleGeometry> sorted = samples;
    std::sort(sorted.begin(), sorted.end(),
              [](const PlanSampleGeometry& lhs,
                 const PlanSampleGeometry& rhs) {
                  if (std::abs(lhs.a - rhs.a) > 1.0e-14) {
                      return lhs.a < rhs.a;
                  }
                  return lhs.s < rhs.s;
              });

    const double tol = std::max(1.0e-10, 1.0e-6 * h);
    std::vector<PlanSampleGeometry> unique;
    unique.reserve(sorted.size());
    for (const auto& sample : sorted) {
        if (!unique.empty() &&
            std::abs(sample.s - unique.back().s) <= tol &&
            std::abs(sample.a - unique.back().a) <= tol) {
            if (sample.r2 < unique.back().r2) unique.back() = sample;
            continue;
        }
        unique.push_back(sample);
    }
    std::sort(unique.begin(), unique.end(),
              [](const PlanSampleGeometry& lhs,
                 const PlanSampleGeometry& rhs) {
                  return lhs.r2 < rhs.r2;
              });
    return unique;
}

double estimateLocalSpacing3D(const std::vector<PlanSampleGeometry>& samples) {
    const int n = (int)samples.size();
    if (n < 2) return 1.0;

    std::vector<double> nearest((size_t)n, std::numeric_limits<double>::max());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double ss = samples[i].s - samples[j].s;
            const double aa = samples[i].a - samples[j].a;
            const double bb = samples[i].b - samples[j].b;
            const double d = ss * ss + aa * aa + bb * bb;
            if (d <= 1.0e-20) continue;
            nearest[(size_t)i] = std::min(nearest[(size_t)i], d);
            nearest[(size_t)j] = std::min(nearest[(size_t)j], d);
        }
    }

    std::vector<double> finiteNearest;
    finiteNearest.reserve((size_t)n);
    for (double d : nearest) {
        if (std::isfinite(d) && d < std::numeric_limits<double>::max()) {
            finiteNearest.push_back(std::sqrt(d));
        }
    }
    if (finiteNearest.empty()) return estimateLocalSpacing2D(samples);
    std::sort(finiteNearest.begin(), finiteNearest.end());
    return std::max(1.0e-8, finiteNearest[finiteNearest.size() / 2]);
}

std::vector<PlanSampleGeometry> uniqueSamples3D(
        const std::vector<PlanSampleGeometry>& samples,
        double h) {
    std::vector<PlanSampleGeometry> sorted = samples;
    std::sort(sorted.begin(), sorted.end(),
              [](const PlanSampleGeometry& lhs,
                 const PlanSampleGeometry& rhs) {
                  if (std::abs(lhs.b - rhs.b) > 1.0e-14) return lhs.b < rhs.b;
                  if (std::abs(lhs.a - rhs.a) > 1.0e-14) return lhs.a < rhs.a;
                  return lhs.s < rhs.s;
              });

    const double tol = std::max(1.0e-10, 1.0e-6 * h);
    std::vector<PlanSampleGeometry> unique;
    unique.reserve(sorted.size());
    for (const auto& sample : sorted) {
        if (!unique.empty() &&
            std::abs(sample.s - unique.back().s) <= tol &&
            std::abs(sample.a - unique.back().a) <= tol &&
            std::abs(sample.b - unique.back().b) <= tol) {
            if (sample.r2 < unique.back().r2) unique.back() = sample;
            continue;
        }
        unique.push_back(sample);
    }
    std::sort(unique.begin(), unique.end(),
              [](const PlanSampleGeometry& lhs,
                 const PlanSampleGeometry& rhs) {
                  return lhs.r2 < rhs.r2;
              });
    return unique;
}

struct Structured2DRow {
    double centerA = 0.0;
    std::vector<PlanSampleGeometry> columns;
};

struct Structured2DTemplate {
    double score = 0.0;
    std::vector<PlanSampleGeometry> samples;
};

struct Structured3DStack {
    int aId = 0;
    int bId = 0;
    double centerA = 0.0;
    double centerB = 0.0;
    std::vector<PlanSampleGeometry> columns;
};

struct Structured3DTemplate {
    double score = 0.0;
    std::vector<PlanSampleGeometry> samples;
};

std::vector<Structured2DRow> buildStructured2DRows(
        const std::vector<PlanSampleGeometry>& samples,
        double h) {
    std::map<int, std::vector<PlanSampleGeometry>> rowBins;
    const double rowSpacing = std::max(1.0e-8, h);
    for (const auto& sample : samples) {
        if (sample.s >= -1.0e-12) continue;
        const int rowId = (int)std::llround(sample.a / rowSpacing);
        rowBins[rowId].push_back(sample);
    }

    std::vector<Structured2DRow> rows;
    rows.reserve(rowBins.size());
    for (auto& rowEntry : rowBins) {
        auto& rowSamples = rowEntry.second;
        if (rowSamples.empty()) continue;

        double centerA = 0.0;
        for (const auto& sample : rowSamples) centerA += sample.a;
        centerA /= (double)rowSamples.size();

        struct ColumnPick {
            PlanSampleGeometry sample;
            double score = std::numeric_limits<double>::max();
        };
        std::map<int, ColumnPick> normalBins;
        for (const auto& sample : rowSamples) {
            const double depth = std::max(0.0, -sample.s);
            const int columnId = std::max(0, (int)std::floor(depth / rowSpacing));
            const double columnCenter = ((double)columnId + 0.5) * rowSpacing;
            const double score =
                std::abs(sample.a - centerA)
                + 0.25 * std::abs(depth - columnCenter);

            auto it = normalBins.find(columnId);
            if (it == normalBins.end() || score < it->second.score) {
                normalBins[columnId] = {sample, score};
            }
        }

        Structured2DRow row;
        row.centerA = centerA;
        row.columns.reserve(normalBins.size());
        for (const auto& columnEntry : normalBins) {
            row.columns.push_back(columnEntry.second.sample);
        }
        if (!row.columns.empty()) rows.push_back(row);
    }

    std::sort(rows.begin(), rows.end(),
              [](const Structured2DRow& lhs,
                 const Structured2DRow& rhs) {
                  return lhs.centerA < rhs.centerA;
              });
    return rows;
}

bool hasDuplicateSamples(const std::vector<PlanSampleGeometry>& samples) {
    std::vector<int> ids;
    ids.reserve(samples.size());
    for (const auto& sample : samples) ids.push_back(sample.idx);
    std::sort(ids.begin(), ids.end());
    return std::adjacent_find(ids.begin(), ids.end()) != ids.end();
}

std::vector<Structured2DTemplate> collectStructured2DTemplates(
        const std::vector<Structured2DRow>& rows,
        int polynomialOrder,
        double h) {
    const int n = polynomialOrder + 1;
    if ((int)rows.size() < n) return {};

    std::vector<Structured2DTemplate> templates;
    const double spacing = std::max(1.0e-8, h);
    for (int start = 0; start + n <= (int)rows.size(); ++start) {
        int minColumns = std::numeric_limits<int>::max();
        double meanA = 0.0;
        double spreadA = 0.0;
        for (int r = 0; r < n; ++r) {
            const auto& row = rows[(size_t)(start + r)];
            minColumns = std::min(minColumns, (int)row.columns.size());
            meanA += row.centerA;
        }
        if (minColumns < n) continue;
        meanA /= (double)n;
        for (int r = 0; r < n; ++r) {
            spreadA = std::max(spreadA,
                               std::abs(rows[(size_t)(start + r)].centerA - meanA));
        }

        const int maxOffset = std::min(2, minColumns - n);
        for (int offset = 0; offset <= maxOffset; ++offset) {
            Structured2DTemplate candidate;
            candidate.samples.reserve((size_t)n * (size_t)n);
            for (int r = 0; r < n; ++r) {
                const auto& row = rows[(size_t)(start + r)];
                for (int c = 0; c < n; ++c) {
                    candidate.samples.push_back(row.columns[(size_t)(offset + c)]);
                }
            }
            if (hasDuplicateSamples(candidate.samples)) continue;

            // Tan 2012 Sec. 2.4 的 E_r 是 (r+1)^2 个结构点。
            // 对IBM曲面，局部网格行可能有多个等价放置；这里按离壁点最近
            // 的切向窗口和法向起始列排序，并把同阶线性权重分摊到这些模板。
            candidate.score = std::abs(meanA) + 0.1 * spreadA
                            + 0.25 * (double)offset * spacing;
            templates.push_back(candidate);
        }
    }

    std::sort(templates.begin(), templates.end(),
              [](const Structured2DTemplate& lhs,
                 const Structured2DTemplate& rhs) {
                  return lhs.score < rhs.score;
              });
    if ((int)templates.size() > kMaxStructured2DVariantsPerOrder) {
        templates.resize(kMaxStructured2DVariantsPerOrder);
    }
    return templates;
}

std::vector<Structured3DStack> buildStructured3DStacks(
        const std::vector<PlanSampleGeometry>& samples,
        double h) {
    struct ColumnPick {
        PlanSampleGeometry sample;
        double score = std::numeric_limits<double>::max();
    };
    struct StackAccumulator {
        double sumA = 0.0;
        double sumB = 0.0;
        int count = 0;
        std::map<int, ColumnPick> normalBins;
    };

    const double spacing = std::max(1.0e-8, h);
    std::map<std::pair<int, int>, StackAccumulator> bins;
    for (const auto& sample : samples) {
        if (sample.s >= -1.0e-12) continue;
        const int aId = (int)std::llround(sample.a / spacing);
        const int bId = (int)std::llround(sample.b / spacing);
        const double depth = std::max(0.0, -sample.s);
        const int columnId = std::max(0, (int)std::floor(depth / spacing));
        const double columnCenter = ((double)columnId + 0.5) * spacing;
        const double score =
            std::abs(sample.a - (double)aId * spacing)
            + std::abs(sample.b - (double)bId * spacing)
            + 0.25 * std::abs(depth - columnCenter);

        auto& stack = bins[{aId, bId}];
        stack.sumA += sample.a;
        stack.sumB += sample.b;
        ++stack.count;
        auto it = stack.normalBins.find(columnId);
        if (it == stack.normalBins.end() || score < it->second.score) {
            stack.normalBins[columnId] = {sample, score};
        }
    }

    std::vector<Structured3DStack> stacks;
    stacks.reserve(bins.size());
    for (auto& entry : bins) {
        const auto& key = entry.first;
        auto& acc = entry.second;
        if (acc.count <= 0 || acc.normalBins.empty()) continue;

        Structured3DStack stack;
        stack.aId = key.first;
        stack.bId = key.second;
        stack.centerA = acc.sumA / (double)acc.count;
        stack.centerB = acc.sumB / (double)acc.count;
        stack.columns.reserve(acc.normalBins.size());
        for (const auto& column : acc.normalBins) {
            stack.columns.push_back(column.second.sample);
        }
        stacks.push_back(stack);
    }

    std::sort(stacks.begin(), stacks.end(),
              [](const Structured3DStack& lhs,
                 const Structured3DStack& rhs) {
                  if (lhs.aId != rhs.aId) return lhs.aId < rhs.aId;
                  return lhs.bId < rhs.bId;
              });
    return stacks;
}

std::vector<int> rankedStructuredWindowStarts(const std::vector<int>& ids,
                                              int n,
                                              double h) {
    std::vector<int> starts;
    if ((int)ids.size() < n) return starts;

    starts.reserve(ids.size());
    for (int start = 0; start + n <= (int)ids.size(); ++start) {
        starts.push_back(start);
    }

    const double spacing = std::max(1.0e-8, h);
    auto score = [&](int start) {
        double mean = 0.0;
        for (int q = 0; q < n; ++q) {
            mean += (double)ids[(size_t)(start + q)];
        }
        mean /= (double)n;

        double spread = 0.0;
        for (int q = 0; q < n; ++q) {
            spread = std::max(
                spread,
                std::abs((double)ids[(size_t)(start + q)] - mean));
        }
        return spacing * (std::abs(mean) + 0.05 * spread);
    };

    std::sort(starts.begin(), starts.end(),
              [&](int lhs, int rhs) {
                  return score(lhs) < score(rhs);
              });
    return starts;
}

std::vector<Structured3DTemplate> collectStructured3DTemplates(
        const std::vector<Structured3DStack>& stacks,
        int polynomialOrder,
        double h) {
    const int n = polynomialOrder + 1;
    if ((int)stacks.size() < n * n) return {};

    std::vector<int> aIds;
    std::vector<int> bIds;
    aIds.reserve(stacks.size());
    bIds.reserve(stacks.size());
    std::map<std::pair<int, int>, const Structured3DStack*> stackById;
    for (const auto& stack : stacks) {
        aIds.push_back(stack.aId);
        bIds.push_back(stack.bId);
        stackById[{stack.aId, stack.bId}] = &stack;
    }
    std::sort(aIds.begin(), aIds.end());
    std::sort(bIds.begin(), bIds.end());
    aIds.erase(std::unique(aIds.begin(), aIds.end()), aIds.end());
    bIds.erase(std::unique(bIds.begin(), bIds.end()), bIds.end());

    const double spacing = std::max(1.0e-8, h);
    const std::vector<int> aStarts =
        rankedStructuredWindowStarts(aIds, n, spacing);
    const std::vector<int> bStarts =
        rankedStructuredWindowStarts(bIds, n, spacing);
    const int localALimit = std::min((int)aStarts.size(),
                                     kMaxStructured3DWindowStartsPerAxis);
    const int localBLimit = std::min((int)bStarts.size(),
                                     kMaxStructured3DWindowStartsPerAxis);

    std::vector<Structured3DTemplate> templates;
    auto collectWindows = [&](int aLimit, int bLimit) {
        for (int aRank = 0; aRank < aLimit; ++aRank) {
            const int aStart = aStarts[(size_t)aRank];
            for (int bRank = 0; bRank < bLimit; ++bRank) {
                const int bStart = bStarts[(size_t)bRank];
                std::vector<const Structured3DStack*> window;
                window.reserve((size_t)n * (size_t)n);
                bool complete = true;
                for (int ia = 0; ia < n && complete; ++ia) {
                    for (int ib = 0; ib < n; ++ib) {
                        auto it = stackById.find({
                            aIds[(size_t)(aStart + ia)],
                            bIds[(size_t)(bStart + ib)]
                        });
                        if (it == stackById.end()) {
                            complete = false;
                            break;
                        }
                        window.push_back(it->second);
                    }
                }
                if (!complete) continue;

                int minColumns = std::numeric_limits<int>::max();
                double meanA = 0.0;
                double meanB = 0.0;
                double spreadA = 0.0;
                double spreadB = 0.0;
                for (const auto* stack : window) {
                    minColumns = std::min(minColumns, (int)stack->columns.size());
                    meanA += stack->centerA;
                    meanB += stack->centerB;
                }
                if (minColumns < n) continue;
                meanA /= (double)window.size();
                meanB /= (double)window.size();
                for (const auto* stack : window) {
                    spreadA = std::max(spreadA,
                                       std::abs(stack->centerA - meanA));
                    spreadB = std::max(spreadB,
                                       std::abs(stack->centerB - meanB));
                }

                const int maxOffset = std::min(2, minColumns - n);
                for (int offset = 0; offset <= maxOffset; ++offset) {
                    Structured3DTemplate candidate;
                    candidate.samples.reserve((size_t)n * (size_t)n * (size_t)n);
                    for (const auto* stack : window) {
                        for (int c = 0; c < n; ++c) {
                            candidate.samples.push_back(
                                stack->columns[(size_t)(offset + c)]);
                        }
                    }
                    if (hasDuplicateSamples(candidate.samples)) continue;

                    // 3D extension of Tan 2012 Sec. 2.4: E_r has (r+1)^3
                    // structured points in local (s,a,b), fitted with a total-degree
                    // polynomial and blended by WENO-type nonlinear weights.
                    candidate.score = std::sqrt(meanA * meanA + meanB * meanB)
                                    + 0.1 * (spreadA + spreadB)
                                    + 0.25 * (double)offset * spacing;
                    templates.push_back(candidate);
                    if ((int)templates.size() >=
                        kMaxStructured3DTemplatesExamined) {
                        return;
                    }
                }
            }
        }
    };

    collectWindows(localALimit, localBLimit);
    if (templates.empty() &&
        (localALimit < (int)aStarts.size() ||
         localBLimit < (int)bStarts.size())) {
        collectWindows((int)aStarts.size(), (int)bStarts.size());
    }

    std::sort(templates.begin(), templates.end(),
              [](const Structured3DTemplate& lhs,
                 const Structured3DTemplate& rhs) {
                  return lhs.score < rhs.score;
              });
    if ((int)templates.size() > kMaxStructured3DVariantsPerOrder) {
        templates.resize(kMaxStructured3DVariantsPerOrder);
    }
    return templates;
}

double tan2DLinearWeight(int polynomialOrder, int targetOrder, double h) {
    const double hw = std::max(1.0e-8, h);
    if (polynomialOrder < targetOrder) {
        return 2.0 * std::pow(hw, targetOrder - polynomialOrder);
    }

    double lowerSum = 0.0;
    for (int r = 0; r < targetOrder; ++r) {
        lowerSum += 2.0 * std::pow(hw, targetOrder - r);
    }
    return std::max(1.0e-12, 1.0 - lowerSum);
}

double evaluate2DPolynomialDerivative(
        const std::vector<Primitive>& coeff,
        int polynomialOrder,
        int ds,
        int da,
        double s,
        double a,
        int v) {
    double value = 0.0;
    int idx = 0;
    for (int total = 0; total <= polynomialOrder; ++total) {
        for (int ps = total; ps >= 0; --ps) {
            const int pa = total - ps;
            if (ps >= ds && pa >= da && idx < (int)coeff.size()) {
                const double denom =
                    Math::Taylor::factorial(ps - ds)
                    * Math::Taylor::factorial(pa - da);
                value += coeff[(size_t)idx][(size_t)v]
                       * Math::Polynomial::powInt(s, ps - ds)
                       * Math::Polynomial::powInt(a, pa - da)
                       / denom;
            }
            ++idx;
        }
    }
    return value;
}

double tan2DSmoothnessIndicator(const IBMILWFitCandidate& fit,
                                const std::vector<Primitive>& values,
                                double h) {
    if (values.empty() || fit.projectionRows.empty()) return 1.0e30;
    if (fit.polynomialOrder <= 0) {
        const double hw = std::max(1.0e-8, h);
        return 2.0 * hw * hw;
    }

    std::vector<Primitive> coeff;
    coeff.reserve(fit.projectionRows.size());
    for (const auto& row : fit.projectionRows) {
        coeff.push_back(Math::WeightedProjection::applyWeights<5>(row, values));
    }

    constexpr double gp[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526
    };
    constexpr double gw[4] = {
        0.3478548451374538, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374538
    };

    const double hw = std::max(1.0e-8, h);
    const double half = 0.5 * hw;
    const double jac = half * half;
    double beta = 0.0;

    for (int total = 1; total <= fit.polynomialOrder; ++total) {
        const double scale = std::pow(hw * hw, total - 1);
        for (int ds = total; ds >= 0; --ds) {
            const int da = total - ds;
            double integral = 0.0;
            for (int qs = 0; qs < 4; ++qs) {
                const double s = half * gp[qs];
                for (int qa = 0; qa < 4; ++qa) {
                    const double a = half * gp[qa];
                    double componentSum = 0.0;
                    for (int v = 0; v < 5; ++v) {
                        const double dv = evaluate2DPolynomialDerivative(
                            coeff, fit.polynomialOrder, ds, da, s, a, v);
                        componentSum += dv * dv;
                    }
                    integral += gw[qs] * gw[qa] * componentSum;
                }
            }
            beta += scale * jac * integral;
        }
    }

    return std::isfinite(beta) ? beta : 1.0e30;
}

double evaluate3DPolynomialDerivative(
        const std::vector<Primitive>& coeff,
        int polynomialOrder,
        int ds,
        int da,
        int db,
        double s,
        double a,
        double b,
        int v) {
    double value = 0.0;
    int idx = 0;
    for (int total = 0; total <= polynomialOrder; ++total) {
        for (int ps = total; ps >= 0; --ps) {
            for (int pa = total - ps; pa >= 0; --pa) {
                const int pb = total - ps - pa;
                if (ps >= ds && pa >= da && pb >= db
                    && idx < (int)coeff.size()) {
                    const double denom =
                        Math::Taylor::factorial(ps - ds)
                        * Math::Taylor::factorial(pa - da)
                        * Math::Taylor::factorial(pb - db);
                    value += coeff[(size_t)idx][(size_t)v]
                           * Math::Polynomial::powInt(s, ps - ds)
                           * Math::Polynomial::powInt(a, pa - da)
                           * Math::Polynomial::powInt(b, pb - db)
                           / denom;
                }
                ++idx;
            }
        }
    }
    return value;
}

double tan3DSmoothnessIndicator(const IBMILWFitCandidate& fit,
                                const std::vector<Primitive>& values,
                                double h) {
    if (values.empty() || fit.projectionRows.empty()) return 1.0e30;
    if (fit.polynomialOrder <= 0) {
        const double hw = std::max(1.0e-8, h);
        return 2.0 * hw * hw;
    }

    std::vector<Primitive> coeff;
    coeff.reserve(fit.projectionRows.size());
    for (const auto& row : fit.projectionRows) {
        coeff.push_back(Math::WeightedProjection::applyWeights<5>(row, values));
    }

    constexpr double gp[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526
    };
    constexpr double gw[4] = {
        0.3478548451374538, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374538
    };

    const double hw = std::max(1.0e-8, h);
    const double half = 0.5 * hw;
    const double jac = half * half * half;
    double beta = 0.0;

    for (int total = 1; total <= fit.polynomialOrder; ++total) {
        const double scale = std::pow(hw * hw, total - 1);
        for (int ds = total; ds >= 0; --ds) {
            const int tangentialOrder = total - ds;
            for (int da = tangentialOrder; da >= 0; --da) {
                const int db = tangentialOrder - da;
                double integral = 0.0;
                for (int qs = 0; qs < 4; ++qs) {
                    const double s = half * gp[qs];
                    for (int qa = 0; qa < 4; ++qa) {
                        const double a = half * gp[qa];
                        for (int qb = 0; qb < 4; ++qb) {
                            const double b = half * gp[qb];
                            double componentSum = 0.0;
                            for (int v = 0; v < 5; ++v) {
                                const double dv = evaluate3DPolynomialDerivative(
                                    coeff, fit.polynomialOrder,
                                    ds, da, db, s, a, b, v);
                                componentSum += dv * dv;
                            }
                            integral += gw[qs] * gw[qa] * gw[qb]
                                      * componentSum;
                        }
                    }
                }
                beta += scale * jac * integral;
            }
        }
    }

    return std::isfinite(beta) ? beta : 1.0e30;
}

bool appendStructured2DWENOFitCandidates(
        const std::vector<PlanSampleGeometry>& samples,
        int derivativeOrder,
        int targetOrder,
        std::vector<IBMILWFitCandidate>& candidates) {
    if (targetOrder < derivativeOrder) return false;

    const double h = estimateLocalSpacing2D(samples);
    const std::vector<PlanSampleGeometry> uniqueSamples =
        uniqueSamples2D(samples, h);
    const std::vector<Structured2DRow> rows =
        buildStructured2DRows(uniqueSamples, h);

    bool appended = false;
    bool hasDerivativeCapableCandidate = false;
    for (int polyOrder = 0; polyOrder <= targetOrder; ++polyOrder) {
        const std::vector<Structured2DTemplate> templates =
            collectStructured2DTemplates(rows, polyOrder, h);

        std::vector<IBMILWFitCandidate> built;
        built.reserve(templates.size());
        for (const auto& structuredTemplate : templates) {
            IBMILWFitCandidate candidate;
            if (!buildFitCandidate(structuredTemplate.samples, false,
                                   polyOrder, derivativeOrder,
                                   1.0, h, candidate)) {
                continue;
            }
            built.push_back(candidate);
        }

        if (built.empty()) {
            const int nBasis =
                Math::Polynomial::tensorBasisSize(false, polyOrder);
            const int structuredSize = (polyOrder + 1) * (polyOrder + 1);
            appendFitCandidatesWithGrowingStencils(
                uniqueSamples, false, polyOrder, derivativeOrder,
                std::max(nBasis, structuredSize),
                1.0, h, built);
        }
        if (built.empty()) continue;

        const double totalLinearWeight =
            tan2DLinearWeight(polyOrder, targetOrder, h);
        const double perTemplateWeight =
            totalLinearWeight / (double)built.size();
        for (auto& candidate : built) {
            candidate.linearWeight = perTemplateWeight;
            candidates.push_back(candidate);
        }
        appended = true;
        if (polyOrder >= derivativeOrder) {
            hasDerivativeCapableCandidate = true;
        }
    }

    return appended && hasDerivativeCapableCandidate;
}

bool appendStructured3DWENOFitCandidates(
        const std::vector<PlanSampleGeometry>& samples,
        int derivativeOrder,
        int targetOrder,
        std::vector<IBMILWFitCandidate>& candidates) {
    if (targetOrder < derivativeOrder) return false;

    const double h = estimateLocalSpacing3D(samples);
    const std::vector<PlanSampleGeometry> uniqueSamples =
        uniqueSamples3D(samples, h);
    const std::vector<Structured3DStack> stacks =
        buildStructured3DStacks(uniqueSamples, h);

    bool appended = false;
    bool hasDerivativeCapableCandidate = false;
    for (int polyOrder = 0; polyOrder <= targetOrder; ++polyOrder) {
        const std::vector<Structured3DTemplate> templates =
            collectStructured3DTemplates(stacks, polyOrder, h);
        if (templates.empty()) continue;

        std::vector<IBMILWFitCandidate> built;
        built.reserve(templates.size());
        for (const auto& structuredTemplate : templates) {
            IBMILWFitCandidate candidate;
            if (!buildFitCandidate(structuredTemplate.samples, true,
                                   polyOrder, derivativeOrder,
                                   1.0, h, candidate)) {
                continue;
            }
            built.push_back(candidate);
        }
        if (built.empty()) continue;

        const double totalLinearWeight =
            tan2DLinearWeight(polyOrder, targetOrder, h);
        const double perTemplateWeight =
            totalLinearWeight / (double)built.size();
        for (auto& candidate : built) {
            candidate.linearWeight = perTemplateWeight;
            candidates.push_back(candidate);
        }
        appended = true;
        if (polyOrder >= derivativeOrder) {
            hasDerivativeCapableCandidate = true;
        }
    }

    return appended && hasDerivativeCapableCandidate;
}

std::vector<PlanSampleGeometry> ranked3DTemplateSamples(
        const std::vector<PlanSampleGeometry>& samples,
        double h,
        int templateSize,
        double tangentBiasA,
        double tangentBiasB) {
    std::vector<PlanSampleGeometry> ranked = samples;
    if (std::abs(tangentBiasA) <= 1.0e-14 &&
        std::abs(tangentBiasB) <= 1.0e-14) {
        if ((int)ranked.size() > templateSize) {
            ranked.resize((size_t)templateSize);
        }
        return ranked;
    }

    const double invH = 1.0 / std::max(1.0e-8, h);
    auto score = [&](const PlanSampleGeometry& sample) {
        const double ss = sample.s * invH;
        const double aa = sample.a * invH;
        const double bb = sample.b * invH;
        const double centered = ss * ss + aa * aa + bb * bb;
        const double directional = tangentBiasA * aa + tangentBiasB * bb;
        return centered - 0.35 * directional + 0.05 * std::abs(ss);
    };

    std::sort(ranked.begin(), ranked.end(),
              [&](const PlanSampleGeometry& lhs,
                  const PlanSampleGeometry& rhs) {
                  return score(lhs) < score(rhs);
              });
    if ((int)ranked.size() > templateSize) {
        ranked.resize((size_t)templateSize);
    }
    return ranked;
}

bool hasTemplateSignature(
        const std::vector<PlanSampleGeometry>& samples,
        const std::vector<std::vector<int>>& usedSignatures) {
    std::vector<int> signature;
    signature.reserve(samples.size());
    for (const auto& sample : samples) signature.push_back(sample.idx);
    std::sort(signature.begin(), signature.end());
    for (const auto& used : usedSignatures) {
        if (used == signature) return true;
    }
    return false;
}

std::vector<int> templateSignature(
        const std::vector<PlanSampleGeometry>& samples) {
    std::vector<int> signature;
    signature.reserve(samples.size());
    for (const auto& sample : samples) signature.push_back(sample.idx);
    std::sort(signature.begin(), signature.end());
    return signature;
}

bool appendRanked3DWENOFitCandidates(
        const std::vector<PlanSampleGeometry>& samples,
        double h,
        int derivativeOrder,
        int targetOrder,
        std::vector<IBMILWFitCandidate>& candidates) {
    if (targetOrder < derivativeOrder) return false;

    const std::vector<PlanSampleGeometry>& uniqueSamples = samples;
    if (uniqueSamples.empty()) return false;

    constexpr double biasDirections[][2] = {
        { 0.0,  0.0},
        { 1.0,  0.0},
        {-1.0,  0.0},
        { 0.0,  1.0},
        { 0.0, -1.0}
    };

    bool appended = false;
    bool hasDerivativeCapableCandidate = false;
    for (int polyOrder = 0;
         polyOrder <= targetOrder;
         ++polyOrder) {
        const int nBasis =
            Math::Polynomial::tensorBasisSize(true, polyOrder);
        const int n = polyOrder + 1;
        const int tensorTemplateSize = n * n * n;
        const int templateSize = std::max(nBasis, tensorTemplateSize);
        if ((int)uniqueSamples.size() < templateSize) continue;

        std::vector<IBMILWFitCandidate> built;
        std::vector<std::vector<int>> usedSignatures;
        for (const auto& dir : biasDirections) {
            if ((int)built.size() >= kMaxRanked3DVariantsPerOrder) break;

            std::vector<PlanSampleGeometry> templateSamples =
                ranked3DTemplateSamples(uniqueSamples, h, templateSize,
                                        dir[0], dir[1]);
            if ((int)templateSamples.size() < templateSize ||
                hasDuplicateSamples(templateSamples) ||
                hasTemplateSignature(templateSamples, usedSignatures)) {
                continue;
            }

            IBMILWFitCandidate candidate;
            if (!buildFitCandidate(templateSamples, true,
                                   polyOrder, derivativeOrder,
                                   1.0, h, candidate)) {
                continue;
            }
            usedSignatures.push_back(templateSignature(templateSamples));
            built.push_back(candidate);
        }
        if (built.empty()) continue;

        const double totalLinearWeight =
            tan2DLinearWeight(polyOrder, targetOrder, h);
        const double perTemplateWeight =
            totalLinearWeight / (double)built.size();
        for (auto& candidate : built) {
            candidate.linearWeight = perTemplateWeight;
            candidates.push_back(candidate);
        }
        appended = true;
        if (polyOrder >= derivativeOrder) {
            hasDerivativeCapableCandidate = true;
        }
    }

    return appended && hasDerivativeCapableCandidate;
}

bool buildPrecomputedPlanForCell(const Field& field,
                                 const Storage& storage,
                                 const IBMGeometry& geometry,
                                 int ibmI, int ibmJ, int ibmK,
                                 int requestedMaxOrder,
                                 const IBMRuntimeConfig& config,
                                 IBMILWPointPlan& plan,
                                 PlanFailure& failure,
                                 bool minimizeWallFitAmplification,
                                 PreprocessStats* timingStats = nullptr) {
    failure = PlanFailure::None;
    plan = IBMILWPointPlan{};

    double wallPoint[3] = {0.0, 0.0, 0.0};
    double exteriorNormal[3] = {1.0, 0.0, 0.0};
    double wallVelocity[3] = {0.0, 0.0, 0.0};
    double targetDistance = 0.0;
    auto phaseStart = Clock::now();
    if (!buildFluidVisibleBoundaryGeometry(field, geometry, ibmI, ibmJ, ibmK,
                                           wallPoint, exteriorNormal,
                                           wallVelocity, targetDistance)) {
        if (timingStats) {
            timingStats->planGeometrySeconds += secondsSince(phaseStart);
        }
        failure = PlanFailure::Fit;
        return false;
    }
    if (timingStats) {
        timingStats->planGeometrySeconds += secondsSince(phaseStart);
    }

    const int inactiveAxis = Math::singleInactiveDirectionIndex();
    const bool useTwoDimensionalILW =
        (Math::activeDimensionCount() == 2 && inactiveAxis >= 0);
    const std::array<double, 3> inactiveNormal =
        Math::inactiveDirectionNormalVector();
    LocalFrame frame = useTwoDimensionalILW
        ? makeLocalFrame2D(exteriorNormal, inactiveNormal)
        : makeLocalFrame(exteriorNormal);
    std::vector<PlanSampleGeometry> samples;
    phaseStart = Clock::now();
    if (!gatherTensorFluidSampleGeometry(field, storage, ibmI, ibmJ, ibmK,
                                         wallPoint, frame, samples)) {
        if (timingStats) {
            timingStats->planSampleSeconds += secondsSince(phaseStart);
        }
        failure = PlanFailure::FluidSamples;
        return false;
    }
    if (timingStats) {
        timingStats->planSampleSeconds += secondsSince(phaseStart);
    }

    phaseStart = Clock::now();
    const double localSpacing = useTwoDimensionalILW
        ? estimateLocalSpacing2D(samples)
        : estimateLocalSpacing3D(samples);
    const std::vector<PlanSampleGeometry> uniqueSamples2DForFit =
        useTwoDimensionalILW ? uniqueSamples2D(samples, localSpacing)
                             : std::vector<PlanSampleGeometry>{};
    const std::vector<PlanSampleGeometry> uniqueSamples3DForFit =
        useTwoDimensionalILW ? std::vector<PlanSampleGeometry>{}
                             : uniqueSamples3D(samples, localSpacing);
    const std::vector<PlanSampleGeometry>& fitSamples =
        useTwoDimensionalILW ? uniqueSamples2DForFit : uniqueSamples3DForFit;
    if (timingStats) {
        timingStats->planSpacingSeconds += secondsSince(phaseStart);
    }

    const int maxOrder =
        std::max(1, std::min(requestedMaxOrder, kMaxSupportedTaylorOrder));
    bool useT2 = false;
    phaseStart = Clock::now();
    if (!useTwoDimensionalILW &&
        hasTangentialSpread(samples, false) &&
        hasTangentialSpread(samples, true)) {
        const int probeOrder = std::min(maxOrder, 2);
        const int nBasis = Math::Polynomial::tensorBasisSize(true, probeOrder);
        if ((int)fitSamples.size() >= nBasis) {
            IBMILWFitCandidate probeFit;
            useT2 = buildFitCandidateWithGrowingStencil(
                fitSamples, true, probeOrder, 0, nBasis + 8,
                1.0, localSpacing, minimizeWallFitAmplification, probeFit);
        }
    }
    if (timingStats) {
        timingStats->planProbeSeconds += secondsSince(phaseStart);
    }

    IBMILWFitCandidate wallFit;
    phaseStart = Clock::now();
    if (!buildFitCandidateWithGrowingStencil(fitSamples, useT2, 1, 0,
                                             Math::Polynomial::tensorBasisSize(useT2, 1) + 8,
                                             1.0, localSpacing,
                                             minimizeWallFitAmplification,
                                             wallFit)) {
        if (timingStats) {
            timingStats->planWallFitSeconds += secondsSince(phaseStart);
        }
        failure = PlanFailure::Fit;
        return false;
    }
    if (timingStats) {
        timingStats->planWallFitSeconds += secondsSince(phaseStart);
    }

    CurvatureTensor curvature;
    phaseStart = Clock::now();
    if (!estimateCurvatureTensor(field, storage, geometry, ibmI, ibmJ, ibmK,
                                 wallPoint, frame, config, useT2,
                                 curvature)) {
        if (timingStats) {
            timingStats->planCurvatureSeconds += secondsSince(phaseStart);
        }
        failure = PlanFailure::Curvature;
        return false;
    }
    if (timingStats) {
        timingStats->planCurvatureSeconds += secondsSince(phaseStart);
    }

    plan.valid = true;
    plan.useT2 = useT2;
    plan.maxOrder = 1;
    plan.wallPoint = {wallPoint[0], wallPoint[1], wallPoint[2]};
    plan.normal = {frame.n[0], frame.n[1], frame.n[2]};
    plan.tangent1 = {frame.t1[0], frame.t1[1], frame.t1[2]};
    plan.tangent2 = {frame.t2[0], frame.t2[1], frame.t2[2]};
    plan.wallVelocity = {wallVelocity[0], wallVelocity[1], wallVelocity[2]};
    plan.targetDistance = targetDistance;
    plan.curvatureK11 = curvature.k11;
    plan.curvatureK12 = curvature.k12;
    plan.curvatureK22 = curvature.k22;
    plan.wallFit = wallFit;

    phaseStart = Clock::now();
    for (int derivativeOrder = 2; derivativeOrder <= maxOrder; ++derivativeOrder) {
        bool hasCandidate = false;
        if (useTwoDimensionalILW) {
            hasCandidate = appendStructured2DWENOFitCandidates(
                fitSamples, derivativeOrder, maxOrder, plan.higherFits);
        } else if (useT2) {
            hasCandidate = appendRanked3DWENOFitCandidates(
                fitSamples, localSpacing, derivativeOrder, maxOrder,
                plan.higherFits);
        } else {
            hasCandidate = appendStructured2DWENOFitCandidates(
                fitSamples, derivativeOrder, maxOrder, plan.higherFits);
        }

        if (!hasCandidate) break;
        plan.maxOrder = derivativeOrder;
    }
    if (timingStats) {
        timingStats->planHigherFitSeconds += secondsSince(phaseStart);
    }

    return true;
}

LocalFrame frameFromPlan(const IBMILWPointPlan& plan) {
    LocalFrame frame;
    for (int c = 0; c < 3; ++c) {
        frame.n[c] = plan.normal[(size_t)c];
        frame.t1[c] = plan.tangent1[(size_t)c];
        frame.t2[c] = plan.tangent2[(size_t)c];
    }
    return frame;
}

bool loadFitValues(const Field& field,
                   const IBMILWFitCandidate& fit,
                   const LocalFrame& frame,
                   double gamma,
                   std::vector<Primitive>& values) {
    values.clear();
    values.reserve(fit.sampleCells.size());
    for (int idx : fit.sampleCells) {
        int i = 0, j = 0, k = 0;
        field.getIJK(idx, i, j, k);
        if (!isIBMReadableCell(field, i, j, k)) return false;
        if (field.CellFlag(i, j, k) != FLUID_CELL) return false;

        double q[5];
        loadConservative(field, i, j, k, q);
        Numerics::requirePhysicalState("ILW primitive", q[0], q[1],
                                       q[2], q[3], q[4], gamma);
        values.push_back(localPrimitiveFromConservative(q, frame, gamma));
    }
    return !values.empty();
}

/// @brief 加载拟合样本并转换为 characteristic-variable 扰动。
bool loadCharacteristicFitValues(const Field& field,
                                 const IBMILWFitCandidate& fit,
                                 const LocalFrame& frame,
                                 const Primitive& reference,
                                 double gamma,
                                 std::vector<Primitive>& values) {
    std::vector<Primitive> primitiveValues;
    if (!loadFitValues(field, fit, frame, gamma, primitiveValues)) {
        return false;
    }

    values.clear();
    values.reserve(primitiveValues.size());
    for (const Primitive& primitive : primitiveValues) {
        Primitive characteristic{};
        if (!primitiveToCharacteristicDifference(primitive, reference,
                                                 gamma,
                                                 characteristic)) {
            return false;
        }
        values.push_back(characteristic);
    }
    return !values.empty();
}

bool solveSlipWallFirstDerivative(const IBMILWPointPlan& plan,
                                  const Primitive& wallState,
                                  const Primitive& wallVelocity,
                                  const Primitive& outgoingCharacteristicD1,
                                  double gamma,
                                  Primitive& derivative) {
    const double c = referenceSoundSpeed(wallState, gamma);
    if (!(c > 0.0) || !std::isfinite(c)) return false;

    const double rho = wallState[LRHO];
    const double relUt1 = wallState[LUT1] - wallVelocity[LUT1];
    const double relUt2 = wallState[LUT2] - wallVelocity[LUT2];
    const double pressureNormal = -rho * (
        plan.curvatureK11 * relUt1 * relUt1
        + 2.0 * plan.curvatureK12 * relUt1 * relUt2
        + plan.curvatureK22 * relUt2 * relUt2);

    derivative[LP] = pressureNormal;
    derivative[LRHO] =
        (pressureNormal - outgoingCharacteristicD1[LRHO]) / (c * c);
    derivative[LUN] =
        (outgoingCharacteristicD1[LUN] - pressureNormal) / (rho * c);
    derivative[LUT1] = outgoingCharacteristicD1[LUT1];
    derivative[LUT2] = outgoingCharacteristicD1[LUT2];

    for (double v : derivative) {
        if (!std::isfinite(v)) return false;
    }
    return true;
}

bool buildWallStateDerivativeFromPlan(const Field& field,
                                      const IBMILWPointPlan& plan,
                                      const LocalFrame& frame,
                                      double gamma,
                                      Primitive& derivative) {
    std::vector<Primitive> wallValues;
    if (!loadFitValues(field, plan.wallFit, frame, gamma, wallValues)) {
        return false;
    }
    if (plan.wallFit.projectionRows.empty()) return false;

    std::vector<double> wallValueWeights = plan.wallFit.projectionRows[0];
    double wallValueWeightSum = 0.0;
    for (double weight : wallValueWeights) wallValueWeightSum += weight;
    const double wallValueCorrection =
        (1.0 - wallValueWeightSum)
        / std::max<size_t>(1, wallValueWeights.size());
    for (double& weight : wallValueWeights) weight += wallValueCorrection;

    derivative = Math::WeightedProjection::applyWeights<5>(
        wallValueWeights, wallValues);
    const double wallVelocity[3] = {
        plan.wallVelocity[0], plan.wallVelocity[1], plan.wallVelocity[2]
    };
    const Primitive wallVelLocal = localWallVelocity(wallVelocity, frame);
    derivative[LUN] = wallVelLocal[LUN];

    if (!std::isfinite(derivative[LRHO]) || derivative[LRHO] <= 0.0 ||
        !std::isfinite(derivative[LP]) || derivative[LP] <= 0.0) {
        double weightL1 = 0.0;
        double weightMin = std::numeric_limits<double>::max();
        double weightMax = -std::numeric_limits<double>::max();
        for (double weight : wallValueWeights) {
            weightL1 += std::abs(weight);
            weightMin = std::min(weightMin, weight);
            weightMax = std::max(weightMax, weight);
        }
        double rhoMin = std::numeric_limits<double>::max();
        double rhoMax = -std::numeric_limits<double>::max();
        double pressureMin = std::numeric_limits<double>::max();
        double pressureMax = -std::numeric_limits<double>::max();
        int communicationHaloSamples = 0;
        for (size_t n = 0; n < wallValues.size(); ++n) {
            rhoMin = std::min(rhoMin, wallValues[n][LRHO]);
            rhoMax = std::max(rhoMax, wallValues[n][LRHO]);
            pressureMin = std::min(pressureMin, wallValues[n][LP]);
            pressureMax = std::max(pressureMax, wallValues[n][LP]);
            int i = 0, j = 0, k = 0;
            field.getIJK(plan.wallFit.sampleCells[n], i, j, k);
            if (field.isCommunicationHalo(i, j, k)) {
                ++communicationHaloSamples;
            }
        }
        std::cerr
            << "[SF ILW] invalid wall-state diagnostic"
            << ": wallPoint=(" << plan.wallPoint[0] << ","
            << plan.wallPoint[1] << "," << plan.wallPoint[2] << ")"
            << ", samples=" << wallValues.size()
            << ", communicationHaloSamples=" << communicationHaloSamples
            << ", weightL1=" << weightL1
            << ", weightMin=" << weightMin
            << ", weightMax=" << weightMax
            << ", sampleRhoRange=[" << rhoMin << "," << rhoMax << "]"
            << ", samplePressureRange=[" << pressureMin << ","
            << pressureMax << "]"
            << ", extrapolatedRho=" << derivative[LRHO]
            << ", extrapolatedPressure=" << derivative[LP]
            << std::endl;
        return ilwFail("precomputed ILW extrapolated invalid wall primitive state");
    }
    return true;
}

bool buildFirstDerivativeFromPlan(const Field& field,
                                  const IBMILWPointPlan& plan,
                                  const LocalFrame& frame,
                                  const Primitive& wallState,
                                  double gamma,
                                  Primitive& derivative) {
    const int d1Index = Math::Polynomial::normalDerivativeIndex(
        plan.useT2, plan.wallFit.polynomialOrder, 1);
    if (d1Index < 0 || d1Index >= (int)plan.wallFit.projectionRows.size()) {
        return false;
    }

    std::vector<Primitive> wallCharacteristicValues;
    if (!loadCharacteristicFitValues(field, plan.wallFit, frame,
                                     wallState, gamma,
                                     wallCharacteristicValues)) {
        return false;
    }

    std::vector<double> d1Weights =
        plan.wallFit.projectionRows[(size_t)d1Index];
    const double invWallFitScale =
        1.0 / std::max(1.0e-12, plan.wallFit.localSpacing);
    for (double& weight : d1Weights) weight *= invWallFitScale;
    double d1WeightSum = 0.0;
    for (double weight : d1Weights) d1WeightSum += weight;
    const double d1Correction =
        d1WeightSum / std::max<size_t>(1, d1Weights.size());
    for (double& weight : d1Weights) weight -= d1Correction;

    Primitive physicalOutgoingCharacteristicD1 =
        Math::WeightedProjection::applyWeights<5>(
            d1Weights, wallCharacteristicValues);
    const double wallVelocity[3] = {
        plan.wallVelocity[0], plan.wallVelocity[1], plan.wallVelocity[2]
    };
    const Primitive wallVelLocal = localWallVelocity(wallVelocity, frame);
    if (!solveSlipWallFirstDerivative(plan, wallState, wallVelLocal,
                                      physicalOutgoingCharacteristicD1,
                                      gamma,
                                      derivative)) {
        return ilwFail("precomputed ILW first-derivative system failed");
    }
    return true;
}

bool buildHigherDerivativeFromFits(const Field& field,
                                   const IBMILWPointPlan& plan,
                                   const LocalFrame& frame,
                                   const Primitive& wallState,
                                   int targetOrder,
                                   int derivativeOrder,
                                   double gamma,
                                   Primitive& derivative) {
    Primitive weightedCharacteristic{};
    double alphaSum = 0.0;

    for (const auto& fit : plan.higherFits) {
        if (fit.derivativeOrder != derivativeOrder) continue;
        if (fit.polynomialOrder > targetOrder) continue;
        if (fit.derivativeWeights.empty()) continue;

        std::vector<Primitive> characteristicValues;
        if (!loadCharacteristicFitValues(field, fit, frame,
                                         wallState, gamma,
                                         characteristicValues)) {
            continue;
        }
        double beta = 0.0;
        if (!plan.useT2) {
            beta = tan2DSmoothnessIndicator(fit, characteristicValues,
                                            fit.localSpacing);
        } else {
            beta = tan3DSmoothnessIndicator(fit, characteristicValues,
                                            fit.localSpacing);
        }
        const double linearWeight = (fit.linearWeight > 0.0)
            ? fit.linearWeight
            : ((fit.polynomialOrder == targetOrder) ? 1.0 : 0.25);
        const double alpha = linearWeight
            / std::pow(kWenoExtrapolationEpsilon + beta,
                       kWenoExtrapolationPower);
        if (!std::isfinite(alpha) || alpha <= 0.0) continue;

        Primitive candidate = Math::WeightedProjection::applyWeights<5>(
            fit.derivativeWeights, characteristicValues);
        for (int v = 0; v < 5; ++v) {
            weightedCharacteristic[(size_t)v] +=
                alpha * candidate[(size_t)v];
        }
        alphaSum += alpha;
    }

    if (alphaSum <= 0.0) return false;

    Primitive characteristicDerivative{};
    for (int v = 0; v < 5; ++v) {
        characteristicDerivative[(size_t)v] =
            weightedCharacteristic[(size_t)v] / alphaSum;
    }

    return characteristicDerivativeToPrimitive(characteristicDerivative,
                                               wallState, gamma, derivative);
}

bool buildTensorILWDerivativesFromPlan(const Field& field,
                                       const IBMILWPointPlan& plan,
                                       int requestedOrder,
                                       double gamma,
                                       std::vector<Primitive>& derivatives) {
    derivatives.clear();
    if (!plan.valid) return false;

    const int targetOrder =
        std::max(1, std::min(requestedOrder, kMaxSupportedTaylorOrder));
    if (plan.maxOrder < targetOrder) return false;

    const LocalFrame frame = frameFromPlan(plan);
    auto derivativeBuilder =
        [&](int derivativeOrder,
            const std::vector<Primitive>& built,
            Primitive& derivative) -> bool {
            if (derivativeOrder == 0) {
                return buildWallStateDerivativeFromPlan(
                    field, plan, frame, gamma, derivative);
            }

            if (built.empty()) return false;
            const Primitive& wallState = built[0];
            if (derivativeOrder == 1) {
                return buildFirstDerivativeFromPlan(
                    field, plan, frame, wallState, gamma, derivative);
            }

            return buildHigherDerivativeFromFits(
                field, plan, frame, wallState,
                targetOrder, derivativeOrder, gamma, derivative);
        };

    return Math::Taylor::buildDerivativeTable<5>(
        targetOrder, derivativeBuilder, derivatives);
}

// ═══════════════════════════════════════════════════════════════
//  Taylor 外推 & Ghost 状态构造
// ═══════════════════════════════════════════════════════════════

bool buildTaylorGhostState(const std::vector<Primitive>& localDerivatives,
                           const LocalFrame& frame,
                           int requestedOrder,
                           double ghostDistance,
                           double gamma,
                           double qGhost[5]) {
    const int order = std::min(requestedOrder,
                               (int)localDerivatives.size() - 1);
    if (order < 0) return false;

    Primitive localGhost = Math::Taylor::evaluate1D<5>(
        localDerivatives, ghostDistance, order);
    if (!std::isfinite(localGhost[LRHO]) || localGhost[LRHO] <= 0.0 ||
        !std::isfinite(localGhost[LP]) || localGhost[LP] <= 0.0) {
        return false;
    }
    Primitive ghostGlobal = toGlobalPrimitive(localGhost, frame);
    conservativeFromPrimitive(ghostGlobal, gamma, qGhost);

    if (!std::isfinite(qGhost[0]) || !std::isfinite(qGhost[1]) ||
        !std::isfinite(qGhost[2]) || !std::isfinite(qGhost[3]) ||
        !std::isfinite(qGhost[4]) || qGhost[0] <= 0.0) {
        return false;
    }
    const double pCheck = Numerics::pressure(
        qGhost[0], qGhost[1], qGhost[2], qGhost[3], qGhost[4], gamma);
    if (!std::isfinite(pCheck) || pCheck <= 0.0) {
        return false;
    }

    return true;
}

void printPrimitiveDebug(const char* label, const Primitive& q) {
    std::cerr << label
              << "[rho=" << q[LRHO]
              << ", un=" << q[LUN]
              << ", ut1=" << q[LUT1]
              << ", ut2=" << q[LUT2]
              << ", p=" << q[LP] << "]";
}

void printTaylorFailure(int ibmI, int ibmJ, int ibmK,
                        const IBMILWPointPlan& plan,
                        int requestedOrder,
                        const std::vector<Primitive>& derivatives) {
    static int emitted = 0;
    if (emitted >= 8) return;
    ++emitted;

    const int order = std::min(requestedOrder,
                               (int)derivatives.size() - 1);
    Primitive localGhost = Math::Taylor::evaluate1D<5>(
        derivatives, plan.targetDistance, order);

    std::cerr << "[SF ILW] Taylor ghost state invalid at ("
              << ibmI << "," << ibmJ << "," << ibmK
              << "), requestedOrder=" << requestedOrder
              << ", usedOrder=" << order
              << ", distance=" << plan.targetDistance
              << ", maxOrder=" << plan.maxOrder
              << ", useT2=" << plan.useT2
              << ", fits=" << plan.higherFits.size()
              << std::endl;
    printPrimitiveDebug("  ghost=", localGhost);
    std::cerr << std::endl;
    for (size_t n = 0; n < derivatives.size(); ++n) {
        std::cerr << "  d" << n << "=";
        printPrimitiveDebug("", derivatives[n]);
        std::cerr << std::endl;
    }
    int fitPrinted = 0;
    for (const auto& fit : plan.higherFits) {
        if (fit.derivativeOrder != order) continue;
        if (fitPrinted >= 12) break;
        std::cerr << "  fit[" << fitPrinted
                  << "]: derivativeOrder=" << fit.derivativeOrder
                  << ", polynomialOrder=" << fit.polynomialOrder
                  << ", samples=" << fit.sampleCells.size()
                  << ", linearWeight=" << fit.linearWeight
                  << ", localSpacing=" << fit.localSpacing
                  << std::endl;
        ++fitPrinted;
    }
}

bool buildFluidVisibleBoundaryGeometry(const Field& field,
                                       const IBMGeometry& geometry,
                                       int ibmI, int ibmJ, int ibmK,
                                       double wallPoint[3],
                                       double exteriorNormal[3],
                                       double wallVelocity[3],
                                       double& targetDistance) {
    if (!geometry.hasGeometry(ibmI, ibmJ, ibmK)) {
        return ilwFail("tensor ILW missing cached IBM wall geometry");
    }

    const Vector3 w = geometry.wallPoint(ibmI, ibmJ, ibmK);
    const Vector3 n = geometry.wallNormal(ibmI, ibmJ, ibmK);
    const Vector3 vel = geometry.wallVelocity(ibmI, ibmJ, ibmK);
    wallPoint[0] = w.x;
    wallPoint[1] = w.y;
    wallPoint[2] = w.z;

    // Field中的IBM法向按普通ghost-cell约定从固体指向流体；
    // Tan局部坐标的外法向从流体域指向固体侧，因此这里取反。
    exteriorNormal[0] = -n.x;
    exteriorNormal[1] = -n.y;
    exteriorNormal[2] = -n.z;
    Math::normalizeArray(exteriorNormal);

    double ghostPoint[3];
    cellPoint(field, ibmI, ibmJ, ibmK, ghostPoint);
    targetDistance = signedDistanceFromWall(ghostPoint, wallPoint, exteriorNormal);
    if (!std::isfinite(targetDistance) || targetDistance <= 1e-12) {
        return ilwFail("tensor ILW cached boundary distance is invalid");
    }

    wallVelocity[0] = vel.x;
    wallVelocity[1] = vel.y;
    wallVelocity[2] = vel.z;

    return true;
}

/// @brief 高阶 ILW ghost 状态构造。
bool highOrderGhostStateForCell(const Field& field,
                                const Storage& storage,
                                const IBMGeometry& geometry,
                                int ibmI, int ibmJ, int ibmK,
                                Math::Dir d,
                                int requestedOrder,
                                const IBMRuntimeConfig& config,
                                double qGhost[5]) {
    (void)d;

    IBMILWPointPlan localPlan;
    const IBMILWPointPlan* planPtr = nullptr;
    const size_t cell = (size_t)field.getIdx(ibmI, ibmJ, ibmK);
    if (storage.hasPlan(cell)) {
        planPtr = &storage.plan(cell);
    } else {
        PlanFailure failure = PlanFailure::None;
        if (!buildPrecomputedPlanForCell(field, storage, geometry,
                                         ibmI, ibmJ, ibmK,
                                         requestedOrder,
                                         config,
                                         localPlan, failure, false)) {
            printRuntimePlanFailure(field, storage, geometry,
                                    ibmI, ibmJ, ibmK,
                                    requestedOrder, failure);
            return false;
        }
        planPtr = &localPlan;
    }
    const IBMILWPointPlan& plan = *planPtr;
    const LocalFrame frame = frameFromPlan(plan);

    std::vector<Primitive> derivatives;
    if (!buildTensorILWDerivativesFromPlan(field, plan,
                                           requestedOrder, config.gamma,
                                           derivatives)) {
        ilwFail("precomputed ILW derivative construction failed");
        return false;
    }

    if (!buildTaylorGhostState(derivatives, frame, requestedOrder,
                               plan.targetDistance, config.gamma, qGhost)) {
        printTaylorFailure(ibmI, ibmJ, ibmK, plan,
                           requestedOrder, derivatives);
        return false;
    }
    return true;
}

} // namespace

// ═══════════════════════════════════════════════════════════════
//  公开 API
// ═══════════════════════════════════════════════════════════════

PreprocessStats preprocessIBMGeometry(const Field& field,
                                      Storage& storage,
                                      const IBMGeometry& geometry,
                                      const IBMRuntimeConfig& config,
                                      int requestedMaxOrder,
                                      bool minimizeWallFitAmplification) {
    PreprocessStats stats;
    storage.clearPlans();

    const int maxOrder =
        std::max(1, std::min(requestedMaxOrder, kMaxSupportedTaylorOrder));
    for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
        for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
            for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                if (field.CellFlag(i, j, k) != IBM_GHOST_CELL) continue;
                if (!geometry.hasGeometry(i, j, k)) continue;

                ++stats.candidates;
                ++stats.ghostCandidates;
                IBMILWPointPlan plan;
                PlanFailure failure = PlanFailure::None;

                if (!buildPrecomputedPlanForCell(field, storage, geometry,
                                                 i, j, k,
                                                 maxOrder, config,
                                                 plan, failure,
                                                 minimizeWallFitAmplification,
                                                 &stats)) {
                    if (failure == PlanFailure::FluidSamples) {
                        ++stats.insufficientFluidSamples;
                    } else if (failure == PlanFailure::Curvature) {
                        ++stats.curvatureFailures;
                    } else {
                        ++stats.fitFailures;
                    }
                    continue;
                }

                if (plan.maxOrder < 2) ++stats.maxOrderLessThanTwo;
                if (plan.maxOrder < 4) ++stats.maxOrderLessThanFour;
                if (plan.maxOrder < 6) ++stats.maxOrderLessThanSix;
                if (plan.maxOrder < 8) ++stats.maxOrderLessThanEight;
                storage.setPlan((size_t)field.getIdx(i, j, k), plan);
                ++stats.directBuiltPlans;
                ++stats.built;
            }
        }
    }

    return stats;
}

bool buildGhostStateForCell(const Field& field,
                            const Storage& storage,
                            const IBMGeometry& geometry,
                            int ibmI, int ibmJ, int ibmK,
                            Math::Dir d,
                            int requestedOrder,
                            const IBMRuntimeConfig& config,
                            double qGhost[5]) {
    return highOrderGhostStateForCell(field, storage, geometry,
                                      ibmI, ibmJ, ibmK,
                                      d, requestedOrder, config, qGhost);
}

} // namespace GhostILW
} // namespace IBM
} // namespace SF
