/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_weightBuilder.h
/// @brief IBM ghost-cell 分类、ILW样本缓存和插值权重构建。

#include "core/mesh/SF_dimension.h"
#include "SF_ibmConfig.h"
#include "SF_weightTypes.h"
#include "SF_field.h"
#include "SF_ilwClosure.h"
#include "SF_ibmTopology.h"
#include "geoProcessing/SF_STLGeometry.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace SF {
namespace IBM {
namespace GhostIBM {

using GeoProcessing::SDFResult;
using GeoProcessing::STLGeometry;
using GeoProcessing::norm2;

/// @brief IBM几何前处理分段耗时统计。
struct BuildTiming {
    double sdfSeconds = 0.0;
    double ghostLayerSeconds = 0.0;
    double classifySeconds = 0.0;
    double ilwSampleCacheSeconds = 0.0;
    double ilwFluidSampleSeconds = 0.0;
    double ilwNormalSampleSeconds = 0.0;
    double ilwPlanSeconds = 0.0;
    double donorWeightSeconds = 0.0;
    double totalSeconds = 0.0;
    int ilwSampleCells = 0;
};

/// @brief 根据STL/SDF生成IBM单元分类和ghost插值权重。
///
/// 该类只负责几何拓扑和插值权重构建，不负责时间推进和通量计算。
/// 分类契约:
/// - FLUID_CELL: 真实流体点，参与RHS/CFL/源项/时间推进；
/// - IBM_GHOST_CELL: 边界闭合点，由ILW/IBM每步重建，不参与推进；
/// - SOLID_CELL: 固体内部点，不参与流体求解。
class WeightBuilder {
public:
    /// @brief 构建IBM分类和ghost-cell权重。
    /// @param field 结构网格场；CellFlag、IBMFluidMask会被更新。
    /// @param geometry 已加载的STL几何和SDF查询器；IBM拓扑几何写入本类geometry()。
    void build(Field& field,
               const STLGeometry& geometry,
               const IBMRuntimeConfig& config) {
        config_ = config;
        using Clock = std::chrono::steady_clock;
        auto secondsSince = [](Clock::time_point start, Clock::time_point end) {
            return std::chrono::duration<double>(end - start).count();
        };

        const auto totalStart = Clock::now();
        weights_.clear();
        counts_ = {};
        timing_ = {};
        ilwStorage_.setup((size_t)field.TotalSize());
        geometry_.setup(field.MX(), field.MY(), field.MZ());

        std::vector<SDFResult> sdfCache((size_t)field.TotalSize());
        std::vector<unsigned char> inside((size_t)field.TotalSize(), 0);
        std::vector<int> ghostLayer((size_t)field.TotalSize(), 0);

        if (Math::activeDimensionCount() == 2 &&
            Math::singleInactiveDirectionIndex() >= 0) {
            broadcast("IBM geometry: ",
                      "using 2D projected STL distance for EMPTY direction");
        }

        auto phaseStart = Clock::now();
        for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
            for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                    SDFResult sdf = signedDistanceForActiveGeometry(
                        field, geometry, i, j, k);
                    sdfCache[(size_t)field.getIdx(i, j, k)] = sdf;
                    inside[(size_t)field.getIdx(i, j, k)] = sdf.inside ? 1 : 0;
                    geometry_.signedDistance(i, j, k) = sdf.signedDistance;
                }
            }
        }
        timing_.sdfSeconds = secondsSince(phaseStart, Clock::now());

        phaseStart = Clock::now();
        const int ghostDepth = config_.requiredGhostLayers;
        for (int layer = 1; layer <= ghostDepth; ++layer) {
            for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
                for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                    for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                        const size_t id = (size_t)field.getIdx(i, j, k);
                        if (!inside[id] || ghostLayer[id] != 0) continue;
                        const bool touchesLayer = (layer == 1)
                            ? touchesFluidNeighbor(field, inside, i, j, k)
                            : touchesGhostLayerNeighbor(field, ghostLayer, layer - 1, i, j, k);
                        if (touchesLayer) ghostLayer[id] = layer;
                    }
                }
            }
        }
        timing_.ghostLayerSeconds = secondsSince(phaseStart, Clock::now());

        phaseStart = Clock::now();
        for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
            for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                    const size_t id = (size_t)field.getIdx(i, j, k);
                    const SDFResult& sdf = sdfCache[id];
                    if (!inside[id]) {
                        field.setIBMFluidMask(i, j, k, true);
                        field.CellFlag(i, j, k) = FLUID_CELL;
                        ++counts_.fluid;
                        continue;
                    }

                    const bool isGhost = ghostLayer[id] > 0;
                    field.setIBMFluidMask(i, j, k, false);
                    field.CellFlag(i, j, k) = isGhost ? IBM_GHOST_CELL : SOLID_CELL;
                    const Point imagePoint = mirrorPoint(sdf);
                    geometry_.setGeometry(
                        i, j, k,
                        Vector3(sdf.closestPoint.x, sdf.closestPoint.y, sdf.closestPoint.z),
                        Vector3(imagePoint.x, imagePoint.y, imagePoint.z),
                        Vector3(sdf.normal.x, sdf.normal.y, sdf.normal.z),
                        Vector3(0.0, 0.0, 0.0));
                    if (isGhost) ++counts_.ghost;
                    else ++counts_.solid;
                }
            }
        }
        timing_.classifySeconds = secondsSince(phaseStart, Clock::now());

        phaseStart = Clock::now();
        buildILWLocalSampleCache(field, geometry_, ilwStorage_, timing_);
        timing_.ilwSampleCacheSeconds = secondsSince(phaseStart, Clock::now());

        ilwStats_ = {};
        if (config_.ilwEnabled) {
            broadcast("IBM ILW preprocessing: ", "start");
            phaseStart = Clock::now();
            ilwStats_ = SF::IBM::GhostILW::preprocessIBMGeometry(
                field, ilwStorage_, geometry_, config_,
                config_.requestedTaylorOrder());
            timing_.ilwPlanSeconds = secondsSince(phaseStart, Clock::now());
            broadcast("IBM ILW preprocessing: ", "complete");
        } else {
            ilwStorage_.clearPlans();
        }

        phaseStart = Clock::now();
        for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
            for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                    if (field.CellFlag(i, j, k) != IBM_GHOST_CELL) continue;

                    const SDFResult& sdf = sdfCache[(size_t)field.getIdx(i, j, k)];

                    GhostCellWeight w;
                    w.ghostI = i;
                    w.ghostJ = j;
                    w.ghostK = k;
                    w.boundaryIntercept = sdf.closestPoint;
                    w.wallNormal = sdf.normal;
                    w.wallDistance = sdf.distance;
                    w.imagePoint = mirrorPoint(sdf);
                    w.donors = interpolationWeights(field, w.imagePoint, i, j, k);
                    if (w.donors.empty()) {
                        DonorWeight donor;
                        nearestFluidCell(field, w.imagePoint, donor.i, donor.j, donor.k);
                        donor.weight = 1.0;
                        w.donors.push_back(donor);
                    }
                    weights_.push_back(w);
                }
            }
        }
        timing_.donorWeightSeconds = secondsSince(phaseStart, Clock::now());
        timing_.totalSeconds = secondsSince(totalStart, Clock::now());
    }

    /// @brief 基于分解前已完成的IBM分类构建当前patch的局部闭合数据。
    ///
    /// 该路径不重新分类物理点，只建立ILW样本、局部投影计划和ghost donor。
    /// donor允许位于已映射的MPI内部halo；找不到donor时直接失败。
    /// @param field 已携带IBM分类和内部halo标记的结构patch。
    /// @return 所有IBM ghost均建立有效donor时返回true。
    bool buildFromPreclassified(Field& field,
                                const IBMGeometry& geometry,
                                const IBMRuntimeConfig& config) {
        config_ = config;
        using Clock = std::chrono::steady_clock;
        auto secondsSince = [](Clock::time_point start, Clock::time_point end) {
            return std::chrono::duration<double>(end - start).count();
        };

        const auto totalStart = Clock::now();
        weights_.clear();
        counts_ = {};
        timing_ = {};
        ilwStorage_.setup((size_t)field.TotalSize());

        for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
            for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                    switch (field.CellFlag(i, j, k)) {
                        case FLUID_CELL: ++counts_.fluid; break;
                        case IBM_GHOST_CELL: ++counts_.ghost; break;
                        case SOLID_CELL: ++counts_.solid; break;
                        default:
                            std::cerr
                                << "[SF FATAL] invalid preclassified IBM cell type "
                                << field.CellFlag(i, j, k) << " at ("
                                << i << "," << j << "," << k << ")"
                                << std::endl;
                            return false;
                    }
                }
            }
        }

        auto phaseStart = Clock::now();
        buildILWLocalSampleCache(field, geometry, ilwStorage_, timing_);
        timing_.ilwSampleCacheSeconds =
            secondsSince(phaseStart, Clock::now());

        ilwStats_ = {};
        if (config_.ilwEnabled) {
            broadcast("IBM ILW preprocessing: ", "start");
            phaseStart = Clock::now();
            ilwStats_ = SF::IBM::GhostILW::preprocessIBMGeometry(
                field, ilwStorage_, geometry, config_,
                config_.requestedTaylorOrder(),
                true);
            timing_.ilwPlanSeconds =
                secondsSince(phaseStart, Clock::now());
            broadcast("IBM ILW preprocessing: ", "complete");
        } else {
            ilwStorage_.clearPlans();
        }

        phaseStart = Clock::now();
        for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
            for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                    if (field.CellFlag(i, j, k) != IBM_GHOST_CELL) continue;
                    if (!geometry.hasGeometry(i, j, k)) {
                        std::cerr
                            << "[SF FATAL] preclassified IBM ghost lacks geometry at ("
                            << i << "," << j << "," << k << ")"
                            << std::endl;
                        return false;
                    }

                    GhostCellWeight w;
                    w.ghostI = i;
                    w.ghostJ = j;
                    w.ghostK = k;
                    const Vector3 wall = geometry.wallPoint(i, j, k);
                    const Vector3 image = geometry.imagePoint(i, j, k);
                    const Vector3 normal = geometry.wallNormal(i, j, k);
                    w.boundaryIntercept = {wall.x, wall.y, wall.z};
                    w.imagePoint = {image.x, image.y, image.z};
                    w.wallNormal = {normal.x, normal.y, normal.z};
                    w.wallDistance = std::abs(geometry.signedDistance(i, j, k));
                    w.donors = interpolationWeights(
                        field, w.imagePoint, i, j, k);
                    if (w.donors.empty()) {
                        std::cerr
                            << "[SF FATAL] no fluid donor for preclassified IBM ghost at ("
                            << i << "," << j << "," << k
                            << "), imagePoint=(" << image.x << ","
                            << image.y << "," << image.z << ")"
                            << std::endl;
                        return false;
                    }
                    weights_.push_back(std::move(w));
                }
            }
        }
        timing_.donorWeightSeconds =
            secondsSince(phaseStart, Clock::now());
        timing_.totalSeconds =
            secondsSince(totalStart, Clock::now());
        return true;
    }

    /// @brief 返回最近一次分类的统计信息。
    /// @return fluid/ghost/solid三类物理点数量。
    const ClassificationCounts& counts() const { return counts_; }

    /// @brief 返回所有IBM ghost点的插值权重。
    /// @return ghost-cell权重数组。
    const std::vector<GhostCellWeight>& weights() const { return weights_; }

    /// @brief 返回ILW几何前处理统计。
    const SF::IBM::GhostILW::PreprocessStats& ilwStats() const { return ilwStats_; }

    /// @brief 返回本patch的ILW专属样本与预计算计划。
    const SF::IBM::GhostILW::Storage& ilwStorage() const { return ilwStorage_; }

    /// @brief 返回最近一次IBM几何前处理分段耗时。
    const BuildTiming& timing() const { return timing_; }

    /// @brief 返回本 patch 的 IBM 拓扑几何（ghost 层、壁面/镜像点、法向）。
    const IBMGeometry& geometry() const { return geometry_; }

private:
    std::vector<GhostCellWeight> weights_;
    ClassificationCounts counts_;
    SF::IBM::GhostILW::PreprocessStats ilwStats_;
    SF::IBM::GhostILW::Storage ilwStorage_;
    IBMRuntimeConfig config_;
    BuildTiming timing_;
    IBMGeometry geometry_;

    struct CachedCandidate {
        int idx = 0;
        double d2 = 0.0;
        double score = 0.0;
        double alignment = 1.0;
        bool aligned = true;
    };

    static CachedCandidate makeDistanceCandidate(int idx, double d2) {
        CachedCandidate c;
        c.idx = idx;
        c.d2 = d2;
        c.score = d2;
        return c;
    }

    static CachedCandidate makeNormalCandidate(int idx,
                                               double d2,
                                               double score,
                                               double alignment,
                                               bool aligned) {
        CachedCandidate c;
        c.idx = idx;
        c.d2 = d2;
        c.score = score;
        c.alignment = alignment;
        c.aligned = aligned;
        return c;
    }

    static Point inactiveDirectionPoint() {
        const auto n = Math::inactiveDirectionNormalVector();
        return Point{n[0], n[1], n[2]};
    }

    static Point activeProjection(const Point& p, const Point& inactive) {
        return p - inactive * dot(p, inactive);
    }

    static Point closestActivePointOnSegment(const Point& p,
                                             const Point& a,
                                             const Point& b,
                                             const Point& inactive,
                                             double& activeDistance2) {
        const Point ap = activeProjection(p - a, inactive);
        const Point ab = activeProjection(b - a, inactive);
        const double ab2 = norm2(ab);
        double t = 0.0;
        if (ab2 > 1.0e-24) {
            t = std::max(0.0, std::min(1.0, dot(ap, ab) / ab2));
        }

        Point c = a + (b - a) * t;
        const double inactiveDelta = dot(p - c, inactive);
        c = c + inactive * inactiveDelta;
        activeDistance2 = norm2(activeProjection(p - c, inactive));
        return c;
    }

    static SDFResult signedDistance2DProjected(const STLGeometry& geometry,
                                               const Point& p) {
        SDFResult out;
        const Point inactive = inactiveDirectionPoint();
        const bool inside = geometry.contains(p);
        double best2 = std::numeric_limits<double>::max();
        Point bestPoint;
        Point bestNormal;
        std::size_t bestTriangle = 0;
        bool found = false;

        const auto& triangles = geometry.triangles();
        for (std::size_t t = 0; t < triangles.size(); ++t) {
            const auto& tri = triangles[t];
            const double inactiveDot = dot(tri.normal, inactive);
            Point activeNormal = tri.normal - inactive * inactiveDot;
            if (norm2(activeNormal) <= 1.0e-12) continue;
            activeNormal = normalize(activeNormal);

            for (int e = 0; e < 3; ++e) {
                const Point& a = tri.v[e];
                const Point& b = tri.v[(e + 1) % 3];
                const Point activeEdge = activeProjection(b - a, inactive);
                if (norm2(activeEdge) <= 1.0e-24) continue;

                double d2 = 0.0;
                Point c = closestActivePointOnSegment(
                    p, a, b, inactive, d2);
                if (d2 >= best2) continue;

                best2 = d2;
                bestPoint = c;
                bestNormal = activeNormal;
                bestTriangle = t;
                found = true;
            }
        }

        if (!found) {
            std::cerr
                << "[SF FATAL] 2D IBM projected geometry found no side "
                << "facets with active-plane normals. Check the STL or use "
                << "a 3D/non-EMPTY setup for this geometry."
                << std::endl;
            std::exit(1);
        }

        out.inside = inside;
        out.distance = std::sqrt(std::max(0.0, best2));
        out.signedDistance = inside ? -out.distance : out.distance;
        out.closestPoint = bestPoint;
        out.normal = bestNormal;
        out.triangleIndex = bestTriangle;

        const Point activeDelta =
            activeProjection(p - out.closestPoint, inactive);
        if (inside && dot(activeDelta, out.normal) > 0.0) {
            out.normal = -1.0 * out.normal;
        }
        if (!inside && dot(activeDelta, out.normal) < 0.0) {
            out.normal = -1.0 * out.normal;
        }
        return out;
    }

    static SDFResult signedDistanceForActiveGeometry(const Field& field,
                                                     const STLGeometry& geometry,
                                                     int i,
                                                     int j,
                                                     int k) {
        const Point p{field.X(i, j, k), field.Y(i, j, k), field.Z(i, j, k)};
        if (Math::activeDimensionCount() == 2 &&
            Math::singleInactiveDirectionIndex() >= 0) {
            return signedDistance2DProjected(geometry, p);
        }
        return geometry.signedDistance(p);
    }

    void buildILWLocalSampleCache(
        const Field& field,
        const IBMGeometry& geometry,
        SF::IBM::GhostILW::Storage& storage,
        BuildTiming& timing) const {
        using Clock = std::chrono::steady_clock;
        auto secondsSince = [](Clock::time_point start, Clock::time_point end) {
            return std::chrono::duration<double>(end - start).count();
        };

        storage.clearSamples();
        for (int k = field.NG(); k < field.NG() + field.NZ(); ++k) {
            for (int j = field.NG(); j < field.NG() + field.NY(); ++j) {
                for (int i = field.NG(); i < field.NG() + field.NX(); ++i) {
                    if (field.CellFlag(i, j, k) != IBM_GHOST_CELL) continue;
                    if (!geometry.hasGeometry(i, j, k)) continue;

                    ++timing.ilwSampleCells;
                    auto sampleStart = Clock::now();
                    const std::vector<int> fluidSamples =
                        collectILWFluidSamples(field, geometry, i, j, k);
                    timing.ilwFluidSampleSeconds +=
                        secondsSince(sampleStart, Clock::now());

                    sampleStart = Clock::now();
                    const std::vector<int> normalSamples =
                        collectILWNormalSamples(field, geometry, i, j, k);
                    timing.ilwNormalSampleSeconds +=
                        secondsSince(sampleStart, Clock::now());

                    const size_t cell = (size_t)field.getIdx(i, j, k);
                    storage.setFluidSamples(cell, fluidSamples);
                    storage.setNormalSamples(cell, normalSamples);
                }
            }
        }
    }

    /// @brief 收集IBM单元的ILW流体侧样本（BFS扩展盒搜索）。
    ///
    /// 对IBM非流体单元，以BFS逐层向外扩展收集法向侧的FLUID_CELL，
    /// 取距离最近的至多256个作为ILW拟合的流体样本。
    static std::vector<int> collectILWFluidSamples(const Field& field,
                                                  const IBMGeometry& geometry,
                                                  int ibmI,
                                                  int ibmJ,
                                                  int ibmK) {
        const Vector3 wall = geometry.wallPoint(ibmI, ibmJ, ibmK);
        const Vector3 normal = geometry.wallNormal(ibmI, ibmJ, ibmK);

        // BFS收集局部候选，排序后保留足够样本供高阶2D/3D模板筛选。
        std::vector<CachedCandidate> candidates =
            collectILWFluidSamplesBFS(field, ibmI, ibmJ, ibmK,
                                      wall, normal, 512);

        std::sort(candidates.begin(), candidates.end(),
                  [](const CachedCandidate& a, const CachedCandidate& b) { return a.d2 < b.d2; });
        if ((int)candidates.size() > 512) candidates.resize(512);

        std::vector<int> samples;
        samples.reserve(candidates.size());
        for (const auto& c : candidates) samples.push_back(c.idx);
        return samples;
    }

    /// @brief 收集IBM单元的ILW法向几何样本（BFS扩展盒搜索）。
    ///
    /// 收集邻近IBM几何点（非流体且有IBM几何），用于曲率张量估计，
    /// 至多保留256个。
    std::vector<int> collectILWNormalSamples(const Field& field,
                                                   const IBMGeometry& geometry,
                                                   int ibmI,
                                                   int ibmJ,
                                                   int ibmK) const {
        const Vector3 wall = geometry.wallPoint(ibmI, ibmJ, ibmK);
        const Vector3 normal = normalize(geometry.wallNormal(ibmI, ibmJ, ibmK));

        std::vector<CachedCandidate> candidates =
            collectILWNormalSamplesBFS(field, geometry, ibmI, ibmJ, ibmK,
                                       wall, normal,
                                       config_.normalSearchTargetCandidates);

        std::sort(candidates.begin(), candidates.end(),
                  [](const CachedCandidate& a, const CachedCandidate& b) {
                      if (a.aligned != b.aligned) return a.aligned;
                      if (a.score != b.score) return a.score < b.score;
                      if (a.d2 != b.d2) return a.d2 < b.d2;
                      return a.idx < b.idx;
                  });
        if ((int)candidates.size() > config_.normalSearchKeepSamples) {
            candidates.resize((size_t)config_.normalSearchKeepSamples);
        }

        std::vector<int> samples;
        samples.reserve(candidates.size());
        for (const auto& c : candidates) samples.push_back(c.idx);
        return samples;
    }

    /// @brief BFS扩展盒搜索：从IBM单元出发，逐层向外扩展收集流体侧样本。
    ///
    /// 替代原先的全局 O(N) 暴力扫描，复杂度降为 O(R³)，
    /// 其中 R 为收集到足够样本所需的扩展半径。
    /// @param field 结构网格场。
    /// @param ibmI,ibmJ,ibmK IBM ghost/solid单元索引。
    /// @param wall IBM壁面截距点。
    /// @param normal 从固体指向流体的壁面单位法向。
    /// @param targetCount 期望收集的最小候选数。
    /// @return 按距离排序的候选样本列表。
    static std::vector<CachedCandidate> collectILWFluidSamplesBFS(const Field& field,
                                                                   int ibmI,
                                                                   int ibmJ,
                                                                   int ibmK,
                                                                   const Vector3& wall,
                                                                   const Vector3& normal,
                                                                   int targetCount) {
        std::vector<CachedCandidate> candidates;
        candidates.reserve((size_t)targetCount);

        const int ng = field.NG();
        const int nx = field.NX();
        const int ny = field.NY();
        const int nz = field.NZ();
        const int maxRadius = std::max({nx, ny, nz});

        for (int r = 1; r <= maxRadius; ++r) {
            const int kMin = std::max(0, ibmK - r);
            const int kMax = std::min(field.MZ() - 1, ibmK + r);
            const int jMin = std::max(0, ibmJ - r);
            const int jMax = std::min(field.MY() - 1, ibmJ + r);
            const int iMin = std::max(0, ibmI - r);
            const int iMax = std::min(field.MX() - 1, ibmI + r);

            for (int kk = kMin; kk <= kMax; ++kk) {
                const bool kShell = (kk == kMin || kk == kMax);
                for (int jj = jMin; jj <= jMax; ++jj) {
                    const bool jShell = (jj == jMin || jj == jMax);
                    for (int ii = iMin; ii <= iMax; ++ii) {
                        // 只检查当前壳层：至少一个分量在边界上
                        if (!kShell && !jShell &&
                            ii != iMin && ii != iMax) continue;
                        if (!isAvailable(field, ii, jj, kk)) continue;
                        if (field.CellFlag(ii, jj, kk) != FLUID_CELL) continue;
                        const Vector3 p(field.X(ii, jj, kk),
                                        field.Y(ii, jj, kk),
                                        field.Z(ii, jj, kk));
                        const Vector3 rvec = p - wall;
                        if (dot(rvec, normal) <= 1e-12) continue;
                        const double d2 = dot(rvec, rvec);
                        candidates.push_back(makeDistanceCandidate(
                            field.getIdx(ii, jj, kk), d2));
                    }
                }
            }

            if ((int)candidates.size() >= targetCount) break;
        }
        return candidates;
    }

    /// @brief BFS扩展盒搜索：从IBM单元出发，逐层向外收集邻近IBM几何点。
    ///
    /// 替代原先的全局 O(N) 暴力扫描。IBM法向样本通常聚集在壁面附近，
    /// BFS 可在很小的半径内完成收集。
    /// @param field 结构网格场。
    /// @param ibmI,ibmJ,ibmK IBM ghost/solid单元索引。
    /// @param wall IBM壁面截距点。
    /// @param targetCount 期望收集的最小候选数。
    /// @return 按距离排序的候选样本列表。
    std::vector<CachedCandidate> collectILWNormalSamplesBFS(const Field& field,
                                                                    const IBMGeometry& geometry,
                                                                    int ibmI,
                                                                    int ibmJ,
                                                                    int ibmK,
                                                                    const Vector3& wall,
                                                                    const Vector3& normal,
                                                                    int targetCount) const {
        std::vector<CachedCandidate> candidates;
        candidates.reserve((size_t)targetCount);

        const int ng = field.NG();
        const int nx = field.NX();
        const int ny = field.NY();
        const int nz = field.NZ();
        const int domainMaxRadius = std::max({nx, ny, nz});
        const int maxRadius = (config_.normalSearchMaxLayers > 0)
            ? std::min(domainMaxRadius, config_.normalSearchMaxLayers)
            : domainMaxRadius;
        const int minRadius = config_.normalSearchMinLayers;
        const double alignmentThreshold =
            config_.normalAlignmentThreshold();

        for (int r = 1; r <= maxRadius; ++r) {
            const int kMin = std::max(0, ibmK - r);
            const int kMax = std::min(field.MZ() - 1, ibmK + r);
            const int jMin = std::max(0, ibmJ - r);
            const int jMax = std::min(field.MY() - 1, ibmJ + r);
            const int iMin = std::max(0, ibmI - r);
            const int iMax = std::min(field.MX() - 1, ibmI + r);

            for (int kk = kMin; kk <= kMax; ++kk) {
                const bool kShell = (kk == kMin || kk == kMax);
                for (int jj = jMin; jj <= jMax; ++jj) {
                    const bool jShell = (jj == jMin || jj == jMax);
                    for (int ii = iMin; ii <= iMax; ++ii) {
                        if (!kShell && !jShell &&
                            ii != iMin && ii != iMax) continue;
                        if (!isAvailable(field, ii, jj, kk)) continue;
                        if (field.CellFlag(ii, jj, kk) == FLUID_CELL) continue;
                        if (!geometry.hasGeometry(ii, jj, kk)) continue;
                        const Vector3 wp = geometry.wallPoint(ii, jj, kk);
                        const Vector3 rvec = wp - wall;
                        const double d2 = dot(rvec, rvec);
                        const Vector3 candidateNormal =
                            normalize(geometry.wallNormal(ii, jj, kk));
                        const double alignment =
                            std::abs(dot(normal, candidateNormal));
                        const bool aligned = alignment >= alignmentThreshold;
                        const double normalOffset = std::abs(dot(rvec, normal));
                        const double misalignment = std::max(0.0, 1.0 - alignment);
                        const double score =
                            d2 * (1.0 + 4.0 * misalignment)
                            + normalOffset * normalOffset;
                        candidates.push_back(makeNormalCandidate(
                            field.getIdx(ii, jj, kk), d2, score,
                            alignment, aligned));
                    }
                }
            }

            if (r >= minRadius && (int)candidates.size() >= targetCount) break;
        }
        return candidates;
    }

    static bool touchesFluidNeighbor(const Field& field,
                                     const std::vector<unsigned char>& inside,
                                     int i, int j, int k) {
        for (int dk = -1; dk <= 1; ++dk) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    if (di == 0 && dj == 0 && dk == 0) continue;
                    int ii = i + di;
                    int jj = j + dj;
                    int kk = k + dk;
                    if (!isAvailable(field, ii, jj, kk)) continue;
                    if (!inside[(size_t)field.getIdx(ii, jj, kk)]) return true;
                }
            }
        }
        return false;
    }

    static bool touchesGhostLayerNeighbor(const Field& field,
                                          const std::vector<int>& ghostLayer,
                                          int targetLayer,
                                          int i, int j, int k) {
        for (int dk = -1; dk <= 1; ++dk) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    if (di == 0 && dj == 0 && dk == 0) continue;
                    int ii = i + di;
                    int jj = j + dj;
                    int kk = k + dk;
                    if (!isAvailable(field, ii, jj, kk)) continue;
                    if (ghostLayer[(size_t)field.getIdx(ii, jj, kk)] == targetLayer) return true;
                }
            }
        }
        return false;
    }

    static Point mirrorPoint(const SDFResult& sdf) {
        constexpr double eps = 1e-10;
        return sdf.closestPoint + sdf.normal * (sdf.distance + eps);
    }

    static std::vector<DonorWeight> interpolationWeights(const Field& field,
                                                         const Point& p,
                                                         int ghostI,
                                                         int ghostJ,
                                                         int ghostK) {
        std::vector<DonorWeight> localDonors = localInverseDistanceWeights(field, p, ghostI, ghostJ, ghostK);
        if (!localDonors.empty()) return localDonors;

        int i0 = 0, j0 = 0, k0 = 0;
        double tx = 0.0, ty = 0.0, tz = 0.0;
        if (!locateAxis(field, 0, p.x, i0, tx)) return {};
        if (!locateAxis(field, 1, p.y, j0, ty)) return {};
        if (!locateAxis(field, 2, p.z, k0, tz)) return {};

        std::vector<DonorWeight> donors;
        donors.reserve(8);
        for (int dk = 0; dk <= 1; ++dk) {
            for (int dj = 0; dj <= 1; ++dj) {
                for (int di = 0; di <= 1; ++di) {
                    int ii = i0 + di;
                    int jj = j0 + dj;
                    int kk = k0 + dk;
                    if (!isAvailable(field, ii, jj, kk)) continue;
                    if (field.CellFlag(ii, jj, kk) != FLUID_CELL) return {};

                    double wx = di ? tx : (1.0 - tx);
                    double wy = dj ? ty : (1.0 - ty);
                    double wz = dk ? tz : (1.0 - tz);
                    double weight = wx * wy * wz;
                    if (weight > 1e-14) donors.push_back({ii, jj, kk, weight});
                }
            }
        }
        normalizeWeights(donors);
        return donors;
    }

    static std::vector<DonorWeight> localInverseDistanceWeights(const Field& field,
                                                                const Point& p,
                                                                int ghostI,
                                                                int ghostJ,
                                                                int ghostK) {
        struct Candidate {
            int i = 0;
            int j = 0;
            int k = 0;
            double d2 = 0.0;
        };

        std::vector<Candidate> candidates;
        candidates.reserve(64);
        constexpr int radius = 4;

        for (int kk = ghostK - radius; kk <= ghostK + radius; ++kk) {
            for (int jj = ghostJ - radius; jj <= ghostJ + radius; ++jj) {
                for (int ii = ghostI - radius; ii <= ghostI + radius; ++ii) {
                    if (!isAvailable(field, ii, jj, kk)) continue;
                    if (field.CellFlag(ii, jj, kk) != FLUID_CELL) continue;

                    Point q{field.X(ii, jj, kk), field.Y(ii, jj, kk), field.Z(ii, jj, kk)};
                    candidates.push_back({ii, jj, kk, norm2(q - p)});
                }
            }
        }

        if (candidates.empty()) return {};
        std::sort(candidates.begin(), candidates.end(),
                  [](const Candidate& a, const Candidate& b) { return a.d2 < b.d2; });

        std::vector<DonorWeight> donors;
        donors.reserve(8);
        const int n = std::min<int>(8, (int)candidates.size());
        for (int m = 0; m < n; ++m) {
            const double weight = 1.0 / (candidates[m].d2 + 1e-20);
            donors.push_back({candidates[m].i, candidates[m].j, candidates[m].k, weight});
        }
        normalizeWeights(donors);
        return donors;
    }

    static bool locateAxis(const Field& field, int axis, double x, int& i0, double& t) {
        int ng = field.NG();
        int n = (axis == 0) ? field.NX() : (axis == 1 ? field.NY() : field.NZ());
        if (n <= 1) {
            i0 = ng;
            t = 0.0;
            return true;
        }

        auto coord = [&](int idx) {
            if (axis == 0) return field.X(idx, ng, ng);
            if (axis == 1) return field.Y(ng, idx, ng);
            return field.Z(ng, ng, idx);
        };

        int first = ng;
        int last = ng + n - 1;
        for (int i = first; i < last; ++i) {
            double a = coord(i);
            double b = coord(i + 1);
            double lo = std::min(a, b) - 1e-12;
            double hi = std::max(a, b) + 1e-12;
            if (x < lo || x > hi) continue;
            double denom = b - a;
            i0 = i;
            t = (std::abs(denom) > 1e-14) ? (x - a) / denom : 0.0;
            t = std::max(0.0, std::min(1.0, t));
            return true;
        }
        return false;
    }

    static bool isInterior(const Field& field, int i, int j, int k) {
        int ng = field.NG();
        return i >= ng && i < ng + field.NX()
            && j >= ng && j < ng + field.NY()
            && k >= ng && k < ng + field.NZ();
    }

    static bool isAvailable(const Field& field, int i, int j, int k) {
        if (i < 0 || i >= field.MX() ||
            j < 0 || j >= field.MY() ||
            k < 0 || k >= field.MZ()) {
            return false;
        }
        return isInterior(field, i, j, k)
            || field.isCommunicationHalo(i, j, k);
    }

    static void normalizeWeights(std::vector<DonorWeight>& donors) {
        double sum = 0.0;
        for (const auto& d : donors) sum += d.weight;
        if (sum <= 1e-14) return;
        for (auto& d : donors) d.weight /= sum;
    }

    /// @brief BFS扩展盒搜索最近的流体单元。
    ///
    /// 先用 locateAxis 沿各轴定位目标点的网格索引，再从该位置
    /// BFS 逐层向外搜索 FLUID_CELL，替代原先的全局 O(N) 暴力扫描。
    /// @param field 结构网格场。
    /// @param p 目标空间坐标。
    /// @param oi,oj,ok 输出最近流体单元的i/j/k索引。
    static void nearestFluidCell(const Field& field, const Point& p, int& oi, int& oj, int& ok) {
        int ng = field.NG();
        oi = ng; oj = ng; ok = ng;

        // 用 locateAxis 沿各轴定位，O(NX+NY+NZ) 替代 O(NX*NY*NZ)
        int i0 = ng, j0 = ng, k0 = ng;
        double tx = 0.0, ty = 0.0, tz = 0.0;
        locateAxis(field, 0, p.x, i0, tx);
        locateAxis(field, 1, p.y, j0, ty);
        locateAxis(field, 2, p.z, k0, tz);

        double bestDist = std::numeric_limits<double>::max();
        const int nx = field.NX();
        const int ny = field.NY();
        const int nz = field.NZ();
        const int maxRadius = std::max({nx, ny, nz});

        for (int r = 0; r <= maxRadius; ++r) {
            bool foundAny = false;
            const int kMin = std::max(ng, k0 - r);
            const int kMax = std::min(ng + nz - 1, k0 + r);
            const int jMin = std::max(ng, j0 - r);
            const int jMax = std::min(ng + ny - 1, j0 + r);
            const int iMin = std::max(ng, i0 - r);
            const int iMax = std::min(ng + nx - 1, i0 + r);

            for (int kk = kMin; kk <= kMax; ++kk) {
                const bool kShell = (r == 0) || (kk == kMin || kk == kMax);
                for (int jj = jMin; jj <= jMax; ++jj) {
                    const bool jShell = (r == 0) || (jj == jMin || jj == jMax);
                    for (int ii = iMin; ii <= iMax; ++ii) {
                        if (r > 0 && !kShell && !jShell &&
                            ii != iMin && ii != iMax) continue;
                        if (field.CellFlag(ii, jj, kk) != FLUID_CELL) continue;
                        foundAny = true;
                        Point q{field.X(ii, jj, kk), field.Y(ii, jj, kk), field.Z(ii, jj, kk)};
                        double d2 = norm2(p - q);
                        if (d2 < bestDist) {
                            bestDist = d2;
                            oi = ii; oj = jj; ok = kk;
                        }
                    }
                }
            }
            if (foundAny) break;
        }
    }
};

} // namespace GhostIBM
} // namespace IBM
} // namespace SF
