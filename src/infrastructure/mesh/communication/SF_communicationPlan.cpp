/// @file SF_communicationPlan.cpp
/// @brief 网格接口通信计划与 donor/receiver 拓扑描述。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_communicationPlan.h"
#include "SF_MultiBlockMesh.h"
#include "core/field/SF_field.h"
#include "methods/numerics/structured/SF_structured.h"
#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace SF {
namespace MeshCommunication {
namespace {

struct Point3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct DonorCell {
    std::array<int, 3> ijk{0, 0, 0};
    std::array<int, 8> indices{};
    std::array<Point3, 8> points{};
    Point3 min;
    Point3 max;
};

struct DonorPoint {
    int blockId = -1;
    int index = -1;
    std::array<int, 3> ijk{0, 0, 0};
    Point3 point;
};

struct DonorCache {
    int blockId = -1;
    const Field* field = nullptr;
    Point3 min;
    Point3 max;
    int binsX = 1;
    int binsY = 1;
    int binsZ = 1;
    std::vector<DonorCell> cells;
    std::vector<std::vector<int>> bins;
};

struct LocatedCell {
    int donorBlock = -1;
    std::array<int, 3> donorCellIJK{0, 0, 0};
    std::array<double, 3> localCoord{0.0, 0.0, 0.0};
    double residual = std::numeric_limits<double>::max();
};

struct GlobalFaceKey {
    std::array<int, 4> ids{-1, -1, -1, -1};

    bool operator<(const GlobalFaceKey& other) const {
        return ids < other.ids;
    }
};

struct FluxParticipantKey {
    int blockId = -1;
    int direction = -1;
    int i = 0;
    int j = 0;
    int k = 0;

    bool operator<(const FluxParticipantKey& other) const {
        if (blockId != other.blockId) return blockId < other.blockId;
        if (direction != other.direction) return direction < other.direction;
        if (i != other.i) return i < other.i;
        if (j != other.j) return j < other.j;
        return k < other.k;
    }

    bool operator==(const FluxParticipantKey& other) const {
        return blockId == other.blockId &&
               direction == other.direction &&
               i == other.i &&
               j == other.j &&
               k == other.k;
    }
};

struct BoundaryFaceCandidate {
    int blockId = -1;
    int faceId = -1;
    std::array<int, 4> pointIds{-1, -1, -1, -1};
    std::array<int, 4> globalPointIds{-1, -1, -1, -1};
};

bool boundarySide(const MeshBlockField& block,
                  const BoundaryFaceCandidate& candidate,
                  int& axis,
                  int& side,
                  std::string& error);

int matchingCorner(const BoundaryFaceCandidate& face,
                   int globalPointId);

struct GhostKey {
    int blockId = -1;
    int index = -1;

    bool operator==(const GhostKey& other) const {
        return blockId == other.blockId && index == other.index;
    }
};

struct GhostKeyHash {
    std::size_t operator()(const GhostKey& key) const {
        std::size_t h = std::hash<int>{}(key.blockId);
        h ^= std::hash<int>{}(key.index) + 0x9e3779b97f4a7c15ULL
           + (h << 6) + (h >> 2);
        return h;
    }
};

using ConformalGhostMap =
    std::unordered_map<GhostKey, HaloCellMapping, GhostKeyHash>;

Point3 pointAt(const Field& f, int i, int j, int k) {
    return {f.X(i, j, k), f.Y(i, j, k), f.Z(i, j, k)};
}

Point3 operator+(const Point3& a, const Point3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Point3 operator-(const Point3& a, const Point3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Point3 operator*(const Point3& a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}

double norm2(const Point3& p) {
    return p.x * p.x + p.y * p.y + p.z * p.z;
}

double distance2(const Point3& a, const Point3& b) {
    return norm2(a - b);
}

double determinant(const Point3& a,
                   const Point3& b,
                   const Point3& c) {
    return a.x * (b.y * c.z - b.z * c.y)
         - a.y * (b.x * c.z - b.z * c.x)
         + a.z * (b.x * c.y - b.y * c.x);
}

double tetraVolume(const MeshPoint& a,
                   const MeshPoint& b,
                   const MeshPoint& c,
                   const MeshPoint& d) {
    const Point3 ab{b.x - a.x, b.y - a.y, b.z - a.z};
    const Point3 ac{c.x - a.x, c.y - a.y, c.z - a.z};
    const Point3 ad{d.x - a.x, d.y - a.y, d.z - a.z};
    return std::abs(determinant(ab, ac, ad)) / 6.0;
}

double hexaVolume(const CellFaceMesh& topology,
                  const MeshCell& cell) {
    std::array<const MeshPoint*, 8> p{};
    for (int n = 0; n < 8; ++n) {
        const int pointId = cell.pointIds[(size_t)n];
        if (pointId < 0 || pointId >= (int)topology.points.size()) {
            return -1.0;
        }
        p[(size_t)n] = &topology.points[(size_t)pointId];
    }
    // 沿 p000--p111 体对角线作六个确定性四面体；相邻 MPI patch 对同一
    // source cell 使用完全相同的点序，因此不会产生分区相关体积。
    return tetraVolume(*p[0], *p[1], *p[2], *p[6])
         + tetraVolume(*p[0], *p[2], *p[3], *p[6])
         + tetraVolume(*p[0], *p[3], *p[7], *p[6])
         + tetraVolume(*p[0], *p[7], *p[4], *p[6])
         + tetraVolume(*p[0], *p[4], *p[5], *p[6])
         + tetraVolume(*p[0], *p[5], *p[1], *p[6]);
}

bool buildDualVolumeFragments(
        const std::vector<MeshBlockField>& blocks,
        std::vector<std::vector<double>>& fragments) {
    fragments.clear();
    fragments.resize(blocks.size());
    for (size_t blockId = 0; blockId < blocks.size(); ++blockId) {
        const auto& topology = blocks[blockId].topology;
        const Field& field = blocks[blockId].field;
        const size_t expected =
            (size_t)field.NX() * field.NY() * field.NZ();
        if (topology.points.size() != expected || topology.cells.empty()) {
            broadcast(
                "Fatal: ",
                "GlobalDof dual-volume assembly requires non-empty hexahedral "
                "cell topology for block " + std::to_string(blockId) + ".");
            return false;
        }
        auto& volume = fragments[blockId];
        volume.assign(topology.points.size(), 0.0);
        for (const MeshCell& cell : topology.cells) {
            const double cellVolume = hexaVolume(topology, cell);
            if (!std::isfinite(cellVolume) || cellVolume <= 1.0e-300) {
                broadcast(
                    "Fatal: ",
                    "GlobalDof dual-volume assembly found a degenerate "
                    "hexahedral cell in block " + std::to_string(blockId)
                    + ".");
                return false;
            }
            for (int pointId : cell.pointIds) {
                volume[(size_t)pointId] += cellVolume / 8.0;
            }
        }
    }
    return true;
}

void includePoint(Point3& bmin, Point3& bmax, const Point3& p) {
    bmin.x = std::min(bmin.x, p.x);
    bmin.y = std::min(bmin.y, p.y);
    bmin.z = std::min(bmin.z, p.z);
    bmax.x = std::max(bmax.x, p.x);
    bmax.y = std::max(bmax.y, p.y);
    bmax.z = std::max(bmax.z, p.z);
}

bool insideBox(const Point3& bmin, const Point3& bmax,
               const Point3& p, double tol) {
    return p.x >= bmin.x - tol && p.x <= bmax.x + tol &&
           p.y >= bmin.y - tol && p.y <= bmax.y + tol &&
           p.z >= bmin.z - tol && p.z <= bmax.z + tol;
}

int clampInt(int v, int lo, int hi) {
    return std::max(lo, std::min(v, hi));
}

double clamp01(double v) {
    return std::max(0.0, std::min(1.0, v));
}

double safeFraction(double value, double lo, double hi) {
    double span = hi - lo;
    if (std::abs(span) < 1.0e-300) return 0.5;
    return clamp01((value - lo) / span);
}

int chooseBins(double extent, double maxExtent, int base, double tol) {
    if (extent <= std::max(tol, 1.0e-300) || maxExtent <= 0.0) return 1;
    int n = (int)std::round(base * extent / maxExtent);
    return clampInt(n, 1, 64);
}

int binCoord(double value, double lo, double hi, int nBins) {
    if (nBins <= 1) return 0;
    double span = hi - lo;
    if (std::abs(span) < 1.0e-300) return 0;
    int b = (int)std::floor((value - lo) / span * nBins);
    return clampInt(b, 0, nBins - 1);
}

int binIndex(const DonorCache& cache, int bx, int by, int bz) {
    return (bz * cache.binsY + by) * cache.binsX + bx;
}

bool solve3x3(const double a[3][3], const double b[3], double x[3]) {
    double det =
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if (std::abs(det) < 1.0e-24) return false;

    auto detReplace = [&](int col) {
        double m[3][3] = {
            {a[0][0], a[0][1], a[0][2]},
            {a[1][0], a[1][1], a[1][2]},
            {a[2][0], a[2][1], a[2][2]}
        };
        m[0][col] = b[0];
        m[1][col] = b[1];
        m[2][col] = b[2];
        return
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    };

    x[0] = detReplace(0) / det;
    x[1] = detReplace(1) / det;
    x[2] = detReplace(2) / det;
    return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}

void trilinearWeights(double u, double v, double w,
                      std::array<double, 8>& weights) {
    double omu = 1.0 - u;
    double omv = 1.0 - v;
    double omw = 1.0 - w;

    weights[0] = omu * omv * omw;
    weights[1] = u   * omv * omw;
    weights[2] = omu * v   * omw;
    weights[3] = u   * v   * omw;
    weights[4] = omu * omv * w;
    weights[5] = u   * omv * w;
    weights[6] = omu * v   * w;
    weights[7] = u   * v   * w;
}

Point3 weightedPoint(const DonorCell& cell,
                     const std::array<double, 8>& weights) {
    Point3 p;
    for (int n = 0; n < 8; ++n) {
        p = p + cell.points[n] * weights[n];
    }
    return p;
}

void trilinearDerivatives(const DonorCell& cell,
                          double u, double v, double w,
                          Point3& du, Point3& dv, Point3& dw) {
    double omu = 1.0 - u;
    double omv = 1.0 - v;
    double omw = 1.0 - w;

    std::array<double, 8> wu = {
        -omv * omw,  omv * omw,
        -v   * omw,  v   * omw,
        -omv * w,    omv * w,
        -v   * w,    v   * w
    };
    std::array<double, 8> wv = {
        -omu * omw, -u * omw,
         omu * omw,  u * omw,
        -omu * w,   -u * w,
         omu * w,    u * w
    };
    std::array<double, 8> ww = {
        -omu * omv, -u * omv,
        -omu * v,   -u * v,
         omu * omv,  u * omv,
         omu * v,    u * v
    };

    du = {};
    dv = {};
    dw = {};
    for (int n = 0; n < 8; ++n) {
        du = du + cell.points[n] * wu[n];
        dv = dv + cell.points[n] * wv[n];
        dw = dw + cell.points[n] * ww[n];
    }
}

double cellDiag(const DonorCell& cell) {
    return std::sqrt(distance2(cell.min, cell.max));
}

bool attemptInvert(const DonorCell& cell, const Point3& target,
                   double startU, double startV, double startW,
                   double tol, std::array<double, 3>& local,
                   std::array<double, 8>& weights,
                   double& residual) {
    double u = startU;
    double v = startV;
    double w = startW;
    double acceptTol = std::max(tol, 1.0e-10 * cellDiag(cell));
    double bestResidual = std::numeric_limits<double>::max();
    double bestU = u, bestV = v, bestW = w;

    for (int iter = 0; iter < 18; ++iter) {
        std::array<double, 8> trialWeights{};
        trilinearWeights(u, v, w, trialWeights);
        Point3 mapped = weightedPoint(cell, trialWeights);
        Point3 r = mapped - target;
        double rn = std::sqrt(norm2(r));
        if (rn < bestResidual) {
            bestResidual = rn;
            bestU = u;
            bestV = v;
            bestW = w;
        }
        if (rn <= acceptTol) break;

        Point3 du, dv, dw;
        trilinearDerivatives(cell, u, v, w, du, dv, dw);
        double a[3][3] = {
            {du.x, dv.x, dw.x},
            {du.y, dv.y, dw.y},
            {du.z, dv.z, dw.z}
        };
        double b[3] = {r.x, r.y, r.z};
        double delta[3] = {0.0, 0.0, 0.0};
        if (!solve3x3(a, b, delta)) break;

        u -= delta[0];
        v -= delta[1];
        w -= delta[2];

        if (!std::isfinite(u) || !std::isfinite(v) || !std::isfinite(w)) break;
        if (u < -0.5 || u > 1.5 || v < -0.5 || v > 1.5 ||
            w < -0.5 || w > 1.5) {
            break;
        }
    }

    double paramTol = 1.0e-6;
    if (bestU < -paramTol || bestU > 1.0 + paramTol ||
        bestV < -paramTol || bestV > 1.0 + paramTol ||
        bestW < -paramTol || bestW > 1.0 + paramTol ||
        bestResidual > acceptTol) {
        return false;
    }

    bestU = clamp01(bestU);
    bestV = clamp01(bestV);
    bestW = clamp01(bestW);
    trilinearWeights(bestU, bestV, bestW, weights);
    local = {bestU, bestV, bestW};
    residual = bestResidual;
    return true;
}

bool invertTrilinear(const DonorCell& cell, const Point3& target,
                     double tol, std::array<double, 3>& local,
                     std::array<double, 8>& weights,
                     double& residual) {
    std::array<std::array<double, 3>, 4> starts = {{
        {safeFraction(target.x, cell.min.x, cell.max.x),
         safeFraction(target.y, cell.min.y, cell.max.y),
         safeFraction(target.z, cell.min.z, cell.max.z)},
        {0.5, 0.5, 0.5},
        {0.25, 0.25, 0.25},
        {0.75, 0.75, 0.75}
    }};

    bool found = false;
    std::array<double, 3> bestLocal{0.0, 0.0, 0.0};
    std::array<double, 8> bestWeights{};
    double bestResidual = std::numeric_limits<double>::max();

    for (const auto& start : starts) {
        std::array<double, 3> trialLocal{0.0, 0.0, 0.0};
        std::array<double, 8> trialWeights{};
        double trialResidual = std::numeric_limits<double>::max();
        if (!attemptInvert(cell, target, start[0], start[1], start[2],
                           tol, trialLocal, trialWeights, trialResidual)) {
            continue;
        }
        if (trialResidual < bestResidual) {
            found = true;
            bestResidual = trialResidual;
            bestLocal = trialLocal;
            bestWeights = trialWeights;
        }
    }

    if (!found) return false;
    local = bestLocal;
    weights = bestWeights;
    residual = bestResidual;
    return true;
}

int interiorAxisSize(const Field& f, int axis) {
    if (axis == 0) return f.NX();
    if (axis == 1) return f.NY();
    return f.NZ();
}

bool lagrangeWeights(double target,
                     int start,
                     int count,
                     std::array<double, kMaxHaloInterpolationStencil>& weights) {
    weights.fill(0.0);
    if (count < 1 || count > kMaxHaloInterpolationStencil) return false;
    if (count == 1) {
        weights[0] = 1.0;
        return true;
    }

    for (int m = 0; m < count; ++m) {
        const double xm = (double)(start + m);
        double w = 1.0;
        for (int q = 0; q < count; ++q) {
            if (q == m) continue;
            const double xq = (double)(start + q);
            const double denom = xm - xq;
            if (std::abs(denom) <= 1.0e-300) return false;
            w *= (target - xq) / denom;
        }
        if (!std::isfinite(w)) return false;
        weights[(size_t)m] = w;
    }
    return true;
}

bool buildTensorInterpolation(const Field& donorField,
                              const LocatedCell& located,
                              HaloCellMapping& mapping,
                              int& tensorPointCount) {
    tensorPointCount = 1;
    const int ng = donorField.NG();
    for (int axis = 0; axis < 3; ++axis) {
        const int axisSize = interiorAxisSize(donorField, axis);
        if (axisSize <= 0) return false;

        const int lowerInterior =
            located.donorCellIJK[(size_t)axis] - ng;
        const double target =
            (double)lowerInterior + located.localCoord[(size_t)axis];
        const int count =
            axisSize <= 1
                ? 1
                : std::min(kMaxHaloInterpolationStencil, axisSize);

        int start = 0;
        if (count > 1) {
            start = (int)std::floor(target) - (count / 2 - 1);
            start = std::max(0, std::min(start, axisSize - count));
        }

        if (start < 0 || start + count > axisSize) return false;
        mapping.donorStencilStart[(size_t)axis] = start;
        mapping.donorStencilSize[(size_t)axis] = count;
        if (!lagrangeWeights(
                target,
                start,
                count,
                mapping.donorStencilWeights[(size_t)axis])) {
            return false;
        }
        tensorPointCount *= count;
    }
    return tensorPointCount > 0;
}

bool isGhostPoint(const Field& f, int i, int j, int k) {
    int ng = f.NG();
    return i < ng || i >= f.NX() + ng ||
           j < ng || j >= f.NY() + ng ||
           k < ng || k >= f.NZ() + ng;
}

int activeExteriorAxisCount(const Field& f, int i, int j, int k) {
    const int ng = f.NG();
    const int lo[3] = {ng, ng, ng};
    const int hi[3] = {
        ng + f.NX() - 1,
        ng + f.NY() - 1,
        ng + f.NZ() - 1
    };
    const int ijk[3] = {i, j, k};
    int count = 0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActive((Math::Dir)axis)) continue;
        if (ijk[axis] < lo[axis] || ijk[axis] > hi[axis]) ++count;
    }
    return count;
}

int activeGhostLayer(const Field& f, int i, int j, int k) {
    const int ng = f.NG();
    const int lo[3] = {ng, ng, ng};
    const int hi[3] = {
        ng + f.NX() - 1,
        ng + f.NY() - 1,
        ng + f.NZ() - 1
    };
    const int ijk[3] = {i, j, k};
    int layer = 0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActive((Math::Dir)axis)) continue;
        if (ijk[axis] < lo[axis]) {
            layer = std::max(layer, lo[axis] - ijk[axis]);
        } else if (ijk[axis] > hi[axis]) {
            layer = std::max(layer, ijk[axis] - hi[axis]);
        }
    }
    return layer;
}

bool localPointIJK(const Field& f,
                   int pointId,
                   std::array<int, 3>& ijk) {
    const int nx = f.NX();
    const int ny = f.NY();
    const int nz = f.NZ();
    const int nPoints = nx * ny * nz;
    if (pointId < 0 || pointId >= nPoints || nx <= 0 || ny <= 0 || nz <= 0) {
        return false;
    }
    ijk[0] = pointId % nx;
    ijk[1] = (pointId / nx) % ny;
    ijk[2] = pointId / (nx * ny);
    return true;
}

int structuredAxisSize(const Field& f, int axis) {
    if (axis == 0) return f.NX();
    if (axis == 1) return f.NY();
    return f.NZ();
}

bool validFluxIndex(const Field& f,
                    int direction,
                    int i,
                    int j,
                    int k) {
    const int ng = f.NG();
    const int nx = f.NX();
    const int ny = f.NY();
    const int nz = f.NZ();
    if (direction == 0) {
        return i >= ng - 1 && i < ng + nx &&
               j >= ng && j < ng + ny &&
               k >= ng && k < ng + nz;
    }
    if (direction == 1) {
        return i >= ng && i < ng + nx &&
               j >= ng - 1 && j < ng + ny &&
               k >= ng && k < ng + nz;
    }
    if (direction == 2) {
        return i >= ng && i < ng + nx &&
               j >= ng && j < ng + ny &&
               k >= ng - 1 && k < ng + nz;
    }
    return false;
}

FluxParticipantKey participantKey(const HaloInterfaceFace& face) {
    return {face.blockId, face.direction, face.i, face.j, face.k};
}

bool boundaryFluxIndexForPoint(const MeshBlockField& block,
                               int blockId,
                               int faceId,
                               int pointId,
                               HaloInterfaceFace& out,
                               std::string& error) {
    const Field& f = block.field;
    const CellFaceMesh& topology = block.topology;
    if (faceId < 0 || faceId >= (int)topology.faces.size()) {
        error = "invalid topology face id";
        return false;
    }

    const MeshFace& face = topology.faces[(size_t)faceId];
    const int direction = face.direction;
    if (direction < 0 || direction > 2) {
        error = "topology face has invalid direction";
        return false;
    }
    if (!Math::isDirectionActive((Math::Dir)direction)) {
        error.clear();
        return false;
    }

    bool containsPoint = false;
    for (int pid : face.pointIds) {
        if (pid == pointId) {
            containsPoint = true;
            break;
        }
    }
    if (!containsPoint) {
        error = "requested point is not on topology face";
        return false;
    }

    std::array<int, 3> firstIJK{};
    if (!localPointIJK(f, face.pointIds[0], firstIJK)) {
        error = "topology face references an invalid local point";
        return false;
    }
    const int normalCoord = firstIJK[(size_t)direction];
    for (int pid : face.pointIds) {
        std::array<int, 3> ijk{};
        if (!localPointIJK(f, pid, ijk)) {
            error = "topology face references an invalid local point";
            return false;
        }
        if (ijk[(size_t)direction] != normalCoord) {
            error = "topology face is not aligned with its declared direction";
            return false;
        }
    }

    const int axisN = structuredAxisSize(f, direction);
    const int ng = f.NG();
    int normalFaceIndex = 0;
    if (normalCoord == 0) {
        normalFaceIndex = ng - 1;
        out.residualSign = -1.0;
    } else if (normalCoord == axisN - 1) {
        normalFaceIndex = ng + axisN - 1;
        out.residualSign = 1.0;
    } else {
        error = "interface topology face is not on a structured patch boundary";
        return false;
    }

    std::array<int, 3> pointIJK{};
    if (!localPointIJK(f, pointId, pointIJK)) {
        error = "topology face references an invalid local point";
        return false;
    }

    out.blockId = blockId;
    out.direction = direction;
    out.orientation = 1.0;
    if (direction == 0) {
        out.i = normalFaceIndex;
        out.j = ng + pointIJK[1];
        out.k = ng + pointIJK[2];
    } else if (direction == 1) {
        out.i = ng + pointIJK[0];
        out.j = normalFaceIndex;
        out.k = ng + pointIJK[2];
    } else {
        out.i = ng + pointIJK[0];
        out.j = ng + pointIJK[1];
        out.k = normalFaceIndex;
    }

    if (!validFluxIndex(f, out.direction, out.i, out.j, out.k)) {
        error = "computed interface flux index is outside Field storage";
        return false;
    }
    return true;
}

Point3 faceCofactorVector(const MeshBlockField& block,
                          const HaloInterfaceFace& face,
                          bool& ok) {
    ok = false;
    if (face.direction < 0 || face.direction > 2) return {};

    double metrics[4] = {0.0, 0.0, 0.0, 0.0};
    Math::faceMetrics(block.field,
                      face.i,
                      face.j,
                      face.k,
                      (Math::Dir)face.direction,
                      metrics);
    Point3 n{metrics[0], metrics[1], metrics[2]};
    const double n2 = norm2(n);
    ok = std::isfinite(n.x) && std::isfinite(n.y) &&
         std::isfinite(n.z) && n2 > 1.0e-300;
    return n;
}

bool buildPatchInterfaceTopology(
    const std::vector<MeshBlockField>& blocks,
    const BoundaryFaceCandidate& ownerFace,
    const BoundaryFaceCandidate& neighbourFace,
    int requiredHaloWidth,
    PatchInterfaceTopology& topology,
    std::string& error) {
    const MeshBlockField& ownerBlock = blocks[(size_t)ownerFace.blockId];
    const MeshBlockField& neighbourBlock =
        blocks[(size_t)neighbourFace.blockId];
    int ownerAxis = -1;
    int ownerSide = -1;
    int neighbourAxis = -1;
    int neighbourSide = -1;
    if (!boundarySide(ownerBlock, ownerFace,
                      ownerAxis, ownerSide, error) ||
        !boundarySide(neighbourBlock, neighbourFace,
                      neighbourAxis, neighbourSide, error)) {
        return false;
    }

    std::array<std::array<int, 3>, 4> ownerCorners{};
    std::array<std::array<int, 3>, 4> neighbourCorners{};
    for (int n = 0; n < 4; ++n) {
        const int neighbourCorner = matchingCorner(
            neighbourFace, ownerFace.globalPointIds[(size_t)n]);
        if (neighbourCorner < 0 ||
            !localPointIJK(ownerBlock.field,
                           ownerFace.pointIds[(size_t)n],
                           ownerCorners[(size_t)n]) ||
            !localPointIJK(neighbourBlock.field,
                           neighbourFace.pointIds[(size_t)neighbourCorner],
                           neighbourCorners[(size_t)n])) {
            error = "failed to pair conformal interface corners";
            return false;
        }
    }

    topology = {};
    topology.ownerBlock = ownerFace.blockId;
    topology.neighbourBlock = neighbourFace.blockId;
    topology.kind = PatchBoundaryKind::CoupledConformal;
    topology.relation = InterfaceRelation::FaceToFace;
    topology.ownerSide.axis = ownerAxis;
    topology.ownerSide.sign = ownerSide == 0 ? -1 : 1;
    topology.neighbourSide.axis = neighbourAxis;
    topology.neighbourSide.sign = neighbourSide == 0 ? -1 : 1;
    topology.haloWidth = requiredHaloWidth;
    topology.canonicalGeometryOwner = ownerFace.blockId;

    topology.ownerSide.begin = ownerCorners[0];
    topology.ownerSide.end = ownerCorners[0];
    topology.neighbourSide.begin = neighbourCorners[0];
    topology.neighbourSide.end = neighbourCorners[0];
    for (int n = 1; n < 4; ++n) {
        for (int axis = 0; axis < 3; ++axis) {
            topology.ownerSide.begin[(size_t)axis] = std::min(
                topology.ownerSide.begin[(size_t)axis],
                ownerCorners[(size_t)n][(size_t)axis]);
            topology.ownerSide.end[(size_t)axis] = std::max(
                topology.ownerSide.end[(size_t)axis],
                ownerCorners[(size_t)n][(size_t)axis]);
            topology.neighbourSide.begin[(size_t)axis] = std::min(
                topology.neighbourSide.begin[(size_t)axis],
                neighbourCorners[(size_t)n][(size_t)axis]);
            topology.neighbourSide.end[(size_t)axis] = std::max(
                topology.neighbourSide.end[(size_t)axis],
                neighbourCorners[(size_t)n][(size_t)axis]);
        }
    }

    auto& transform = topology.indexTransform;
    transform.neighbourAxisForOwner.fill(-1);
    transform.orientation.fill(0);
    transform.neighbourAxisForOwner[(size_t)ownerAxis] = neighbourAxis;
    transform.orientation[(size_t)ownerAxis] =
        ownerSide == neighbourSide ? -1 : 1;

    std::array<bool, 3> neighbourAxisUsed{false, false, false};
    neighbourAxisUsed[(size_t)neighbourAxis] = true;
    for (int ownerTangential = 0; ownerTangential < 3;
         ++ownerTangential) {
        if (ownerTangential == ownerAxis) continue;
        bool found = false;
        for (int a = 0; a < 4 && !found; ++a) {
            for (int b = a + 1; b < 4 && !found; ++b) {
                int ownerChangedAxes = 0;
                for (int axis = 0; axis < 3; ++axis) {
                    ownerChangedAxes +=
                        ownerCorners[(size_t)a][(size_t)axis] !=
                        ownerCorners[(size_t)b][(size_t)axis];
                }
                if (ownerChangedAxes != 1 ||
                    ownerCorners[(size_t)a][(size_t)ownerTangential] ==
                    ownerCorners[(size_t)b][(size_t)ownerTangential]) {
                    continue;
                }

                int mappedAxis = -1;
                for (int axis = 0; axis < 3; ++axis) {
                    if (neighbourCorners[(size_t)a][(size_t)axis] ==
                        neighbourCorners[(size_t)b][(size_t)axis]) {
                        continue;
                    }
                    if (mappedAxis >= 0) {
                        error = "conformal interface edge is not aligned with one neighbour index axis";
                        return false;
                    }
                    mappedAxis = axis;
                }
                if (mappedAxis < 0 || neighbourAxisUsed[(size_t)mappedAxis]) {
                    error = "conformal interface index transform is singular";
                    return false;
                }
                const int ownerDelta =
                    ownerCorners[(size_t)b][(size_t)ownerTangential] -
                    ownerCorners[(size_t)a][(size_t)ownerTangential];
                const int neighbourDelta =
                    neighbourCorners[(size_t)b][(size_t)mappedAxis] -
                    neighbourCorners[(size_t)a][(size_t)mappedAxis];
                transform.neighbourAxisForOwner[(size_t)ownerTangential] =
                    mappedAxis;
                transform.orientation[(size_t)ownerTangential] =
                    ownerDelta * neighbourDelta > 0 ? 1 : -1;
                neighbourAxisUsed[(size_t)mappedAxis] = true;
                found = true;
            }
        }

        if (!found) {
            // 单层 inactive 方向没有边长；以剩余轴完成置换，仍保留显式映射。
            for (int axis = 0; axis < 3; ++axis) {
                if (neighbourAxisUsed[(size_t)axis]) continue;
                transform.neighbourAxisForOwner[(size_t)ownerTangential] = axis;
                transform.orientation[(size_t)ownerTangential] = 1;
                neighbourAxisUsed[(size_t)axis] = true;
                found = true;
                break;
            }
        }
        if (!found) {
            error = "conformal interface index transform is incomplete";
            return false;
        }
    }

    HaloInterfaceFace canonicalFace;
    if (!boundaryFluxIndexForPoint(ownerBlock,
                                   ownerFace.blockId,
                                   ownerFace.faceId,
                                   ownerFace.pointIds[0],
                                   canonicalFace,
                                   error)) {
        return false;
    }
    bool metricOk = false;
    const Point3 canonical =
        faceCofactorVector(ownerBlock, canonicalFace, metricOk);
    if (!metricOk) {
        error = "canonical interface face geometry is invalid";
        return false;
    }
    topology.canonicalCofactor = {canonical.x, canonical.y, canonical.z};
    return true;
}

GlobalFaceKey makeGlobalFaceKey(std::array<int, 4> ids) {
    std::sort(ids.begin(), ids.end());
    return {ids};
}

bool boundarySide(const MeshBlockField& block,
                  const BoundaryFaceCandidate& candidate,
                  int& axis,
                  int& side,
                  std::string& error) {
    axis = -1;
    side = -1;
    const Field& field = block.field;
    const CellFaceMesh& topology = block.topology;
    if (candidate.faceId < 0 ||
        candidate.faceId >= (int)topology.faces.size()) {
        error = "invalid conformal interface face id";
        return false;
    }
    const MeshFace& face = topology.faces[(size_t)candidate.faceId];
    axis = face.direction;
    if (axis < 0 || axis > 2) {
        error = "conformal interface face has invalid direction";
        return false;
    }
    const int axisN = structuredAxisSize(field, axis);
    if (axisN <= 1) {
        error = "conformal interface axis has fewer than two points";
        return false;
    }

    int normalCoord = -1;
    for (int pointId : candidate.pointIds) {
        std::array<int, 3> ijk{};
        if (!localPointIJK(field, pointId, ijk)) {
            error = "conformal interface face references invalid point";
            return false;
        }
        if (normalCoord < 0) {
            normalCoord = ijk[(size_t)axis];
        } else if (normalCoord != ijk[(size_t)axis]) {
            error = "conformal interface face is not index-aligned";
            return false;
        }
    }

    if (normalCoord == 0) {
        side = 0;
        return true;
    }
    if (normalCoord == axisN - 1) {
        side = 1;
        return true;
    }
    error = "conformal interface face is not on a patch boundary";
    return false;
}

int matchingCorner(const BoundaryFaceCandidate& face,
                   int globalPointId) {
    for (int n = 0; n < 4; ++n) {
        if (face.globalPointIds[(size_t)n] == globalPointId) return n;
    }
    return -1;
}

bool addConformalGhostLayerMappings(
    const std::vector<MeshBlockField>& blocks,
    const BoundaryFaceCandidate& ownerFace,
    const BoundaryFaceCandidate& donorFace,
    int requiredHaloWidth,
    ConformalGhostMap& mappings,
    std::string& error) {
    const MeshBlockField& ownerBlock = blocks[(size_t)ownerFace.blockId];
    const MeshBlockField& donorBlock = blocks[(size_t)donorFace.blockId];
    const Field& owner = ownerBlock.field;
    const Field& donor = donorBlock.field;

    int ownerAxis = -1;
    int ownerSide = -1;
    int donorAxis = -1;
    int donorSide = -1;
    if (!boundarySide(ownerBlock, ownerFace, ownerAxis, ownerSide, error) ||
        !boundarySide(donorBlock, donorFace, donorAxis, donorSide, error)) {
        return false;
    }

    const int ng = owner.NG();
    if (donor.NG() != ng) {
        error = "conformal interface blocks have different ghost depths";
        return false;
    }
    const int ownerAxisN = structuredAxisSize(owner, ownerAxis);
    const int donorAxisN = structuredAxisSize(donor, donorAxis);
    if (requiredHaloWidth <= 0 || requiredHaloWidth > ng) {
        error = "numerics-required halo width is outside allocated ghost depth";
        return false;
    }
    if (ownerAxisN <= requiredHaloWidth ||
        donorAxisN <= requiredHaloWidth) {
        error = "conformal interface block is too thin for requested ghost depth";
        return false;
    }

    for (int n = 0; n < 4; ++n) {
        const int donorCorner =
            matchingCorner(donorFace,
                           ownerFace.globalPointIds[(size_t)n]);
        if (donorCorner < 0) {
            error = "conformal interface faces do not share all corner ids";
            return false;
        }

        std::array<int, 3> ownerPoint{};
        std::array<int, 3> donorPoint{};
        if (!localPointIJK(owner, ownerFace.pointIds[(size_t)n], ownerPoint) ||
            !localPointIJK(donor,
                           donorFace.pointIds[(size_t)donorCorner],
                           donorPoint)) {
            error = "conformal interface corner point is invalid";
            return false;
        }

        for (int layer = 1; layer <= requiredHaloWidth; ++layer) {
            std::array<int, 3> ownerIJK{
                ng + ownerPoint[0],
                ng + ownerPoint[1],
                ng + ownerPoint[2]
            };
            ownerIJK[(size_t)ownerAxis] =
                ownerSide == 0
                    ? ng - layer
                    : ng + ownerAxisN - 1 + layer;
            if (!isGhostPoint(owner, ownerIJK[0], ownerIJK[1],
                              ownerIJK[2])) {
                error = "conformal owner layer did not produce a ghost point";
                return false;
            }

            std::array<int, 3> donorIJK{
                ng + donorPoint[0],
                ng + donorPoint[1],
                ng + donorPoint[2]
            };
            donorIJK[(size_t)donorAxis] =
                donorSide == 0
                    ? ng + layer
                    : ng + donorAxisN - 1 - layer;
            if (isGhostPoint(donor, donorIJK[0], donorIJK[1],
                             donorIJK[2])) {
                error = "conformal donor layer is outside donor interior";
                return false;
            }

            HaloCellMapping mapping;
            mapping.kind = HaloMappingKind::DirectCopy;
            mapping.ownerBlock = ownerFace.blockId;
            mapping.ownerIJK = ownerIJK;
            mapping.ownerIndex =
                owner.getIdx(ownerIJK[0], ownerIJK[1], ownerIJK[2]);
            mapping.donorBlock = donorFace.blockId;
            mapping.donorIndex =
                donor.getIdx(donorIJK[0], donorIJK[1], donorIJK[2]);
            mapping.donorCellIJK = donorIJK;

            const GhostKey key{mapping.ownerBlock, mapping.ownerIndex};
            auto inserted = mappings.emplace(key, mapping);
            if (!inserted.second) {
                const HaloCellMapping& old = inserted.first->second;
                if (old.donorBlock != mapping.donorBlock ||
                    old.donorIndex != mapping.donorIndex) {
                    error = "conformal interface gives ambiguous ghost donor";
                    return false;
                }
            }
        }
    }
    return true;
}

bool buildConformalGhostMap(const std::vector<MeshBlockField>& blocks,
                            ConformalGhostMap& mappings,
                            HaloExchangePlan& plan,
                            int requiredHaloWidth) {
    mappings.clear();
    std::map<GlobalFaceKey, std::vector<BoundaryFaceCandidate>> faceGroups;

    for (size_t blockId = 0; blockId < blocks.size(); ++blockId) {
        const MeshBlockField& block = blocks[blockId];
        const Field& field = block.field;
        const size_t expected =
            (size_t)field.NX() * field.NY() * field.NZ();
        if (block.topology.empty() || block.globalPointIds.size() != expected) {
            continue;
        }

        for (int faceId = 0; faceId < (int)block.topology.faces.size();
             ++faceId) {
            const MeshFace& face = block.topology.faces[(size_t)faceId];
            if (face.neighbourCell >= 0) continue;
            if (face.direction < 0 || face.direction > 2) continue;
            if (!Math::isDirectionActive((Math::Dir)face.direction)) continue;

            BoundaryFaceCandidate candidate;
            candidate.blockId = (int)blockId;
            candidate.faceId = faceId;
            candidate.pointIds = face.pointIds;
            bool valid = true;
            for (int n = 0; n < 4; ++n) {
                const int pointId = face.pointIds[(size_t)n];
                if (pointId < 0 || pointId >= (int)expected) {
                    valid = false;
                    break;
                }
                candidate.globalPointIds[(size_t)n] =
                    block.globalPointIds[(size_t)pointId];
                if (candidate.globalPointIds[(size_t)n] < 0) {
                    valid = false;
                    break;
                }
            }
            if (!valid) continue;
            faceGroups[makeGlobalFaceKey(candidate.globalPointIds)]
                .push_back(candidate);
        }
    }

    int facePairs = 0;
    for (const auto& [faceKey, candidates] : faceGroups) {
        (void)faceKey;
        if (candidates.size() < 2) continue;

        std::vector<BoundaryFaceCandidate> unique = candidates;
        std::sort(unique.begin(), unique.end(),
                  [](const BoundaryFaceCandidate& a,
                     const BoundaryFaceCandidate& b) {
                      if (a.blockId != b.blockId) return a.blockId < b.blockId;
                      return a.faceId < b.faceId;
                  });
        unique.erase(
            std::unique(unique.begin(), unique.end(),
                        [](const BoundaryFaceCandidate& a,
                           const BoundaryFaceCandidate& b) {
                            return a.blockId == b.blockId &&
                                   a.faceId == b.faceId;
                        }),
            unique.end());

        bool spansBlocks = false;
        for (size_t n = 1; n < unique.size(); ++n) {
            if (unique[n].blockId != unique.front().blockId) {
                spansBlocks = true;
                break;
            }
        }
        if (!spansBlocks) continue;
        if (unique.size() != 2) {
            broadcast("Fatal: ",
                      "conformal halo requires exactly two faces per shared "
                      "interface face; found "
                      + std::to_string(unique.size()));
            return false;
        }

        std::string error;
        if (!addConformalGhostLayerMappings(
                blocks, unique[0], unique[1], requiredHaloWidth,
                mappings, error) ||
            !addConformalGhostLayerMappings(
                blocks, unique[1], unique[0], requiredHaloWidth,
                mappings, error)) {
            broadcast("Fatal: ",
                      "failed to build conformal interface ghost mapping: "
                      + error);
            return false;
        }
        PatchInterfaceTopology topology;
        if (!buildPatchInterfaceTopology(
                blocks, unique[0], unique[1], requiredHaloWidth,
                topology, error)) {
            broadcast("Fatal: ",
                      "failed to build structured interface topology: "
                      + error);
            return false;
        }
        plan.interfaces.push_back(std::move(topology));
        ++facePairs;
    }

    broadcast("haloExchange conformal ghost map: ",
              std::to_string(facePairs) + " face pairs, "
              + std::to_string(mappings.size()) + " direct ghost points");
    return true;
}

bool looksLikeHaloSet(const std::string& name) {
    std::string u = name;
    std::transform(u.begin(), u.end(), u.begin(),
                   [](unsigned char c) { return (char)std::toupper(c); });
    return u.find("HALO") != std::string::npos ||
           u.find("OVERLAP") != std::string::npos ||
           u.find("INTERFACE") != std::string::npos ||
           u.find("PROCESSOR") != std::string::npos ||
           u.rfind("PROC", 0) == 0;
}

bool isPhysicalPoint(const Field& f, int i, int j, int k) {
    int ng = f.NG();
    return i >= ng && i < ng + f.NX() &&
           j >= ng && j < ng + f.NY() &&
           k >= ng && k < ng + f.NZ();
}

bool isGhostAlongAxis(const Field& f, int i, int j, int k, int axis) {
    int ng = f.NG();
    if (axis == 0) return i < ng || i >= f.NX() + ng;
    if (axis == 1) return j < ng || j >= f.NY() + ng;
    if (axis == 2) return k < ng || k >= f.NZ() + ng;
    return false;
}

int axisSize(const Field& f, int axis) {
    if (axis == 0) return f.NX();
    if (axis == 1) return f.NY();
    return f.NZ();
}

int coordinateForAxis(int i, int j, int k, int axis) {
    if (axis == 0) return i;
    if (axis == 1) return j;
    return k;
}

int configuredBoundaryAxis(const Field& f,
                           const std::vector<int>& indices) {
    std::vector<int> candidates;
    int ng = f.NG();
    for (int axis = 0; axis < 3; ++axis) {
        bool allOnBoundary = true;
        int nPhysical = 0;
        int lo = ng;
        int hi = ng + axisSize(f, axis) - 1;
        for (int idx : indices) {
            int i, j, k;
            f.getIJK(idx, i, j, k);
            if (!isPhysicalPoint(f, i, j, k)) continue;
            ++nPhysical;
            int c = coordinateForAxis(i, j, k, axis);
            if (c != lo && c != hi) {
                allOnBoundary = false;
                break;
            }
        }
        if (nPhysical > 0 && allOnBoundary) candidates.push_back(axis);
    }
    return candidates.size() == 1 ? candidates.front() : -1;
}

std::unordered_set<int> configuredPhysicalBoundaryGhosts(
        const Field& f,
        const std::vector<std::string>& names) {
    std::unordered_set<int> ids;
    const auto& sets = f.getAllSets();
    for (const std::string& name : names) {
        auto it = sets.find(name);
        if (it == sets.end()) continue;
        int axis = configuredBoundaryAxis(f, it->second);
        if (axis < 0) continue;

        for (int idx : it->second) {
            int i, j, k;
            f.getIJK(idx, i, j, k);
            if (isGhostAlongAxis(f, i, j, k, axis)) {
                ids.insert(idx);
            }
        }
    }
    return ids;
}

std::unordered_set<int> declaredHaloGhosts(
        const Field& f,
        const std::vector<std::string>& physicalBoundaryNames,
        bool& hasDeclared) {
    std::unordered_set<int> ids;
    hasDeclared = false;
    for (const auto& kv : f.getAllSets()) {
        if (!looksLikeHaloSet(kv.first)) continue;
        hasDeclared = true;
        for (int idx : kv.second) ids.insert(idx);
    }

    const std::unordered_set<int> physicalGhosts =
        configuredPhysicalBoundaryGhosts(f, physicalBoundaryNames);
    for (int idx : physicalGhosts) ids.erase(idx);
    return ids;
}

DonorCell makeCell(const Field& f, int i, int j, int k) {
    DonorCell cell;
    cell.ijk = {i, j, k};

    std::array<std::array<int, 3>, 8> ijk = {{
        {i,     j,     k},
        {i + 1, j,     k},
        {i,     j + 1, k},
        {i + 1, j + 1, k},
        {i,     j,     k + 1},
        {i + 1, j,     k + 1},
        {i,     j + 1, k + 1},
        {i + 1, j + 1, k + 1}
    }};

    cell.min = {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };
    cell.max = {
        -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max()
    };

    for (int n = 0; n < 8; ++n) {
        int ci = ijk[n][0];
        int cj = ijk[n][1];
        int ck = ijk[n][2];
        cell.indices[n] = f.getIdx(ci, cj, ck);
        cell.points[n] = pointAt(f, ci, cj, ck);
        includePoint(cell.min, cell.max, cell.points[n]);
    }
    return cell;
}

DonorCache buildDonorCache(const MeshBlockField& block, int blockId,
                           double tol) {
    DonorCache cache;
    cache.blockId = blockId;
    cache.field = &block.field;

    const Field& f = block.field;
    cache.min = {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };
    cache.max = {
        -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max()
    };

    int ng = f.NG();
    for (int k = ng; k < f.NZ() + ng; ++k) {
        for (int j = ng; j < f.NY() + ng; ++j) {
            for (int i = ng; i < f.NX() + ng; ++i) {
                includePoint(cache.min, cache.max, pointAt(f, i, j, k));
            }
        }
    }

    if (f.NX() < 2 || f.NY() < 2 || f.NZ() < 2) {
        return cache;
    }

    int nCells = (f.NX() - 1) * (f.NY() - 1) * (f.NZ() - 1);
    int base = clampInt((int)std::round(std::cbrt((double)std::max(nCells, 1))), 4, 64);
    double ex = cache.max.x - cache.min.x;
    double ey = cache.max.y - cache.min.y;
    double ez = cache.max.z - cache.min.z;
    double maxExtent = std::max(ex, std::max(ey, ez));
    cache.binsX = chooseBins(ex, maxExtent, base, tol);
    cache.binsY = chooseBins(ey, maxExtent, base, tol);
    cache.binsZ = chooseBins(ez, maxExtent, base, tol);
    cache.bins.assign((size_t)cache.binsX * cache.binsY * cache.binsZ, {});
    cache.cells.reserve((size_t)nCells);

    for (int k = ng; k < f.NZ() + ng - 1; ++k) {
        for (int j = ng; j < f.NY() + ng - 1; ++j) {
            for (int i = ng; i < f.NX() + ng - 1; ++i) {
                int cellId = (int)cache.cells.size();
                cache.cells.push_back(makeCell(f, i, j, k));
                const DonorCell& cell = cache.cells.back();

                int bx0 = binCoord(cell.min.x - tol, cache.min.x, cache.max.x, cache.binsX);
                int bx1 = binCoord(cell.max.x + tol, cache.min.x, cache.max.x, cache.binsX);
                int by0 = binCoord(cell.min.y - tol, cache.min.y, cache.max.y, cache.binsY);
                int by1 = binCoord(cell.max.y + tol, cache.min.y, cache.max.y, cache.binsY);
                int bz0 = binCoord(cell.min.z - tol, cache.min.z, cache.max.z, cache.binsZ);
                int bz1 = binCoord(cell.max.z + tol, cache.min.z, cache.max.z, cache.binsZ);

                for (int bz = bz0; bz <= bz1; ++bz) {
                    for (int by = by0; by <= by1; ++by) {
                        for (int bx = bx0; bx <= bx1; ++bx) {
                            cache.bins[binIndex(cache, bx, by, bz)].push_back(cellId);
                        }
                    }
                }
            }
        }
    }

    return cache;
}

bool findInDonor(const DonorCache& cache, const Point3& p,
                 double tol, LocatedCell& located) {
    if (cache.cells.empty()) return false;

    double blockTol = std::max(tol, 1.0e-12 * std::sqrt(distance2(cache.min, cache.max)));
    if (!insideBox(cache.min, cache.max, p, blockTol)) return false;

    int bx = binCoord(p.x, cache.min.x, cache.max.x, cache.binsX);
    int by = binCoord(p.y, cache.min.y, cache.max.y, cache.binsY);
    int bz = binCoord(p.z, cache.min.z, cache.max.z, cache.binsZ);

    bool found = false;
    LocatedCell best;
    for (int dz = -1; dz <= 1; ++dz) {
        int qz = bz + dz;
        if (qz < 0 || qz >= cache.binsZ) continue;
        for (int dy = -1; dy <= 1; ++dy) {
            int qy = by + dy;
            if (qy < 0 || qy >= cache.binsY) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                int qx = bx + dx;
                if (qx < 0 || qx >= cache.binsX) continue;
                const auto& bucket = cache.bins[binIndex(cache, qx, qy, qz)];
                for (int cellId : bucket) {
                    const DonorCell& cell = cache.cells[cellId];
                    double cellTol = std::max(tol, 1.0e-12 * cellDiag(cell));
                    if (!insideBox(cell.min, cell.max, p, cellTol)) continue;

                    std::array<double, 3> local{0.0, 0.0, 0.0};
                    std::array<double, 8> weights{};
                    double residual = std::numeric_limits<double>::max();
                    if (!invertTrilinear(cell, p, cellTol, local, weights, residual)) {
                        continue;
                    }

                    if (residual < best.residual) {
                        found = true;
                        best.donorBlock = cache.blockId;
                        best.donorCellIJK = cell.ijk;
                        best.localCoord = local;
                        best.residual = residual;
                    }
                }
            }
        }
    }

    if (!found) return false;
    located = best;
    return true;
}

bool locateGhostPoint(const std::vector<DonorCache>& caches,
                      int ownerBlock, const Point3& p,
                      double tol, LocatedCell& located) {
    bool found = false;
    LocatedCell best;
    for (const DonorCache& cache : caches) {
        if (cache.blockId == ownerBlock) continue;
        LocatedCell trial;
        if (!findInDonor(cache, p, tol, trial)) continue;
        if (trial.residual < best.residual) {
            found = true;
            best = trial;
        }
    }

    if (!found) return false;
    located = best;
    return true;
}

} // namespace

/// @brief 将donor field的全域平铺索引(含ghost padding)转换为interior-only索引。
static int interiorIndexFromFull(const Field& donorField, int fullIndex) {
    int mx = donorField.MX();
    int my = donorField.MY();
    int i = fullIndex % mx;
    int j = (fullIndex / mx) % my;
    int k = fullIndex / (mx * my);
    int ng = donorField.NG();
    int nx = donorField.NX();
    int ny = donorField.NY();
    return (k - ng) * ny * nx + (j - ng) * nx + (i - ng);
}

/// @brief 为所有映射补充 interior donor index 并收集 donor block info。
static void finalizePlanInteriorIndices(
    HaloExchangePlan& plan,
    const std::vector<MeshBlockField>& blocks) {
    plan.blockInfos.clear();
    plan.blockInfos.reserve(blocks.size());
    for (size_t blockId = 0; blockId < blocks.size(); ++blockId) {
        const Field& f = blocks[blockId].field;
        DonorBlockInfo info;
        info.blockId = (int)blockId;
        info.nx = f.NX();
        info.ny = f.NY();
        info.nz = f.NZ();
        info.ng = f.NG();
        info.interiorPointCount = info.nx * info.ny * info.nz;
        info.variableCount = f.NVar();
        info.interiorDataCount = info.interiorPointCount * info.variableCount;
        plan.blockInfos.push_back(info);
    }

    for (size_t b = 0; b < plan.blockPlans.size(); ++b) {
        HaloBlockPlan& blockPlan = plan.blockPlans[b];
        // 收集本块所有donor block的info
        std::unordered_map<int, const Field*> donorFieldMap;
        for (const HaloCellMapping& mapping : blockPlan.cells) {
            if (donorFieldMap.count(mapping.donorBlock) == 0) {
                donorFieldMap[mapping.donorBlock] =
                    &blocks[(size_t)mapping.donorBlock].field;
            }
        }
        // 加上自己的info（exchange时自己也需要pack interior）
        donorFieldMap[blockPlan.blockId] = &blocks[b].field;

        for (auto& [blkId, f] : donorFieldMap) {
            DonorBlockInfo info;
            info.blockId = blkId;
            info.nx = f->NX();
            info.ny = f->NY();
            info.nz = f->NZ();
            info.ng = f->NG();
            info.interiorPointCount = info.nx * info.ny * info.nz;
            info.variableCount = f->NVar();
            info.interiorDataCount = info.interiorPointCount * info.variableCount;
            blockPlan.donorBlockInfos.push_back(info);
        }

        // 为每个映射计算 interior donor index
        for (HaloCellMapping& mapping : blockPlan.cells) {
            const Field* donorF = donorFieldMap.count(mapping.donorBlock) > 0
                ? donorFieldMap[mapping.donorBlock] : nullptr;
            if (!donorF) continue;

            if (mapping.kind == HaloMappingKind::DirectCopy) {
                mapping.donorInteriorIndex =
                    interiorIndexFromFull(*donorF, mapping.donorIndex);
            } else {
                for (int axis = 0; axis < 3; ++axis) {
                    const int axisSize =
                        axis == 0 ? donorF->NX()
                        : (axis == 1 ? donorF->NY() : donorF->NZ());
                    const int start =
                        mapping.donorStencilStart[(size_t)axis];
                    const int size =
                        mapping.donorStencilSize[(size_t)axis];
                    if (start < 0 || size <= 0 ||
                        size > kMaxHaloInterpolationStencil ||
                        start + size > axisSize) {
                        broadcast("Fatal: ",
                                  "haloExchange tensor interpolation stencil "
                                  "is invalid during plan finalization.");
                        std::exit(1);
                    }
                }
            }
        }
    }
}

enum class GlobalPointSyncResult {
    Unavailable,
    Success,
    Invalid
};

static GlobalPointSyncResult addGlobalPhysicalPointSync(
    std::vector<MeshBlockField>& blocks,
    HaloExchangePlan& plan) {
    std::unordered_map<int, std::vector<DonorPoint>> groups;
    bool hasGlobalTopology = false;

    for (size_t blockId = 0; blockId < blocks.size(); ++blockId) {
        const MeshBlockField& block = blocks[blockId];
        const Field& field = block.field;
        const size_t expected =
            (size_t)field.NX() * field.NY() * field.NZ();
        if (block.globalPointIds.empty()) continue;
        if (block.globalPointIds.size() != expected) {
            broadcast("Fatal: ",
                      "global point topology size does not match structured patch "
                      + std::to_string(blockId));
            return GlobalPointSyncResult::Invalid;
        }
        hasGlobalTopology = true;

        const int ng = field.NG();
        size_t localId = 0;
        for (int k = 0; k < field.NZ(); ++k) {
            for (int j = 0; j < field.NY(); ++j) {
                for (int i = 0; i < field.NX(); ++i, ++localId) {
                    const int globalId = block.globalPointIds[localId];
                    if (globalId < 0) continue;

                    DonorPoint point;
                    point.blockId = (int)blockId;
                    point.ijk = {i + ng, j + ng, k + ng};
                    point.index = field.getIdx(point.ijk[0],
                                               point.ijk[1],
                                               point.ijk[2]);
                    point.point = pointAt(field,
                                          point.ijk[0],
                                          point.ijk[1],
                                          point.ijk[2]);
                    groups[globalId].push_back(point);
                }
            }
        }
    }

    if (!hasGlobalTopology) return GlobalPointSyncResult::Unavailable;

    bool hasSharedGlobalDof = false;
    for (const auto& entry : groups) {
        if (entry.second.size() > 1) {
            hasSharedGlobalDof = true;
            break;
        }
    }

    // 单块/无重复 GlobalDof 时不需要对偶体积。这样二维单层网格不会仅因
    // Runtime 支持 GlobalDof 装配而被无关的三维体积检查拒绝。
    std::vector<std::vector<double>> dualVolumes;
    if (hasSharedGlobalDof
        && !buildDualVolumeFragments(blocks, dualVolumes)) {
        return GlobalPointSyncResult::Invalid;
    }

    int syncedPoints = 0;
    int sharedGroups = 0;
    int maxMultiplicity = 1;
    for (auto& [globalId, points] : groups) {
        std::sort(points.begin(), points.end(),
                  [](const DonorPoint& a, const DonorPoint& b) {
                      if (a.blockId != b.blockId) return a.blockId < b.blockId;
                      return a.index < b.index;
                  });
        points.erase(
            std::unique(points.begin(), points.end(),
                        [](const DonorPoint& a, const DonorPoint& b) {
                            return a.blockId == b.blockId && a.index == b.index;
                        }),
            points.end());

        if (points.empty()) continue;
        const DonorPoint& ownerPoint = points.front();
        if (ownerPoint.blockId < 0 ||
            ownerPoint.blockId >= (int)blocks.size()) {
            broadcast("Fatal: ",
                      "GlobalDof canonical owner references an invalid block.");
            return GlobalPointSyncResult::Invalid;
        }
        const int ownerRank =
            blocks[(size_t)ownerPoint.blockId].ownerRank;
        for (const DonorPoint& point : points) {
            if (point.blockId < 0 ||
                point.blockId >= (int)blocks.size()) {
                broadcast("Fatal: ",
                          "GlobalDof replica references an invalid block.");
                return GlobalPointSyncResult::Invalid;
            }
            const double mismatch2 = distance2(point.point, ownerPoint.point);
            if (!std::isfinite(mismatch2) ||
                mismatch2 > plan.tolerance * plan.tolerance) {
                broadcast(
                    "Fatal: ",
                    "one GlobalDof maps to geometrically different points; "
                    "coordinate-based averaging/correction is forbidden "
                    "(GlobalDof=" + std::to_string(globalId) + ").");
                return GlobalPointSyncResult::Invalid;
            }
            Field& field = blocks[(size_t)point.blockId].field;
            field.setGlobalDofOwnership(
                point.ijk[0], point.ijk[1], point.ijk[2],
                globalId, ownerRank,
                point.blockId == ownerPoint.blockId &&
                point.index == ownerPoint.index);
        }

        if (points.size() < 2) continue;

        ++sharedGroups;
        maxMultiplicity = std::max(maxMultiplicity, (int)points.size());
        HaloInterfaceSyncGroup group;
        group.globalDofId = globalId;
        group.canonicalOwner = 0;
        group.ownerRank = ownerRank;
        group.points.reserve(points.size());
        for (const DonorPoint& point : points) {
            HaloInterfacePoint syncPoint;
            syncPoint.blockId = point.blockId;
            syncPoint.interiorIndex = interiorIndexFromFull(
                blocks[(size_t)point.blockId].field,
                point.index);
            if (syncPoint.interiorIndex < 0
                || syncPoint.interiorIndex >=
                    (int)dualVolumes[(size_t)point.blockId].size()) {
                broadcast(
                    "Fatal: ",
                    "GlobalDof dual-volume point index is invalid.");
                return GlobalPointSyncResult::Invalid;
            }
            syncPoint.dualVolume =
                dualVolumes[(size_t)point.blockId]
                           [(size_t)syncPoint.interiorIndex];
            if (!std::isfinite(syncPoint.dualVolume)
                || syncPoint.dualVolume <= 1.0e-300) {
                broadcast(
                    "Fatal: ",
                    "GlobalDof has a non-positive dual-volume fragment "
                    "(GlobalDof=" + std::to_string(globalId)
                    + ", block=" + std::to_string(point.blockId) + ").");
                return GlobalPointSyncResult::Invalid;
            }
            group.points.push_back(syncPoint);
        }
        syncedPoints += (int)group.points.size() - 1;
        plan.synchronizedPhysicalPoints +=
            (int)group.points.size() - 1;
        plan.interfaceSyncGroups.push_back(std::move(group));
    }

    std::ostringstream msg;
    msg << sharedGroups << " global point groups, "
        << syncedPoints << " duplicate references, "
        << "max multiplicity=" << maxMultiplicity;
    broadcast("haloExchange global point sync: ", msg.str());
    return GlobalPointSyncResult::Success;
}

enum class GlobalFluxSyncResult {
    Success,
    Invalid
};

struct GlobalDofFaceKey {
    int sourceZoneId = -1;
    int sourceDirection = -1;
    int lower = -1;
    int upper = -1;

    bool operator<(const GlobalDofFaceKey& other) const {
        return std::tie(sourceZoneId, sourceDirection, lower, upper) <
               std::tie(other.sourceZoneId, other.sourceDirection,
                        other.lower, other.upper);
    }
};

int interiorGlobalDof(const MeshBlockField& block, int fullIndex) {
    const Field& field = block.field;
    int i = 0, j = 0, k = 0;
    field.getIJK(fullIndex, i, j, k);
    const int ng = field.NG();
    if (i < ng || i >= ng + field.NX() ||
        j < ng || j >= ng + field.NY() ||
        k < ng || k >= ng + field.NZ()) {
        return -1;
    }
    const int interior =
        ((k - ng) * field.NY() + (j - ng)) * field.NX()
        + (i - ng);
    if (interior < 0 ||
        interior >= (int)block.globalPointIds.size()) {
        return -1;
    }
    return block.globalPointIds[(size_t)interior];
}

int faceEndpointGlobalDof(
        const std::vector<MeshBlockField>& blocks,
        const ConformalGhostMap& conformalGhosts,
        int blockId, int fullIndex) {
    if (blockId < 0 || blockId >= (int)blocks.size()) return -1;
    const int local = interiorGlobalDof(blocks[(size_t)blockId], fullIndex);
    if (local >= 0) return local;

    const auto mapping = conformalGhosts.find({blockId, fullIndex});
    if (mapping == conformalGhosts.end() ||
        mapping->second.kind != HaloMappingKind::DirectCopy) {
        return -1;
    }
    const HaloCellMapping& donor = mapping->second;
    if (donor.donorBlock < 0 ||
        donor.donorBlock >= (int)blocks.size()) {
        return -1;
    }
    return interiorGlobalDof(
        blocks[(size_t)donor.donorBlock], donor.donorIndex);
}

GlobalFluxSyncResult addGlobalDofFaceFluxSync(
        std::vector<MeshBlockField>& blocks,
        const ConformalGhostMap& conformalGhosts,
        HaloExchangePlan& plan) {
    std::map<GlobalDofFaceKey, std::vector<HaloInterfaceFace>> groups;

    for (int blockId = 0; blockId < (int)blocks.size(); ++blockId) {
        Field& field = blocks[(size_t)blockId].field;
        for (int direction = 0; direction < 3; ++direction) {
            const Math::Dir dir = (Math::Dir)direction;
            if (!Math::isDirectionActive(dir)) continue;
            Math::forFaces(field, dir, [&](int i, int j, int k) {
                int di = 0, dj = 0, dk = 0;
                Math::dirOffset(dir, di, dj, dk);
                const int lowerIndex = field.getIdx(i, j, k);
                const int upperIndex = field.getIdx(i + di, j + dj, k + dk);
                const int localLower = faceEndpointGlobalDof(
                    blocks, conformalGhosts, blockId, lowerIndex);
                const int localUpper = faceEndpointGlobalDof(
                    blocks, conformalGhosts, blockId, upperIndex);
                if (localLower < 0 || localUpper < 0 ||
                    localLower == localUpper) {
                    return;
                }

                GlobalDofFaceKey key{
                    blocks[(size_t)blockId].sourceZoneId,
                    direction,
                    std::min(localLower, localUpper),
                    std::max(localLower, localUpper)};
                HaloInterfaceFace face;
                face.blockId = blockId;
                face.direction = direction;
                face.i = i;
                face.j = j;
                face.k = k;
                face.orientation = localLower == key.lower ? 1.0 : -1.0;
                groups[key].push_back(face);
            });
        }
    }

    int globalFaceId = 0;
    int sharedGroups = 0;
    int duplicateReferences = 0;
    int maxMultiplicity = 1;
    double maxGeometryMismatch = 0.0;
    for (auto& [key, faces] : groups) {
        (void)key;
        std::sort(faces.begin(), faces.end(),
                  [](const HaloInterfaceFace& a,
                     const HaloInterfaceFace& b) {
                      return participantKey(a) < participantKey(b);
                  });
        faces.erase(std::unique(
                        faces.begin(), faces.end(),
                        [](const HaloInterfaceFace& a,
                           const HaloInterfaceFace& b) {
                            return participantKey(a) == participantKey(b);
                        }),
                    faces.end());
        const int thisGlobalFaceId = globalFaceId++;
        if (faces.size() < 2) continue;

        int canonicalOwner = -1;
        for (int n = 0; n < (int)faces.size(); ++n) {
            const HaloInterfaceFace& candidate = faces[(size_t)n];
            const Field& field =
                blocks[(size_t)candidate.blockId].field;
            int di = 0, dj = 0, dk = 0;
            Math::dirOffset((Math::Dir)candidate.direction,
                            di, dj, dk);
            if (isPhysicalPoint(field,
                                candidate.i,
                                candidate.j,
                                candidate.k) &&
                isPhysicalPoint(field,
                                candidate.i + di,
                                candidate.j + dj,
                                candidate.k + dk)) {
                canonicalOwner = n;
                break;
            }
        }
        if (canonicalOwner < 0) {
            // 原始 source-zone 边界面在另一坐标方向被 MPI 切分时，所有副本
            // 都可能是外侧半面。几何已经在切分前由 source SCMM 唯一生成，
            // 因而此处按稳定 participant 顺序选择 face owner，不重算/平均度规。
            canonicalOwner = 0;
        }
        const HaloInterfaceFace& owner =
            faces[(size_t)canonicalOwner];
        bool ownerOk = false;
        const Point3 ownerLocal = faceCofactorVector(
            blocks[(size_t)owner.blockId], owner, ownerOk);
        if (!ownerOk) {
            broadcast("Fatal: ",
                      "GlobalFace canonical metric is invalid.");
            return GlobalFluxSyncResult::Invalid;
        }
        const Point3 canonical{
            owner.orientation * ownerLocal.x,
            owner.orientation * ownerLocal.y,
            owner.orientation * ownerLocal.z};
        double ownerMetrics[4]{0.0, 0.0, 0.0, 0.0};
        Math::faceMetrics(
            blocks[(size_t)owner.blockId].field,
            owner.i, owner.j, owner.k,
            (Math::Dir)owner.direction, ownerMetrics);
        if (!std::isfinite(ownerMetrics[3]) ||
            std::abs(ownerMetrics[3]) <= 1.0e-300) {
            broadcast("Fatal: ",
                      "GlobalFace canonical inverse Jacobian is invalid.");
            return GlobalFluxSyncResult::Invalid;
        }

        for (const HaloInterfaceFace& face : faces) {
            bool localOk = false;
            const Point3 local = faceCofactorVector(
                blocks[(size_t)face.blockId], face, localOk);
            if (!localOk) {
                broadcast("Fatal: ",
                          "GlobalFace replica metric is invalid.");
                return GlobalFluxSyncResult::Invalid;
            }
            const Point3 oriented{
                face.orientation * local.x,
                face.orientation * local.y,
                face.orientation * local.z};
            const double scale = std::max(
                std::sqrt(norm2(canonical)), std::sqrt(norm2(oriented)));
            const double mismatch = std::sqrt(
                (canonical.x - oriented.x) * (canonical.x - oriented.x)
              + (canonical.y - oriented.y) * (canonical.y - oriented.y)
              + (canonical.z - oriented.z) * (canonical.z - oriented.z));
            if (!std::isfinite(scale) || scale <= 1.0e-300 ||
                !std::isfinite(mismatch) ||
                mismatch > 64.0 * plan.tolerance * scale) {
                broadcast(
                    "Fatal: ",
                    "fully matching GlobalFace replicas have incompatible "
                    "cofactor geometry (GlobalFace="
                    + std::to_string(thisGlobalFaceId)
                    + ", sourceZone=" + std::to_string(key.sourceZoneId)
                    + ", sourceDirection="
                    + std::to_string(key.sourceDirection)
                    + ", dofs=" + std::to_string(key.lower) + ":"
                    + std::to_string(key.upper)
                    + ", ownerBlock=" + std::to_string(owner.blockId)
                    + ", replicaBlock=" + std::to_string(face.blockId)
                    + ", relativeMismatch="
                    + std::to_string(mismatch / scale)
                    + ", canonical=(" + std::to_string(canonical.x)
                    + "," + std::to_string(canonical.y) + ","
                    + std::to_string(canonical.z) + "), replica=("
                    + std::to_string(oriented.x) + ","
                    + std::to_string(oriented.y) + ","
                    + std::to_string(oriented.z) + ")).");
                return GlobalFluxSyncResult::Invalid;
            }
            maxGeometryMismatch = std::max(
                maxGeometryMismatch, mismatch / scale);
        }

        for (int n = 0; n < (int)faces.size(); ++n) {
            const HaloInterfaceFace& face = faces[(size_t)n];
            Field& field = blocks[(size_t)face.blockId].field;
            field.setCanonicalFaceMetrics(
                face.direction, face.i, face.j, face.k,
                {face.orientation * canonical.x,
                 face.orientation * canonical.y,
                 face.orientation * canonical.z,
                 ownerMetrics[3]});
            field.setCanonicalFaceOwner(
                face.direction, face.i, face.j, face.k,
                n == canonicalOwner);
        }

        HaloInterfaceFluxSyncGroup group;
        group.faces = std::move(faces);
        group.canonicalOwner = canonicalOwner;
        group.globalFaceId = thisGlobalFaceId;
        group.canonicalCofactor = {
            canonical.x, canonical.y, canonical.z};
        group.canonicalInverseJacobian = ownerMetrics[3];
        maxMultiplicity = std::max(
            maxMultiplicity, (int)group.faces.size());
        duplicateReferences += (int)group.faces.size() - 1;
        plan.synchronizedInterfaceFluxes +=
            (int)group.faces.size() - 1;
        plan.interfaceFluxSyncGroups.push_back(std::move(group));
        ++sharedGroups;
    }

    std::ostringstream msg;
    msg << globalFaceId << " GlobalFaces, "
        << sharedGroups << " shared groups, "
        << duplicateReferences << " duplicate references, "
        << "max multiplicity=" << maxMultiplicity
        << ", max relative metric mismatch=" << maxGeometryMismatch;
    broadcast("GlobalFace preprocess: ", msg.str());
    return GlobalFluxSyncResult::Success;
}

bool buildHaloExchangePlan(std::vector<MeshBlockField>& blocks,
                           const HaloPreprocessOptions& options,
                           HaloExchangePlan& plan) {
    plan.clear();
    if (!std::isfinite(options.tolerance) || options.tolerance <= 0.0) {
        broadcast("Fatal: ",
                  "halo/interface tolerance must be explicitly positive and finite.");
        return false;
    }
    plan.tolerance = options.tolerance;
    plan.haloWidth = options.requiredHaloWidth;

    if (plan.haloWidth <= 0) {
        broadcast("Fatal: ",
                  "numerics-required halo width must be positive.");
        return false;
    }
    for (size_t b = 0; b < blocks.size(); ++b) {
        if (plan.haloWidth > blocks[b].field.NG()) {
            broadcast(
                "Fatal: ",
                "numerics requires " + std::to_string(plan.haloWidth)
                + " halo layers but patch " + std::to_string(b)
                + " allocated only "
                + std::to_string(blocks[b].field.NG()) + ".");
            return false;
        }
    }

    if (blocks.size() < 2) return true;

    std::vector<DonorCache> caches;
    caches.reserve(blocks.size());
    for (size_t b = 0; b < blocks.size(); ++b) {
        caches.push_back(buildDonorCache(blocks[b], (int)b, plan.tolerance));
    }
    ConformalGhostMap conformalGhosts;
    if (!buildConformalGhostMap(
            blocks, conformalGhosts, plan, plan.haloWidth)) {
        return false;
    }

    plan.blockPlans.reserve(blocks.size());
    for (size_t owner = 0; owner < blocks.size(); ++owner) {
        const Field& f = blocks[owner].field;
        HaloBlockPlan blockPlan;
        blockPlan.blockId = (int)owner;

        bool hasDeclaredHalo = false;
        std::unordered_set<int> declaredIds;
        const std::unordered_set<int> physicalGhostIds =
            configuredPhysicalBoundaryGhosts(
                f, options.physicalBoundaryNames);
        if (options.useDeclaredHaloSets) {
            declaredIds = declaredHaloGhosts(
                f, options.physicalBoundaryNames, hasDeclaredHalo);
        }

        for (int k = 0; k < f.MZ(); ++k) {
            for (int j = 0; j < f.MY(); ++j) {
                for (int i = 0; i < f.MX(); ++i) {
                    if (!isGhostPoint(f, i, j, k)) continue;
                    int ownerIndex = f.getIdx(i, j, k);
                    const int ghostLayer = activeGhostLayer(f, i, j, k);
                    if (ghostLayer == 0 || ghostLayer > plan.haloWidth) {
                        continue;
                    }
                    if (physicalGhostIds.count(ownerIndex) != 0) continue;
                    if (hasDeclaredHalo && declaredIds.count(ownerIndex) == 0) {
                        continue;
                    }

                    ++plan.searchedGhostCells;
                    const GhostKey ghostKey{(int)owner, ownerIndex};
                    auto conformal = conformalGhosts.find(ghostKey);
                    if (conformal != conformalGhosts.end()) {
                        blockPlan.cells.push_back(conformal->second);
                        ++plan.mappedGhostCells;
                        ++plan.directMappedGhostCells;
                        continue;
                    }

                    Point3 p = pointAt(f, i, j, k);
                    LocatedCell located;
                    if (locateGhostPoint(caches, (int)owner, p,
                                         plan.tolerance, located)) {
                        if (activeExteriorAxisCount(f, i, j, k) != 1) {
                            ++plan.unmappedGhostCells;
                            continue;
                        }
                        if (!options.allowNonMatchingInterfaces) {
                            broadcast(
                                "Fatal: ",
                                "non-matching/interpolated interface detected "
                                "at ownerBlock=" + std::to_string(owner)
                                + ", ownerIJK=(" + std::to_string(i) + ","
                                + std::to_string(j) + ","
                                + std::to_string(k)
                                + "). Current structured-FDM MPI path only "
                                  "accepts fully matching interfaces.");
                            return false;
                        }
                        HaloCellMapping mapping;
                        mapping.kind = HaloMappingKind::Interpolated;
                        mapping.ownerBlock = (int)owner;
                        mapping.ownerIndex = ownerIndex;
                        mapping.ownerIJK = {i, j, k};
                        mapping.donorBlock = located.donorBlock;
                        mapping.donorCellIJK = located.donorCellIJK;
                        mapping.localCoord = located.localCoord;
                        int tensorPoints = 0;
                        if (located.donorBlock < 0 ||
                            located.donorBlock >= (int)blocks.size() ||
                            !buildTensorInterpolation(
                                blocks[(size_t)located.donorBlock].field,
                                located,
                                mapping,
                                tensorPoints)) {
                            broadcast(
                                "Fatal: ",
                                "failed to build high-order non-matching halo "
                                "interpolation at ownerBlock="
                                + std::to_string(owner) + ", ownerIJK=("
                                + std::to_string(i) + ","
                                + std::to_string(j) + ","
                                + std::to_string(k) + ").");
                            return false;
                        }
                        blockPlan.cells.push_back(mapping);
                        ++plan.mappedGhostCells;
                        ++plan.interpolatedMappedGhostCells;
                        if (plan.minInterpolationTensorPoints == 0) {
                            plan.minInterpolationTensorPoints = tensorPoints;
                        } else {
                            plan.minInterpolationTensorPoints = std::min(
                                plan.minInterpolationTensorPoints,
                                tensorPoints);
                        }
                        plan.maxInterpolationTensorPoints = std::max(
                            plan.maxInterpolationTensorPoints,
                            tensorPoints);
                        continue;
                    } else {
                        ++plan.unmappedGhostCells;
                        continue;
                    }
                }
            }
        }

        plan.blockPlans.push_back(std::move(blockPlan));
    }

    const GlobalPointSyncResult globalSync =
        addGlobalPhysicalPointSync(blocks, plan);
    if (globalSync == GlobalPointSyncResult::Invalid) return false;
    if (globalSync == GlobalPointSyncResult::Unavailable) {
        broadcast(
            "Fatal: ",
            "GlobalDof ownership requires globalPointIds on every patch; "
            "coordinate-coincidence grouping and averaging are forbidden.");
        return false;
    }

    if (addGlobalDofFaceFluxSync(blocks, conformalGhosts, plan)
        == GlobalFluxSyncResult::Invalid) {
        return false;
    }

    std::ostringstream msg;
    msg << plan.mappedGhostCells << " mapped / "
        << plan.searchedGhostCells
        << " ghost points, direct="
        << plan.directMappedGhostCells
        << ", interpolated="
        << plan.interpolatedMappedGhostCells
        << ", tensorPoints="
        << plan.minInterpolationTensorPoints
        << ".."
        << plan.maxInterpolationTensorPoints
        << ", physicalSync=" << plan.synchronizedPhysicalPoints
        << ", fluxSync=" << plan.synchronizedInterfaceFluxes
        << ", unmappedGhost=" << plan.unmappedGhostCells
        << ", haloWidth=" << plan.haloWidth
        << ", tolerance="
        << plan.tolerance;
    broadcast("haloExchange preprocess: ", msg.str());

    finalizePlanInteriorIndices(plan, blocks);
    return true;
}

} // namespace MeshCommunication
} // namespace SF
