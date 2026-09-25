/// @file SF_STLGeometry.h
/// @brief STL 读取、三角形查询、内外分类和最近点几何工具。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "core/interfaces/SF_log.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace SF {
namespace IBM {
namespace GeoProcessing {

struct Point {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Point operator+(const Point& a, const Point& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Point operator-(const Point& a, const Point& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Point operator*(const Point& p, double s) { return {p.x * s, p.y * s, p.z * s}; }
inline Point operator*(double s, const Point& p) { return p * s; }
inline Point operator/(const Point& p, double s) { return {p.x / s, p.y / s, p.z / s}; }
inline double dot(const Point& a, const Point& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Point cross(const Point& a, const Point& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm2(const Point& p) { return dot(p, p); }
inline double norm(const Point& p) { return std::sqrt(norm2(p)); }
inline Point normalize(const Point& p) {
    double n = norm(p);
    return (n > 1e-14) ? p / n : Point{0.0, 0.0, 0.0};
}

struct Triangle {
    Point normal;
    Point v[3];
};

struct SDFResult {
    double signedDistance = std::numeric_limits<double>::max();
    double distance = std::numeric_limits<double>::max();
    bool inside = false;
    Point closestPoint;
    Point normal;
    std::size_t triangleIndex = 0;
};

struct WatertightReport {
    bool watertight = false;
    int boundaryEdges = 0;
    int nonManifoldEdges = 0;
    int duplicateEdges = 0;
};

class STLGeometry {
public:
    bool loadSTLFiles(const std::vector<std::string>& files, const std::string& caseDir) {
        triangles_.clear();
        bool ok = true;
        for (const auto& file : files) {
            std::string path = resolvePath(file, caseDir);
            bool loaded = loadASCIISTL(path);
            if (!loaded) loaded = loadBinarySTL(path);
            ok = ok && loaded;
            if (!loaded) broadcast("IBM warning: failed to read STL ", path);
        }
        finalize();
        return ok && !triangles_.empty();
    }

    const std::vector<Triangle>& triangles() const { return triangles_; }
    const WatertightReport& watertightReport() const { return watertight_; }
    bool watertight() const { return watertight_.watertight; }
    bool normalsOutward() const { return normalsOutward_; }

    SDFResult signedDistance(const Point& p) const {
        SDFResult out;
        if (triangles_.empty()) return out;

        nearestTriangle(0, p, out);

        out.inside = contains(p);
        out.signedDistance = out.inside ? -out.distance : out.distance;
        if (norm2(out.normal) < 1e-24) out.normal = normalize(p - out.closestPoint);
        if (out.inside && dot(p - out.closestPoint, out.normal) > 0.0) out.normal = -1.0 * out.normal;
        if (!out.inside && dot(p - out.closestPoint, out.normal) < 0.0) out.normal = -1.0 * out.normal;
        return out;
    }

    bool contains(const Point& p) const {
        if (triangles_.empty() || !insideBounds(p)) return false;

        static constexpr Point dirs[3] = {
            {1.0, 0.1732050808, 0.0977197538},
            {0.131557, 1.0, 0.213791},
            {0.271828, 0.314159, 1.0}
        };

        int votes = 0;
        for (const Point& dir : dirs) {
            int hits = countRayHits(0, p, dir);
            if ((hits % 2) == 1) ++votes;
        }
        return votes >= 2;
    }

private:
    std::vector<Triangle> triangles_;
    std::vector<int> triIndices_;
    struct BVHNode {
        Point min;
        Point max;
        int left = -1;
        int right = -1;
        int start = 0;
        int count = 0;
    };
    std::vector<BVHNode> nodes_;
    Point min_, max_;
    WatertightReport watertight_;
    bool normalsOutward_ = true;

    static std::string resolvePath(const std::string& file, const std::string& caseDir) {
        if (!file.empty() && (file[0] == '/' || file[0] == '\\')) return file;
        return caseDir + "/" + file;
    }

    bool loadASCIISTL(const std::string& path) {
        std::ifstream in(path);
        if (!in) return false;

        std::string first;
        std::getline(in, first);
        if (first.find("solid") == std::string::npos) return false;
        in.clear();
        in.seekg(0);

        std::string word;
        std::size_t before = triangles_.size();
        while (in >> word) {
            if (word != "facet") continue;
            Triangle tri;
            std::string normalWord;
            in >> normalWord >> tri.normal.x >> tri.normal.y >> tri.normal.z;
            in >> word >> word;
            for (int v = 0; v < 3; ++v) in >> word >> tri.v[v].x >> tri.v[v].y >> tri.v[v].z;
            in >> word;
            in >> word;
            normalizeTriangle(tri);
            triangles_.push_back(tri);
        }
        return triangles_.size() > before;
    }

    bool loadBinarySTL(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) return false;

        char header[80];
        uint32_t count = 0;
        in.read(header, 80);
        in.read(reinterpret_cast<char*>(&count), sizeof(uint32_t));
        if (!in || count == 0) return false;

        std::size_t before = triangles_.size();
        for (uint32_t n = 0; n < count; ++n) {
            float data[12];
            uint16_t attr = 0;
            in.read(reinterpret_cast<char*>(data), sizeof(data));
            in.read(reinterpret_cast<char*>(&attr), sizeof(attr));
            if (!in) break;

            Triangle tri;
            tri.normal = {data[0], data[1], data[2]};
            tri.v[0] = {data[3], data[4], data[5]};
            tri.v[1] = {data[6], data[7], data[8]};
            tri.v[2] = {data[9], data[10], data[11]};
            normalizeTriangle(tri);
            triangles_.push_back(tri);
        }
        return triangles_.size() > before;
    }

    static void normalizeTriangle(Triangle& tri) {
        Point geomN = normalize(cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0]));
        tri.normal = (norm2(tri.normal) < 1e-24) ? geomN : normalize(tri.normal);
        if (dot(tri.normal, geomN) < 0.0) tri.normal = -1.0 * tri.normal;
    }

    void finalize() {
        if (triangles_.empty()) return;
        min_ = { std::numeric_limits<double>::max(),  std::numeric_limits<double>::max(),  std::numeric_limits<double>::max()};
        max_ = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};

        std::map<std::array<long long, 6>, int> edgeCount;
        double volume6 = 0.0;

        auto keyPoint = [](const Point& p) {
            constexpr double scale = 1e9;
            return std::array<long long, 3>{
                (long long)std::llround(p.x * scale),
                (long long)std::llround(p.y * scale),
                (long long)std::llround(p.z * scale)
            };
        };
        auto addEdge = [&](const Point& a, const Point& b) {
            auto ka = keyPoint(a);
            auto kb = keyPoint(b);
            if (kb < ka) std::swap(ka, kb);
            edgeCount[{ka[0], ka[1], ka[2], kb[0], kb[1], kb[2]}]++;
        };

        for (const auto& tri : triangles_) {
            for (const auto& v : tri.v) {
                min_.x = std::min(min_.x, v.x); min_.y = std::min(min_.y, v.y); min_.z = std::min(min_.z, v.z);
                max_.x = std::max(max_.x, v.x); max_.y = std::max(max_.y, v.y); max_.z = std::max(max_.z, v.z);
            }
            addEdge(tri.v[0], tri.v[1]);
            addEdge(tri.v[1], tri.v[2]);
            addEdge(tri.v[2], tri.v[0]);
            volume6 += dot(tri.v[0], cross(tri.v[1], tri.v[2]));
        }

        watertight_ = {};
        for (const auto& kv : edgeCount) {
            if (kv.second == 1) ++watertight_.boundaryEdges;
            else if (kv.second > 2) ++watertight_.nonManifoldEdges;
            if (kv.second > 1) watertight_.duplicateEdges += kv.second - 1;
        }
        watertight_.watertight = !edgeCount.empty()
                              && watertight_.boundaryEdges == 0
                              && watertight_.nonManifoldEdges == 0;
        normalsOutward_ = volume6 >= 0.0;

        std::ostringstream msg;
        msg << triangles_.size()
            << " triangle(s), watertight="
            << (watertight_.watertight ? "true" : "false")
            << ", normals="
            << (normalsOutward_ ? "outward" : "inward");
        broadcast("IBM STL loaded: ", msg.str());

        buildBVH();
    }

    void buildBVH() {
        nodes_.clear();
        triIndices_.resize(triangles_.size());
        for (std::size_t i = 0; i < triangles_.size(); ++i) triIndices_[i] = (int)i;
        if (!triIndices_.empty()) buildNode(0, (int)triIndices_.size());
    }

    int buildNode(int start, int count) {
        int id = (int)nodes_.size();
        nodes_.push_back({});
        nodes_[id].start = start;
        nodes_[id].count = count;
        nodes_[id].min = { std::numeric_limits<double>::max(),  std::numeric_limits<double>::max(),  std::numeric_limits<double>::max()};
        nodes_[id].max = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};

        Point cmin = nodes_[id].min;
        Point cmax = nodes_[id].max;
        for (int n = 0; n < count; ++n) {
            const Triangle& tri = triangles_[triIndices_[start + n]];
            Point c = (tri.v[0] + tri.v[1] + tri.v[2]) / 3.0;
            cmin.x = std::min(cmin.x, c.x); cmin.y = std::min(cmin.y, c.y); cmin.z = std::min(cmin.z, c.z);
            cmax.x = std::max(cmax.x, c.x); cmax.y = std::max(cmax.y, c.y); cmax.z = std::max(cmax.z, c.z);
            for (const auto& v : tri.v) {
                nodes_[id].min.x = std::min(nodes_[id].min.x, v.x); nodes_[id].min.y = std::min(nodes_[id].min.y, v.y); nodes_[id].min.z = std::min(nodes_[id].min.z, v.z);
                nodes_[id].max.x = std::max(nodes_[id].max.x, v.x); nodes_[id].max.y = std::max(nodes_[id].max.y, v.y); nodes_[id].max.z = std::max(nodes_[id].max.z, v.z);
            }
        }

        if (count <= 12) return id;

        Point extent = cmax - cmin;
        int axis = 0;
        if (extent.y > extent.x && extent.y >= extent.z) axis = 1;
        else if (extent.z > extent.x && extent.z >= extent.y) axis = 2;

        auto centroidCoord = [&](int triId) {
            const Triangle& tri = triangles_[triId];
            Point c = (tri.v[0] + tri.v[1] + tri.v[2]) / 3.0;
            return axis == 0 ? c.x : (axis == 1 ? c.y : c.z);
        };

        int mid = start + count / 2;
        std::nth_element(triIndices_.begin() + start, triIndices_.begin() + mid,
                         triIndices_.begin() + start + count,
                         [&](int a, int b) { return centroidCoord(a) < centroidCoord(b); });

        int left = buildNode(start, mid - start);
        int right = buildNode(mid, start + count - mid);
        nodes_[id].left = left;
        nodes_[id].right = right;
        nodes_[id].count = 0;
        return id;
    }

    bool insideBounds(const Point& p) const {
        constexpr double pad = 1e-12;
        return p.x >= min_.x - pad && p.x <= max_.x + pad
            && p.y >= min_.y - pad && p.y <= max_.y + pad
            && p.z >= min_.z - pad && p.z <= max_.z + pad;
    }

    static Point closestPointOnTriangle(const Point& p, const Point& a, const Point& b, const Point& c) {
        Point ab = b - a;
        Point ac = c - a;
        Point ap = p - a;
        double d1 = dot(ab, ap);
        double d2 = dot(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0) return a;

        Point bp = p - b;
        double d3 = dot(ab, bp);
        double d4 = dot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3) return b;

        double vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
            double v = d1 / (d1 - d3);
            return a + ab * v;
        }

        Point cp = p - c;
        double d5 = dot(ab, cp);
        double d6 = dot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6) return c;

        double vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
            double w = d2 / (d2 - d6);
            return a + ac * w;
        }

        double va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
            double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + (c - b) * w;
        }

        double denom = 1.0 / (va + vb + vc);
        double v = vb * denom;
        double w = vc * denom;
        return a + ab * v + ac * w;
    }

    void nearestTriangle(int nodeId, const Point& p, SDFResult& out) const {
        if (nodeId < 0 || nodeId >= (int)nodes_.size()) return;
        const BVHNode& node = nodes_[nodeId];
        double best2 = out.distance * out.distance;
        if (distance2ToBox(p, node.min, node.max) > best2) return;

        if (node.left < 0 && node.right < 0) {
            for (int n = 0; n < node.count; ++n) {
                int triId = triIndices_[node.start + n];
                const Triangle& tri = triangles_[triId];
                Point c = closestPointOnTriangle(p, tri.v[0], tri.v[1], tri.v[2]);
                double d2 = norm2(p - c);
                if (d2 < out.distance * out.distance) {
                    out.distance = std::sqrt(d2);
                    out.closestPoint = c;
                    out.normal = tri.normal;
                    out.triangleIndex = (std::size_t)triId;
                }
            }
            return;
        }

        double dl = node.left >= 0 ? distance2ToBox(p, nodes_[node.left].min, nodes_[node.left].max)
                                   : std::numeric_limits<double>::max();
        double dr = node.right >= 0 ? distance2ToBox(p, nodes_[node.right].min, nodes_[node.right].max)
                                    : std::numeric_limits<double>::max();
        if (dl < dr) {
            nearestTriangle(node.left, p, out);
            nearestTriangle(node.right, p, out);
        } else {
            nearestTriangle(node.right, p, out);
            nearestTriangle(node.left, p, out);
        }
    }

    int countRayHits(int nodeId, const Point& p, const Point& dir) const {
        if (nodeId < 0 || nodeId >= (int)nodes_.size()) return 0;
        const BVHNode& node = nodes_[nodeId];
        if (!rayIntersectsBox(p, dir, node.min, node.max)) return 0;
        if (node.left < 0 && node.right < 0) {
            int hits = 0;
            for (int n = 0; n < node.count; ++n) {
                double t = 0.0;
                if (rayIntersectsTriangle(p, dir, triangles_[triIndices_[node.start + n]], t) && t > 1e-12) ++hits;
            }
            return hits;
        }
        return countRayHits(node.left, p, dir) + countRayHits(node.right, p, dir);
    }

    static double distance2ToBox(const Point& p, const Point& bmin, const Point& bmax) {
        double d2 = 0.0;
        auto axis = [&](double x, double lo, double hi) {
            if (x < lo) return (lo - x) * (lo - x);
            if (x > hi) return (x - hi) * (x - hi);
            return 0.0;
        };
        d2 += axis(p.x, bmin.x, bmax.x);
        d2 += axis(p.y, bmin.y, bmax.y);
        d2 += axis(p.z, bmin.z, bmax.z);
        return d2;
    }

    static bool rayIntersectsBox(const Point& p, const Point& dir, const Point& bmin, const Point& bmax) {
        double tmin = 0.0;
        double tmax = std::numeric_limits<double>::max();
        auto slab = [&](double origin, double d, double lo, double hi) {
            if (std::abs(d) < 1e-14) return origin >= lo && origin <= hi;
            double inv = 1.0 / d;
            double t0 = (lo - origin) * inv;
            double t1 = (hi - origin) * inv;
            if (t0 > t1) std::swap(t0, t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            return tmax >= tmin;
        };
        return slab(p.x, dir.x, bmin.x, bmax.x)
            && slab(p.y, dir.y, bmin.y, bmax.y)
            && slab(p.z, dir.z, bmin.z, bmax.z);
    }

    static bool rayIntersectsTriangle(const Point& origin, const Point& dir,
                                      const Triangle& tri, double& t) {
        constexpr double eps = 1e-12;
        Point e1 = tri.v[1] - tri.v[0];
        Point e2 = tri.v[2] - tri.v[0];
        Point h = cross(dir, e2);
        double a = dot(e1, h);
        if (std::abs(a) < eps) return false;
        double f = 1.0 / a;
        Point s = origin - tri.v[0];
        double u = f * dot(s, h);
        if (u < -eps || u > 1.0 + eps) return false;
        Point q = cross(s, e1);
        double v = f * dot(dir, q);
        if (v < -eps || u + v > 1.0 + eps) return false;
        t = f * dot(e2, q);
        return t > eps;
    }
};

} // namespace GeoProcessing
} // namespace IBM
} // namespace SF
