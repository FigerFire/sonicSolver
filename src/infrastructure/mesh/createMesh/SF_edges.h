/// @file SF_edges.h
/// @brief 结构网格 block、edge 与生成流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include <cmath>
#include <stdexcept>
#include "SF_valueTypes.h"
#include "SF_meshGen.h"

namespace SF {

// ============================================================
//  ArcGeometry — 圆弧参数 (供 ArcEdge / 内联函数使用)
// ============================================================
struct ArcGeometry {
    Vector3 center{0,0,0};
    double  radius = 0;
    Vector3 axis{0,0,1};
    double  angle = 0;
    Vector3 startDir{1,0,0};
    Vector3 bitangent{0,1,0};
};

// ============================================================
//  Edge 抽象基类 (OOP 风格, 供 SF_block.cpp 使用)
// ============================================================
class Edge {
public:
    virtual ~Edge() = default;
    virtual Vector3 evaluate(double t) const = 0;   // t∈[0,1]
};

class LineEdge : public Edge {
    Vector3 v0, v1;
public:
    LineEdge(const Vector3& a, const Vector3& b) : v0(a), v1(b) {}
    Vector3 evaluate(double t) const override {
        return { v0.x + t*(v1.x - v0.x),
                 v0.y + t*(v1.y - v0.y),
                 v0.z + t*(v1.z - v0.z) };
    }
};

class ArcEdge : public Edge {
    ArcGeometry geom;
public:
    ArcEdge() = default;
    ArcEdge(const Vector3& v0, const Vector3& v1, const Vector3& mid) {
        buildFromThreePoints(v0, v1, mid);
    }
    ArcEdge(const Vector3& v0, const Vector3& v1, const Vector3& origin, bool useCenter) {
        if (useCenter) buildFromCenter(v0, v1, origin);
    }

    void buildFromThreePoints(const Vector3& v0, const Vector3& v1, const Vector3& mid);
    void buildFromCenter(const Vector3& v0, const Vector3& v1, const Vector3& origin);
    void buildFromAngle(const Vector3& v0, const Vector3& v1, double angle, const Vector3& axis);

    Vector3 evaluate(double t) const override {
        double theta = t * geom.angle;
        double c = std::cos(theta), s = std::sin(theta);
        return { geom.center.x + geom.radius*(c*geom.startDir.x + s*geom.bitangent.x),
                 geom.center.y + geom.radius*(c*geom.startDir.y + s*geom.bitangent.y),
                 geom.center.z + geom.radius*(c*geom.startDir.z + s*geom.bitangent.z) };
    }
};

// ============================================================
//  内联工具函数 (供 SF_meshGen.cpp 的 generateTFIBlock 直接使用)
// ============================================================

/// 过三点 (v0, v1, mid) 计算圆弧几何参数
/// @return radius, angle, center, axis, startDir, bitangent
inline void arcGeometry(const Vector3& v0, const Vector3& v1, const Vector3& mid,
                         Vector3& center, double& radius,
                         Vector3& axis, double& angle,
                         Vector3& startDir, Vector3& bitangent) {
    Vector3 a = v1 - v0;
    Vector3 b = mid - v0;
    axis = cross(a, b);
    double len = std::sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
    if (len < 1e-20) { center = v0; radius = 0; angle = 0; startDir = {1,0,0}; bitangent = {0,1,0}; return; }
    axis.x /= len; axis.y /= len; axis.z /= len;

    Vector3 mAB = (v0 + v1) * 0.5;
    Vector3 mAC = (v0 + mid) * 0.5;
    Vector3 dAB = v1 - v0;
    Vector3 dAC = mid - v0;
    Vector3 nAB = cross(axis, dAB);
    Vector3 nAC = cross(axis, dAC);

    // Solve midAB + t1*nAB = midAC + t2*nAC → use yz plane
    double det = nAB.y * nAC.z - nAB.z * nAC.y;
    if (std::abs(det) < 1e-12) det = nAB.x * nAC.z - nAB.z * nAC.x;
    if (std::abs(det) < 1e-12) det = nAB.x * nAC.y - nAB.y * nAC.x;
    double t1 = 0;
    if (std::abs(det) > 1e-20) {
        double dx = mAC.x - mAB.x, dy = mAC.y - mAB.y, dz = mAC.z - mAB.z;
        if (std::abs(nAB.y * nAC.z - nAB.z * nAC.y) > 1e-12)
            t1 = (dy * nAC.z - dz * nAC.y) / det;
        else
            t1 = (dx * nAC.y - dy * nAC.x) / (det + 1e-20);
    }
    center = { mAB.x + t1 * nAB.x, mAB.y + t1 * nAB.y, mAB.z + t1 * nAB.z };

    startDir = v0 - center;
    radius = std::sqrt(startDir.x*startDir.x + startDir.y*startDir.y + startDir.z*startDir.z);
    if (radius < 1e-20) { angle = 0; startDir = {1,0,0}; bitangent = {0,1,0}; return; }
    startDir.x /= radius; startDir.y /= radius; startDir.z /= radius;

    Vector3 endDir = v1 - center;
    double endLen = std::sqrt(endDir.x*endDir.x + endDir.y*endDir.y + endDir.z*endDir.z);
    endDir.x /= (endLen + 1e-20); endDir.y /= (endLen + 1e-20); endDir.z /= (endLen + 1e-20);

    double dotProd = startDir.x*endDir.x + startDir.y*endDir.y + startDir.z*endDir.z;
    if (dotProd > 1.0) dotProd = 1.0; if (dotProd < -1.0) dotProd = -1.0;
    double shortAngle = std::acos(dotProd);

    // 判断 mid 点是否在从 startDir→endDir 的短弧上 (绕 axis 方向)
    Vector3 crSE = cross(startDir, endDir);
    double signSE = crSE.x*axis.x + crSE.y*axis.y + crSE.z*axis.z;

    Vector3 rMid = mid - center;
    double rMidLen = std::sqrt(rMid.x*rMid.x + rMid.y*rMid.y + rMid.z*rMid.z);
    Vector3 midDir = {rMid.x / (rMidLen + 1e-20), rMid.y / (rMidLen + 1e-20), rMid.z / (rMidLen + 1e-20)};
    Vector3 crSM = cross(startDir, midDir);
    double signSM = crSM.x*axis.x + crSM.y*axis.y + crSM.z*axis.z;

    // mid 和 endDir 在 axis 同侧 → mid 在短弧上；否则在长弧上
    bool midOnShortArc = (signSE * signSM >= 0);

    if (midOnShortArc) {
        angle = shortAngle;
        // 确保 bitangent 指向 endDir 的方向
        if (signSE < 0) {
            axis.x = -axis.x; axis.y = -axis.y; axis.z = -axis.z;
        }
    } else {
        angle = 2.0 * M_PI - shortAngle;
        // 确保 bitangent 指向 mid 的方向 (先经过 mid 再到达 endDir)
        if (signSM < 0) {
            axis.x = -axis.x; axis.y = -axis.y; axis.z = -axis.z;
        }
    }

    bitangent = cross(axis, startDir);
}

/// 沿弧线求点 P(θ) = center + R*(cosθ*startDir + sinθ*bitangent)
inline Vector3 arcPoint(const Vector3& center, double radius,
                         const Vector3& startDir, const Vector3& bitangent, double theta) {
    double c = std::cos(theta), s = std::sin(theta);
    return { center.x + radius*(c*startDir.x + s*bitangent.x),
             center.y + radius*(c*startDir.y + s*bitangent.y),
             center.z + radius*(c*startDir.z + s*bitangent.z) };
}

/// 便捷函数: 沿三点圆弧求 t∈[0,1] 处的点
inline Vector3 arcEval(const Vector3& v0, const Vector3& v1, const Vector3& mid, double t) {
    Vector3 center, axis, startDir, bitangent; double radius, angle;
    arcGeometry(v0, v1, mid, center, radius, axis, angle, startDir, bitangent);
    return arcPoint(center, radius, startDir, bitangent, t * angle);
}

/// 沿边求点: 直线 或 弧线, t∈[0,1]
inline Vector3 edgeEval(const Vector3& v0, const Vector3& v1, const MeshEdge& e, double t) {
    if (e.type == EdgeType::ARC)
        return arcEval(v0, v1, e.arcP, t);
    // LINE: linear interpolation
    return { v0.x + t*(v1.x - v0.x), v0.y + t*(v1.y - v0.y), v0.z + t*(v1.z - v0.z) };
}

} // namespace SF
