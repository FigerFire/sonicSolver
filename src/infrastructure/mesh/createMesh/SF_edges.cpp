/// @file SF_edges.cpp
/// @brief 结构网格 block、edge 与生成流程实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_edges.h"
#include <cmath>

namespace SF {

void ArcEdge::buildFromThreePoints(const Vector3& v0, const Vector3& v1, const Vector3& mid) {
    Vector3 v01 = v1 - v0, v0mid = mid - v0;
    Vector3 normal = cross(v01, v0mid);
    double len = norm(normal) + 1e-30;
    geom.axis = (1.0 / len) * normal;

    Vector3 m01 = (v0 + v1) * 0.5, m0m = (v0 + mid) * 0.5;
    Vector3 n01 = cross(geom.axis, v01);
    Vector3 n0m = cross(geom.axis, v0mid);

    double det = dot(cross(n01, n0m), geom.axis);
    double t1 = 0;
    if (std::abs(det) > 1e-20) {
        Vector3 delta = m0m - m01;
        t1 = dot(cross(delta, n0m), geom.axis) / det;
    }
    geom.center = { m01.x + t1*n01.x, m01.y + t1*n01.y, m01.z + t1*n01.z };

    Vector3 r0 = v0 - geom.center;
    geom.radius = norm(r0);
    geom.startDir = normalize(r0);
    geom.bitangent = cross(geom.axis, geom.startDir);

    Vector3 r1 = v1 - geom.center;
    Vector3 endDir = normalize(r1);
    double cosA = dot(geom.startDir, endDir);
    if (cosA > 1.0) cosA = 1.0; if (cosA < -1.0) cosA = -1.0;
    double fullAngle = std::acos(cosA);
    if (dot(cross(geom.startDir, endDir), geom.axis) < 0) fullAngle = 2.0*M_PI - fullAngle;

    Vector3 rMid = normalize(mid - geom.center);
    double sinMid = dot(cross(geom.startDir, rMid), geom.axis);
    double cosMid = dot(geom.startDir, rMid);
    double midAngle = std::atan2(sinMid, cosMid);
    if (midAngle < 0) midAngle += 2.0*M_PI;

    geom.angle = (midAngle < fullAngle) ? fullAngle : (fullAngle - 2.0*M_PI);
    if (geom.angle < 0) geom.angle += 2.0*M_PI;
}

void ArcEdge::buildFromCenter(const Vector3& v0, const Vector3& v1, const Vector3& origin) {
    geom.center = origin;
    Vector3 r0 = v0 - origin, r1 = v1 - origin;
    double r0Len = norm(r0), r1Len = norm(r1);
    geom.radius = 0.5 * (r0Len + r1Len);
    geom.startDir = (1.0 / r0Len) * r0;
    geom.axis = normalize(cross(r0, r1));
    geom.bitangent = cross(geom.axis, geom.startDir);

    double cosA = dot(geom.startDir, r1) / r1Len;
    double sinA = dot(cross(geom.startDir, r1), geom.axis);
    geom.angle = std::atan2(sinA, cosA);
    if (std::abs(geom.angle) < 1e-8) geom.angle = 2.0*M_PI;
    if (geom.angle < 0) geom.angle += 2.0*M_PI;
}

void ArcEdge::buildFromAngle(const Vector3& v0, const Vector3& v1, double angle, const Vector3& axis) {
    geom.axis = normalize(axis);
    Vector3 chord = v1 - v0;
    double d = norm(chord);
    double half = 0.5 * angle;
    double R = 0.5 * d / std::sin(half);
    double h = R * std::cos(half);
    Vector3 chordDir = chord * (1.0/d);
    Vector3 perpDir = normalize(cross(geom.axis, chordDir));
    Vector3 M = (v0 + v1) * 0.5;

    geom.center = M + h * perpDir;
    geom.radius = norm(v0 - geom.center);
    geom.startDir = normalize(v0 - geom.center);
    geom.bitangent = cross(geom.axis, geom.startDir);
    geom.angle = angle;
}

} // namespace SF
