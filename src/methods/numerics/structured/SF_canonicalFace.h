#pragma once

/// @file SF_canonicalFace.h
/// @brief 结构网格唯一有向面的几何、索引和 owner-neighbour 符号契约。

#include "SF_field.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace SF::Numerics::CanonicalFace {

/// @brief 面方向总是从 lower 结构点指向 upper 结构点。
enum class CellSide {
    LowerOwner,
    UpperNeighbour
};

/// @brief 一个 canonical 面的守恒几何。
struct Geometry {
    std::array<double,3> cofactor{0.0,0.0,0.0};
    std::array<double,3> unitNormal{0.0,0.0,0.0};
    double area=0.0;
    double inverseJacobian=0.0;
};

inline double metricCofactor(double metric,double inverseJacobian) {
    if (!std::isfinite(metric)||!std::isfinite(inverseJacobian)
        ||std::abs(inverseJacobian)<1.0e-30) {
        throw std::runtime_error(
            "Canonical face received invalid metric/Jacobian.");
    }
    return metric/inverseJacobian;
}

/// @brief 返回 lower 点所拥有的唯一面索引。
inline std::size_t index(
        const Field& field,int direction,int i,int j,int k) {
    if (direction<0||direction>=3
        ||i<0||i>=field.MX()
        ||j<0||j>=field.MY()
        ||k<0||k>=field.MZ()) {
        throw std::runtime_error(
            "Canonical face descriptor is outside Field storage.");
    }
    return static_cast<std::size_t>(direction)
        *static_cast<std::size_t>(field.TotalSize())
        +static_cast<std::size_t>(field.getIdx(i,j,k));
}

/// @brief 构造 lower->upper 的面共因子、面积和单位法向。
inline Geometry geometry(
        const Field& field,int direction,int i,int j,int k) {
    (void)index(field,direction,i,j,k);
    int di=0,dj=0,dk=0;
    if(direction==0)di=1;
    else if(direction==1)dj=1;
    else dk=1;
    if(i+di>=field.MX()||j+dj>=field.MY()||k+dk>=field.MZ()) {
        throw std::runtime_error(
            "Canonical face upper point is outside Field storage.");
    }

    double stored[4]{0.0,0.0,0.0,0.0};
    if(field.hasCanonicalFaceMetrics(direction,i,j,k)) {
        field.canonicalFaceMetrics(direction,i,j,k,stored);
    } else {
        auto component=[&](int axis,int ii,int jj,int kk) {
            if(direction==0) {
                if(axis==0)return field.XiX(ii,jj,kk);
                if(axis==1)return field.XiY(ii,jj,kk);
                return field.XiZ(ii,jj,kk);
            }
            if(direction==1) {
                if(axis==0)return field.EtX(ii,jj,kk);
                if(axis==1)return field.EtY(ii,jj,kk);
                return field.EtZ(ii,jj,kk);
            }
            if(axis==0)return field.ZeX(ii,jj,kk);
            if(axis==1)return field.ZeY(ii,jj,kk);
            return field.ZeZ(ii,jj,kk);
        };
        for(int axis=0;axis<3;++axis) {
            stored[axis]=0.5*(
                metricCofactor(
                    component(axis,i,j,k),field.Jac(i,j,k))
                +metricCofactor(
                    component(axis,i+di,j+dj,k+dk),
                    field.Jac(i+di,j+dj,k+dk)));
        }
        stored[3]=0.5*(
            field.Jac(i,j,k)+field.Jac(i+di,j+dj,k+dk));
    }

    Geometry result;
    result.cofactor={stored[0],stored[1],stored[2]};
    result.inverseJacobian=stored[3];
    if(!std::isfinite(result.inverseJacobian)
        ||std::abs(result.inverseJacobian)<1.0e-30) {
        throw std::runtime_error(
            "Canonical face has invalid inverse Jacobian.");
    }
    result.area=std::sqrt(
        stored[0]*stored[0]+stored[1]*stored[1]+stored[2]*stored[2]);
    if(!std::isfinite(result.area)||result.area<=0.0) {
        throw std::runtime_error(
            "Canonical face has non-positive area.");
    }
    for(int axis=0;axis<3;++axis) {
        result.unitNormal[(size_t)axis]=
            result.cofactor[(size_t)axis]/result.area;
    }
    return result;
}

/// @brief 将 lower->upper 通量转换为指定相邻单元的外向通量。
inline double outward(double canonicalFlux,CellSide side) {
    if(!std::isfinite(canonicalFlux)) {
        throw std::runtime_error(
            "Canonical face flux is non-finite.");
    }
    return side==CellSide::LowerOwner
        ?canonicalFlux:-canonicalFlux;
}

/// @brief 将 canonical owner 的通量映射到局部面方向。
inline double orient(double canonicalFlux,double orientation) {
    if(!std::isfinite(orientation)
        ||std::abs(std::abs(orientation)-1.0)>1.0e-12) {
        throw std::runtime_error(
            "Canonical face orientation must be +1 or -1.");
    }
    return orientation*canonicalFlux;
}

/// @brief 检查一对 owner-neighbour 外向通量是否严格相消。
inline void requireConservativePair(
        double lowerOutward,double upperOutward,
        double tolerance=1.0e-12) {
    const double scale=std::max({
        std::abs(lowerOutward),std::abs(upperOutward),1.0});
    if(!std::isfinite(lowerOutward)||!std::isfinite(upperOutward)
        ||std::abs(lowerOutward+upperOutward)>tolerance*scale) {
        throw std::runtime_error(
            "Canonical owner-neighbour flux pair is not conservative.");
    }
}

} // namespace SF::Numerics::CanonicalFace
