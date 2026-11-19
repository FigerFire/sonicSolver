#pragma once
/// @file SF_scalarUpwind.h
/// @brief 节点面速度/质量通量的标量 Upwind1；不选择状态版本或时间算法。
#include "solver/discretization/SF_formulaOperator.h"
namespace SF::Discretization {
/// @brief 原metric节点面上风形式；保持每方向right-left及Jac归一化顺序。
inline System::FormulaValueKernel scalarMetricUpwind1(std::size_t velocity,std::size_t scalar) {
    return [velocity,scalar](const auto& v,int cell,int) {
        const auto face=[&](int lower,int axis){
            const double un=v.metricFaceVelocity(velocity,lower,axis);
            const int upper=v.spatialNeighbor(lower,axis,1);
            const double left=v.boundReads.at(scalar)(lower,0),right=v.boundReads.at(scalar)(upper,0);
            return un>=0?un*left:un*right;
        };
        double divergence=0;
        for(int axis=0;axis<3;++axis)if(v.spacing(axis)!=0)
            divergence+=face(cell,axis)-face(v.spatialNeighbor(cell,axis,-1),axis);
        return v.inverseVolume(cell)*divergence;
    };
}
/// A single face expression is reused with +/- divergence signs. This is an
/// explicit scalar Upwind1 recipe, not a replacement for the NS flux recipe.
inline double conservativeScalarAdvection(const System::FormulaValues& v,int cell) {
    double divergence=0;
    for(int axis=0;axis<3;++axis) {
        const double h=v.spacing(axis);if(h==0)continue;
        const double rho=v.boundReads[1](cell,0),m=v.boundReads[2](cell,axis);
        if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Scalar transport density is not positive.");
        const double c=v.boundReads[0](cell,0)/rho;
        for(int sign:{-1,1}) {
            const double rn=v.boundNeighbor(1,cell,axis,sign,0),mn=v.boundNeighbor(2,cell,axis,sign,axis);
            if(!std::isfinite(rn) || rn<=0)throw std::runtime_error("Scalar transport neighbor density is not positive.");
            const double cn=v.boundNeighbor(0,cell,axis,sign,0)/rn;
            const double mass=.5*(m+mn);
            const double left=sign>0?c:cn,right=sign>0?cn:c;
            divergence+=sign*mass*(mass>=0?left:right)/h;
        }
    }return divergence;
}

/// @brief 与原中心面速 Upwind1 相同的节点面表达；输入是显式绑定速度，不读取 packed Q。
inline System::FormulaValueKernel scalarVelocityUpwind1(std::size_t velocity,std::size_t scalar) {
    return [velocity,scalar](const auto& v,int cell,int) {
        double divergence=0;
        for(int axis=0;axis<3;++axis) {
            const double h=v.spacing(axis);if(h==0)continue;
            const double u=v.boundReads.at(velocity)(cell,axis),a=v.boundReads.at(scalar)(cell,0);
            for(int sign:{-1,1}) {
                const double un=v.boundNeighbor(velocity,cell,axis,sign,axis);
                const double an=v.boundNeighbor(scalar,cell,axis,sign,0);
                const double face=.5*(u+un),left=sign>0?a:an,right=sign>0?an:a;
                divergence+=sign*face*(face>=0?left:right)/h;
            }
        }return divergence;
    };
}

/// @brief 同一校正体积通量的上风标量散度，Jac 归一化仅一次。
inline System::FormulaValueKernel scalarVolumeUpwind1(std::size_t scalar) {
    return [scalar](const auto& v,int cell,int) {
        double net=0;
        for(int axis=0;axis<3;++axis) {
            if(v.spacing(axis)==0)continue;
            for(int sign:{-1,1}) {
                const int neighbor=v.spatialNeighbor(cell,axis,sign);
                const int lower=sign>0?cell:neighbor;
                const double phi=v.volumeFlux(axis,lower);
                const int upwind=phi>=0?lower:v.spatialNeighbor(lower,axis,1);
                net+=sign*phi*v.boundReads.at(scalar)(upwind,0);
            }
        }return net*v.inverseVolume(cell);
    };
}

}
