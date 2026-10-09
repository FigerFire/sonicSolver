#pragma once
/// @file SF_viscousWall.h
/// @brief 静止无滑移、绝热壁面的约束二次拟合；计划只保存几何权重。
#include "core/field/SF_field.h"
#include "SF_ibmConfig.h"
#include "SF_ilwStorage.h"
#include "SF_ibmTopology.h"
#include "SF_localFrame.h"
#include "methods/math/discrete/SF_polynomial.h"
#include "core/mesh/SF_dimension.h"
#include <cmath>
#include <stdexcept>

namespace SF::IBM::GhostIBM {
/// @brief 同一固定支撑上，Dirichlet 速度和 Neumann 热力学量的 ghost 权重。
struct ViscousWallPlan {
    int ghost=0;
    std::vector<int> cells;
    std::vector<double> dirichlet,neumann,normalDerivative;
    int firstFluid=0;
    double firstDistance=0;
};
inline std::vector<double> wallBasis(double s,double t,double b,bool threeD,bool dirichlet) {
    if (dirichlet) return threeD?std::vector<double>{s,s*s,s*t,s*b}:std::vector<double>{s,s*s,s*t};
    return threeD?std::vector<double>{1,t,b,s*s,t*t,b*b,t*b}:std::vector<double>{1,t,s*s,t*t};
}
/// @brief 投影固定边界条件下的多项式；不满秩支撑直接失败。
inline std::vector<double> wallProjection(const std::vector<std::vector<double>>& rows,
                                          const std::vector<double>& atGhost) {
    const auto n=atGhost.size();
    std::vector<std::vector<double>> gram(n,std::vector<double>(n,0)),inverse;
    for (const auto& row:rows) for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<n;++j)
        gram[i][j]+=row[i]*row[j];
    if (!Math::Polynomial::invertNormalMatrix(gram,inverse))
        throw std::runtime_error("stationaryNoSlipAdiabatic requires full-rank quadratic fluid support; no order downgrade.");
    std::vector<double> weights(rows.size(),0);
    for (std::size_t q=0;q<rows.size();++q) for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<n;++j)
        weights[q]+=atGhost[i]*inverse[i][j]*rows[q][j];
    for (std::size_t j=0;j<n;++j) {
        double residual=-atGhost[j];
        for (std::size_t q=0;q<rows.size();++q) residual+=weights[q]*rows[q][j];
        if (!std::isfinite(residual) || std::abs(residual)>1e-9)
            throw std::runtime_error("Viscous wall projection failed polynomial reproduction.");
    }
    return weights;
}
inline ViscousWallPlan buildViscousWallPlan(const Field& field,const GhostILW::Storage& storage,
                                          const IBMGeometry& geometry,int i,int j,int k) {
    ViscousWallPlan result;result.ghost=field.getIdx(i,j,k);
    const auto wall=geometry.wallPoint(i,j,k),normal=geometry.wallNormal(i,j,k);
    const double normalArray[3]={normal.x,normal.y,normal.z},wallArray[3]={wall.x,wall.y,wall.z};
    const bool threeD=Math::activeDimensionCount()==3;
    const auto inactive=Math::inactiveDirectionNormalVector();
    const auto frame=threeD?Math::makeLocalFrame(normalArray):Math::makeLocalFrame2D(normalArray,inactive);
    const double targetDistance=std::abs(geometry.signedDistance(i,j,k));
    if (!std::isfinite(targetDistance)) throw std::runtime_error("Viscous wall requires finite geometric distance.");
    const double position[3]={field.X(i,j,k),field.Y(i,j,k),field.Z(i,j,k)};
    double ghostS,ghostT,ghostB;
    Math::localCoordinates(position,wallArray,frame,ghostS,ghostT,ghostB);
    const double surfaceTolerance=32*std::numeric_limits<double>::epsilon()
        *std::max({1.0,std::abs(position[0]),std::abs(position[1]),std::abs(position[2])});
    if (!std::isfinite(ghostS) || !std::isfinite(ghostT) || !std::isfinite(ghostB)
        || ghostS>surfaceTolerance)
        throw std::runtime_error("Viscous ghost is not on the solid side of its wall frame; inconsistent geometry, no reflection fallback.");
    double scale=0;
    std::vector<std::array<double,3>> points;
    std::vector<std::vector<double>> velocityRows,thermalRows;
    result.firstDistance=std::numeric_limits<double>::max();
    const int count=std::min(48,storage.fluidCount(result.ghost));
    for (int q=0;q<count;++q) {
        const int cell=storage.fluidSample(result.ghost,q);int ii,jj,kk;field.getIJK(cell,ii,jj,kk);
        if (field.CellFlag(ii,jj,kk)!=FLUID_CELL) throw std::runtime_error("Viscous wall donor is not fluid.");
        const double position[3]={field.X(ii,jj,kk),field.Y(ii,jj,kk),field.Z(ii,jj,kk)};
        double s,t,b;Math::localCoordinates(position,wallArray,frame,s,t,b);
        if (s<=0) continue;
        result.cells.push_back(cell);
        points.push_back({s,t,b});scale=std::max(scale,s);
        const double distance=field.wallDistance(ii,jj,kk);
        if (distance>0 && distance<result.firstDistance) {result.firstDistance=distance;result.firstFluid=cell;}
    }
    if (result.cells.empty() || result.firstDistance==std::numeric_limits<double>::max())
        throw std::runtime_error("Viscous wall has no valid first fluid distance.");
    for (const auto& point:points) {
        velocityRows.push_back(wallBasis(point[0]/scale,point[1]/scale,point[2]/scale,threeD,true));
        thermalRows.push_back(wallBasis(point[0]/scale,point[1]/scale,point[2]/scale,threeD,false));
    }
    // Evaluate at the actual node, including its tangential offset from the wall anchor.
    result.dirichlet=wallProjection(velocityRows,wallBasis(ghostS/scale,ghostT/scale,ghostB/scale,threeD,true));
    result.neumann=wallProjection(thermalRows,wallBasis(ghostS/scale,ghostT/scale,ghostB/scale,threeD,false));
    auto derivative=wallBasis(0,0,0,threeD,true);derivative[0]=1/scale;
    result.normalDerivative=wallProjection(velocityRows,derivative);
    return result;
}
/// @brief rho/p 零法向梯度与 U=0；PerfectGas 下等价于绝热 T 条件。
inline void applyViscousWall(Field& field,const std::vector<ViscousWallPlan>& plans,double gamma) {
    if (!std::isfinite(gamma) || gamma<=1 || field.NVar()!=5)
        throw std::runtime_error("Viscous immersed closure requires five-component PerfectGas state and gamma > 1.");
    for (const auto& plan:plans) {
        double rho=0,p=0,u[3]={0,0,0};
        for (std::size_t q=0;q<plan.cells.size();++q) {
            int i,j,k;field.getIJK(plan.cells[q],i,j,k);
            const double density=field(i,j,k,RHO);
            if (!std::isfinite(density) || density<=0)
                throw std::runtime_error("Viscous wall requires positive donor density.");
            const double velocity[3]={field(i,j,k,RU)/density,field(i,j,k,RV)/density,field(i,j,k,RW)/density};
            const double pressure=(gamma-1)*(field(i,j,k,E)-0.5*density*(velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]));
            if (!std::isfinite(density) || density<=0 || !std::isfinite(pressure) || pressure<=0)
                throw std::runtime_error("Viscous wall requires physical PerfectGas donors.");
            rho+=plan.neumann[q]*density;p+=plan.neumann[q]*pressure;
            for (int d=0;d<3;++d) u[d]+=plan.dirichlet[q]*velocity[d];
        }
        if (!std::isfinite(rho) || rho<=0 || !std::isfinite(p) || p<=0)
            throw std::runtime_error("Viscous wall quadratic closure produced nonphysical rho/p; no clamp or fallback.");
        int i,j,k;field.getIJK(plan.ghost,i,j,k);
        field(i,j,k,RHO)=rho;for (int d=0;d<3;++d) field(i,j,k,RU+d)=rho*u[d];
        field(i,j,k,E)=p/(gamma-1)+0.5*rho*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    }
}
/// @brief SST 低 Reynolds 壁面：k=0、omega=60 nu/(beta1 d1²)、mu_t=0。
inline void applySSTWall(const Field& field,const std::vector<ViscousWallPlan>& plans,double laminarMu,
        std::vector<double>& k,std::vector<double>& omega,std::vector<double>& muT) {
    if (!std::isfinite(laminarMu) || laminarMu<=0 || k.size()!=static_cast<std::size_t>(field.TotalSize())
        || omega.size()!=k.size() || muT.size()!=k.size())
        throw std::runtime_error("SST immersed boundary requires positive molecular viscosity and matching original storage.");
    for (const auto& plan:plans) {
        int i,j,l;field.getIJK(plan.firstFluid,i,j,l);
        const double rho=field(i,j,l,RHO);
        if (!std::isfinite(rho) || rho<=0) throw std::runtime_error("SST wall density is nonphysical.");
        const double wallOmega=60*laminarMu/(rho*0.075*plan.firstDistance*plan.firstDistance);
        double kGhost=0,omegaGhost=wallOmega;
        for (std::size_t q=0;q<plan.cells.size();++q) {
            kGhost+=plan.dirichlet[q]*k[plan.cells[q]];
            omegaGhost+=plan.dirichlet[q]*(omega[plan.cells[q]]-wallOmega);
        }
        if (!std::isfinite(kGhost) || !std::isfinite(omegaGhost))
            throw std::runtime_error("SST immersed polynomial boundary produced invalid omega; no fallback.");
        k[plan.ghost]=kGhost;omega[plan.ghost]=omegaGhost;muT[plan.ghost]=0;
    }
}

}
