#pragma once
/// @brief 可选的只读 publication 诊断；关闭时不遍历或复制数组。
#include "SF_state.h"
#include "core/residual/SF_residual.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <unistd.h>
namespace SF::Physics::Multiphase {
inline void tracePublication(const char* operation,const Field& field,
        const LevelSetField& state,const std::vector<double>* rhs=nullptr,const Residual* source=nullptr) {
    const char* directory=std::getenv("SF_LEVEL_SET_CHECKPOINT_DIR");
    if (!directory) return;
    static std::map<const Field*,unsigned long long> sequences;
    const auto sequence=sequences[&field]++;
    const char* rank=std::getenv("OMPI_COMM_WORLD_RANK");
    std::ofstream out(std::string(directory)+"/levelSet-"+std::to_string(getpid())+"-"
        +std::to_string(reinterpret_cast<std::uintptr_t>(&field))+"-"+std::to_string(sequence)+"-"+operation+".csv");
    if (!out) throw std::runtime_error("Cannot write level-set publication checkpoint.");
    out<<std::setprecision(17)<<"sequence,operation,rank,id,ownerRank,owner,i,j,k,x,y,z,interior,halo,boundary,flag,phi,nx,ny,nz,curvature,mask,density,viscosity,rhs,fx,fy,fz,power\n";
    const int ng=field.NG();
    for(int k=ng;k<ng+field.NZ();++k) for(int j=0;j<field.MY();++j) for(int i=0;i<field.MX();++i) {
        const auto cell=field.getIdx(i,j,k);const auto& n=state.normal(i,j,k);
        out<<sequence<<','<<operation<<','<<(rank?rank:"0")<<','<<field.globalDofId(i,j,k)<<','<<field.globalDofOwnerRank(i,j,k)<<','<<field.isGlobalDofOwner(i,j,k)<<','<<i<<','<<j<<','<<k<<','<<field.X(i,j,k)<<','<<field.Y(i,j,k)<<','<<field.Z(i,j,k)<<','
            <<(i>=ng && i<ng+field.NX() && j>=ng && j<ng+field.NY() && k>=ng && k<ng+field.NZ())<<','<<field.isCommunicationHalo(i,j,k)<<','<<field.isSolverBoundaryPoint(i,j,k)<<','<<field.CellFlag(i,j,k)<<','<<state.phi(i,j,k)<<','<<n.x<<','<<n.y<<','<<n.z<<','<<state.curvature(i,j,k)<<','<<state.isInterfaceCell(i,j,k)<<','<<state.density(i,j,k)<<','<<state.viscosity(i,j,k)<<','<<(rhs?rhs->at(cell):0)<<','<<(source?source->source(i,j,k,RU):0)<<','<<(source?source->source(i,j,k,RV):0)<<','<<(source?source->source(i,j,k,RW):0)<<','<<(source?source->source(i,j,k,E):0)<<'\n';
    }
}
}
