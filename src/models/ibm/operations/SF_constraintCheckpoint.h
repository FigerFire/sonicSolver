#pragma once

/// @file SF_constraintCheckpoint.h
/// @brief 可选的只读 IBM 首差异取证；不持有 state，不改变求解顺序。
#include "operations/SF_fieldOps.h"
#include "core/interfaces/SF_executionRuntime.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <unistd.h>

namespace SF::IBM::Checkpoint {
inline bool active(double time) {
    const char* directory=std::getenv("SF_IBM_CHECKPOINT_DIR");
    const char* limit=std::getenv("SF_IBM_CHECKPOINT_TIME");
    if(!directory) return false;
    const double end=limit?std::stod(limit):.001;
    if(!std::isfinite(end)||end<0) throw std::runtime_error("Invalid SF_IBM_CHECKPOINT_TIME.");
    const char* startValue=std::getenv("SF_IBM_CHECKPOINT_START_TIME");
    const double start=startValue?std::stod(startValue):0.;
    if (!std::isfinite(start) || start<0. || start>end)
        throw std::runtime_error("Invalid SF_IBM_CHECKPOINT_START_TIME.");
    return time>=start && time<=end;
}
inline std::ofstream output(const char* label,double time) {
    static std::size_t sequence=0;
    const std::string path=std::string(std::getenv("SF_IBM_CHECKPOINT_DIR"))+"/ibm-"
        +std::to_string(getpid())+"-"+std::to_string(sequence++)+"-"+label+".csv";
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Cannot write IBM checkpoint: "+path);
    out<<std::setprecision(17)<<"time,"<<time<<'\n';
    return out;
}
inline void state(const char* label,const Field& field,double time) {
    if(!active(time)) return;
    auto out=output(label,time);out<<"x,y,z,id,owner,flag,jac,mask";
    for(int v=0;v<field.NVar();++v)out<<",Q"<<v;out<<'\n';
    const int ng=field.NG();
    for(int k=ng;k<ng+field.NZ();++k)for(int j=ng;j<ng+field.NY();++j)for(int i=ng;i<ng+field.NX();++i) {
        out<<field.X(i,j,k)<<','<<field.Y(i,j,k)<<','<<field.Z(i,j,k)<<','<<field.globalDofId(i,j,k)<<','
           <<field.isGlobalDofOwner(i,j,k)<<','<<field.CellFlag(i,j,k)<<','<<field.Jac(i,j,k)<<','<<field.isSolverBoundaryPoint(i,j,k);
        for(int v=0;v<field.NVar();++v)out<<','<<field(i,j,k,v);out<<'\n';
    }
}
inline void markers(const char* label,const FDM::ImmersedSurfaceSystem& surface,
                    const std::vector<Vector3>& values,double time) {
    if(!active(time))return;
    if(values.size()!=surface.points.size())throw std::runtime_error("IBM checkpoint marker count mismatch.");
    auto out=output(label,time);out<<"marker,x,y,z,measure,UbX,UbY,UbZ,xValue,yValue,zValue\n";
    for(std::size_t n=0;n<values.size();++n) {
        const auto& point=surface.points[n];const auto& value=values[n];
        out<<point.globalConstraintId.value()<<','<<point.position.x<<','<<point.position.y<<','<<point.position.z<<','<<point.measure<<','
            <<point.prescribedVelocity.x<<','<<point.prescribedVelocity.y<<','<<point.prescribedVelocity.z<<','<<value.x<<','<<value.y<<','<<value.z<<'\n';
    }
}
inline void body(const Field& field,const FDM::ImmersedBodySystem& system,double time) {
    if(!active(time))return;
    auto out=output("body-points",time);out<<"x,y,z,id,owner,volume,rho\n";
    for(const auto& point:system.points) {
        const auto index=FieldOps::index3(field,point.localCell);
        out<<point.position.x<<','<<point.position.y<<','<<point.position.z<<','<<field.globalDofId(index.i,index.j,index.k)<<','
            <<field.isGlobalDofOwner(index.i,index.j,index.k)<<','<<point.dualVolume<<','<<FieldOps::density(field,index)<<'\n';
    }
}
inline void diagonal(const char* label,const FDM::ImmersedSurfaceSystem& surface,
                     const std::vector<double>& values,double time) {
    if(!active(time))return;
    std::vector<Vector3> columns;for(double value:values)columns.push_back({value,0.,0.});
    markers(label,surface,columns,time);
}
inline void result(const FDM::ImmersedConstraintResult& value,double time,double dt) {
    if(!active(time))return;
    auto out=output("constraint-result",time);
    out<<"cells,residual,stationarity,forceX,forceY,forceZ,torqueX,torqueY,torqueZ,power,dt\n"
       <<value.constrainedCells<<','<<value.maximumVelocityResidual<<','<<value.maximumStationarityResidual<<','
       <<value.forceOnBody.x<<','<<value.forceOnBody.y<<','<<value.forceOnBody.z<<','
       <<value.torqueOnBody.x<<','<<value.torqueOnBody.y<<','<<value.torqueOnBody.z<<','<<value.fluidMechanicalPower<<','<<dt<<'\n';
}
inline void force(const char* label,const Field& field,const std::vector<Vector3>& values,double time) {
    if(!active(time))return;
    if(values.size()!=static_cast<std::size_t>(field.TotalSize()))
        throw std::runtime_error("IBM force checkpoint storage mismatch.");
    auto out=output(label,time);out<<"x,y,z,owner,volume,Fx,Fy,Fz\n";
    const int ng=field.NG();
    for(int k=ng;k<ng+field.NZ();++k)for(int j=ng;j<ng+field.NY();++j)for(int i=ng;i<ng+field.NX();++i) {
        const auto& value=values[static_cast<std::size_t>(field.getIdx(i,j,k))];
        out<<field.X(i,j,k)<<','<<field.Y(i,j,k)<<','<<field.Z(i,j,k)<<','<<field.isGlobalDofOwner(i,j,k)<<','
           <<FieldOps::volume(field,{i,j,k})<<','<<value.x<<','<<value.y<<','<<value.z<<'\n';
    }
}
inline void graph(const Field& field,const FDM::ImmersedSurfaceSystem& surface,
                  FDM::IExecutionRuntime* runtime,double time) {
    if(!active(time))return;
    auto out=output("edges",time);out<<"marker,markerOwner,x,y,z,eulerian,weight,volume,rho\n";
    std::vector<double> unity(surface.points.size(),0.0),ju(surface.points.size(),0.0);
    double eulerianPair=0.0,markerPair=0.0;
    for(std::size_t n=0;n<surface.points.size();++n) {
        const auto& point=surface.points[n];
        // Independent deterministic fields, with physical M_L and M_E weights.
        const double lambda=1+.125*static_cast<double>(n);
        for(const auto& edge:point.interpolation) {
            const auto cell=FieldOps::index3(field,edge.cell);
            const double volume=FieldOps::volume(field,cell);
            const double testVelocity=1+.1*field.X(cell.i,cell.j,cell.k)+.2*field.Y(cell.i,cell.j,cell.k);
            unity[n]+=edge.value;ju[n]+=edge.value*testVelocity;
            const double spreadDensity=edge.value*point.measure*lambda/volume;
            eulerianPair+=volume*testVelocity*spreadDensity;
            out<<point.globalConstraintId.value()<<','<<(!runtime||runtime->ownsCanonicalEntity(point.globalConstraintId.value()))<<','
               <<field.X(cell.i,cell.j,cell.k)<<','<<field.Y(cell.i,cell.j,cell.k)<<','<<field.Z(cell.i,cell.j,cell.k)<<','
               <<edge.globalEulerianDofId<<','<<edge.value<<','<<volume<<','<<FieldOps::density(field,cell)<<'\n';
        }
    }
    if(runtime){runtime->globalSum(unity);runtime->globalSum(ju);}
    double unityError=0.;
    for(std::size_t n=0;n<surface.points.size();++n) {
        unityError=std::max(unityError,std::abs(unity[n]-1));
        if(!runtime||runtime->ownsCanonicalEntity(surface.points[n].globalConstraintId.value()))
            markerPair+=surface.points[n].measure*(1+.125*static_cast<double>(n))*ju[n];
    }
    std::vector<double> pairs{markerPair,eulerianPair};if(runtime)runtime->globalSum(pairs);
    auto metrics=output("weighted-adjoint",time);
    metrics<<"unityError,markerPair,eulerianPair,absoluteError\n"<<unityError<<','<<pairs[0]<<','<<pairs[1]<<','<<std::abs(pairs[0]-pairs[1])<<'\n';
}
}
