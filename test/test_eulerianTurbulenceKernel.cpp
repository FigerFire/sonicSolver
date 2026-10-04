#include "solver/algorithm/eulerian/equation/SF_transport.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Independent old transport-kernel harness: no composer, HOW or numerical-provider routing.
int main() {
    SF::Field field;field.setup(5,5,5,1,5);
    for (int cell=0;cell<field.TotalSize();++cell) {
        int i,j,k;field.getIJK(cell,i,j,k);
        field.X(i,j,k)=i;field.Y(i,j,k)=j;field.Z(i,j,k)=k;
        field.Jac(i,j,k)=1;field.XiX(i,j,k)=1;field.EtY(i,j,k)=1;field.ZeZ(i,j,k)=1;
    }
    const int center=field.getIdx(3,3,3);
    SF::LinearAlgebra::DistributedRowMap map;map.global.assign(field.TotalSize(),-1);
    map.localCells={center};map.first=map.last=0;map.total=7;map.global[center]=0;
    int column=1;
    for (int axis=0;axis<3;++axis) for (int sign:{-1,1}) {
        int di=0,dj=0,dk=0;SF::EulerianEulerian::Ops::offset(axis,sign,di,dj,dk);
        map.global[field.getIdx(3+di,3+dj,3+dk)]=column++;
    }
    SF::ScalarField phi,old,mass,diffusion,source,sink,diagonal;
    phi.setupLike(field,"k",.06);old.setupLike(field,"previous.mass*k",.15);
    mass.setupLike(field,"phaseMass",2.5);diffusion.setupLike(field,"diffusion",.25);
    source.setupLike(field,"production",.8);sink.setupLike(field,"implicitSink",.6);
    diagonal.setupLike(field,"diagonal");
    std::vector<double> flux(3*field.TotalSize());
    const double value[]={2,-3,4};
    for (int axis=0;axis<3;++axis) for (int cell=0;cell<field.TotalSize();++cell) {
        int i,j,k;field.getIJK(cell,i,j,k);
        flux[SF::EulerianEulerian::PhaseFaceFlux::index(field,axis,i,j,k)]=value[axis];
    }
    const auto matrix=SF::EulerianEulerian::assembleTransportEquation(field,map,phi,old,mass,
        diffusion,source,flux,.02,false,0,diagonal,&sink);
    const auto& row=matrix.rows.at(0);
    double expected[]={136.1,-2.25,-.25,-.25,-3.25,-4.25,-.25};
    if (row.columns.size()!=7 || std::abs(matrix.rhs.at(0)-8.3)>1e-13)
        throw std::runtime_error("Phase-mass/old-time/source contribution changed.");
    for (std::size_t i=0;i<row.columns.size();++i)
        if (std::abs(row.values[i]-expected[row.columns[i]])>1e-13)
            throw std::runtime_error("Canonical upwind/diffusion/sink sign changed.");
    const auto explicitMatrix=SF::EulerianEulerian::assembleTransportEquation(field,map,phi,old,mass,
        diffusion,source,flux,.02,false,0,diagonal,nullptr);
    const auto findDiagonal=[](const auto& row) {
        for (std::size_t i=0;i<row.columns.size();++i) if (row.columns[i]==0) return row.values[i];
        throw std::runtime_error("Missing implicit diagonal.");
    };
    if (std::abs(findDiagonal(row)-findDiagonal(explicitMatrix.rows[0])-.6)>1e-13
        || explicitMatrix.rhs!=matrix.rhs || phi.values()[center]!=.06 || old.values()[center]!=.15)
        throw std::runtime_error("Sink became an explicit RHS or assembly committed primary/old-time state.");
    std::cout<<std::setprecision(17);
    for (std::size_t i=0;i<row.columns.size();++i) std::cout<<row.columns[i]<<':'<<row.values[i]<<'\n';
    std::cout<<"rhs:"<<matrix.rhs[0]<<" initial:"<<matrix.initialGuess[0]<<'\n';
}
