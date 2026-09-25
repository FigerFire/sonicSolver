/// @file SF_pressureOperators.cpp
/// @brief Serial structured constant-density pressure operations.

#include "solver/algorithm/pressureBased/SF_pressureOperators.h"

#include "methods/numerics/structured/SF_canonicalFace.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"
#include "solver/discretization/pressure/SF_rhieChow.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::PressureBased {
namespace {

std::array<double,3> components(const Vector3& value) {
    return {value.x,value.y,value.z};
}

int inwardCell(const Field& field, int cell) {
    int i=0,j=0,k=0;
    field.getIJK(cell,i,j,k);
    const int ng=field.NG();
    if (i==ng) ++i;
    else if (i==ng+field.NX()-1) --i;
    else if (j==ng) ++j;
    else if (j==ng+field.NY()-1) --j;
    else if (k==ng) ++k;
    else if (k==ng+field.NZ()-1) --k;
    return field.getIdx(i,j,k);
}

} // namespace

PressureOperators::PressureOperators(
        const FDM::SolverConfig& config,
        const System::CompiledNumericalSystem& numerics)
    : config_(config), numerics_(numerics),
      pressureSolver_(config.pressure.linear.pressure) {}

void PressureOperators::bind(const System::StateRealization& state) {
    const auto& u=state.at("U");
    const auto& p=state.at("p");
    const auto& rho=state.at("rho");
    if (!u.descriptor || !p.descriptor || !rho.descriptor
        || u.fields.size()!=1 || p.fields.size()!=1
        || !u.fields[0] || !p.fields[0]
        || !rho.descriptor->constantValue
        || !std::isfinite(*rho.descriptor->constantValue)
        || *rho.descriptor->constantValue<=0.0
        || u.descriptor->components!=3 || p.descriptor->components!=1) {
        throw std::runtime_error(
            "Pressure operations require one realized U[3], p[1], "
            "and positive derived rhoConst.");
    }
    velocity_=u.fields[0];
    pressure_=p.fields[0];
    geometry_=velocity_->geometry;
    density_=*rho.descriptor->constantValue;
    if (!geometry_ || pressure_->geometry!=geometry_)
        throw std::runtime_error("U/p pressure realization geometry differs.");
    if (numerics_.pressureFaceCoupling
        !=System::PressureFaceCoupling::RhieChow)
        throw std::runtime_error(
            "Constant-density pressure requires compiled RhieChow face coupling.");
    if (velocity_->blockId!=pressure_->blockId)
        throw std::runtime_error("U/p pressure realization block differs.");
    const auto convection=std::find_if(
        numerics_.terms.begin(),numerics_.terms.end(),
        [](const System::BoundTerm& term) {
            return term.equationId=="E_MOMENTUM"
                && term.kind==Equation::TermKind::Divergence;
        });
    if (convection==numerics_.terms.end()
        || convection->recipe.id()!=FDM::TermRecipeId::PrimitiveUpwind1
        || convection->primary!="momentumFlux") {
        throw std::runtime_error(
            "Constant-density Momentum requires compiled "
            "div(momentumFlux) -> primitiveUpwind1; no implicit "
            "convection replacement is available.");
    }
    for (const auto& term:numerics_.terms) {
        if (term.equationId!="E_MOMENTUM") continue;
        if (term.kind==Equation::TermKind::Diffusion) {
            if (term.recipe.id()!=FDM::TermRecipeId::Central2Explicit
                || term.primary!="U" || term.secondary!="nu") {
                throw std::runtime_error(
                    "Constant-density Momentum diffusion requires "
                    "compiled diffusion(nu,U) -> central2Explicit.");
            }
        } else if (term.kind==Equation::TermKind::Source) {
            throw std::runtime_error(
                "No constant-density numerical source provider is bound "
                "for Momentum term '"+term.primary+"'.");
        }
    }
    const int total=geometry_->TotalSize();
    cells_.clear();
    rowOfCell_.assign((size_t)total,-1);
    boundaryCell_.assign((size_t)total,0);
    pressureDirichlet_.assign((size_t)total,0);
    for (const auto& bc:config_.boundaries.velocity) {
        if (bc.type==EMPTY) continue;
        if (bc.type!=FIXED_VALUE && bc.type!=ZERO_GRADIENT)
            throw std::runtime_error(
                "Constant-density U boundary supports fixedValue/zeroGradient only.");
        for (int cell:geometry_->getSet(bc.name)) {
            if (cell<0 || cell>=total) {
                throw std::runtime_error("Velocity boundary cell is outside Field.");
            }
            boundaryCell_[(size_t)cell]=1;
        }
    }
    for (const auto& bc:config_.boundaries.energyFromPressure) {
        if (bc.type!=EMPTY && bc.type!=FIXED_VALUE
            && bc.type!=ZERO_GRADIENT)
            throw std::runtime_error(
                "Constant-density p boundary supports fixedValue/zeroGradient only.");
        if (bc.type!=FIXED_VALUE) continue;
        for (int cell:geometry_->getSet(bc.name)) {
            if (cell<0 || cell>=total)
                throw std::runtime_error("Pressure boundary cell is outside Field.");
            pressureDirichlet_[(size_t)cell]=1;
        }
    }
    const int ng=geometry_->NG();
    for (int k=ng;k<ng+geometry_->NZ();++k)
        for (int j=ng;j<ng+geometry_->NY();++j)
            for (int i=ng;i<ng+geometry_->NX();++i) {
                const int cell=geometry_->getIdx(i,j,k);
                if (boundaryCell_[(size_t)cell]
                    || geometry_->CellFlag(i,j,k)!=FLUID_CELL) continue;
                rowOfCell_[(size_t)cell]=static_cast<int>(cells_.size());
                cells_.push_back(cell);
            }
    if (cells_.empty())
        throw std::runtime_error("Pressure equation has no interior U/p cells.");
    for (int cell:cells_) {
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign:{-1,1}) {
                const int lower=sign>0?cell:adjacent(cell,axis,-1);
                const auto area=faceArea(axis,lower);
                int i=0,j=0,k=0, ni=0,nj=0,nk=0;
                geometry_->getIJK(lower,i,j,k);
                geometry_->getIJK(adjacent(lower,axis,1),ni,nj,nk);
                const std::array<double,3> delta{
                    geometry_->X(ni,nj,nk)-geometry_->X(i,j,k),
                    geometry_->Y(ni,nj,nk)-geometry_->Y(i,j,k),
                    geometry_->Z(ni,nj,nk)-geometry_->Z(i,j,k)};
                double dot=0.0,area2=0.0,delta2=0.0;
                for (int c=0;c<3;++c) {
                    dot+=area[(size_t)c]*delta[(size_t)c];
                    area2+=area[(size_t)c]*area[(size_t)c];
                    delta2+=delta[(size_t)c]*delta[(size_t)c];
                }
                if (!(area2>0.0) || !(delta2>0.0)
                    || std::abs(dot)/std::sqrt(area2*delta2)<1.0-1.0e-8) {
                    throw std::runtime_error(
                        "Constant-density pressure response currently requires "
                        "an orthogonal structured grid.");
                }
            }
        }
    }
    if (config_.pressure.reference.referenceCell<0
        || config_.pressure.reference.referenceCell
            >=static_cast<int>(cells_.size())) {
        throw std::runtime_error("Pressure referenceCell is outside pressure rows.");
    }
    hByA_.assign((size_t)total*3,0.0);
    rAU_.assign((size_t)total,0.0);
    correction_.assign((size_t)total,0.0);
    predictedFlux_.assign((size_t)total*3,0.0);
    faceResponse_.assign((size_t)total*3,0.0);
    correctedFlux_.assign((size_t)total*3,0.0);
}

int PressureOperators::adjacent(int cell,int axis,int sign) const {
    int i=0,j=0,k=0;
    geometry_->getIJK(cell,i,j,k);
    if (axis==0) i+=sign;
    else if (axis==1) j+=sign;
    else k+=sign;
    if (i<0 || i>=geometry_->MX() || j<0 || j>=geometry_->MY()
        || k<0 || k>=geometry_->MZ())
        throw std::runtime_error("Pressure stencil exceeds allocated halo.");
    return geometry_->getIdx(i,j,k);
}

bool PressureOperators::solved(int cell) const {
    return cell>=0 && cell<geometry_->TotalSize()
        && rowOfCell_[(size_t)cell]>=0;
}

std::size_t PressureOperators::faceIndex(int axis,int lowerCell) const {
    return (size_t)axis*(size_t)geometry_->TotalSize()+(size_t)lowerCell;
}

std::array<double,3> PressureOperators::faceArea(
        int axis,int lowerCell) const {
    int i=0,j=0,k=0;
    geometry_->getIJK(lowerCell,i,j,k);
    return Numerics::CanonicalFace::geometry(*geometry_,axis,i,j,k).cofactor;
}

double PressureOperators::spacing(int axis,int lowerCell) const {
    const int upperCell=adjacent(lowerCell,axis,1);
    int i=0,j=0,k=0, ni=0,nj=0,nk=0;
    geometry_->getIJK(lowerCell,i,j,k);
    geometry_->getIJK(upperCell,ni,nj,nk);
    const double dx=geometry_->X(ni,nj,nk)-geometry_->X(i,j,k);
    const double dy=geometry_->Y(ni,nj,nk)-geometry_->Y(i,j,k);
    const double dz=geometry_->Z(ni,nj,nk)-geometry_->Z(i,j,k);
    const double h=std::sqrt(dx*dx+dy*dy+dz*dz);
    if (!std::isfinite(h) || h<=0.0)
        throw std::runtime_error("Pressure operator found degenerate face spacing.");
    return h;
}

void PressureOperators::extendHalos(
        State::DistributedFieldView& view,int offset,int components) {
    const int ng=geometry_->NG();
    for (int k=0;k<geometry_->MZ();++k)
        for (int j=0;j<geometry_->MY();++j)
            for (int i=0;i<geometry_->MX();++i) {
                const int ci=std::clamp(i,ng,ng+geometry_->NX()-1);
                const int cj=std::clamp(j,ng,ng+geometry_->NY()-1);
                const int ck=std::clamp(k,ng,ng+geometry_->NZ()-1);
                if (i==ci && j==cj && k==ck) continue;
                const int cell=geometry_->getIdx(i,j,k);
                const int donor=geometry_->getIdx(ci,cj,ck);
                for (int c=0;c<components;++c)
                    view.write(cell,offset+c,view.read(donor,offset+c));
            }
}

void PressureOperators::applyVelocityBoundary() {
    const int offset=0;
    for (BCType type:{ZERO_GRADIENT,FIXED_VALUE}) {
        for (const auto& bc:config_.boundaries.velocity) {
            if (bc.type==EMPTY || bc.type!=type) continue;
            for (int cell:geometry_->getSet(bc.name)) {
                const int donor=inwardCell(*geometry_,cell);
                const auto value=components(bc.value);
                for (int c=0;c<3;++c) velocity_->write(
                    cell,offset+c,bc.type==FIXED_VALUE
                        ? value[(size_t)c]
                        : velocity_->read(donor,offset+c));
            }
        }
    }
    extendHalos(*velocity_,0,3);
}

void PressureOperators::applyPressureBoundary() {
    for (BCType type:{ZERO_GRADIENT,FIXED_VALUE}) {
        for (const auto& bc:config_.boundaries.energyFromPressure) {
            if (bc.type==EMPTY || bc.type!=type) continue;
            for (int cell:geometry_->getSet(bc.name)) {
                const int donor=inwardCell(*geometry_,cell);
                pressure_->write(cell,0,bc.type==FIXED_VALUE
                    ? bc.value : pressure_->read(donor,0));
            }
        }
    }
    extendHalos(*pressure_,0,1);
}

std::array<double,3> PressureOperators::gradient(
        const std::vector<double>& values,int cell) const {
    int i=0,j=0,k=0;
    geometry_->getIJK(cell,i,j,k);
    const auto g=CENTRAL2::gradSampled(*geometry_,[&](int ii,int jj,int kk) {
        return values.at((size_t)geometry_->getIdx(ii,jj,kk));
    },i,j,k);
    return {g.x,g.y,g.z};
}

std::array<double,3> PressureOperators::pressureGradient(int cell) const {
    int i=0,j=0,k=0;
    geometry_->getIJK(cell,i,j,k);
    const auto g=CENTRAL2::gradSampled(*geometry_,[&](int ii,int jj,int kk) {
        return pressure_->read(geometry_->getIdx(ii,jj,kk),0);
    },i,j,k);
    return {g.x,g.y,g.z};
}

double PressureOperators::faceFlux(
        const std::vector<double>& flux,int axis,int lowerCell) const {
    return flux.at(faceIndex(axis,lowerCell));
}

double PressureOperators::divergence(
        const std::vector<double>& flux,int cell) const {
    double net=0.0;
    for (int axis=0;axis<3;++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        net+=faceFlux(flux,axis,cell)
            -faceFlux(flux,axis,adjacent(cell,axis,-1));
    }
    int i=0,j=0,k=0;
    geometry_->getIJK(cell,i,j,k);
    return net*geometry_->Jac(i,j,k);
}

double PressureOperators::maxDivergence(
        const std::vector<double>& flux) const {
    double maximum=0.0;
    for (int cell:cells_)
        maximum=std::max(maximum,std::abs(divergence(flux,cell)));
    return maximum;
}

double PressureOperators::prepare(double maximumTimeStep,double) {
    if (!velocity_ || !pressure_ || !std::isfinite(maximumTimeStep)
        || maximumTimeStep<=0.0)
        throw std::runtime_error("Pressure prepare requires bound state and positive dt limit.");
    applyVelocityBoundary();
    applyPressureBoundary();
    const double nu=config_.numerics.dynamicViscosity/density_;
    double limit=std::min(maximumTimeStep,numerics_.dt.maxDeltaT);
    for (int cell:cells_) {
        double speed2=0.0;
        for (int c=0;c<3;++c) {
            const double u=velocity_->read(cell,c);
            if (!std::isfinite(u))
                throw std::runtime_error("Constant-density U is non-finite.");
            speed2+=u*u;
        }
        const double speed=std::sqrt(speed2);
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            const double h=spacing(axis,cell);
            if (speed>0.0)
                limit=std::min(limit,numerics_.dt.cfl*h/speed);
            if (nu>0.0)
                limit=std::min(limit,
                    numerics_.dt.cfl*h*h/(2.0*Math::activeDimensionCount()*nu));
        }
    }
    if (!std::isfinite(limit) || limit<=0.0)
        throw std::runtime_error("Constant-density dt is invalid.");
    dt_=limit;
    summary_={};
    return dt_;
}

void PressureOperators::assembleMomentum() {
    if (!(dt_>0.0)) throw std::runtime_error("Momentum assemble requires dt.");
    const int total=geometry_->TotalSize();
    const double nu=config_.numerics.dynamicViscosity/density_;
    const bool diffusion=std::any_of(
        numerics_.terms.begin(),numerics_.terms.end(),
        [](const System::BoundTerm& term) {
            return term.equationId=="E_MOMENTUM"
                && term.kind==Equation::TermKind::Diffusion;
        });
    std::vector<double> oldFlux((size_t)total*3,0.0);
    for (int cell:cells_) {
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign:{-1,1}) {
                const int lower=sign>0?cell:adjacent(cell,axis,-1);
                const int upper=adjacent(lower,axis,1);
                const auto area=faceArea(axis,lower);
                double flux=0.0;
                for (int c=0;c<3;++c)
                    flux+=0.5*(velocity_->read(lower,c)
                               +velocity_->read(upper,c))*area[(size_t)c];
                oldFlux[faceIndex(axis,lower)]=flux;
            }
        }
    }
    for (int cell:cells_) {
        int i=0,j=0,k=0;
        geometry_->getIJK(cell,i,j,k);
        const double inverseVolume=geometry_->Jac(i,j,k);
        if (!std::isfinite(inverseVolume) || inverseVolume<=0.0)
            throw std::runtime_error("Momentum cell has invalid volume.");
        rAU_[(size_t)cell]=dt_/density_;
        for (int c=0;c<3;++c) {
            double advective=0.0;
            double laplacian=0.0;
            for (int axis=0;axis<3;++axis) {
                if (!Math::isDirectionActiveIndex(axis)) continue;
                for (int sign:{-1,1}) {
                    const int lower=sign>0?cell:adjacent(cell,axis,-1);
                    const int upper=adjacent(lower,axis,1);
                    const double flux=oldFlux[faceIndex(axis,lower)];
                    const double upwind=velocity_->read(
                        flux>=0.0?lower:upper,c);
                    advective+=sign*flux*upwind;
                    if (diffusion) {
                        const auto area=faceArea(axis,lower);
                        const double magnitude=std::sqrt(
                            area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
                        const int neighbour=adjacent(cell,axis,sign);
                        laplacian+=magnitude
                            *(velocity_->read(neighbour,c)-velocity_->read(cell,c))
                            /spacing(axis,lower);
                    }
                }
            }
            hByA_[(size_t)c*(size_t)total+(size_t)cell]=velocity_->read(cell,c)
                +dt_*inverseVolume*(-advective+nu*laplacian);
        }
    }
}

void PressureOperators::solveMomentum() {
    const int total=geometry_->TotalSize();
    for (int cell:cells_) {
        const auto grad=pressureGradient(cell);
        for (int c=0;c<3;++c) {
            const double predicted=hByA_[(size_t)c*(size_t)total+(size_t)cell]
                -rAU_[(size_t)cell]*grad[(size_t)c];
            if (!std::isfinite(predicted))
                throw std::runtime_error("Momentum predictor produced non-finite U.");
            velocity_->write(cell,c,predicted);
        }
    }
    applyVelocityBoundary();
    std::fill(predictedFlux_.begin(),predictedFlux_.end(),0.0);
    std::fill(faceResponse_.begin(),faceResponse_.end(),0.0);
    for (int cell:cells_) {
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign:{-1,1}) {
                const int lower=sign>0?cell:adjacent(cell,axis,-1);
                const int upper=adjacent(lower,axis,1);
                const auto area=faceArea(axis,lower);
                const auto index=faceIndex(axis,lower);
                const bool both=solved(lower)&&solved(upper);
                double flux=0.0;
                const int boundary=both?-1:(solved(lower)?upper:lower);
                if (both || (boundary>=0
                             && pressureDirichlet_[(size_t)boundary])) {
                    const auto lowerGrad=pressureGradient(lower);
                    const auto upperGrad=pressureGradient(upper);
                    std::array<double,3> lowerVelocity{},upperVelocity{};
                    for (int c=0;c<3;++c) {
                        lowerVelocity[(size_t)c]=velocity_->read(lower,c);
                        upperVelocity[(size_t)c]=velocity_->read(upper,c);
                    }
                    const double response=both ? 0.5
                        *(rAU_[(size_t)lower]+rAU_[(size_t)upper])
                        : rAU_[(size_t)(solved(lower)?lower:upper)];
                    const auto face=RhieChow::predict(
                        lowerVelocity,upperVelocity,lowerGrad,upperGrad,area,
                        pressure_->read(lower,0),pressure_->read(upper,0),
                        response,spacing(axis,lower));
                    flux=face.flux;
                    faceResponse_[index]=face.correctionCoefficient;
                } else {
                    for (int c=0;c<3;++c)
                        flux+=velocity_->read(boundary,c)*area[(size_t)c];
                }
                predictedFlux_[index]=flux;
            }
        }
    }
    correctedFlux_=predictedFlux_;
}

void PressureOperators::preparePressureBoundary() {
    applyVelocityBoundary();
    applyPressureBoundary();
}

void PressureOperators::assemblePressure() {
    summary_.maxDivergenceBefore=maxDivergence(correctedFlux_);
    const bool hasDirichlet=std::any_of(
        pressureDirichlet_.begin(),pressureDirichlet_.end(),
        [](unsigned char value) { return value!=0; });
    pressureMatrix_={};
    pressureMatrix_.firstRow=0;
    pressureMatrix_.lastRow=static_cast<std::int64_t>(cells_.size())-1;
    pressureMatrix_.globalSize=static_cast<std::int64_t>(cells_.size());
    pressureMatrix_.rows.reserve(cells_.size());
    pressureMatrix_.rhs.reserve(cells_.size());
    pressureMatrix_.initialGuess.assign(cells_.size(),0.0);
    for (size_t index=0;index<cells_.size();++index) {
        const int cell=cells_[index];
        LinearAlgebra::SparseRow row;
        row.globalRow=static_cast<std::int64_t>(index);
        if (!hasDirichlet
            && static_cast<int>(index)==config_.pressure.reference.referenceCell) {
            row.columns={row.globalRow};
            row.values={1.0};
            pressureMatrix_.rows.push_back(std::move(row));
            pressureMatrix_.rhs.push_back(0.0);
            continue;
        }
        double diagonal=0.0;
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign:{-1,1}) {
                const int neighbour=adjacent(cell,axis,sign);
                const int lower=sign>0?cell:neighbour;
                const double coefficient=faceResponse_[faceIndex(axis,lower)];
                if (!solved(neighbour) && coefficient==0.0) continue;
                if (!std::isfinite(coefficient) || coefficient<=0.0)
                    throw std::runtime_error("Pressure face response is invalid.");
                diagonal+=coefficient;
                if (solved(neighbour)) {
                    row.columns.push_back(rowOfCell_[(size_t)neighbour]);
                    row.values.push_back(-coefficient);
                }
            }
        }
        if (!(diagonal>0.0))
            throw std::runtime_error("Pressure row has no connected faces.");
        row.columns.push_back(row.globalRow);
        row.values.push_back(diagonal);
        double net=0.0;
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            net+=faceFlux(correctedFlux_,axis,cell)
                -faceFlux(correctedFlux_,axis,adjacent(cell,axis,-1));
        }
        pressureMatrix_.rows.push_back(std::move(row));
        pressureMatrix_.rhs.push_back(-net);
    }
}

void PressureOperators::solvePressure() {
    pressureResult_=pressureSolver_.solve(pressureMatrix_);
    if (!pressureResult_.converged
        || pressureResult_.solution.size()!=cells_.size())
        throw std::runtime_error("Constant-density pressure solve did not converge.");
    std::fill(correction_.begin(),correction_.end(),0.0);
    for (size_t row=0;row<cells_.size();++row)
        correction_[(size_t)cells_[row]]=pressureResult_.solution[row];
    // Zero-gradient correction at velocity boundaries preserves the
    // pressure-null-space condition while supplying a centered cell gradient.
    for (size_t cell=0;cell<boundaryCell_.size();++cell) {
        if (!boundaryCell_[cell]) continue;
        correction_[cell]=pressureDirichlet_[cell]
            ? 0.0 : correction_[(size_t)inwardCell(*geometry_,(int)cell)];
    }
    summary_.iterations=pressureResult_.iterations;
    summary_.relativeResidual=pressureResult_.relativeResidual;
}

void PressureOperators::preparePressureUpdate() {
    for (int cell:cells_) {
        const double updated=pressure_->read(cell,0)+correction_[(size_t)cell];
        if (!std::isfinite(updated))
            throw std::runtime_error("Pressure correction produced non-finite p.");
        pressure_->write(cell,0,updated);
    }
    const bool hasDirichlet=std::any_of(
        pressureDirichlet_.begin(),pressureDirichlet_.end(),
        [](unsigned char value) { return value!=0; });
    if (!hasDirichlet) {
        // The pressure correction removes the null space; its absolute gauge
        // is the configured reference pressure, independent of initial p.
        const int reference=cells_[(size_t)config_.pressure.reference.referenceCell];
        const double shift=config_.pressure.reference.referencePressure
            -pressure_->read(reference,0);
        for (int cell:cells_)
            pressure_->write(cell,0,pressure_->read(cell,0)+shift);
    }
    applyPressureBoundary();
}

void PressureOperators::correctVelocity() {
    for (int cell:cells_) {
        const auto grad=gradient(correction_,cell);
        for (int c=0;c<3;++c) {
            const double corrected=velocity_->read(cell,c)
                -rAU_[(size_t)cell]*grad[(size_t)c];
            if (!std::isfinite(corrected))
                throw std::runtime_error("Velocity correction produced non-finite U.");
            velocity_->write(cell,c,corrected);
        }
    }
    applyVelocityBoundary();
}

void PressureOperators::correctFlux() {
    for (int cell:cells_) {
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign:{-1,1}) {
                const int lower=sign>0?cell:adjacent(cell,axis,-1);
                const int upper=adjacent(lower,axis,1);
                if (sign<0 && solved(lower)) continue;
                const auto index=faceIndex(axis,lower);
                if (faceResponse_[index]==0.0) continue;
                correctedFlux_[index]=RhieChow::correct(
                    correctedFlux_[index],faceResponse_[index],
                    correction_[(size_t)lower],correction_[(size_t)upper]);
            }
        }
    }
    summary_.maxDivergenceAfter=maxDivergence(correctedFlux_);
}

PressureOperationSummary PressureOperators::commitCorrection() {
    if (!std::isfinite(summary_.maxDivergenceAfter))
        throw std::runtime_error("Corrected continuity defect is non-finite.");
    return summary_;
}

void PressureOperators::commitStep() {
    applyVelocityBoundary();
    applyPressureBoundary();
}

} // namespace SF::PressureBased
