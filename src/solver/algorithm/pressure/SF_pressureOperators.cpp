/// @file SF_pressureOperators.cpp
/// @brief Structured constant-density pressure operations.

#include "solver/algorithm/pressure/SF_pressureOperators.h"
#include "solver/algorithm/pressure/SF_fixedTimeMath.h"

#include "methods/numerics/structured/SF_canonicalFace.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"
#include "solver/discretization/pressure/SF_rhieChow.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <stdexcept>

namespace SF::Pressure {
namespace {

std::array<double,3> components(const Vector3& value) {
    return {value.x,value.y,value.z};
}

int inwardCell(const Field& field, int cell, const std::array<bool,6>& sides) {
    int i=0,j=0,k=0;
    field.getIJK(cell,i,j,k);
    const int ng=field.NG();
    if (i==ng && sides[0]) ++i;
    else if (i==ng+field.NX()-1 && sides[1]) --i;
    else if (j==ng && sides[2]) ++j;
    else if (j==ng+field.NY()-1 && sides[3]) --j;
    else if (k==ng && sides[4]) ++k;
    else if (k==ng+field.NZ()-1 && sides[5]) --k;
    return field.getIdx(i,j,k);
}

} // namespace

PressureOperators::PressureOperators(
        PressureOperatorConfig config,
        const System::CompiledNumericalSystem& numerics)
    : config_(std::move(config)), numerics_(numerics) {}

void PressureOperators::bind(System::StateRealization& state,
                             FDM::IExecutionRuntime& runtime) {
    runtime_=&runtime;
    stateViews_=&state;
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
    physicalVelocity_=velocity_;
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
    int convectionCount=0;
    int diffusionCount=0;
    momentumTerms_.clear();
    for (const auto& term:numerics_.operators) {
        if (term.output!="U" || term.equationMethod!="PressureMomentum") continue;
        if (term.operatorName=="div") {
            ++convectionCount;
            if (!term.recipe
                || term.recipe->id()!=FDM::TermRecipeId::PrimitiveUpwind1
                || term.primary!="momentumFlux"
                || !term.primitiveSpatial) {
                throw std::runtime_error(
                    "Constant-density Momentum requires compiled "
                    "div(momentumFlux) -> primitiveUpwind1.");
            }
        } else if (term.operatorName=="diffusion") {
            ++diffusionCount;
            if (!term.recipe
                || term.recipe->id()!=FDM::TermRecipeId::Central2Explicit
                || term.primary!="U" || term.secondary!="nu"
                || !term.primitiveSpatial) {
                throw std::runtime_error(
                    "Constant-density Momentum diffusion requires "
                    "compiled diffusion(nu,U) -> central2Explicit.");
            }
        } else if (term.operatorName=="source") {
            if (term.side!=System::CompiledSpatialBinding::Side::Right
                || !term.primitiveSource) {
                throw std::runtime_error(
                    "Constant-density Momentum source '"+term.primary
                    +"' has no right-hand primitive-velocity provider.");
            }
        } else {
            throw std::runtime_error(
                "No constant-density numerical provider is bound for a "
                "compiled Momentum term.");
        }
        momentumTerms_.push_back(term);
    }
    if (convectionCount!=1 || diffusionCount>1)
        throw std::runtime_error(
            "Constant-density Momentum requires exactly one compiled "
            "primitiveUpwind1 convection and at most one central2 diffusion term.");
    const int total=geometry_->TotalSize();
    cells_.clear();
    rowOfCell_.assign((size_t)total,-1);
    boundaryCell_.assign((size_t)total,0);
    pressureDirichlet_.assign((size_t)total,0);
    for (const auto& bc:config_.velocityBoundary) {
        if (bc.type==EMPTY) continue;
        if (bc.type!=FIXED_VALUE && bc.type!=ZERO_GRADIENT)
            throw std::runtime_error(
                "Constant-density U boundary supports fixedValue/zeroGradient only.");
        const auto& points=geometry_->getSet(bc.name);
        for (int cell:points) {
            if (cell<0 || cell>=total) {
                throw std::runtime_error("Velocity boundary cell is outside Field.");
            }
            boundaryCell_[(size_t)cell]=1;
        }
    }
    // A named boundary may group several physical sides (e.g. cavity walls).
    // A processor plane has interior fluid points and must never choose a BC donor.
    for (int axis=0;axis<3;++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        for (int side=0;side<2;++side) {
            const int ng=geometry_->NG();
            const int extent=axis==0?geometry_->NX():axis==1?geometry_->NY():geometry_->NZ();
            const int plane=ng+(side?extent-1:0);
            bool physical=true;
            for (int k=ng;k<ng+geometry_->NZ();++k)
                for (int j=ng;j<ng+geometry_->NY();++j)
                    for (int i=ng;i<ng+geometry_->NX();++i)
                        if ((axis==0?i:axis==1?j:k)==plane)
                            physical=physical && boundaryCell_[(size_t)geometry_->getIdx(i,j,k)];
            physicalSides_[(size_t)(2*axis+side)]=physical;
        }
    }
    for (const auto& bc:config_.pressureBoundary) {
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
    hByA_.assign((size_t)total*3,0.0);
    rAU_.assign((size_t)total,0.0);
    correction_.assign((size_t)total,0.0);
    predictedFlux_.assign((size_t)total*3,0.0);
    faceResponse_.assign((size_t)total*3,0.0);
    correctedFlux_.assign((size_t)total*3,0.0);
    workingVelocity_.resize((size_t)total*3);
    workingVelocityView_=State::workspaceView("pressureMomentumWorkspace",physicalVelocity_->blockId,
        *geometry_,workingVelocity_,3,physicalVelocity_->haloDepth,physicalVelocity_->stages);
    correctionView_=State::workspaceView("pressureCorrection",pressure_->blockId,
        *geometry_,correction_,1,pressure_->haloDepth,pressure_->stages);
    fluxView_.name="pressureFaceFlux";fluxView_.blockId=physicalVelocity_->blockId;
    fluxView_.geometry=geometry_;fluxView_.location=State::FieldLocation::Face;
    fluxView_.exchange=State::ExchangeKind::CanonicalFaceFlux;fluxView_.components=1;
    fluxView_.read=[this](int face,int) { return correctedFlux_.at((size_t)face); };
    fluxView_.write=[this](int face,int,double value) { correctedFlux_.at((size_t)face)=value; };
    if (stateViews_->requestsView("U",System::StateViewKind::Working))
        stateViews_->bindView("U",System::StateViewKind::Working,workingVelocityView_);
    stateViews_->bindView("p",System::StateViewKind::Correction,correctionView_);
    stateViews_->bindView("phi",System::StateViewKind::Physical,fluxView_);
    initializeNumbering();
}

void PressureOperators::requireCollectively(bool valid,const char* reason) const {
    if (!runtime_ || runtime_->globalMinimum(valid?1.0:0.0)<0.5)
        throw std::runtime_error(reason);
}

void PressureOperators::synchronize(State::DistributedFieldView& view) {
    bool valid=true;
    for (int cell=0;cell<view.geometry->TotalSize();++cell)
        for (int component=0;component<view.components;++component)
            valid=std::isfinite(view.read(cell,component)) && valid;
    requireCollectively(valid,"Non-finite pressure state/workspace before halo or canonical COPY.");
    runtime_->synchronizeTransient({view});
}

void PressureOperators::synchronize(std::vector<double>& values,int count,
                                     State::ExchangeKind kind) {
    auto view=State::workspaceView("pressureWorkspace",velocity_->blockId,
        *geometry_,values,count,kind);
    synchronize(view);
}

LinearAlgebra::GlobalDofId PressureOperators::dof(int cell) const {
    return LinearAlgebra::GlobalDofId::make(LinearAlgebra::GlobalDofSpace::Pressure,
        entityOfCell_.at((size_t)cell));
}

void PressureOperators::initializeNumbering() {
    candidateCells_=cells_;
    const int total=geometry_->TotalSize();
    entityOfCell_.assign((size_t)total,-1);
    backendRow_.assign((size_t)total,-1);
    cells_.clear();
    bool valid=true;
    for (int cell:candidateCells_) {
        int i=0,j=0,k=0;
        geometry_->getIJK(cell,i,j,k);
        const int entity=runtime_->distributed()
            ? geometry_->globalDofId(i,j,k) : cell;
        valid=valid && entity>=0;
        entityOfCell_[(size_t)cell]=entity;
        if (!runtime_->distributed() || geometry_->isGlobalDofOwner(i,j,k))
            cells_.push_back(cell);
    }
    if (runtime_->globalMinimum(valid && !cells_.empty() ? 1.0 : 0.0)<0.5)
        throw std::runtime_error("Distributed pressure requires valid GlobalDof identities and nonempty owned rows on every execution partition (empty-row partitions Unsupported).");
    const auto range=runtime_->allocateDistributedIndices((std::int64_t)cells_.size());
    for (size_t row=0;row<cells_.size();++row)
        backendRow_[(size_t)cells_[row]]=range.first+(std::int64_t)row;
    runtime_->synchronizeIdentifiers(*geometry_,velocity_->blockId,entityOfCell_);
    runtime_->synchronizeIdentifiers(*geometry_,velocity_->blockId,backendRow_);
    std::vector<std::pair<LinearAlgebra::GlobalDofId,std::int64_t>> entries;
    for (int cell=0;cell<total;++cell)
        if (backendRow_[(size_t)cell]>=0)
            entries.push_back({dof(cell),backendRow_[(size_t)cell]});
    numbering_=std::make_unique<LinearAlgebra::StaticDistributedNumbering>(
        range.first,range.last,range.total,entries);
    pressureSolver_=std::make_unique<LinearAlgebra::DistributedLinearSystem>(
        config_.pressureLinear,*numbering_);
    hasDirichlet_=runtime_->globalMaximum(std::any_of(
        pressureDirichlet_.begin(),pressureDirichlet_.end(),
        [](unsigned char value) { return value!=0; }) ? 1.0 : 0.0)>0.5;
    if (config_.reference.referenceCell<0 || config_.reference.referenceCell>=range.total)
        throw std::runtime_error("Pressure referenceCell is outside global pressure rows.");
    // Global mesh identities preserve source point order; backend row numbers do not.
    std::int64_t low=std::numeric_limits<std::int64_t>::max(),high=-1;
    for (int cell:cells_) {
        low=std::min(low,entityOfCell_[(size_t)cell]);
        high=std::max(high,entityOfCell_[(size_t)cell]);
    }
    auto left=runtime_->globalMinimum(low);
    auto right=runtime_->globalMaximum(high);
    while (left<right) {
        const auto mid=left+(right-left)/2;
        std::int64_t count=0;
        for (int cell:cells_) count+=entityOfCell_[(size_t)cell]<=mid;
        if (runtime_->globalSum(count)>config_.reference.referenceCell) right=mid;
        else left=mid+1;
    }
    referenceEntity_=left;
}

void PressureOperators::restoreGauge() {
    if (hasDirichlet_) return;
    double reference=0.0;
    for (int cell:cells_)
        if (entityOfCell_[(size_t)cell]==referenceEntity_)
            reference=pressure_->read(cell,0);
    const double shift=config_.reference.referencePressure-runtime_->globalSum(reference);
    for (int cell:cells_)
        pressure_->write(cell,0,pressure_->read(cell,0)+shift);
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
        && (backendRow_.empty() ? rowOfCell_[(size_t)cell]>=0
                               : backendRow_[(size_t)cell]>=0);
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
                if (geometry_->isCommunicationHalo(i,j,k)) continue;
                const int cell=geometry_->getIdx(i,j,k);
                const int donor=geometry_->getIdx(ci,cj,ck);
                for (int c=0;c<components;++c)
                    view.write(cell,offset+c,view.read(donor,offset+c));
            }
}

void PressureOperators::applyVelocityBoundary() {
    const int offset=0;
    for (BCType type:{ZERO_GRADIENT,FIXED_VALUE}) {
        for (const auto& bc:config_.velocityBoundary) {
            if (bc.type==EMPTY || bc.type!=type) continue;
            for (int cell:geometry_->getSet(bc.name)) {
                const int donor=inwardCell(*geometry_,cell,physicalSides_);
                const auto value=components(bc.value);
                for (int c=0;c<3;++c) velocity_->write(
                    cell,offset+c,bc.type==FIXED_VALUE
                        ? value[(size_t)c]
                        : velocity_->read(donor,offset+c));
            }
        }
    }
    extendHalos(*velocity_,0,3);
    synchronize(*velocity_);
}

void PressureOperators::applyPressureBoundary() {
    for (BCType type:{ZERO_GRADIENT,FIXED_VALUE}) {
        for (const auto& bc:config_.pressureBoundary) {
            if (bc.type==EMPTY || bc.type!=type) continue;
            for (int cell:geometry_->getSet(bc.name)) {
                const int donor=inwardCell(*geometry_,cell,physicalSides_);
                pressure_->write(cell,0,bc.type==FIXED_VALUE
                    ? bc.value : pressure_->read(donor,0));
            }
        }
    }
    extendHalos(*pressure_,0,1);
    synchronize(*pressure_);
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
    bool valid=true;
    for (int cell:cells_) {
        const double value=divergence(flux,cell);
        valid=std::isfinite(value) && valid;
        maximum=std::max(maximum,std::abs(value));
    }
    requireCollectively(valid,"Non-finite pressure divergence before globalMaximum.");
    return runtime_->globalMaximum(maximum);
}

double PressureOperators::prepare(double maximumTimeStep) {
    requireCollectively(velocity_ && pressure_ && std::isfinite(maximumTimeStep)
        && maximumTimeStep>0.0,"Pressure prepare requires bound state and positive dt limit.");
    applyVelocityBoundary();
    applyPressureBoundary();
    const double nu=config_.dynamicViscosity/density_;
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
    dt_=runtime_->globalMinimum(limit);
    summary_={};
    fixedTimeActive_=false;
    iterationActive_=false;
    iterationIndex_=0;
    return dt_;
}

void PressureOperators::beginFixedTimeStep() {
    if (!(dt_>0.0) || !velocity_ || !pressure_)
        throw std::runtime_error("Fixed-time pressure step requires prepared state.");
    const int total=geometry_->TotalSize();
    baseVelocity_.resize((size_t)total*3);
    iterationVelocity_.resize((size_t)total*3);
    iterationPressure_.resize((size_t)total);
    iterationFlux_.resize(correctedFlux_.size());
    for (int cell=0;cell<geometry_->TotalSize();++cell)
        if (solved(cell)) rAU_[(size_t)cell]=dt_/density_;
    // The first F_k must match U_n/p_n at this timestep's pressure response,
    // rather than inheriting an uninitialized or previous-dt face buffer.
    reconstructFaceFlux(correctedFlux_);
    for (int cell=0;cell<total;++cell)
        for (int c=0;c<3;++c)
            baseVelocity_[(size_t)c*(size_t)total+(size_t)cell]
                =velocity_->read(cell,c);
    fixedTimeActive_=true;
    iterationIndex_=0;
    converged_=false;
    candidateContinuityDefect_=0.0;
}

void PressureOperators::beginIteration() {
    if (!fixedTimeActive_ || iterationActive_)
        throw std::runtime_error("Pressure iteration requires one fixed-time step.");
    const int total=geometry_->TotalSize();
    for (int cell=0;cell<total;++cell) {
        iterationPressure_[(size_t)cell]=pressure_->read(cell,0);
        for (int c=0;c<3;++c)
            iterationVelocity_[(size_t)c*(size_t)total+(size_t)cell]
                =velocity_->read(cell,c);
    }
    iterationFlux_=correctedFlux_;
    iterationActive_=true;
    ++iterationIndex_;
    converged_=false;
}

void PressureOperators::assembleMomentum() {
    if (!(dt_>0.0)) throw std::runtime_error("Momentum assemble requires dt.");
    const int total=geometry_->TotalSize();
    const double nu=config_.dynamicViscosity/density_;
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
    const Discretization::PressureMomentum::Stencil stencil{
        *velocity_,*geometry_,oldFlux};
    for (int cell:cells_) {
        int i=0,j=0,k=0;
        geometry_->getIJK(cell,i,j,k);
        const double inverseVolume=geometry_->Jac(i,j,k);
        if (!std::isfinite(inverseVolume) || inverseVolume<=0.0)
            throw std::runtime_error("Momentum cell has invalid volume.");
        rAU_[(size_t)cell]=dt_/density_;
        const Vector3 velocity(velocity_->read(cell,0),
                               velocity_->read(cell,1),
                               velocity_->read(cell,2));
        Vector3 sourceAcceleration;
        for (const auto& term:momentumTerms_)
            if (term.operatorName=="source")
                sourceAcceleration=sourceAcceleration+term.primitiveSource(
                    *geometry_,i,j,k,velocity);
        const auto source=components(sourceAcceleration);
        for (int c=0;c<3;++c) {
            double advective=0.0;
            double laplacian=0.0;
            for (const auto& term:momentumTerms_) {
                if (term.operatorName=="div")
                    advective=term.primitiveSpatial(stencil,cell,c);
                else if (term.operatorName=="diffusion")
                    laplacian=term.primitiveSpatial(stencil,cell,c);
            }
            const size_t index=(size_t)c*(size_t)total+(size_t)cell;
            if (fixedTimeActive_ && iterationIndex_>1) {
                hByA_[index]=fixedTimeCandidate(baseVelocity_[index],
                    dt_*inverseVolume*(-advective+nu*laplacian));
            } else {
                // At the first iterate U_k=U_n; retain the Phase 26 PISO
                // evaluation order exactly for the outer=1 degeneration.
                hByA_[index]=velocity_->read(cell,c)
                    +dt_*inverseVolume*(-advective+nu*laplacian);
            }
            hByA_[index]+=dt_*source[(size_t)c];
        }
    }
}

void PressureOperators::solveMomentum(System::TargetKind target) {
    if (target!=System::TargetKind::Physical && target!=System::TargetKind::Working)
        throw std::runtime_error("Momentum provider requires physical or working target storage.");
    const int total=geometry_->TotalSize();
    workingVelocity_.resize((size_t)total*3);
    for (int cell=0;cell<total;++cell)
        for (int c=0;c<3;++c)
            workingVelocity_[(size_t)c*total+cell]=velocity_->read(cell,c);
    workingVelocityView_=State::workspaceView("pressureMomentumWorkspace",physicalVelocity_->blockId,
        *geometry_,workingVelocity_,3,physicalVelocity_->haloDepth,physicalVelocity_->stages);
    velocity_=&workingVelocityView_;
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
    synchronize(rAU_,1);
    reconstructFaceFlux(predictedFlux_);
    correctedFlux_=predictedFlux_;
    if (target==System::TargetKind::Physical) publishVelocity();
}

void PressureOperators::publishVelocity() {
    if (velocity_==physicalVelocity_) return;
    // COPY the already boundary/halo-ready result. Only StateBundle owns the
    // physical field; the predictor buffer belongs to this numerical provider.
    for (int cell=0;cell<geometry_->TotalSize();++cell)
        for (int c=0;c<3;++c)
            physicalVelocity_->write(cell,c,velocity_->read(cell,c));
    velocity_=physicalVelocity_;
}

void PressureOperators::reconstructFaceFlux(std::vector<double>& fluxField) {
    std::fill(fluxField.begin(),fluxField.end(),0.0);
    std::fill(faceResponse_.begin(),faceResponse_.end(),0.0);
    for (int cell:candidateCells_) {
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
                fluxField[index]=flux;
            }
        }
    }
    synchronize(faceResponse_,3,State::ExchangeKind::CanonicalFaceFlux);
    synchronize(fluxField,3,State::ExchangeKind::CanonicalFaceFlux);
}

void PressureOperators::preparePressureBoundary() {
    applyVelocityBoundary();
    applyPressureBoundary();
}

void PressureOperators::assemblePressure() {
    summary_.maxDivergenceBefore=maxDivergence(correctedFlux_);
    pressureMatrix_.rows.clear();
    pressureMatrix_.rows.reserve(cells_.size());
    for (size_t index=0;index<cells_.size();++index) {
        const int cell=cells_[index];
        LinearAlgebra::GlobalDofRow row(dof(cell));
        if (!hasDirichlet_ && entityOfCell_[(size_t)cell]==referenceEntity_) {
            row.add(row.row(),1.0);
            row.setRightHandSide(0.0);
            pressureMatrix_.rows.push_back(std::move(row));
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
                    row.add(dof(neighbour),-coefficient);
                }
            }
        }
        if (!(diagonal>0.0))
            throw std::runtime_error("Pressure row has no connected faces.");
        row.add(row.row(),diagonal);
        double net=0.0;
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            net+=faceFlux(correctedFlux_,axis,cell)
                -faceFlux(correctedFlux_,axis,adjacent(cell,axis,-1));
        }
        row.setRightHandSide(-net);
        pressureMatrix_.rows.push_back(std::move(row));
    }
    traceAssembly();
}

void PressureOperators::traceAssembly() const {
    const char* directory=std::getenv("SF_PRESSURE_TRACE_DIR");
    if (!directory || !*directory || traceAssemblyIndex_++!=0) return;
    std::filesystem::create_directories(directory);
    std::ofstream output(std::filesystem::path(directory)/(
        "matrix-0-patch-"
        +std::to_string(velocity_->blockId)+".txt"));
    output << std::setprecision(17);
    std::unordered_map<std::int64_t,int> cells;
    for (int cell=0;cell<geometry_->TotalSize();++cell)
        if (solved(cell)) cells.emplace(dof(cell).value(),cell);
    const auto position=[&](LinearAlgebra::GlobalDofId id) {
        int i=0,j=0,k=0;
        geometry_->getIJK(cells.at(id.value()),i,j,k);
        output << geometry_->X(i,j,k) << ' ' << geometry_->Y(i,j,k)
               << ' ' << geometry_->Z(i,j,k);
    };
    for (const auto& row:pressureMatrix_.rows) {
        output << "row "; position(row.row());
        output << ' ' << row.rightHandSide() << '\n';
        for (const auto& entry:row.coefficients()) {
            output << "entry "; position(entry.column);
            output << ' ' << entry.value << '\n';
        }
    }
    if (!output) throw std::runtime_error("Pressure assembly trace write failed.");
}

void PressureOperators::solvePressure() {
    pressureResult_=pressureSolver_->solve(pressureMatrix_);
    if (!pressureResult_.converged
        || pressureResult_.solution.size()!=cells_.size())
        throw std::runtime_error("Constant-density pressure solve did not converge.");
    std::fill(correction_.begin(),correction_.end(),0.0);
    for (size_t row=0;row<cells_.size();++row)
        correction_[(size_t)cells_[row]]=pressureResult_.solution[row];
    synchronize(correction_,1);
    // Zero-gradient correction at velocity boundaries preserves the
    // pressure-null-space condition while supplying a centered cell gradient.
    for (size_t cell=0;cell<boundaryCell_.size();++cell) {
        if (!boundaryCell_[cell]) continue;
        correction_[cell]=pressureDirichlet_[cell]
            ? 0.0 : correction_[(size_t)inwardCell(*geometry_,(int)cell,physicalSides_)];
    }
    synchronize(correction_,1);
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
    restoreGauge();
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
    publishVelocity();
}

void PressureOperators::correctFlux() {
    for (int cell:candidateCells_) {
        for (int axis=0;axis<3;++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign:{-1,1}) {
                const int lower=sign>0?cell:adjacent(cell,axis,-1);
                const int upper=adjacent(lower,axis,1);
                if (sign<0 && rowOfCell_[(size_t)lower]>=0) continue;
                const auto index=faceIndex(axis,lower);
                if (faceResponse_[index]==0.0) continue;
                correctedFlux_[index]=RhieChow::correct(
                    correctedFlux_[index],faceResponse_[index],
                    correction_[(size_t)lower],correction_[(size_t)upper]);
            }
        }
    }
    synchronize(correctedFlux_,3,State::ExchangeKind::CanonicalFaceFlux);
    summary_.maxDivergenceAfter=maxDivergence(correctedFlux_);
}

PressureOperationSummary PressureOperators::commitCorrection() {
    if (!std::isfinite(summary_.maxDivergenceAfter))
        throw std::runtime_error("Corrected continuity defect is non-finite.");
    return summary_;
}

void PressureOperators::applyRelaxation() {
    if (!iterationActive_)
        throw std::runtime_error("Pressure relaxation requires active iteration.");
    candidateContinuityDefect_=maxDivergence(correctedFlux_);
    const int total=geometry_->TotalSize();
    for (int cell:cells_) {
        const double oldPressure=iterationPressure_[(size_t)cell];
        const double relaxedPressure=relaxSolution(oldPressure,
            pressure_->read(cell,0),config_.pressureRelaxation);
        if (!std::isfinite(relaxedPressure))
            throw std::runtime_error("Pressure relaxation produced non-finite p.");
        pressure_->write(cell,0,relaxedPressure);
        for (int c=0;c<3;++c) {
            const double oldVelocity=
                iterationVelocity_[(size_t)c*(size_t)total+(size_t)cell];
            const double relaxedVelocity=relaxSolution(oldVelocity,
                velocity_->read(cell,c),config_.velocityRelaxation);
            if (!std::isfinite(relaxedVelocity))
                throw std::runtime_error("Velocity relaxation produced non-finite U.");
            velocity_->write(cell,c,relaxedVelocity);
        }
    }
    restoreGauge();
    applyPressureBoundary();
    applyVelocityBoundary();
}

void PressureOperators::restoreFluxConsistency() {
    if (!iterationActive_)
        throw std::runtime_error("Flux consistency requires active iteration.");
    if (config_.velocityRelaxation!=1.0
        || config_.pressureRelaxation!=1.0) {
        // A relaxed U/p iterate needs its own Rhie-Chow face state; the
        // unrelaxed candidate flux cannot remain authoritative.
        reconstructFaceFlux(correctedFlux_);
    }
    summary_.maxDivergenceAfter=maxDivergence(correctedFlux_);
}

void PressureOperators::evaluateConvergence() {
    if (!iterationActive_)
        throw std::runtime_error("Pressure convergence requires active iteration.");
    const int total=geometry_->TotalSize();
    bool valid=true;
    for (int cell:cells_) {
        valid=std::isfinite(pressure_->read(cell,0)-iterationPressure_[(size_t)cell]) && valid;
        for (int c=0;c<3;++c)
            valid=std::isfinite(velocity_->read(cell,c)-iterationVelocity_[(size_t)c*total+cell]) && valid;
    }
    for (size_t n=0;n<correctedFlux_.size();++n)
        valid=std::isfinite(correctedFlux_[n]-iterationFlux_[n]) && valid;
    requireCollectively(valid,"Non-finite fixed-point difference before globalMaximum.");
    velocityDelta_=0.0;
    pressureDelta_=0.0;
    fluxDelta_=0.0;
    double velocityScale=1.0;
    double pressureScale=1.0;
    for (int cell:cells_) {
        pressureScale=std::max(pressureScale,
            std::abs(iterationPressure_[(size_t)cell]));
        pressureDelta_=std::max(pressureDelta_,std::abs(pressure_->read(cell,0)
            -iterationPressure_[(size_t)cell]));
        for (int c=0;c<3;++c) {
            velocityScale=std::max(velocityScale,std::abs(
                iterationVelocity_[(size_t)c*(size_t)total+(size_t)cell]));
            velocityDelta_=std::max(velocityDelta_,std::abs(velocity_->read(cell,c)
                -iterationVelocity_[(size_t)c*(size_t)total+(size_t)cell]));
        }
    }
    velocityDelta_=runtime_->globalMaximum(velocityDelta_);
    pressureDelta_=runtime_->globalMaximum(pressureDelta_);
    velocityScale=runtime_->globalMaximum(velocityScale);
    pressureScale=runtime_->globalMaximum(pressureScale);
    iterationDelta_=std::max(
        velocityDelta_/velocityScale,pressureDelta_/pressureScale);
    for (std::size_t index=0;index<correctedFlux_.size();++index)
        fluxDelta_=std::max(fluxDelta_,std::abs(
            correctedFlux_[index]-iterationFlux_[index]));
    fluxDelta_=runtime_->globalMaximum(fluxDelta_);
    summary_.maxDivergenceAfter=maxDivergence(correctedFlux_);
    if (!std::isfinite(iterationDelta_)
        || !std::isfinite(fluxDelta_)
        || !std::isfinite(summary_.maxDivergenceAfter))
        throw std::runtime_error("Pressure fixed-point residual is non-finite.");
    converged_=iterationDelta_<=config_.relativeTolerance
        && summary_.maxDivergenceAfter<=config_.absoluteTolerance;
}

void PressureOperators::endIteration() {
    if (!iterationActive_)
        throw std::runtime_error("Pressure iteration end has no active iteration.");
    iterationActive_=false;
}

void PressureOperators::commitStep() {
    if (iterationActive_)
        throw std::runtime_error("Cannot commit an incomplete pressure iteration.");
    applyVelocityBoundary();
    applyPressureBoundary();
    fixedTimeActive_=false;
}

} // namespace SF::Pressure
