/// @file SF_pressureOperators.cpp
/// @brief Structured constant-density pressure operations.

#include "solver/algorithm/pressure/SF_pressureOperators.h"
#include "solver/algorithm/SF_iteration.h"

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
#include <set>
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
        SolverAlgorithm::FaceFluxExecution& faces)
    : config_(std::move(config)), faces_(faces),predictor_(faces.predictor()),rAU_(predictor_.response()),
      hByA_(predictor_.candidate()),physicalSides_(predictor_.sides()),boundaryCell_(predictor_.boundary()),
      pressureDirichlet_(faces.pressureDirichlet()),iterationVelocity_(predictor_.previousIterate()),predictedFlux_(faces.predicted()),
      faceResponse_(faces.response()),correctedFlux_(faces.corrected()) {}

void PressureOperators::bind(System::StateRealization& state,FDM::IExecutionRuntime& runtime) {
    runtime_=&runtime;
    velocity_=physicalVelocity_=&predictor_.physical();pressure_=&faces_.pressure();geometry_=velocity_->geometry;
    const auto& rho=state.at("rho");
    if(!rho.descriptor || !rho.descriptor->constantValue || !std::isfinite(*rho.descriptor->constantValue)
        || *rho.descriptor->constantValue<=0)throw std::runtime_error("Pressure correction needs positive declared rhoConst");
    rowOfCell_.assign(geometry_->TotalSize(),-1);cells_=predictor_.candidates();
    for(std::size_t n=0;n<cells_.size();++n)rowOfCell_[cells_[n]]=static_cast<int>(n);
    correction_.assign(geometry_->TotalSize(),0);
    correctionView_=State::workspaceView("pressureCorrection",pressure_->blockId,*geometry_,correction_,1,pressure_->haloDepth,pressure_->stages);
    state.bindView("deltaP",System::StateViewKind::Correction,correctionView_);
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

void PressureOperators::applyVelocityBoundary() {predictor_.prepareOutput(*velocity_,true);}

void PressureOperators::applyPressureBoundary() {
    predictor_.prepareInput(faces_.contract().numericalInputs.at(1));
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
    faces_.invalidateCorrected();
    requireCollectively(velocity_ && pressure_ && std::isfinite(maximumTimeStep)
        && maximumTimeStep>0.0,"Pressure prepare requires bound state and positive dt limit.");
    applyPressureBoundary();
    dt_=maximumTimeStep;
    faces_.beginEpoch();
    summary_={};
    fixedTimeActive_=false;
    iterationActive_=false;
    return dt_;
}

void PressureOperators::beginFixedTimeStep() {
    if (!(dt_>0.0) || !velocity_ || !pressure_)
        throw std::runtime_error("Fixed-time pressure step requires prepared state.");
    const int total=geometry_->TotalSize();
    iterationPressure_.resize((size_t)total);
    iterationFlux_.resize(correctedFlux_.size());
    std::vector<int> responseCells;for(int cell=0;cell<geometry_->TotalSize();++cell)if(solved(cell))responseCells.push_back(cell);
    predictor_.setResponse(dt_,&responseCells);
    // The first F_k must match U_n/p_n at this timestep's pressure response,
    // rather than inheriting an uninitialized or previous-dt face buffer.
    faces_.reconstruct(correctedFlux_,*velocity_);
    if(!predictor_.hasFrozenBase())throw std::runtime_error("Pressure fixed-time relation requires the ordinary entry snapshot");
    fixedTimeActive_=true;
    candidateContinuityDefect_=0.0;
}

void PressureOperators::beginIteration() {
    if (!fixedTimeActive_ || iterationActive_)
        throw std::runtime_error("Pressure iteration requires one fixed-time step.");
    const int total=geometry_->TotalSize();
    for (int cell=0;cell<total;++cell) {
        iterationPressure_[(size_t)cell]=pressure_->read(cell,0);
    }
    iterationFlux_=correctedFlux_;
    iterationActive_=true;
    predictor_.beginFixedIteration();
}

void PressureOperators::publishVelocity() {
    if(velocity_!=physicalVelocity_) {
        predictor_.copyToPhysical();velocity_=physicalVelocity_;
    }
}

void PressureOperators::preparePressureBoundary() {
    faces_.requirePrediction();velocity_=&predictor_.output();
    applyVelocityBoundary();
    applyPressureBoundary();
}

void PressureOperators::assemblePressure() {
    faces_.requirePrediction();
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

void PressureOperators::traceState(const char* phase,int step,double time) const {
    const char* directory=std::getenv("SF_PRESSURE_TRACE_DIR");
    const char* selected=std::getenv("SF_PRESSURE_TRACE_STEP");
    if(!directory || !*directory || !selected || !*selected)return;
    std::size_t parsed=0;const int wanted=std::stoi(selected,&parsed);
    if(selected[parsed]!='\0' || wanted<0)throw std::runtime_error("Invalid pressure trace step");
    if(step!=0 && step!=wanted)return;
    if(runtime_->distributed())throw std::runtime_error("Pressure verification snapshots require serial execution; distributed snapshot ownership is not implemented");
    std::filesystem::create_directories(directory);
    const auto prefix=std::filesystem::path(directory)/("state-"+std::to_string(step)+"-"+std::to_string(traceStateIndex_++)+"-"+phase);
    std::ofstream cells(prefix.string()+".csv");
    cells<<std::setprecision(17)<<"cell,i,j,k,x,y,z,solved,inverseVolume,dt,U0,U1,U2,p,H0,H1,H2,rAU,deltaP\n";
    for(int cell=0;cell<geometry_->TotalSize();++cell) {
        int i,j,k;geometry_->getIJK(cell,i,j,k);
        cells<<cell<<','<<i<<','<<j<<','<<k<<','<<geometry_->X(i,j,k)<<','<<geometry_->Y(i,j,k)<<','<<geometry_->Z(i,j,k)<<','<<solved(cell)<<','<<geometry_->Jac(i,j,k)<<','<<dt_;
        for(int c=0;c<3;++c)cells<<','<<velocity_->read(cell,c);
        cells<<','<<pressure_->read(cell,0);
        for(int c=0;c<3;++c)cells<<','<<hByA_[(size_t)c*geometry_->TotalSize()+cell];
        cells<<','<<rAU_[cell]<<','<<correction_[cell]<<'\n';
    }
    std::ofstream faces(prefix.string()+"-faces.csv");
    faces<<std::setprecision(17)<<"axis,lower,upper,Ax,Ay,Az,response,flux,predicted,deltaL,deltaR\n";
    std::set<std::pair<int,int>> visited;
    for(int cell:cells_)for(int axis=0;axis<3;++axis)if(Math::isDirectionActiveIndex(axis))for(int sign:{-1,1}) {
        const int lower=sign>0?cell:adjacent(cell,axis,-1),upper=adjacent(lower,axis,1);
        if(!visited.emplace(axis,lower).second)continue;
        const auto index=faceIndex(axis,lower);const auto area=faceArea(axis,lower);
        faces<<axis<<','<<lower<<','<<upper<<','<<area[0]<<','<<area[1]<<','<<area[2]<<','<<faceResponse_[index]<<','<<correctedFlux_[index]<<','<<predictedFlux_[index]<<','<<correction_[lower]<<','<<correction_[upper]<<'\n';
    }
    if(std::string_view(phase)=="accepted") {
        std::ofstream matrix(prefix.string()+"-matrix.csv");matrix<<std::setprecision(17)<<"row,column,coefficient,rhs,solution\n";
        for(size_t n=0;n<pressureMatrix_.rows.size();++n) {
            const auto& row=pressureMatrix_.rows[n];
            for(const auto& entry:row.coefficients())matrix<<row.row().value()<<','<<entry.column.value()<<','<<entry.value<<','<<row.rightHandSide()<<','<<pressureResult_.solution.at(n)<<'\n';
        }
        if(!matrix)throw std::runtime_error("Pressure trace matrix write failed");
    }
    if(!cells || !faces)throw std::runtime_error("Pressure trace state write failed");
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
    predictor_.invalidateWorking();
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
    faces_.invalidateCorrected();
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
    faces_.synchronize(correctedFlux_,3,State::ExchangeKind::CanonicalFaceFlux);
    summary_.maxDivergenceAfter=maxDivergence(correctedFlux_);
}

PressureOperationSummary PressureOperators::commitCorrection() {
    if (!std::isfinite(summary_.maxDivergenceAfter))
        throw std::runtime_error("Corrected continuity defect is non-finite.");
    faces_.publishCorrected();
    return summary_;
}

void PressureOperators::applyRelaxation() {
    predictor_.invalidateWorking();
    faces_.invalidateCorrected();
    if (!iterationActive_)
        throw std::runtime_error("Pressure relaxation requires active iteration.");
    candidateContinuityDefect_=maxDivergence(correctedFlux_);
    const int total=geometry_->TotalSize();
    for (int cell:cells_) {
        const double oldPressure=iterationPressure_[(size_t)cell];
        const double relaxedPressure=SolverAlgorithm::relaxSolution(oldPressure,
            pressure_->read(cell,0),config_.pressureRelaxation);
        if (!std::isfinite(relaxedPressure))
            throw std::runtime_error("Pressure relaxation produced non-finite p.");
        pressure_->write(cell,0,relaxedPressure);
        for (int c=0;c<3;++c) {
            const double oldVelocity=
                iterationVelocity_[(size_t)c*(size_t)total+(size_t)cell];
            const double relaxedVelocity=SolverAlgorithm::relaxSolution(oldVelocity,
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
        faces_.reconstruct(correctedFlux_,*velocity_);
    }
    summary_.maxDivergenceAfter=maxDivergence(correctedFlux_);
    faces_.publishCorrected();
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
}

void PressureOperators::endIteration() {
    if (!iterationActive_)
        throw std::runtime_error("Pressure iteration end has no active iteration.");
    iterationActive_=false;
}

void PressureOperators::commitStep() {
    predictor_.invalidateWorking();
    if (iterationActive_)
        throw std::runtime_error("Cannot commit an incomplete pressure iteration.");
    applyVelocityBoundary();
    applyPressureBoundary();
    fixedTimeActive_=false;
}

} // namespace SF::Pressure
