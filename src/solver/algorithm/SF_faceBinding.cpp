#include "SF_faceBinding.h"
#include "core/system/SF_operationIds.h"
#include "core/mesh/SF_dimension.h"
#include <algorithm>
#include <cmath>
namespace SF::SolverAlgorithm {
FaceFluxExecution::FaceFluxExecution(const System::CompiledEquationCall& call,std::shared_ptr<CellEquationExecution> predictor,
        State::StateBundle& state,FDM::IExecutionRuntime& runtime,System::StateRealization& realized)
    :call_(call),contract_(*call.faceFlux),predictor_(std::move(predictor)),state_(state),runtime_(runtime) {
    const auto& contract=contract_;
    auto& u=realized.at(contract.velocity);auto& p=realized.at(contract.numericalInputs.at(1));
    if(u.fields.size()!=1 || p.fields.size()!=1 || !u.fields[0] || !p.fields[0] || u.fields[0]->components!=3 || p.fields[0]->components!=1
        || !p.fields[0]->read || u.fields[0]->geometry!=p.fields[0]->geometry || u.fields[0]->blockId!=p.fields[0]->blockId
        || predictor_->physical().name!=u.fields[0]->name || predictor_->call().source.occurrence!=contract.response)
        throw std::runtime_error("RhieChow input shape/geometry/response backing mismatch: "+contract.response);
    physical_=u.fields[0];pressure_=p.fields[0];auto& geometry=*physical_->geometry;const int total=geometry.TotalSize();
    predicted_.assign((size_t)total*3,0);response_.assign((size_t)total*3,0);
    if(realized.requestsView(contract.symbol,System::StateViewKind::Physical))corrected_.assign((size_t)total*3,0);
    pressureDirichlet_.assign(total,0);
    for(const auto& bc:predictor_->inputData(contract.numericalInputs.at(1))->cellBoundary)if(bc.kind==FDM::ScalarBoundaryKind::FixedValue)
        for(int cell:bc.cells)pressureDirichlet_.at(cell)=1;
    // This is the old pressure active-row mask, without creating a matrix or row numbering.
    std::vector<std::int64_t> active(total,-1);
    for(int cell:predictor_->candidates()) {
        int i,j,k;geometry.getIJK(cell,i,j,k);
        if(!runtime_.distributed() || geometry.isGlobalDofOwner(i,j,k))active[cell]=1;
    }
    runtime_.synchronizeIdentifiers(geometry,physical_->blockId,active);
    const auto solved=[&](int cell){return cell>=0 && cell<total && active[cell]>=0;};
    for(int cell:predictor_->candidates())for(int axis=0;axis<3;++axis)if(Math::isDirectionActiveIndex(axis))for(int sign:{-1,1}) {
        const int lower=sign>0?cell:Discretization::CellEquation::detail::adjacent(geometry,cell,axis,-1);
        const int upper=Discretization::CellEquation::detail::adjacent(geometry,lower,axis,1);
        const bool both=solved(lower)&&solved(upper);const int boundary=both?-1:(solved(lower)?upper:lower);
        const bool open=boundary>=0 && pressureDirichlet_[boundary];
        faces_.push_back({lower,upper,boundary,(size_t)axis*total+lower,
            Discretization::CellEquation::detail::faceArea(geometry,axis,lower),
            both || open?Discretization::CellEquation::detail::spacing(geometry,axis,lower):0.0,both,open});
    }
    predictedView_.name=contract.storage;predictedView_.blockId=physical_->blockId;predictedView_.geometry=&geometry;
    predictedView_.location=State::FieldLocation::Face;predictedView_.exchange=State::ExchangeKind::CanonicalFaceFlux;predictedView_.components=1;
    predictedView_.write=[](int,int,double){throw std::runtime_error("Face output is writable only by its selected numerical operation");};
    predictedView_.requireReadable=[this]{requirePrediction();};
    predictedView_.revision=[this]{return predictionGeneration_;};
    predictedView_.read=[this](int face,int){requirePrediction();return predicted_.at(face);};
    realized.bindView(contract.symbol,System::StateViewKind::Working,predictedView_);
    if(realized.requestsView(contract.symbol,System::StateViewKind::Physical)) {
        correctedView_=predictedView_;correctedView_.name=realized.at(contract.symbol).descriptor->storageKey;
        correctedView_.requireReadable=[this]{requireInterval();if(!correctedReady_ || generation_<=epochGeneration_)throw std::runtime_error("Unpublished corrected face flux: "+contract_.symbol);};
        correctedView_.revision=[this]{return generation_;};
        correctedView_.read=[this](int face,int){correctedView_.requireReadable();return corrected_.at(face);};
        realized.bindView(contract.symbol,System::StateViewKind::Physical,correctedView_);
    }
}
void FaceFluxExecution::requireInterval() const {
    interval_.require(state_,"step="+call_.source.occurrence+" equation="+call_.source.equation+" target="+contract_.symbol);
}
void FaceFluxExecution::beginEpoch() {
    interval_.open(state_,contract_.symbol);predictedReady_=correctedReady_=false;epochGeneration_=generation_;
}
void FaceFluxExecution::synchronize(std::vector<double>& data,int count,State::ExchangeKind kind) {
    bool finite=std::all_of(data.begin(),data.end(),[](double value){return std::isfinite(value);});
    if(runtime_.globalMinimum(finite?1.0:0.0)<.5)throw std::runtime_error("Nonfinite face resource before canonical COPY");
    auto view=State::workspaceView("face."+contract_.symbol,physical_->blockId,*physical_->geometry,data,count,kind);
    runtime_.synchronizeTransient({view});
}
void FaceFluxExecution::reconstruct(std::vector<double>& destination,const State::DistributedFieldView& velocity) {
    Pressure::RhieChow::construct(velocity,*pressure_,predictor_->response(),faces_,
        [this](int cell){return Discretization::CellEquation::centralGradient(*physical_->geometry,*pressure_,cell);},destination,response_);
    synchronize(response_,3,State::ExchangeKind::CanonicalFaceFlux);
    synchronize(destination,3,State::ExchangeKind::CanonicalFaceFlux);
}
void FaceFluxExecution::predict() {
    if(interval_.step!=state_.step)beginEpoch();
    requireInterval();predictor_->requireResponse();
    if(consumed_==predictor_->generation())throw std::runtime_error("RhieChow response already consumed: "+contract_.response);
    auto& input=contract_.velocityView==System::TargetKind::Working?predictor_->working():predictor_->physical();
    if(input.requireReadable)input.requireReadable();
    // Cell completion has already prepared these exact inputs; do not refresh mathematics here.
    synchronize(predictor_->response(),1,State::ExchangeKind::State);reconstruct(predicted_,input);
    if(!corrected_.empty())corrected_=predicted_;
    consumed_=predictor_->generation();++predictionGeneration_;predictedReady_=true;correctedReady_=false;
}
void FaceFluxExecution::requirePrediction() const {requireInterval();if(!predictedReady_)throw std::runtime_error("Unpublished predicted face flux: "+contract_.symbol);}
void FaceFluxExecution::publishCorrected() {requireInterval();if(corrected_.empty())throw std::runtime_error("No corrected face output requested: "+contract_.symbol);if(!correctedReady_)++generation_;correctedReady_=true;}
double FaceFluxExecution::published(std::string_view symbol,int axis,int lower) const {
    if(symbol!=contract_.symbol || axis<0 || axis>2)throw std::runtime_error("Missing published face resource: "+std::string(symbol));
    if(!correctedView_.requireReadable)throw std::runtime_error("No corrected publication is declared for face resource: "+contract_.symbol);
    correctedView_.requireReadable();const auto value=corrected_.at((size_t)axis*physical_->geometry->TotalSize()+lower);
    if(!std::isfinite(value))throw std::runtime_error("Nonfinite published face flux");return value;
}
}
