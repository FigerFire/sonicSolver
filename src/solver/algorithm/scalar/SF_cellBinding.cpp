#include "SF_cellBinding.h"
#include "solver/algorithm/SF_iteration.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <cmath>
namespace SF::SolverAlgorithm {
namespace {
int inward(const Field& field,int cell,const std::array<bool,6>& sides) {
    int i,j,k;field.getIJK(cell,i,j,k);const int g=field.NG();
    // Preserve the original ordered side donor, including grouped-wall corners.
    if(i==g && sides[0])++i;
    else if(i==g+field.NX()-1 && sides[1])--i;
    else if(j==g && sides[2])++j;
    else if(j==g+field.NY()-1 && sides[3])--j;
    else if(k==g && sides[4])++k;
    else if(k==g+field.NZ()-1 && sides[5])--k;
    return field.getIdx(i,j,k);
}
}
CellEquationExecution::CellEquationExecution(const System::CompiledEquationCall& call,const System::CompiledNumericalSystem& numerical,
        State::StateBundle& state,const FDM::SolverServices& services,System::StateRealization& realized,bool responseRequired)
    :call_(call),contract_(*System::cellEquation(call)),numerical_(numerical),state_(state),realized_(realized),runtime_(services.executionRuntime),mesh_(state.singlePatch()) {
    if(call.temporalResidual && responseRequired)throw std::runtime_error("Recipe participant cannot provide a Step Euler face response: "+call.source.occurrence);
    if(call.temporalResidual && runtime_ && runtime_->distributed())throw std::runtime_error("Cell recipe boundary backing currently requires one serial patch: "+call.source.occurrence);
    if(state.patches.size()!=1)throw std::runtime_error("Cell equation needs one local patch: "+call.source.occurrence);
    for(const auto& need:call.physicalResources) {
        if(!need.fields || !need.fields->cellVolumeBoundary)throw std::runtime_error("Cell equation lacks its boundary resource contract");
        const auto& instance=equationData(services,need.id);auto data=instance.data.lock();
        if(!data || data->target!=need.fields->target)throw std::runtime_error("Missing cell field resource: "+need.id);
        auto& backing=equationField(realized,{data->target},mesh_,true,true);
        FieldPort port{&backing,data,std::vector<unsigned char>(mesh_.TotalSize(),0)};
        for(const auto& bc:data->cellBoundary)for(int cell:bc.cells) {
            if(cell<0 || cell>=mesh_.TotalSize())throw std::runtime_error("Cell boundary outside geometry");
            port.boundary[cell]=1;
        }
        const int g=mesh_.NG(),n[3]={mesh_.NX(),mesh_.NY(),mesh_.NZ()};
        for(int axis=0;axis<3;++axis)if(data->activeDimensions[axis] && n[axis]>1)for(int side=0;side<2;++side) {
            bool physical=true;
            for(int k=g;k<g+n[2];++k)for(int j=g;j<g+n[1];++j)for(int i=g;i<g+n[0];++i)
                if((axis==0?i:axis==1?j:k)==g+(side?n[axis]-1:0))physical=physical && port.boundary[mesh_.getIdx(i,j,k)];
            port.sides[2*axis+side]=physical;
            if(need.fields->boundaryRequired && (!runtime_ || !runtime_->distributed()) && !physical)throw std::runtime_error("Missing complete cell boundary: "+data->target);
        }
        fields_.push_back(std::move(port));
    }
    physical_=fields_.front().view;
    if(physical_->components!=3)throw std::runtime_error("CellExplicit requires actual three-component backing");
    const int g=mesh_.NG();for(int k=g;k<g+mesh_.NZ();++k)for(int j=g;j<g+mesh_.NY();++j)for(int i=g;i<g+mesh_.NX();++i) {
        const int cell=mesh_.getIdx(i,j,k);if(!boundary()[cell] && mesh_.CellFlag(i,j,k)==FLUID_CELL)cells_.push_back(cell);
    }
    candidates_=cells_;
    if(runtime_)cells_=runtime_->canonicalOwnerCells(mesh_,cells_);
    const auto targetData=fields_.front().data.lock();
    for(int cell:cells_)for(int axis=0;axis<3;++axis)if(call.physicalResources.front().fields->boundaryRequired && targetData->activeDimensions[axis] && (axis==0?mesh_.NX():axis==1?mesh_.NY():mesh_.NZ())>1)for(int sign:{-1,1}) {
        const int lower=sign>0?cell:Discretization::CellEquation::detail::adjacent(mesh_,cell,axis,-1);
        const int upper=Discretization::CellEquation::detail::adjacent(mesh_,lower,axis,1);
        const auto area=Discretization::CellEquation::detail::faceArea(mesh_,axis,lower);
        int i,j,k,ni,nj,nk;mesh_.getIJK(lower,i,j,k);mesh_.getIJK(upper,ni,nj,nk);
        const double d[3]={mesh_.X(ni,nj,nk)-mesh_.X(i,j,k),mesh_.Y(ni,nj,nk)-mesh_.Y(i,j,k),mesh_.Z(ni,nj,nk)-mesh_.Z(i,j,k)};
        double dot=0,a2=0,d2=0;for(int c=0;c<3;++c){dot+=area[c]*d[c];a2+=area[c]*area[c];d2+=d[c]*d[c];}
        if(!(a2>0) || !(d2>0) || std::abs(dot)/std::sqrt(a2*d2)<1.0-1.0e-8)throw std::runtime_error("Cell stencil requires an orthogonal structured grid");
    }
    if(cells_.empty())throw std::runtime_error("Cell equation has no owned interior cells");
    if(responseRequired)candidate_.resize((size_t)mesh_.TotalSize()*3);
    if(responseRequired){if(contract_.potential.empty())throw std::runtime_error("This cell equation has no scalar gradient response: "+call.source.occurrence);response_.resize(mesh_.TotalSize());}
    for(const auto index:contract_.terms) {
        const auto& term=numerical.operators.at(index);
        terms_.push_back({&term,term.primitiveSource});
    }
    if(runtime_ && runtime_->distributed() && std::any_of(terms_.begin(),terms_.end(),[](const auto& term){return term.compiled->primitiveSpatial==&Discretization::CellEquation::stressArithmeticCentral2;}))
        throw std::runtime_error("Cell stress arithmetic-face publication currently needs one serial patch; canonical MPI traction exchange is not implemented: "+call.source.occurrence);
    retainWorking_=realized.requestsView(call.target.symbol,System::StateViewKind::Working);
    if(retainWorking_) {
        working_.resize((size_t)mesh_.TotalSize()*3);
        const auto storage=call.target.kind==System::TargetKind::Working?call.target.workspace:call.viewOutputs.at(0).storage;
        workingView_=State::workspaceView(storage,physical_->blockId,mesh_,working_,3,physical_->haloDepth,physical_->stages);
        workingView_.requireReadable=[this]{interval_.require(state_,call_.source.occurrence);if(!ready_)throw std::runtime_error("Unpublished/expired Working view: "+call_.source.occurrence);};
        workingView_.revision=[this]{return revision_;};
        workingView_.write=[this](int cell,int component,double value){working_.at((size_t)component*mesh_.TotalSize()+cell)=value;++revision_;};
        workingView_.prepareSpatial=[this](double time){workingView_.requireReadable();if(time!=state_.time)throw std::runtime_error("Working cell input has a different interval");prepareOutput(workingView_);};
        realized.bindView(call.target.symbol,System::StateViewKind::Working,workingView_);
    }
    if(!call.temporalResidual)for(std::size_t n=0;n<fields_.size();++n)fields_[n].view->prepareSpatial=[this,n](double time){
        if(time!=state_.time)throw std::runtime_error("Cell boundary preparation time differs from current interval");prepareField(fields_[n],*fields_[n].view);
    };
    if(call.temporalResidual) {
        old_.resize((size_t)mesh_.TotalSize()*3);stage_.resize(old_.size());
        for(auto& increment:increments_)increment.resize(old_.size());
        oldView_=State::workspaceView(call.oldTimeWorkspace,physical_->blockId,mesh_,old_,3);
        stageView_=State::workspaceView(call.stageWorkspace,physical_->blockId,mesh_,stage_,3);
        stageView_.requireReadable=[this]{realized_.requireStageRead(nextStage_,evaluationTime_);};
        stageView_.revision=[this]{return revision_;};
        stageView_.prepareSpatial=[this](double time){stageView_.requireReadable();if(time!=evaluationTime_)throw std::runtime_error("Cell Stage boundary time mismatch");prepareField(fields_.front(),stageView_);};
        if(realized.requestsView(call.target.symbol,System::StateViewKind::OldTime))realized.bindView(call.target.symbol,System::StateViewKind::OldTime,oldView_);
        for(int n=0;n<numerical.time.recipe.stageCount();++n)realized.bindView(call.target.symbol,System::StateViewKind::Stage,stageView_,n);
    }

}
void CellEquationExecution::prepareField(FieldPort& port,State::DistributedFieldView& view) {
    const auto data=port.data.lock();if(!data)throw std::runtime_error("Expired cell boundary resource");
    const bool working=&view==&workingView_ || &view==&stageView_;
    const auto revision=working?revision_:realized_.revision(data->target);
    if(working && preparedRevision_==revision && preparedTime_==evaluationTime_)return;
    if(!working && port.prepared==revision && port.time==evaluationTime_)return;
    for(auto kind:{FDM::ScalarBoundaryKind::ZeroGradient,FDM::ScalarBoundaryKind::FixedValue})for(const auto& bc:data->cellBoundary)if(bc.kind==kind)
        for(int cell:bc.cells) {
            const int donor=inward(mesh_,cell,port.sides);
            for(int c=0;c<view.components;++c)view.write(cell,c,kind==FDM::ScalarBoundaryKind::FixedValue?bc.value.at(c):view.read(donor,c));
        }
    const int g=mesh_.NG();for(int k=0;k<mesh_.MZ();++k)for(int j=0;j<mesh_.MY();++j)for(int i=0;i<mesh_.MX();++i) {
        int ci=std::clamp(i,g,g+mesh_.NX()-1),cj=std::clamp(j,g,g+mesh_.NY()-1),ck=std::clamp(k,g,g+mesh_.NZ()-1);
        if((i==ci && j==cj && k==ck) || mesh_.isCommunicationHalo(i,j,k))continue;
        for(int c=0;c<view.components;++c)view.write(mesh_.getIdx(i,j,k),c,view.read(mesh_.getIdx(ci,cj,ck),c));
    }
    bool finite=true;for(int cell=0;cell<mesh_.TotalSize();++cell)for(int c=0;c<view.components;++c)finite=std::isfinite(view.read(cell,c)) && finite;
    if(runtime_)finite=runtime_->globalMinimum(finite?1.0:0.0)>.5;
    if(!finite)throw std::runtime_error("Nonfinite cell field before canonical COPY");
    if(runtime_)runtime_->synchronizeTransient({view});
    if(working){preparedRevision_=revision_;preparedTime_=evaluationTime_;}
    else {port.prepared=realized_.revision(data->target);port.time=evaluationTime_;}
}
void CellEquationExecution::prepareOutput(State::DistributedFieldView& view,bool force) {
    if(force){if(&view==&workingView_)preparedRevision_=~0ULL;else fields_.front().prepared=~0ULL;}prepareField(fields_.front(),view);
}
std::shared_ptr<const FDM::ScalarEquationData> CellEquationExecution::inputData(std::string_view symbol) const {
    for(const auto& field:fields_)if(auto data=field.data.lock();data && data->target==symbol)return data;
    throw std::runtime_error("Missing bound cell input field: "+std::string(symbol));
}
void CellEquationExecution::prepareInput(std::string_view symbol) {
    for(auto& field:fields_)if(auto data=field.data.lock();data && data->target==symbol) {
        field.prepared=~0ULL;prepareField(field,*field.view);return;
    }
    throw std::runtime_error("Missing bound cell boundary input: "+std::string(symbol));
}
void CellEquationExecution::bindInputs() {
    if(inputsBound_)throw std::runtime_error("Cell inputs already bound: "+call_.source.occurrence);
    const auto useOf=[&](std::string_view symbol)->const System::StateUse& {
        const auto found=std::find_if(call_.stateUses.begin(),call_.stateUses.end(),[&](const auto& input){return input.symbol==symbol;});
        if(found==call_.stateUses.end())throw std::runtime_error("Missing cell STATE input: "+std::string(symbol));
        return *found;
    };
    for(auto& bound:terms_) {
        const auto& term=*bound.compiled;
        if(term.fieldSource.bind) {
            std::vector<System::FieldSource::Reader> reads;
            for(const auto& use:term.fieldSource.inputs)reads.push_back(realized_.reader(use,numerical_.time.recipe.stageCount()));
            bound.source=[value=term.fieldSource.bind(std::move(reads))](const Field& field,int i,int j,int k,const Vector3&) {
                const int cell=field.getIdx(i,j,k);
                return Vector3(value(cell,0),value(cell,1),value(cell,2));
            };
        }
        if(!term.coefficientState.empty()) {
            const auto& use=useOf(term.coefficientState);
            bound.coefficient=bindInputView(use);
            const auto& symbol=*realized_.at(use.symbol).descriptor;
            bound.prepareCoefficient=symbol.storageBinding!=System::StorageBinding::TransientWorkspace
                && (use.version==System::StateVersion::Stage || use.kind==System::TargetKind::Working
                    || (call_.temporalResidual && use.version==System::StateVersion::Current));
        }
    }
    for(const auto& gradient:contract_.gradients)
        gradients_.push_back(bindInputView(useOf(gradient.input)));
    inputsBound_=true;
}
void CellEquationExecution::begin() {if(preparedStep_==state_.step)throw std::runtime_error("Repeated physical interval preparation: "+call_.source.occurrence);preparedStep_=state_.step;evaluationTime_=state_.time;for(auto& field:fields_)field.prepared=~0ULL;interval_={};updated_=-1;ready_=assembled_=false;base_=nullptr;}
double CellEquationExecution::stepSize() {
    evaluationTime_=state_.time;
    prepareOutput(*physical_);
    double limit=numerical_.dt.maxDeltaT;
    for(int cell:cells_) {
        double speed2=0;for(int c=0;c<3;++c){double u=physical_->read(cell,c);speed2+=u*u;}
        double speed=std::sqrt(speed2)*contract_.advectionWeight;
        const auto data=fields_.front().data.lock();if(!data)throw std::runtime_error("Expired cell field port");
        for(int axis=0;axis<3;++axis)if(data->activeDimensions[axis] && (axis==0?mesh_.NX():axis==1?mesh_.NY():mesh_.NZ())>1) {
            const double h=Discretization::CellEquation::detail::spacing(mesh_,axis,cell);
            if(speed>0)limit=std::min(limit,numerical_.dt.cfl*h/speed);
            if(contract_.diffusivity>0)limit=std::min(limit,numerical_.dt.cfl*h*h/(2.0*std::count(data->activeDimensions.begin(),data->activeDimensions.end(),true)*contract_.diffusivity));
        }
    }
    if(runtime_)limit=runtime_->globalMinimum(limit);
    if(!std::isfinite(limit) || limit<=0)throw std::runtime_error("Invalid cell time-step bound");return limit;
}
void CellEquationExecution::setResponse(double dt,const std::vector<int>* region) {if(response_.empty())return;for(int cell:region?*region:cells_)response_[cell]=contract_.potential.empty()?0:dt/contract_.capacityValue;}
void CellEquationExecution::assemble() {
    if(!inputsBound_ || preparedStep_!=state_.step || !(state_.dt>0))throw std::runtime_error("Cell RHS without prepared physical interval");
    if(interval_.step<0)interval_.open(state_,call_.source.occurrence);
    interval_.require(state_,call_.source.occurrence);
    evaluate(*physical_);
    if(!candidate_.empty())for(int cell:cells_)for(int c=0;c<3;++c){const auto at=(size_t)c*mesh_.TotalSize()+cell;candidate_[at]=fixedTimeCandidate(base_?base_->at(at):physical_->read(cell,c),state_.dt*withoutGradient_[at]);}
    setResponse(state_.dt);assembled_=true;
}
CellEquationExecution::ViewReader CellEquationExecution::bindInputView(const System::StateUse& use) {
    if(use.version==System::StateVersion::Stage) {
        std::vector<State::DistributedFieldView*> stages;
        for(int stage=0;stage<numerical_.time.recipe.stageCount();++stage) {
            const auto& view=realized_.view(use.symbol,System::StateViewKind::Stage,stage);
            if(view.fields.size()!=1 || view.fields.front()->geometry!=&mesh_
                || view.fields.front()->components!=view.descriptor.components || view.descriptor.componentOffset)
                throw std::runtime_error("Cell Stage input requires unsliced backing on the same patch: "+use.symbol);
            stages.push_back(view.fields.front());
        }
        return [this,stages=std::move(stages)]()->State::DistributedFieldView& {return *stages.at(nextStage_);};
    }
    if(use.version==System::StateVersion::Frozen || use.version==System::StateVersion::OldTime) {
        const auto kind=use.version==System::StateVersion::Frozen?System::StateViewKind::Frozen:System::StateViewKind::OldTime;
        return [view=realized_.view(use.symbol,kind).fields.at(0)]()->State::DistributedFieldView& {return *view;};
    }
    return [view=&equationField(realized_,use,mesh_)]()->State::DistributedFieldView& {return *view;};
}
void CellEquationExecution::evaluate(State::DistributedFieldView& input) {
    if(!inputsBound_)throw std::runtime_error("Cell input readers not bound: "+call_.source.occurrence);
    if(input.requireReadable)input.requireReadable();
    if(!call_.temporalResidual)for(auto& field:fields_)prepareField(field,*field.view);
    const auto data=fields_.front().data.lock();if(!data)throw std::runtime_error("Expired cell field port");
    auto flux=contract_.advectionWeight>0?Discretization::CellEquation::candidateFlux(input,mesh_,cells_,data->activeDimensions):std::vector<double>{};
    Discretization::CellEquation::Stencil stencil{input,mesh_,flux,data->activeDimensions};
    std::vector<Discretization::CellEquation::Contribution> terms;
    for(const auto& bound:terms_) {
        const auto& term=*bound.compiled;
        State::DistributedFieldView* coefficient=nullptr;
        if(bound.coefficient) {
            coefficient=&bound.coefficient();
            if(coefficient->requireReadable)coefficient->requireReadable();
            if(coefficient->geometry!=&mesh_ || coefficient->components!=1)
                throw std::runtime_error("Stress coefficient view differs from the target geometry or scalar shape: "+term.coefficientState);
            if(bound.prepareCoefficient) {
                if(!coefficient->prepareSpatial)throw std::runtime_error("Stress coefficient lacks actual-version boundary preparation: "+term.coefficientState);
                coefficient->prepareSpatial(evaluationTime_);
            }
        }
        terms.push_back({term.primitiveSpatial,bound.source,term.coefficient,coefficient,term.constantCoefficient});
    }
    std::vector<Discretization::CellEquation::GradientContribution> gradients;
    for(std::size_t index=0;index<contract_.gradients.size();++index) {
        const auto& term=contract_.gradients[index];
        auto& view=gradients_[index]();
        if(view.requireReadable)view.requireReadable();
        if(call_.temporalResidual){if(!view.prepareSpatial)throw std::runtime_error("Missing real Stage gradient boundary: "+term.input);view.prepareSpatial(evaluationTime_);}
        gradients.push_back({&view,term.kernel,term.weight,term.divisor});
    }
    Discretization::CellEquation::evaluate(stencil,cells_,terms,gradients,rhs_,candidate_.empty()?nullptr:&withoutGradient_);
}
void CellEquationExecution::advance() {
    interval_.require(state_,call_.source.occurrence);
    if(!assembled_)throw std::runtime_error("Cell advance before RHS: "+call_.source.occurrence);
    auto* output=physical_;
    if(retainWorking_) {
        for(int cell=0;cell<mesh_.TotalSize();++cell)for(int c=0;c<3;++c)working_[(size_t)c*mesh_.TotalSize()+cell]=physical_->read(cell,c);
        output=&workingView_;
    }
    for(int cell:cells_)for(int c=0;c<3;++c){const auto at=(size_t)c*mesh_.TotalSize()+cell;output->write(cell,c,fixedTimeCandidate(base_?base_->at(at):physical_->read(cell,c),state_.dt*rhs_[at]));}
    if(output==physical_)fields_.front().prepared=~0ULL;
    prepareOutput(*output);
    if(retainWorking_){++revision_;preparedRevision_=revision_;}
    ready_=retainWorking_;updated_=state_.step;assembled_=false;
    if(call_.target.kind==System::TargetKind::Physical){copyToPhysical();state_.distributed.markModified(physical_->name,physical_->blockId);}
    predictionRevision_=retainWorking_?revision_:realized_.revision(call_.target.symbol);
    ++generation_;
    if(!contract_.potential.empty())potentialRevision_=realized_.revision(contract_.potential);
}
void CellEquationExecution::copyToPhysical() {
    if(retainWorking_)for(int cell=0;cell<mesh_.TotalSize();++cell)for(int c=0;c<3;++c)physical_->write(cell,c,workingView_.read(cell,c));
}
void CellEquationExecution::requireResponse() const {
    interval_.require(state_,call_.source.occurrence);
    if(response_.empty() || updated_!=state_.step || predictionRevision_!=(retainWorking_?revision_:realized_.revision(call_.target.symbol))
        || (!contract_.potential.empty() && potentialRevision_!=realized_.revision(contract_.potential)))
        throw std::runtime_error("Expired/unpublished response occurrence: "+call_.source.occurrence);
}
void CellEquationExecution::accept() const {interval_.require(state_,call_.source.occurrence);if(updated_!=state_.step)throw std::runtime_error("Cell result not produced in this physical step: "+call_.source.occurrence);}
void CellEquationExecution::freezeBase() {
    if(call_.temporalResidual || frozenStep_==state_.step || preparedStep_!=state_.step)
        throw std::runtime_error("Invalid/repeated fixed-time snapshot: "+call_.source.occurrence);
    if(interval_.step<0)interval_.open(state_,call_.source.occurrence);
    snapshotValues(*physical_,old_);frozenStep_=state_.step;fixedIteration_=0;
}
void CellEquationExecution::beginFixedIteration() {
    interval_.require(state_,call_.source.occurrence);
    if(frozenStep_!=state_.step)throw std::runtime_error("Missing fixed-time entry snapshot: "+call_.source.occurrence);
    snapshotValues(*physical_,previousIterate_);++fixedIteration_;
    base_=fixedIteration_>1?&old_:nullptr;
}
void CellEquationExecution::snapshot() {
    begin();interval_.open(state_,call_.source.occurrence);nextStage_=0;
    snapshotValues(*physical_,old_);
    stage_=old_;++revision_;preparedRevision_=~0ULL;
}
void CellEquationExecution::prepareStage(const Run::ExecutionContext& context) {
    numerical_.time.recipe.requireProviderStages(context.stageCount);
    if(context.stageIndex!=nextStage_)throw std::runtime_error("Cell stage order mismatch: "+call_.source.occurrence);
    evaluationTime_=state_.time+numerical_.time.recipe.stage(nextStage_).abscissa*state_.dt;
    prepareField(fields_.front(),stageView_);
}
void CellEquationExecution::evaluateStage(const Run::ExecutionContext& context) {
    interval_.require(state_,call_.source.occurrence);realized_.requireStageRead(context.stageIndex,evaluationTime_);
    evaluate(stageView_);increments_.at(context.stageIndex)=rhs_;
}
void CellEquationExecution::advanceStage(const Run::ExecutionContext& context) {
    realized_.requireStageAdvance(context.stageIndex,evaluationTime_);
    for(int cell:cells_)for(int c=0;c<3;++c){const auto at=(size_t)c*mesh_.TotalSize()+cell;std::array<double,4> increments{};for(int n=0;n<numerical_.time.recipe.stageCount();++n)increments[n]=increments_[n][at];stage_[at]=explicitValue(numerical_.time.recipe,context.stageIndex,state_.dt,old_[at],stage_[at],rhs_[at],increments);}
    ++revision_;++nextStage_;
}
void CellEquationExecution::validatePublish() {
    if(nextStage_!=numerical_.time.recipe.stageCount())throw std::runtime_error("Incomplete cell recipe: "+call_.source.occurrence);
    evaluationTime_=state_.time+state_.dt;prepareField(fields_.front(),stageView_);
    for(double value:stage_)if(!std::isfinite(value))throw std::runtime_error("Nonfinite cell final state");
}
void CellEquationExecution::publish() {
    for(int cell=0;cell<mesh_.TotalSize();++cell)for(int c=0;c<3;++c)physical_->write(cell,c,stage_[(size_t)c*mesh_.TotalSize()+cell]);
    state_.distributed.markModified(physical_->name,physical_->blockId);nextStage_=-1;
}
BoundEquation bindCell(const EquationBindingInput& input,bool responseRequired,bool fixedBaseRequired,
        std::shared_ptr<CellEquationExecution>& instance) {
    const auto& call=input.call;
    instance=std::make_shared<CellEquationExecution>(call,input.numerical,input.state,input.services,input.views,responseRequired);
    auto owner=instance;
    BoundEquation bound;
    bound.connect=[owner](const FluxReaders&){owner->bindInputs();};
    if(call.temporalResidual) {
        auto& cb=bound.stage;
        cb.identity=call.source.occurrence;
        cb.stepSize=[owner]{return owner->stepSize();};
        cb.snapshot=[owner]{owner->snapshot();};
        cb.prepare=[owner](const auto& c){owner->prepareStage(c);};
        cb.rhs=[owner](const auto& c){owner->evaluateStage(c);};
        cb.advance=[owner](const auto& c){owner->advanceStage(c);};
        cb.validatePublish=[owner]{owner->validatePublish();};
        cb.publish=[owner]{owner->publish();};
        return bound;
    }
    bound.interval.prepare=[owner]{owner->begin();};
    bound.interval.stepSize=[owner]{return owner->stepSize();};
    bound.interval.accept=[owner]{owner->accept();};
    if(fixedBaseRequired)bound.interval.snapshot=[owner]{owner->freezeBase();};
    bound.install=[owner,step=&call](Run::OpRegistry& ops) {
        const auto& assemble=step->fragment.children.at(0);
        const auto& advance=step->fragment.children.at(1);
        ops.bind(assemble.operation,assemble.provider,[owner]{owner->assemble();});
        ops.bind(advance.operation,advance.provider,[owner]{owner->advance();});
    };
    return bound;
}
}
