#include "SF_scalarBinding.h"
#include "solver/system/SF_scalarMethod.h"
#include "core/system/SF_operationIds.h"
#include "core/mesh/SF_dimension.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
namespace SF::SolverAlgorithm {
namespace {
/// Owns only numerical workspace. Physical writes use the realized STATE view.
class ScalarOperations {
    State::StateBundle& state_;
    Field& geometry_;
    const FDM::ScalarEquationData& data_;
    const System::CompiledTimeRecipe& recipe_;
    const double& limit_;
    double fixedDt_;
    const System::CompiledScalarEquation contract_;
    State::DistributedFieldView* physical_;
    std::array<int,3> shape_;
    std::array<bool,3> active_{};
    std::array<double,3> spacing_{};
    std::vector<double> old_,stage_,rhs_;
    std::array<std::vector<double>,4> increments_;
    State::DistributedFieldView oldView_,stageView_,rhsView_;
    System::StateRealization& realized_;
    std::shared_ptr<const FDM::ScalarEquationData> modelData_;
    std::weak_ptr<const FDM::ScalarEquationData> lifetime_;
    bool checkedLifetime_=false;
    void requireData() const {
        if (checkedLifetime_ && lifetime_.expired()) throw std::runtime_error("Scalar instance data lifetime expired.");
    }
    System::FormulaValues values_;
    double stageTime_=0;
    int nextStage_=-1;
    std::function<StageMassFluxView(StageFluxIdentity)> massFlux_;

    std::array<int,3> index(int cell) const {
        int i,j,k;geometry_.getIJK(cell,i,j,k);return {i,j,k};
    }
    double coordinate(int axis,const std::array<int,3>& at) const {
        return axis==0?geometry_.X(at[0],at[1],at[2])
            :axis==1?geometry_.Y(at[0],at[1],at[2]):geometry_.Z(at[0],at[1],at[2]);
    }
    double known(const FDM::ScalarKnownValue& value,int cell,double time) const {
        const auto at=index(cell);
        const double result=value(geometry_.X(at[0],at[1],at[2]),geometry_.Y(at[0],at[1],at[2]),geometry_.Z(at[0],at[1],at[2]),time);
        if (!std::isfinite(result)) throw std::runtime_error("Scalar known boundary/source value is non-finite.");
        return result;
    }
    template<class F> void owned(F&& f) const {
        const int g=geometry_.NG();
        for (int k=g;k<g+shape_[2];++k) for (int j=g;j<g+shape_[1];++j)
            for (int i=g;i<g+shape_[0];++i) f(geometry_.getIdx(i,j,k));
    }
    bool fixed(int cell) const {
        const auto at=index(cell);const int g=geometry_.NG();
        for (int axis=0;axis<3;++axis) if (active_[axis])
            for (int side=0;side<2;++side)
                if (at[axis]==g+(side?shape_[axis]-1:0)
                    && data_.boundary[2*axis+side].kind==FDM::ScalarBoundaryKind::FixedValue) return true;
        return false;
    }
    void boundary(double time) {
        const int g=geometry_.NG();
        // Boundary nodes are part of the nodal physical array. Dirichlet is
        // imposed at each true stage time; conflicting corner data are errors.
        owned([&](int cell) {
            const auto at=index(cell);bool assigned=false;double value=0;
            for (int axis=0;axis<3;++axis) if (active_[axis])
                for (int side=0;side<2;++side) {
                    const auto& bc=data_.boundary[2*axis+side];
                    if (at[axis]!=g+(side?shape_[axis]-1:0) || bc.kind!=FDM::ScalarBoundaryKind::FixedValue) continue;
                    const double next=known(bc.value,cell,time);
                    if (assigned && std::abs(value-next)>1e-12*std::max({1.0,std::abs(value),std::abs(next)}))
                        throw std::runtime_error("Conflicting scalar fixedValue corner boundaries.");
                    value=next;assigned=true;
                }
            if (assigned) stage_[cell]=value*(contract_.primitiveStorage?geometry_(at[0],at[1],at[2],0):1);
        });
        // Legacy scalar retains its nodal reflection. Shared transport publishes
        // concentration ghosts against the prepared NS density before Ready;
        // evaluation of any participant RHS then performs no scalar Stage writes.
        for (int axis=0;axis<3;++axis) if (active_[axis])
            owned([&](int cell) {
                auto at=index(cell);
                for (int side=0;side<2;++side) {
                    const int wall=g+(side?shape_[axis]-1:0);
                    if (at[axis]!=wall) continue;
                    auto inside=at,outside=at;inside[axis]+=side?-1:1;outside[axis]+=side?1:-1;
                    const double interior=stage_[geometry_.getIdx(inside[0],inside[1],inside[2])];
                    const auto& bc=data_.boundary[2*axis+side];
                    const int ghost=geometry_.getIdx(outside[0],outside[1],outside[2]);
                    if(contract_.sharedMassFlux) {
                        const double rho=geometry_(at[0],at[1],at[2],0),ri=geometry_(inside[0],inside[1],inside[2],0),rg=geometry_(outside[0],outside[1],outside[2],0);
                        if(!std::isfinite(rho) || !std::isfinite(ri) || !std::isfinite(rg) || rho<=0 || ri<=0 || rg<=0)
                            throw std::runtime_error("Primitive scalar boundary density must be positive and finite.");
                        const double wallConcentration=bc.kind==FDM::ScalarBoundaryKind::FixedValue?known(bc.value,cell,time)/(contract_.primitiveStorage?1:rho):0;
                        stage_[ghost]=(bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?stage_[cell]/rho:2*wallConcentration-interior/ri)*rg;
                    } else stage_[ghost]=bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?interior:2*stage_[cell]-interior;
                }
            });
    }
public:
    void ownBoundaryData(std::shared_ptr<const FDM::ScalarEquationData> data) {modelData_=std::move(data);}
    void trackData(std::weak_ptr<const FDM::ScalarEquationData> lifetime) {
        lifetime_=std::move(lifetime);checkedLifetime_=true;
    }
    ScalarOperations(const System::CompiledNumericalSystem& numerical,const System::CompiledEquationCall& call,
        State::StateBundle& state,const FDM::ScalarEquationData& data,System::StateRealization& realized,const double& limit,
        std::function<StageMassFluxView(StageFluxIdentity)> massFlux)
        :state_(state),geometry_(state.singlePatch()),data_(data),recipe_(numerical.time.recipe),limit_(limit),
         fixedDt_(numerical.dt.maxDeltaT),contract_(std::any_cast<System::CompiledScalarEquation>(call.providerContract)),
         physical_(realized.at(call.target.symbol).fields.at(0)),shape_{geometry_.NX(),geometry_.NY(),geometry_.NZ()},realized_(realized),massFlux_(std::move(massFlux)) {
        for(int axis=0;axis<3;++axis)active_[axis]=shape_[axis]>1 && (!contract_.sharedMassFlux || Math::isDirectionActiveIndex(axis));
        if(contract_.sharedMassFlux && !massFlux_)throw std::runtime_error("Shared scalar requires a bound canonical NS mass flux producer.");
        for (const auto& name:contract_.stageReads)
            if(std::none_of(call.stateUses.begin(),call.stateUses.end(),[&](const auto& use){return use.symbol==name && use.version==System::StateVersion::Stage;}))
                throw std::runtime_error("Scalar compiled cross read lacks declared Stage binding: "+name);
        if(physical_->geometry!=&geometry_ || !physical_->read || !physical_->write)
            throw std::runtime_error("Scalar storage geometry/read/write binding mismatch.");
        if (data.target!=contract_.target) throw std::runtime_error("Scalar boundary/source port targets a different STATE.");
        if (geometry_.NG()<1) throw std::runtime_error("Scalar Central2 requires one ghost layer.");
        if (physical_->components!=1) throw std::runtime_error("Scalar physical storage requires exactly one component.");
        const std::array<int,3> origin{geometry_.NG(),geometry_.NG(),geometry_.NG()};
        for (int axis=0;axis<3;++axis) if (active_[axis]) {
            if (shape_[axis]<3) throw std::runtime_error("Scalar active axis requires at least three nodes.");
            auto neighbor=origin;++neighbor[axis];
            spacing_[axis]=coordinate(axis,neighbor)-coordinate(axis,origin);
            if (!std::isfinite(spacing_[axis]) || spacing_[axis]<=0)
                throw std::runtime_error("Unsupported scalar geometry: positive axis-aligned uniform spacing required.");
            for (int side=0;side<2;++side) {
                const auto& bc=data.boundary[2*axis+side];
                if (bc.kind==FDM::ScalarBoundaryKind::Unspecified
                    || (bc.kind==FDM::ScalarBoundaryKind::FixedValue && !bc.value))
                    throw std::runtime_error("Missing scalar boundary contract on active axis.");
            }
        }
        if (spacing_==std::array<double,3>{}) throw std::runtime_error("Scalar diffusion requires an active spatial axis.");
        owned([&](int cell) {
            const auto at=index(cell);
            for (int axis=0;axis<3;++axis) if(active_[axis]) {
                const double expected=coordinate(axis,origin)+(at[axis]-origin[axis])*spacing_[axis];
                const double actual=coordinate(axis,at);
                if (!std::isfinite(actual) || std::abs(expected-actual)>1e-11*std::max({1.0,std::abs(expected),std::abs(actual)}))
                    throw std::runtime_error("Unsupported scalar geometry: nonuniform/nonorthogonal coordinates.");
            }
        });
        for (const auto& name:contract_.sources)
            if (!data.sources.count(name) || !data.sources.at(name))
                throw std::runtime_error("Missing scalar known-source binding: "+name);
        const auto total=geometry_.TotalSize();old_.resize(total);stage_.resize(total);rhs_.resize(total);
        for (int i=0;i<recipe_.stageCount();++i) increments_[i].resize(total);
        oldView_=State::workspaceView(call.oldTimeWorkspace,0,geometry_,old_,1);
        stageView_=State::workspaceView(call.stageWorkspace,0,geometry_,stage_,1);
        rhsView_=State::workspaceView(call.residualWorkspace,0,geometry_,rhs_,1);
        realized.bindView(contract_.target,System::StateViewKind::Workspace,rhsView_);
        if(contract_.primitiveStorage) {
            oldView_.read=[this](int cell,int) {
                const auto& density=realized_.view("rho",System::StateViewKind::OldTime);
                const double rho=density.fields.at(0)->read(cell,0);
                if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Primitive OldTime transport density invalid.");
                return old_[cell]/rho;
            };
            stageView_.read=[this](int cell,int) {int i,j,k;geometry_.getIJK(cell,i,j,k);double rho=geometry_(i,j,k,0);
                if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Primitive stage transport density invalid.");return stage_[cell]/rho;};
        }
        if (realized.requestsView(contract_.target,System::StateViewKind::OldTime))
            realized.bindView(contract_.target,System::StateViewKind::OldTime,oldView_);
        for (int i=0;i<recipe_.stageCount();++i)
            realized.bindView(contract_.target,System::StateViewKind::Stage,stageView_,i);
        values_.spacing=[this](int axis) {return spacing_.at(axis);};
        values_.neighbor=[this](int cell,int axis,int offset) {
            auto at=index(cell);at.at(axis)+=offset;return stage_.at(geometry_.getIdx(at[0],at[1],at[2]));
        };
    }
    void bindInputs() {
        for (const auto& name:contract_.stageReads)
            values_.boundReads.push_back(realized_.stageReader(name,recipe_.stageCount()));
        values_.boundNeighbor=[this](std::size_t slot,int cell,int axis,int offset,int component) {
            auto at=index(cell);at.at(axis)+=offset;
            return values_.boundReads.at(slot)(geometry_.getIdx(at[0],at[1],at[2]),component);
        };
        for (const auto& name:contract_.sources) {
            auto callback=data_.sources.at(name);
            values_.boundReads.push_back([this,callback=std::move(callback)](int cell,int) {return known(callback,cell,stageTime_);});
        }
    }
    double stepSize() const {
        requireData();
        if (!std::isfinite(fixedDt_) || fixedDt_<=0 || !std::isfinite(limit_) || limit_<=0)
            throw std::runtime_error("Scalar fixed dt/driver limit must be finite and positive.");
        return fixedDt_;
    }
    void begin() {
        requireData();
        for (int cell=0;cell<geometry_.TotalSize();++cell) {
            int i,j,k;geometry_.getIJK(cell,i,j,k);
            stage_[cell]=physical_->read(cell,0)*(contract_.primitiveStorage?geometry_(i,j,k,0):1);
            if (!std::isfinite(stage_[cell])) throw std::runtime_error("Scalar physical STATE is non-finite.");
        }
        boundary(state_.time);old_=stage_;nextStage_=0;
    }
    void prepareStage(const Run::ExecutionContext& context) {
        requireData();
        recipe_.requireProviderStages(context.stageCount);
        const int n=context.stageIndex;
        if (n!=nextStage_ || n<0 || n>=recipe_.stageCount()) throw std::runtime_error("Scalar stage order does not match frozen recipe.");
        stageTime_=state_.time+recipe_.stage(n).abscissa*state_.dt;
        boundary(stageTime_);
    }
    void evaluateRHS(const Run::ExecutionContext& context) {
        requireData();
        const int n=context.stageIndex;
        std::optional<StageMassFluxView> shared;
        if(contract_.sharedMassFlux) {
            shared.emplace(massFlux_({state_.step,n,stageTime_}));

        }
        const auto advection=[&](int cell) {
            if(!shared)return Discretization::conservativeScalarAdvection(values_,cell);
            double divergence=0;auto center=index(cell);
            const auto concentration=[&](const std::array<int,3>& at) {
                int node=geometry_.getIdx(at[0],at[1],at[2]);
                double rho=values_.boundReads.at(1)(node,0);
                if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Shared scalar Stage density must be positive.");
                return values_.boundReads.at(0)(node,0)/(contract_.primitiveStorage?1:rho);
            };
            for(int axis=0;axis<3;++axis)if(active_[axis]) {
                for(int sign:{-1,1}) {
                    auto left=center,right=center;
                    if(sign<0)--left[axis];else ++right[axis];
                    double mass=shared->at(axis,left[0],left[1],left[2]);
                    divergence+=sign*mass*concentration(mass>=0?left:right);
                }
            }
            return geometry_.Jac(center[0],center[1],center[2])*divergence;
        };
        owned([&](int cell) {
            rhs_[cell]=fixed(cell)?0:contract_.rhs(values_,cell,0)
                -(contract_.conservativeTransport?advection(cell):0);
            if (!std::isfinite(rhs_[cell])) throw std::runtime_error("Scalar RHS is non-finite.");
            increments_[n][cell]=rhs_[cell];
        });
    }
    void advanceStage(const Run::ExecutionContext& context) {
        const int n=context.stageIndex;
        if (n+1<recipe_.stageCount() || recipe_.stageCount()==1) {
            owned([&](int cell) {stage_[cell]=old_[cell]+recipe_.stage(n).incrementWeight*state_.dt*rhs_[cell];});
        } else {
            const auto weights=recipe_.finalWeights();
            owned([&](int cell) {
                double increment=0;
                for (int i=0;i<recipe_.stageCount();++i) increment+=weights[i]*increments_[i][cell];
                stage_[cell]=old_[cell]+state_.dt*increment/recipe_.finalDivisor();
            });
        }
        ++nextStage_;
    }
    void validatePublish() {
        requireData();
        if (nextStage_!=recipe_.stageCount()) throw std::runtime_error("Scalar cannot commit incomplete stages.");
        boundary(state_.time+state_.dt);
        for (int cell=0;cell<geometry_.TotalSize();++cell) {
            if (!std::isfinite(stage_[cell])) throw std::runtime_error("Scalar update is non-finite.");
        }
        if(contract_.primitiveStorage)owned([&](int cell){if(stage_[cell]<=0)throw std::runtime_error("Native stage scalar physical update nonpositive; no repair.");});
    }
    void publish() {
        for (int cell=0;cell<geometry_.TotalSize();++cell) {
            int i,j,k;geometry_.getIJK(cell,i,j,k);
            const double value=stage_[cell]/(contract_.primitiveStorage?geometry_(i,j,k,0):1);
            if(contract_.primitiveStorage && (!std::isfinite(value)))throw std::runtime_error("Stage turbulence update is nonpositive; no repair.");
            physical_->write(cell,0,value);
        }
        state_.distributed.markModified(physical_->name);
        nextStage_=-1;
    }
};
class StageTransportClosure {
    State::StateBundle& state_;
    Field& field_;
    FDM::ITransportModel& model_;
    System::StateRealization& realized_;
    int stages_;
    std::vector<std::string> targets_;
    std::array<std::vector<double>,7> values_;
    std::vector<double> viscosity_;
    std::array<State::DistributedFieldView,7> views_;
    std::array<std::function<double(int,int)>,2> reads_;
public:
    StageTransportClosure(State::StateBundle& state,Field& field,FDM::ITransportModel& model,System::StateRealization& realized,int stages,
            std::vector<std::string> targets):state_(state),field_(field),model_(model),realized_(realized),stages_(stages),targets_(std::move(targets)) {
        if(model_.stageScalarTargets()!=targets_)throw std::runtime_error("Stage scalar closure port target/model mismatch.");
        const std::array<std::string,7> names={"SST.nu.k","SST.nu.omega","SST.production.k","SST.production.omega","SST.sink.k","SST.sink.omega","SST.cross"};
        viscosity_.resize(field.TotalSize());
        for(int stage=0;stage<stages_;++stage)if(realized_.requestsView("mu_t",System::StateViewKind::Stage,stage))realized_.bindView("mu_t",System::StateViewKind::Stage,*realized_.at("mu_t").fields.at(0),stage);
        for(int slot=0;slot<7;++slot) {
            values_[slot].resize(field.TotalSize());views_[slot]=State::workspaceView(names[slot],0,field,values_[slot],1);
            for(int stage=0;stage<stages_;++stage)if(realized_.requestsView(names[slot],System::StateViewKind::Stage,stage))realized_.bindView(names[slot],System::StateViewKind::Stage,views_[slot],stage);
        }
    }
    void bindInputs() {for(int s=0;s<2;++s)reads_[s]=realized_.stageReader(targets_[s],stages_);}
    void refresh(bool physical) {
        FDM::TransportStageSamples samples;
        for(int s=0;s<2;++s) {
            const auto read=physical?realized_.at(targets_[s]).fields.at(0)->read:reads_[s];
            samples[s]=[this,read](int i,int j,int k){const int g=field_.NG();
                // Active axes consume the declared boundary ghost. Inactive axes
                // have no scalar gradient; their sole physical node is aliased.
                const int ii=field_.NX()>1 && Math::isDirectionActiveIndex(0)?std::clamp(i,g-1,g+field_.NX()):g;
                const int jj=field_.NY()>1 && Math::isDirectionActiveIndex(1)?std::clamp(j,g-1,g+field_.NY()):g;
                const int kk=field_.NZ()>1 && Math::isDirectionActiveIndex(2)?std::clamp(k,g-1,g+field_.NZ()):g;
                return read(field_.getIdx(ii,jj,kk),0);};
        }
        int g=field_.NG();
        for(int k=g;k<g+field_.NZ();++k)for(int j=g;j<g+field_.NY();++j)for(int i=g;i<g+field_.NX();++i) {
            int cell=field_.getIdx(i,j,k);const auto c=model_.stageCoefficients(field_,samples,i,j,k);
            std::array<double,7> value={c.diffusivity[0],c.diffusivity[1],c.production[0],c.production[1],c.sink[0],c.sink[1],c.crossSource};
            for(int s=0;s<7;++s){if(!std::isfinite(value[s]))throw std::runtime_error("Stage coefficient/source is non-finite.");values_[s][cell]=value[s];}
            viscosity_[cell]=c.dynamicViscosity;
        }
        for(int cell=0;cell<field_.TotalSize();++cell) {
            int i,j,k;field_.getIJK(cell,i,j,k);int interior=field_.getIdx(std::clamp(i,g,g+field_.NX()-1),std::clamp(j,g,g+field_.NY()-1),std::clamp(k,g,g+field_.NZ()-1));
            for(int s=0;s<7;++s)values_[s][cell]=values_[s][interior];viscosity_[cell]=viscosity_[interior];
        }
        model_.publishStageViscosity(field_,viscosity_);
        state_.distributed.markModified(realized_.at("mu_t").fields.at(0)->name);
    }
};

}
std::vector<Time::TemporalCallbacks> scalarParticipants(const System::CompiledNumericalSystem& numerical,
        const System::CompiledSolvePlan& plan,State::StateBundle& state,const FDM::SolverServices& services,
        System::StateRealization& realized,const double& limit,
        std::function<StageMassFluxView(StageFluxIdentity)> massFlux) {
    if (state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
        throw std::runtime_error("Unsupported scalar execution: serial single patch required.");
    std::vector<const System::CompiledEquationCall*> calls;
    for (const auto& call:plan.compiledProgram.steps) if (call.backendProvider==System::ScalarOps::Provider || call.backendProvider==System::ScalarOps::TransportProvider || call.backendProvider==System::ScalarOps::SharedFluxProvider || call.backendProvider==System::ScalarOps::StageTransportProvider) calls.push_back(&call);
    if (services.scalarEquation && (!services.scalarInstances.empty() || calls.size()!=1))
        throw std::runtime_error("Single scalar data port cannot bind multiple occurrences or coexist with instance table.");
    for (std::size_t i=0;i<services.scalarInstances.size();++i) {
        const auto& item=services.scalarInstances[i];
        if (item.data.expired() || std::count_if(services.scalarInstances.begin(),services.scalarInstances.end(),
                [&](const auto& entry){return entry.occurrence==item.occurrence;})!=1
            || std::none_of(calls.begin(),calls.end(),[&](const auto* call){return call->source.occurrence==item.occurrence;}))
            throw std::runtime_error("Invalid/duplicate/unconsumed scalar instance data key.");
    }
    std::shared_ptr<StageTransportClosure> closure;
    const auto native=std::count_if(calls.begin(),calls.end(),[](const auto* call){return call->backendProvider==System::ScalarOps::StageTransportProvider;});
    if(native) {
        if(native!=2 || !services.transportModel || services.immersed.boundary || services.immersed.constraint || services.immersed.system)throw std::runtime_error("Native stage turbulence requires complete pair, stage coefficient port, and no immersed stage treatment.");
        auto contract=std::any_cast<System::CompiledScalarEquation>((*std::find_if(calls.begin(),calls.end(),[](const auto* c){return c->backendProvider==System::ScalarOps::StageTransportProvider;}))->providerContract);
        closure=std::make_shared<StageTransportClosure>(state,state.singlePatch(),*services.transportModel,realized,numerical.time.recipe.stageCount(),contract.modelTargets);
    }
    std::vector<Time::TemporalCallbacks> callbacks;
    std::vector<std::shared_ptr<ScalarOperations>> owners;
    for (const auto* call:calls) {
        const auto checked=[call,&state](const char* phase,int stage,auto&& operation)->decltype(auto) {
            try {return operation();}
            catch(const std::exception& error) {
                throw std::runtime_error("Scalar operation="+std::string(phase)+" occurrence="+call->source.occurrence
                    +" equation="+call->source.equation+" target="+call->target.symbol+" provider="+call->backendProvider
                    +" step="+std::to_string(state.step)+" stage="+std::to_string(stage)
                    +" patch=0: "+error.what());
            }
        };
        const auto found=std::find_if(services.scalarInstances.begin(),services.scalarInstances.end(),
            [&](const auto& item){return item.occurrence==call->source.occurrence;});
        const auto locked=found==services.scalarInstances.end()?std::shared_ptr<const FDM::ScalarEquationData>{}:found->data.lock();
        std::shared_ptr<FDM::ScalarEquationData> modelData;
        if(call->backendProvider==System::ScalarOps::StageTransportProvider) {
            if(locked || services.scalarEquation)throw std::runtime_error("Stage model boundary port cannot coexist with duplicate scalar data.");
            modelData=checked("boundary binding",-1,[&]{return std::make_shared<FDM::ScalarEquationData>(services.transportModel->stageScalarData(state.singlePatch(),call->target.symbol));});
        }
        const auto* data=modelData?modelData.get():services.scalarEquation?services.scalarEquation:locked.get();
        if (!data) throw std::runtime_error("Missing scalar instance data: "+call->source.occurrence);
        auto owner=checked("STATE binding",-1,[&]{return std::make_shared<ScalarOperations>(numerical,*call,state,*data,realized,limit,massFlux);});
        if (locked) owner->trackData(found->data);
        if(modelData)owner->ownBoundaryData(modelData);
        owners.push_back(owner);
        Time::TemporalCallbacks cb;cb.identity=call->source.occurrence;
        if(!call->coefficientPrepare.empty()){
            cb.coefficients=[closure,checked](const auto& c){checked("coefficients",c.stageIndex,[&]{closure->refresh(false);});};
            cb.publishCoefficients=[closure,checked]{checked("coefficient publication",-1,[&]{closure->refresh(true);});};
        }
        cb.stepSize=[owner,checked]{return checked("dt",-1,[&]{return owner->stepSize();});};
        cb.snapshot=[owner,checked]{checked("snapshot",-1,[&]{owner->begin();});};
        cb.prepare=[owner,checked](const auto& c){checked("prepare",c.stageIndex,[&]{owner->prepareStage(c);});};
        cb.rhs=[owner,checked](const auto& c){checked("rhs",c.stageIndex,[&]{owner->evaluateRHS(c);});};
        cb.advance=[owner,checked](const auto& c){checked("advance",c.stageIndex,[&]{owner->advanceStage(c);});};
        cb.validatePublish=[owner,checked]{checked("validate publish",-1,[&]{owner->validatePublish();});};
        cb.publish=[owner,checked]{checked("publish",-1,[&]{owner->publish();});};
        callbacks.push_back(std::move(cb));
    }
    for (const auto& owner:owners) owner->bindInputs();
    if(closure)closure->bindInputs();
    return callbacks;
}
}
