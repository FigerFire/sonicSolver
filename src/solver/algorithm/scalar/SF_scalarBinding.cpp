#include "SF_scalarBinding.h"
#include "solver/system/SF_executionContract.h"
#include "solver/system/SF_scalarMethod.h"
#include "core/system/SF_operationIds.h"
#include "core/mesh/SF_dimension.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
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
    std::weak_ptr<const FDM::ScalarEquationData> lifetime_;
    bool checkedLifetime_=false,inputsBound_=false;
    void requireData() const {
        if(!inputsBound_)throw std::runtime_error("Scalar input readers have not completed startup binding");
        if (checkedLifetime_ && lifetime_.expired()) throw std::runtime_error("Scalar instance data lifetime expired.");
    }
    System::FormulaValues values_;
    double stageTime_=0;
    int nextStage_=-1;
    std::function<StageMassFluxView(StageFluxIdentity)> massFlux_;

    int entry(int cell,int component) const {return component*geometry_.TotalSize()+cell;}
    template<class F> void ownedValues(F&& f) const {
        for(int component=0;component<contract_.components;++component)
            owned([&](int cell){f(cell,component);});
    }
    double boundaryValue(const FDM::ScalarBoundaryCondition& bc,int cell,double time,int component) const {
        if(contract_.components==1)return known(bc.value,cell,time);
        const auto at=index(cell);
        const double v=bc.componentValue(geometry_.X(at[0],at[1],at[2]),geometry_.Y(at[0],at[1],at[2]),geometry_.Z(at[0],at[1],at[2]),time,component);
        if(!std::isfinite(v))throw std::runtime_error("Non-finite component boundary value.");
        return v;
    }
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
        ownedValues([&](int cell,int component) {
            const auto at=index(cell);bool assigned=false;double value=0;
            for (int axis=0;axis<3;++axis) if (active_[axis])
                for (int side=0;side<2;++side) {
                    const auto& bc=data_.boundary[2*axis+side];
                    if (at[axis]!=g+(side?shape_[axis]-1:0) || bc.kind!=FDM::ScalarBoundaryKind::FixedValue) continue;
                    const double next=boundaryValue(bc,cell,time,component);
                    if (assigned && std::abs(value-next)>1e-12*std::max({1.0,std::abs(value),std::abs(next)}))
                        throw std::runtime_error("Conflicting scalar fixedValue corner boundaries.");
                    value=next;assigned=true;
                }
            if (assigned) stage_[entry(cell,component)]=value*(contract_.primitiveStorage?geometry_(at[0],at[1],at[2],0):1);
        });
        // Legacy scalar retains its nodal reflection. Shared transport publishes
        // concentration ghosts against the prepared NS density before Ready;
        // evaluation of any participant RHS then performs no scalar Stage writes.
        for (int axis=0;axis<3;++axis) if (active_[axis])
            ownedValues([&](int cell,int component) {
                auto at=index(cell);
                for (int side=0;side<2;++side) {
                    const int wall=g+(side?shape_[axis]-1:0);
                    if (at[axis]!=wall) continue;
                    auto inside=at,outside=at;inside[axis]+=side?-1:1;outside[axis]+=side?1:-1;
                    const double interior=stage_[entry(geometry_.getIdx(inside[0],inside[1],inside[2]),component)];
                    const auto& bc=data_.boundary[2*axis+side];
                    const int ghost=geometry_.getIdx(outside[0],outside[1],outside[2]);
                    if(contract_.sharedMassFlux) {
                        const double rho=geometry_(at[0],at[1],at[2],0),ri=geometry_(inside[0],inside[1],inside[2],0),rg=geometry_(outside[0],outside[1],outside[2],0);
                        if(!std::isfinite(rho) || !std::isfinite(ri) || !std::isfinite(rg) || rho<=0 || ri<=0 || rg<=0)
                            throw std::runtime_error("Primitive scalar boundary density must be positive and finite.");
                        const double wallConcentration=bc.kind==FDM::ScalarBoundaryKind::FixedValue?boundaryValue(bc,cell,time,component)/(contract_.primitiveStorage?1:rho):0;
                        stage_[entry(ghost,component)]=(bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?stage_[entry(cell,component)]/rho:2*wallConcentration-interior/ri)*rg;
                    } else stage_[entry(ghost,component)]=bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?interior:2*stage_[entry(cell,component)]-interior;
                }
            });
    }
public:
    void trackData(std::weak_ptr<const FDM::ScalarEquationData> lifetime) {
        lifetime_=std::move(lifetime);checkedLifetime_=true;
    }
    ScalarOperations(const System::CompiledNumericalSystem& numerical,const System::CompiledEquationCall& call,
        State::StateBundle& state,const FDM::ScalarEquationData& data,System::StateRealization& realized,const double& limit,
        std::function<StageMassFluxView(StageFluxIdentity)> massFlux,bool stencil)
        :state_(state),geometry_(state.singlePatch()),data_(data),recipe_(numerical.time.recipe),limit_(limit),
         fixedDt_(numerical.dt.maxDeltaT),contract_(std::any_cast<System::CompiledScalarEquation>(call.providerContract)),
         physical_(realized.at(call.target.symbol).fields.at(0)),shape_{geometry_.NX(),geometry_.NY(),geometry_.NZ()},realized_(realized),massFlux_(std::move(massFlux)) {
        for(int axis=0;axis<3;++axis)active_[axis]=stencil && shape_[axis]>1 && (!contract_.sharedMassFlux || Math::isDirectionActiveIndex(axis));
        if(contract_.sharedMassFlux && !massFlux_)throw std::runtime_error("Shared scalar requires a bound canonical NS mass flux producer.");
        for (const auto& read:contract_.valueReads)
            if(std::none_of(call.stateUses.begin(),call.stateUses.end(),[&](const auto& use){return use.symbol==read.symbol && use.version==read.version;}))
                throw std::runtime_error("Scalar compiled input lacks declared version binding: "+read.symbol);
        if(physical_->geometry!=&geometry_ || !physical_->read || !physical_->write)
            throw std::runtime_error("Scalar storage geometry/read/write binding mismatch.");
        if (data.target!=contract_.target) throw std::runtime_error("Scalar boundary/source port targets a different STATE.");
        if (stencil && geometry_.NG()<1) throw std::runtime_error("Scalar Central2 requires one ghost layer.");
        if (physical_->components!=contract_.components) throw std::runtime_error("Physical storage component count differs from compiled target.");
        const std::array<int,3> origin{geometry_.NG(),geometry_.NG(),geometry_.NG()};
        for (int axis=0;axis<3;++axis) if (active_[axis]) {
            if (shape_[axis]<3) throw std::runtime_error("Scalar active axis requires at least three nodes.");
            auto neighbor=origin;++neighbor[axis];
            spacing_[axis]=coordinate(axis,neighbor)-coordinate(axis,origin);
            if (!std::isfinite(spacing_[axis]) || spacing_[axis]<=0)
                throw std::runtime_error("Unsupported scalar geometry: positive axis-aligned uniform spacing required.");
            for (int side=0;side<2;++side) {
                const auto& bc=data.boundary[2*axis+side];
                if(contract_.arithmeticFaces && bc.kind!=FDM::ScalarBoundaryKind::FixedValue && bc.kind!=FDM::ScalarBoundaryKind::ZeroGradient)
                    throw std::runtime_error("Unsupported arithmetic diffusion boundary: existing fixedValue/zeroGradient ghost required");
                if (bc.kind==FDM::ScalarBoundaryKind::Unspecified
                    || (bc.kind==FDM::ScalarBoundaryKind::FixedValue && (contract_.components==1?!bc.value:!bc.componentValue)))
                    throw std::runtime_error("Missing scalar boundary contract on active axis.");
            }
        }
        if (stencil && spacing_==std::array<double,3>{}) throw std::runtime_error("Scalar diffusion requires an active spatial axis.");
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
        const auto total=geometry_.TotalSize()*contract_.components;old_.resize(total);stage_.resize(total);rhs_.resize(total);
        for (int i=0;i<recipe_.stageCount();++i) increments_[i].resize(total);
        if(std::getenv("SF_SCALAR_TRACE")) {
            std::ostringstream trace;trace<<"[SF SCALAR storage] target="<<contract_.target<<" old="<<static_cast<const void*>(old_.data())<<" stage="<<static_cast<const void*>(stage_.data())<<" rhs="<<static_cast<const void*>(rhs_.data());
            std::cout<<trace.str()<<'\n';
        }
        oldView_=State::workspaceView(call.oldTimeWorkspace,0,geometry_,old_,contract_.components);
        stageView_=State::workspaceView(call.stageWorkspace,0,geometry_,stage_,contract_.components);
        rhsView_=State::workspaceView(call.residualWorkspace,0,geometry_,rhs_,contract_.components);
        realized.bindView(contract_.target,System::StateViewKind::Workspace,rhsView_);
        if(contract_.primitiveStorage) {
            oldView_.read=[this](int cell,int) {
                const auto& density=realized_.view("rho",System::StateViewKind::OldTime);
                const double rho=density.fields.at(0)->read(cell,0);
                if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Primitive OldTime transport density invalid.");
                return old_[cell]/rho;
            };
            stageView_.read=[this](int cell,int) {int i,j,k;geometry_.getIJK(cell,i,j,k);double rho=geometry_(i,j,k,0);
                if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Capacity stage transport density invalid.");return stage_[cell]/rho;};
        }
        if (realized.requestsView(contract_.target,System::StateViewKind::OldTime))
            realized.bindView(contract_.target,System::StateViewKind::OldTime,oldView_);
        for (int i=0;i<recipe_.stageCount();++i)
            realized.bindView(contract_.target,System::StateViewKind::Stage,stageView_,i);
        values_.spacing=[this](int axis) {return spacing_.at(axis);};

    }
    void bindInputs() {
        if(inputsBound_)throw std::runtime_error("Scalar input readers were already bound");
        for (const auto& input:contract_.valueReads) {
            if(input.version==System::StateVersion::Stage)values_.boundReads.push_back(realized_.stageReader(input.symbol,recipe_.stageCount()));
            else if(input.version==System::StateVersion::OldTime)values_.boundReads.push_back(realized_.oldTimeReader(input.symbol));
            else if(input.version==System::StateVersion::Current) {
                const auto& physical=realized_.at(input.symbol);if(physical.fields.size()!=1)throw std::runtime_error("Missing Current input resource: "+input.symbol);
                const auto view=*physical.fields.at(0);const int offset=physical.descriptor->componentOffset;
                values_.boundReads.push_back([view,offset](int cell,int component){if(view.requireReadable)view.requireReadable();return view.read(cell,offset+component);});
            } else throw std::runtime_error("Unsupported runtime coefficient version: "+input.symbol);
        }
        values_.boundNeighbor=[this](std::size_t slot,int cell,int axis,int offset,int component) {
            auto at=index(cell);at.at(axis)+=offset;
            return values_.boundReads.at(slot)(geometry_.getIdx(at[0],at[1],at[2]),component);
        };
        values_.coefficientNeighbor=[this](std::size_t slot,int cell,int axis,int offset) {
            auto at=index(cell);at.at(axis)+=offset;const int lo=geometry_.NG(),hi=lo+shape_.at(axis)-1;
            // Match the existing nodal scalar ghost reflection; the coefficient owner is read, never extended or overwritten.
            if(at[axis]<lo)at[axis]=2*lo-at[axis];else if(at[axis]>hi)at[axis]=2*hi-at[axis];
            return values_.boundReads.at(slot)(geometry_.getIdx(at[0],at[1],at[2]),0);
        };
        for (const auto& name:contract_.sources) {
            auto callback=data_.sources.at(name);
            values_.boundReads.push_back([this,callback=std::move(callback)](int cell,int) {return known(callback,cell,stageTime_);});
        }
        inputsBound_=true;
    }
    double stepSize() const {
        requireData();
        if (!std::isfinite(fixedDt_) || fixedDt_<=0 || !std::isfinite(limit_) || limit_<=0)
            throw std::runtime_error("Scalar fixed dt/driver limit must be finite and positive.");
        return fixedDt_;
    }
    void begin() {
        requireData();
        for(int component=0;component<contract_.components;++component) for (int cell=0;cell<geometry_.TotalSize();++cell) {
            int i,j,k;geometry_.getIJK(cell,i,j,k);
            stage_[entry(cell,component)]=physical_->read(cell,component)*(contract_.primitiveStorage?geometry_(i,j,k,0):1);
            if (!std::isfinite(stage_[entry(cell,component)])) throw std::runtime_error("Scalar physical STATE is non-finite.");
        }
        boundary(state_.time);old_=stage_;nextStage_=0;
    }
    void prepareStage(const Run::ExecutionContext& context) {
        requireData();
        recipe_.requireProviderStages(context.stageCount);
        const int n=context.stageIndex;
        if (n!=nextStage_ || n<0 || n>=recipe_.stageCount()) throw std::runtime_error("Scalar stage order does not match frozen recipe.");
        stageTime_=state_.time+recipe_.stage(n).abscissa*state_.dt;
        if(std::getenv("SF_SCALAR_TRACE")) {
            std::ostringstream trace;trace<<std::setprecision(17)<<"[SF SCALAR stage] target="<<contract_.target<<" step="<<state_.step<<" stage="<<n<<" time="<<stageTime_<<" dt="<<state_.dt;
            std::cout<<trace.str()<<'\n';
        }
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
        ownedValues([&](int cell,int component) {
            rhs_[entry(cell,component)]=fixed(cell)?0:contract_.rhs(values_,cell,component)
                -(contract_.conservativeTransport?advection(cell):0);
            if (!std::isfinite(rhs_[entry(cell,component)])) throw std::runtime_error("Scalar RHS is non-finite.");
            increments_[n][entry(cell,component)]=rhs_[entry(cell,component)];
        });
    }
    void advanceStage(const Run::ExecutionContext& context) {
        const int n=context.stageIndex;
        if (n+1<recipe_.stageCount() || recipe_.stageCount()==1) {
            ownedValues([&](int cell,int component) {stage_[entry(cell,component)]=old_[entry(cell,component)]+recipe_.stage(n).incrementWeight*state_.dt*rhs_[entry(cell,component)];});
        } else {
            const auto weights=recipe_.finalWeights();
            ownedValues([&](int cell,int component) {
                double increment=0;
                for (int i=0;i<recipe_.stageCount();++i) increment+=weights[i]*increments_[i][entry(cell,component)];
                stage_[entry(cell,component)]=old_[entry(cell,component)]+state_.dt*increment/recipe_.finalDivisor();
            });
        }
        ++nextStage_;
    }
    void validatePublish() {
        requireData();
        if (nextStage_!=recipe_.stageCount()) throw std::runtime_error("Scalar cannot commit incomplete stages.");
        boundary(state_.time+state_.dt);
        for(int component=0;component<contract_.components;++component) for (int cell=0;cell<geometry_.TotalSize();++cell) {
            if (!std::isfinite(stage_[entry(cell,component)])) throw std::runtime_error("Scalar update is non-finite.");
        }
        if(contract_.primitiveStorage)ownedValues([&](int cell,int component){if(stage_[entry(cell,component)]<=0)throw std::runtime_error("Native stage scalar physical update nonpositive; no repair.");});
    }
    void publish() {
        for(int component=0;component<contract_.components;++component) for (int cell=0;cell<geometry_.TotalSize();++cell) {
            int i,j,k;geometry_.getIJK(cell,i,j,k);
            const double value=stage_[entry(cell,component)]/(contract_.primitiveStorage?geometry_(i,j,k,0):1);
            if(contract_.primitiveStorage && (!std::isfinite(value)))throw std::runtime_error("Capacity transport update is nonpositive; no repair.");
            physical_->write(cell,component,value);
        }
        state_.distributed.markModified(physical_->name);
        nextStage_=-1;
    }
};
// This is the original final-coefficient extension, not evaluation on ghost inputs.
void extendNearestInterior(const Field& mesh,std::vector<double>& values) {
    const int g=mesh.NG();
    for(int cell=0;cell<mesh.TotalSize();++cell) {
        int i,j,k;mesh.getIJK(cell,i,j,k);
        const int interior=mesh.getIdx(std::clamp(i,g,g+mesh.NX()-1),
            std::clamp(j,g,g+mesh.NY()-1),std::clamp(k,g,g+mesh.NZ()-1));
        values.at(cell)=values.at(interior);
    }
}

class SstClosureOperations {
    State::StateBundle& state_;
    Field& field_;
    FDM::ITransportModel& model_;
    System::StateRealization& realized_;
    int stages_;
    double molecularViscosity_;
    std::vector<std::string> targets_;
    std::array<std::vector<double>,4> values_;
    std::vector<double> viscosity_;
    std::array<State::DistributedFieldView,4> physicalViews_;
    std::vector<State::DistributedFieldView> stageViews_;
    std::array<std::function<double(int,int)>,2> reads_;
    struct Publication {
        int step=-1,stage=-1;double time=0,dt=0;bool physical=false;
        void require(const State::StateBundle& state,bool current,const std::string& input,int index=-1) const {
            const bool committed=current && physical && step+1==state.step && time+dt==state.time;
            if(!committed && (step!=state.step || physical!=current || (current?time!=state.time:stage!=index)))
                throw std::runtime_error("Unpublished closure input="+input+" requested="+std::string(current?"Current":"Stage")+" step="+std::to_string(state.step)+" stage="+std::to_string(index));
        }
    };
    std::shared_ptr<Publication> publication_=std::make_shared<Publication>();
public:
    SstClosureOperations(State::StateBundle& state,Field& field,FDM::ITransportModel& model,System::StateRealization& realized,int stages,
            std::vector<std::string> targets,const System::CompiledEquationCall& publisher):state_(state),field_(field),model_(model),realized_(realized),stages_(stages),targets_(std::move(targets)) {
        if(model_.stageScalarTargets()!=targets_)throw std::runtime_error("Stage scalar closure port target/model mismatch.");
        const std::array<std::string,4> names={"SST.F1","SST.production.k","SST.production.omega","SST.cross"};
        std::vector<std::string> provided{"mu_t"};provided.insert(provided.end(),names.begin(),names.end());
        if(publisher.coefficientVersion!=System::StateVersion::Stage || publisher.coefficientStates.size()!=provided.size()
            || std::any_of(provided.begin(),provided.end(),[&](const auto& name){return std::count(publisher.coefficientStates.begin(),publisher.coefficientStates.end(),name)!=1;}))
            throw std::runtime_error("Stage coefficient publication set does not match the implemented SST workspace layout: "+publisher.source.occurrence);
        const std::vector<std::string> requiredInputs={"rho","rhoU",targets_.at(0),targets_.at(1)};
        if(publisher.coefficientInputs.size()!=requiredInputs.size() || std::any_of(requiredInputs.begin(),requiredInputs.end(),[&](const auto& name){
            return std::count_if(publisher.coefficientInputs.begin(),publisher.coefficientInputs.end(),[&](const auto& use){return use.symbol==name && use.version==System::StateVersion::Stage && use.kind==System::TargetKind::Physical;})!=1;
        }))throw std::runtime_error("Stage closure inputs differ from the implemented rho/momentum/transport sampling: "+publisher.source.occurrence);
        const auto mu=realized_.at("kOmegaSST.laminarMu").descriptor->constantValue;
        if(!mu || !std::isfinite(*mu) || *mu<0)throw std::runtime_error("Stage SST requires the actual configured molecular-viscosity constant");
        molecularViscosity_=*mu;
        viscosity_.resize(field.TotalSize());stageViews_.reserve(5*stages_);
        auto* physicalMu=realized_.at("mu_t").fields.at(0);
        if(physicalMu->geometry!=&field || physicalMu->components!=1 || !physicalMu->read)
            throw std::runtime_error("Stage SST mu_t backing/layout/lifetime is missing");
        const auto originalMu=*physicalMu;
        // Only the read guard is decorated; the original model array remains the owner.
        physicalMu->requireReadable=[original=originalMu.requireReadable,stamp=publication_,&state] {
            if(original)original();stamp->require(state,true,"mu_t");
        };
        physicalMu->read=[original=originalMu.read,guard=physicalMu->requireReadable](int cell,int component) {
            guard();return original(cell,component);
        };
        const auto stageView=[&](const std::string& name,const State::DistributedFieldView& backing,int stage) {
            auto view=backing;
            view.write=[name](int,int,double){throw std::runtime_error("Read-only prepared closure input="+name);};
            view.read=[this,backing,stage,name](int cell,int component) {
                realized_.requireStageRead(stage,publication_->time);
                publication_->require(state_,false,name,stage);
                if(backing.requireReadable)backing.requireReadable();
                return backing.read(cell,component);
            };
            stageViews_.push_back(std::move(view));
            realized_.bindView(name,System::StateViewKind::Stage,stageViews_.back(),stage);
        };
        for(int stage=0;stage<stages_;++stage)stageView("mu_t",originalMu,stage);
        for(int slot=0;slot<4;++slot) {
            values_[slot].resize(field.TotalSize());
            auto view=State::workspaceView(names[slot],0,field,values_[slot],1);
            // Stage aliases have their own guard; do not inherit the Current guard.
            for(int stage=0;stage<stages_;++stage)stageView(names[slot],view,stage);
            view.write=[name=names[slot]](int,int,double){throw std::runtime_error("Read-only prepared closure input="+name);};
            view.requireReadable=[stamp=publication_,&state,name=names[slot]]{stamp->require(state,true,name);};
            view.read=[original=view.read,guard=view.requireReadable](int cell,int component){guard();return original(cell,component);};
            physicalViews_[slot]=std::move(view);
            realized_.bindView(names[slot],System::StateViewKind::Physical,physicalViews_[slot]);
        }
    }
    void bindInputs() {
        for(int s=0;s<2;++s)reads_[s]=realized_.stageReader(targets_[s],stages_);
    }
    void refresh(bool physical,int stage=-1,double time=0) {
        if(!physical && publication_->step==state_.step && !publication_->physical && publication_->stage==stage)
            throw std::runtime_error("Duplicate nonlocal closure evaluation in one stage; iterative closure treatment is not implemented");
        publication_->step=-1;
        FDM::TransportStageSamples samples;
        for(int s=0;s<2;++s) {
            const auto read=physical?realized_.at(targets_[s]).fields.at(0)->read:reads_[s];
            samples[s]=[this,read](int i,int j,int k){const int g=field_.NG();
                // Keep the old active-axis stencil and inactive-axis alias exactly.
                const int ii=field_.NX()>1 && Math::isDirectionActiveIndex(0)?std::clamp(i,g-1,g+field_.NX()):g;
                const int jj=field_.NY()>1 && Math::isDirectionActiveIndex(1)?std::clamp(j,g-1,g+field_.NY()):g;
                const int kk=field_.NZ()>1 && Math::isDirectionActiveIndex(2)?std::clamp(k,g-1,g+field_.NZ()):g;
                return read(field_.getIdx(ii,jj,kk),0);};
        }
        const int g=field_.NG();
        for(int k=g;k<g+field_.NZ();++k)for(int j=g;j<g+field_.NY();++j)for(int i=g;i<g+field_.NX();++i) {
            const int cell=field_.getIdx(i,j,k);const auto c=model_.stageCoefficients(field_,samples,i,j,k);
            if(c.molecularViscosity!=molecularViscosity_)
                throw std::runtime_error("Stage SST molecular-viscosity backing disagrees with WHAT constant kOmegaSST.laminarMu");
            const std::array<double,4> value={c.blending1,c.production[0],c.production[1],c.crossSource};
            for(int s=0;s<4;++s){if(!std::isfinite(value[s]))throw std::runtime_error("Stage coefficient/source is non-finite.");values_[s][cell]=value[s];}
            if(!std::isfinite(c.dynamicViscosity))throw std::runtime_error("Stage mu_t is non-finite");
            viscosity_[cell]=c.dynamicViscosity;
        }
        for(auto& values:values_)extendNearestInterior(field_,values);
        extendNearestInterior(field_,viscosity_);
        model_.publishStageViscosity(field_,viscosity_);
        state_.distributed.markModified(realized_.at("mu_t").fields.at(0)->name);
        *publication_={state_.step,stage,physical?state_.time:time,state_.dt,physical};
    }
};


const std::string& fieldPort(const System::CompiledEquationCall& call) {
    const System::ResourceViewRequirement* port=nullptr;
    for(const auto& need:call.physicalResources)if((need.owner=="native.scalar" || need.owner=="model.scalar") && need.fields) {
        if(port)throw std::runtime_error("Multiple field ports for one ordinary equation: "+call.source.occurrence);
        port=&need;
    }
    if(!port)throw std::runtime_error("Missing compiled field port: "+call.source.occurrence);
    return port->id;
}

void validateScalarInstances(const System::CompiledSolvePlan& plan,const FDM::SolverServices& services) {
    for (std::size_t i=0;i<services.scalarInstances.size();++i) {
        const auto& item=services.scalarInstances[i];
        if(item.data.expired())throw std::runtime_error("Expired scalar resource owner: "+item.occurrence);
        if (std::count_if(services.scalarInstances.begin(),services.scalarInstances.end(),
                [&](const auto& entry){return entry.occurrence==item.occurrence;})!=1
            || std::none_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[&](const auto& call){return std::any_of(call.physicalResources.begin(),call.physicalResources.end(),[&](const auto& resource){return (resource.owner=="native.scalar" || resource.owner=="model.scalar") && resource.id==item.occurrence;});}))
            throw std::runtime_error("Invalid/duplicate/unconsumed scalar instance data key.");
    }
}

/// Owns matrix/solve workspace only. Current STATE remains in its original storage.
class BlockOperations {
    const System::CompiledEquationCall& call_;
    const System::CompiledLinearBlock contract_;
    State::StateBundle& state_;
    System::StateRealization& realized_;
    Field& mesh_;
    std::array<int,3> shape_;
    std::array<double,3> spacing_{};
    std::vector<int> cells_,local_;
    std::vector<State::DistributedFieldView*> physical_;
    std::vector<std::weak_ptr<const FDM::ScalarEquationData>> ports_;
    System::FormulaValues values_;
    LinearAlgebra::GlobalDofSystem matrix_;
    std::unique_ptr<LinearAlgebra::StaticDistributedNumbering> numbering_;
    std::unique_ptr<LinearAlgebra::DistributedLinearSystem> linear_;
    std::vector<double> pending_;
    bool assembled_=false,solved_=false,bound_=false;
    std::array<int,3> index(int cell)const{int i,j,k;mesh_.getIJK(cell,i,j,k);return {i,j,k};}
    std::shared_ptr<const FDM::ScalarEquationData> port(std::size_t slot)const{auto p=ports_.at(slot).lock();if(!p)throw std::runtime_error("Block resource owner expired: "+contract_.targets.at(slot));return p;}
    double known(const FDM::ScalarKnownValue& f,int cell)const{if(!f)throw std::runtime_error("Missing block boundary value");auto at=index(cell);double v=f(mesh_.X(at[0],at[1],at[2]),mesh_.Y(at[0],at[1],at[2]),mesh_.Z(at[0],at[1],at[2]),state_.time);if(!std::isfinite(v))throw std::runtime_error("Nonfinite block boundary");return v;}
    std::optional<double> fixed(std::size_t slot,int cell)const{
        auto data=port(slot);auto at=index(cell);std::optional<double> result;
        for(int axis=0;axis<3;++axis)if(spacing_[axis])for(int side=0;side<2;++side)if(at[axis]==mesh_.NG()+(side?shape_[axis]-1:0) && data->boundary[2*axis+side].kind==FDM::ScalarBoundaryKind::FixedValue){double next=known(data->boundary[2*axis+side].value,cell);if(result && std::abs(*result-next)>1e-12*std::max({1.,std::abs(*result),std::abs(next)}))throw std::runtime_error("Conflicting block fixedValue corner");result=next;}
        return result;
    }
public:
    BlockOperations(const System::CompiledEquationCall& call,State::StateBundle& state,System::StateRealization& realized,const FDM::SolverServices& services)
        :call_(call),contract_(std::any_cast<System::CompiledLinearBlock>(call.providerContract)),state_(state),realized_(realized),mesh_(state.singlePatch()),shape_{mesh_.NX(),mesh_.NY(),mesh_.NZ()} {
        if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))throw std::runtime_error("Unsupported linear block backend: serial single patch required");
        int g=mesh_.NG();local_.assign(mesh_.TotalSize(),-1);
        for(int k=g;k<g+shape_[2];++k)for(int j=g;j<g+shape_[1];++j)for(int i=g;i<g+shape_[0];++i){int cell=mesh_.getIdx(i,j,k);local_[cell]=cells_.size();cells_.push_back(cell);}
        const auto coordinate=[&](int axis,int cell){auto at=index(cell);return axis==0?mesh_.X(at[0],at[1],at[2]):axis==1?mesh_.Y(at[0],at[1],at[2]):mesh_.Z(at[0],at[1],at[2]);};
        int origin=cells_.front();auto start=index(origin);
        for(int axis=0;axis<3;++axis)if(shape_[axis]>1){if(shape_[axis]<3)throw std::runtime_error("Block Central2 requires at least three nodes per active direction");auto next=start;++next[axis];double h=coordinate(axis,mesh_.getIdx(next[0],next[1],next[2]))-coordinate(axis,origin);if(!std::isfinite(h) || h<=0)throw std::runtime_error("Unsupported block Cartesian spacing");spacing_[axis]=h;for(int cell:cells_){double expected=coordinate(axis,origin)+(index(cell)[axis]-g)*h;if(std::abs(coordinate(axis,cell)-expected)>1e-11*std::max(1.,std::abs(expected)))throw std::runtime_error("Unsupported nonuniform/nonorthogonal block geometry");}}
        for(std::size_t n=0;n<contract_.targets.size();++n){const auto& target=contract_.targets[n];auto& field=realized_.at(target);if(field.fields.size()!=1)throw std::runtime_error("Unbound block target: "+target);auto* view=field.fields[0];if(view->geometry!=&mesh_ || view->components!=1 || !view->read || !view->write || (field.descriptor->storageBinding==System::StorageBinding::ProviderDistributed && !view->requireReadable))throw std::runtime_error("Block target layout/access/lifetime capability missing: "+target);physical_.push_back(view);const auto& id=call_.physicalResources.at(n).id;auto instance=std::find_if(services.scalarInstances.begin(),services.scalarInstances.end(),[&](const auto& p){return p.occurrence==id;});if(instance==services.scalarInstances.end())throw std::runtime_error("Missing block boundary resource: "+id);ports_.push_back(instance->data);auto data=port(n);if(data->target!=target)throw std::runtime_error("Block boundary port target mismatch");for(int axis=0;axis<3;++axis)if(spacing_[axis])for(int side=0;side<2;++side)if(data->boundary[2*axis+side].kind!=FDM::ScalarBoundaryKind::FixedValue && data->boundary[2*axis+side].kind!=FDM::ScalarBoundaryKind::ZeroGradient)throw std::runtime_error("Unsupported block boundary: existing fixedValue/zeroGradient required");}
        values_.spacing=[this](int axis){return spacing_.at(axis);};
        values_.neighborCell=[this](int cell,int axis,int sign){auto at=index(cell);at[axis]+=sign;int lo=mesh_.NG(),hi=lo+shape_.at(axis)-1;if(at[axis]<lo)at[axis]=2*lo-at[axis];else if(at[axis]>hi)at[axis]=2*hi-at[axis];int next=mesh_.getIdx(at[0],at[1],at[2]);if(local_.at(next)<0)throw std::runtime_error("Block stencil escaped unhandled boundary");return next;};
        values_.dof=[this](const std::string& target,int cell,int component){auto n=std::find(contract_.targets.begin(),contract_.targets.end(),target);if(n==contract_.targets.end() || component!=0 || local_.at(cell)<0)throw std::runtime_error("Missing block column unknown/layout");return LinearAlgebra::GlobalDofId::make(LinearAlgebra::GlobalDofSpace::ScalarTransport,(n-contract_.targets.begin())*cells_.size()+local_[cell]);};
        std::vector<std::pair<LinearAlgebra::GlobalDofId,std::int64_t>> entries;for(const auto& target:contract_.targets)for(int cell:cells_){auto dof=values_.dof(target,cell,0);entries.push_back({dof,static_cast<std::int64_t>(entries.size())});}
        numbering_=std::make_unique<LinearAlgebra::StaticDistributedNumbering>(0,entries.size()-1,entries.size(),entries);linear_=std::make_unique<LinearAlgebra::DistributedLinearSystem>(contract_.linear,*numbering_);
    }
    void bindReads(){if(bound_)throw std::runtime_error("Block input readers already bound");for(const auto& read:contract_.valueReads){const auto& input=realized_.at(read.symbol);if(input.fields.size()!=1)throw std::runtime_error("Missing Current block input resource: "+read.symbol);auto view=*input.fields[0];if(view.geometry!=&mesh_ || view.components!=1 || !view.read || (input.descriptor->storageBinding==System::StorageBinding::ProviderDistributed && !view.requireReadable))throw std::runtime_error("Block input layout/lifetime capability missing: "+read.symbol);values_.boundReads.push_back([view](int cell,int component){if(view.requireReadable)view.requireReadable();double v=view.read(cell,component);if(!std::isfinite(v))throw std::runtime_error("Nonfinite block STATE read");return v;});}bound_=true;}
    void assemble(){
        assembled_=solved_=false;matrix_.rows.clear();pending_.clear();if(!bound_)throw std::runtime_error("Unbound block inputs");for(const auto& view:physical_)if(view->requireReadable)view->requireReadable();
        for(std::size_t target=0;target<contract_.targets.size();++target){auto row=std::find(contract_.rowTargets.begin(),contract_.rowTargets.end(),contract_.targets[target]);if(row==contract_.rowTargets.end())throw std::runtime_error("Missing block equation row");std::size_t eq=row-contract_.rowTargets.begin();
            // Principal diffusion identifies the physical boundary row; HOW list order only numbers columns.
            for(int cell:cells_){auto dof=values_.dof(contract_.targets[target],cell,0);LinearAlgebra::GlobalDofRow output(dof);if(auto prescribed=fixed(target,cell)){output.add(dof,1);output.setRightHandSide(*prescribed);}else{auto term=contract_.rows.at(eq)(values_,cell,0);if(!std::isfinite(term.constant))throw std::runtime_error("Nonfinite block RHS");for(const auto& a:term.coefficients)output.add(a.column,a.value);output.setRightHandSide(-term.constant);}double initial=physical_[target]->read(cell,0);if(!std::isfinite(initial))throw std::runtime_error("Nonfinite block initial target");output.setInitialGuess(initial);matrix_.rows.push_back(std::move(output));}
        }
        if(std::getenv("SF_BLOCK_TRACE"))for(const auto& row:matrix_.rows){std::cout<<std::setprecision(17)<<"[SF BLOCK ROW] "<<numbering_->row(row.row())<<" rhs="<<row.rightHandSide();for(const auto& a:row.coefficients())std::cout<<" "<<numbering_->row(a.column)<<":"<<a.value;std::cout<<'\n';}
        assembled_=true;
    }
    void solve(){
        solved_=false;if(!assembled_)throw std::runtime_error("Block solve before assembly");auto result=linear_->solve(matrix_);if(!result.converged || result.solution.size()!=matrix_.rows.size())throw std::runtime_error("Linear block solve failed/nonconverged; physical targets were not published");double squared=0,rhs=0;
        for(std::size_t n=0;n<matrix_.rows.size();++n){const auto& row=matrix_.rows[n];double r=-row.rightHandSide();for(const auto& a:row.coefficients())r+=a.value*result.solution.at(numbering_->row(a.column));squared+=r*r;rhs+=row.rightHandSide()*row.rightHandSide();}
        if(!std::isfinite(squared) || !std::isfinite(rhs))throw std::runtime_error("Nonfinite block residual norm; no publication");double residual=std::sqrt(squared),limit=std::max(contract_.linear.absoluteTolerance,contract_.linear.relativeTolerance*std::sqrt(rhs));if(!std::isfinite(residual) || residual>limit || std::any_of(result.solution.begin(),result.solution.end(),[](double v){return !std::isfinite(v);}))throw std::runtime_error("Linear block residual acceptance failed; physical targets were not published");
        std::cout<<std::setprecision(17)<<"[SF BLOCK] occurrence="<<call_.source.occurrence<<" rows="<<matrix_.rows.size()<<" iterations="<<result.iterations<<" residual="<<residual<<" limit="<<limit<<'\n';pending_=std::move(result.solution);solved_=true;
    }
    void publish(){if(!solved_)throw std::runtime_error("Block publication before accepted whole solution");for(std::size_t n=0;n<physical_.size();++n){port(n);if(physical_[n]->requireReadable)physical_[n]->requireReadable();}for(std::size_t n=0;n<physical_.size();++n){for(std::size_t cell=0;cell<cells_.size();++cell)physical_[n]->write(cells_[cell],0,pending_.at(n*cells_.size()+cell));state_.distributed.markModified(physical_[n]->name,physical_[n]->blockId);}assembled_=solved_=false;}
};

class AffineOperations {
    const System::CompiledEquationCall& call_;
    const System::CompiledTimeRecipe& recipe_;
    System::CompiledAffineEquation contract_;
    State::StateBundle& state_;
    System::StateRealization& realized_;
    Field& mesh_;
    State::DistributedFieldView* output_;
    State::DistributedFieldView materialized_;
    bool workspace_=false,physicalReady_=false;
    int publishedStep_=-1;
    std::vector<std::pair<std::string,unsigned long long>> inputRevisions_;
    std::shared_ptr<std::set<std::string>> publications_;
    std::weak_ptr<const FDM::ScalarEquationData> port_;
    std::vector<double> stage_,pending_,materializedValues_;
    std::vector<State::DistributedFieldView> views_;
    System::FormulaValues values_,physicalValues_;
    bool inputsBound_=false;
    int readyStage_=-1,readyStep_=-1;
    double readyTime_=0;
    template<class F>void owned(F f){int g=mesh_.NG();for(int k=g;k<g+mesh_.NZ();++k)for(int j=g;j<g+mesh_.NY();++j)for(int i=g;i<g+mesh_.NX();++i)f(mesh_.getIdx(i,j,k));}
    bool interior(int cell) const {
        int i,j,k;mesh_.getIJK(cell,i,j,k);const int g=mesh_.NG();
        return i>=g && i<g+mesh_.NX() && j>=g && j<g+mesh_.NY() && k>=g && k<g+mesh_.NZ();
    }
    void requireRegion(int cell) const {
        if(contract_.boundaryPreparation==System::WorkspaceBoundary::None && !interior(cell))
            throw std::runtime_error("Unprepared coefficient read region: step="+call_.source.occurrence+" equation="+call_.source.equation+" target="+contract_.target+" cell="+std::to_string(cell));
    }
    void calculate(const System::FormulaValues& values,std::vector<double>& output){
        if(!inputsBound_)throw std::runtime_error("Affine inputs have not completed startup binding");
        if(!workspace_ && port_.expired())throw std::runtime_error("Expired algebraic resource port");
        owned([&](int cell){double A=contract_.coefficient(values,cell,0),R=contract_.remainder(values,cell,0);if(!std::isfinite(A) || !std::isfinite(R) || A==0)throw std::runtime_error("LocalAffine singular/non-finite coefficient: equation="+call_.source.equation+" target="+contract_.target);double x=-R/A;if(!std::isfinite(x))throw std::runtime_error("Non-finite LocalAffine solution");output.at(cell)=x;});
        if(contract_.boundaryPreparation==System::WorkspaceBoundary::NearestInterior)extendNearestInterior(mesh_,output);
    }
    void write(){
        if(!workspace_ && output_->requireReadable)output_->requireReadable();
        if(contract_.boundaryPreparation==System::WorkspaceBoundary::NearestInterior)
            for(int cell=0;cell<mesh_.TotalSize();++cell)output_->write(cell,0,pending_.at(cell));
        else owned([&](int cell){output_->write(cell,0,pending_.at(cell));});
        physicalReady_=true;
        publishedStep_=state_.step;inputRevisions_.clear();
        for(const auto& input:contract_.valueReads)if(input.version==System::StateVersion::Current)
            inputRevisions_.push_back({input.symbol,realized_.revision(input.symbol)});
        if(workspace_)realized_.publishWorkspace(contract_.target);
        if(call_.target.kind==System::TargetKind::Physical && !workspace_)state_.distributed.markModified(output_->name,output_->blockId);
        else if(call_.target.kind!=System::TargetKind::Physical)publications_->insert(call_.target.workspace);
    }
public:
    AffineOperations(const System::CompiledEquationCall& call,const System::CompiledTimeRecipe& recipe,State::StateBundle& state,System::StateRealization& realized,const FDM::ScalarInstanceData* instance,std::shared_ptr<std::set<std::string>> publications)
        :call_(call),recipe_(recipe),contract_(std::any_cast<System::CompiledAffineEquation>(call.providerContract)),state_(state),realized_(realized),mesh_(state.singlePatch()),publications_(std::move(publications)){
        workspace_=realized.at(contract_.target).descriptor->storageBinding==System::StorageBinding::TransientWorkspace;
        if(contract_.boundaryPreparation!=realized.at(contract_.target).descriptor->workspaceBoundary
            || (contract_.boundaryPreparation!=System::WorkspaceBoundary::None && contract_.boundaryPreparation!=System::WorkspaceBoundary::NearestInterior))
            throw std::runtime_error("Affine output read-region contract differs from its STATE declaration: "+contract_.target);
        pending_.resize(mesh_.TotalSize());
        if(workspace_) {
            if(instance)throw std::runtime_error("Duplicate resource authority for affine workspace: "+contract_.target);
            materializedValues_.resize(mesh_.TotalSize());
            materialized_=State::workspaceView(realized.at(contract_.target).descriptor->storageKey,0,mesh_,materializedValues_,1);
            materialized_.read=[this](int cell,int){requireRegion(cell);return materializedValues_.at(cell);};
            materialized_.requireReadable=[this]{
                if(!physicalReady_ || publishedStep_!=state_.step)throw std::runtime_error("Unpublished materialized STATE in this step: "+contract_.target+"; execute its HOW relation first");
                for(const auto& input:inputRevisions_) {
                    const auto& backing=realized_.at(input.first);
                    if(backing.fields.size()!=1)throw std::runtime_error("Missing materialized dependency backing: "+input.first);
                    if(backing.fields.front()->requireReadable)backing.fields.front()->requireReadable();
                    if(realized_.revision(input.first)!=input.second)throw std::runtime_error("Stale materialized dependency: "+contract_.target+" reads "+input.first+"; HOW must recompute the relation");
                }
            };
            realized.bindView(contract_.target,System::StateViewKind::Physical,materialized_);
        } else {
            if(!instance)throw std::runtime_error("Missing affine physical resource port: "+contract_.target);
            port_=instance->data;
        }
        output_=call.target.kind==System::TargetKind::Physical?realized.at(contract_.target).fields.at(0)
            :realized.view(contract_.target,call.target.kind==System::TargetKind::Working?System::StateViewKind::Working:System::StateViewKind::Correction).fields.at(0);
        const auto port=port_.lock();if((!workspace_ && (!port || port->target!=contract_.target)) || output_->geometry!=&mesh_ || output_->components!=1 || !output_->write)throw std::runtime_error("Invalid affine STATE resource/target binding");
        if(call_.stageLocal){stage_.resize(mesh_.TotalSize());views_.reserve(recipe_.stageCount());for(int stage=0;stage<recipe_.stageCount();++stage){views_.push_back(State::workspaceView(call_.stageWorkspace,0,mesh_,stage_,1));views_.back().read=[this,stage](int cell,int){double time=state_.time+recipe_.stage(stage).abscissa*state_.dt;realized_.requireStageRead(stage,time);if(readyStage_!=stage || readyStep_!=state_.step || readyTime_!=time)throw std::runtime_error("Unpublished algebraic Stage: "+contract_.target);requireRegion(cell);return stage_.at(cell);};realized_.bindView(contract_.target,System::StateViewKind::Stage,views_.back(),stage);}}
    }
    void bindReads(){
        if(inputsBound_)throw std::runtime_error("Affine input readers were already bound");
        for(const auto& read:contract_.valueReads){
            auto current=read;current.kind=System::TargetKind::Physical;current.version=System::StateVersion::Current;
            auto physical=realized_.reader(current,recipe_.stageCount());
            physicalValues_.boundReads.push_back(read.version==System::StateVersion::OldTime?realized_.oldTimeReader(read.symbol):physical);
            if(read.kind!=System::TargetKind::Physical) {
                const auto kind=read.kind==System::TargetKind::Working?System::StateViewKind::Working:System::StateViewKind::Correction;
                const auto& resource=realized_.view(read.symbol,kind);
                if(resource.fields.size()!=1)throw std::runtime_error("Missing qualified affine STATE view");
                auto work=*resource.fields.front();
                values_.boundReads.push_back([work,key=resource.descriptor.storage,ready=publications_](int cell,int c){
                    if(!ready->count(key))throw std::runtime_error("Unpublished Working/Correction STATE: "+key);
                    if(work.requireReadable)work.requireReadable();return work.read(cell,c);
                });
            } else values_.boundReads.push_back(read.version==System::StateVersion::Stage?realized_.stageReader(read.symbol,recipe_.stageCount()):read.version==System::StateVersion::OldTime?realized_.oldTimeReader(read.symbol):physical);
        }
        inputsBound_=true;
    }
    void solve(const Run::ExecutionContext& context){if(call_.stageLocal){if(context.stageCount!=recipe_.stageCount())throw std::runtime_error("Affine Stage recipe/context mismatch");double time=state_.time+recipe_.stage(context.stageIndex).abscissa*state_.dt;realized_.requireStageRead(context.stageIndex,time);readyStage_=-1;calculate(values_,stage_);readyTime_=time;readyStep_=state_.step;readyStage_=context.stageIndex;}else{calculate(values_,pending_);write();}}
    /// Final relation evaluation is the selected LocalAffine lifecycle, not a COPY of its last stage.
    void publish(){if(readyStep_!=state_.step || readyStage_!=recipe_.stageCount()-1)throw std::runtime_error("Algebraic physical publication before complete stages");calculate(physicalValues_,pending_);write();readyStage_=-1;}
    std::vector<double> samples(System::ResidualKind kind,const Run::ExecutionContext& context){
        if(!workspace_ && port_.expired())throw std::runtime_error("Expired residual STATE resource");
        if(call_.stageLocal && (context.stageIndex<0 || context.stageCount!=recipe_.stageCount()))throw std::runtime_error("Residual Stage context mismatch");
        if(kind!=System::ResidualKind::StateChange && kind!=System::ResidualKind::DiscreteEquation)throw std::runtime_error("LocalAffine has no linear-solver residual sampler");
        bool stageReady=false;
        if(call_.stageLocal){
            const double time=state_.time+recipe_.stage(context.stageIndex).abscissa*state_.dt;
            realized_.requireStageRead(context.stageIndex,time);
            stageReady=readyStep_==state_.step && readyStage_==context.stageIndex && readyTime_==time;
            if(!stageReady && kind!=System::ResidualKind::StateChange)throw std::runtime_error("Residual before algebraic Stage publication");
        }
        if(!stageReady && output_->requireReadable)output_->requireReadable();
        std::vector<double> result;
        owned([&](int cell){
            // An unpublished stage has the physical initial guess, never a prior stage's workspace.
            double x=stageReady?stage_.at(cell):output_->read(cell,0);
            double value=kind==System::ResidualKind::StateChange?x:contract_.coefficient(values_,cell,0)*x+contract_.remainder(values_,cell,0);
            if(!std::isfinite(value))throw std::runtime_error("Nonfinite LocalAffine residual sample");
            result.push_back(value);
        });return result;
    }
};
}
std::vector<Time::TemporalCallbacks> scalarParticipants(const System::CompiledNumericalSystem& numerical,
        const System::CompiledSolvePlan& plan,State::StateBundle& state,const FDM::SolverServices& services,
        System::StateRealization& realized,const double& limit,
        std::function<StageMassFluxView(StageFluxIdentity)> massFlux) {
    if (state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
        throw std::runtime_error("Unsupported scalar execution: serial single patch required.");
    System::validateCoefficientContracts(plan.compiledProgram);
    std::vector<const System::CompiledEquationCall*> calls;
    for (const auto& call:plan.compiledProgram.steps) if(std::any_cast<System::CompiledScalarEquation>(&call.providerContract))calls.push_back(&call);
    if (services.scalarEquation && (!services.scalarInstances.empty() || calls.size()!=1))
        throw std::runtime_error("Single scalar data port cannot bind multiple occurrences or coexist with instance table.");
    validateScalarInstances(plan,services);
    for(const auto* consumer:calls) {
        const auto& contract=std::any_cast<const System::CompiledScalarEquation&>(consumer->providerContract);
        if(contract.diffusivityState.empty())continue;
        const auto& input=*realized.at(contract.diffusivityState).descriptor;
        if(input.storageBinding!=System::StorageBinding::TransientWorkspace)continue;
        const auto count=std::count_if(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[&](const auto& call) {
            const auto* affine=std::any_cast<System::CompiledAffineEquation>(&call.providerContract);
            return affine && affine->target==input.id && call.stageLocal
                && affine->boundaryPreparation==System::WorkspaceBoundary::NearestInterior;
        });
        if(count!=1 || input.workspaceBoundary!=System::WorkspaceBoundary::NearestInterior)
            throw std::runtime_error("Missing unique prepared spatial coefficient publisher: step="+consumer->source.occurrence+" equation="+consumer->source.equation+" target="+consumer->target.symbol+" input="+input.id);
    }
    std::vector<Time::TemporalCallbacks> callbacks;
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
        const auto found=std::find_if(services.scalarInstances.begin(),services.scalarInstances.end(),[&](const auto& item){return item.occurrence==fieldPort(*call);});
        const auto locked=found==services.scalarInstances.end()?std::shared_ptr<const FDM::ScalarEquationData>{}:found->data.lock();
        const auto* data=services.scalarEquation?services.scalarEquation:locked.get();
        if (!data) throw std::runtime_error("Missing scalar instance data: "+call->source.occurrence);
        auto owner=checked("STATE binding",-1,[&]{return std::make_shared<ScalarOperations>(numerical,*call,state,*data,realized,limit,massFlux,
            std::any_of(call->physicalResources.begin(),call->physicalResources.end(),[](const auto& need){return need.fields && need.fields->boundaryRequired;}));});
        if (locked) owner->trackData(found->data);
        Time::TemporalCallbacks cb;cb.identity=call->source.occurrence;
        cb.bindReads=[owner]{owner->bindInputs();};
        cb.stepSize=[owner,checked]{return checked("dt",-1,[&]{return owner->stepSize();});};
        cb.snapshot=[owner,checked]{checked("snapshot",-1,[&]{owner->begin();});};
        cb.prepare=[owner,checked](const auto& c){checked("prepare",c.stageIndex,[&]{owner->prepareStage(c);});};
        cb.rhs=[owner,checked](const auto& c){checked("rhs",c.stageIndex,[&]{owner->evaluateRHS(c);});};
        cb.advance=[owner,checked](const auto& c){checked("advance",c.stageIndex,[&]{owner->advanceStage(c);});};
        cb.validatePublish=[owner,checked]{checked("validate publish",-1,[&]{owner->validatePublish();});};
        cb.publish=[owner,checked]{checked("publish",-1,[&]{owner->publish();});};
        callbacks.push_back(std::move(cb));
    }
    return callbacks;
}
std::vector<std::function<void()>> bindAlgebraicOperations(Run::OpRegistry& ops,const System::CompiledNumericalSystem& numerical,const System::CompiledSolvePlan& plan,State::StateBundle& state,const FDM::SolverServices& services,System::StateRealization& realized){
    validateScalarInstances(plan,services);
    std::vector<std::function<void()>> inputs;
    auto publications=std::make_shared<std::set<std::string>>();
    for(const auto& call:plan.compiledProgram.steps)if(const auto* contract=std::any_cast<System::CompiledStageClosure>(&call.providerContract)) {
        if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()) || services.immersed.boundary || services.immersed.constraint || services.immersed.system)
            throw std::runtime_error("Unsupported closure stage treatment: step="+call.source.occurrence+" requires serial single patch without IBM");
        if(!services.transportModel)throw std::runtime_error("Missing SST closure implementation resource: "+call.source.occurrence);
        auto owner=std::make_shared<SstClosureOperations>(state,state.singlePatch(),*services.transportModel,realized,numerical.time.recipe.stageCount(),contract->targets,call);
        ops.bind(call.backendOperation,call.backendProvider,[owner,&state,&numerical](const Run::ExecutionContext& c){owner->refresh(false,c.stageIndex,state.time+numerical.time.recipe.stage(c.stageIndex).abscissa*state.dt);});
        ops.bind(call.temporalPublish,call.backendProvider,[owner]{owner->refresh(true);});
        inputs.push_back([owner]{owner->bindInputs();});
    }

    for(const auto& call:plan.compiledProgram.steps)if(std::any_cast<System::CompiledAffineEquation>(&call.providerContract)){
        if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
            throw std::runtime_error("Unsupported LocalAffine layout: one serial patch required: "+call.source.occurrence);
        const FDM::ScalarInstanceData* instance=nullptr;
        if(!call.physicalResources.empty()) {
            const auto found=std::find_if(services.scalarInstances.begin(),services.scalarInstances.end(),[&](const auto& p){return p.occurrence==fieldPort(call);});
            if(found==services.scalarInstances.end())throw std::runtime_error("Missing affine resource instance: "+call.source.occurrence);
            instance=&*found;
        }
        auto owner=std::make_shared<AffineOperations>(call,numerical.time.recipe,state,realized,instance,publications);
        const auto checked=[label="LocalAffine step="+call.source.occurrence+" equation="+call.source.equation+" target="+call.target.symbol+": "](auto&& action) {
            try{action();}catch(const std::exception& e){throw std::runtime_error(label+e.what());}
        };
        inputs.push_back([owner,checked]{checked([&]{owner->bindReads();});});
        ops.bind(call.backendOperation,call.backendProvider,[owner,checked](const Run::ExecutionContext& c){checked([&]{owner->solve(c);});});
        if(call.stageLocal)ops.bind(call.temporalPublish,call.backendProvider,[owner,checked]{checked([&]{owner->publish();});});
        for(const auto& residual:plan.compiledProgram.residuals)if(residual.occurrence==call.source.occurrence)
            ops.bindResidual(residual.id,residual.provider,[owner,kind=residual.kind](const Run::ExecutionContext& c){return owner->samples(kind,c);});
    }
    for(const auto& call:plan.compiledProgram.steps)if(std::any_cast<System::CompiledLinearBlock>(&call.providerContract)){
        auto owner=std::make_shared<BlockOperations>(call,state,realized,services);inputs.push_back([owner]{owner->bindReads();});
        std::string label="Block occurrence="+call.source.occurrence+" provider="+call.equationMethod+" equations/targets=[";for(const auto& c:call.calls)label+=" "+c.equation+"->"+c.target;label+=" ]: ";
        const auto checked=[owner,label](auto operation){return [owner,label,operation]{try{((*owner).*operation)();}catch(const std::exception& e){throw std::runtime_error(label+e.what());}};};
        ops.bind(call.operations.at(0).operation,call.backendProvider,checked(&BlockOperations::assemble));ops.bind(call.operations.at(1).operation,call.backendProvider,checked(&BlockOperations::solve));ops.bind(call.operations.at(2).operation,call.backendProvider,checked(&BlockOperations::publish));
    }
    const auto bindCommit=[&](const auto& self,const System::ExecutionScope& scope)->void {
        if(!scope.commitSources.empty()) {
            struct Copy {State::DistributedFieldView source,destination;std::string key;int components;};
            std::vector<Copy> copies;
            for(std::size_t i=0;i<scope.commitSources.size();++i) {
                const auto& src=scope.commitSources[i];const auto& dst=scope.commitTargets.at(i);
                const auto kind=src.kind==System::TargetKind::Working?System::StateViewKind::Working:System::StateViewKind::Correction;
                const auto& view=realized.view(src.symbol,kind);const auto& physical=realized.at(dst.symbol);
                if(view.fields.size()!=1 || physical.fields.size()!=1 || !physical.fields[0]->write
                    || view.fields[0]->geometry!=physical.fields[0]->geometry || view.fields[0]->components!=physical.descriptor->components)
                    throw std::runtime_error("Variable Commit resource/layout missing: "+scope.id);
                copies.push_back({*view.fields[0],*physical.fields[0],view.descriptor.storage,physical.descriptor->components});
            }
            ops.bind(System::instanceOperation("state.commit",scope.id),"execution.state",[copies,publications,&state]{
                auto& mesh=state.singlePatch();const int g=mesh.NG();std::vector<std::vector<double>> pending(copies.size());
                for(std::size_t n=0;n<copies.size();++n) {
                    const auto& copy=copies[n];
                    if(!publications->count(copy.key))throw std::runtime_error("Variable Commit before accepted Working publication: "+copy.key);
                    if(copy.source.requireReadable)copy.source.requireReadable();if(copy.destination.requireReadable)copy.destination.requireReadable();
                    for(int k=g;k<g+mesh.NZ();++k)for(int j=g;j<g+mesh.NY();++j)for(int i=g;i<g+mesh.NX();++i)for(int c=0;c<copy.components;++c){double value=copy.source.read(mesh.getIdx(i,j,k),c);if(!std::isfinite(value))throw std::runtime_error("Nonfinite variable Commit; physical targets unmodified");pending[n].push_back(value);}
                }
                for(std::size_t n=0;n<copies.size();++n){const auto& copy=copies[n];std::size_t at=0;for(int k=g;k<g+mesh.NZ();++k)for(int j=g;j<g+mesh.NY();++j)for(int i=g;i<g+mesh.NX();++i)for(int c=0;c<copy.components;++c)copy.destination.write(mesh.getIdx(i,j,k),c,pending[n].at(at++));state.distributed.markModified(copy.destination.name,copy.destination.blockId);}
            });
        }
        for(const auto& child:scope.children)self(self,child);
    };bindCommit(bindCommit,plan.compiledProgram.root);
    return inputs;
}
}
