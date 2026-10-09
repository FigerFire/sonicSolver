#include "SF_scalarBinding.h"
#include "solver/system/SF_scalarMethod.h"
#include "core/system/SF_operationIds.h"
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
    std::array<double,3> spacing_{};
    std::vector<double> old_,stage_,rhs_;
    std::array<std::vector<double>,4> increments_;
    State::DistributedFieldView oldView_,stageView_;
    System::StateRealization& realized_;
    std::weak_ptr<const FDM::ScalarEquationData> lifetime_;
    bool checkedLifetime_=false;
    void requireData() const {
        if (checkedLifetime_ && lifetime_.expired()) throw std::runtime_error("Scalar instance data lifetime expired.");
    }
    System::FormulaValues values_;
    double stageTime_=0;
    int nextStage_=-1;

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
        for (int axis=0;axis<3;++axis) if (shape_[axis]>1)
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
            for (int axis=0;axis<3;++axis) if (shape_[axis]>1)
                for (int side=0;side<2;++side) {
                    const auto& bc=data_.boundary[2*axis+side];
                    if (at[axis]!=g+(side?shape_[axis]-1:0) || bc.kind!=FDM::ScalarBoundaryKind::FixedValue) continue;
                    const double next=known(bc.value,cell,time);
                    if (assigned && std::abs(value-next)>1e-12*std::max({1.0,std::abs(value),std::abs(next)}))
                        throw std::runtime_error("Conflicting scalar fixedValue corner boundaries.");
                    value=next;assigned=true;
                }
            if (assigned) stage_[cell]=value;
        });
        // Central2 needs exactly one exterior node. Reflection about the wall
        // gives zero normal derivative; fixedValue uses odd reflection.
        for (int axis=0;axis<3;++axis) if (shape_[axis]>1)
            owned([&](int cell) {
                auto at=index(cell);
                for (int side=0;side<2;++side) {
                    const int wall=g+(side?shape_[axis]-1:0);
                    if (at[axis]!=wall) continue;
                    auto inside=at,outside=at;inside[axis]+=side?-1:1;outside[axis]+=side?1:-1;
                    const double interior=stage_[geometry_.getIdx(inside[0],inside[1],inside[2])];
                    const auto& bc=data_.boundary[2*axis+side];
                    stage_[geometry_.getIdx(outside[0],outside[1],outside[2])]=
                        bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?interior:2*stage_[cell]-interior;
                }
            });
    }
public:
    void trackData(std::weak_ptr<const FDM::ScalarEquationData> lifetime) {
        lifetime_=std::move(lifetime);checkedLifetime_=true;
    }
    ScalarOperations(const System::CompiledNumericalSystem& numerical,const System::CompiledEquationCall& call,
        State::StateBundle& state,const FDM::ScalarEquationData& data,System::StateRealization& realized,const double& limit)
        :state_(state),geometry_(state.singlePatch()),data_(data),recipe_(numerical.time.recipe),limit_(limit),
         fixedDt_(numerical.dt.maxDeltaT),contract_(std::any_cast<System::CompiledScalarEquation>(call.providerContract)),
         physical_(realized.at(call.target.symbol).fields.at(0)),shape_{geometry_.NX(),geometry_.NY(),geometry_.NZ()},realized_(realized) {
        for (const auto& name:contract_.stageReads)
            if(std::none_of(call.stateUses.begin(),call.stateUses.end(),[&](const auto& use){return use.symbol==name && use.version==System::StateVersion::Stage;}))
                throw std::runtime_error("Scalar compiled cross read lacks declared Stage binding: "+name);
        if(physical_->geometry!=&geometry_ || !physical_->read || !physical_->write)
            throw std::runtime_error("Scalar storage geometry/read/write binding mismatch.");
        if (data.target!=contract_.target) throw std::runtime_error("Scalar boundary/source port targets a different STATE.");
        if (geometry_.NG()<1) throw std::runtime_error("Scalar Central2 requires one ghost layer.");
        if (physical_->components!=1) throw std::runtime_error("Scalar physical storage requires exactly one component.");
        const std::array<int,3> origin{geometry_.NG(),geometry_.NG(),geometry_.NG()};
        for (int axis=0;axis<3;++axis) if (shape_[axis]>1) {
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
            for (int axis=0;axis<3;++axis) {
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
            stage_[cell]=physical_->read(cell,0);
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
        owned([&](int cell) {
            rhs_[cell]=fixed(cell)?0:contract_.rhs(values_,cell,0)
                -(contract_.conservativeTransport?Discretization::conservativeScalarAdvection(values_,cell):0);
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
    }
    void publish() {
        for (int cell=0;cell<geometry_.TotalSize();++cell) {
            physical_->write(cell,0,stage_[cell]);
        }
        state_.distributed.markModified(physical_->name);
        nextStage_=-1;
    }
};
}
std::vector<Time::TemporalCallbacks> scalarParticipants(const System::CompiledNumericalSystem& numerical,
        const System::CompiledSolvePlan& plan,State::StateBundle& state,const FDM::SolverServices& services,
        System::StateRealization& realized,const double& limit) {
    if (state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
        throw std::runtime_error("Unsupported scalar execution: serial single patch required.");
    std::vector<const System::CompiledEquationCall*> calls;
    for (const auto& call:plan.compiledProgram.steps) if (call.backendProvider==System::ScalarOps::Provider || call.backendProvider==System::ScalarOps::TransportProvider) calls.push_back(&call);
    if (services.scalarEquation && (!services.scalarInstances.empty() || calls.size()!=1))
        throw std::runtime_error("Single scalar data port cannot bind multiple occurrences or coexist with instance table.");
    for (std::size_t i=0;i<services.scalarInstances.size();++i) {
        const auto& item=services.scalarInstances[i];
        if (item.data.expired() || std::count_if(services.scalarInstances.begin(),services.scalarInstances.end(),
                [&](const auto& entry){return entry.occurrence==item.occurrence;})!=1
            || std::none_of(calls.begin(),calls.end(),[&](const auto* call){return call->source.occurrence==item.occurrence;}))
            throw std::runtime_error("Invalid/duplicate/unconsumed scalar instance data key.");
    }
    std::vector<Time::TemporalCallbacks> callbacks;
    std::vector<std::shared_ptr<ScalarOperations>> owners;
    for (const auto* call:calls) {
        const auto found=std::find_if(services.scalarInstances.begin(),services.scalarInstances.end(),
            [&](const auto& item){return item.occurrence==call->source.occurrence;});
        const auto locked=found==services.scalarInstances.end()?std::shared_ptr<const FDM::ScalarEquationData>{}:found->data.lock();
        const auto* data=services.scalarEquation?services.scalarEquation:locked.get();
        if (!data) throw std::runtime_error("Missing scalar instance data: "+call->source.occurrence);
        auto owner=std::make_shared<ScalarOperations>(numerical,*call,state,*data,realized,limit);
        if (locked) owner->trackData(found->data);
        owners.push_back(owner);
        Time::TemporalCallbacks cb;cb.identity=call->source.occurrence;
        cb.stepSize=[owner]{return owner->stepSize();};cb.snapshot=[owner]{owner->begin();};
        cb.prepare=[owner](const auto& c){owner->prepareStage(c);};
        cb.rhs=[owner](const auto& c){owner->evaluateRHS(c);};
        cb.advance=[owner](const auto& c){owner->advanceStage(c);};
        cb.validatePublish=[owner]{owner->validatePublish();};cb.publish=[owner]{owner->publish();};
        callbacks.push_back(std::move(cb));
    }
    for (const auto& owner:owners) owner->bindInputs();
    return callbacks;
}
}
