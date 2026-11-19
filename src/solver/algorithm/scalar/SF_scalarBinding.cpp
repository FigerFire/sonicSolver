#include "solver/system/SF_sstMethod.h"
#include "SF_scalarBinding.h"
#include "methods/numerics/scalar/SF_scalarTransport.h"
#include "solver/algorithm/SF_iteration.h"
#include "solver/system/SF_affine.h"
#include "solver/system/SF_block.h"
#include "solver/system/SF_scalarMethod.h"
#include "core/system/SF_operationIds.h"
#include "core/mesh/SF_dimension.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "solver/discretization/convection/SF_scalarUpwind.h"
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
#include <fstream>
#include <filesystem>
namespace SF::SolverAlgorithm {
namespace {
/// Owns only numerical workspace. Physical writes use the realized STATE view.
class ScalarOperations {
    State::StateBundle& state_;
    Field& geometry_;
    const FDM::ScalarEquationData& data_;
    const System::CompiledTimeRecipe& recipe_;
    const double& limit_;
    double fixedDt_;int fixedStep_=-1,fixedUpdated_=-1;EquationInterval fixedInterval_;
    const System::CompiledScalarEquation& contract_;
    State::DistributedFieldView* physical_;
    std::array<int,3> shape_;
    std::array<bool,3> active_{};
    std::array<double,3> spacing_{};
    std::vector<double> old_,stage_,rhs_;
    std::map<std::string,unsigned long long> frozenInputs_;
    // Metric FE uses its boundary/bounds-prepared base; published OldTime stays
    // the interval-entry snapshot for every equation sharing this time group.
    std::vector<double> updateBase_;
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
    std::shared_ptr<std::ofstream> faceTrace_;
    unsigned long scratchRevision_=0,boundaryRevision_=~0UL;
    double boundaryTime_=0;
    std::function<StageMassFluxView(StageFluxIdentity)> massFlux_;
public:
    void bindFluxes(const FluxReaders& fluxes) {
        if(contract_.sharedMassFlux && !fluxes.mass)
            throw std::runtime_error("Shared scalar requires a bound canonical NS mass flux producer.");
        if(contract_.sharedMassFlux)massFlux_=fluxes.mass;
        bindVolumeFlux(fluxes.volume);
    }
private:
    void bindVolumeFlux(std::function<double(std::string_view,int,int)> reader) {
        if(!contract_.volumeFlux.empty() && !reader)throw std::runtime_error("Missing published volume flux binding");
        if(!reader)return;
        values_.volumeFlux=[this,reader=std::move(reader)](int axis,int lower) {
            const double flux=reader(contract_.volumeFlux,axis,lower);
            if(faceTrace_ && *faceTrace_) {
                const int upper=values_.spatialNeighbor(lower,axis,1),upwind=flux>=0?lower:upper;
                *faceTrace_<<axis<<','<<lower<<','<<upper<<','<<flux<<','<<values_.boundReads.at(0)(upwind,0)<<'\n';
            }
            return flux;
        };
    }
private:

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
        if(contract_.metric)return false;
        const auto at=index(cell);const int g=geometry_.NG();
        for (int axis=0;axis<3;++axis) if (active_[axis])
            for (int side=0;side<2;++side)
                if (at[axis]==g+(side?shape_[axis]-1:0)
                    && data_.boundary[2*axis+side].kind==FDM::ScalarBoundaryKind::FixedValue) return true;
        return false;
    }
    void metricBoundary(bool after) {
        ScalarArrayView scratch{geometry_,stage_,contract_.target};const auto& options=*contract_.metric;
        FDM::Scalar::applyBoundaryConditions(geometry_,scratch,data_.setBoundary,options.ilw,options.order);
        const auto report=FDM::Scalar::applyBounds(geometry_,scratch,options.bounds);
        if(after && report.violated() && options.bounds.mode==FDM::Scalar::BoundMode::Project)
            FDM::Scalar::applyBoundaryConditions(geometry_,scratch,data_.setBoundary,options.ilw,options.order);
        boundaryRevision_=scratchRevision_;boundaryTime_=stageTime_;
    }
    template<class Read,class Write> void nodalBoundary(double time,const Read& read,const Write& write) {
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
            if (assigned) write(cell,component,value*(contract_.primitiveStorage?geometry_(at[0],at[1],at[2],0):1));
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
                    const double interior=read(geometry_.getIdx(inside[0],inside[1],inside[2]),component);
                    const auto& bc=data_.boundary[2*axis+side];
                    const int ghost=geometry_.getIdx(outside[0],outside[1],outside[2]);
                    if(contract_.sharedMassFlux) {
                        const double rho=geometry_(at[0],at[1],at[2],0),ri=geometry_(inside[0],inside[1],inside[2],0),rg=geometry_(outside[0],outside[1],outside[2],0);
                        if(!std::isfinite(rho) || !std::isfinite(ri) || !std::isfinite(rg) || rho<=0 || ri<=0 || rg<=0)
                            throw std::runtime_error("Primitive scalar boundary density must be positive and finite.");
                        const double wallConcentration=bc.kind==FDM::ScalarBoundaryKind::FixedValue?boundaryValue(bc,cell,time,component)/(contract_.primitiveStorage?1:rho):0;
                        write(ghost,component,(bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?read(cell,component)/rho:2*wallConcentration-interior/ri)*rg);
                    } else write(ghost,component,bc.kind==FDM::ScalarBoundaryKind::ZeroGradient?interior:2*read(cell,component)-interior);
                }
            });
    }
    void boundary(double time) {
        if(contract_.metric){metricBoundary(false);return;}
        if(boundaryRevision_==scratchRevision_ && boundaryTime_==time)return;
        nodalBoundary(time,[this](int cell,int c){return stage_[entry(cell,c)];},
            [this](int cell,int c,double value){stage_[entry(cell,c)]=value;});
        boundaryRevision_=scratchRevision_;boundaryTime_=time;
    }
public:
    void trackData(std::weak_ptr<const FDM::ScalarEquationData> lifetime) {
        lifetime_=std::move(lifetime);checkedLifetime_=true;
    }
    ScalarOperations(const System::CompiledNumericalSystem& numerical,const System::CompiledEquationCall& call,
        State::StateBundle& state,const FDM::ScalarEquationData& data,System::StateRealization& realized,const double& limit,bool stencil)
        :state_(state),geometry_(state.singlePatch()),data_(data),recipe_(numerical.time.recipe),limit_(limit),
         fixedDt_(numerical.dt.maxDeltaT),contract_(std::any_cast<const System::CompiledScalarEquation&>(call.providerContract)),
         physical_(realized.at(call.target.symbol).fields.at(0)),shape_{geometry_.NX(),geometry_.NY(),geometry_.NZ()},realized_(realized) {
        for(int axis=0;axis<3;++axis)active_[axis]=stencil && shape_[axis]>1 && data_.activeDimensions[axis] && Math::isDirectionActiveIndex(axis);
        for (const auto& read:contract_.valueReads)
            if(std::none_of(call.stateUses.begin(),call.stateUses.end(),[&](const auto& use){return use.symbol==read.symbol && use.version==read.version;}))
                throw std::runtime_error("Scalar compiled input lacks declared version binding: "+read.symbol);
        if(physical_->geometry!=&geometry_ || !physical_->read || !physical_->write)
            throw std::runtime_error("Scalar storage geometry/read/write binding mismatch.");
        if (data.target!=contract_.target) throw std::runtime_error("Scalar boundary/source port targets a different STATE.");
        if (stencil && geometry_.NG()<1) throw std::runtime_error("Scalar Central2 requires one ghost layer.");
        if (physical_->components!=contract_.components) throw std::runtime_error("Physical storage component count differs from compiled target.");
        const std::array<int,3> origin{geometry_.NG(),geometry_.NG(),geometry_.NG()};
        if(contract_.metric && contract_.metric->ilw && Math::activeDimensionCount()<2)throw std::runtime_error("Metric scalar ILW requires an active tangential direction");
        for (int axis=0;axis<3;++axis) if (active_[axis]) {
            if(contract_.metric){spacing_[axis]=1;continue;}
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
        if(!contract_.metric)owned([&](int cell) {
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
        const auto total=geometry_.TotalSize()*contract_.components;old_.resize(total);stage_.resize(total);rhs_.resize(total);if(contract_.metric)updateBase_.resize(total);
        for (int i=0;i<recipe_.stageCount();++i) increments_[i].resize(total);
        if(std::getenv("SF_SCALAR_TRACE")) {
            std::ostringstream trace;trace<<"[SF SCALAR storage] target="<<contract_.target<<" old="<<static_cast<const void*>(old_.data())<<" stage="<<static_cast<const void*>(stage_.data())<<" rhs="<<static_cast<const void*>(rhs_.data());
            std::cout<<trace.str()<<'\n';
        }
        oldView_=State::workspaceView(call.oldTimeWorkspace,0,geometry_,old_,contract_.components);
        stageView_=State::workspaceView(call.stageWorkspace,0,geometry_,stage_,contract_.components);
        stageView_.prepareSpatial=[this](double time){realized_.requireStageRead(nextStage_,stageTime_);if(time!=stageTime_)throw std::runtime_error("Scalar Stage spatial preparation time mismatch");boundary(time);};
        rhsView_=State::workspaceView(call.residualWorkspace,0,geometry_,rhs_,contract_.components);
        if(!contract_.sequential)realized.bindView(contract_.target,System::StateViewKind::Workspace,rhsView_);
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
        if ((!contract_.sequential || contract_.fixedTime) && realized.requestsView(contract_.target,System::StateViewKind::OldTime))
            realized.bindView(contract_.target,System::StateViewKind::OldTime,oldView_);
        if(!contract_.sequential)for (int i=0;i<recipe_.stageCount();++i)
            realized.bindView(contract_.target,System::StateViewKind::Stage,stageView_,i);
        values_.metricFaceVelocity=[this](std::size_t slot,int cell,int axis) {
            const auto at=index(cell);return FDM::Scalar::metricFaceVelocity(geometry_,at[0],at[1],at[2],static_cast<Math::Dir>(axis),
                [&](int i,int j,int k,int c){return values_.boundReads.at(slot)(geometry_.getIdx(i,j,k),c);});
        };
        values_.metricDiffusion=[this](std::size_t slot,int cell,double d) {
            const auto at=index(cell);return FDM::Scalar::metricDiffusion(geometry_,[&](int i,int j,int k){return values_.boundReads.at(slot)(geometry_.getIdx(i,j,k),0);},contract_.target,at[0],at[1],at[2],d);
        };
        values_.spacing=[this](int axis) {return spacing_.at(axis);};
        values_.spatialNeighbor=[this](int cell,int axis,int sign){auto at=index(cell);at[axis]+=sign;return geometry_.getIdx(at[0],at[1],at[2]);};
        values_.inverseVolume=[this](int cell){auto at=index(cell);return geometry_.Jac(at[0],at[1],at[2]);};
        if(contract_.sequential && !contract_.metric && !contract_.conservativeTransport && !physical_->prepareSpatial)
            physical_->prepareSpatial=[this](double time) {
                requireData();
                nodalBoundary(time,physical_->read,physical_->write);
                state_.distributed.markModified(physical_->name,physical_->blockId);
            };
        if(!contract_.volumeFlux.empty())for(int axis=0;axis<3;++axis)if(active_[axis])for(int side=0;side<2;++side)
            if(data_.boundary[2*axis+side].kind!=FDM::ScalarBoundaryKind::FixedValue)
                throw std::runtime_error("Shared volume transport currently requires fixedValue scalar boundary donor nodes on each active face");

    }
    void bindInputs() {
        if(inputsBound_)throw std::runtime_error("Scalar input readers were already bound");
        for (const auto& input:contract_.valueReads) {
            if(contract_.sequential && input.symbol==contract_.target && input.version==System::StateVersion::Current) {
                values_.boundReads.push_back([this](int cell,int c){return stage_.at(entry(cell,c));});continue;
            }
            if(input.version==System::StateVersion::Stage)values_.boundReads.push_back(realized_.stageReader(input.symbol,recipe_.stageCount()));
            else if(input.version==System::StateVersion::Frozen) {
                auto reader=realized_.reader(input,recipe_.stageCount());
                values_.boundReads.push_back([this,input,reader=std::move(reader)](int cell,int c) {
                    const auto at=frozenInputs_.find(input.symbol);
                    if(at==frozenInputs_.end() || at->second!=realized_.revision(input))
                        throw std::runtime_error("Frozen publication overwritten/not captured: target="+contract_.target+" input="+input.symbol);
                    return reader(cell,c);
                });
            }
            else if(input.version==System::StateVersion::OldTime)values_.boundReads.push_back(realized_.oldTimeReader(input.symbol));
            else if(input.version==System::StateVersion::Current) {
                equationField(realized_,input,geometry_);
                values_.boundReads.push_back(realized_.reader(input,recipe_.stageCount()));
            } else throw std::runtime_error("Unsupported runtime coefficient version: "+input.symbol);
        }
        values_.boundNeighbor=[this](std::size_t slot,int cell,int axis,int offset,int component) {
            auto at=index(cell);at.at(axis)+=offset;
            return values_.boundReads.at(slot)(geometry_.getIdx(at[0],at[1],at[2]),component);
        };
        values_.coefficientNeighbor=[this](std::size_t slot,int cell,int axis,int offset) {
            auto at=index(cell);at.at(axis)+=offset;const int lo=geometry_.NG(),hi=lo+shape_.at(axis)-1;
            // Match the existing nodal scalar ghost reflection; the coefficient owner is read, never extended or overwritten.
            const auto& read=contract_.valueReads.at(slot);
            const auto& symbol=*realized_.at(read.symbol).descriptor;
            if(symbol.storageBinding!=System::StorageBinding::TransientWorkspace) {
                if(at[axis]<lo)at[axis]=2*lo-at[axis];else if(at[axis]>hi)at[axis]=2*hi-at[axis];
            }
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
        ++scratchRevision_;if(!contract_.metric)boundary(state_.time);old_=stage_;nextStage_=0;
    }
    void prepareStage(const Run::ExecutionContext& context) {
        requireData();
        recipe_.requireProviderStages(context.stageCount);
        const int n=context.stageIndex;
        if (n!=nextStage_ || n<0 || n>=recipe_.stageCount()) throw std::runtime_error("Scalar stage order does not match frozen recipe.");
        stageTime_=contract_.fixedTime?state_.time+state_.dt:state_.time+recipe_.stage(n).abscissa*state_.dt;
        if(std::getenv("SF_SCALAR_TRACE")) {
            std::ostringstream trace;trace<<std::setprecision(17)<<"[SF SCALAR stage] target="<<contract_.target<<" step="<<state_.step<<" stage="<<n<<" time="<<stageTime_<<" dt="<<state_.dt;
            std::cout<<trace.str()<<'\n';
        }
        if(contract_.sequential) {
            for(const auto& read:contract_.stencilReads) {
                if(read.symbol==contract_.target || read.version!=System::StateVersion::Current)continue;
                const auto& fields=read.kind==System::TargetKind::Working
                    ?realized_.view(read.symbol,System::StateViewKind::Working).fields:realized_.at(read.symbol).fields;
                const auto& symbol=*realized_.at(read.symbol).descriptor;
                if(symbol.storageBinding==System::StorageBinding::TransientWorkspace) {
                    // The producing call owns final NearestInterior extension and publication.
                    if(fields.size()!=1 || !fields.front()->requireReadable)throw std::runtime_error("Missing published coefficient neighborhood: "+read.symbol);
                    fields.front()->requireReadable();continue;
                }
                if(fields.size()!=1 || !fields.front()->prepareSpatial)throw std::runtime_error("Missing selected boundary preparation for spatial input "+read.symbol);
                if(fields.front()->requireReadable)fields.front()->requireReadable();
                fields.front()->prepareSpatial(stageTime_);
            }
        }
        if(contract_.metric)for(const auto& read:contract_.stencilReads)if(read.version==System::StateVersion::Current && read.symbol!=contract_.target) {
            const auto& actual=realized_.at(read.symbol);const auto& fields=actual.fields;
            if(!actual.descriptor->producer && actual.descriptor->storageBinding==System::StorageBinding::NamedDistributed
                && std::any_of(fields.begin(),fields.end(),[](const auto* field){return !field->prepareSpatial;}))
                throw std::runtime_error("Metric spatial input lacks a Current boundary port: "+read.symbol);
            // Existing flow preparation owns packed/producer boundaries; named
            // input ports prepare their actual backing, never a substitute Q.
            for(auto* field:fields)if(field->prepareSpatial)field->prepareSpatial(stageTime_);
        }
        boundary(stageTime_);
        if(contract_.metric)updateBase_=stage_;
    }
    void evaluateRHS(const Run::ExecutionContext& context) {
        requireData();
        const int n=context.stageIndex;
        if(boundaryRevision_!=scratchRevision_ || boundaryTime_!=stageTime_)throw std::runtime_error("RHS before matching version/time boundary preparation");
        if(!contract_.volumeFlux.empty())if(const auto* dir=std::getenv("SF_SCALAR_FLUX_TRACE_DIR")) {
            std::filesystem::create_directories(dir);
            faceTrace_=std::make_shared<std::ofstream>(std::filesystem::path(dir)/("faces-"+std::to_string(state_.step)+".csv"));
            *faceTrace_<<std::setprecision(17)<<"axis,lower,upper,phi,upwind\n";
        }
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
            const auto at=index(cell);
            rhs_[entry(cell,component)]=(fixed(cell) || (contract_.metric && geometry_.CellFlag(at[0],at[1],at[2])!=FLUID_CELL))?0:contract_.rhs(values_,cell,component)
                -(contract_.conservativeTransport?advection(cell):0);
            if (!std::isfinite(rhs_[entry(cell,component)])) throw std::runtime_error("Scalar RHS is non-finite.");
            increments_[n][entry(cell,component)]=rhs_[entry(cell,component)];
        });
    }
    void advanceStage(const Run::ExecutionContext& context) {
        const int n=context.stageIndex;
        ownedValues([&](int cell,int component) {
            const auto location=index(cell);if(contract_.metric && geometry_.CellFlag(location[0],location[1],location[2])!=FLUID_CELL)return;
            const auto at=entry(cell,component);std::array<double,4> increments{};
            for(int i=0;i<recipe_.stageCount();++i)increments[i]=increments_[i][at];
            stage_[at]=explicitValue(recipe_,n,state_.dt,contract_.metric?updateBase_[at]:old_[at],stage_[at],rhs_[at],increments);
        });
        ++scratchRevision_;++nextStage_;if(contract_.metric)metricBoundary(true);
    }
    void fixedSnapshot() {fixedInterval_.open(state_,contract_.target);begin();fixedStep_=state_.step;fixedUpdated_=-1;}
    void loadIterate() {
        fixedInterval_.require(state_,contract_.target);
        if(fixedStep_!=state_.step)throw std::runtime_error("Missing/expired fixed-time snapshot: "+contract_.target);
        for(int cell=0;cell<geometry_.TotalSize();++cell)for(int c=0;c<contract_.components;++c)stage_[entry(cell,c)]=physical_->read(cell,c);
        ++scratchRevision_;nextStage_=0;Run::ExecutionContext context;context.stageCount=1;context.stageIndex=0;prepareStage(context);
    }
    void iterate() {
        frozenInputs_.clear();
        for(const auto& read:contract_.valueReads)if(read.version==System::StateVersion::Frozen) {
            const auto& actual=realized_.view(read.symbol,System::StateViewKind::Frozen);
            actual.fields.at(0)->requireReadable();frozenInputs_[read.symbol]=realized_.revision(read);
        }
        loadIterate();Run::ExecutionContext context;context.stageCount=1;context.stageIndex=0;evaluateRHS(context);
        ownedValues([&](int cell,int c){if(fixed(cell))return;auto at=entry(cell,c);stage_[at]=relaxSolution(stage_[at],fixedTimeCandidate(old_[at],state_.dt*rhs_[at]),contract_.relaxation);});
        ++scratchRevision_;nextStage_=1;validatePublish();publish();fixedUpdated_=state_.step;
    }
    void acceptFixed() {fixedInterval_.require(state_,contract_.target);if(fixedUpdated_!=state_.step)throw std::runtime_error("Fixed-time result not produced: "+contract_.target);}
    bool fixedSampleAvailable(System::ResidualKind kind) const {
        if(fixedStep_!=state_.step)return false;
        if(kind==System::ResidualKind::StateChange)return true;
        if(fixedUpdated_!=state_.step)return false;
        for(const auto& [symbol,revision]:frozenInputs_)if(realized_.revision(symbol)!=revision)return false;
        return true;
    }
    std::vector<double> fixedSamples(System::ResidualKind kind) {
        fixedInterval_.require(state_,contract_.target);
        if(kind!=System::ResidualKind::StateChange && kind!=System::ResidualKind::DiscreteEquation)throw std::runtime_error("Unsupported fixed-time residual kind");
        if(kind==System::ResidualKind::DiscreteEquation){loadIterate();Run::ExecutionContext c;c.stageCount=1;c.stageIndex=0;evaluateRHS(c);}
        std::vector<double> result;ownedValues([&](int cell,int c){if(fixed(cell))return;const auto at=entry(cell,c);const double a=physical_->read(cell,c);result.push_back(kind==System::ResidualKind::StateChange?a:a-old_[at]-state_.dt*rhs_[at]);});return result;
    }
    void validatePublish() {
        requireData();
        if (nextStage_!=recipe_.stageCount()) throw std::runtime_error("Scalar cannot commit incomplete stages.");
        if(!contract_.metric)boundary(state_.time+state_.dt);
        for(int component=0;component<contract_.components;++component) for (int cell=0;cell<geometry_.TotalSize();++cell) {
            if (!std::isfinite(stage_[entry(cell,component)])) throw std::runtime_error("Scalar update is non-finite.");
        }
        if(contract_.primitiveStorage)ownedValues([&](int cell,int component){if(stage_[entry(cell,component)]<=0)throw std::runtime_error("Native stage scalar physical update nonpositive; no repair.");});
    }
    void publish() {
        if(faceTrace_)faceTrace_->flush();
        if(!contract_.volumeFlux.empty())if(const auto* dir=std::getenv("SF_SCALAR_FLUX_TRACE_DIR")) {
            std::ofstream out(std::filesystem::path(dir)/("cells-"+std::to_string(state_.step)+".csv"));
            out<<std::setprecision(17)<<"cell,old,rhs,new,volume,dt\n";
            owned([&](int cell){if(fixed(cell))return;auto at=index(cell);double volume=1/geometry_.Jac(at[0],at[1],at[2]);
                for(int axis=0;axis<3;++axis)if(!active_[axis] && shape_[axis]>1 && (at[axis]==geometry_.NG() || at[axis]==geometry_.NG()+shape_[axis]-1))volume*=.5;
                out<<cell<<','<<old_[cell]<<','<<rhs_[cell]<<','<<stage_[cell]<<','<<volume<<','<<state_.dt<<'\n';});
        }
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
    System::StateRealization& realized_;
    int stages_;
    const System::CompiledEquationCall& call_;
    const System::CompiledSstClosure& contract_;
    std::vector<std::vector<double>> values_;
    std::vector<State::DistributedFieldView> physical_,staged_,frozen_,owners_;
    std::vector<std::function<double(int,int)>> reads_;
    std::vector<unsigned long long> revisions_;
    struct Publication {
        int step=-1,stage=-1;double time=0,dt=0;bool physical=false;
        void require(const State::StateBundle& state,bool current,const std::string& input,int index=-1) const {
            const bool committed=current && physical && step+1==state.step && time+dt==state.time;
            if(!committed && (step!=state.step || dt!=state.dt || physical!=current || (current?time!=state.time:stage!=index)))
                throw std::runtime_error("Unpublished closure input="+input+" view="+(current?"Current":"Stage"));
        }
    } publication_;
    void requireInputs() const {
        for(std::size_t i=0;i<revisions_.size();++i)if(realized_.revision(call_.stateUses[i])!=revisions_[i])
            throw std::runtime_error("Closure input changed since publication: "+call_.source.occurrence+" / "+call_.stateUses[i].symbol);
    }
public:
    SstClosureOperations(State::StateBundle& state,Field& field,System::StateRealization& realized,int stages,
            const System::CompiledSstClosure& contract,const System::CompiledEquationCall& call)
        :state_(state),field_(field),realized_(realized),stages_(stages),call_(call),contract_(contract) {
        if(contract.packedVelocity && field.NVar()<5)throw std::runtime_error("SST packed primitive gradient needs its declared five-component backing");
        for(const auto& id:{contract.density,contract.momentum}) {
            if(realized.at(id).descriptor->constantValue)continue;
            const auto& backing=realized.at(id).fields;
            if(backing.size()!=1 || backing.front()->geometry!=&field)throw std::runtime_error("SST input geometry differs: "+id);
        }
        const auto count=contract.outputs.size();values_.resize(count);physical_.resize(count);owners_.resize(count);frozen_.resize(count);staged_.reserve(count*stages);
        for(std::size_t n=0;n<count;++n) {
            const auto& name=contract.outputs[n];values_[n].resize(field.TotalSize());
            auto view=State::workspaceView(name,0,field,values_[n],1);
            if(realized.at(name).descriptor->storageBinding!=System::StorageBinding::TransientWorkspace) {
                const auto& backing=realized.at(name).fields;
                if(backing.size()!=1 || backing.front()->geometry!=&field || backing.front()->components!=1 || !backing.front()->write)
                    throw std::runtime_error("Closure output lacks actual writable scalar backing: "+name);
                owners_[n]=*backing.front();view=owners_[n];
            }
            for(int stage=0;call.source.stageScoped && stage<stages;++stage) {
                auto stageView=view;
                stageView.requireReadable=[this,name,stage]{realized_.requireStageRead(stage,publication_.time);publication_.require(state_,false,name,stage);};
                stageView.read=[read=view.read,guard=stageView.requireReadable](int c,int k){guard();return read(c,k);};
                stageView.write=[](int,int,double){throw std::runtime_error("Closure Stage result is read-only");};
                staged_.push_back(std::move(stageView));realized.bindView(name,System::StateViewKind::Stage,staged_.back(),stage);
            }
            if(realized.requestsView(name,System::StateViewKind::Frozen)) {
                auto frozen=view;frozen.requireReadable=[this,name]{publication_.require(state_,true,name);};
                frozen.read=[read=view.read,guard=frozen.requireReadable](int c,int k){guard();return read(c,k);};
                frozen.write=[](int,int,double){throw std::runtime_error("Frozen closure result is read-only");};
                frozen.revision=[this,name]{return realized_.revision(name);};
                frozen_[n]=std::move(frozen);realized.bindView(name,System::StateViewKind::Frozen,frozen_[n]);
            }
            view.requireReadable=[this,name]{publication_.require(state_,true,name);if(call_.source.stepScoped)requireInputs();};
            view.read=[read=view.read,guard=view.requireReadable](int c,int k){guard();return read(c,k);};
            view.write=[](int,int,double){throw std::runtime_error("Closure result is read-only outside its producing call");};
            physical_[n]=std::move(view);
            if(owners_[n].write)*realized.at(name).fields.front()=physical_[n];
            else realized.bindView(name,System::StateViewKind::Physical,physical_[n]);
        }
    }
    void bindInputs() {
        for(const auto& use:call_.stateUses)reads_.push_back(realized_.reader(use,stages_));
    }
    void refresh(bool physical,int stage=-1,double time=0) {
        if(!physical && publication_.step==state_.step && !publication_.physical && publication_.stage==stage)
            throw std::runtime_error("Repeated closure Stage evaluation requires explicit iterative stage treatment: "+call_.source.occurrence);
        publication_.step=-1;
        const bool final=physical && call_.source.stageScoped;
        const double evaluationTime=final?state_.time+state_.dt:physical?state_.time:time;
        auto inputs=reads_;
        if(final)for(std::size_t n=0;n<inputs.size();++n)
            inputs[n]=realized_.reader({call_.stateUses[n].symbol,System::StateVersion::Current},stages_);
        if(!contract_.packedVelocity)for(std::size_t n=1;n<4;++n) {
            const auto& use=call_.stateUses[n];
            if(!physical && use.version==System::StateVersion::Stage)continue;
            if(final && use.version==System::StateVersion::Stage && realized_.at(use.symbol).descriptor->storageBinding==System::StorageBinding::ProviderDistributed)continue;
            auto preparedUse=use;if(final){preparedUse.version=System::StateVersion::Current;preparedUse.kind=System::TargetKind::Physical;}
            realized_.prepareSpatial(preparedUse,evaluationTime,stage);
        }
        Discretization::SstSamples samples;
        for(int n=0;n<2;++n)samples[n]=[this,read=inputs[n+2]](int i,int j,int k) {
            const int g=field_.NG();
            const int x=field_.NX()>1 && Math::isDirectionActiveIndex(0)?std::clamp(i,g-1,g+field_.NX()):g;
            const int y=field_.NY()>1 && Math::isDirectionActiveIndex(1)?std::clamp(j,g-1,g+field_.NY()):g;
            const int z=field_.NZ()>1 && Math::isDirectionActiveIndex(2)?std::clamp(k,g-1,g+field_.NZ()):g;
            return read(field_.getIdx(x,y,z),0);
        };
        const auto density=[&](int i,int j,int k){return inputs[0](field_.getIdx(i,j,k),0);};
        const auto gradient=[&](int i,int j,int k) {
            if(contract_.packedVelocity)return Viscous::velocityGradientAt(field_,i,j,k,Viscous::CentralOrder::SECOND);
            Tensor3 result;for(int c=0;c<3;++c)result[c]=CENTRAL2::gradSampled(field_,[&](int x,int y,int z){return inputs[1](field_.getIdx(x,y,z),c);},i,j,k);return result;
        };
        const int g=field_.NG();
        for(int k=g;k<g+field_.NZ();++k)for(int j=g;j<g+field_.NY();++j)for(int i=g;i<g+field_.NX();++i) {
            const auto distance=[this](int i,int j,int k){const double d=field_.wallDistance(i,j,k);return std::isfinite(d) && d>1e-14?d:1e20;};
            const int cell=field_.getIdx(i,j,k);
            std::array<double,5> output{};
            if(contract_.eddyViscosity)output[1]=Discretization::sstEddyViscosity(density(i,j,k),samples[0](i,j,k),samples[1](i,j,k),distance(i,j,k),gradient(i,j,k),contract_.parameters);
            else {
                const double viscosity=inputs[4](cell,0);
                if(!std::isfinite(viscosity) || viscosity<0)throw std::runtime_error("SST production reads invalid published viscosity");
                const auto c=Discretization::sstCentral2(field_,samples,contract_.parameters,density,gradient,distance,viscosity,i,j,k);
                output={c.blending1,viscosity,c.production[0],c.production[1],c.crossSource};
            }
            for(std::size_t n=0;n<values_.size();++n) {
                const double value=output[contract_.outputRoles[n]];
                if(!std::isfinite(value) || (contract_.eddyViscosity && value<0))throw std::runtime_error("Nonfinite closure output: "+contract_.outputs[n]);
                values_[n][cell]=value;
            }
        }
        for(std::size_t n=0;n<values_.size();++n) {
            extendNearestInterior(field_,values_[n]);
            if(owners_[n].write) {
                for(int cell=0;cell<field_.TotalSize();++cell)owners_[n].write(cell,0,values_[n][cell]);
                state_.distributed.markModified(owners_[n].name);
            } else realized_.publishWorkspace(contract_.outputs[n]);
        }
        revisions_.clear();for(const auto& use:call_.stateUses)revisions_.push_back(realized_.revision(use));
        publication_={state_.step,stage,physical?state_.time:time,state_.dt,physical};
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

/// Owns matrix/solve workspace only. Current STATE remains in its original storage.
class BlockOperations {
    const System::CompiledEquationCall& call_;
    const System::CompiledLinearBlock& contract_;
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
        :call_(call),contract_(std::any_cast<const System::CompiledLinearBlock&>(call.providerContract)),state_(state),realized_(realized),mesh_(state.singlePatch()),shape_{mesh_.NX(),mesh_.NY(),mesh_.NZ()} {
        if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))throw std::runtime_error("Unsupported linear block backend: serial single patch required");
        int g=mesh_.NG();local_.assign(mesh_.TotalSize(),-1);
        for(int k=g;k<g+shape_[2];++k)for(int j=g;j<g+shape_[1];++j)for(int i=g;i<g+shape_[0];++i){int cell=mesh_.getIdx(i,j,k);local_[cell]=cells_.size();cells_.push_back(cell);}
        const auto coordinate=[&](int axis,int cell){auto at=index(cell);return axis==0?mesh_.X(at[0],at[1],at[2]):axis==1?mesh_.Y(at[0],at[1],at[2]):mesh_.Z(at[0],at[1],at[2]);};
        int origin=cells_.front();auto start=index(origin);
        for(int axis=0;axis<3;++axis)if(shape_[axis]>1){if(shape_[axis]<3)throw std::runtime_error("Block Central2 requires at least three nodes per active direction");auto next=start;++next[axis];double h=coordinate(axis,mesh_.getIdx(next[0],next[1],next[2]))-coordinate(axis,origin);if(!std::isfinite(h) || h<=0)throw std::runtime_error("Unsupported block Cartesian spacing");spacing_[axis]=h;for(int cell:cells_){double expected=coordinate(axis,origin)+(index(cell)[axis]-g)*h;if(std::abs(coordinate(axis,cell)-expected)>1e-11*std::max(1.,std::abs(expected)))throw std::runtime_error("Unsupported nonuniform/nonorthogonal block geometry");}}
        for(std::size_t n=0;n<contract_.targets.size();++n){const auto& target=contract_.targets[n];auto& field=realized_.at(target);if(field.fields.size()!=1)throw std::runtime_error("Unbound block target: "+target);auto* view=field.fields[0];if(view->geometry!=&mesh_ || view->components!=1 || !view->read || !view->write || (field.descriptor->storageBinding==System::StorageBinding::ProviderDistributed && !view->requireReadable))throw std::runtime_error("Block target layout/access/lifetime capability missing: "+target);physical_.push_back(view);const auto& id=call_.physicalResources.at(n).id;ports_.push_back(equationData(services,id).data);auto data=port(n);if(data->target!=target)throw std::runtime_error("Block boundary port target mismatch");for(int axis=0;axis<3;++axis)if(spacing_[axis])for(int side=0;side<2;++side)if(data->boundary[2*axis+side].kind!=FDM::ScalarBoundaryKind::FixedValue && data->boundary[2*axis+side].kind!=FDM::ScalarBoundaryKind::ZeroGradient)throw std::runtime_error("Unsupported block boundary: existing fixedValue/zeroGradient required");}
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
    const System::CompiledAffineEquation& contract_;
    State::StateBundle& state_;
    System::StateRealization& realized_;
    Field& mesh_;
    State::DistributedFieldView* output_;
    State::DistributedFieldView materialized_,frozen_;
    bool workspace_=false,publicationRequired_=false,physicalReady_=false;
    int publishedStep_=-1;
    std::vector<std::pair<System::StateUse,unsigned long long>> inputRevisions_;
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
            throw std::runtime_error("Unprepared coefficient read region: step="+call_.source.occurrence+" equation="+call_.source.equation+" target="+call_.target.symbol+" cell="+std::to_string(cell));
    }
    void requirePublication() const {
        if(!physicalReady_ || publishedStep_!=state_.step)
            throw std::runtime_error("Unpublished materialized STATE in this step: "+call_.target.symbol+"; execute its HOW relation first");
        for(const auto& input:inputRevisions_) {
            const auto& use=input.first;
            if(use.kind==System::TargetKind::Physical) {
                const auto& backing=realized_.at(use.symbol);
                if(backing.fields.size()!=1)throw std::runtime_error("Missing materialized dependency backing: "+use.symbol);
                if(backing.fields.front()->requireReadable)backing.fields.front()->requireReadable();
            }
            if(realized_.revision(use)!=input.second)
                throw std::runtime_error("Stale materialized dependency: "+call_.target.symbol+" reads "+use.symbol+" in its declared view; HOW must recompute the relation");
        }
    }
    void calculate(const System::FormulaValues& values,std::vector<double>& output){
        if(!inputsBound_)throw std::runtime_error("Affine inputs have not completed startup binding");
        if(!workspace_ && port_.expired())throw std::runtime_error("Expired algebraic resource port");
        owned([&](int cell){double A=contract_.coefficient(values,cell,0),R=contract_.remainder(values,cell,0);if(!std::isfinite(A) || !std::isfinite(R) || A==0)throw std::runtime_error("LocalAffine singular/non-finite coefficient: equation="+call_.source.equation+" target="+call_.target.symbol);double x=-R/A;if(!std::isfinite(x))throw std::runtime_error("Non-finite LocalAffine solution");output.at(cell)=x;});
        if(contract_.boundaryPreparation==System::WorkspaceBoundary::NearestInterior)extendNearestInterior(mesh_,output);
    }
    void write(){
        if(!workspace_ && output_->requireReadable)output_->requireReadable();
        if(contract_.boundaryPreparation==System::WorkspaceBoundary::NearestInterior)
            for(int cell=0;cell<mesh_.TotalSize();++cell)output_->write(cell,0,pending_.at(cell));
        else owned([&](int cell){output_->write(cell,0,pending_.at(cell));});
        physicalReady_=true;
        publishedStep_=state_.step;inputRevisions_.clear();
        for(const auto& input:call_.stateUses)if(input.version==System::StateVersion::Current)
            inputRevisions_.push_back({input,realized_.revision(input)});
        if(workspace_)realized_.publishWorkspace(call_.target.symbol);
        if(call_.target.kind==System::TargetKind::Physical && !workspace_)state_.distributed.markModified(output_->name,output_->blockId);
        else if(call_.target.kind!=System::TargetKind::Physical) {
            publications_->insert(call_.target.workspace);
            realized_.publishWorkspace(call_.target.symbol,call_.target.kind==System::TargetKind::Working?System::StateViewKind::Working:System::StateViewKind::Correction);
        }
    }
public:
    AffineOperations(const System::CompiledEquationCall& call,const System::CompiledTimeRecipe& recipe,State::StateBundle& state,System::StateRealization& realized,const FDM::ScalarInstanceData* instance,std::shared_ptr<std::set<std::string>> publications)
        :call_(call),recipe_(recipe),contract_(std::any_cast<const System::CompiledAffineEquation&>(call.providerContract)),state_(state),realized_(realized),mesh_(state.singlePatch()),publications_(std::move(publications)){
        const auto& target=*realized.at(call_.target.symbol).descriptor;
        workspace_=target.storageBinding==System::StorageBinding::TransientWorkspace;
        publicationRequired_=workspace_ || (call.target.kind==System::TargetKind::Physical
            && target.storageBinding==System::StorageBinding::ProviderDistributed
            && target.evaluation==System::StateEvaluation::Materialized);
        if(contract_.boundaryPreparation!=realized.at(call_.target.symbol).descriptor->workspaceBoundary
            || (contract_.boundaryPreparation!=System::WorkspaceBoundary::None && contract_.boundaryPreparation!=System::WorkspaceBoundary::NearestInterior))
            throw std::runtime_error("Affine output read-region contract differs from its STATE declaration: "+call_.target.symbol);
        pending_.resize(mesh_.TotalSize());
        if(workspace_) {
            if(instance)throw std::runtime_error("Duplicate resource authority for affine workspace: "+call_.target.symbol);
            materializedValues_.resize(mesh_.TotalSize());
            materialized_=State::workspaceView(realized.at(call_.target.symbol).descriptor->storageKey,0,mesh_,materializedValues_,1);
            materialized_.read=[this](int cell,int){requireRegion(cell);return materializedValues_.at(cell);};
            materialized_.requireReadable=[this]{requirePublication();};
            realized.bindView(call_.target.symbol,System::StateViewKind::Physical,materialized_);
            if(realized.requestsView(call_.target.symbol,System::StateViewKind::Frozen)) {
                frozen_=materialized_;
                frozen_.requireReadable=[this]{
                    if(!physicalReady_ || publishedStep_!=state_.step)throw std::runtime_error("Unpublished Frozen relation: "+call_.target.symbol);
                };
                frozen_.revision=[this]{return realized_.revision(call_.target.symbol);};
                frozen_.write=[](int,int,double){throw std::runtime_error("Frozen published relation is read-only");};
                realized.bindView(call_.target.symbol,System::StateViewKind::Frozen,frozen_);
            }
        } else {
            if(!instance)throw std::runtime_error("Missing affine physical resource port: "+call_.target.symbol);
            port_=instance->data;
        }
        output_=call.target.kind==System::TargetKind::Physical?realized.at(call_.target.symbol).fields.at(0)
            :realized.view(call_.target.symbol,call.target.kind==System::TargetKind::Working?System::StateViewKind::Working:System::StateViewKind::Correction).fields.at(0);
        const auto port=port_.lock();if((!workspace_ && (!port || port->target!=call_.target.symbol)) || output_->geometry!=&mesh_ || output_->components!=1 || !output_->write)throw std::runtime_error("Invalid affine STATE resource/target binding");
        if(publicationRequired_ && !workspace_) {
            // Keep the original write port while consumers see the published
            // view of the same owner array.
            materialized_=*output_;
            output_->requireReadable=[this,guard=materialized_.requireReadable] {
                if(guard)guard();
                requirePublication();
            };
            output_->read=[this,read=materialized_.read,guard=output_->requireReadable](int cell,int component) {
                guard();
                requireRegion(cell);
                return read(cell,component);
            };
            output_=&materialized_;
        }
        if(call_.stageLocal){stage_.resize(mesh_.TotalSize());views_.reserve(recipe_.stageCount());for(int stage=0;stage<recipe_.stageCount();++stage){views_.push_back(State::workspaceView(call_.stageWorkspace,0,mesh_,stage_,1));views_.back().read=[this,stage](int cell,int){double time=state_.time+recipe_.stage(stage).abscissa*state_.dt;realized_.requireStageRead(stage,time);if(readyStage_!=stage || readyStep_!=state_.step || readyTime_!=time)throw std::runtime_error("Unpublished algebraic Stage: "+call_.target.symbol);requireRegion(cell);return stage_.at(cell);};realized_.bindView(call_.target.symbol,System::StateViewKind::Stage,views_.back(),stage);}}
    }
    void bindReads(){
        if(inputsBound_)throw std::runtime_error("Affine input readers were already bound");
        for(const auto& read:call_.stateUses){
            if(read.kind==System::TargetKind::Working && realized_.view(read.symbol,System::StateViewKind::Working).descriptor.owner==System::StateViewOwner::NumericalProvider) {
                auto reader=realized_.reader(read,recipe_.stageCount());
                values_.boundReads.push_back(reader);physicalValues_.boundReads.push_back(std::move(reader));continue;
            }
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
    bool sampleAvailable() const {return !publicationRequired_ || (physicalReady_ && publishedStep_==state_.step);}
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
        if(!stageReady) {
            if(publicationRequired_)requirePublication();
            else if(output_->requireReadable)output_->requireReadable();
        }
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
void EquationInterval::open(const State::StateBundle& state,const std::string& label) {
    if(step==state.step || !std::isfinite(state.dt) || state.dt<=0)throw std::runtime_error("Missing/new physical interval: "+label);
    step=state.step;time=state.time;dt=state.dt;
}
void EquationInterval::require(const State::StateBundle& state,const std::string& label) const {
    if(step!=state.step || time!=state.time || dt!=state.dt || !(dt>0))throw std::runtime_error("Expired physical interval: "+label);
}
State::DistributedFieldView& equationField(System::StateRealization& realized,const System::StateUse& use,
        Field& geometry,bool writable,bool unsliced) {
    const auto& base=realized.at(use.symbol);const auto& descriptor=*base.descriptor;
    if(use.version!=System::StateVersion::Current || (use.kind!=System::TargetKind::Physical && use.kind!=System::TargetKind::Working))
        throw std::runtime_error("Ordinary field port requires implemented Current Physical/Working backing: "+use.symbol);
    const auto& fields=use.kind==System::TargetKind::Physical?base.fields:realized.view(use.symbol,System::StateViewKind::Working).fields;
    const int offset=use.kind==System::TargetKind::Physical?descriptor.componentOffset:0;
    if(fields.size()!=1 || !fields[0] || fields[0]->geometry!=&geometry || !fields[0]->read || (writable && !fields[0]->write)
        || fields[0]->components<offset+descriptor.components || (unsliced && (offset || fields[0]->components!=descriptor.components)))
        throw std::runtime_error("Ordinary field geometry/shape/access/offset mismatch: "+use.symbol);
    return *fields[0];
}
const FDM::ScalarInstanceData& equationData(const FDM::SolverServices& services,const std::string& id) {
    const FDM::ScalarInstanceData* result=nullptr;
    for(const auto& entry:services.scalarInstances)if(entry.occurrence==id) {
        if(result)throw std::runtime_error("Ambiguous ordinary field resource: "+id);result=&entry;
    }
    if(!result || result->data.expired())throw std::runtime_error("Missing/expired ordinary field resource: "+id);
    return *result;
}
BoundEquation bindScalar(const EquationBindingInput& input) {
    const auto& call=input.call;
    const auto& numerical=input.numerical;
    auto& state=input.state;
    const auto& services=input.services;
    auto& realized=input.views;
    const auto& limit=input.limit;
    const auto& residuals=input.residuals;
    const auto& scalar=std::any_cast<const System::CompiledScalarEquation&>(call.providerContract);
    if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
        throw std::runtime_error("Scalar execution requires one serial patch: "+call.source.occurrence);
    const auto* port=services.scalarEquation && !scalar.sequential?nullptr:&equationData(services,fieldPort(call));
    const auto data=port?port->data.lock():std::shared_ptr<const FDM::ScalarEquationData>{};
    if(port && !data)throw std::runtime_error("Missing scalar instance data: "+call.source.occurrence);
    const auto checked=[step=&call,&state](const char* phase,int stage,auto&& operation)->decltype(auto) {
        try {return operation();}
        catch(const std::exception& error) {
            throw std::runtime_error("Scalar operation="+std::string(phase)+" occurrence="+step->source.occurrence
                +" equation="+step->source.equation+" target="+step->target.symbol+" provider="+step->backendProvider
                +" step="+std::to_string(state.step)+" stage="+std::to_string(stage)+" patch=0: "+error.what());
        }
    };
    const auto owner=checked("STATE binding",-1,[&] {
        return std::make_shared<ScalarOperations>(numerical,call,state,port?*data:*services.scalarEquation,realized,limit,
            std::any_of(call.physicalResources.begin(),call.physicalResources.end(),[](const auto& resource) {
                return resource.fields && resource.fields->boundaryRequired;
            }));
    });
    if(port)owner->trackData(port->data);
    BoundEquation bound;
    if(!port)bound.singleScalarPort=services.scalarEquation;
    bound.connect=[owner](const FluxReaders& fluxes){owner->bindFluxes(fluxes);owner->bindInputs();};
    if(!scalar.sequential) {
        auto& cb=bound.stage;
        cb.identity=call.source.occurrence;
        cb.stepSize=[owner,checked]{return checked("dt",-1,[&]{return owner->stepSize();});};
        cb.snapshot=[owner,checked]{checked("snapshot",-1,[&]{owner->begin();});};
        cb.prepare=[owner,checked](const auto& c){checked("prepare",c.stageIndex,[&]{owner->prepareStage(c);});};
        cb.rhs=[owner,checked](const auto& c){checked("rhs",c.stageIndex,[&]{owner->evaluateRHS(c);});};
        cb.advance=[owner,checked](const auto& c){checked("advance",c.stageIndex,[&]{owner->advanceStage(c);});};
        cb.validatePublish=[owner,checked]{checked("validate publish",-1,[&]{owner->validatePublish();});};
        cb.publish=[owner,checked]{checked("publish",-1,[&]{owner->publish();});};
        return bound;
    }
    bound.interval.stepSize=[owner]{return owner->stepSize();};
    if(scalar.fixedTime) {
        bound.interval.snapshot=[owner]{owner->fixedSnapshot();};
        bound.interval.accept=[owner]{owner->acceptFixed();};
        bound.interval.publishesOldTime=true;
        bound.install=[owner,step=&call,residuals](Run::OpRegistry& ops) {
            ops.bind(step->backendOperation,step->backendProvider,[owner]{owner->iterate();});
            for(const auto* residual:residuals)
                ops.bindResidual(residual->id,residual->provider,[owner,kind=residual->kind](const auto&){return owner->fixedSamples(kind);},
                    [owner,kind=residual->kind](const auto&){return owner->fixedSampleAvailable(kind);});
        };
        return bound;
    }
    bound.install=[owner,step=&call,&state](Run::OpRegistry& ops) {
        auto interval=std::make_shared<EquationInterval>();
        auto next=std::make_shared<int>(0);
        for(int phase=0;phase<5;++phase) {
            const auto& leaf=step->fragment.children.at(phase);
            ops.bind(leaf.operation,leaf.provider,[owner,interval,next,&state,phase,label=step->source.occurrence] {
                if(phase==0){interval->open(state,label);*next=0;}
                interval->require(state,label);
                if(*next!=phase)throw std::runtime_error("Out-of-order scalar operation: "+label);
                Run::ExecutionContext context;context.stageIndex=0;context.stageCount=1;
                switch(phase) {
                case 0:owner->begin();break;
                case 1:owner->prepareStage(context);break;
                case 2:owner->evaluateRHS(context);break;
                case 3:owner->advanceStage(context);break;
                case 4:owner->validatePublish();owner->publish();break;
                }
                ++*next;
            });
        }
    };
    return bound;
}

BoundEquation bindSst(const EquationBindingInput& input) {
    const auto& call=input.call;
    const auto& numerical=input.numerical;
    auto& state=input.state;
    const auto& services=input.services;
    auto& realized=input.views;
    if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed())
        || services.immersed.boundary || services.immersed.constraint || services.immersed.system)
        throw std::runtime_error("Unsupported closure stage treatment: step="+call.source.occurrence+" requires serial single patch without IBM");
    const auto& recipe=numerical.time.recipe;
    auto owner=std::make_shared<SstClosureOperations>(state,state.singlePatch(),realized,recipe.stageCount(),
        std::any_cast<const System::CompiledSstClosure&>(call.providerContract),call);
    BoundEquation bound;
    bound.connect=[owner](const FluxReaders&){owner->bindInputs();};
    bound.install=[owner,step=&call,&state,&recipe](Run::OpRegistry& ops) {
        ops.bind(step->backendOperation,step->backendProvider,[owner,&state,&recipe,physical=step->source.stepScoped](const Run::ExecutionContext& c) {
            if(physical)owner->refresh(true);
            else owner->refresh(false,c.stageIndex,state.time+recipe.stage(c.stageIndex).abscissa*state.dt);
        });
        if(!step->temporalPublish.empty())ops.bind(step->temporalPublish,step->backendProvider,[owner]{owner->refresh(true);});
    };
    return bound;
}

BoundEquation bindAffine(const EquationBindingInput& input,const std::shared_ptr<std::set<std::string>>& publications) {
    const auto& call=input.call;
    const auto& numerical=input.numerical;
    auto& state=input.state;
    const auto& services=input.services;
    auto& realized=input.views;
    const auto& residuals=input.residuals;
    if(state.patches.size()!=1 || (services.executionRuntime && services.executionRuntime->distributed()))
        throw std::runtime_error("Unsupported LocalAffine layout: one serial patch required: "+call.source.occurrence);
    const auto* instance=call.physicalResources.empty()?nullptr:&equationData(services,fieldPort(call));
    auto owner=std::make_shared<AffineOperations>(call,numerical.time.recipe,state,realized,instance,publications);
    const auto checked=[label="LocalAffine step="+call.source.occurrence+" equation="+call.source.equation+" target="+call.target.symbol+": "](auto&& action) {
        try{action();}catch(const std::exception& e){throw std::runtime_error(label+e.what());}
    };
    BoundEquation bound;
    bound.connect=[owner,checked](const FluxReaders&){checked([&]{owner->bindReads();});};
    bound.install=[owner,checked,step=&call,residuals](Run::OpRegistry& ops) {
        ops.bind(step->backendOperation,step->backendProvider,[owner,checked](const Run::ExecutionContext& c){checked([&]{owner->solve(c);});});
        if(step->stageLocal)ops.bind(step->temporalPublish,step->backendProvider,[owner,checked]{checked([&]{owner->publish();});});
        for(const auto* residual:residuals)
            ops.bindResidual(residual->id,residual->provider,[owner,kind=residual->kind](const Run::ExecutionContext& c){return owner->samples(kind,c);},
                [owner](const auto&){return owner->sampleAvailable();});
    };
    return bound;
}

BoundEquation bindBlock(const EquationBindingInput& input) {
    const auto& call=input.call;
    auto owner=std::make_shared<BlockOperations>(call,input.state,input.views,input.services);
    std::string label="Block occurrence="+call.source.occurrence+" provider="+call.equationMethod+" equations/targets=[";
    for(const auto& relation:call.calls)label+=" "+relation.equation+"->"+relation.target;
    label+=" ]: ";
    const auto checked=[owner,label](auto operation) {
        return [owner,label,operation]{try{((*owner).*operation)();}catch(const std::exception& e){throw std::runtime_error(label+e.what());}};
    };
    BoundEquation bound;
    bound.connect=[owner](const FluxReaders&){owner->bindReads();};
    bound.install=[checked,step=&call](Run::OpRegistry& ops) {
        ops.bind(step->operations.at(0).operation,step->backendProvider,checked(&BlockOperations::assemble));
        ops.bind(step->operations.at(1).operation,step->backendProvider,checked(&BlockOperations::solve));
        ops.bind(step->operations.at(2).operation,step->backendProvider,checked(&BlockOperations::publish));
    };
    return bound;
}

void bindVariableCommits(Run::OpRegistry& ops,const System::ExecutionScope& root,
        State::StateBundle& state,System::StateRealization& realized) {
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
            ops.bind(System::instanceOperation("state.commit",scope.id),"execution.state",[copies,&state]{
                auto& mesh=state.singlePatch();const int g=mesh.NG();std::vector<std::vector<double>> pending(copies.size());
                for(std::size_t n=0;n<copies.size();++n) {
                    const auto& copy=copies[n];
                    if(copy.source.requireReadable)copy.source.requireReadable();if(copy.destination.requireReadable)copy.destination.requireReadable();
                    for(int k=g;k<g+mesh.NZ();++k)for(int j=g;j<g+mesh.NY();++j)for(int i=g;i<g+mesh.NX();++i)for(int c=0;c<copy.components;++c){double value=copy.source.read(mesh.getIdx(i,j,k),c);if(!std::isfinite(value))throw std::runtime_error("Nonfinite variable Commit; physical targets unmodified");pending[n].push_back(value);}
                }
                for(std::size_t n=0;n<copies.size();++n){const auto& copy=copies[n];std::size_t at=0;for(int k=g;k<g+mesh.NZ();++k)for(int j=g;j<g+mesh.NY();++j)for(int i=g;i<g+mesh.NX();++i)for(int c=0;c<copy.components;++c)copy.destination.write(mesh.getIdx(i,j,k),c,pending[n].at(at++));state.distributed.markModified(copy.destination.name,copy.destination.blockId);}
            });
        }
        for(const auto& child:scope.children)self(self,child);
    };bindCommit(bindCommit,root);
}

}
