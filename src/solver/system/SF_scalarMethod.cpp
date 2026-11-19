#include "SF_builtinState.h"
#include "SF_scalarMethod.h"
#include "SF_formulaCompiler.h"
#include "SF_numericalSystem.h"
#include "SF_stateRealization.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "solver/discretization/convection/SF_scalarUpwind.h"
#include "core/system/SF_operationIds.h"
#include <cmath>
#include <stdexcept>
namespace SF::System {
namespace {
FDM::ScalarBoundaryOptions metricOptions(const NumericalBinding& binding) {
    FDM::ScalarBoundaryOptions options;
    for(const auto& [key,value]:binding.parameters) {
        if(!std::isfinite(value))throw std::runtime_error("Non-finite metric scalar parameter: "+key);
        if(key=="bounded" || key=="ilw") {
            if(value!=0 && value!=1)throw std::runtime_error("Metric scalar flag must be 0/1: "+key);
            if(key=="bounded")options.bounds.enabled=value;
            else options.ilw=value;
        } else if(key=="lower")options.bounds.lower=value;
        else if(key=="upper")options.bounds.upper=value;
        else if(key=="boundMode") {
            if(value!=0 && value!=1 && value!=2)throw std::runtime_error("Unknown metric scalar bound mode");
            options.bounds.mode=static_cast<FDM::Scalar::BoundMode>((int)value);
        } else if(key=="ilwOrder") {
            if(value!=3 && value!=5 && value!=7 && value!=9 && value!=0)
                throw std::runtime_error("Metric scalar ILW order must be 3/5/7/9, or 0 when disabled");
            options.order=(int)value;
        } else throw std::runtime_error("Unconsumed metric scalar parameter: "+key);
    }
    if(options.bounds.lower>options.bounds.upper || (options.ilw && !options.order))
        throw std::runtime_error("Invalid metric scalar bounds/ILW choice");
    return options;
}
double fixedRelaxation(const NumericalBinding& binding,const std::string& occurrence) {
    double relaxation=1.0;
    for(const auto& [key,value]:binding.parameters) {
        if(key!="relaxation" || !std::isfinite(value))
            throw std::runtime_error("Fixed-point method requires a finite relaxation parameter: "+occurrence);
        relaxation=value;
    }
    return relaxation;
}
enum class ScalarImplementation {
    Diffusion, ArithmeticDiffusion, VelocityTransport, VolumeTransport,
    ConservativeTransport, SharedTransport, CapacityTransport,
    FixedPoint, FixedPointVolume, MetricTransport
};
class ScalarMethod final {
    ScalarImplementation implementation_;
public:
    explicit ScalarMethod(ScalarImplementation implementation):implementation_(implementation) {}
    std::string_view id() const {
        switch(implementation_) {
        case ScalarImplementation::Diffusion:return "ScalarDiffusionCentral2";
        case ScalarImplementation::ArithmeticDiffusion:return "ExplicitCentral2Arithmetic";
        case ScalarImplementation::VelocityTransport:return "ScalarAdvectionCentral2";
        case ScalarImplementation::VolumeTransport:return "ScalarVolumeFluxUpwind1";
        case ScalarImplementation::ConservativeTransport:return "ScalarTransportUpwindCentral2";
        case ScalarImplementation::SharedTransport:return "ScalarTransportSharedFluxCentral2";
        case ScalarImplementation::CapacityTransport:return "CapacityTransportSharedFluxCentral2";
        case ScalarImplementation::FixedPoint:return "FixedPointCentral2";
        case ScalarImplementation::FixedPointVolume:return "FixedPointVolumeFluxUpwind1";
        case ScalarImplementation::MetricTransport:return "ScalarMetricUpwindCentral2";
        }throw std::runtime_error("Unknown scalar implementation");
    }
    std::string_view runtimeProvider() const {
        switch(implementation_) {
        case ScalarImplementation::CapacityTransport:return ScalarOps::CapacityProvider;
        case ScalarImplementation::SharedTransport:return ScalarOps::SharedFluxProvider;
        case ScalarImplementation::ConservativeTransport:return ScalarOps::TransportProvider;
        default:return ScalarOps::Provider;
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
            const NumericalBinding& binding) const {
        using M=ScalarImplementation;
        const bool arithmeticFaces=implementation_==M::ArithmeticDiffusion
            || ((implementation_==M::FixedPoint || implementation_==M::FixedPointVolume)
                && std::any_of(binding.operators.begin(),binding.operators.end(),[](const auto& op){return op.provider=="scalar.central2.arithmetic-faces";}));
        const bool volumeFlux=implementation_==M::VolumeTransport || implementation_==M::FixedPointVolume;
        const bool metric=implementation_==M::MetricTransport;
        const bool fixedPoint=implementation_==M::FixedPoint || implementation_==M::FixedPointVolume;
        const bool velocityTransport=implementation_==M::VelocityTransport || volumeFlux || metric;
        const bool primitive=implementation_==M::CapacityTransport;
        const bool shared=primitive || implementation_==M::SharedTransport;
        const bool transport=shared || implementation_==M::ConservativeTransport;
        if(fixedPoint && (!call.stepScoped || call.stageScoped))throw std::runtime_error("Fixed-point equation needs explicit Step scope: "+call.occurrence);
        if ((!volumeFlux && !binding.inputs.empty()) || (!fixedPoint && !metric && !binding.parameters.empty()))
            throw std::runtime_error(std::string(id())+": unsupported extra inputs/parameters; coefficients belong to WHAT, method options must be explicitly supported.");
        const auto& state=system.state.at(call.target.symbol);
        if (call.target.kind!=TargetKind::Physical || (state.shape!=ValueShape::Scalar && (transport || state.shape!=ValueShape::Vector || state.components!=3))
            || (state.shape==ValueShape::Scalar && state.components!=1) || state.componentOffset!=0
            || state.location!=VariableLocation::EulerianCell
            || (state.storageBinding!=StorageBinding::NamedDistributed && (!(primitive || fixedPoint) || state.storageBinding!=StorageBinding::ProviderDistributed))
            || state.constantValue || state.readOnly())
            throw std::runtime_error("Unsupported scalar target: requires writable cell named storage, or a supported borrowed scalar owner for fixed/capacity transport.");
        const auto& equation=system.registry.at(call.equation);
        using K=FormulaExpr::Kind;
        const auto target=state.expression();
        auto concentration=target;
        if(transport) {
            const auto density=BuiltinStateCatalog().expression("rho");
            concentration=primitive?target:FormulaExpr::divide(target,density);
            const auto transient=FormulaExpr::op(Math::Op::Ddt,{primitive?FormulaExpr::multiply(density,target):target});
            const auto lhs=FormulaExpr::add(transient,FormulaExpr::op(Math::Op::Div,{BuiltinStateCatalog().expression("rhoU"),concentration}));
            if(!sameFormulaExpression(equation.lhs,lhs))
                throw std::runtime_error("Selected conservative transport implementation requires its declared capacity and mass-flux LHS; generic explicit mass inversion is not implemented.");
        }
        double transportD=0;
        CompiledScalarEquation contract;contract.components=state.components;contract.arithmeticFaces=arithmeticFaces;contract.target=call.target.symbol;contract.sequential=call.stepScoped;contract.fixedTime=fixedPoint;contract.conservativeTransport=transport;contract.sharedMassFlux=shared;contract.primitiveStorage=primitive;
        if(metric) {
            if(contract.components!=1 || call.stepScoped)throw std::runtime_error("Metric scalar requires a scalar recipe participant; Step/fixed-point treatment is not implemented: "+call.occurrence);
            contract.metric=metricOptions(binding);
        }
        if(fixedPoint) {
            if(contract.components!=1)throw std::runtime_error("Fixed-point ordinary implementation currently requires scalar target");
            contract.relaxation=fixedRelaxation(binding,call.occurrence);
        }
        if(volumeFlux) {
            if(!call.stepScoped || binding.inputs.size()!=1)throw std::runtime_error("Volume-flux transport requires Step and one explicitly selected face STATE input");
            const auto& flux=system.state.at(binding.inputs.front());
            Math::requireCompatible(flux.signature,Math::signature(0,{0,3,-1,0,0,0,0}),"volume face flux");
            if(flux.location!=VariableLocation::EulerianFace || flux.components!=1 || flux.ownership!=OwnershipKind::CanonicalFace)
                throw std::runtime_error("Volume flux requires scalar canonical face backing");
            contract.volumeFlux=flux.id;
            if(std::count_if(call.inputs.begin(),call.inputs.end(),[&](const auto& use){return use.symbol==flux.id && use.version==StateVersion::Current;})!=1)throw std::runtime_error("Volume flux requires explicit Current STATE input");
        }
        contract.valueReads={{contract.target,call.stepScoped?StateVersion::Current:StateVersion::Stage}};
        if (transport) {
            const auto& rho=system.state.at("rho");const auto& momentum=system.state.at("rhoU");
            if (rho.components!=1 || momentum.components!=3 || rho.storageBinding!=StorageBinding::PackedDistributed
                || momentum.storageBinding!=StorageBinding::PackedDistributed)
                throw std::runtime_error("Scalar transport requires packed Stage density/momentum.");
            contract.valueReads.push_back({"rho",StateVersion::Stage});contract.valueReads.push_back({"rhoU",StateVersion::Stage});
        }
        const auto appendRead=[&](const std::string& name,StateVersion version,TargetKind kind=TargetKind::Physical) {
            const auto found=std::find_if(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==name;});
            if(found!=contract.valueReads.end()) {
                if(name==contract.target && version==StateVersion::OldTime && !transport) {found->version=version;return;}
                if(found->version!=version || found->kind!=kind)throw std::runtime_error("Conflicting value read versions: "+name);
            } else {
                StateUse use{name,version};use.kind=kind;
                const auto authored=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& input){return input.symbol==name;});
                if(authored!=call.inputs.end()){use.halo=authored->halo;use.snapshot=authored->snapshot;use.perIteration=authored->perIteration;}
                contract.valueReads.push_back(std::move(use));
            }
        };
        const auto readSymbol=[&](const auto& self,const std::string& name,std::optional<StateVersion> inherited=std::nullopt)->void {
            const auto& partner=system.state.at(name);
            if(volumeFlux && name==contract.volumeFlux)return;
            const auto explicitRead=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==name;});
            const auto kind=explicitRead==call.inputs.end()?TargetKind::Physical:explicitRead->kind;
            if(kind!=TargetKind::Physical && (kind!=TargetKind::Working || !call.stepScoped || transport || name==contract.target || partner.producer))
                throw std::runtime_error("Unsupported Working input treatment: step="+call.occurrence+" input="+name+"; requires ordinary Step Current backing independent of target");
            const bool evolved=std::find(system.state.solutionVariables().begin(),system.state.solutionVariables().end(),name)!=system.state.solutionVariables().end();
            if(!transport && !inherited && !evolved && explicitRead==call.inputs.end())throw std::runtime_error("Missing explicit input version: "+name);
            const auto version=inherited.value_or(explicitRead==call.inputs.end()?StateVersion::Stage:explicitRead->version);
            if(inherited && explicitRead!=call.inputs.end() && explicitRead->version!=version)throw std::runtime_error("Producer dependency version conflict: "+name);
            if(version!=StateVersion::Current && version!=StateVersion::Stage && version!=StateVersion::OldTime && !(fixedPoint && version==StateVersion::Frozen && explicitRead!=call.inputs.end() && explicitRead->perIteration))throw std::runtime_error("Unsupported input version: "+name+" / "+toString(version)+" requires a real snapshot binding");
            if(call.stepScoped && version!=StateVersion::Current && !(fixedPoint && (version==StateVersion::OldTime || version==StateVersion::Frozen)))throw std::runtime_error("Sequential Step input requires Current backing: "+name);
            if(version==StateVersion::OldTime && (transport || primitive || !evolved || partner.producer))throw std::runtime_error("Unsupported OldTime input: requires an existing ordinary temporal snapshot: "+name);
            // The complete program validates actual Stage publishers; membership in
            // the evolved unknown set does not decide whether a coefficient is readable.
            if((partner.shape!=ValueShape::Scalar && (transport || partner.shape!=ValueShape::Vector || partner.components!=3)) || partner.componentOffset!=0
                || partner.location!=VariableLocation::EulerianCell || partner.constantValue)
                throw std::runtime_error("Unsupported scalar input shape/storage: "+name);
            if(partner.producer) {
                if(!(arithmeticFaces || metric) || partner.evaluation!=StateEvaluation::Lazy)throw std::runtime_error("Unsupported coefficient producer: explicit materialized HOW equation/publication required: "+name);
                const bool visited=std::any_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==name;});
                appendRead(name,version,kind);if(visited)return;
                for(const auto& dep:partner.dependencies)self(self,dep,version);
                return;
            }
            if(partner.storageBinding!=StorageBinding::NamedDistributed && (!(primitive || arithmeticFaces || velocityTransport || fixedPoint) || (partner.storageBinding!=StorageBinding::ProviderDistributed && partner.storageBinding!=StorageBinding::TransientWorkspace && partner.storageBinding!=StorageBinding::PackedDistributed)))
                throw std::runtime_error("Unsupported scalar input shape/storage: "+name);
            appendRead(name,version,kind);
        };
        const auto stageSymbol=[&](const std::string& name){readSymbol(readSymbol,name);};
        const auto rhs=transport?equation.rhs:explicitFormulaRhs(equation,contract.target);
        const auto structure=inspectFormula(rhs,contract.target);
        int diffusion=0;std::vector<std::string> velocities;
        for(const auto& item:structure.nodes) {
            const auto& node=*item.expression;
            if(node.kind()==K::Operator && node.operation()==Math::Op::Div && velocityTransport) {
                if(node.arguments().size()!=2 || node.arguments()[0].kind()!=K::Symbol
                    || node.arguments()[1].kind()!=K::Symbol || node.arguments()[1].name()!=contract.target
                    || system.state.at(node.arguments()[0].name()).components!=3 || contract.components!=1)
                    throw std::runtime_error("Upwind1 requires div(registered velocity vector, target scalar); no compound flux lowering");
                const auto& velocity=node.arguments()[0];
                Math::requireCompatible(velocity.signature(),Math::signature(1,{0,1,-1,0,0,0,0}),"Upwind1 velocity");
                if(volumeFlux)contract.fluxVelocity=velocity.name();
                stageSymbol(velocity.name());velocities.push_back(velocity.name());continue;
            }
            if (node.kind()==K::Operator && node.operation()==Math::Op::Diffusion) {
                if(node.arguments().size()!=2 || !(transport?sameFormulaExpression(node.arguments()[1],concentration)
                    :node.arguments()[1].kind()==K::Symbol))
                    throw std::runtime_error("Unsupported: Central2 diffusion requires a declared Stage symbol/concentration operand");
                if(metric && node.arguments()[1].name()!=contract.target)throw std::runtime_error("Metric diffusion currently requires its own scalar target; cross-field spatial preparation is not implemented: "+call.occurrence);
                if (++diffusion!=1 && transport) throw std::runtime_error("Scalar transport provider requires one diffusion term.");
                if(!transport)stageSymbol(node.arguments()[1].name());
                std::string coefficientState;
                const FormulaExpr* coefficientPtr=&node.arguments()[0];
                if (transport) {
                    const auto& coefficient=*coefficientPtr;
                    if (coefficient.kind()!=K::Multiply || coefficient.arguments().size()!=2
                        || !sameFormulaExpression(coefficient.arguments()[0],BuiltinStateCatalog().expression("rho")))
                        throw std::runtime_error("Scalar transport diffusion coefficient requires rho*constantD.");
                    coefficientPtr=&coefficient.arguments()[1];
                }
                const auto& coefficient=*coefficientPtr;
                if (coefficient.kind()==K::Constant) transportD=coefficient.constant();
                else if (coefficient.kind()==K::Symbol && system.state.at(coefficient.name()).constantValue)
                    transportD=*system.state.at(coefficient.name()).constantValue;
                else if(arithmeticFaces && coefficient.kind()==K::Symbol) {
                    const auto occurrence=node.occurrence().empty()?item.path:node.occurrence();
                    if(fixedPoint && std::none_of(binding.operators.begin(),binding.operators.end(),[&](const auto& op) {
                        return (op.formula.empty() || op.formula==equation.id) && op.occurrence==occurrence
                            && op.provider=="scalar.central2.arithmetic-faces";
                    }))throw std::runtime_error("Variable diffusion needs explicit WHICH arithmetic face interpolation: equation="+equation.id+" occurrence="+occurrence);
                    stageSymbol(coefficient.name());coefficientState=coefficient.name();
                    const auto& field=system.state.at(coefficient.name());
                    if(field.storageBinding==StorageBinding::TransientWorkspace && field.workspaceBoundary!=WorkspaceBoundary::NearestInterior)
                        throw std::runtime_error("Coefficient lacks published spatial read region: "+field.id);
                }
                else if(shared && coefficient.kind()==K::Symbol) {
                    stageSymbol(coefficient.name());contract.diffusivityState=coefficient.name();
                    const auto& field=system.state.at(coefficient.name());
                    if(field.storageBinding==StorageBinding::TransientWorkspace && field.workspaceBoundary!=WorkspaceBoundary::NearestInterior)
                        throw std::runtime_error("Missing coefficient read-region preparation: step="+call.occurrence+" equation="+call.equation+" target="+call.target.symbol+" input="+field.id);
                } else throw std::runtime_error("Unsupported scalar diffusivity: select a supported constant, explicit arithmetic-face occurrence, or shared mass-flux coefficient.");
                if (contract.diffusivityState.empty() && (!std::isfinite(transportD) || (transport && (shared?transportD<0:transportD<=0))))
                    throw std::runtime_error("Diffusion coefficient must be finite; transport requires positive/nonnegative diffusivity.");
                contract.diffusions.push_back({node.occurrence().empty()?item.path:node.occurrence(),
                    transport?contract.target:node.arguments()[1].name(),transportD,node.arguments()[1].signature().rank==0?1:3,coefficientState});
                continue;
            }
            if (node.kind()==K::Operator && node.operation()==Math::Op::Source) {
                if(contract.components!=1)throw std::runtime_error("Unsupported: vector known-source callbacks need a component-aware source port.");
                if(node.arguments().size()!=1)throw std::runtime_error("Source requires one complete local expression");
                if(node.arguments()[0].kind()!=K::Symbol || system.state.contains(node.arguments()[0].name()))continue;
                if(std::find(contract.sources.begin(),contract.sources.end(),node.arguments()[0].name())==contract.sources.end())contract.sources.push_back(node.arguments()[0].name());
            }
        }
        for(const auto& name:structure.symbols) {
            if(std::any_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==name;}))continue;
            if(std::find(contract.sources.begin(),contract.sources.end(),name)!=contract.sources.end())continue;
            if(system.state.contains(name) && system.state.at(name).constantValue)continue;
            stageSymbol(name);
        }
        if (transport && diffusion!=1) throw std::runtime_error("Scalar transport provider requires one diffusion term.");
        if(arithmeticFaces)for(const auto& d:contract.diffusions) {
            const auto input=std::find_if(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==d.input;});
            if(input==contract.valueReads.end() || input->version!=(call.stepScoped?StateVersion::Current:StateVersion::Stage))
                throw std::runtime_error("Unsupported diffusion operand version: selected Current (Step) or Stage publisher required: "+d.input);
        }
        for(const auto& read:call.inputs) {
            if((read.kind!=TargetKind::Physical && !(read.kind==TargetKind::Working && call.stepScoped && read.version==StateVersion::Current && !transport))
                || read.halo || !read.snapshot.empty() || (read.perIteration && !(fixedPoint && read.version==StateVersion::Frozen)))throw std::runtime_error("Unsupported authored scalar read qualifier: "+read.symbol);
            if(std::count_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==read.symbol;})!=1)throw std::runtime_error("Duplicate input read declaration: "+read.symbol);
            if(!(volumeFlux && read.symbol==contract.volumeFlux && read.version==StateVersion::Current) && std::none_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==read.symbol && r.version==read.version && r.kind==read.kind;}))throw std::runtime_error("Unconsumed input version declaration: "+read.symbol);
        }
        std::vector<std::string> slots;for(const auto& read:contract.valueReads)slots.push_back(read.symbol);
        FormulaOperatorCatalog operators;
        if(velocityTransport) {
            FormulaOperatorProvider advection;advection.id=metric?"scalar.metric-upwind1":volumeFlux?"scalar.volume-upwind1":"scalar.velocity-upwind1";advection.mathematicalOperator=Math::Op::Div;
            advection.matches=[velocities,target=contract.target](const auto& e){return e.arguments().size()==2 && e.arguments()[0].kind()==K::Symbol && e.arguments()[1].kind()==K::Symbol && e.arguments()[1].name()==target && std::find(velocities.begin(),velocities.end(),e.arguments()[0].name())!=velocities.end();};
            advection.compileValue=[slots,target=contract.target,volume=volumeFlux,metric=metric](const auto& e){
                const auto velocity=std::find(slots.begin(),slots.end(),e.arguments()[0].name())-slots.begin();
                const auto scalar=std::find(slots.begin(),slots.end(),target)-slots.begin();
                return metric?Discretization::scalarMetricUpwind1(velocity,scalar):volume?Discretization::scalarVolumeUpwind1(scalar):Discretization::scalarVelocityUpwind1(velocity,scalar);
            };operators.add(std::move(advection));
        }
        if(metric) {
            FormulaOperatorProvider diffusion;diffusion.id="scalar.metric-central2";diffusion.mathematicalOperator=Math::Op::Diffusion;
            diffusion.matches=[](const auto& e){return e.arguments().size()==2 && e.arguments()[1].kind()==K::Symbol;};
            diffusion.compileValue=[slots,&system](const auto& e) {
                const auto& coefficient=e.arguments()[0];
                const double d=coefficient.kind()==K::Constant?coefficient.constant():*system.state.at(coefficient.name()).constantValue;
                if(!std::isfinite(d) || d<0)throw std::runtime_error("Metric diffusion requires a finite nonnegative constant coefficient");
                const auto slot=std::find(slots.begin(),slots.end(),e.arguments()[1].name())-slots.begin();
                return FormulaValueKernel([slot,d](const auto& values,int cell,int){return values.metricDiffusion(slot,cell,d);});
            };operators.add(std::move(diffusion));
        } else if(shared) {
            auto slot=std::find(slots.begin(),slots.end(),contract.diffusivityState)-slots.begin();
            operators.add(Discretization::conservativeScalarCentral2(transportD,contract.diffusivityState.empty()?-1:static_cast<int>(slot),primitive));
        } else if(transport)operators.add(Discretization::conservativeScalarCentral2(transportD));
        else operators.add(Discretization::central2ScalarDiffusion(system.state,slots,arithmeticFaces));
        FormulaOperatorProvider source;source.id="scalar.known-source";source.mathematicalOperator=Math::Op::Source;
        source.matches=[](const auto&) {return true;};
        source.compileValue=[offset=slots.size(),names=contract.sources,slots,&equation,&system](const FormulaExpr& expression) -> FormulaValueKernel {
            const auto& argument=expression.arguments()[0];
            if(argument.kind()!=K::Symbol || std::find(names.begin(),names.end(),argument.name())==names.end()) {
                std::vector<std::string> used;
                return compileFormulaValue(equation,argument,{},used,slots,&system.state);
            }
            const auto slot=offset+static_cast<std::size_t>(std::find(names.begin(),names.end(),argument.name())-names.begin());
            return [slot](const FormulaValues& values,int cell,int) {
                return values.boundReads.at(slot)(cell,0);
            };
        };
        operators.add(std::move(source));
        CompiledEquationCall result;
        contract.nonlinear=structure.targetDegree>1;
        contract.normalized=formulaText(transport?equation:Equation{equation.id,FormulaExpr::op(Math::Op::Ddt,{target}),rhs,equation.origin});
        contract.rhs=compileFormulaValue(equation,rhs,operators,result.operatorBindings,slots,&system.state,&contract.lowering,binding.operators);
        if(transport)result.operatorBindings.push_back(shared?"scalar.canonical-mass-flux.upwind1":"scalar.independent-center-rhoU.upwind1");
        if(shared)result.stageRequires={"canonicalMassFlux(Stage)"};
        for(const auto& read:contract.valueReads) {
            const bool neighbor=(!volumeFlux && std::find(velocities.begin(),velocities.end(),read.symbol)!=velocities.end())
                || std::any_of(contract.diffusions.begin(),contract.diffusions.end(),[&](const auto& d){return d.input==read.symbol || d.coefficientState==read.symbol;});
            if(neighbor)contract.stencilReads.push_back(read);
        }
        if(fixedPoint)for(const auto& read:contract.stencilReads)if(read.version==StateVersion::OldTime)
            throw std::runtime_error("Fixed-time step="+call.occurrence+" input="+read.symbol+": OldTime is supported for local reads, not endpoint stencil preparation");
        const auto valueReads=contract.valueReads;
        {
            std::vector<std::string> boundaries;
            if(diffusion || transport || !velocities.empty())boundaries.push_back(call.target.symbol);
            for(const auto& d:contract.diffusions) {
                boundaries.push_back(d.input);
                if(!d.coefficientState.empty() && system.state.at(d.coefficientState).storageBinding!=StorageBinding::TransientWorkspace)
                    boundaries.push_back(d.coefficientState);
            }
            if(metric)for(const auto& velocity:velocities)boundaries.push_back(velocity);
            auto resource=equationFieldResource(call.occurrence,call.target.symbol,contract.valueReads,std::move(boundaries),contract.sources);
            if(metric)resource.fields->setBoundaryOptions=contract.metric;
            if(state.storageBinding==StorageBinding::ProviderDistributed)resource.owner="model.scalar";
            result.physicalResources.push_back(std::move(resource));
        }
        result.providerContract=std::move(contract);
        result.source=call;result.target.symbol=call.target.symbol;result.target.kind=call.target.kind;
        result.equationMethod=id();result.temporalResidual=true;result.synchronousStages=true;
        result.temporalStepSize=ScalarOps::Dt;
        result.temporalSnapshot=instanceOperation(ScalarOps::Begin,call.occurrence);
        result.temporalStagePrepare=instanceOperation("scalar.stage.prepare",call.occurrence);
        result.temporalRhs=instanceOperation("scalar.stage.rhs",call.occurrence);
        result.temporalAdvance=instanceOperation("scalar.stage.advance",call.occurrence);
        result.temporalPublish=instanceOperation(ScalarOps::Commit,call.occurrence);
        result.backendOperation=ScalarOps::Stage;
        result.oldTimeWorkspace="scalar."+call.target.symbol+".old";
        result.stageWorkspace="scalar."+call.target.symbol+".stage";
        result.residualWorkspace="scalar."+call.target.symbol+".rhs";
        result.requiredHaloWidth=(metric || diffusion || !velocities.empty())?1:0;
        result.numericalWorkspaces={result.oldTimeWorkspace,result.stageWorkspace,result.residualWorkspace};
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4};
        if(!transport)result.temporalCapabilities.push_back(FDM::TimeRecipeId::SSPRK3);
        if(metric)result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        result.calls={{call.equation,call.target.symbol}};
        for(const auto& read:valueReads)result.reads.push_back(read.symbol);result.writes={call.target.symbol};
        result.stateUses=valueReads;
        if(volumeFlux) {result.stateUses.push_back({binding.inputs.front(),StateVersion::Current});result.workspaceRequires.push_back(faceFluxPublication(system.state.at(binding.inputs.front()).storageKey));}
        result.stateEffects={{result.residualWorkspace,true,false,false}};
        result.requirements={"serial single patch","common group dt","component-major scalar/vector workspace","stage-time known source"};
        if(diffusion)result.requirements.push_back(arithmeticFaces
            ? "Cartesian Central2 arithmetic faces; existing nodal fixedValue/zeroGradient ghost; coefficient neighborhood reads mirrored physical nodes in the selected version"
            : "Cartesian Central2; each stencil input requires its own fixedValue/zeroGradient Stage boundary");
        if(metric)result.requirements={"serial single patch; ForwardEuler; metric Upwind1/Central2; explicit set-boundary/ILW/bounds; no EOS feedback", "velocity reads its declared version; no packed-state identity requirement"};
        if(primitive)result.requirements.push_back("primitive physical alias; rho-weighted conservative stage workspace; declared coefficient publishers");
        if(transport)result.requirements.push_back("conservative Upwind1 advection and rho-weighted Central2 concentration diffusion; boundary values refer to the declared physical scalar target");
        for (const auto& op:{result.temporalSnapshot,result.temporalStagePrepare,result.temporalRhs,
                result.temporalAdvance,result.temporalPublish})
            result.operations.push_back({op,op,op==result.temporalPublish?OperationStage::StepCommit:OperationStage::Prepare,
                {OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        result.backendProvider=runtimeProvider();
        if(call.stepScoped) {
            if(transport)throw std::runtime_error("Sequential capacity/packed transport has no implemented Step treatment");
            result.temporalResidual=false;result.synchronousStages=false;
            result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
            result.temporalStepSize.clear();result.temporalSnapshot.clear();result.temporalStagePrepare.clear();
            result.temporalRhs.clear();result.temporalAdvance.clear();result.temporalPublish.clear();
            result.backendOperation.clear();result.operations.clear();
            result.fragment.kind=PlanNodeKind::Sequence;result.fragment.id="scalar.interval."+call.occurrence;
            for(const auto* phase:{"snapshot","prepare","rhs","advance","accept"}) {
                const auto operation=instanceOperation(std::string("scalar.step.")+phase,call.occurrence);
                SolvePlanNode node;node.kind=PlanNodeKind::Update;node.id=node.name=node.operation=operation;
                node.provider=result.backendProvider;node.equationCalls=result.calls;result.fragment.children.push_back(std::move(node));
                result.operations.push_back({operation,phase,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
            }
            result.stateEffects={{call.target.symbol,true,false,true}};
            result.requirements.push_back("explicit Step: ForwardEuler only; consume host dt; once outside repeated scopes; no clock ownership");
        }
        if(fixedPoint) {
            result.temporalCapabilities={FDM::TimeRecipeId::PartitionedFixedEuler};
            result.stepInputs={{OpIds::CellPrepare,"equation.cell-explicit",true}};
            result.physicalPublication={OpIds::CellAccept,"equation.cell-explicit"};
            result.backendOperation=instanceOperation("scalar.iterate",call.occurrence);result.fragment={};result.operations.clear();
            result.residualCapabilities={ResidualKind::StateChange,ResidualKind::DiscreteEquation};
            result.requirements.back()="fixed physical interval: old base, endpoint RHS and explicit relaxation; repetition and stopping belong to HOW";
            CompiledStateView old;old.symbol=state.id;old.kind=StateViewKind::OldTime;old.storage="scalar."+call.occurrence+".old";
            old.components=state.components;old.owner=StateViewOwner::NumericalProvider;result.viewOutputs.push_back(old);result.oldTimeWorkspace=old.storage;
            for(const auto& operation:{result.backendOperation,std::string(OpIds::CellPrepare),std::string(OpIds::CellAccept)})
                result.operations.push_back({operation,operation,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        }
        return result;
    }
};
} // namespace
std::optional<CompiledEquationCall> compileScalarEquation(const ExecutableEquationSystem& system,const EquationCall& call,const NumericalBinding& binding) {
    using M=ScalarImplementation;
    static const ScalarMethod implementations[]={
        ScalarMethod(M::Diffusion),ScalarMethod(M::ArithmeticDiffusion),
        ScalarMethod(M::VelocityTransport),ScalarMethod(M::VolumeTransport),
        ScalarMethod(M::ConservativeTransport),ScalarMethod(M::SharedTransport),
        ScalarMethod(M::CapacityTransport),ScalarMethod(M::FixedPoint),
        ScalarMethod(M::FixedPointVolume),ScalarMethod(M::MetricTransport)};
    if(binding.method=="ThermalEulerTransport")throw std::runtime_error("ThermalEulerTransport was retired; select ScalarMetricUpwindCentral2 and declare all mathematics in WHAT");
    for(const auto& implementation:implementations)
        if(implementation.id()==binding.method)return implementation.compile(system,call,binding);
    return std::nullopt;
}

} // namespace SF::System
