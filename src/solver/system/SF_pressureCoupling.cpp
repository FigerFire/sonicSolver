#include "core/system/SF_operationIds.h"
/// @file SF_pressureCoupling.cpp
/// @brief pressure-constraint coupling preset 的匹配与贡献实现。

#include "SF_pressureCoupling.h"
#include "SF_equationContribution.h"
#include "SF_affine.h"
#include "SF_cell.h"
#include "SF_scalarMethod.h"
#include "SF_solvePlan.h"
#include "SF_numericalSystem.h"

#include <algorithm>
#include <stdexcept>
#include <set>

namespace SF::System {
namespace {

/// @brief 约束声明的 pressure 未知量是否存在（role 不限：Derived closure 与
///        Multiplier 都是合法压力表示）。
const ConstraintDescriptor* pressureConstraintOf(
        const RawEquationSystem& raw,const PressureCouplingInputs* inputs=nullptr) {
    for (const ConstraintDescriptor& constraint : raw.constraints) {
        if (inputs ? constraint.id==inputs->constraint : (constraint.id == "C_INCOMPRESSIBILITY"
            || constraint.id == "C_SHARED_PRESSURE")) {
            return &constraint;
        }
    }
    return nullptr;
}

bool rawHasUnknown(const RawEquationSystem& raw, const std::string& id) {
    return std::any_of(
        raw.state.symbols().begin(),raw.state.symbols().end(),
        [&](const StateSymbol& unknown) { return unknown.id == id; });
}

bool rawHasEquation(const RawEquationSystem& raw,const std::string& id) {
    return raw.registry.contains(id);
}

bool hasMomentumEquation(const RawEquationSystem& raw) {
    return std::any_of(raw.registry.entries().begin(),raw.registry.entries().end(),
        [&](const Equation& equation) {
            const auto vectorDerivative=[&](const auto& self,const FormulaExpr& expression)->bool {
                if (expression.kind()==FormulaExpr::Kind::Operator && expression.operation()==Math::Op::Ddt
                    && expression.arguments().size()==1 && expression.arguments()[0].kind()==FormulaExpr::Kind::Symbol) {
                    const auto& id=expression.arguments()[0].name();
                    return raw.state.contains(id) && raw.state.at(id).components>1;
                }
                return std::any_of(expression.arguments().begin(),expression.arguments().end(),[&](const auto& child){return self(self,child);});
            };
            return vectorDerivative(vectorDerivative,equation.lhs);
        });
}

} // namespace

std::string validatePressureProgram(const ExecutableEquationSystem& equations,
        const CompiledNumericalSystem& numerics,const CompiledSolvePlan& plan) {
    const auto& program=plan.compiledProgram;
    if(!numerics.facePrediction || numerics.facePrediction->response.empty() || !numerics.pressureOperator
        || !numerics.facePrediction)
        return "Pressure correction needs a compiled cell response, numerical inputs and persistent Rhie-Chow flux";
    if(equations.constraints.size()!=1 || equations.constraints.front().id!="C_INCOMPRESSIBILITY")
        return "Constant-density pressure requires its incompressibility constraint; additional constraint treatment is not implemented";
    if(plan.compiledProgram.hasTemporalRoot || !program.temporalParticipants.empty())
        return "Pressure correction has no coupled RK participant/stage treatment";
    for(const auto& boundary:equations.boundaryClosures)
        return "Pressure spatial boundary port is not implemented for declared closure "+boundary.id+" (owner="+boundary.provider+")";
    bool fixed=false;
    const auto response=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& call){return call.source.occurrence==numerics.facePrediction->response;});
    const bool working=response!=program.steps.end() && response->target.kind==TargetKind::Working;
    std::set<std::string> protectedStorage;
    const auto protect=[&](const std::string& id) {
        protectedStorage.insert(id);
        if(equations.state.contains(id))protectedStorage.insert(equations.state.at(id).storageKey);
    };
    for(const auto& call:program.steps)if(!cellEquation(call) && !std::any_cast<CompiledAffineEquation>(&call.providerContract) && !std::any_cast<CompiledScalarEquation>(&call.providerContract) && call.coefficientProvides.empty()) {
        const auto context="Pressure HOW step="+call.source.occurrence+" equation="+call.source.equation+" target="+call.target.symbol+": ";
        const std::set<std::string> pressureActions{
            OpIds::CellPrepare,OpIds::CellAccept,OpIds::FaceFluxPredict,OpIds::PressureBoundaryPrepare,
            OpIds::PressureAssemble,OpIds::PressureSolve,OpIds::PressureUpdatePrepare,
            OpIds::VelocityCorrect,OpIds::FluxCorrect,OpIds::PressureCorrectionCommit,
            OpIds::PressureRelaxationApply,OpIds::PressureFluxConsistencyRestore,OpIds::PressureConvergenceEvaluate};
        bool found=false,valid=true;
        const auto action=[&](const std::string& id){if(!id.empty()){found=true;valid=valid && pressureActions.count(id);}};
        action(call.backendOperation);
        const auto fragment=[&](const auto& self,const SolvePlanNode& node)->void {
            action(node.operation);for(const auto& child:node.children)self(self,child);
        };fragment(fragment,call.fragment);
        if(!found || !valid)return context+"selected operations lack an implemented fixed-time pressure resource treatment";
        fixed|=std::any_of(call.stepInputs.begin(),call.stepInputs.end(),[](const auto& input) {
            return input.operation==OpIds::PressureStepBegin;
        });
        for(const auto* ids:{&call.writes,&call.workspaceRequires,&call.workspaceProvides,&call.workspaceInvalidates})
            for(const auto& id:*ids)protect(id);
        protect(call.target.workspace);
        for(const auto& resource:call.target.resources)protect(resource.storage);
    }
    for(const auto* id:{"U","p","rho","phi","deltaP"}) {
        if(!equations.state.contains(id))return "Pressure resource missing STATE "+std::string(id);
        protect(id);
    }
    for(const auto& layout:std::vector<std::pair<std::string,int>>{{"U",3},{"p",1},{"phi",1}}) {
        const auto& state=equations.state.at(layout.first);
        if(state.components!=layout.second || state.readOnly())
            return "Incompatible writable pressure STATE layout: "+state.id;
    }
    if(working!=fixed)return "Pressure Working predictor and fixed-time relaxation must be selected together";
    std::set<std::string> localOwners;
    for(const auto& call:program.steps)if(std::any_cast<CompiledAffineEquation>(&call.providerContract)) {
        const auto context="Pressure local step="+call.source.occurrence+" equation="+call.source.equation+" target="+call.target.symbol+": ";
        const auto& output=equations.state.at(call.target.symbol);
        if(protectedStorage.count(output.id) || protectedStorage.count(output.storageKey))
            return context+"aliases protected predictor/pressure storage="+output.storageKey+"; no explicit consistency treatment for this writer";
        const bool independentBacking=output.storageBinding==StorageBinding::TransientWorkspace
            || (output.evaluation==StateEvaluation::Materialized
                && (output.storageBinding==StorageBinding::NamedDistributed
                    || output.storageBinding==StorageBinding::ProviderDistributed));
        if(call.stageLocal || call.temporalResidual || call.target.kind!=TargetKind::Physical
            || !independentBacking || output.storageKey.empty())
            return context+"requires independent materialized scalar storage; pressure state/flux writes need an implemented consistency treatment";
        if(!localOwners.insert(output.storageKey).second)
            return context+"storage has multiple local publication owners";
        for(const auto& read:call.stateUses)
            if(read.version!=StateVersion::Current || (read.kind!=TargetKind::Physical && read.kind!=TargetKind::Working) || read.halo || !read.snapshot.empty() || read.perIteration)
                return context+"requires Current Physical/Working input; no Stage/OldTime/snapshot pressure treatment is bound";
    }
    for(const auto& call:program.steps)if(!call.coefficientProvides.empty()) {
        const auto context="Pressure coefficient step="+call.source.occurrence+": ";
        if(call.temporalResidual || call.stageLocal || call.coefficientVersion!=StateVersion::Current)
            return context+"requires a nonintegrated Step publication; no pressure RK coefficient treatment";
        for(const auto& id:call.writes) {
            const auto& output=equations.state.at(id);
            if(protectedStorage.count(id) || protectedStorage.count(output.storageKey))return context+"aliases pressure/predictor storage="+output.storageKey;
        }
        for(const auto& use:call.stateUses)if(use.version!=StateVersion::Current || use.halo || use.perIteration || !use.snapshot.empty())
            return context+"requires actual Current Physical/Working inputs";
    }
    for(const auto& call:program.steps)if(const auto* scalar=std::any_cast<CompiledScalarEquation>(&call.providerContract)) {
        const auto& output=equations.state.at(scalar->target);
        if(!scalar->sequential || protectedStorage.count(output.id) || protectedStorage.count(output.storageKey)
            || (output.storageBinding!=StorageBinding::NamedDistributed && output.storageBinding!=StorageBinding::ProviderDistributed))
            return "Pressure scalar step="+call.source.occurrence+" requires independent owned scalar backing and explicit Step treatment";
    }
    const auto method=[&](const ExecutionScope& node)->std::string_view {
        const auto call=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& c){return c.source.occurrence==node.step.occurrence;});
        return call==program.steps.end()?std::string_view{}:call->equationMethod;
    };
    // The prediction consumes the published working sample and its original
    // response. Independent reads can intervene; writes cannot silently change
    // that tuple. Compare actual backing slices, including aliases.
    struct Slice {std::string storage;int offset,components;VariableLocation location;};
    const auto compiled=[&](const ExecutionScope& node)->const CompiledEquationCall* {
        if(node.kind!=ExecutionKind::EquationCall)return nullptr;
        const auto at=std::find_if(program.steps.begin(),program.steps.end(),[&](const auto& c){return c.source.occurrence==node.step.occurrence;});
        return at==program.steps.end()?nullptr:&*at;
    };
    const auto stateSlice=[&](const StateUse& use)->Slice {
        const auto& symbol=equations.state.at(use.symbol);
        if(use.kind==TargetKind::Physical)return {symbol.storageKey,symbol.componentOffset,symbol.components,symbol.location};
        const auto kind=use.kind==TargetKind::Working?StateViewKind::Working:StateViewKind::Correction;
        const auto view=std::find_if(program.stateViews.begin(),program.stateViews.end(),[&](const auto& v){return v.symbol==use.symbol && v.kind==kind;});
        if(view==program.stateViews.end())throw std::runtime_error("Missing prediction input view: "+use.symbol);
        return {view->storage,view->componentOffset,view->components,view->location};
    };
    std::string predictionError;
    const auto predictionInputs=[&](const auto& self,const ExecutionScope& scope)->void {
        std::set<std::size_t> consumed;
        for(std::size_t i=0;i<scope.children.size();++i) {
            const auto* predictor=compiled(scope.children[i]);
            if(!predictor || predictor->source.occurrence!=numerics.facePrediction->response)continue;
            std::size_t end=i+1;
            for(;end<scope.children.size();++end) {
                const auto* call=compiled(scope.children[end]);
                if(!call || call->source.occurrence==numerics.facePrediction->response)break;
                if(call->faceFlux && !call->faceFlux->corrected)break;
            }
            const auto* flux=end<scope.children.size()?compiled(scope.children[end]):nullptr;
            if(!flux || !flux->faceFlux || flux->faceFlux->corrected) {
                predictionError="Producer="+predictor->source.occurrence+" view=U.Working: missing surfaceFlux consumer in the same scope; no automatic snapshot/reordering";return;
            }
            if(flux->faceFlux->velocityView==TargetKind::Physical && predictor->target.kind!=TargetKind::Physical) {
                predictionError="Consumer="+flux->source.occurrence+" requests Physical velocity not published by producer="+predictor->source.occurrence;return;
            }
            consumed.insert(end);
            std::vector<Slice> guarded;
            for(const auto& use:flux->stateUses)guarded.push_back(stateSlice(use));
            guarded.push_back(stateSlice(StateUse{"U"}));
            std::set<std::string> workspaces(predictor->workspaceProvides.begin(),predictor->workspaceProvides.end());
            for(const auto& output:predictor->viewOutputs)workspaces.insert(output.storage);
            for(std::size_t j=i+1;j<end;++j) {
                const auto& call=*compiled(scope.children[j]);
                const auto fail=[&](const std::string& storage) {
                    predictionError="Producer="+predictor->source.occurrence+" consumer="+flux->source.occurrence
                        +" view=U.Working/p/response conflicts with writer="+call.source.occurrence
                        +" equation="+call.source.equation+" target="+call.target.symbol+" storage="+storage;
                };
                const auto writes=[&](const Slice& slice) {
                    if(workspaces.count(slice.storage)){fail(slice.storage);return;}
                    for(const auto& input:guarded)if(!slice.storage.empty() && slice.storage==input.storage && slice.location==input.location
                        && slice.offset<input.offset+input.components && input.offset<slice.offset+slice.components)fail(slice.storage);
                };
                for(const auto& resource:call.target.resources)if(resource.access!=ResourceAccessMode::Read)
                    writes({resource.storage,resource.componentOffset,resource.components,equations.state.contains(call.target.symbol)?equations.state.at(call.target.symbol).location:call.target.workspaceLocation});
                for(const auto& effect:call.stateEffects) {
                    if(equations.state.contains(effect.symbol))writes(stateSlice(StateUse{effect.symbol}));
                    else if(workspaces.count(effect.symbol))fail(effect.symbol);
                }
                for(const auto& output:call.viewOutputs)writes({output.storage,output.componentOffset,output.components,output.location});
                for(const auto* tokens:{&call.workspaceProvides,&call.workspaceInvalidates})for(const auto& token:*tokens)
                    if(workspaces.count(token))fail(token);
                if(!predictionError.empty())return;
            }
        }
        for(std::size_t i=0;i<scope.children.size();++i) {
            const auto* call=compiled(scope.children[i]);
            if(call && call->faceFlux && !call->faceFlux->corrected && !consumed.count(i)) {
                predictionError="Consumer="+call->source.occurrence+" requires one published predictor tuple in the same scope";return;
            }
            self(self,scope.children[i]);if(!predictionError.empty())return;
        }
    };predictionInputs(predictionInputs,program.root);
    if(!predictionError.empty())return predictionError;
    std::string error;
    const auto scopes=[&](const auto& self,const ExecutionScope& scope)->unsigned {
        if(scope.kind==ExecutionKind::StageLoop || scope.kind==ExecutionKind::BlockSolve || !scope.commitSources.empty())
            error="Pressure HOW scope="+scope.id+": stage/block or explicit state COPY has no pressure consistency treatment";
        // Local arithmetic does not turn an uncommitted pressure-solve loop
        // into a complete correction loop with velocity/pressure/flux updates.
        constexpr unsigned correction=1,other=2,predictor=4;
        const auto implementation=scope.kind==ExecutionKind::EquationCall?method(scope):std::string_view{};
        unsigned members=implementation=="PressureCorrection"?correction
            :scope.kind==ExecutionKind::EquationCall && scope.step.occurrence==numerics.facePrediction->response?other|predictor
            :implementation.empty() || implementation=="LocalAffine" || implementation=="CellExplicit"?0:other;
        for(const auto& child:scope.children)members|=self(self,child);
        if(scope.kind==ExecutionKind::Sequence && (members&predictor) && scope.repetitions!=1)
            error="Pressure HOW scope="+scope.id+": predictor repetition requires the explicit fixed-time iteration Loop contract";
        if(scope.kind==ExecutionKind::Loop && (members&predictor)) {
            if(!working && scope.repetitions!=1)
                error="Pressure HOW scope="+scope.id+": repeated physical predictor would advance a full dt; select the implemented fixed-time Working treatment";
            if(working && !scope.terminationSignal.empty() && scope.terminationSignal!=kPressureOuterConvergedSignal)
                error="Pressure HOW scope="+scope.id+": selected convergence kernel does not provide signal "+scope.terminationSignal;
        }
        if(members==correction && scope.repetitions!=1)
            error="Pressure HOW scope="+scope.id+": repeated uncommitted pressure solves require a non-orthogonal implementation";
        return members;
    };
    scopes(scopes,program.root);
    return error;
}

std::vector<Equation> conservativePressureRelations() {
    using Expr=FormulaExpr;
    const auto symbol=[](const char* id) { return PressureMathematics::symbol(id,true); };
    const Provenance origin{OriginKind::Generated,"conservative pressure correction"};
    // Corrector::assemble uses reaction + face-averaged inverse-density diffusion.
    const auto compliance=Expr::multiply(symbol("rho"),Expr::multiply(
        Expr::multiply(symbol("soundSpeed"),symbol("soundSpeed")),
        Expr::multiply(symbol("dt"),symbol("dt"))));
    std::vector<Equation> equations{
        {"pSimple",Expr::subtract(Expr::divide(symbol("deltaP"),compliance),
            Expr::op(Math::Op::Div,{Expr::multiply(
                PressureMathematics::op("faceAverage",{Expr::divide(Expr::constantValue(1.0),symbol("rho"))}),
                Expr::op(Math::Op::Grad,{symbol("deltaP")}))})),
            Expr::negate(Expr::divide(Expr::op(Math::Op::Div,{symbol("U")}),symbol("dt"))),origin},
        {"correctU",symbol("rhoU"),Expr::subtract(symbol("rhoU"),
            Expr::multiply(Expr::multiply(symbol("momentumRelaxation"),symbol("dt")),
                Expr::op(Math::Op::Grad,{symbol("deltaP")}))),origin},
        {"correctP",symbol("p"),Expr::add(symbol("p"),
            Expr::multiply(symbol("pressureRelaxation"),symbol("deltaP"))),origin},
        {"correctFluxp",symbol("fluxValidity"),
            PressureMathematics::op("invalidateDerivedFlux",{symbol("rhoU")}),origin},
        {"publishPressure",symbol("rhoE"),PressureMathematics::op("totalEnergyFromPressure",
            {symbol("rho"),symbol("rhoU"),symbol("p")}),origin},
        {"relaxIterate",symbol("iterate"),PressureMathematics::op("relax",
            {symbol("iterate"),symbol("laggedIterate")}),origin},
        {"restoreFlux",symbol("fluxValidity"),PressureMathematics::op("consistentFlux",
            {symbol("rhoU"),symbol("p")}),origin},
        {"checkConvergence",symbol("converged"),PressureMathematics::op("residualConvergence",
            {symbol("iterate"),symbol("laggedIterate")}),origin}};
    return equations;
}

const char* toString(CouplingStatus status) {
    switch (status) {
        case CouplingStatus::Active: return "active";
        case CouplingStatus::Inactive: return "inactive";
        case CouplingStatus::Invalid: return "invalid";
        case CouplingStatus::Unsupported: return "unsupported";
    }
    return "inactive";
}

CouplingPresetRequest couplingRequestFrom(
        const FDM::PressureCorrectionConfig& pressure,
        bool explicitlyRegistered) {
    const auto& coupling=pressure.coupling;
    CouplingPresetRequest request;
    request.iterationTolerance=pressure.relativeTolerance;request.continuityTolerance=pressure.absoluteTolerance;
    request.presetKind = coupling.preset;
    request.preset = FDM::toString(coupling.preset);
    request.outerCorrectors = coupling.outerCorrectors;
    request.pressureCorrectors = coupling.pressureCorrectors;
    request.nonOrthogonalCorrectors = coupling.nonOrthogonalCorrectors;
    request.explicitlyRegistered = explicitlyRegistered;
    request.origin = {OriginKind::BuiltinPreset,"coupling preset"};
    return request;
}

PressureCouplingInputs pressureCouplingInputs(EquationCall predictor,PressurePredictorKind representation) {
    const bool conservative=representation==PressurePredictorKind::ConservativeMomentum;
    PressureCouplingInputs input;
    input.predictor=std::move(predictor);input.representation=representation;
    input.predictorTarget={input.predictor.target.symbol,conservative ? TargetKind::Physical : TargetKind::Working};
    input.predictorOccurrence="predictor";input.constraint="C_INCOMPRESSIBILITY";
    input.correction={"pSimple",{"deltaP",TargetKind::Correction}};
    input.velocityUpdate={"correctU",{input.predictor.target.symbol}};
    input.pressureUpdate={"correctP",{"p",conservative ? TargetKind::Working : TargetKind::Physical}};
    input.fluxUpdate={"correctFluxp",{conservative ? "fluxValidity" : "phi",conservative ? TargetKind::Workspace : TargetKind::Physical}};
    if (!conservative) {
        input.fluxPrediction={"predictFlux",{"phi",TargetKind::Working}};
    }
    if (conservative) {
        input.fusedBalances={{"continuity",{"rho"}},{"energy",{"rhoE"}}};
        input.replacedCommit="physicalStep.commit";
        input.publication={"publishPressure",{"rhoE"}};
    }
    input.relaxation={"relaxIterate",{"iterate",TargetKind::Workspace}};
    input.fluxConsistency={"restoreFlux",{conservative ? "fluxValidity" : "phi",TargetKind::Workspace}};
    input.convergence={"checkConvergence",{"converged",TargetKind::Workspace}};
    return input;
}

void applyPressureExecution(ExecutionProgram& program,
                            const CouplingPresetRequest& request,
                            const PressureCouplingInputs& input) {
    auto& entries=program.root.children;
    const auto matches=[](const ExecutionScope& node,const EquationCall& call) {
        return node.kind==ExecutionKind::EquationCall && node.step.equation==call.equation
            && node.step.target.symbol==call.target.symbol && node.step.target.kind==call.target.kind
            && node.step.occurrence==call.occurrence;
    };
    if (input.predictorOccurrence.empty() || input.predictorTarget.symbol!=input.predictor.target.symbol)
        throw std::runtime_error("Pressure method requires an explicit predictor occurrence and matching target.");
    const auto count=std::count_if(entries.begin(),entries.end(),[&](const auto& node){return matches(node,input.predictor);});
    if (count!=1) throw std::runtime_error("Pressure HOW requires one referenced predictor: "+input.predictor.equation+" / "+input.predictor.occurrence);
    for (const auto& balance:input.fusedBalances)
        if (std::count_if(entries.begin(),entries.end(),[&](const auto& node){return matches(node,balance);})!=1)
            throw std::runtime_error("Pressure HOW requires one referenced fused balance: "+balance.equation);
    const auto found=std::find_if(entries.begin(),entries.end(),[&](const auto& node){return matches(node,input.predictor);});
    ExecutionScope predictor=*found;
    entries.erase(found);
    predictor.step.target=input.predictorTarget;
    // The single PISO predictor publishes U to the physical target. RhieChow can
    // still borrow the same occurrence's Working view for its face response.
    if(request.presetKind==FDM::PressureCouplingPreset::PISO
        && input.representation==PressurePredictorKind::Velocity)
        predictor.step.target.kind=TargetKind::Physical;
    predictor.step.occurrence=input.predictorOccurrence;
    predictor.order=10;predictor.slots={"predictor"};
    predictor.origin={OriginKind::Generated,request.preset+" transforms "+input.predictor.equation+" occurrence="+input.predictor.occurrence+" -> "+input.predictor.target.symbol+" @ 10; source="+predictor.origin.source};
    entries.erase(std::remove_if(entries.begin(),entries.end(),[&](const auto& node) {
        return (!input.replacedCommit.empty() && node.kind==ExecutionKind::Commit && node.id==input.replacedCommit)
            || std::any_of(input.fusedBalances.begin(),input.fusedBalances.end(),[&](const auto& call){return matches(node,call);});
    }),entries.end());
    if(request.presetKind!=FDM::PressureCouplingPreset::PISO
        && request.presetKind!=FDM::PressureCouplingPreset::SIMPLE
        && request.presetKind!=FDM::PressureCouplingPreset::PIMPLE)
        throw std::runtime_error("No selected built-in pressure algorithm.");
    const Provenance origin{OriginKind::BuiltinPreset,"builtin "+request.preset};
    const auto call=[&](const EquationCall& reference,Order order) {
        if(reference.equation.empty() || reference.target.symbol.empty())
            throw std::runtime_error("Pressure contribution has an empty equation/target reference.");
        ExecutionScope node;node.kind=ExecutionKind::EquationCall;
        node.step=reference;node.order=order;node.origin=origin;
        return node;
    };
    const auto loop=[&](const char* id,int repetitions,Order order) {
        ExecutionScope node;node.kind=ExecutionKind::Loop;
        node.id=id;node.repetitions=repetitions;node.order=order;node.origin=origin;
        return node;
    };
    const auto pressureCorrection=[&] {
        ExecutionScope correction;correction.id="pressureCorrection";correction.origin=origin;
        auto nonOrthogonal=loop("nonOrthogonal",request.nonOrthogonalCorrectors+1,10);
        nonOrthogonal.children.push_back(call(input.correction,10));
        correction.children.push_back(std::move(nonOrthogonal));
        correction.children.push_back(call(input.velocityUpdate,20));
        correction.children.push_back(call(input.pressureUpdate,30));
        correction.children.push_back(call(input.fluxUpdate,40));
        if(!input.publication.equation.empty())
            correction.children.push_back(call(input.publication,50));
        return correction;
    };
    const auto pressureLoop=[&] {
        auto pressure=loop("pressure",request.pressureCorrectors,20);
        pressure.children.push_back(pressureCorrection());
        return pressure;
    };
    if(request.presetKind==FDM::PressureCouplingPreset::PISO) {
        entries.push_back(std::move(predictor));
        if(!input.fluxPrediction.equation.empty())entries.push_back(call(input.fluxPrediction,15));
        auto correction=pressureLoop();correction.slots={"transport"};entries.push_back(std::move(correction));
        program.root.slots.push_back("coupling");
    } else {
        auto outer=loop("outer",request.outerCorrectors,20);outer.slots={"coupling"};
        if(input.representation==PressurePredictorKind::ConservativeMomentum)outer.terminationSignal=kPressureOuterConvergedSignal;
        else {
            outer.kind=ExecutionKind::While;outer.condition="res(pressure.outer)";
            auto& convergence=outer;convergence.minimumIterations=1;
            const std::array<std::pair<const char*,double>,2> criteria={{{"pressure.iterate",request.iterationTolerance},{"pressure.continuity",request.continuityTolerance}}};
            for(int component=0;component<2;++component) {
                ResidualDefinition residual;residual.id=criteria[component].first;
                residual.occurrence=input.convergence.occurrence.empty()?"pressure.convergence":input.convergence.occurrence;
                residual.target="p";residual.kind=ResidualKind::DiscreteEquation;residual.component=component;
                residual.normalization=ResidualNormalization::None;residual.absoluteTolerance=criteria[component].second;residual.relativeTolerance=0;
                residual.origin=origin;program.residuals.push_back(std::move(residual));
            }
            ResidualDefinition group;group.id="pressure.outer";group.members={"pressure.iterate","pressure.continuity"};group.origin=origin;
            program.residuals.push_back(std::move(group));
        }
        outer.children.push_back(std::move(predictor));
        if(!input.fluxPrediction.equation.empty())outer.children.push_back(call(input.fluxPrediction,15));
        if(request.presetKind==FDM::PressureCouplingPreset::SIMPLE) {
            auto correction=pressureCorrection();correction.order=20;
            outer.children.push_back(std::move(correction));
        } else outer.children.push_back(pressureLoop());
        outer.children.push_back(call(input.relaxation,30));
        auto consistency=call(input.fluxConsistency,40);consistency.slots={"transport"};outer.children.push_back(std::move(consistency));
        auto convergence=call(input.convergence,50);convergence.slots={"convergence"};
        if(input.representation!=PressurePredictorKind::ConservativeMomentum && convergence.step.occurrence.empty())convergence.step.occurrence="pressure.convergence";
        outer.children.push_back(std::move(convergence));
        entries.push_back(std::move(outer));
    }
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;
    commit.id="physicalStep.commit";commit.order=1000000;commit.origin=origin;
    entries.push_back(std::move(commit));

}

std::vector<NumericalBinding> pressureNumerics(const CouplingPresetRequest& request,
        const PressureCouplingInputs& input) {
    const bool conservative=input.representation==PressurePredictorKind::ConservativeMomentum;
    const auto binding=[&](const EquationCall& call,const char* method) {
        return NumericalBinding{call.equation,method,{},call.occurrence};
    };
    auto predictor=binding(input.predictor,conservative ? "ConservativePressureMomentum" : "CellExplicit");
    predictor.occurrence=input.predictorOccurrence;
    for (const auto& balance:input.fusedBalances) predictor.inputs.push_back(balance.equation);
    std::vector<NumericalBinding> bindings{predictor,
        binding(input.correction,conservative ? "ConservativePressureCorrection" : "PressureCorrection"),
        binding(input.velocityUpdate,conservative ? "ConservativeVelocityCorrection" : "VelocityCorrection"),
        binding(input.pressureUpdate,conservative ? "ConservativePressureUpdate" : "PressureUpdate"),
        binding(input.fluxUpdate,conservative ? "ConservativeFluxCorrection" : "FluxCorrection")};
    if (!input.fluxPrediction.equation.empty()) {
        auto flux=binding(input.fluxPrediction,"RhieChow");flux.response=input.predictorOccurrence;
        bindings.push_back(std::move(flux));
    }
    if (!input.publication.equation.empty()) bindings.push_back(binding(input.publication,"ConservativePressurePublication"));
    if (request.presetKind!=FDM::PressureCouplingPreset::PISO) {
        bindings.push_back(binding(input.relaxation,conservative ? "ConservativeFixedTimeRelaxation" : "FixedTimeRelaxation"));
        bindings.push_back(binding(input.fluxConsistency,conservative ? "ConservativeFluxConsistency" : "FluxConsistency"));
        bindings.push_back(binding(input.convergence,conservative ? "ConservativeResidualConvergence" : "ResidualConvergence"));
    }
    return bindings;
}

CouplingReport matchPressureCoupling(
        const CouplingPresetRequest& request,
        const RawEquationSystem& raw,const PressureCouplingInputs* inputs) {
    CouplingReport report;
    report.id = "coupling."+request.preset;
    report.preset = request.preset;
    report.requirements = {
        "momentum equation",
        "incompressibility / pressure-multiplier constraint",
        "pressure multiplier unknown"};

    if (inputs ? !rawHasEquation(raw,inputs->predictor.equation) : !hasMomentumEquation(raw)) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "momentum equation not present in the resolved equation system";
        return report;
    }
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw,inputs);
    if (!pressureConstraint) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "incompressibility / pressure-multiplier constraint not present";
        return report;
    }
    const std::string multiplier = pressureConstraint->multiplierUnknown;
    const bool pressureBound = multiplier.empty()
        ? rawHasUnknown(raw,"p")
        : rawHasUnknown(raw,multiplier);
    if (!pressureBound) {
        report.status = request.explicitlyRegistered
            ? CouplingStatus::Invalid : CouplingStatus::Inactive;
        report.reason =
            "the pressure constraint is not bound to a declared pressure "
            "unknown";
        return report;
    }
    // pressure-constraint formulation 需要生成算法压力方程；它由 transformer
    // 产生（pSimple / E_SHARED_PRESSURE），不是用户输入。
    const bool shared =
        pressureConstraint->id == "C_SHARED_PRESSURE";
    const auto pressureEquation=inputs ? inputs->correction.equation : (shared ? "E_SHARED_PRESSURE" : "pSimple");
    const bool hasPressureEquation=rawHasEquation(raw,pressureEquation);
    if (!hasPressureEquation) {
        report.derivedEquations = {
            pressureEquation};
    }
    if (request.outerCorrectors <= 0 || request.pressureCorrectors <= 0
        || request.nonOrthogonalCorrectors < 0) {
        report.status = CouplingStatus::Invalid;
        report.reason =
            "registered coupling preset has invalid corrector counts";
        return report;
    }
    if (request.presetKind == FDM::PressureCouplingPreset::PISO
        && request.outerCorrectors != 1) {
        report.status = CouplingStatus::Invalid;
        report.reason = "PISO does not define an outer fixed-point iteration; "
            "use PIMPLE for outerCorrectors > 1.";
        return report;
    }
    if (request.presetKind == FDM::PressureCouplingPreset::SIMPLE
        && request.pressureCorrectors != 1) {
        report.status = CouplingStatus::Invalid;
        report.reason = "SIMPLE uses one pressure correction per outer "
            "iteration; use PIMPLE for multiple inner pressure corrections.";
        return report;
    }
    // Final operations are reported from the compiled owned plan.
    report.status = CouplingStatus::Active;
    report.reason = "resolved equation/constraint system satisfies the "
                    "pressure-coupling requirements";
    return report;
}

CouplingReport contributePressureCoupling(
        SystemCompositionBuilder& system,
        const CouplingPresetRequest& request,
        const RawEquationSystem& raw,const PressureCouplingInputs* inputs) {
    CouplingReport report = matchPressureCoupling(request,raw,inputs);
    if (report.status != CouplingStatus::Active) return report;

    // 1. formulation / transformation contribution。
    //    Equation source 可能已经声明了同一 formulation（例如 Eulerian 模板
    //    本身要求 shared-pressure formulation），此时不重复请求。
    const ConstraintDescriptor* pressureConstraint =
        pressureConstraintOf(raw,inputs);
    const bool shared = pressureConstraint
        && pressureConstraint->id == "C_SHARED_PRESSURE";
    if (shared) {
        if (!system.requestsTransformation("sharedPressureConstraint")) {
            system.requestTransformation({
                "sharedPressureConstraint",
                request.preset+" shared-pressure transformation",
                110,true,request.origin});
        }
    } else if (!system.requestsTransformation("pressureConstraint")) {
        system.requestTransformation({
            "pressureConstraint",
            request.preset+" pressure transformation",
            100,true,request.origin});
    }

    // Coupling counts are authored only in native HOW.
    return report;
}

} // namespace SF::System
