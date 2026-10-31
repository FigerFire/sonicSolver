#include "SF_numericalSelection.h"
#include "SF_buildRequest.h"
#include "SF_runtimeRequirements.h"
#include "core/system/SF_operationIds.h"
#include <cmath>
#include <algorithm>
#include <set>
#include "SF_formulaCompiler.h"
namespace SF::System {
namespace {

void overlay(NumericalBinding& destination,const NumericalBinding& source) {
    std::set<std::pair<std::string,std::string>> sites;
    for(const auto& op:source.operators)if(!sites.emplace(op.formula,op.occurrence).second)
        throw std::runtime_error("Duplicate numerical operator selection in one WHICH declaration: "+op.formula+" / "+op.occurrence);
    if(!source.method.empty())destination.method=source.method;
    if(source.inputsDeclared || !source.inputs.empty()){destination.inputs=source.inputs;destination.inputsDeclared=true;}
    for(const auto& parameter:source.parameters)destination.parameters[parameter.first]=parameter.second;
    if(source.time)destination.time=source.time;
    for(const auto& op:source.operators) {
        const auto same=[&](const auto& previous){return previous.formula==op.formula && previous.occurrence==op.occurrence && previous.operation==op.operation;};
        destination.operators.erase(std::remove_if(destination.operators.begin(),destination.operators.end(),same),destination.operators.end());
        destination.operators.push_back(op);
    }
    destination.origin=source.origin;
}
}
ResolvedBinding resolveNumerics(const EquationCall& call,const std::vector<NumericalBinding>& bindings,bool block,const ExecutableEquationSystem* system) {
    ResolvedBinding result;
    const auto apply=[&](NumericalBinding& destination,const NumericalBinding& source) {
        auto normalized=source;
        for(auto& op:normalized.operators)if(op.formula.empty() || op.formula=="*")op.formula=call.equation;
        if(system)resolveOperatorReferences(*system,call,normalized);
        overlay(destination,normalized);
    };
    for(int specificity=0;specificity<3;++specificity) {
        std::vector<const NumericalBinding*> defaults,authored;
        for(const auto& binding:bindings) {
            if(binding.equation!="*" && binding.equation!=call.equation)continue;
            if(!binding.occurrence.empty() && binding.occurrence!=call.occurrence)continue;
            if(block && binding.occurrence.empty())continue;
            const int rank=!binding.occurrence.empty()?2:binding.equation!="*"?1:0;
            if(rank!=specificity)continue;
            (binding.explicitChoice?authored:defaults).push_back(&binding);
        }
        if(authored.size()>1)throw std::runtime_error("Conflicting highest-priority numerical bindings for occurrence: "+call.occurrence);
        if(defaults.size()>1) {
            NumericalBinding covering;
            if(!authored.empty())apply(covering,*authored.front());
            for(const auto& higher:bindings) {
                const int rank=!higher.occurrence.empty()?2:higher.equation!="*"?1:0;
                if(rank<=specificity || (higher.equation!="*" && higher.equation!=call.equation)
                    || (!higher.occurrence.empty() && higher.occurrence!=call.occurrence))continue;
                apply(covering,higher);
            }
            const auto uncovered=[&](const NumericalBinding& a,const NumericalBinding& b) {
                if(!a.method.empty() && !b.method.empty() && covering.method.empty())return true;
                if(a.inputs!=b.inputs && !covering.inputsDeclared && covering.inputs.empty())return true;
                if(a.time!=b.time && !covering.time)return true;
                auto parameters=a.parameters;for(const auto& item:b.parameters)parameters.emplace(item);
                for(const auto& item:parameters)if(!covering.parameters.count(item.first)
                    && (a.parameters.count(item.first)!=b.parameters.count(item.first)
                        || (a.parameters.count(item.first) && a.parameters.at(item.first)!=b.parameters.at(item.first))))return true;
                std::map<std::pair<std::string,std::string>,std::string> left,right;
                for(const auto& op:a.operators)left[{op.formula,op.occurrence}]=op.provider;
                for(const auto& op:b.operators)right[{op.formula,op.occurrence}]=op.provider;
                for(const auto& op:covering.operators){left.erase({op.formula,op.occurrence});right.erase({op.formula,op.occurrence});}
                return left!=right;
            };
            for(std::size_t i=1;i<defaults.size();++i)if(uncovered(*defaults.front(),*defaults[i]))
                throw std::runtime_error("Ambiguous default WHICH for occurrence: "+call.occurrence+" from "+defaults.front()->origin.source+" and "+defaults[i]->origin.source);
        }
        if(!defaults.empty()){apply(result.value,*defaults.front());result.origins.push_back(defaults.front()->origin);}
        if(!authored.empty()) {apply(result.value,*authored.front());result.origins.push_back(authored.front()->origin);}
    }
    if(result.value.method.empty())throw std::runtime_error("Implementation capability missing: no provider for equation occurrence "+call.equation+" / "+call.occurrence+(block?"; whole block requires an explicit Step binding":""));
    result.value.equation=call.equation;result.value.occurrence=call.occurrence;
    return result;
}
void resolveOperatorReferences(const ExecutableEquationSystem& system,const EquationCall& call,NumericalBinding& binding) {
    std::set<std::pair<std::string,std::string>> used;
    for(auto& op:binding.operators) {
        if(op.formula.empty() || op.formula=="*")op.formula=call.equation;
        if(op.formula!=call.equation)throw std::runtime_error("WHICH operator belongs to another equation at Step "+call.occurrence+": "+op.formula);
        const auto& equation=system.registry.at(op.formula);
        std::vector<const FormulaExpr*> exact,abbreviated;
        const auto inspect=[&](const auto& self,const FormulaExpr& expression)->void {
            if(expression.kind()==FormulaExpr::Kind::Operator) {
                if(expression.occurrence()==op.occurrence)exact.push_back(&expression);
                if(Math::name(expression.operation())==op.occurrence)abbreviated.push_back(&expression);
            }
            for(const auto& child:expression.arguments())self(self,child);
        };
        inspect(inspect,equation.lhs);inspect(inspect,equation.rhs);
        const auto& matches=exact.empty()?abbreviated:exact;
        if(matches.size()!=1)throw std::runtime_error("Missing/ambiguous operator occurrence: "+call.occurrence+" / "+op.formula+" / "+op.occurrence);
        op.occurrence=matches.front()->occurrence();
        if(op.occurrence.empty() || !used.emplace(op.formula,op.occurrence).second)throw std::runtime_error("Duplicate/unnamed numerical operator binding: "+call.occurrence);
    }
}

const PressureNumericalConfig* requiredPressureInputs(const Program& program) {
    const auto predictor=std::find_if(program.solve.compiledProgram.steps.begin(),program.solve.compiledProgram.steps.end(),
        [](const auto& call){return call.equationMethod=="PressureMomentum";});
    if(predictor==program.solve.compiledProgram.steps.end())return nullptr;
    const auto context="Pressure resources step="+predictor->source.occurrence+" equation="+predictor->source.equation+": ";
    if(!program.numerics.pressureMomentum || program.numerics.pressureMomentum->equation!=predictor->source.equation
        || program.numerics.pressureMomentum->occurrence!=predictor->source.occurrence)
        throw std::runtime_error(context+"missing/mismatched frozen momentum term contract.");
    if(!program.numerics.pressureOperator || program.numerics.pressureFaceCoupling!=PressureFaceCoupling::RhieChow)
        throw std::runtime_error(context+"missing frozen pressure inputs or persistent Rhie-Chow face-flux treatment.");
    for(const char* symbol:{"U","p","rho","phi"})if(!program.state.contains(symbol))
        throw std::runtime_error(context+"missing required STATE: "+symbol);
    const auto& velocity=program.state.at("U");const auto& pressure=program.state.at("p");
    const auto& density=program.state.at("rho");const auto& flux=program.state.at("phi");
    if(predictor->target.symbol!=velocity.id || velocity.storageBinding!=StorageBinding::PackedDistributed
        || velocity.storageKey!="velocity" || velocity.components!=3 || velocity.componentOffset!=0 || velocity.readOnly()
        || pressure.storageBinding!=StorageBinding::NamedDistributed || pressure.storageKey!="pressure"
        || pressure.components!=1 || pressure.componentOffset!=0 || pressure.readOnly()
        || density.storageKey!="rhoConst" || !density.constantValue
        || !std::isfinite(*density.constantValue) || *density.constantValue<=0.0)
        throw std::runtime_error(context+"requires writable U[3]/velocity and p[1]/pressure, and positive declared rhoConst; conservative Q and phase layouts are distinct.");
    if(flux.location!=VariableLocation::EulerianFace || flux.ownership!=OwnershipKind::CanonicalFace
        || flux.storageKey!="pressureFaceFlux")
        throw std::runtime_error(context+"requires canonical pressureFaceFlux owned by the pressure workspace.");
    const auto hasView=[&](std::string_view symbol,StateViewKind kind) {
        return std::any_of(program.solve.compiledProgram.stateViews.begin(),program.solve.compiledProgram.stateViews.end(),
            [&](const auto& view){return view.symbol==symbol && view.kind==kind && view.owner==StateViewOwner::NumericalProvider;});
    };
    if(!hasView("deltaP",StateViewKind::Correction) || !hasView(flux.id,StateViewKind::Physical))
        throw std::runtime_error(context+"requires numerical-owner deltaP Correction and persistent physical face-flux views.");
    for(const char* id:{OpIds::PressurePrepare,OpIds::MomentumAssemble,OpIds::MomentumSolve,
            OpIds::PressureBoundaryPrepare,OpIds::PressureAssemble,OpIds::PressureSolve,OpIds::PressureUpdatePrepare,
            OpIds::VelocityCorrect,OpIds::FluxCorrect,OpIds::PressureCorrectionCommit,OpIds::PressureStepCommit,OpIds::TimeCommit}) {
        const auto binding=std::find_if(program.runtime.operationBindings.begin(),program.runtime.operationBindings.end(),
            [&](const auto& item){return item.operation==id;});
        if(!requiresOperation(program.runtime,id) || binding==program.runtime.operationBindings.end())
            throw std::runtime_error(context+"missing pressure operation requirement/binding: "+id);
        // Labels verify the frozen implementation, never select an application family.
        const auto expected=std::string_view(id)==OpIds::FluxCorrect?"flow.rhie-chow":"flow.pressure-operators";
        if(binding->status==BindingStatus::Resolved && binding->provider!=expected)
            throw std::runtime_error(context+"conflicting operation/face-flux owner: "+id);
    }
    return &*program.numerics.pressureOperator;
}

NumericalSelection selectNumerics(const FDM::SolverConfig& config,const BuildRequest& request,
        const ExecutableEquationSystem& equations,std::vector<NumericalBinding> bindings) {
    NumericalSelection selection;
    selection.bindings=std::move(bindings);
    selection.recipes=config.numerics.recipes;
    selection.recipes.time=config.numerics.timeRecipe;
    if (!selection.recipes.diffusion && !config.numerics.termRecipesDeclared)
        selection.defaultDiffusion=FDM::builtInDiffusionRecipe(config.numerics.viscous);
    selection.dt.cfl=config.numerics.cfl;
    selection.dt.maxDeltaT=config.numerics.maxDeltaT;
    selection.phaseTransport.convection=config.pressure.phaseTransport.convection;
    selection.phaseTransport.sourceCfl=config.pressure.phaseTransport.sourceCfl;
    selection.termProviders=request.termProviders;
    selection.operatorBindings=request.operatorBindings;
    if (std::any_of(selection.bindings.begin(),selection.bindings.end(),[](const auto& binding) {
            return binding.method=="PressureMomentum";
        })) {
        selection.pressureFaceCoupling=PressureFaceCoupling::RhieChow;
        PressureNumericalConfig pressure;
        pressure.dynamicViscosity=config.numerics.dynamicViscosity;
        pressure.velocityRelaxation=config.pressure.coupling.momentumRelaxation;
        pressure.pressureRelaxation=config.pressure.coupling.pressureRelaxation;
        pressure.relativeTolerance=config.pressure.relativeTolerance;
        pressure.absoluteTolerance=config.pressure.absoluteTolerance;
        pressure.reference=config.pressure.reference;
        pressure.pressureLinear=config.pressure.linear.pressure;
        pressure.velocityBoundary=config.boundaries.velocity;
        pressure.pressureBoundary=config.boundaries.pressure;
        selection.pressureOperator=std::move(pressure);
    }
    if(std::any_of(selection.bindings.begin(),selection.bindings.end(),[](const auto& binding) {
        return binding.method=="EulerianPhaseContinuity";
    })) {
        FDM::validatePressureCorrectionConfig(config.pressure);
        selection.phaseOperator=PhaseNumericalConfig{config.pressure.coupling.momentumRelaxation,
            config.pressure.coupling.pressureRelaxation,config.pressure.reference,config.pressure.linear,
            config.boundaries.pressure};
    }
    return selection;
}
} // namespace SF::System
