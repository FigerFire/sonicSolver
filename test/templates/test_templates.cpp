#include "solver/system/SF_declarations.h"
#include "core/system/SF_scheduleIds.h"
#include "solver/system/SF_builtinState.h"
#include "solver/system/SF_singleFluidPreset.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/system/SF_presets.h"
#include <iostream>
using namespace SF;
using namespace SF::System;
namespace {
void check(bool ok,const char* reason) {if(!ok)throw std::runtime_error(reason);}
template<class F> void rejects(F fn) {
    try {fn();}catch(const std::exception&){return;}
    throw std::runtime_error("Expected a rejected declaration");
}
// Independent typed declarations: the reference does not call the built-in constructor.
std::vector<Equation> balances(bool viscous) {
    const BuiltinStateCatalog state;
    const auto temporal=[&](const char* q,const char* id){return FormulaExpr::op(Math::Op::Ddt,{state.expression(q)},id);};
    const auto div=[&](const char* q,const char* id){return FormulaExpr::op(Math::Op::Div,{state.expression(q)},id);};
    auto momentum=temporal("rhoU","momentum.time")+div("momentumFlux","momentum.convection");
    auto energy=temporal("rhoE","energy.time")+div("energyFlux","energy.convection");
    if(viscous) {
        momentum=momentum+FormulaExpr::op(Math::Op::Diffusion,{state.expression("mu"),state.expression("U")},"momentum.diffusion");
        energy=energy+FormulaExpr::op(Math::Op::Diffusion,{state.expression("conductivity"),state.expression("T")},"energy.diffusion");
    }
    return {{"continuity",temporal("rho","mass.time")+div("massFlux","mass.convection"),FormulaExpr::constantValue(0)},
        {"momentum",momentum,FormulaExpr::constantValue(0)},{"energy",energy,FormulaExpr::constantValue(0)}};
}
bool sameOperators(const FormulaExpr& a,const FormulaExpr& b) {
    if(a.kind()!=b.kind() || a.arguments().size()!=b.arguments().size())return false;
    if(a.kind()==FormulaExpr::Kind::Operator && (a.operation()!=b.operation() || a.occurrence()!=b.occurrence()))return false;
    for(std::size_t i=0;i<a.arguments().size();++i)
        if(!sameOperators(a.arguments()[i],b.arguments()[i]))return false;
    return true;
}
bool same(const ExecutionScope& a,const ExecutionScope& b) {
    if(a.kind!=b.kind || a.id!=b.id || a.order!=b.order || a.frequency!=b.frequency
        || a.repetitions!=b.repetitions || a.minimumIterations!=b.minimumIterations
        || a.terminationSignal!=b.terminationSignal || a.children.size()!=b.children.size()
        || a.step.equation!=b.step.equation || a.step.target.symbol!=b.step.target.symbol
        || a.step.target.kind!=b.step.target.kind || a.step.occurrence!=b.step.occurrence)return false;
    for(std::size_t i=0;i<a.children.size();++i)if(!same(a.children[i],b.children[i]))return false;
    return true;
}
ExecutionScope call(const EquationCall& ref,int order) {
    ExecutionScope node;node.kind=ExecutionKind::EquationCall;node.step=ref;node.order=order;return node;
}
ExecutionScope loop(const char* id,int n,int order) {
    ExecutionScope node;node.kind=ExecutionKind::Loop;node.id=id;node.repetitions=n;node.order=order;return node;
}
ExecutionProgram reference(const CouplingPresetRequest& req,const PressureCouplingInputs& in) {
    ExecutionProgram out;
    auto predictor=call(in.predictor,10);predictor.step.target=in.predictorTarget;predictor.step.occurrence=in.predictorOccurrence;
    const auto correction=[&] {
        ExecutionScope body;body.id="pressureCorrection";
        auto nonOrthogonal=loop("nonOrthogonal",req.nonOrthogonalCorrectors+1,10);
        nonOrthogonal.children.push_back(call(in.correction,10));body.children.push_back(nonOrthogonal);
        body.children.push_back(call(in.velocityUpdate,20));body.children.push_back(call(in.pressureUpdate,30));
        body.children.push_back(call(in.fluxUpdate,40));
        if(!in.publication.equation.empty())body.children.push_back(call(in.publication,50));
        return body;
    };
    if(req.presetKind==FDM::PressureCouplingPreset::PISO) {
        out.root.children.push_back(predictor);
        auto pressure=loop("pressure",req.pressureCorrectors,20);pressure.children.push_back(correction());out.root.children.push_back(pressure);
    } else {
        auto outer=loop("outer",req.outerCorrectors,20);outer.terminationSignal=kPressureOuterConvergedSignal;
        outer.children.push_back(predictor);
        if(req.presetKind==FDM::PressureCouplingPreset::SIMPLE) {
            auto body=correction();body.order=20;outer.children.push_back(body);
        } else {
            auto pressure=loop("pressure",req.pressureCorrectors,20);pressure.children.push_back(correction());outer.children.push_back(pressure);
        }
        outer.children.push_back(call(in.relaxation,30));outer.children.push_back(call(in.fluxConsistency,40));
        outer.children.push_back(call(in.convergence,50));out.root.children.push_back(outer);
    }
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="physicalStep.commit";commit.order=1000000;
    out.root.children.push_back(commit);return out;
}
void pressure() {
    for(const auto representation:{PressurePredictorKind::Velocity,PressurePredictorKind::ConservativeMomentum})
    for(const auto kind:{FDM::PressureCouplingPreset::SIMPLE,FDM::PressureCouplingPreset::PISO,FDM::PressureCouplingPreset::PIMPLE}) {
        auto in=pressureCouplingInputs({"renamedBalance",{"velocity"},"renamedOccurrence"},representation);
        in.predictorOccurrence="renamedPredictor";
        CouplingPresetRequest req;req.presetKind=kind;req.preset=FDM::toString(kind);
        req.outerCorrectors=3;req.pressureCorrectors=2;req.nonOrthogonalCorrectors=1;
        ExecutionProgram actual;actual.root.children.push_back(call(in.predictor,20));
        for(const auto& balance:in.fusedBalances)actual.root.children.push_back(call(balance,10));
        if(!in.replacedCommit.empty()) {
            ExecutionScope node;node.kind=ExecutionKind::Commit;node.id=in.replacedCommit;actual.root.children.push_back(node);
        }
        applyPressureExecution(actual,req,in);
        check(same(actual.root,reference(req,in).root),"Pressure schedule differs from independent typed reference");
    }
}
}
int main() {
    try {
        for(const bool viscous:{false,true}) {
            const auto actual=Preset::singleFluidEquations({viscous});
            const auto expected=balances(viscous);
            check(actual.registeredEquations.size()==expected.size(),"Missing balance");
            for(std::size_t i=0;i<expected.size();++i) {
                const auto& equation=actual.registeredEquations[i];
                check(equation.id==expected[i].id && canonicalFormula(equation)==canonicalFormula(expected[i]),"WHAT differs from manual declarations");
                check(sameOperators(equation.lhs,expected[i].lhs),"Physical operator occurrence changed");
            }
            check(actual.requiredStates.size()==(viscous?5:3),"Wrong property requests");
        }
        ExecutionScope custom;custom.kind=ExecutionKind::EquationCall;
        custom.step={"customEquation",{"a"},"customStep"};custom.frequency=ExecutionFrequency::Stage;
        const auto program=Compose::explicitProgram({custom});
        check(program.root.children.size()==2 && program.root.children.front().step.occurrence=="customStep","Built-in changes authored occurrence");
        check(program.root.children.front().frequency==ExecutionFrequency::Stage,"Built-in loses explicit frequency");
        check(program.root.children.back().kind==ExecutionKind::Commit,"Built-in loses physical commit");
        check(program.recipeScopes,"Built-in loses recipe scope expansion");
        const auto defaults=Preset::singleFluidDefaults();
        check(defaults.solutionVariables==std::vector<std::string>{"rho","rhoU","rhoE"}
            && defaults.equations==std::vector<std::string>{"Continuity","Momentum","Energy"}
            && defaults.algorithm=="Explicit","Built-in package defaults changed");
        pressure();
        using P=Model::Parameters;
        const Provenance origin{OriginKind::User,"native declarations"};
        rejects([&]{Declarations::execution({{"template","explicit"}},2,origin);});
        rejects([&]{Declarations::execution({{"root",P::array()}},2,origin);});
        rejects([&]{Declarations::execution({{"root",P::array({{{"$scope","predictor"}}})}},2,origin);});
        rejects([&]{Declarations::execution({{"root",P::array({{{"equation","E"},{"target","a"},{"frequency","Unknown"}}})}},2,origin);});
        rejects([&]{Declarations::numerics({{"bindings",P::array({{{"equation","E"},{"time","unknown"}}})}},origin);});
        std::cout<<"Built-in declaration contract passed\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
