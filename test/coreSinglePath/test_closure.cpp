#include "solver/system/SF_numericalCompiler.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/system/SF_builtinState.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "models/physics/gravity/SF_accelerationProvider.h"
#include "models/ibm/SF_ibmSystemContribution.h"
#include <iostream>
using namespace SF;
using namespace SF::System;
namespace {
int checks=0;
void require(bool ok,const std::string& why){++checks;if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(const std::string& message,F f) {
    try{f();}catch(const std::exception& e){require(std::string(e.what()).find(message)!=std::string::npos,e.what());return;}
    throw std::runtime_error("Expected rejection: "+message);
}
void sourceIdentity() {
    using E=FormulaExpr;
    BuiltinStateCatalog states;
    const auto force=E::symbol("gravity",Math::signature(1,{1,-2,-2,0,0,0,0}));
    CompiledExecutionProgram program;CompiledEquationCall call;
    call.source.equation="force";call.source.target.symbol="rhoU";
    call.equationMethod="ConservativeResidual";call.spatialTerms=true;call.calls={{"force","rhoU"}};
    program.steps.push_back(call);
    auto providers=TermProviderCatalog::builtIn();ZoneVectorSetting g;g.value=Vector3(.01,0,0);
    for(auto p:Physics::Gravity::termProviders({g}))providers.addSource(std::move(p));
    FDM::NumericalRecipeSet recipes;recipes.time=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);
    const auto time=builtinTemporalMethods().at(recipes.time.id()).compile(recipes.time);
    std::string canonical,provider;
    for(auto origin:{OriginKind::BuiltinDefault,OriginKind::BuiltinPreset,OriginKind::Model,OriginKind::User,OriginKind::Generated}) {
        const auto lhs=E::op(Math::Op::Ddt,{states.expression("rhoU")});
        Equation equation{"force",lhs,E::op(Math::Op::Source,{force},"force.source"),{origin,"different provenance"}};
        ExecutableEquationSystem system;system.registry.add(equation);
        const auto numerical=NumericalCompiler::compile(system,program,recipes,providers,time);
        require(numerical.operators.size()==1 && numerical.operators[0].conservativeSource,"Explicit source lost execution");
        if(canonical.empty()){canonical=canonicalFormula(equation);provider=numerical.operators[0].provider;}
        require(canonicalFormula(equation)==canonical && numerical.operators[0].provider==provider,"Provenance changed mathematics/provider");
        ExecutableEquationSystem naked;naked.registry.add({"force",lhs,force,{origin,"same origin"}});
        rejects("explicit Source AST",[&]{NumericalCompiler::compile(naked,program,recipes,providers,time);});
    }
}
struct CountingModel:Physics::FluidStateModel::Model {
    std::shared_ptr<const Physics::FluidStateModel::Model> original=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,1e-5,.72);
    mutable int calls=0;
    auto family()const->Physics::FluidStateModel::Family override{return original->family();}
    const State::StateLayout& layout()const override{return original->layout();}
    const std::vector<Physics::FluidStateModel::Variable>& primaryVariables()const override{return original->primaryVariables();}
    int densityVariableCount()const override{return original->densityVariableCount();}
    int momentumIndex(int c)const override{return original->momentumIndex(c);}
    int energyIndex()const override{return original->energyIndex();}
    Physics::FluidStateModel::ThermodynamicState close(const double* q,int n)const override{++calls;return original->close(q,n);}
};
void cachedProjections() {
    Field field;field.setup(1,1,1,0,5);auto model=std::make_shared<CountingModel>();field.setStateModel(model);
    field(0,0,0,0)=2;field(0,0,0,1)=6;field(0,0,0,2)=field(0,0,0,3)=0;field(0,0,0,4)=250009;
    BuiltinStateCatalog catalog;StateRegistry symbols;
    for(auto id:{"rho","rhoU","rhoE"})symbols.add(catalog.solution(id));
    for(auto id:{"U","p","T","h","mu","nu","conductivity"})symbols.add(catalog.at(id));
    for(auto id:{"U","p","T"})require(!catalog.solution(id).producer && !catalog.solution(id).readOnly(),"Primary view retained closure authority");
    State::StateBundle bundle;bundle.patches={&field};bundle.distributed.add(State::conservativeView("conservative",0,field));
    std::vector<CompiledStateView> views;
    for(auto id:{"rho","rhoU","rhoE","p","T","h","mu","nu","conductivity"}) {
        const auto& s=symbols.at(id);CompiledStateView v;v.symbol=id;v.kind=StateViewKind::Stage;v.stage=0;
        v.components=s.components;v.componentOffset=s.componentOffset;v.owner=StateViewOwner::NumericalProvider;views.push_back(v);
    }
    auto realized=realizeState(symbols,{},bundle,views);require(model->calls==0,"Binding closed EOS");
    auto expected=model->original->close(std::array<double,5>{2,6,0,0,250009}.data(),5);
    const auto verify=[&] {
        for(auto pair:std::vector<std::pair<std::string,double>>{{"p",expected.pressure},{"T",expected.temperature},
                {"h",(expected.internalEnergyDensity+expected.pressure)/expected.density},{"mu",expected.dynamicViscosity},
                {"nu",expected.dynamicViscosity/expected.density},{"conductivity",expected.thermalConductivity}}) {
            require(realized.at(pair.first).fields[0]->read(0,0)==pair.second,"Changed original cache projection: "+pair.first);
            rejects("read-only",[&]{realized.at(pair.first).fields[0]->write(0,0,0);});
        }
    };
    verify();require(model->calls==1,"Projections created additional EOS caches");
    std::vector<State::DistributedFieldView> stages;
    stages.reserve(3);
    for(auto id:{"rho","rhoU","rhoE"}) {
        const auto s=symbols.at(id);auto v=State::conservativeView(id,0,field);v.components=s.components;
        v.read=[&field,s](int,int c){return field(0,0,0,s.componentOffset+c);};
        stages.push_back(v);realized.bindView(id,StateViewKind::Stage,stages.back(),0);
    }
    realized.bindDerivedViews();auto stage=realized.view("T",StateViewKind::Stage,0).fields[0]->read;
    rejects("wrong Stage",[&]{stage(0,0);});
    field(0,0,0,4)=500009;field.invalidateThermodynamicCache();
    expected=model->original->close(std::array<double,5>{2,6,0,0,500009}.data(),5);
    require(model->calls==1,"Q publication closed EOS eagerly");
    realized.beginStage(0,1,.1);realized.stageReady();
    require(stage(0,0)==expected.temperature && model->calls==2,"Stage producer lost original cache invalidation");
    verify();require(model->calls==2,"Physical/Stage acquired separate caches");
    realized.rhsReady();realized.finishStage();rejects("wrong Stage",[&]{stage(0,0);});
    auto expired=realized.at("p").fields[0]->read;
    realized=StateRealization{};rejects("expired",[&]{expired(0,0);});
    StateRegistry transportOnly;transportOnly.add(catalog.at("nu"));
    auto coefficient=realizeState(transportOnly,{},bundle);
    require(transportOnly.symbols().size()==1 && coefficient.at("nu").fields[0]->read(0,0)==5e-6,"Cache coefficient invented conservative STATE");
}
void wallGeometry() {
    FDM::ImmersedAlgorithmDescriptor descriptor;descriptor.enforcement=FDM::IBMEnforcement::GhostCell;
    descriptor.wallClosure=FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic;
    SystemContribution contribution;IBM::SystemContribution::contribute(contribution,descriptor);
    StateRegistry symbols;for(const auto& s:contribution.states)symbols.add(s);
    Field field;field.setupGeometry(1,1,1,0);field.setWallDistance(0,0,0,.125);
    State::StateBundle bundle;bundle.patches={&field};auto state=realizeState(symbols,{},bundle);
    require(state.at("immersed.wallDistance").fields[0]->read(0,0)==.125,"Geometry producer requires EOS or copied wall distance");
    field.setWallDistance(0,0,0,.25);
    require(state.at("immersed.wallDistance").fields[0]->read(0,0)==.25,"Geometry producer stopped aliasing owner");
}
}
int main(){try{sourceIdentity();cachedProjections();wallGeometry();std::cout<<"Core single-path checks passed: "<<checks<<'\n';}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
