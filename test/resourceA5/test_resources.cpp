#define main standaloneReferenceMain
#include "test/test_standaloneScalar.cpp"
#undef main
#include "app/application/adapters/SF_scalarResources.h"
#include "infrastructure/io/SF_caseConfig.h"
#include "solver/system/SF_numericalCompiler.h"

namespace {
struct LocalLowering final : System::IProvider {
    mutable int calls=0;
    std::string_view id()const override{return "test.localLowering";}
    System::CompiledEquationCall compile(const System::ExecutableEquationSystem&,const System::EquationCall&,
        const System::NumericalBinding&)const override{throw std::runtime_error("not used");}
    void lowerSpatial(const System::ExecutableEquationSystem& system,const System::CompiledEquationCall& call,
        const FDM::NumericalRecipeSet&,const System::TermProviderCatalog&,
        const std::vector<System::FormulaOperatorBinding>&,System::CompiledNumericalSystem& numerical)const override {
        const auto& equation=system.registry.at(call.source.equation);
        require(equation.lhs.kind()==System::FormulaExpr::Kind::Negate,"AST changed before selected implementation");
        ++calls;numerical.providerRequirements.push_back(std::string(id()));
    }
};
void delegatedLowering() {
    using E=System::FormulaExpr;
    System::ExecutableEquationSystem system;
    const auto target=E::symbol("a",Math::signature(0,{}));
    system.registry.add({"signed",E::negate(E::op(Math::Op::Ddt,{target})),E::constantValue(0)});
    System::CompiledExecutionProgram program;System::CompiledEquationCall call;
    call.source.equation="signed";call.equationMethod="test.localLowering";call.spatialTerms=true;program.steps.push_back(call);
    LocalLowering method;System::ProviderRegistry providers;providers.add(method);
    FDM::NumericalRecipeSet recipes;const auto time=System::builtinTemporalMethods().at(recipes.time.id()).compile(recipes.time);
    const auto lowered=System::NumericalCompiler::compile(system,program,recipes,{},time,providers);
    require(method.calls==1 && lowered.providerRequirements==std::vector<std::string>{"test.localLowering"},"Public compiler bypassed selected lowering");
    rejects([&]{System::NumericalCompiler::compile(system,program,recipes,{},time,System::builtinProviders());});
}
System::ResolvedSimulationSystem pair(FDM::SolverConfig& config) {
    auto a=request("a",.03),b=request("b",.08);
    a.userContributions.front().records.front().id="user.scalar.a";
    b.userContributions.front().records.front().id="user.scalar.b";
    a.userContributions.front().execution.pop_back();
    a.composition.solutionVariables.push_back("b");
    a.userContributions.push_back(b.userContributions.front());
    return System::build(config,a);
}
void resourceBinding() {
    FDM::SolverConfig config;config.numerics.maxDeltaT=.0005;
    config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto system=pair(config);
    Field geometry;geometry.setupGeometry(9,1,1,1);
    for(int k=0;k<geometry.MZ();++k)for(int j=0;j<geometry.MY();++j)for(int i=0;i<geometry.MX();++i)geometry.X(i,j,k)=(i-1)/8.;
    geometry.getAllSets()["Left"]={geometry.getIdx(1,1,1)};geometry.getAllSets()["Right"]={geometry.getIdx(9,1,1)};
    CaseConfig input;auto description=std::make_shared<Model::Description>();
    for(auto name:{"a","b"}) {
        Model::FieldDescriptor f;f.name=f.symbol=name;f.initial=name==std::string("a")?1.:.5;
        f.boundaries={{"Left",{{"type","fixedValue"},{"value",0.}}},{"Right",{{"type","fixedValue"},{"value",0.}}}};
        description->fields.push_back(f);
    }
    input.modelDescription=description;
    State::StateBundle bound;bound.patches={&geometry};
    auto resources=Application::Execution::Detail::bindScalarResources(geometry,bound,system,input);
    require(resources.arrays.size()==2 && resources.ports.size()==2 && resources.instances.size()==2,"Wrong number of physical owners/ports");
    require(geometry.NVar()==0 && !geometry.stateModel(),"Binder allocated a hidden Q/EOS");
    // Independent resources bypass the application binder, but execute the same frozen equations.
    ScalarField a,b;a.setupLike(geometry,"a",1);b.setupLike(geometry,"b",.5);
    State::StateBundle explicitState;explicitState.patches={&geometry};
    explicitState.distributed.add(State::scalarView("a",0,geometry,a));explicitState.distributed.add(State::scalarView("b",0,geometry,b));
    auto da=std::make_shared<FDM::ScalarEquationData>(),db=std::make_shared<FDM::ScalarEquationData>();
    da->target="a";db->target="b";
    for(auto* d:{da.get(),db.get()}) for(int side=0;side<2;++side) {
        d->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;
        d->boundary[side].value=[](double,double,double,double){return 0.;};
    }
    FDM::SolverServices boundServices,explicitServices;boundServices.scalarInstances=resources.instances;
    explicitServices.scalarInstances={{"advance.a",da},{"advance.b",db}};
    Time::RunControl control;control.endTime=.005;control.writeByStep=true;control.writeIntervalSteps=100;
    Application::Execution::executeEquations(config,system,bound,boundServices,control);
    Application::Execution::executeEquations(config,system,explicitState,explicitServices,control);
    require(bound.time==explicitState.time && bound.dt==explicitState.dt && bound.step==10,"Resource binding changed clock");
    require(resources.arrays[0]->values()==a.values() && resources.arrays[1]->values()==b.values(),"Binder vs manual physical storage changed numerical result");
    rejects([&]{Application::Execution::Detail::bindScalarResources(geometry,bound,system,input);});
    for(int mode=0;mode<5;++mode) {
        auto wrong=system;auto& requirement=wrong.solvePlan.compiledProgram.steps.front().physicalResources.front();
        if(mode==0)requirement.owner="other.owner";
        if(mode==1)requirement.type=typeid(double);
        if(mode==2)requirement.writable=true;
        if(mode==3)requirement.id="wrong.occurrence";
        if(mode==4)wrong.solvePlan.compiledProgram.steps.front().physicalResources.push_back(requirement);
        State::StateBundle candidate;candidate.patches={&geometry};
        rejects([&]{auto result=Application::Execution::Detail::bindScalarResources(geometry,candidate,wrong,input);
            System::validatePhysicalResources(wrong.solvePlan.compiledProgram,result.views);});
        require(candidate.distributed.select("a",State::HaloSyncStage::None).empty() && candidate.distributed.select("b",State::HaloSyncStage::None).empty(),"Failed resource binding mutated active STATE");
    }
    for(int mode=0;mode<3;++mode) {
        auto missing=input;auto fields=std::make_shared<Model::Description>(*description);missing.modelDescription=fields;
        if(mode==0)fields->fields.clear();
        if(mode==1)fields->fields.front().initial=nullptr;
        if(mode==2)fields->fields.front().boundaries.erase("Left");
        State::StateBundle candidate;candidate.patches={&geometry};
        rejects([&]{Application::Execution::Detail::bindScalarResources(geometry,candidate,system,missing);});
    }
    const auto weak=resources.instances.front().data;resources.ports.clear();
    require(weak.expired(),"Services secretly own boundary-port lifetime");
}
}
int main() {
    try{delegatedLowering();resourceBinding();std::cout<<"A5 selected lowering and resource contracts passed\n";return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
