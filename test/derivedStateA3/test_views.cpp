#include "solver/system/SF_equationContribution.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/system/SF_stateRealization.h"
#include "solver/system/SF_builtinState.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/algorithm/SF_patchWorkspace.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include <fstream>
#include <filesystem>
#include <iostream>
using namespace SF;
using namespace SF::System;
namespace {
int checks=0;
void require(bool ok,const std::string& why) {++checks;if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(const std::string& expected,F action) {
    try {action();} catch(const std::exception& e) {
        require(std::string(e.what()).find(expected)!=std::string::npos,"Wrong rejection: "+std::string(e.what()));
        std::cout<<"Rejected ["<<expected<<"] "<<e.what()<<'\n';return;
    }
    throw std::runtime_error("Expected rejection: "+expected);
}
StateSymbol scalar(std::string id) {
    StateSymbol s;s.id=s.name=s.storageKey=std::move(id);s.signature=Math::signature(0,{0,1,0,0,0,0,0});
    s.storageBinding=StorageBinding::NamedDistributed;return s;
}
StateSymbol derived(std::string id="q") {
    auto s=scalar(std::move(id));s.role=StateRole::Derived;s.evaluation=StateEvaluation::Materialized;
    s.initializationRequired=s.boundaryRequired=s.restartEligible=false;s.dependencies={"a","b"};
    auto p=std::make_shared<StateProducer>();p->id="test.affine";
    p->expression=FormulaExpr::constantValue(2)*scalar("a").expression()+scalar("b").expression();
    p->versions={StateViewKind::Physical,StateViewKind::Stage};
    p->bind=[](Field& geometry,const auto& inputs,State::DistributedFieldView& output)->std::function<void()> {
        require(inputs.size()==2,"Missing input bindings");
        const auto a=inputs[0],b=inputs[1],out=output;
        return [a,b,out,field=&geometry] {
            const int g=field->NG();
            for(int k=g;k<g+field->NZ();++k)for(int j=g;j<g+field->NY();++j)for(int i=g;i<g+field->NX();++i) {
                const int cell=field->getIdx(i,j,k);out.write(cell,0,2*a.read(cell,0)+b.read(cell,0));
            }
        };
    };
    s.producer=std::move(p);return s;
}
struct UpdateContract {ResourceViewRequirement residual{"residual","patch",typeid(Residual),true};};
class UpdateMethod final:public IProvider {
public:
    std::string_view id()const override{return "test.affine.update";}
    std::string_view runtimeProvider()const override{return "test.affine";}
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,const NumericalBinding&)const override {
        const auto& s=system.state.at(call.target.symbol);const auto& equation=system.registry.at(call.equation);
        if(!s.producer || canonicalFormula(equation)!=canonicalFormula(Equation{equation.id,s.expression(),s.producer->expression}))
            throw std::runtime_error("Unsupported affine mathematics");
        CompiledEquationCall result;result.source=call;result.target.symbol=call.target.symbol;
        result.target.kind=call.target.kind;result.target.producer=s.producer;
        result.equationMethod=id();result.backendOperation="derive."+call.occurrence;
        result.calls={{equation.id,s.id}};result.reads=s.dependencies;result.writes={s.id};
        result.providerContract=UpdateContract{};return result;
    }
};
CompiledStateView demand(std::string id,std::string storage,StateViewKind kind,int stage=-1) {
    CompiledStateView v;v.symbol=std::move(id);v.storage=std::move(storage);v.kind=kind;v.stage=stage;
    v.owner=StateViewOwner::NumericalProvider;return v;
}
void extension() {
    SystemContribution contribution;contribution.addState(scalar("a"));contribution.addState(scalar("b"));
    contribution.addState(scalar("unrelated"));auto q=derived();contribution.addState(q);
    contribution.addEquation({"affine",q.expression(),q.producer->expression});
    RawEquationSystem raw;std::vector<TransformationDescriptor> transformations;
    SystemCompositionBuilder builder(raw,transformations,{OriginKind::User,"A3 internal extension"});builder.applyContribution(contribution);
    ExecutableEquationSystem system;system.state=raw.state;system.registry=raw.registry;
    ProviderRegistry providers;UpdateMethod method;providers.add(method);
    ExecutionProgram how{{"affine",{"q"},"current"}};how.explicitOrder=true;
    auto compiled=compileExecutionProgram(system,how,{{"affine","test.affine.update"}},providers);
    for(int stage=0;stage<2;++stage)for(const auto& id:{"a","b","q"})
        compiled.stateViews.push_back(demand(id,std::string(id)+".stage",StateViewKind::Stage,stage));
    compileDerivedViews(system.state,compiled.stateViews);
    Field geometry;geometry.setupGeometry(3,1,1,1);const int cell=geometry.getIdx(1,1,1);
    State::StateBundle bundle;bundle.patches={&geometry};
    ScalarField a,b,as,bs,qc,qs,other;
    for(auto pair:std::vector<std::pair<std::string,ScalarField*>>{{"a",&a},{"b",&b},{"a.stage",&as},{"b.stage",&bs},{"q",&qc},{"q.stage",&qs},{"unrelated",&other}}) {
        pair.second->setupLike(geometry,pair.first);bundle.distributed.add(State::scalarView(pair.first,0,geometry,*pair.second));
    }
    std::fill(a.values().begin(),a.values().end(),1);std::fill(b.values().begin(),b.values().end(),3);
    std::fill(as.values().begin(),as.values().end(),4);std::fill(bs.values().begin(),bs.values().end(),5);
    auto state=realizeState(system.state,{},bundle,compiled.stateViews);
    auto av=State::scalarView("a.stage",0,geometry,as),bv=State::scalarView("b.stage",0,geometry,bs);
    for(int stage=0;stage<2;++stage) {state.bindView("a",StateViewKind::Stage,av,stage);state.bindView("b",StateViewKind::Stage,bv,stage);}
    rejects("cannot bind",[&]{state.bindView("q",StateViewKind::Stage,av,0);});
    state.bindDerivedViews();
    const auto read=[&] {return state.view("q",StateViewKind::Physical).fields[0]->read(cell,0);};
    rejects("unready",[&]{read();});
    auto owner=std::make_shared<SolverAlgorithm::PatchWorkspace>();owner->residual.setup(geometry.MX(),geometry.MY(),geometry.MZ(),1);
    const auto resource=std::any_cast<UpdateContract>(compiled.steps[0].providerContract).residual;
    auto residual=ResourceView::bind(resource,owner,owner->residual,"patch");
    const auto* address=&owner->residual;
    Run::OpRegistry ops;ops.bind("derive.current","test.affine",[&]{state.updateView("q",StateViewKind::Physical);residual.write<Residual>()->source(1,1,1,0)=read();});
    CompiledSolvePlan plan;plan.compiledProgram=compiled;plan.root=compileMethodProgram(compiled,&providers);
    Run::PlanExecutor::validateBindings(plan,ops);Run::PlanExecutor::execute(plan,ops);
    require(read()==5 && owner->residual.source(1,1,1,0)==5,"Compiled current update/reference failed");
    require(residual.read<Residual>().get()==address,"Resource was copied");
    bundle.distributed.markModified("unrelated");require(read()==5,"Unrelated write invalidated q");
    a.values()[cell]=2;bundle.distributed.markModified("a");rejects("stale",[&]{read();});
    Run::PlanExecutor::execute(plan,ops);require(read()==7,"Current update read wrong inputs");
    rejects("read-only",[&]{state.view("q",StateViewKind::Physical).fields[0]->write(cell,0,99);});
    auto stagePlan=plan;stagePlan.root.kind=PlanNodeKind::StageLoop;stagePlan.root.repetitions=2;
    stagePlan.root.children.clear();
    for(const auto* op:{"open","update","ready","consume","rhs","close"}) {
        SolvePlanNode leaf;leaf.kind=PlanNodeKind::Update;leaf.operation=op;leaf.provider="test.affine";stagePlan.root.children.push_back(leaf);
    }
    const auto stageReader=state.stageReader("q",2);
    std::vector<double> observed;
    ops.bind("open","test.affine",[&](const Run::ExecutionContext& c){state.beginStage(c.stageIndex,c.stageCount,c.stageIndex*.1);});
    ops.bind("update","test.affine",[&](const Run::ExecutionContext& c){
        rejects("wrong Stage",[&]{state.view("q",StateViewKind::Stage,c.stageIndex).fields[0]->read(cell,0);});
        state.updateView("q",StateViewKind::Stage,c.stageIndex);
    });
    ops.bind("ready","test.affine",[&]{state.stageReady();});
    ops.bind("consume","test.affine",[&](const Run::ExecutionContext& c){
        const auto& view=state.view("q",StateViewKind::Stage,c.stageIndex);
        observed.push_back(stageReader(cell,0));require(read()==7,"Stage overwrote Current");
        if(c.stageIndex==0) {as.values()[cell]=6;bundle.distributed.markModified("a.stage");rejects("stale",[&]{view.fields[0]->read(cell,0);});}
    });
    ops.bind("rhs","test.affine",[&]{state.rhsReady();});ops.bind("close","test.affine",[&]{state.finishStage();});
    Run::PlanExecutor::validateBindings(stagePlan,ops);Run::PlanExecutor::execute(stagePlan,ops);
    require(observed==std::vector<double>{13,17},"Stage version/reference mismatch");
    rejects("wrong Stage",[&]{state.view("q",StateViewKind::Stage,0).fields[0]->read(cell,0);});
    state.beginStage(0,2,1);state.stageReady();rejects("unready",[&]{state.view("q",StateViewKind::Stage,0).fields[0]->read(cell,0);});state.rhsReady();state.finishStage();
    auto halo=compiled.stateViews;
    for(auto& v:halo)if(v.symbol=="q" && v.kind==StateViewKind::Physical)v.haloRequired=true;
    auto published=realizeState(system.state,{},bundle,halo);published.updateView("q",StateViewKind::Physical);
    rejects("unpublished",[&]{published.view("q",StateViewKind::Physical).fields[0]->read(cell,0);});
    for(int k=0;k<geometry.MZ();++k)for(int j=0;j<geometry.MY();++j)for(int i=0;i<geometry.MX();++i)
        qc(i,j,k)=qc(std::clamp(i,1,3),1,1); // Explicit serial zeroGradient replica copy.
    bundle.distributed.markSynchronized("q",1);require(published.view("q",StateViewKind::Physical).fields[0]->read(cell,0)==7,"Explicit publication failed");
    auto readonly=resource;readonly.writable=false;auto ro=ResourceView::bind(readonly,owner,owner->residual,"patch");
    rejects("access",[&]{ro.write<Residual>();});rejects("type",[&]{residual.read<FluxField>();});
    rejects("type/owner",[&]{ResourceView::bind(resource,owner,owner->residual,"wrong");});
    rejects("type/owner",[&]{ResourceView::bind(resource,owner,owner->convectiveFlux,"patch");});
    const auto size=owner->residual.xSize();auto* existing=&owner->residual;
    auto checked=ResourceView::bind(resource,owner,owner->residual,"patch",[existing,size] {
        if(existing->xSize()!=size)throw std::runtime_error("Execution contract: residual layout lifetime ended");
    });
    owner->residual.setup(1,1,1,1);rejects("layout lifetime",[&]{checked.read<Residual>();});
    owner.reset();rejects("expired",[&]{residual.read<Residual>();});
    StateRegistry missingBindings;missingBindings.add(scalar("a"));missingBindings.add(scalar("b"));missingBindings.add(derived());
    auto absent=demand("q","q.stage",StateViewKind::Stage,0);
    auto needsStage=realizeState(missingBindings,{},bundle,{absent});
    rejects("Unrequested STATE",[&]{needsStage.bindDerivedViews();});
    auto handle=state.view("q",StateViewKind::Physical).fields[0]->read;
    {auto another=realizeState(system.state,{},bundle,compiled.stateViews);another.updateView("q",StateViewKind::Physical);handle=another.at("q").fields[0]->read;}
    rejects("expired",[&]{handle(cell,0);});
    std::ofstream("result/affine.csv")<<"current,stage0,stage1\n7,13,17\n";
}

void chained() {
    Field geometry;geometry.setupGeometry(1,1,1,0);
    State::StateBundle bundle;bundle.patches={&geometry};
    ScalarField a,b,q,r;
    for (auto item:std::vector<std::pair<std::string,ScalarField*>>{{"a",&a},{"b",&b},{"q",&q},{"r",&r}}) {
        item.second->setupLike(geometry,item.first);bundle.distributed.add(State::scalarView(item.first,0,geometry,*item.second));
    }
    a.values()[0]=1;b.values()[0]=3;
    auto second=derived("r");second.dependencies={"q","b"};
    auto binding=std::make_shared<StateProducer>(*second.producer);
    binding->expression=FormulaExpr::constantValue(2)*scalar("q").expression()+scalar("b").expression();second.producer=binding;
    StateRegistry symbols;symbols.add(second);symbols.add(derived());symbols.add(scalar("b"));symbols.add(scalar("a"));
    auto state=realizeState(symbols,{},bundle);
    rejects("unready",[&]{state.updateView("r",StateViewKind::Physical);});
    state.updateView("q",StateViewKind::Physical);state.updateView("r",StateViewKind::Physical);
    const auto read=state.at("r").fields[0]->read;
    require(read(0,0)==13,"Reverse registration order changed dependent bindings");
    a.values()[0]=2;bundle.distributed.markModified("a");rejects("stale",[&]{read(0,0);});
    state.updateView("q",StateViewKind::Physical);rejects("stale",[&]{read(0,0);});
    state.updateView("r",StateViewKind::Physical);require(read(0,0)==17,"Explicit dependency refresh failed");
    bundle.distributed.markModified("r");rejects("overwritten",[&]{read(0,0);});
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
void velocity() {
    Field field;field.setup(1,1,1,0,5);auto model=std::make_shared<CountingModel>();field.setStateModel(model);
    field(0,0,0,0)=2;field(0,0,0,1)=6;field(0,0,0,2)=field(0,0,0,3)=0;field(0,0,0,4)=250009;
    StateRegistry symbols;BuiltinStateCatalog catalog;
    for(const auto* id:{"rho","rhoU","rhoE"})symbols.add(catalog.solution(id));symbols.add(catalog.at("U"));
    State::StateBundle bundle;bundle.patches={&field};bundle.distributed.add(State::conservativeView("conservative",0,field));
    auto state=realizeState(symbols,{},bundle);
    require(model->calls==0,"Binding eagerly computed EOS");
    const auto reader=state.at("U").fields[0]->read;
    require(reader(0,0)==3 && reader(0,1)==0 && model->calls==1,"Existing cache miss count changed");
    require(field.thermodynamicState(0,0,0).velocity[0]==3 && model->calls==1,"Original cache not reused");
    state.beginStage(0,2,.1);field(0,0,0,1)=10;field.invalidateThermodynamicCache();state.stageReady();
    require(model->calls==1 && reader(0,0)==5 && model->calls==2,"Stage publication closure timing changed");state.rhsReady();state.finishStage();
    state.beginStage(1,2,.2);field(0,0,0,0)=4;field.invalidateThermodynamicCache();state.stageReady();
    require(reader(0,0)==2.5 && model->calls==3,"Second stage invalidation changed");state.rhsReady();state.finishStage();
    rejects("read-only",[&]{state.at("U").fields[0]->write(0,0,0);});
    auto working=demand("U","working.U",StateViewKind::Working);working.owner=StateViewOwner::CompilerWorkspace;
    auto snapshot=realizeState(symbols,{},bundle,{working});
    field(0,0,0,1)=12;field.invalidateThermodynamicCache();
    require(snapshot.view("U",StateViewKind::Working).fields[0]->read(0,0)==2.5 && snapshot.at("U").fields[0]->read(0,0)==3,"Working snapshot/physical alias lost identity");
    std::cout<<"Velocity: original cache miss timing preserved; Current plus two stage publications; Working snapshot unchanged\n";
}
void negatives() {
    const auto validate=[](const StateRegistry& s,std::vector<CompiledStateView> views={}){compileDerivedViews(s,views);};
    Field first,second;first.setupGeometry(1,1,1,0);second.setupGeometry(1,1,1,0);
    ScalarField a0,a1;a0.setupLike(first,"a");a1.setupLike(second,"a");State::DistributedFieldRegistry distributed;
    distributed.add(State::scalarView("a",0,first,a0));distributed.add(State::scalarView("a",1,second,a1));
    distributed.markModified("a",1);require(distributed.generation("a",0)==1 && distributed.generation("a",1)==2,"Per-patch publication dirtied an unrelated patch");
    distributed.markModified("a");require(distributed.generation("a",0)==2 && distributed.generation("a",1)==3,"Existing all-patch dirty semantics changed");
    StateRegistry s;s.add(scalar("a"));s.add(scalar("b"));s.add(derived());validate(s);
    rejects("Duplicate",[&]{s.add(derived());});
    StateRegistry missing;missing.add(scalar("a"));missing.add(derived());rejects("missing derived dependency",[&]{validate(missing);});
    StateRegistry cycle;auto a=derived("a");a.dependencies={"q","b"};auto p=std::make_shared<StateProducer>(*a.producer);
    p->expression=scalar("q").expression()+scalar("b").expression();a.producer=p;cycle.add(a);cycle.add(scalar("b"));cycle.add(derived());rejects("cycle",[&]{validate(cycle);});
    StateRegistry wrong;wrong.add(scalar("a"));wrong.add(scalar("b"));auto q=derived();q.signature=Math::signature(0,{});wrong.add(q);rejects("derived output",[&]{validate(wrong);});
    StateRegistry typo;auto bad=scalar("a");bad.signature=Math::signature(0,{});typo.add(bad);typo.add(scalar("b"));typo.add(derived());rejects("STATE reference",[&]{validate(typo);});
    auto view=demand("q","q.old",StateViewKind::OldTime);rejects("Unsupported",[&]{validate(s,{view});});
    view=demand("q","q.stage",StateViewKind::Stage,0);view.components=3;rejects("layout",[&]{validate(s,{view});});
    view=demand("q","q",StateViewKind::Physical);view.producer=std::make_shared<StateProducer>(*s.at("q").producer);rejects("conflicting",[&]{validate(s,{view});});
    StateRegistry noImplementation;noImplementation.add(scalar("a"));noImplementation.add(scalar("b"));
    auto unavailable=derived();auto absent=std::make_shared<StateProducer>(*unavailable.producer);absent->bind={};unavailable.producer=absent;noImplementation.add(unavailable);
    rejects("Unsupported missing STATE producer",[&]{validate(noImplementation);});
    view=demand("q","q.stage",StateViewKind::Stage,-1);rejects("Invalid STATE view stage",[&]{validate(s,{view});});
    CompiledTarget target;target.symbol="q";rejects("read-only",[&]{realizeTarget(s,target);});
}
}
int main(){std::filesystem::create_directories("result");extension();chained();velocity();negatives();std::cout<<"A3 STATE/resource checks passed: "<<checks<<'\n';}
