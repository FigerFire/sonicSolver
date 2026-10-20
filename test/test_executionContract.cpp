#include "solver/system/SF_executionContract.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_builtinState.h"
#include "solver/system/SF_formulaStorage.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include <cmath>
#include <functional>
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
void fails(const std::function<void()>& test,const std::string& reason) {
    try {test();} catch(const std::runtime_error& error) {
        require(std::string(error.what()).find(reason)!=std::string::npos,error.what());return;
    }
    throw std::runtime_error("Expected rejection: "+reason);
}
StateSymbol symbol(std::string id,std::string storage,int offset=0,int components=1) {
    StateSymbol result;result.id=id;result.name=id;result.storageKey=storage;
    result.storageBinding=StorageBinding::PackedDistributed;result.componentOffset=offset;result.components=components;
    result.shape=components==1?ValueShape::Scalar:ValueShape::Vector;return result;
}
ExecutionScope call(std::string id) {
    ExecutionScope result;result.kind=ExecutionKind::EquationCall;result.id=id;result.step={id,{"m"},id};return result;
}
CompiledEquationCall method(std::string id,std::vector<StateUse> reads,std::vector<StateEffect> writes) {
    CompiledEquationCall result;result.source=call(id).step;result.stateUses=std::move(reads);result.stateEffects=std::move(writes);return result;
}
CompiledExecutionProgram program(std::initializer_list<CompiledEquationCall> calls) {
    CompiledExecutionProgram result;
    for(auto& compiled:calls) {result.root.children.push_back(call(compiled.source.occurrence));result.steps.push_back(compiled);}
    return result;
}
}
int main() {
    StateRegistry state;
    state.add(symbol("m","Q",1,3));state.add(symbol("e","Q",4));state.add(symbol("packed","Q",0,5));
    auto prior=symbol("laggedM","laggedM");prior.availableVersion=StateVersion::Lagged;state.add(prior);
    auto u=symbol("U","u");u.role=StateRole::Derived;u.evaluation=StateEvaluation::Materialized;u.dependencies={"m"};state.add(u);
    auto lazy=symbol("instantU","instant");lazy.evaluation=StateEvaluation::Lazy;lazy.dependencies={"m"};state.add(lazy);
    auto materialized=symbol("closure","closure");materialized.evaluation=StateEvaluation::Materialized;materialized.dependencies={"instantU"};state.add(materialized);
    const auto write=method("writer",{},{{"m"}}),refresh=method("refresh",{},{{"U",true}}),consumer=method("consumer",{{"U"}},{});
    fails([&]{validateDataFlow(state,program({write,consumer}));},"invalidating writer 'writer'");
    validateDataFlow(state,program({write,refresh,consumer}));
    validateDataFlow(state,program({method("unrelated",{},{{"e"}}),consumer}));
    fails([&]{validateDataFlow(state,program({method("packedWrite",{},{{"packed"}}),consumer}));},"packedWrite");
    fails([&]{validateDataFlow(state,program({write,method("indirectRead",{{"closure"}}, {})}));},"invalidating writer 'writer'");
    validateDataFlow(state,program({write,method("lazyRead",{{"instantU"}}, {})}));
    fails([&]{validateDataFlow(state,program({write,method("haloRead",{{"m",StateVersion::Current,true}}, {})}));},"halo publication");
    validateDataFlow(state,program({method("synchronizedWrite",{},{{"m",false,true}}),method("haloRead",{{"m",StateVersion::Current,true}}, {})}));
    auto frozen=program({method("snapshot",{},{{"frozenM",true}}),write,method("frozenRead",{{"frozenM",StateVersion::Frozen}}, {})});
    validateDataFlow(state,frozen);
    validateDataFlow(state,program({write,method("laggedRead",{{"laggedM",StateVersion::Lagged}}, {})}));
    fails([&]{validateDataFlow(state,program({method("falseLag",{{"m",StateVersion::Lagged}}, {})}));},"Unavailable lagged");
    fails([&]{validateDataFlow(state,program({method("oldRead",{{"m",StateVersion::OldTime,false,"oldM"}}, {})}));},"Unavailable old-time");
    validateDataFlow(state,program({method("oldSnapshot",{},{{"oldM",true}}),write,method("oldRead",{{"m",StateVersion::OldTime,false,"oldM"}}, {})}));
    fails([&]{validateDataFlow(state,program({method("untrackedWriter",{},{{"m",false,false,false}}),method("lazyRead",{{"instantU"}}, {})}));},"Lazy STATE invalidation missing");
    auto loop=program({refresh,consumer,write});loop.root.kind=ExecutionKind::Loop;loop.root.repetitions=3;
    // Refresh each loop iteration precedes the consumer; no requirement to use newest m after its last write.
    validateDataFlow(state,loop);
    auto invalidLoop=program({consumer,write});invalidLoop.root.kind=ExecutionKind::Loop;invalidLoop.root.repetitions=3;
    fails([&]{validateDataFlow(state,invalidLoop);},"Stale STATE");
    auto acrossStep=program({consumer,write});fails([&]{validateDataFlow(state,acrossStep);},"Stale STATE");
    auto frozenLoop=program({method("snapshot",{},{{"frozenM",true}}),method("frozenRead",{{"frozenM",StateVersion::Frozen}}, {})});
    frozenLoop.root.kind=ExecutionKind::Loop;frozenLoop.root.repetitions=3;
    fails([&]{validateDataFlow(state,frozenLoop);},"Frozen STATE overwritten");
    auto legalOuter=frozenLoop;
    legalOuter.steps[1].stateUses[0].perIteration=true;
    validateDataFlow(state,legalOuter);
    auto escape=program({method("snapshot",{},{{"frozenM",true}}),method("frozenRead",{{"frozenM",StateVersion::Frozen}}, {})});
    ExecutionScope scope;scope.kind=ExecutionKind::Loop;scope.id="local";scope.children={escape.root.children[0]};
    escape.root.children[0]=scope;fails([&]{validateDataFlow(state,escape);},"Unavailable frozen");
    auto persistentLoop=escape;
    persistentLoop.steps[0].stepWorkspaces={"frozenM"};
    validateDataFlow(state,persistentLoop);
    StateRegistry cycle;auto a=symbol("a","a"),b=symbol("b","b");a.dependencies={"b"};b.dependencies={"a"};cycle.add(a);cycle.add(b);
    fails([&]{validateDataFlow(cycle,{});},"a -> b -> a");
    ExecutionProgram ordering;ordering.root.children={call("A"),call("B")};ordering.requirements={{{},"B",{"A"},{}}};
    validatePlacement(ordering,ordering.root);
    ordering.requirements.push_back({{},"A",{"B"},{}});
    fails([&]{validatePlacement(ordering,ordering.root);},"dependency cycle");
    ordering.requirements.pop_back();std::swap(ordering.root.children[0],ordering.root.children[1]);
    fails([&]{validatePlacement(ordering,ordering.root);},"placement violation");
    ExecutionProgram defaults;defaults.root.children={call("flow"),call("finish"),call("module")};
    defaults.requirements={{{},"module",{"flow"},{"finish"}}};
    composeDefaultPlacement(defaults,defaults.root);validatePlacement(defaults,defaults.root);
    require(defaults.root.children[1].id=="module","Relative placement depended on absolute order.");
    auto ambiguous=defaults;ambiguous.root.children.push_back(call("other"));
    ambiguous.requirements.push_back({{},"other",{"flow"},{"finish"}});
    fails([&]{composeDefaultPlacement(ambiguous,ambiguous.root);},"Ambiguous default HOW");
    auto severalSlots=defaults;
    severalSlots.root.children.insert(severalSlots.root.children.begin()+1,call("intermediate"));
    fails([&]{composeDefaultPlacement(severalSlots,severalSlots.root);},"Ambiguous default HOW");
    auto wrongScope=defaults;wrongScope.requirements[0].scope="absent";
    fails([&]{validatePlacement(wrongScope,wrongScope.root);},"Missing placement scope");
    auto cyclic=defaults;cyclic.requirements.push_back({{},"flow",{"finish"},{}});
    fails([&]{composeDefaultPlacement(cyclic,cyclic.root);},"dependency cycle");
    auto capability=program({method("providerA",{},{}),method("providerB",{}, {})});
    capability.steps[0].capabilities={"velocity.corrected"};
    capability.steps[1].capabilityRequirements={{"velocity.corrected","boundary.response","explicit response is missing"}};
    require(missingCapability(capability,{},capability.steps[1]).find("boundary.response")!=std::string::npos,"Missing capability was hidden.");
    capability.steps[0].capabilities.push_back("boundary.response");require(missingCapability(capability,{},capability.steps[1]).empty(),"Provider name influenced capability validation.");
    // Algebraic unknowns may form a coupled block: only STATE derivation graphs are checked for cycles.
    StateRegistry algebraic;algebraic.add(symbol("a","a"));algebraic.add(symbol("b","b"));validateDataFlow(algebraic,{});
    auto stage=program({method("stageRead",{{"m",StateVersion::Stage}}, {})});
    fails([&]{validateDataFlow(state,stage);},"outside declared StageLoop");
    stage.root.kind=ExecutionKind::StageLoop;stage.root.repetitions=2;validateDataFlow(state,stage);
    SF::Field field;field.setup(1,1,1,1,5);
    field.setStateModel(SF::Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0.0,0.72));
    field(1,1,1,0)=2;field(1,1,1,1)=6;field(1,1,1,4)=30;
    SF::State::StateBundle bundle;bundle.patches={&field};bundle.registerConservativeState();
    StateRegistry physical;const BuiltinStateCatalog catalog;
    for (const auto* id:{"rho","rhoU","rhoE"}) physical.add(catalog.solution(id));
    catalog.require(physical,"U");catalog.require(physical,"p");catalog.require(physical,"T");
    const auto realized=realizeState(physical,RuntimeRequirements{},bundle);
    const auto cell=field.getIdx(1,1,1);
    const auto velocity=realized.at("U").fields[0];
    require(velocity->read(cell,0)==3.0,"Initial lazy velocity does not alias packed Q.");
    auto values=bindFormulaCellValues(realized,CompiledTarget{"rho"},
        [](const std::string&,int,int,int,int){return 0.0;});
    values.write("rho",0,0,4);
    require(field(1,1,1,0)==4 && velocity->read(cell,0)==1.5,
        "Packed STATE write failed to invalidate the actual lazy thermodynamic view.");
    const auto temperature=realized.at("T").fields[0]->read(cell,0);
    auto energy=bindFormulaCellValues(realized,CompiledTarget{"rhoE"},
        [](const std::string&,int,int,int,int){return 0.0;});
    energy.write("rhoE",0,0,40);
    require(realized.at("T").fields[0]->read(cell,0)>temperature && velocity->read(cell,0)==1.5,
        "Total-energy write did not refresh temperature or altered independent velocity.");

    auto scoped=capability;
    scoped.steps[1].capabilityRequirements[0]={"velocity.corrected","boundary.response","wrong wall binding","momentum","rhoU","wallA"};
    scoped.steps[0].boundCapabilities={{"boundary.response","momentum","rhoU","wallB"}};
    require(!missingCapability(scoped,{},scoped.steps[1]).empty(),"Wrong immersed wall satisfied the consumer.");
    scoped.steps[0].boundCapabilities[0].boundary="wallA";
    scoped.steps[0].boundCapabilities[0].fluidPort="rhoU.other";
    require(!missingCapability(scoped,{},scoped.steps[1]).empty(),"Wrong fluid port satisfied the consumer.");
    scoped.steps[0].boundCapabilities[0].fluidPort="rhoU";
    scoped.steps[0].boundCapabilities[0].equation="energy";
    require(!missingCapability(scoped,{},scoped.steps[1]).empty(),"Wrong equation satisfied the consumer.");
    scoped.steps[0].boundCapabilities[0].equation="momentum";
    require(missingCapability(scoped,{},scoped.steps[1]).empty(),"Matching wall binding was rejected.");

}
