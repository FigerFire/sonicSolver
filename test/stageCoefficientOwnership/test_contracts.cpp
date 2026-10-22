#include "test/typedSymbols.h"
#define SF_STAGE_COEFFICIENT_LIBRARY
#include "test_evidence.cpp"
#undef SF_STAGE_COEFFICIENT_LIBRARY
#include "solver/system/SF_executionContract.h"
#include <iostream>
namespace {
using namespace SF::System;
void failure(const std::function<void()>& run,const std::string& expected) {
    bool rejected=false;
    try{run();}catch(const std::runtime_error& e){rejected=std::string(e.what()).find(expected)!=std::string::npos;std::cout<<e.what()<<'\n';}
    require(rejected,"Missing/wrong coefficient diagnostic");
}
StateRegistry symbols() {
    StateRegistry state;
    for(const auto* id:{"a","b"}){StateSymbol s;s.id=s.name=s.storageKey=id;s.role=StateRole::Derived;state.add(SF::Test::state(s));}
    return state;
}
CompiledExecutionProgram declarations(StateVersion provided=StateVersion::Stage) {
    CompiledExecutionProgram p;
    CompiledEquationCall pub;pub.source.occurrence="publisher";pub.backendProvider="test.coefficients";
    pub.coefficientProvides={"coefficient-set"};pub.coefficientStates={"a","b"};pub.coefficientVersion=provided;
    if(provided==StateVersion::Stage){pub.temporalResidual=true;pub.coefficientPrepare="prepare.coefficients";pub.coefficientPublish="publish.coefficients";}
    else pub.stateEffects={{"a",true},{"b",true},{"coefficient-set",true}};
    CompiledEquationCall read;read.source.occurrence="consumer";read.coefficientUses={{"a",{StateVersion::Stage,StateVersion::Frozen}}};
    p.steps={pub,read};return p;
}
void staticContracts() {
    const auto state=symbols();
    for(const auto treatment:{StateVersion::Stage,StateVersion::Frozen}) {
        auto p=declarations(treatment);bindCoefficientContracts(state,p);validateCoefficientContracts(state,p);
        const auto& use=p.steps[1].coefficientUses[0];require(use.version==treatment && use.producer=="publisher" && use.resource=="coefficient-set","Selected numerical treatment did not freeze its publisher");
    }
    require(state.size()==2 && state.at("a").storageKey=="a" && state.at("a").dependencies.empty(),"Binding rewrote base STATE");
    auto p=declarations();p.steps.erase(p.steps.begin());failure([&]{bindCoefficientContracts(state,p);},"Missing coefficient publisher");
    p=declarations();auto another=p.steps[0];another.source.occurrence="competing";p.steps.push_back(another);failure([&]{bindCoefficientContracts(state,p);},"Conflicting coefficient publishers");
    p=declarations();another=p.steps[0];another.temporalOwner="publisher";another.coefficientProvides={"another-bundle"};another.coefficientStates={"b"};p.steps[0].coefficientStates={"a"};p.steps.push_back(another);failure([&]{bindCoefficientContracts(state,p);},"Conflicting coefficient publishers");
    p=declarations(StateVersion::Frozen);p.steps[1].coefficientUses[0].acceptedVersions={StateVersion::Stage};failure([&]{bindCoefficientContracts(state,p);},"Incompatible coefficient treatment");
    p=declarations();p.steps[1].coefficientUses[0].acceptedVersions={StateVersion::Current};failure([&]{bindCoefficientContracts(state,p);},"Incompatible coefficient treatment");
    p=declarations();p.steps[0].coefficientStates={"b"};failure([&]{bindCoefficientContracts(state,p);},"Missing coefficient publisher");
    p=declarations();p.steps[0].coefficientStates.push_back("a");failure([&]{bindCoefficientContracts(state,p);},"Duplicate coefficient state");
    p=declarations();p.steps[1].stateUses={{"a",StateVersion::Stage,true}};failure([&]{bindCoefficientContracts(state,p);},"halo incompatible");
    p=declarations();bindCoefficientContracts(state,p);p.steps[1].coefficientUses.clear();failure([&]{validateCoefficientContracts(state,p);},"lacks frozen producer binding");
    p=declarations();bindCoefficientContracts(state,p);p.steps[1].coefficientUses[0].producer="wrong";failure([&]{validateCoefficientContracts(state,p);},"Missing/incompatible frozen coefficient binding");
    p=declarations();bindCoefficientContracts(state,p);p.steps[1].stateUses[0].version=StateVersion::Frozen;failure([&]{validateCoefficientContracts(state,p);},"stateUses/version/halo");
    p=declarations();bindCoefficientContracts(state,p);TemporalParticipant participant;participant.identity="publisher";participant.coefficientProvides=p.steps[0].coefficientProvides;
    participant.coefficientStates={"a"};participant.coefficientPrepare=p.steps[0].coefficientPrepare;participant.coefficientPublish=p.steps[0].coefficientPublish;p.temporalParticipants={participant};
    failure([&]{validateCoefficientContracts(state,p);},"publication differs from compiled provider");
}
}
int main(){
    staticContracts();
    const auto base=std::filesystem::current_path();
    for(const auto* fault:{"unpublished","staleStage","staleStep","wrongProvidedSet"}) {
        auto path=base/fault;std::filesystem::create_directories(path);std::filesystem::current_path(path);
        sst(FDM::TimeRecipeId::ClassicalRK4,true,true,.0005,.001,0,fault);std::filesystem::current_path(base);
    }
    std::cout<<"coefficient static binding and runtime epoch contracts passed\n";
}
