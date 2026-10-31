#include "solver/system/SF_stateRealizer.h"
#include "solver/system/SF_numericalSelection.h"
#include "core/field/SF_scalarField.h"
#include <iostream>
#include <filesystem>
#include "infrastructure/io/case/SF_case.h"
using namespace SF;using namespace SF::System;
namespace {
void check(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F>void rejects(const char* expected,F action){try{action();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,e.what());return;}throw std::runtime_error("Expected rejection");}
void history(){
 StateRegistry registry;StateSymbol a;a.id=a.name=a.storageKey="a";a.signature=Math::signature(0,{});a.storageBinding=StorageBinding::NamedDistributed;registry.add(a);
 Field geometry;geometry.setupGeometry(3,1,1,1);const int cell=geometry.getIdx(1,1,1);State::StateBundle bundle;bundle.patches={&geometry};bundle.step=7;bundle.time=3.;
 ScalarField physical,old;physical.setupLike(geometry,"a");old.setupLike(geometry,"original.q0");physical.values()[cell]=2;old.values()[cell]=1;
 bundle.distributed.add(State::scalarView("a",0,geometry,physical));
 CompiledStateView demand;demand.symbol="a";demand.storage="original.q0";demand.kind=StateViewKind::OldTime;demand.owner=StateViewOwner::NumericalProvider;
 auto state=realizeState(registry,{},bundle,{demand});auto borrowed=State::scalarView("original.q0",0,geometry,old);state.bindView("a",StateViewKind::OldTime,borrowed);
 auto read=state.oldTimeReader("a");rejects("not published",[&]{read(cell,0);});state.publishOldTime();check(read(cell,0)==1,"OldTime fell back to physical");
 physical.values()[cell]=3;check(read(cell,0)==1,"Physical mutation overwrote q0");rejects("repeated",[&]{state.publishOldTime();});
 ++bundle.step;bundle.time=3.1;rejects("not published",[&]{read(cell,0);});state.beginOldTime();old.values()[cell]=3;state.publishOldTime();check(read(cell,0)==3,"Next-step history not rebound");
 check(physical.values()[cell]==3 && bundle.step==8 && bundle.time==3.1,"Reader modified state/clock");
 auto moved=std::move(state);check(read(cell,0)==3,"Move broke a valid borrowed history reader");
 std::function<double(int,int)> expired;
 {auto temporary=realizeState(registry,{},bundle,{demand});temporary.bindView("a",StateViewKind::OldTime,borrowed);temporary.publishOldTime();expired=temporary.oldTimeReader("a");}
 rejects("Expired OldTime",[&]{expired(cell,0);});
}
void choices(){
 EquationCall call;call.equation="E";call.occurrence="S";
 NumericalBinding defaults;defaults.equation="E";defaults.method="method";defaults.parameters={{"a",1},{"b",2}};defaults.origin={OriginKind::BuiltinDefault,"defaults"};
 NumericalBinding override;override.equation="*";override.occurrence="S";override.parameters={{"b",4}};override.origin={OriginKind::User,"authored"};override.explicitChoice=true;
 auto result=resolveNumerics(call,{defaults,override});check(result.value.method=="method" && result.value.parameters.at("a")==1 && result.value.parameters.at("b")==4,"Partial override discarded defaults");
 rejects("Conflicting",[&]{resolveNumerics(call,{defaults,override,override});});rejects("Ambiguous",[&]{resolveNumerics(call,{defaults,defaults});});auto ambiguous=defaults;ambiguous.method="other";rejects("Ambiguous",[&]{resolveNumerics(call,{defaults,ambiguous});});
 override.method="explicit";check(resolveNumerics(call,{defaults,ambiguous,override}).value.method=="explicit","Explicit choice failed to cover defaults");
 rejects("whole block",[&]{resolveNumerics(call,{defaults},true);});
}
}
int main(int argc,char** argv){try{
 if(argc==3){
  auto input=CaseIO::read(argv[1]);CaseIO::write(input,argv[2]);
  std::filesystem::copy_file(std::filesystem::path(argv[1])/"mesh/mesh.sfm",std::filesystem::path(argv[2])/"mesh/mesh.sfm",std::filesystem::copy_options::overwrite_existing);
  auto reread=CaseIO::read(argv[2]);check(reread.formatVersion==input.formatVersion,"Contract version lost during roundtrip");
  for(const auto& original:input.objects)if(original.category=="solver") {
   const auto found=std::find_if(reread.objects.begin(),reread.objects.end(),[&](const auto& value){return value.category==original.category && value.name==original.name;});
   check(found!=reread.objects.end() && found->parameters==original.parameters,"Four-axis roundtrip changed source declarations");
  }
 }
 history();choices();std::cout<<"OldTime borrowed snapshot epochs and WHICH coverage passed\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
