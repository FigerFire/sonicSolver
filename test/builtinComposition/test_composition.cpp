#include "test/typedSymbols.h"
#include "app/application/model/SF_model.h"
#include "app/application/model/SF_output.h"
#include "app/application/SF_inspection.h"
#include "app/application/SF_environment.h"
#include "app/application/SF_application.h"
#include "app/application/execution/SF_execution.h"
#include "infrastructure/io/case/SF_case.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemPrinter.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_builtinState.h"
#include "app/application/SF_hostCapabilities.h"
#include "models/physics/gravity/SF_accelerationProvider.h"
#include <filesystem>
#include <fstream>
#include <iostream>
using namespace SF;
using namespace SF::System;
namespace {
void check(bool value,const std::string& why) {if(!value)throw std::runtime_error(why);}
void prepare(const std::string& source,const std::string& output,const std::string& kind) {
    auto m=CaseIO::read(source);m.caseDir=output;m.name=kind;
    if(kind.rfind("negative.",0)==0) {
        const auto find=[&](const std::string& type)->Model::ObjectDescriptor& {for(auto& o:m.objects)if(o.type==type)return o;throw std::runtime_error("Missing mutation object: "+type);};
        if(kind=="negative.ambiguity")m.objects.push_back({"otherComposition","solver","compositionPreset",{},{{"preset","singleFluidLaminar"}}});
        else if(kind=="negative.conflict") {auto second=find("turbulence");second.name="otherTurbulence";second.type="Turbulence";second.parameters["enabled"]=false;m.objects.push_back(second);}
        else if(kind=="negative.time")m.numerics["time"]["default"]="BDF2";
        else if(kind=="negative.initial") {for(auto& f:m.fields)if(f.symbol=="k")f.initial=nullptr;}
        else if(kind=="negative.boundary") {for(auto& f:m.fields)if(f.symbol=="omega")f.boundaries=Model::Parameters::object();}
        else if(kind=="negative.field")m.fields.erase(std::remove_if(m.fields.begin(),m.fields.end(),[](const auto& f){return f.symbol=="k";}),m.fields.end());
        else if(kind=="negative.output") {for(auto& f:m.fields)if(f.symbol=="omega"){f.output=true;f.outputDeclared=true;}}
        else if(kind=="negative.model")find("compositionPreset").parameters["preset"]="unknown";
        else if(kind=="negative.enable")find("gravity").parameters["enable"]=true;
        else if(kind=="negative.sourceAlias")find("gravity").parameters["g"]=Model::Parameters::array({1,0,0});
        else if(kind=="negative.source")find("gravity").parameters.erase("value");
        else if(kind=="negative.mpi")m.runtime["parallel"]={{"enabled",true},{"nProcs",2}};
        else if(kind=="negative.binding")m.objects.push_back({"nativeHow","solver","executionProgram",{},{{"root",{{"kind","Sequence"},{"children",Model::Parameters::array({{{"kind","Commit"},{"id","commit"}}})}}}}});
        else throw std::runtime_error("Unknown mutation "+kind);
        CaseIO::write(m,output);return;
    }
    m.mesh["files"]=Model::Parameters::array({"mesh/mesh.sfm"});
    std::ofstream mesh(output+"/mesh/mesh.sfm");mesh<<"#Information\n101 21 2 4\n\n#Point\n";
    for(int k=0;k<2;++k)for(int j=0;j<21;++j)for(int i=0;i<101;++i)mesh<<i<<" "<<j<<" "<<.1*k<<"\n";
    for(const auto& set:{"nativeLeft","nativeRight","nativeSides","frontBack"}) {
        mesh<<"\n#"<<set<<"\n";
        for(int k=0;k<2;++k)for(int j=0;j<21;++j)for(int i=0;i<101;++i)
            if((std::string(set)=="nativeLeft" && i==0) || (std::string(set)=="nativeRight" && i==100)
                || (std::string(set)=="nativeSides" && (j==0 || j==20)) || std::string(set)=="frontBack")mesh<<(k*21+j)*101+i<<" ";
        mesh<<"\n";
    }
    mesh.close();
    for(auto& field:m.fields) {
        field.sets=Model::Parameters::object();field.initial=0.;
        if(field.symbol=="rho")field.initial=1.2;
        if(field.symbol=="p")field.initial=100000.;
        if(field.symbol=="U")field.initial=Model::Parameters::array({.3,0,0});
        if(field.symbol=="k")field.initial=.01;
        if(field.symbol=="omega")field.initial=1.;
        field.boundaries={{"nativeRight",{{"type","zeroGradient"}}},{"nativeSides",{{"type","zeroGradient"}}},{"frontBack",{{"type","empty"}}}};
        if(field.symbol=="k" || field.symbol=="omega")field.boundaries["nativeLeft"]={{"type","zeroGradient"}};
        if(field.output)field.outputDeclared=false;
    }
    for(auto& object:m.objects)if(object.type=="thermoDynamics")object.parameters={
        {"thermoDynamics",{{"equationOfState","perfectGas"},{"thermo","hConst"},{"transport","const"}}},
        {"properties",{{"equationOfState",{{"gamma",1.4},{"R",287.05}}},{"transport",{{"mu",1e-5},{"Pr",.72}}}}}};

    const bool preset=kind.find("Preset")!=std::string::npos;
    m.objects.erase(std::remove_if(m.objects.begin(),m.objects.end(),[&](const auto& o){
        return o.type=="turbulence" || o.type=="gravity" || o.type=="MRF" || o.type=="IBM"
            || (preset && (o.type=="stateRegistry" || o.type=="equationRegistry" || o.type=="algorithmRegistry" || o.type=="equationSystem"));}),m.objects.end());
    m.runtime["createMesh"]=false;m.runtime["endTime"]=.001;m.runtime["maxDeltaT"]=.0001;
    m.runtime["writeControl"]="timeStep";m.runtime["writeInterval"]=1;m.runtime["writeInitial"]=true;
    m.runtime["parallel"]={{"enabled",false},{"nProcs",1},{"split","auto"},{"haloTolerance",1e-5}};
    m.numerics={{"time",{{"default","classicalRK4"}}},{"terms",{{"convection","weno5Rusanov"},{"diffusion","central2Explicit"}}}};
    const auto add=[&](std::string id,std::string type,Model::Parameters p,std::string category="models") {m.objects.push_back({id,category,type,{},p});};
    if(preset) {
        add("composition","compositionPreset",{{"preset",kind=="laminarPreset"?"singleFluidLaminar":"singleFluidSST"}},"solver");
        m.numerics=Model::Parameters::object();
    }
    if(kind=="laminarAbsent" || kind=="laminarExplicit")m.numerics["time"]["default"]="forwardEuler";
    if(kind=="laminarExplicit")add("turbulence","turbulence",{{"simulationType","laminar"}});
    if(kind=="sstExplicit" || kind=="sstDefaultOn") {
        Model::Parameters t={{"simulationType","RAS"},{"RAS",{{"model","kOmegaSST"}}},{"transportMethod","stageConservative"}};
        if(kind=="sstExplicit")t["enabled"]=true;
        add("turbulence","turbulence",t);
    }
    if(kind=="resolvedDefaultsPreset") {
        add("otherComposition","compositionPreset",{{"preset","singleFluidLaminar"}},"solver");
        add("turbulence","turbulence",{{"simulationType","RAS"},{"RAS",{{"model","kOmegaSST"}}},{"transportMethod","stageConservative"}});
        m.numerics["time"]["default"]="forwardEuler";
    }
    if(kind=="laminarOverridePreset")add("turbulence","turbulence",{{"simulationType","laminar"}});
    if(kind=="sstDisabledPreset")add("turbulence","turbulence",{{"enabled",false}});
    if(kind=="overridePreset")m.numerics={{"time",{{"default","forwardEuler"}}},{"terms",{{"convection","weno7Rusanov"}}}};
    if(kind=="gravityOn" || kind=="gravityOff")add("gravity","gravity",{{"enabled",kind=="gravityOn"},{"value",Model::Parameters::array({9.8,0,0})}});
    if(kind=="gravityZone")add("gravity","gravity",{{"all",{{"value",Model::Parameters::array({9.8,0,0})}}}});
    if(kind=="gravityDefaultOn")add("gravity","gravity",{{"value",Model::Parameters::array({9.8,0,0})}});
    CaseIO::write(m,output);
}
ResolvedSimulationSystem manual(const CaseConfig& config) {
    using E=FormulaExpr;
    BuildRequest request;request.includeDefaultFluidPreset=false;request.composition=config.composition;
    request.capabilities=Application::detectBuildCapabilities();
    request.termProviders=Physics::Gravity::termProviders(config.solver.sources.gravity);
    SystemContribution c;c.recordContribution("test.manual","independent explicit contribution");
    BuiltinStateCatalog catalog;
    for(const auto* id:{"rho","rhoU","rhoE"})c.addState(catalog.solution(id));
    for(const auto* id:{"U","p","T","mu","conductivity"})c.requireState(id);
    const auto ddt=[](const char* value,const char* id){return E::op("ddt",{SF::Test::symbol(value)},id);};
    const auto div=[](const char* value,const char* id){return E::op("div",{SF::Test::symbol(value)},id);};
    c.addEquation({"continuity",E::add(ddt("rho","mass.time"),div("massFlux","mass.convection")),E::constantValue(0)});
    c.addEquation({"momentum",E::add(E::add(ddt("rhoU","momentum.time"),div("momentumFlux","momentum.convection")),E::op("diffusion",{SF::Test::symbol("mu"),SF::Test::symbol("U")},"momentum.diffusion")),E::constantValue(0)});
    c.addEquation({"energy",E::add(E::add(ddt("rhoE","energy.time"),div("energyFlux","energy.convection")),E::op("diffusion",{SF::Test::symbol("conductivity"),SF::Test::symbol("T")},"energy.diffusion")),E::constantValue(0)});
    int n=0;for(const auto* id:{"continuity","momentum","energy"}) {
        ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.order=n==0?10:n==1?20:60;
        call.step={id,{n==0?"rho":n==1?"rhoU":"rhoE"}};++n;c.addExecution(call);c.bindNumerics({id,"ConservativeResidual"});
    }
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="physicalStep.commit";commit.order=1000000;c.addExecution(commit);
    if(config.solver.turbulence.enabled) {
        for(const auto& id:{"k","omega","mu_t"}) {
            StateSymbol s;s.id=s.name=s.storageKey=id;s.nameSpace="turbulence";s.storageBinding=StorageBinding::ProviderDistributed;
            s.role=std::string(id)=="mu_t"?StateRole::Derived:StateRole::Transported;s.restartEligible=false;
            s.initializationRequired=s.boundaryRequired=std::string(id)!="mu_t";
            if(std::string(id)=="mu_t"){s.evaluation=StateEvaluation::Materialized;s.dependencies={"k","omega","rho","U"};}c.addState(SF::Test::state(s));
        }
        for(const auto* id:{"SST.nu.k","SST.nu.omega","SST.production.k","SST.production.omega","SST.sink.k","SST.sink.omega","SST.cross"}) {
            StateSymbol s;s.id=s.name=s.storageKey=id;s.role=StateRole::Derived;s.storageBinding=StorageBinding::TransientWorkspace;
            s.runtimeStorageRequired=s.initializationRequired=s.boundaryRequired=s.restartEligible=s.outputEligible=false;c.addState(SF::Test::state(s));
        }
        for(const auto& id:std::vector<std::string>{"k","omega"}) {
            const auto symbol=[](const std::string& x){return SF::Test::symbol(x);};
            auto rhs=E::add(E::op("diffusion",{E::multiply(symbol("rho"),symbol("SST.nu."+id)),symbol(id)}),E::subtract(symbol("SST.production."+id),E::multiply(symbol("SST.sink."+id),symbol(id))));
            if(id=="omega")rhs=E::add(rhs,symbol("SST.cross"));
            c.addEquation({id,E::add(E::op("ddt",{E::multiply(symbol("rho"),symbol(id))}),E::op("div",{symbol("rhoU"),symbol(id)})),rhs});
            ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.step={id,{id}};c.addExecution(call);c.bindNumerics({id,"SSTStageTransport",{"k","omega"}});
        }
        c.placement={{{},"equation:k",{"equation:energy"},{"kind:Commit"}},{{},"equation:omega",{"equation:k"},{"kind:Commit"}}};
        c.requireProvider("closure.turbulence","bind turbulence closure and transported equations");
        c.addClosure("mu_t from kOmegaSST",{"k","omega"},{"mu_t"});
    }
    if(!config.solver.sources.gravity.empty())c.extendMathematics({"momentum"},E::op("source",{SF::Test::symbol("gravity")},"gravity"));
    request.userContributions.push_back(c);return build(config.solver,request);
}
std::string signature(const ResolvedSimulationSystem& system) {
    Model::Parameters result=Model::Parameters::object();
    for(const auto& e:system.executableSystem.registry.entries())result["math"][e.id]=canonicalFormula(e);
    for(const auto& s:system.executableSystem.state.symbols())result["state"][s.id]={int(s.role),int(s.storageBinding),s.storageKey,s.componentOffset,s.components,s.dependencies,s.initializationRequired,s.boundaryRequired,s.runtimeStorageRequired};
    for(const auto& c:system.solvePlan.compiledProgram.steps) {
        auto& call=result["calls"][c.source.occurrence];call={{"equation",c.source.equation},{"target",targetText(c.source.target)},{"method",c.equationMethod},{"provider",c.backendProvider},{"old",c.oldTimeWorkspace},{"stage",c.stageWorkspace},{"rhs",c.residualWorkspace}};
        for(const auto& use:c.stateUses)call["reads"].push_back({use.symbol,int(use.version),use.halo});
        for(const auto& effect:c.stateEffects)call["writes"].push_back({effect.symbol,effect.publish,effect.halo});
    }
    for(const auto& binding:system.numericalSystem.operators)result["operators"].push_back({binding.formulaId,binding.occurrence,binding.output,binding.operatorName(),binding.primary,binding.secondary,binding.provider,binding.recipe?FDM::toString(binding.recipe->id()):"provider-defined"});
    for(const auto& c:system.solvePlan.compiledProgram.steps) {
        auto& call=result["calls"][c.source.occurrence];
        call["coefficientStates"]=c.coefficientStates;call["coefficientProvides"]=c.coefficientProvides;
        for(const auto& use:c.coefficientUses)call["coefficientUses"].push_back({use.symbol,use.producer,use.resource,int(use.version)});
    }
    const auto plan=[&](const auto& self,const SolvePlanNode& p)->Model::Parameters {
        Model::Parameters node={{"kind",int(p.kind)},{"id",p.id},{"op",p.operation},{"provider",p.provider},{"occurrence",p.occurrence},{"repeat",p.repetitions},{"minimum",p.minimumIterations},{"termination",p.terminationSignal}};
        for(const auto& c:p.children)node["children"].push_back(self(self,c));return node;
    };result["plan"]=plan(plan,system.solvePlan.root);return nlohmann::json::parse(result.dump()).dump(2);
}
}
int main(int argc,char** argv) {
    if(argc==5 && std::string(argv[1])=="prepare"){prepare(argv[2],argv[3],argv[4]);return 0;}
    check(argc==4,"CASE OUTPUT mode (reference/inspect)");
    auto inspected=Application::inspectCase(argv[1]);auto config=inspected.config;
    std::filesystem::create_directories(argv[2]);
    std::ofstream(std::string(argv[2])+"/signature.json")<<signature(inspected.system);
    if(std::string(argv[3])=="inspect") {
        const auto model=*config.modelDescription;
        const auto first=signature(inspected.system);
        auto inactive=model;inactive.objects.erase(std::remove_if(inactive.objects.begin(),inactive.objects.end(),[](const auto& o){return o.type=="turbulence";}),inactive.objects.end());
        inactive.objects.push_back({"turbulence","models","turbulence",{},{{"enabled",false},{"simulationType","laminar"}}});
        const auto off=Application::inspectCase(Application::ModelLoader::build(inactive));
        check(!off.system.executableSystem.state.contains("k") && !requiresProvider(off.system,"flow.turbulence"),"Consecutive disabled build retained turbulence resources");
        check(off.system.executableSystem.state.contains("mu"),"Shared flow diffusion resource was removed with turbulence");
        if(!inspected.system.executableSystem.state.contains("k"))check(!requiresProvider(inspected.system,"flow.turbulence"),"Inactive model retained transport owner");
        for(int i=0;i<2;++i) {
            auto next=Application::ModelLoader::build(model);
            check(next.modelDescription->selectionOrigins==model.selectionOrigins,"Consecutive builds changed selection provenance");
            check(signature(Application::inspectCase(next).system)==first,"Consecutive builds leaked composition state");
        }
        std::cout<<Application::describeCase(inspected);return 0;
    }
    auto reference=manual(config);validate(reference);
    std::ofstream(std::string(argv[2])+"/manual-signature.json")<<signature(reference);
    check(signature(reference)==signature(inspected.system),"Independent explicit contribution differs: see signatures");
    check(reference.rawSystem.contributions.back().origin.kind==OriginKind::User,"Manual contribution lost independent provenance");
    inspected.system=std::move(reference);
    config.output.caseDir=argv[2];std::filesystem::create_directories(config.output.caseDir+"/result");ResultWriter writer(config.output);Application::ModelLoader::configureOutput(writer,config,&inspected.system.executableSystem.state);
    Application::ExecutionEnvironment env;int count=argc;char** args=argv;
    check(Application::Environment::build(env,config,inspected,inspected.system.runtime,count,args,writer),"Reference resource binding failed");
    check(Application::Environment::validate(env),"Reference environment invalid");
    Application::RunRequest request;request.casePath=argv[1];
    return Application::Execution::execute(inspected.system,inspected.system.solvePlan,config,env,request);
}
