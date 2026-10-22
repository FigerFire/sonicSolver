#include "test/typedSymbols.h"
#include "app/application/SF_inspection.h"
#include "app/application/execution/SF_flowLoop.h"
#include "app/application/output/SF_fields.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemPrinter.h"
#include "solver/system/SF_scalarMethod.h"
#include "core/field/SF_scalarField.h"
#include "infrastructure/execution/SF_localExecutionRuntime.h"
#include "infrastructure/mesh/SF_mesh.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace SF;
int main(int argc,char** argv) {
    if(argc<6 || (argc-3)%3)throw std::runtime_error("reference CASE OUTPUT symbol equation D ...");
    const auto inspected=Application::inspectCase(argv[1]);const auto& config=inspected.config;
    using E=System::FormulaExpr;System::BuildRequest request;request.includeDefaultFluidPreset=false;
    request.composition.stateDeclared=true;
    System::SystemContribution contribution;System::ExecutionProgram how;how.root.id="native.diffusion";
    for(int n=3;n<argc;n+=3) {
        const std::string symbol=argv[n],equation=argv[n+1];
        std::istringstream parameters(argv[n+2]);std::string token;std::getline(parameters,token,':');const double D=std::stod(token);
        double reaction=0,source=0;
        if(std::getline(parameters,token,':'))reaction=std::stod(token);
        if(std::getline(parameters,token,':'))source=std::stod(token);
        request.composition.solutionVariables.push_back(symbol);
        System::StateSymbol state;state.id=state.name=state.storageKey=symbol;state.storageBinding=System::StorageBinding::NamedDistributed;state.restartEligible=false;
        contribution.addState(([&]{state.signature=SF::Math::signature(0,{});return state;})());
        auto rhs=E::op("diffusion",{E::constantValue(D,SF::Math::signature(0,{0,2,-1,0,0,0,0})),E::symbol(symbol,SF::Math::signature(0,{}))});
        if(reaction)rhs=E::add(std::move(rhs),E::multiply(E::constantValue(reaction,SF::Math::signature(0,{0,0,-1,0,0,0,0})),E::symbol(symbol,SF::Math::signature(0,{}))));
        if(source)rhs=E::add(std::move(rhs),E::constantValue(source,SF::Math::signature(0,{0,0,-1,0,0,0,0})));
        contribution.addEquation(System::Equation{equation,E::op("ddt",{E::symbol(symbol,SF::Math::signature(0,{}))}),std::move(rhs)});
        System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;call.step={equation,{symbol,System::TargetKind::Physical},"advance."+symbol};how.root.children.push_back(call);
        request.authoredNumerics.push_back({equation,"ScalarDiffusionCentral2"});
    }
    System::ExecutionScope commit;commit.kind=System::ExecutionKind::Commit;commit.id="physical.commit";how.root.children.push_back(commit);
    request.authoredExecution=how;request.userContributions.push_back(contribution);
    auto reference=System::build(config.solver,request);
    const auto& native=inspected.system;
    if(reference.executableSystem.registry.entries().size()!=native.executableSystem.registry.entries().size())throw std::runtime_error("Native/C++ equation count mismatch");
    for(const auto& e:reference.executableSystem.registry.entries())
        if(System::canonicalFormula(e)!=System::canonicalFormula(native.executableSystem.registry.at(e.id)))throw std::runtime_error("Native/C++ mathematical IR mismatch");
    const auto& a=reference.solvePlan.compiledProgram.steps;const auto& b=native.solvePlan.compiledProgram.steps;
    if(a.size()!=b.size())throw std::runtime_error("Native/C++ call count mismatch");
    for(std::size_t i=0;i<a.size();++i)
        if(a[i].source.equation!=b[i].source.equation || a[i].source.occurrence!=b[i].source.occurrence || a[i].target.symbol!=b[i].target.symbol || a[i].backendProvider!=b[i].backendProvider || a[i].stageWorkspace!=b[i].stageWorkspace || a[i].residualWorkspace!=b[i].residualWorkspace || a[i].oldTimeWorkspace!=b[i].oldTimeWorkspace)
            throw std::runtime_error("Native/C++ compiled binding mismatch");
    Field geometry;if(!Mesh::setupComplexMesh(geometry,config,inspected.mesh))throw std::runtime_error("Reference mesh loading failed");
    if(geometry.NVar() || geometry.stateModel())throw std::runtime_error("Reference acquired fake Q/EOS");
    State::StateBundle state;state.patches={&geometry};state.time=config.time.startTime;
    std::vector<std::unique_ptr<ScalarField>> arrays;std::vector<std::shared_ptr<FDM::ScalarEquationData>> ports;
    FDM::SolverServices services;Execution::LocalRuntime runtime;services.executionRuntime=&runtime;
    for(int n=3;n<argc;n+=3) {
        const int slot=(n-3)/3;const std::string name=argv[n];
        auto physical=std::make_unique<ScalarField>();physical->setupLike(geometry,name);
        for(int k=0;k<geometry.MZ();++k)for(int j=0;j<geometry.MY();++j)for(int i=0;i<geometry.MX();++i) {
            const bool left=i<=geometry.NG()+16;
            (*physical)(i,j,k)=slot==0?(left?1:0):slot==1?(left?.5:2):(left?.2:.75);
        }
        state.distributed.add(State::scalarView(name,0,geometry,*physical));arrays.push_back(std::move(physical));
        auto data=std::make_shared<FDM::ScalarEquationData>();data->target=name;
        for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
        if(slot==0)for(int side=0;side<2;++side){data->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;data->boundary[side].value=[side](double,double,double,double){return side?0:1;};}
        services.scalarInstances.push_back({"advance."+name,data});ports.push_back(std::move(data));
    }
    auto output=Application::Output::registeredStateFields(reference.executableSystem.state,state,geometry);
    std::sort(output.begin(),output.end(),[](const auto& x,const auto& y){return x.name<y.name;});
    auto writeConfig=config.output;writeConfig.caseDir=argv[2];std::filesystem::create_directories(writeConfig.caseDir+"/result");ResultWriter writer(writeConfig);
    auto save=[&](int step){writer.saveStep(geometry,step,state.time,output);};save(0);
    Time::DriverCallbacks callbacks;callbacks.saveStep=save;
    Application::Execution::executeEquations(config.solver,reference,state,services,config.time,callbacks);
    if(state.step!=10 || state.time!=.005)throw std::runtime_error("Reference clock mismatch");
    std::ofstream(std::string(argv[2])+"/explain.txt")<<System::describe(reference);
    std::cout<<"Native/C++ canonical WHAT and compiled bindings identical; ten steps completed without Q/EOS\n";
}
