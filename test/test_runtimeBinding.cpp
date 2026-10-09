#include "solver/algorithm/transport/SF_transportBinding.h"
#include "solver/algorithm/SF_singleFluidStepper.h"
#include <type_traits>
#include "solver/algorithm/interface/SF_levelSetBinding.h"
#include "solver/system/SF_providerResolver.h"
#include "solver/system/SF_systemBuilder.h"
#include "core/system/SF_operationIds.h"
#include <map>
#include <stdexcept>
#include <vector>
using namespace SF;
namespace {
void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
struct Transport : FDM::ITransportModel {
    std::vector<std::string> events;
    std::vector<double> steps;
    void applyBoundary(const Field&) override { events.push_back("boundary"); }
    void correct(const Field&,double dt) override { events.push_back("correct");steps.push_back(dt); }
    double dynamicViscosity(const Field&,int,int,int,double mu) const override { return mu; }
};
struct Interface : FDM::ILevelSetOperations {
    Field& geometry;
    double phi=2,reference=0;
    std::vector<double> references,physicalDt;
    int stages=0;
    explicit Interface(Field& field):geometry(field) {}
    void validateAdvection(int,double,double,bool,bool,double,double) const override {}
    void beginReinitialization(int order,double dt,double,double,double) override {
        require(order==5 && dt==.01,"Frozen pseudo-time recipe changed");
        reference=phi;references.push_back(reference);
    }
    State::DistributedFieldView referenceView() override {
        State::DistributedFieldView view;view.name="phi0";view.geometry=&geometry;view.components=1;view.blockId=0;
        view.read=[this](int,int) { return reference; };
        view.write=[](int,int,double) { throw std::runtime_error("Frozen reference is read-only"); };
        return view;
    }
    void reinitializeStage() override { ++stages;phi+=1; }
    void publishGeometry(double dt) override { physicalDt.push_back(dt); }
};
System::SolvePlanNode leaf(const char* operation,const char* provider) {
    System::SolvePlanNode node;node.kind=System::PlanNodeKind::Update;node.operation=operation;node.provider=provider;return node;
}
}
int main() {
    static_assert(!std::is_copy_constructible_v<SolverAlgorithm::SingleFluidStepper>);
    static_assert(!std::is_move_constructible_v<SolverAlgorithm::SingleFluidStepper>);
    Field field;field.setup(1,1,1,0,5);
    State::StateBundle state;state.patches={&field};state.time=12;state.step=99;state.dt=.1;
    FDM::SolverServices services;Transport model;services.transportModel=&model;
    System::RuntimeRequirements requirements;
    requirements.operationBindings.push_back({System::OpIds::TurbulenceAdvance,"flow.turbulence",System::BindingStatus::Resolved,{}});
    Run::OpRegistry operations;
    SolverAlgorithm::bindTransportOperations(operations,requirements,state,services);
    Run::ExecutionContext execution;
    operations.invoke(System::OpIds::TurbulenceAdvance,execution,"flow.turbulence");
    state.dt=.2;
    operations.invoke(System::OpIds::TurbulenceAdvance,execution,"flow.turbulence");
    require(model.steps==std::vector<double>{.1,.2},"Long-lived binding captured a stale timestep");
    require(model.events==std::vector<std::string>{"boundary","correct","boundary","boundary","correct","boundary"},"Transport BC/correction order changed");
    require(state.time==12 && state.step==99,"Transport gained clock authority");
    auto missing=services;missing.transportModel=nullptr;
    bool rejected=false;
    try {Run::OpRegistry invalid;SolverAlgorithm::bindTransportOperations(invalid,requirements,state,missing);}
    catch(const std::runtime_error& e) { rejected=std::string(e.what()).find("flow.turbulence")!=std::string::npos && std::string(e.what()).find(System::OpIds::TurbulenceAdvance)!=std::string::npos; }
    require(rejected,"Missing provider port did not identify its operation");

    Interface interface(field);services.levelSet=&interface;
    System::CompiledSolvePlan plan;plan.root.kind=System::PlanNodeKind::Sequence;
    plan.root.children.push_back(leaf(System::OpIds::LevelSetReference,"interface.level-set"));
    System::SolvePlanNode loop;loop.kind=System::PlanNodeKind::Loop;loop.id="pseudoTime";loop.repetitions=3;
    loop.children.push_back(leaf(System::OpIds::LevelSetPseudoStage,"interface.level-set"));
    plan.root.children.push_back(loop);plan.root.children.push_back(leaf(System::OpIds::LevelSetGeometry,"interface.level-set"));
    for (const auto* id:{System::OpIds::LevelSetReference,System::OpIds::LevelSetPseudoStage,System::OpIds::LevelSetGeometry}) {
        System::CompiledEquationCall call;call.backendProvider="interface.level-set";call.backendOperation=id;
        if (id==System::OpIds::LevelSetReference)
            call.providerContract=std::map<std::string,double>{{"order",5},{"pseudoDt",.01},{"epsilon",1e-6},{"power",2},{"signFactor",1}};
        plan.compiledProgram.steps.push_back(call);
    }
    System::CompiledStateView demand;demand.symbol="phi0";demand.storage="phi0";
    demand.kind=System::StateViewKind::Workspace;demand.owner=System::StateViewOwner::NumericalProvider;
    System::StateRegistry symbols;
    auto realized=System::realizeState(symbols,{},state,{demand});
    Run::OpRegistry interfaceOps;
    SolverAlgorithm::bindLevelSetOperations(interfaceOps,plan,state,services,realized);
    Run::PlanExecutor::validateBindings(plan,interfaceOps);
    Run::PlanExecutor::execute(plan,interfaceOps);
    require(interface.stages==3 && interface.references==std::vector<double>{2},"Binder hid a pseudo-time loop or froze repeatedly");
    require(realized.view("phi0",System::StateViewKind::Workspace).fields.front()->read(0,0)==2,"phi0 lost its original frozen backing");
    state.dt=.3;Run::PlanExecutor::execute(plan,interfaceOps);
    require(interface.stages==6 && interface.references==std::vector<double>{2,5},"Next physical step reused the old reference");
    require(interface.physicalDt==std::vector<double>{.2,.3},"Geometry captured stale physical dt");
    require(state.time==12 && state.step==99,"Pseudo-time changed the physical clock");
    auto wrong=plan;wrong.root.children.front().provider="wrong.port";rejected=false;
    try {Run::PlanExecutor::validateBindings(wrong,interfaceOps);} catch(const std::runtime_error& e) {
        rejected=std::string(e.what()).find("wrong.port")!=std::string::npos;
    }
    require(rejected,"An available operation incorrectly satisfied a different frozen provider");

    FDM::SolverConfig config;System::BuildRequest request;
    request.composition.stateDeclared=true;request.composition.solutionVariables={"rho","rhoU","rhoE"};
    request.singleFluidPreset=System::SingleFluidPresetSpec{false};
    auto compiled=System::build(config,request);
    require(System::validateCompiledConservativeStage(compiled.executableSystem,compiled.solvePlan).empty(),"Valid compiled stage failed shared validation");
    const auto corrupt=[&](const auto& self,System::SolvePlanNode& node)->void {
        if (node.operation==System::OpIds::ExplicitStageExecute) node.equationCalls.front().target="badTarget";
        for (auto& child:node.children) self(self,child);
    };
    corrupt(corrupt,compiled.solvePlan.root);
    require(!System::validateCompiledConservativeStage(compiled.executableSystem,compiled.solvePlan).empty(),"Runtime/compiler shared guard missed corrupted fusion");
    auto bindings=System::compileOperationBindings(compiled.executableSystem,compiled.numericalSystem,compiled.solvePlan,{},false);
    require(System::reportOperationBindings(compiled.solvePlan,bindings).status==System::RuntimeStatus::Unsupported,"Static report falsely advertised corrupted fusion as Runnable");
}
