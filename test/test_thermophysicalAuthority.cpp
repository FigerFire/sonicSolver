#include "app/application/model/SF_model.h"
#include "infrastructure/io/case/SF_case.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "solver/boundary/SF_applicator.h"
#include "methods/numerics/viscous/SF_viscous.h"
#include "app/application/model/SF_runtimeConfig.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include "app/application/SF_inspection.h"
#include "solver/system/SF_systemPrinter.h"
#include "models/initial/SF_init.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <filesystem>
using namespace SF;
namespace {
void require(bool x,const char* m){if(!x)throw std::runtime_error(m);}
Model::Description input() {
    auto m=CaseIO::read(std::string(SF_TEST_SOURCE_DIR)+"/test/Sod/sodCase_weno7_t0p2");
    for(auto& o:m.objects)if(o.type=="thermoDynamics") {
        o.parameters={{"thermoDynamics",{{"equationOfState","perfectGas"},{"thermo","hConst"},{"transport","const"}}},
            {"properties",{{"equationOfState",{{"gamma",1.32},{"R",310.}}},{"transport",{{"mu",2.1e-5},{"Pr",.81}}}}}};
    }
    m.numerics.erase("transport");return m;
}
Model::Parameters& properties(Model::Description& m) {
    for(auto& o:m.objects)if(o.type=="thermoDynamics"||o.type=="thermophysical"||o.type=="Thermophysical")return o.parameters;
    throw std::runtime_error("missing thermophysical test object");
}
void rejected(Model::Description m,const std::string& token) {
    std::cout<<"negative: "<<token<<std::endl;
    try{(void)Application::ModelLoader::build(m);}catch(const std::exception& e){
        require(std::string(e.what()).find(token)!=std::string::npos,e.what());return;}
    throw std::runtime_error("Invalid thermophysical input accepted: "+token);
}
void ioChecks(const Model::Description& m) {
    namespace fs=std::filesystem;
    const auto base=fs::path(SF_TEST_SOURCE_DIR)/"test/t/thermophysical-io";
    fs::create_directories(base);
    auto sameContract=[](const auto& a,const auto& b) {
        require(a.gamma==b.gamma&&a.gasConstant==b.gasConstant&&a.dynamicViscosity==b.dynamicViscosity
            &&a.prandtl==b.prandtl&&a.selection.equationOfState==b.selection.equationOfState
            &&a.selection.thermo==b.selection.thermo&&a.selection.transport==b.selection.transport
            &&a.native==b.native,"IO changed frozen parameters or selection provenance");
    };
    CaseIO::write(m,(base/"roundtrip").string());
    fs::copy_file(fs::path(m.caseDir)/"mesh/mesh.sfm",base/"roundtrip/mesh/mesh.sfm",fs::copy_options::overwrite_existing);
    const auto a=Application::ModelLoader::build(m);
    const auto b=Application::ModelLoader::read((base/"roundtrip").string());
    sameContract(*a.solver.thermophysical,*b.solver.thermophysical);
    const auto ra=Application::inspectCase(a),rb=Application::inspectCase(b);
    require(ra.system.runtime.report.requiredOperations==rb.system.runtime.report.requiredOperations,"IO changed execution operations");
    const auto explain=System::describe(rb.system);
    require(explain.find("gamma=1.32 R=310")!=std::string::npos&&explain.find("mu=2.1e-05 Pr=0.81")!=std::string::npos,"explain omitted actual parameters");
    auto rejectWrite=[&](Model::Description bad,const std::string& token) {
        std::cout<<"IO negative: "<<token<<std::endl;
        const auto output=base/"invalid";fs::create_directories(output);
        {std::ofstream file(output/"case.yaml");file<<"sentinel\n";}
        bool rejected=false;try{CaseIO::write(bad,output.string());}catch(const std::runtime_error& e){rejected=std::string(e.what()).find(token)!=std::string::npos;}
        require(rejected,"invalid output paths accepted");
        std::ifstream file(output/"case.yaml");std::string line;std::getline(file,line);require(line=="sentinel","preflight overwrote case.yaml");
        require(!fs::exists(output/"solvers/runtime.yaml"),"preflight partially wrote files");
    };
    auto bad=m;bad.objects.push_back({"a.b","models","custom",{},Model::Parameters::object(),{}});
    bad.objects.push_back({"a/b","models","custom",{},Model::Parameters::object(),{}});rejectWrite(bad,"collision");
    bad=m;bad.objects.push_back({"Name","models","custom",{},Model::Parameters::object(),{}});
    bad.objects.push_back({"name","models","custom",{},Model::Parameters::object(),{}});rejectWrite(bad,"collision");
    bad=m;bad.objects.push_back({"runtime","solver","custom",{},Model::Parameters::object(),{}});rejectWrite(bad,"reserved");
    bad=m;bad.objects.push_back({"models","models","custom",{},Model::Parameters::object(),{}});rejectWrite(bad,"collision");
    bad=m;bad.objects.push_back(m.objects.front());rejectWrite(bad,"Duplicate");
    auto separate=m;separate.objects.push_back({"same","models","custom",{},Model::Parameters::object(),{}});
    separate.objects.push_back({"same","solver","custom",{},Model::Parameters::object(),{}});
    CaseIO::write(separate,(base/"separate").string());
    require(fs::exists(base/"separate/models/same.yaml")&&fs::exists(base/"separate/solvers/same.yaml"),"same names in different categories collided");
    require(CaseIO::read((base/"separate").string()).objects.size()==separate.objects.size(),"separate objects lost");
    auto rejectRead=[&](const std::string& token) {
        bool rejected=false;try{(void)CaseIO::read((base/"roundtrip").string());}catch(const std::runtime_error& e){rejected=std::string(e.what()).find(token)!=std::string::npos;}
        require(rejected,"invalid referenced file accepted");
    };
    const auto regPath=base/"roundtrip/models/models.yaml";
    const auto reg=Serialization::readSonicFile(regPath.string());auto changed=reg;
    const auto entry=changed.body.begin().value();changed.body["duplicate"]=entry;
    Serialization::writeSonicFile(changed,regPath.string());rejectRead("duplicate referenced file");
    changed=reg;changed.body.begin().value()["file"]="models/missing.yaml";
    Serialization::writeSonicFile(changed,regPath.string());rejectRead("Cannot read");
    changed=reg;changed.body.begin().value()["type"]="wrongType";
    Serialization::writeSonicFile(changed,regPath.string());rejectRead("expected object");
    Serialization::writeSonicFile(reg,regPath.string());
    const auto modelPath=base/"roundtrip"/reg.body.begin().value().at("file").get<std::string>();
    auto doc=Serialization::readSonicFile(modelPath.string());doc.object="solver";
    Serialization::writeSonicFile(doc,modelPath.string());rejectRead("expected object");
    fs::remove_all(base);
}

}
int main(){try{
    const auto m=input();const auto c=Application::ModelLoader::build(m);const auto& t=*c.solver.thermophysical;
    require(t.native&&t.gamma==1.32&&t.gasConstant==310.&&t.dynamicViscosity==2.1e-5&&t.prandtl==.81,"native parameters lost");
    require(std::abs(t.cp()-1278.75)<1e-10&&std::abs(t.cv()-968.75)<1e-10,"caloric relations failed");
    require(c.solver.numerics.idealGasGamma==t.gamma&&c.solver.numerics.dynamicViscosity==t.dynamicViscosity
        &&c.solver.turbulence.laminarDynamicViscosity==t.dynamicViscosity&&c.solver.boundaries.thermalPrandtl==t.prandtl,"projection differs");
    auto boundaryOnly=m;
    Model::FieldDescriptor temperature;
    temperature.name="T";
    temperature.boundaries={{"Right",{{"type","fixedTemperature"},{"value",300.}}}};
    boundaryOnly.fields.push_back(temperature);
    const auto boundaryOnlyConfig=Application::ModelLoader::build(boundaryOnly);
    require(boundaryOnlyConfig.solver.initial.temperature.empty()
        &&boundaryOnlyConfig.solver.boundaries.thermal.size()==1
        &&boundaryOnlyConfig.solver.boundaries.thermal.front().value==300.,
        "temperature boundary law was discarded when initialization used pressure");
    Field fromP,fromT;fromP.setup(2,2,2,1,5);fromT.setup(2,2,2,1,5);
    const int point=fromP.getIdx(1,1,1);
    fromP.setBoundarySets({{"all",{point}}});fromT.setBoundarySets({{"all",{point}}});
    FDM::InitialConditionConfig initial;
    initial.density={{"all",FIXED_VALUE,1.225}};initial.velocity={{"all",FIXED_VALUE,{2.36643,0.,0.}}};
    initial.pressure={{"all",FIXED_VALUE,113925.}};
    Init::setupAll(fromP,initial,t.gamma,t.gasConstant);
    initial.temperature={{"all",FIXED_VALUE,300.}};
    Init::setupAll(fromT,initial,t.gamma,t.gasConstant);
    std::cout<<std::setprecision(17)<<"initial raw Q: p->rhoE="<<fromP(1,1,1,E)
        <<" T->rhoE="<<fromT(1,1,1,E)<<" delta="<<fromP(1,1,1,E)-fromT(1,1,1,E)<<std::endl;
    require(std::abs(fromP(1,1,1,E)-fromT(1,1,1,E))<1.e-9,"p/T EOS initialization differs");
    const auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(t);
    Math::resetActiveDirections();Math::deactivateDirection(2);
    Field f;f.setup(15,15,2,3,5);f.setStateModel(model);
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
        f.X(i,j,k)=i;f.Y(i,j,k)=j;f.Z(i,j,k)=k;
        f.Jac(i,j,k)=1.;f.XiX(i,j,k)=1.;f.EtY(i,j,k)=1.;f.ZeZ(i,j,k)=1.;
        f(i,j,k,RHO)=1.2;f(i,j,k,RU)=0.;f(i,j,k,RV)=0.;f(i,j,k,RW)=0.;
        Boundary::setEnergyFromTemperature(f,i,j,k,300.+2.*i);
    }
    auto equal=[](double a,double b){require(std::abs(a-b)<1.e-10*std::max(1.,std::abs(b)),"actual thermophysical kernel differs");};
    const auto state=f.thermodynamicState(10,10,3);
    equal(state.temperature,320.);equal(state.pressure,1.2*310.*320.);
    equal(state.soundSpeed,std::sqrt(1.32*310.*320.));
    equal(state.internalEnergyDensity,1.2*t.cv()*320.);
    equal(state.dynamicViscosity,t.dynamicViscosity);equal(state.thermalConductivity,t.conductivity());
    Viscous::TransportProperties transport{t.dynamicViscosity,t.prandtl,t.gamma,t.gasConstant};
    std::vector<double> flux;
    Viscous::viscousFluxAtFace(f,10,10,3,Math::XI,Viscous::CentralOrder::SECOND,transport,flux);
    equal(flux[E],2.*t.conductivity());
    f.setBoundarySets({{"Top",{f.getIdx(10,17,3)}}});
    Boundary::updateEnergyFromThermalBoundary(f,{{"Top",ThermalBCType::HeatFlux,2.*t.conductivity()}},0.,.72,true,3);
    f.invalidateThermodynamicCache();equal(Boundary::temperatureAt(f,10,17,3),322.);
    Boundary::updateEnergyFromPressure(f,{{"Top",FIXED_VALUE,101000.}},true,3);
    f.invalidateThermodynamicCache();equal(Boundary::pressureAt(f,10,17,3),101000.);
    equal(f(10,17,3,E),101000./.32);
    const auto mesh=Application::Runtime::makeMeshRuntimeConfig(c,c.solver,3);
    equal(mesh.idealGasGamma,t.gamma);equal(mesh.idealGasConstant,t.gasConstant);
    // A shear flux must consume the same molecular mu, independently of heat flux.
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
        f(i,j,k,RV)=1.2*(2.*i);
        Boundary::setEnergyFromTemperature(f,i,j,k,300.+2.*i);
    }
    f.invalidateThermodynamicCache();
    Viscous::viscousFluxAtFace(f,10,10,3,Math::XI,Viscous::CentralOrder::SECOND,transport,flux);
    equal(flux[RV],2.*t.dynamicViscosity);
    equal(flux[E],21.*2.*t.dynamicViscosity+2.*t.conductivity());
    Physics::FluidStateModel::initializeSingleFluidRestart(f,model);
    equal(f.thermodynamicState(10,10,3).temperature,320.);
    auto drift=c.solver;drift.numerics.idealGasGamma=1.4;
    bool driftRejected=false;try{drift.validateThermophysicalProjection();}catch(const std::runtime_error& e){driftRejected=std::string(e.what()).find("projection drift")!=std::string::npos;}
    require(driftRejected,"changed kernel projection accepted");
    Field unbound;unbound.setup(2,2,2,1,5);
    bool unboundRejected=false;try{Boundary::Applicator(c.solver.boundaries).apply(unbound);}catch(const std::runtime_error& e){unboundRejected=std::string(e.what()).find("attached thermodynamic")!=std::string::npos;}
    require(unboundRejected,"native boundary accepted fixed-gamma fallback");
    for(auto alias:{"thermophysical","Thermophysical"}){auto a=m;for(auto& o:a.objects)if(o.type=="thermoDynamics")o.type=alias;
        require(Application::ModelLoader::build(a).solver.thermophysical->gamma==1.32,"alias bypassed parameters");}
    auto a=m;properties(a)["properties"]["equationOfState"].erase("gamma");rejected(a,"gamma");
    a=m;properties(a)["properties"].erase("transport");rejected(a,"properties.transport");
    for(double gamma:{1.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
        a=m;properties(a)["properties"]["equationOfState"]["gamma"]=gamma;rejected(a,gamma<=1?"gamma>1":"nonfinite");}
    a=m;properties(a)["properties"]["equationOfState"]["R"]=0.;rejected(a,"R>0");
    a=m;properties(a)["properties"]["transport"]["mu"]=-1.;rejected(a,"mu>=0");
    a=m;properties(a)["properties"]["transport"]["Pr"]=0.;rejected(a,"Pr>0");
    a=m;properties(a)["properties"]["thermo"]={{"cp",1000.},{"cv",700.}};rejected(a,"inconsistent cp/cv");
    a=m;properties(a)["thermoDynamics"]["thermo"]="janaf";rejected(a,"Unsupported caloric provider");
    a=m;properties(a)["thermoDynamics"]["transport"]="sutherland";rejected(a,"Unsupported transport provider");
    a=m;properties(a)["thermoDynamics"]["equationOfState"]="madeUp";rejected(a,"unknown equationOfState");
    a=m;a.numerics["transport"]={{"mu",1e-3}};rejected(a,"parameter conflict: mu");
    a=m;for(const auto& o:m.objects)if(o.type=="thermoDynamics"){auto copy=o;copy.name="duplicate";a.objects.push_back(copy);}rejected(a,"Duplicate thermophysical authority");
    a=m;properties(a)["thermoDynamics"]["energy"]="rhoE";rejected(a,"energy is WHAT");
    a=m;properties(a)["thermoDynamics"]={{"equationOfState","rhoConst"},{"transport","const"}};
    properties(a)["properties"]["equationOfState"]={{"rho",1.2}};
    auto r=Application::ModelLoader::build(a);require(r.solver.thermophysical->selection.constantDensity==1.2
        &&std::isnan(r.solver.thermophysical->gamma),"rhoConst invented PerfectGas");
    a=m;properties(a)["thermoDynamics"].erase("transport");properties(a)["properties"].erase("transport");
    a.numerics["terms"]["diffusion"]="central2Explicit";rejected(a,"selected transport provider");
    a=m;for(auto& o:a.objects)if(o.type=="stateRegistry")o.parameters["use"]=Model::Parameters::array({"U","p"});
    bool unsupported=false;try{(void)Application::inspectCase(Application::ModelLoader::build(a));}catch(const std::runtime_error&){unsupported=true;}
    require(unsupported,"unsupported energy/state combination accepted");
    ioChecks(m);
    std::cout<<"Thermophysical typed authority and negative cases passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
